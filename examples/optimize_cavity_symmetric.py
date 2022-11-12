"""
This example optimizes a cavity from a given starting geometry.
"""

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.utils import BBox
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.engine import LumericalEngine
import os
from scipy.optimize import minimize
import numpy as np
import sividl.sividl_devices as sivp


#  Initialize Lumerical File Locations
FDTDLoc = '/n/home08/eknall/sw_ENK/lumerical-2021-R2-2717-7bf43e7149_seas'
FDTDexeLoc = os.path.join(FDTDLoc,'bin/fdtd-solutions')
FDTDmpiLoc = os.path.join(FDTDLoc,'bin/fdtd-engine-ompi-lcl')
print(FDTDexeLoc)

def fitness(cavity_params):
    aL, aR, hx, hy, maxDef, beam_w = cavity_params
    print(cavity_params)
    pcc_params = {
        'layer'               : 2,
        'aL'                  : aL,
        'aR'                  : aR,
        'hxL'                 : hx,
        'hyL'                 : hy,
        'hxR'                 : hx,
        'hyR'                 : hy,
        'maxDef'              : maxDef,
        'nholesLMirror'       : 4,
        'nholesRMirror'       : 4,
        'nholes_wvg-mirr_trans_L': 5,
        'nholes_wvg-mirr_trans_R': 5,
        'nholes_defect'       : 7,
        'min_hole_dim'        : 0.05,
        'effective_index'     : 1.6,
        'resonance_wavelength': 0.737
    }

    # beam_w = 0.482e-6
    #center_cell = 16
    n_beam = 2.4028
    apex_half_angle = 50*np.pi/180


    cavity = sivp.AirholePCC_PeriodOnly(pcc_params)
    all_hx, all_hy, all_a = cavity.compute_geometry()
    all_a = np.append(all_a,all_a[-1])

    # convert from microns to meters
    all_hx *= 1e-6
    all_hy *= 1e-6
    all_a *= 1e-6
    beam_w *= 1e-6
    beam_h = (beam_w / 2) * np.tan(np.pi/2 - apex_half_angle)

    # Use level 1 automeshing accuracy, and don't show the Lumerical GUI while running simulations
    engine = LumericalEngine(mesh_accuracy=1, hide=False, lumerical_path=FDTDLoc,save_fsp=True,working_path='./fsps')

    cavity_cells = []
    for hx, hy, a in zip(all_hx, all_hy, all_a):
        print(hx,hy,a)
        cell_size = Vec3(a,beam_w,beam_h)
        cell_box = TriStructure(Vec3(0), Vec3(beam_w, apex_half_angle, a), 
                                DielectricMaterial(2.4028, order=2, color="blue"), 
                                rot_angles=(np.pi/2, np.pi/2, 0))

        #cell_box = BoxStructure(Vec3(0), cell_size, DielectricMaterial(n_beam, order=2, color="blue"))
        # offset the hole to respect the way we define the relevant lattice constant
        cell_hole = CylinderStructure(Vec3(-a/2,0,0), beam_h, hx/2, DielectricMaterial(1, order=1, color="red"), radius2=hy/2)
        unit_cell = UnitCell(structures = [cell_box, cell_hole], size=cell_size, engine=engine)
        cavity_cells += [unit_cell]

    # The length of the cavity beam
    beam_length = 10e-6
    # The target resonance frequency, in Hz
    target_frequency = 406.7e12

    cavity = Cavity1D(
      unit_cells=cavity_cells,
      structures=[TriStructure(Vec3(0), Vec3(beam_w, apex_half_angle, beam_length), 
                  DielectricMaterial(2.4028, order=2, color="gray"), rot_angles=(np.pi/2, np.pi/2, 0))], engine=engine)

    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity_name = np.array2string(cavity_params,prefix='',separator='_')
    cavity_name = cavity_name[1:]
    cavity_name = cavity_name[:-1]
    cavity.save(log_name[:-4]+cavity_name+".obj")

    
    man_mesh = MeshRegion(BBox(Vec3(0),Vec3(7e-6,0.7e-6,0.4e-6)), 12e-9, dy=None, dz=None)

    #r3 = cavity.simulate("guidedness", target_freq=target_frequency)
    r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, 
                        analyze_time=1000e-15, mesh_regions = [man_mesh], sim_size=Vec3(2, 4, 8))

    qx = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    qy = 1/(2/r1["qymax"])
    qz = 1/(1/r1["qzmin"] + 1/r1["qzmax"])
    qscat = 1/((1/qy)+(1/qz))
    qtot = 1/(1/qx + 1/qy + 1/qz)
    vmode = r1["vmode"]
    vmode = 0.3 if vmode < 0.3 else vmode
    purcell = qtot/vmode
    F = r1["freq"]

    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
      r1["freq"], r1["vmode"],
      1/(1/r1["qxmin"] + 1/r1["qxmax"]),
      1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))
    print(purcell)
    wavelen_pen = np.exp(-((target_frequency-F)/4e12)**2)
    witness = -1*purcell*wavelen_pen
    with open(log_name, "ab") as f:
        f.write(b"\n")
        step_info = np.append(cavity_params,np.array([witness, wavelen_pen,purcell,r1["qxmin"],r1["qxmax"],qscat,qtot,vmode,F]))
        np.savetxt(f, step_info.reshape(1, step_info.shape[0]), fmt='%.4f')
    #r1["xyprofile"].show()
    #r1["yzprofile"].show()


    return witness
    
log_name = "optimal_sym-tri_090722_00.txt"
#p0 = np.array([0.270,0.250,0.114,0.161,0.1392,0.482])
p0 = np.array([0.270,0.270,0.152,0.162,0.098264,0.474])
bounds = ((0.200,0.300),(0.200,0.300),(0.100,0.200),(0.100,0.200),(0.08,0.18),(0.25,0.52))


with open(log_name, "ab") as f:
    f.write(b"aL    aR    hx    hy    maxDef    beam_w    fitness    wavelen_pen    purcell    qxmin    qxmax    qscat    qtot    vmode    freq")


popt = minimize(fitness,p0,bounds = bounds,method='Nelder-Mead')
print(popt)
# Print the reults and plot the electric field profiles
#print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
#  r1["freq"], r1["vmode"],
#  1/(1/r1["qxmin"] + 1/r1["qxmax"]),
#  1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
#))



#print("Guidedness: %s" % r3)
