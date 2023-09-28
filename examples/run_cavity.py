"""
This example optimizes a cavity from a given starting geometry.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sividl.sividl_devices as sivp

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.utils import BBox
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.engine import LumericalEngine


#  Initialize Lumerical File Locations
FDTDLoc = 'D:\\Program Files\\Lumerical\\v232'  # '/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas'
fsps_dir = './fsps'

nmirrsL = 7
nmirrsR = 5
ndefs = 5
target_qx = 6000

# Target resonance frequency (Hz)
target_frequency = 406.7e12
# Frequency of excitation light in simulation (Hz)
source_frequency = 406.7e12
# Threshold for how far the frequency can be off before we rerun
rerun_thresh = 0.955

n_beam = 2.4028
iter_count = 0
hide = False

def gauss(x, x0, sigma):
    return np.exp(-((x - x0) / sigma) ** 2)


def plot_geom(cavity, file_name, hide):
    plt.close()
    unit_cells = cavity.get_unit_cells()

    a_s = np.zeros(len(unit_cells))
    hx = np.zeros(len(unit_cells))
    hy = np.zeros(len(unit_cells))
    for i,unit_cell in enumerate(unit_cells):
        size = unit_cell.get_size()
        structures = unit_cell.get_structures()
        hx[i] = structures[1].radius*2
        hy[i] = structures[1].radius2*2
        a_s[i] = size.x

    x_coords = np.copy(a_s)
    x_coords = x_coords[:-1]
    x_coords = np.insert(x_coords,0,0)
    x_coord = 0
    for i,a in enumerate(x_coords):
        x_coord += a
        x_coords[i] = x_coord

    a_coords = np.array([(x+x_coords[i+1])/2 for i, x in enumerate(x_coords[:-1])])

    plt.close()

    plt.vlines(a_coords, 0.18, 1.01, color='gray', alpha=0.3)
    plt.hlines(1, np.min(x_coords), np.max(x_coords), colors='red', linestyle='dashed')
    plt.plot(a_coords, a_s[:-1] / a_s[5], 'ko-', label='a')
    plt.plot(x_coords, hx/a_s[5], 'go-', label='hx')
    plt.plot(x_coords, hy/a_s[5], 'bo-', label='hy')

    plt.ylabel(r"$\mathregular{a}_{\mathregular{nominal}}$",fontsize=14)
    plt.xlabel(r"Unit Cell Position ($\mu$m)",fontsize=12)
    plt.legend(loc=8)

    fig = plt.gcf()
    fig.set_size_inches(6, 3)
    plt.tight_layout()
    plt.savefig(file_name)
    if not hide:
        plt.show()
    plt.close()

def build_cavity(cavity_params):
    global iter_count
    iter_count += 1
    maxDef, beam_w, hxL, hyL, hxR, hyR, aL, aR = cavity_params
    print(cavity_params)
    pcc_params = {
        'layer'               : 2,
        'aL'                  : aL,
        'aR'                  : aR,
        'hxL'                 : hxL,
        'hyL'                 : hyL,
        'hxR'                 : hxR,
        'hyR'                 : hyR,
        'maxDef'              : maxDef,
        'nholesLMirror'       : nmirrsL,
        'nholesRMirror'       : nmirrsR,
        'nholes_wvg-mirr_trans_L': 5,
        'nholes_wvg-mirr_trans_R': 5,
        'nholes_defect'       : ndefs,
        'min_hole_dim'        : 0.05,
        'effective_index'     : 1.6,
        'resonance_wavelength': 0.737
    }


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

    # specify working_path to save solved fsp file
    engine = LumericalEngine(mesh_accuracy=5, hide=hide, working_path=fsps_dir,lumerical_path=FDTDLoc,save_fsp=True)

    cavity_cells = []
    for hx, hy, a in zip(all_hx, all_hy, all_a):
        print(hx,hy,a)
        cell_size = Vec3(a,beam_w,beam_h)
        cell_box = TriStructure(Vec3(0), Vec3(beam_w, apex_half_angle, a), 
                                DielectricMaterial(n_beam, order=2, color="blue"), 
                                rot_angles=(np.pi/2, np.pi/2, 0))

        # offset the hole to respect the way we define the relevant lattice constant
        cell_hole = CylinderStructure(Vec3(-a/2,0,0), beam_h, hx/2, DielectricMaterial(1, order=1, color="red"), radius2=hy/2)
        unit_cell = UnitCell(structures = [cell_box, cell_hole], size=cell_size, engine=engine)
        cavity_cells += [unit_cell]

    # The length of the cavity beam
    # TODO: why not calculated?
    beam_length = 15e-6
    
    # shift the cavity so that the source is centered in the dielectric
    # shift = cavity_cells[int(len(cavity_cells)/2)].get_size().x/2
    shift = all_a[len(all_a)//2 - 1] / 2
    print("Center cavity a is ", 2 * shift)

    cavity = Cavity1D(
      unit_cells=cavity_cells,
      structures=[TriStructure(Vec3(0), Vec3(beam_w, apex_half_angle, beam_length), 
                  DielectricMaterial(n_beam, order=2, color="gray"), rot_angles=(np.pi/2, np.pi/2, 0))], engine=engine, center_shift = shift)

    # By setting the save path here, the cavity will save itself after each simulation to this file
    print(cavity_params)
    cavity_name = "_".join([str(n) for n in cavity_params])
    print(cavity_name)

    file_name = log_name[:-4]+"_"+cavity_name+"_"+str(iter_count)
    cavity.save(file_name+".obj")
    plot_geom(cavity,file_name+"_geom.png",hide)
    return cavity, file_name

def fitness(cavity_params):
    global source_frequency
    
    cavity, file_name = build_cavity(cavity_params)
    
    man_mesh = MeshRegion(BBox(Vec3(0), Vec3(10e-6, 0.7e-6, 0.4e-6)), 12e-9, dy=None, dz=None)

    r1 = cavity.simulate("resonance", target_freq=source_frequency, source_pulselength=60e-15, 
                        analyze_time=600e-15, mesh_regions = [man_mesh], sim_size=Vec3(1.25, 3, 5.5))

    qx = 1 / (1 / r1["qxmin"] + 1 / r1["qxmax"])
    qx1 = r1["qxmin"]
    qx2 = r1["qxmax"]

    qy = 1 / (2 / r1["qymax"])
    qz = 1 / (1 / r1["qzmin"] + 1 / r1["qzmax"])

    qscat = 1 / ((1 / qy)+(1 / qz))
    qtot = 1 / (1 / qx + 1 / qy + 1 / qz)

    vmode = r1["vmode"]
    vmode_copy = vmode
    vmode = 1e6 if vmode < 0.48 else vmode

    qtot_max = 300000
    purcell = min(qtot, qtot_max) / vmode

    F = r1["freq"]
    wavelen = (2.99e8 / F) * 1e9 # nm

    print(f"F: {r1['freq']}, Vmode: {r1['vmode']}, "
          f"Qwvg: {1/(1/r1['qxmin'] + 1/r1['qxmax'])}, "
          f"Qsc: {1/(2/r1['qymax'] + 1/r1['qzmin'] + 1/r1['qzmax'])}")
    print(purcell)

    wavelen_pen = gauss(F, target_frequency, 4e12)

    # TODO: check if / 4 is correct
    qx_pen = (gauss(qx, target_qx, 120000) + gauss(qx, target_qx, 60000) + 
              gauss(qx, target_qx, 30000) + gauss(qx, target_qx, 2000) / 4)

    qscat_max = 500000
    guidedness = min(qscat, qscat_max) / qx

    # FIGURE OF MERIT
    witness = -1 * purcell * wavelen_pen * guidedness * qx_pen

    with open(log_name, "ab") as f:
        f.write(b"\n")
        step_info = np.append(cavity_params, np.array([witness, wavelen_pen, purcell, 
                                                       r1["qxmin"], r1["qxmax"], 
                                                       qscat, qtot, vmode, vmode_copy, F]))
        np.savetxt(f, step_info.reshape(1, step_info.shape[0]), fmt='%.6f')

    r1["xyprofile"].save(file_name + "_xy.png", 
                         title=f"Q = {qtot:.0f} \nQ_scat = {qscat:.04} Qx = {qx:.0f}\nV = {vmode_copy:.3f}")
    r1["yzprofile"].save(file_name + "_yz.png",
                         title=f"Q = {qtot:.0f} Q_scat = {qscat:.04}\n Qx1 = {qx1:.0f} Qx2 = {qx2:.0f}\nV = {vmode_copy:.3f} "+r"$\lambda$"+f" = {wavelen:.1f}")

    # Plot the quasipotential
    # r2 = cavity.simulate("quasipotential", target_freq=target_frequency, freqs=(0.25e15, 0.7e15, 100000), window_pos = 0)
    # r2.show()

    # Second condition ensures that we only rerun once
    if((wavelen_pen < rerun_thresh) and (source_frequency == target_frequency)):
        # Shift source frequency to cavity resonance and rerun simulation.
        # (This should help avoid non cavities with artificially low mode volumes)
        source_frequency = F
        witness_rerun = fitness(cavity_params)
        print("rerun. Fitness when source is recentered:",witness_rerun)
        source_frequency = target_frequency
        return witness_rerun

    return witness

log_name = f"cavity_run-{nmirrsL}-{ndefs}-{nmirrsR}_111722_00.txt"

# maxDef, beam_w, hxL, hyL, hxR, hyR, aL, aR
p0 = np.array([0.14276, 0.483143,
               0.108918, 0.164633,
               0.108918, 0.164633,
               0.273535, 0.251509])

# p0 = np.array([0.1392,    0.482,
#                0.1135849, 0.1605274,
#                0.1135849, 0.1605274,
#                0.2717,    0.2502])



with open(log_name, "ab") as f:
    f.write(b"maxDef    beam_w    hx    hy    a    fitness    wavelen_pen    purcell    qxmin    qxmax    qscat    qtot    vmode    vmode_copy    freq")

fitness_result = fitness(p0)
print(fitness_result)


