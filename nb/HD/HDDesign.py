import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from IPython.display import Markdown as md
module_dir = os.path.abspath('../../pynextsw/')  # Adjust path as needed
sys.path.append(module_dir)

from pynext.system_of_units import *
import pynext.pynext_types as pn
from pynext.CylinderGeomEff import barrel_detection_efficiency
from HDF import simulate_cylinder_barrel_hits, simulate_cylinder_barrel_sectors, simulate_photons_in_cylinder

def nhits(N, A, r):
    return N * A /(4 * np.pi * r**2)
    
def adcrange(nbits):
    return 2**nbits -1

def maxpes(nbits, adctopes):
    return adcrange(nbits)/adctopes


def barrel_s2_eff(Rg,zmin, zmax, Npoints = 100, Ndir = 1000, Nsectors=500): 
    d = 2 * Rg   # mm, diameter
    L = zmax - zmin   # mm, cylinder length

    xs, total_hits, max_sector_hits = simulate_cylinder_barrel_sectors(d, L, Npoints, Ndir, Nsectors, 
                                                                       seed=123)
    return xs, total_hits/Ndir, max_sector_hits/Ndir


def ng_s2_dtp(hdtp, eblob, wblob, ts_DTP = 0.5*mus, dd=2*mm): # dd distance between anode ans SiPM plane
    area_si = hdtp.dtp.sipm.area
    dsi    =  hdtp.tpcel.d + dd
    ng_blob_DTP = int(hdtp.tpcel.el_photons(eblob))
    ns_DTP = wblob/ts_DTP
    ng_xs_DTP = int(ng_blob_DTP/ns_DTP)
    ng_xs_plane_DTP = ng_xs_DTP * hdtp.efficiency_s2 
    ng_xs_si_DTP = int(nhits(ng_xs_plane_DTP, area_si, dsi))
    
    return ns_DTP,ng_blob_DTP, ng_xs_DTP, ng_xs_si_DTP 


def power(dtp):
    nsipm = int(2 * dtp.n_sipm)
    nasic = int(nsipm/64)+1
    n_readout_tile = int(nasic/8)+ 1 
    p_readout_tile = 14 * watt
    tpw = p_readout_tile * n_readout_tile
    return nsipm, nasic, n_readout_tile, p_readout_tile, tpw 
    
 
def cost_dtp(hdtp, sipm_cost_mm2, asic_cost_channel):
    cost_si =hdtp.dtp.area/mm2 * sipm_cost_mm2 /1e+6
    cost_asic = hdtp.dtp.n_sipm * asic_cost_channel/1e+6
    return cost_si, cost_asic
    
    
def detector_parameters_md(Ps, Ts, hdxe):
    f =f"""## Detector parameters at P={Ps} bar, T= {Ts}C\n
- R = {hdxe.radius/m} mm;  L = {hdxe.length/m} m 
- Area barrel ={hdxe.area_barrel/m2:.2f} m2
- Area end-cap ={hdxe.area_endcap/m2:.2f} m2
- Volume= {hdxe.volume/m3:.2f} m3 
- $\\rho$ ={hdxe.density/(kg/m3):.2f} kg/m3
- mass = {hdxe.mass/kg:.1f} kg"""
    return md(f)
    

def md_el(Ps, Ts, hdel):
    f0 =f"## Detector EL parameters at P={Ps} bar, T= {Ts}C\n"
    kv_cm_b = kilovolt / (cm * bar)
    kv_cm   = kilovolt / cm
    cm_bar = 1 / (cm * bar)
    s= """
- E/P = %7.2f kV * cm$^{-1} \\times$ bar$^{-1}$
- Drift voltage = %7.2f kV $\\times$ cm$^{-1}$
- EL grid gap = %7.2f mm
- Drift lenght =%7.2f m
- Grid voltage = %7.2f kV
- Cathode voltage = %7.2f kV
- Yield =  %7.2e photons/e

"""%(hdel.e_over_p / kv_cm_b,
     hdel.drift_voltage / kv_cm,
     hdel.grid_distance / mm,
     hdel.drift_length / m,
     hdel.grid_voltage / kilovolt,
     hdel.cathode_voltage /kilovolt,
     hdel.optical_gain)
    

    return md(f0+s)

def kr_qbb_s1s2(tpcel):
    s=f"""## S1 and S2 photons for Kr and $Q_{{\\beta\\beta}}$ \n
- Primary scintillation Krypton              = {tpcel.scintillation_photons(41.5 * keV):.2e}
- EL photons Krypton                         = {tpcel.el_photons(41.5 * keV):.2e}
- Primary scintillation $Q_{{\\beta\\beta}}$ = {tpcel.scintillation_photons(2458 * keV):.2e}
- EL photons $Q_{{\\beta\\beta}}$            = {tpcel.el_photons(2458 * keV):.2e}
"""
    return md(s)
    

def fiber_wls(fwls):
    f0 =f"## WLS fibers\n"
    s= f"""
    diameter ={fwls.d/mm} mm, Q = {fwls.qfib}
    ncore = {fwls.ncore}, nclad1 ={fwls.nclad1}, nclad2 ={fwls.nclad2}
    Absoprtion prob at 450 nm     = {fwls.blue_absorption(450*nm)}
    Trapping efficieny            = {fwls.trapping_efficiency}
    Fiber coated with WLS         = {fwls.wls.name}
    WLS QE                        = {fwls.wls.quantum_efficiency}

    """
    return md(f0+s)
    

def fiber_detector_s1(fdHD):
    ren = fdHD.total_efficiency/ fdHD.n_fibers
    
    s =f"""## Light efficiency collection S1\n
- Fibers efficiency = {100 * fdHD.fiber_efficiency:.2f} %
    - Transport = {100 * fdHD.transport:.2f} %
    - Attenuation = {100* fdHD.attenuation:.2f} %
- SiPM PDE          = {100* fdHD.sipm.PDE:.2f} %
    
- Total detection efficiency = {100* fdHD.total_efficiency:.1f} %
    
- Sampling S1       = {fdHD.sampling/ns:.1f} ns
- Operating temperature ({fdHD.operating_temperatureC:.1f} C)
    
- Primary scintillation Krypton          = {fdHD.tpcel.scintillation_photons(41.5 * keV):.2e}
- Primary scintillation Krypton detected = {fdHD.tpcel.scintillation_photons(41.5*keV)*fdHD.total_efficiency:.1f} 
- Number of DCR photons in the detector = {fdHD.dcr_sipm_per_time(fdHD.tempC, fdHD.sampling):.2f}
    """
    return md(s)
    

def bfd(fdHD):
    f0 =f"## BFD\n"
    
    ren = fdHD.total_efficiency/ fdHD.n_fibers
    
    s =f"""
    Total detection efficiency = {100*fdHD.total_efficiency:.2f} %
    SiPM PDE          = {fdHD.sipm.PDE:.1f}
    fiber size        = {fdHD.fwls.diameter/mm:.1f} mm
    Sampling S1       = {fdHD.sampling/ns:.1f} ns
    number of fibers  = {fdHD.n_fibers:d}

    Primary scintillation Krypton          = {fdHD.tpcel.scintillation_photons(41.5 * keV):.2e}
    EL photons Krypton                     = {fdHD.tpcel.el_photons(41.5 * keV):.2e}
    Primary scintillation Qbb              = {fdHD.tpcel.scintillation_photons(2458 * keV):.2e}
    EL photons Qbb                         = {fdHD.tpcel.el_photons(2458 * keV):.2e}

    Primary scintillation Krypton detected = {fdHD.tpcel.scintillation_photons(41.5 * keV)* fdHD.total_efficiency:.2e}
    EL photons Krypton detected            = {fdHD.tpcel.el_photons(41.5 * keV)* fdHD.total_efficiency:.2e}
    Primary scintillation Qbb detected     = {fdHD.tpcel.scintillation_photons(2458 * keV)* fdHD.total_efficiency:.2e}
    EL photons Qbb detected                = {fdHD.tpcel.el_photons(2458 * keV)* fdHD.total_efficiency:.2e}

    Number of DCR photons in the detector for:
     --operating temperature ({fdHD.operating_temperatureC:.2f} C)
     --sampling time of {fdHD.sampling/ns:.2f} 
     -- nDCR = {fdHD.dcr_sipm_per_time(fdHD.tempC, fdHD.sampling):.2f} pes
    """
    return md(f0+s)


def dr_fiber_detector(eblob,wblob,ts, nxs, ngblob, ngblobxs, ngx, ngx2, rf):
     
    ff = f"### Dynamic range for fiber detector\n"
    
    s=f"""
- Assume $E_b \\sim$ {eblob} keV
- Minimum diffusion ~{wblob/mus} mus.
- Sampling: {ts} ns.
- Thus: number of samples = {int(nxs)}
- N$_\\gamma$ (blob) = {ngblob}
- N$_\\gamma$ (sample) = {ngblobxs}
- N$_\\gamma$ (max/sector) = {ngx}
- Cutting at Rf = {rf:.2f}
- N$_\\gamma$ (Rf) = {ngx2}
- PES for 13 bit ADC = {int(maxpes(13, 15))}
    """
    return md(ff+s)


def dr_dtp(eblob,wblob,ts, nxs, ngblob, ngblobxs, ngx):
     
    ff = f"## Dynamic range for DTP detector\n"
    
    s=f"""
- Assume $E_b \\sim$ {eblob/keV} keV
- Minimum diffusion ~{wblob/mus} mus.
- Sampling: {ts/mus} mus.
- Thus: number of samples = {int(nxs)}
- N$_\\gamma$ (blob) = {ngblob}
- N$_\\gamma$ (sample) = {ngblobxs}
- N$_\\gamma$ (SiPM) = {ngx}
- PES for 13 bit ADC = {int(maxpes(13, 11))}
    """
    return md(ff+s)


def hdtp_dc(hdtp, Ps, Ts):
    elsqd = hdtp.tpcel.el_photons(2458 * keV)* hdtp.efficiency_s2 * hdtp.dtp.ff
    u = hertz / mm2
    s=f"""## DTP\n
- Pressure = {Ps} bar
- Temperature = {Ts} C
- Filling factor = {hdtp.dtp.ff}
- area = {hdtp.dtp.area/mm2:.2g} mm$^2$, 
- number of SiPMs = {hdtp.dtp.n_sipm:.3g}
- Light efficiency (s1) = {100*hdtp.efficiency_s1:.2f} %
- Light efficiency (s2) = {100*hdtp.efficiency_s2:.2f} %
- SiPM (efective) PDE   = {hdtp.dtp.sipm.pde:.2f}
- SiPM size         = {hdtp.dtp.sipm.xsize/mm:.2f} mm
- Sampling S1       = {hdtp.sampling/ns:.2f} ns
- Primary scintillation Qbb keV detected = {hdtp.tpcel.scintillation_photons(2480 * keV)* hdtp.efficiency_s1:.2e}
- EL photons $Q_{{\\beta\\beta}}$ detected                = {elsqd:.2e}
- nDCR = {2* hdtp.dtp.dcr_sipm_per_time(hdtp.operating_temperatureC, hdtp.sampling):.2e}
    """
    return md(s)


def dc_cooling(nsipm, nasic, n_readout_tile, p_readout_tile, tpw):
    s=f""" ## Power dissipation in DTP\n
- number of SiPMs in fDTP ={nsipm} 
- number of ASICS in fDTP ={nasic} 
- number of readout tiles = {n_readout_tile} 
- power per readout tiles = {p_readout_tile/watt} watt
- Total Power to dissipate = {tpw/watt:.2e} watt
    """
    return md(s)

def cost_dtp_md(cost_si, cost_asic):
    f = f"""## Cost of DTP\n
- cost of SiPMs = {cost_si:.2f} M€
- cost of electronics = {cost_asic:.2f} M€
- Total =  {cost_asic + cost_si:.2f} M€
    """
    return md(f)

def display_png(fn, title=" ", figsize=(8, 6)):
    plt.figure(figsize=figsize)
    img = mpimg.imread(fn)
    plt.imshow(img)
    plt.title(title)
    plt.axis('off')  # Hide axes
    plt.show()
    

def display_BFD(figsize=(8, 6)):
    plt.figure(figsize=figsize)
    img = mpimg.imread('BFD.png')
    plt.imshow(img)
    #plt.title("BFD")
    plt.axis('off')  # Hide axes
    plt.show()
    

def display_HD(figsize=(8, 6)):
    plt.figure(figsize=figsize)
    img = mpimg.imread('HDB.png')
    plt.imshow(img)
    plt.axis('off')  # Hide axes
    plt.show()

# Plot data and fitted line
def plot_dc(figsize=(8, 6)):
    xx = np.arange(-100, 20, 0.1)
    yy = pn.logdc_vs_t(xx)
    fig, axs = plt.subplots(1, 1, figsize=figsize)
    axs.scatter(pn.tC, pn.ldR, label='Data')
    axs.plot(pn.tC, pn.y_fit, 'r-', label='Fitted line')
    axs.plot(xx, yy, 'g-', label='log(dc)')
    axs.legend()
    axs.grid()
    axs.set_title("Dependence of DC with Temperature (Hamamatsu S13360 series)")
    axs.set_xlabel("T (C)")
    axs.set_ylabel("log(dC) in kHz")
    fig.tight_layout()


def display_dtp(dtp):
        u = hertz / mm2
        s =f"""
        - DTP: Fill factor = {dtp.ff}
        - area = {dtp.area/mm2:.2g} mm2, number of SiPMs = {dtp.n_sipm:.3g}
        - SiPMs: 
            {dtp.sipm}
        - DCR  (tC = 20C) = {dtp.dcr_sipm(20)/u:.2g} Hz;
        - DCR  time window (tC = 20C, t=1 mus) = {dtp.dcr_sipm_per_time(20, 1*mus):.2g} counts;
        
        """
        
        return md(s)




def dcr_bfd_kr(hdel, hdxe, fwls, s1mm, figsize=(6,4)):
    DCR = []
    tt = np.arange(-20, 20, 0.5)
    
    for t in tt:
        fdHD = pn.FiberDetector(hdel, hdxe, fwls, s1mm, eff_t=0.85, sampling= 100 * ns, tempC=t)
        DCR.append(fdHD.dcr_sipm_per_time(fdHD.tempC, fdHD.sampling))

    nkr = fdHD.tpcel.scintillation_photons(41.5 * keV)* fdHD.total_efficiency
    yy = np.array(DCR)
    fig, axs = plt.subplots(1, 1, figsize=figsize)
    axs.plot(tt, yy)
    plt.axhline(y=nkr, color='r', linestyle='--', linewidth=2)
    axs.grid()
    axs.set_title("Dependence of DC with Temperature (BFD, Kr)")
    axs.set_xlabel("T (C)")
    axs.set_ylabel("dcr")


def plot_s2_dr(xs, total_hits, max_sector_hits, ngblobxs, Rg, Rc, mbeff):
    xx = np.arange(-Rg+1, Rg)
    fig, axs = plt.subplots(1, 2, figsize=(12,4))
    axs[0].plot(xs, max_sector_hits)
    axs[0].grid()
    axs[0].set_title("Max eff = (R)")
    axs[0].set_xlabel(" R (mm)")
    axs[0].set_ylabel("eff")
    axs[1].plot(xx, mbeff(xx)* ngblobxs)
    plt.axvline(x=-Rc, color='r', linestyle='--', linewidth=2)
    plt.axvline(x=Rc, color='r', linestyle='--', linewidth=2)
    plt.axhline(y=maxpes(13, 15), color='g', linestyle='--', linewidth=2)
    plt.text(-600, 600, "13 bits ADC", fontsize=16, color='blue')
    axs[1].grid()
    axs[1].set_title("Nphe = (R)")
    axs[1].set_xlabel(" R (mm)")
    axs[1].set_ylabel("eff")
    fig.tight_layout()


def plot_dtp_dc(tt, dc, label, figsize=(6,4)):
    
    fig, axs = plt.subplots(1, 1, figsize=figsize)
    axs.plot(tt, dc, 'r', label=label)
    axs.legend()
    axs.grid()
    axs.set_title("Total DC/sample (1 mus) in DTP")
    axs.set_xlabel("T(C)")
    axs.set_ylabel("DC (pes)")
    fig.tight_layout()



def dcr_dtp(vhdel, vhdxe, vdtp, eff_grid=0.9, sampling=0.5*mus, figsize=(6,4)):
    DCR = []
    tt = np.arange(-40, -10, 0.5)
    
    for t in tt:
        vhdtp = pn.DTPDetector(vhdel, vhdxe, vdtp, eff_grid, sampling, tempC=t)
        ndcr = 2* vhdtp.dtp.dcr_sipm_per_time(vhdtp.operating_temperatureC, vhdtp.sampling)
        DCR.append(ndcr)

    yy = np.array(DCR)
    
    nkr = vhdtp.tpcel.scintillation_photons(41.5 * keV)* vhdtp.efficiency_s1
    n511 = vhdtp.tpcel.scintillation_photons(511 * keV)* vhdtp.efficiency_s1
    nqbb = vhdtp.tpcel.scintillation_photons(2458 * keV)* vhdtp.efficiency_s1
    #print(nkr, n511, nqbb)
   
    fig, axs = plt.subplots(1, 1, figsize=figsize)
    axs.plot(tt, yy)
    plt.axhline(y=n511, color='r', linestyle='--', linewidth=2)
    plt.text(-40, 250, "511 keV", fontsize=14, color='blue')
    plt.axhline(y=nqbb, color='g', linestyle='--', linewidth=2)
    plt.text(-40, 850, "Qbb", fontsize=14, color='blue')
    axs.grid()
    axs.set_title("Dependence of DC with Temperature (DTP, VUV)")
    axs.set_xlabel("T (C)")
    axs.set_ylabel("dcr")

    
def set_fonts(ax, fontsize=20):
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)

def plotxy(xx, yy, xl, yl, tit, color='r', label="", figsize=(6,4)):
    fig, axs = plt.subplots(1, 1, figsize=figsize)
    axs.plot(xx, yy, color, label=label)
    axs.legend()
    axs.grid()
    axs.set_title(tit)
    axs.set_xlabel(xl)
    axs.set_ylabel(yl)
    fig.tight_layout()






