#-*- coding: utf-8 -*-

### DISCRETIZATION and VISUALIZATION jobs
# Depending on the job_types (OPT, FREQ, TD, etc.) different pictures will be generated

import os
import math
import numpy as np
from PIL import Image, ImageOps 

from quchemreport.processing import calc_orb, TD2UVvis
from quchemreport.visualization import visu_mayavi, visu_plots, visu_txt


FWHM = 3000 # 3000 cm-1 for full width at half maximum of gaussian band. 


def autocrop(imgfile):
    image=Image.open(imgfile)
    # The automatic crop only works with RGB mode
    if image.mode != 'RGB':
        RGB_im=image.convert("RGB")
    else :
        RGB_im=image
    # Invert image to get a black border
    invert_im = ImageOps.invert(RGB_im)
    # get box search for outer blanck pixels in pictures
    imageBox = invert_im.getbbox()
    cropped=image.crop(imageBox)
    cropped.save(imgfile, fomat='PNG', dpi=(300,300))


def jobsSI(jf, job_types, nres_noES, charges, charge_ref, discret_proc, mo_viz_done, data_for_discretization, nproc, restart):
# Parameters list of logfiles and JSON data, list of job-types and key data obtained by the conformity tests.

    # SETUP for ORBKIT
    ## Get grid parameters and initialize grid
    ## Oversizing the grid aroung the molecule in Bohr radii
    ## ORBKIT uses 5 by default, tune this as required
    extand = 5
    ## Spacing of cubic voxels in the grid 
    step = 0.1
    # MO list initialization
    MO_list = []
    # electronic transitions   
    Elec_T = []  


    for i, jt in enumerate(job_types):	  
        if job_types[i] == ['OPT'] or job_types[i] == ['FREQ', 'OPT'] or job_types[i] == ['FREQ', 'OPT', 'TD'] : 
            # Test if Ground state by geometry, if atomic positions present: generate picture
            if jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"] == nres_noES[0] \
                              and jf[i]["results"]["geometry"]["elements_3D_coords_converged"] != 'N/A': 
                visu_mayavi.topo(jf[i],file_name="img")
                print("Picture generated for the topology of the ground state")
            # Test if Ground state by geometry and by charge and if MO coeffs presents 
            # So generate HOMO and LUMO pictures for the ground state at the reference oxydation state
            if jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"] == nres_noES[0] \
                               and int(jf[i]["molecule"]["charge"]) == charge_ref \
                               and (discret_proc is True) and (mo_viz_done is False):
                # SETUP for the HOMO and LUMO (HOMO+1) index list. Orbkit starts MO_list counter from zero
                # Test if calculation is unrestricted (alpha and beta spin electrons)
                HO_ind = jf[i]["results"]["wavefunction"]["homo_indexes"]
                if len(HO_ind) == 2:
                    # Unrestricted calculation: treat the alpha orbitals first
                    print("Unrestricted calculation detected")
                    MO_list_alpha = [HO_ind[0]-1, HO_ind[0], HO_ind[0]+1, HO_ind[0]+2]
                    MO_labels_alpha = ['homo-1_alpha','homo_alpha', 'lumo_alpha', 'lumo+1_alpha'] 
                    if (restart == 1) and (os.path.isfile("img-MO-homo-1_alpha.png")) and (os.path.isfile("img-MO-homo_alpha.png")) \
                                      and (os.path.isfile("img-MO-lumo_alpha.png")) and (os.path.isfile("img-MO-lumo+1_alpha.png")):
                        print("Alpha Molecular orbitals pictures already done!")
                    else:                                   
                        ## Calculations of MO
                        out, X, Y, Z = calc_orb.MO(data_for_discretization, MO_list_alpha, spin="alpha", grid_step=step, nproc=nproc)
                        ## Visulation of the MO
                        visu_mayavi.viz_MO(out, X, Y, Z, data_for_discretization, file_name="img", labels=MO_labels_alpha)
                        mo_viz_done = True
                    # Unrestricted calculation: now treat the beta orbitals
                    MO_list_beta = [HO_ind[1]-1, HO_ind[1], HO_ind[1]+1, HO_ind[1]+2] 
                    MO_labels_beta = ['homo-1_beta','homo_beta', 'lumo_beta', 'lumo+1_beta'] 
                    if (restart == 1) and (os.path.isfile("img-MO-homo-1_beta.png")) and (os.path.isfile("img-MO-homo_beta.png")) \
                                      and (os.path.isfile("img-MO-lumo_beta.png")) and (os.path.isfile("img-MO-lumo+1_beta.png")):
                        print("Beta Molecular orbitals pictures already done!")
                    else:                                   
                        ## Calculations of MO
                        out, X, Y, Z = calc_orb.MO(data_for_discretization, MO_list_beta, spin="beta", grid_step=step, nproc=nproc)
                        ## Visulation of the MO
                        visu_mayavi.viz_MO(out, X, Y, Z, data_for_discretization, file_name="img", labels=MO_labels_beta)
                        mo_viz_done = True

                else:
                    MO_list = ['homo', 'lumo', 'lumo+1']
                    MO_labels = MO_list
                    if (restart == 1) and (os.path.isfile("img-MO-lumo.png")) and (os.path.isfile("img-MO-homo.png")) :
                        print("Molecular orbitals pictures already done!")
                    else:                                   
                        ## Calculations of MO
                        out, X, Y, Z = calc_orb.MO(data_for_discretization, MO_list, spin="none", grid_step=step, nproc=nproc)
                        ## Visulation of the MO
                        visu_mayavi.viz_MO(out, X, Y, Z, data_for_discretization, file_name="img", labels=MO_labels)
                        mo_viz_done = True

              

    
        #if  job_types[i] == ['FREQ']  or job_types[i] == ['FREQ', 'OPT'] :
        # For now, no spectrum of IR. 
        
        if job_types[i] == ['TD'] or job_types[i] == ['FREQ', 'OPT', 'TD'] :
            if jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"] == nres_noES[0] \
                and int(jf[i]["molecule"]["charge"]) == charge_ref :                
                print("Electronic transitions detected. Calculating the UV absorption spectrum.")                     
                if jf[i]["results"]["excited_states"]["et_oscs"].count(0.0) == len(jf[i]["results"]["excited_states"]["et_oscs"]):
                    print("The oscillators strengths are always 0.0. The absorption spectrum won't be plotted.")
                else :
                    # Spectrum width should min cm-1 - 9000 (min 5000), max cm-1 + 9000 (max 100000 cm-1)  
                    td_start = round(min(jf[i]["results"]["excited_states"]["et_energies"]) - 9000, -3)
                    td_end = round(max(jf[i]["results"]["excited_states"]["et_energies"]) + 9000, -3)
                    if td_start < 5000:
                        td_start = 5000
                    if td_end > 100000 :
                        td_end = 100000
                    numpts = int(np.floor((td_end - td_start)/20))
                    et_energies = jf[i]["results"]["excited_states"]["et_energies"]
                    et_oscs = jf[i]["results"]["excited_states"]["et_oscs"]
                    if ((len(et_energies)) > 0 ) and  ((len(et_oscs)) > 0 ):
                        # Calculating the Gaussian broadening based on wavenumbers
                        heights = [[x*2.174e8/FWHM for x in et_oscs]]
                        xvalues, spectrum = TD2UVvis.GaussianSpectrum(td_start,td_end,numpts,et_energies,heights,FWHM)
                        # If Rotational strength in calculation treat Circular Dichroism       
                        try:
                            et_rotats=jf[i]["results"]["excited_states"]["et_rot"]
                        except KeyError:
                            et_rotats = []
                        if ((len(et_rotats)) == 0 ) and (et_rotats.count(0.0) == len(et_rotats)):
                            print("The rotational strengths are always 0.0. The circular dichroism spectrum won't be plotted.")
                            CDspectrum = []
                        else : 
                            heights = TD2UVvis.CDheights(et_energies, et_rotats, FWHM)
                            xvalues, CDspectrum = TD2UVvis.GaussianSpectrum(td_start,td_end,numpts,et_energies,[heights], FWHM)

                        # Output calculated spectra in text files and figures
                        if not ((len(spectrum)) == 0 ):
                            if not ((len(CDspectrum)) == 0 ):
                                # Write both Absorption and Circular dichroism spectra in text file to compare with experimental data
                                visu_txt.textfile_CD(xvalues, spectrum, CDspectrum)                                 
                                # IF CD plot CD spectrum 
                                visu_plots.absoCD(et_energies, et_rotats, xvalues, CDspectrum)
                                print("Circular dichroism spectrum done.")
                            else :
                                # If no CD write only Absorption spectrum in text file to compare with experimental data
                                visu_txt.textfile_UV(xvalues, spectrum)                                                                                                   
                            # Plotting UV spectrum 
                            visu_plots.absoUV(et_energies, et_oscs, xvalues, spectrum)
                            print("Absorption spectrum done.")

                # Calculating Electronic density differences discretization and transition properties (Tozer lambda and charge tranfer)
                # Even if all transition are silent. Because of triplet states
                if discret_proc is True :
                    print("Calculating the electronic density differences.")
                    # Set the transition list to process = T_list. By default all transtions are taken into account.
                    T_list = []
                    et_transitions = jf[i]["results"]["excited_states"]["et_transitions"]
                    et_sym = jf[i]["results"]["excited_states"]["et_sym"]
                    jf[i]["results"]["excited_states"]["Tozer_lambda"] = ["N/A"] * len(et_energies)
                    jf[i]["results"]["excited_states"]["d_ct"] = ["N/A"]* len(et_energies) 
                    jf[i]["results"]["excited_states"]["q_ct"] = ["N/A"] * len(et_energies)
                    jf[i]["results"]["excited_states"]["mu_ct"] = ["N/A"] * len(et_energies)
                    jf[i]["results"]["excited_states"]["e-_barycenter"] = ["N/A"] * len(et_energies)
                    jf[i]["results"]["excited_states"]["hole_barycenter"] = ["N/A"] * len(et_energies)
                    ## Discretization of all MO used in the transitions
                    out, X, Y, Z = calc_orb.TD(data_for_discretization, et_transitions, grid_step=step, nproc=nproc)
                    for k, transitions in enumerate(et_transitions):
                        if restart == 0 :
                            print("EDD visualization in progress for the transition(s):", k)
                            visu_mayavi.viz_EDD([out[k][0]], X, Y, Z, data_for_discretization, et_sym[k], 
                                                                 file_name="img", labels=[k+1]) 
                            if (len(et_rotats) > 0):
                                visu_mayavi.viz_BARY([out[k][2][3:]], data_for_discretization, et_sym[k],
                                                          file_name="img", labels=[k+1]) 
                        ## Returns the calculated values of the tozer_lambda, d_CT, Q_CT, Mu_CT and e- barycenter and hole barycenter to the json   
                        jf[i]["results"]["excited_states"]["Tozer_lambda"][k] = out[k][1]
                        jf[i]["results"]["excited_states"]["d_ct"][k] = out[k][2][0]
                        jf[i]["results"]["excited_states"]["q_ct"][k] = out[k][2][1]
                        jf[i]["results"]["excited_states"]["mu_ct"][k] = out[k][2][2]
                        jf[i]["results"]["excited_states"]["e-_barycenter"][k] = out[k][2][3].tolist()
                        jf[i]["results"]["excited_states"]["hole_barycenter"][k] = out[k][2][4].tolist()
                    print("EDD done")
    
    
        if  jf[i]["comp_details"]["general"]["job_type"] == ['OPT_ES'] :
            print("Optimization of an excitated state detected. Calculating the UV emission spectrum")
            # TODO : get optimized excited state number 
            # Nothing in the log file in Gaussian allow us to clearly identify the Excited state number that is optimized
            # Workaround based on log files names !!! Expected log names OPT_ESX_step_Y.inp or OPT_ETX_step_Y.inp with X the excited state number            
            emi_state = [int(s) for s in jf[i]["metadata"]["log_file"] if s.isdigit()][0] 
            emi_index = emi_state-1
            if (emi_index) < 0:
                print("Incoherent excited state optimization detected. Problem with root number")
            else:
                emi_transition = jf[i]["results"]["excited_states"]["et_transitions"][emi_index]
                emi_sym = jf[i]["results"]["excited_states"]["et_sym"][emi_index]
                emi_energy = jf[i]["results"]["excited_states"]["et_energies"][emi_index]
                emi_osc = jf[i]["results"]["excited_states"]["et_oscs"][emi_index]
                try:
                    emi_rotat=jf[i]["results"]["excited_states"]["et_rot"][emi_index]
                except KeyError:
                    emi_rotat = 0.0
                jf[i]["results"]["excited_states"]["Tozer_lambda"] = ["N/A"] 
                jf[i]["results"]["excited_states"]["d_ct"] = ["N/A"] 
                jf[i]["results"]["excited_states"]["q_ct"] = ["N/A"] 
                jf[i]["results"]["excited_states"]["mu_ct"] = ["N/A"] 
                jf[i]["results"]["excited_states"]["e-_barycenter"] = ["N/A"] 
                jf[i]["results"]["excited_states"]["hole_barycenter"] = ["N/A"] 
                # Calculating Electronic density difference discretization and transition property (Tozer lambda and charge tranfer)
                if discret_proc is True :
                    print("Calculating the emission electronic density difference.")
                    ## Discretization of all MO used in the transition
                    out, X, Y, Z = calc_orb.TD(data_for_discretization, [emi_transition], grid_step=step, nproc=nproc)
                    print("EDD visualization in progress for the transition :", emi_state)
                    visu_mayavi.viz_EDD([out[0][0]], X, Y, Z, data_for_discretization, emi_sym, 
                                                             file_name="img-emi", labels=[emi_state]) 
                    if (emi_rotat != 0.0):
                        visu_mayavi.viz_BARY([out[0][2][3:]], data_for_discretization, emi_sym,
                                                      file_name="img-emi", labels=[emi_state]) 
                    ## Returns the calculated values of the tozer_lambda, d_CT, Q_CT, Mu_CT and e- barycenter and hole barycenter to the json   
                    jf[i]["results"]["excited_states"]["Tozer_lambda"][0] = out[0][1]
                    jf[i]["results"]["excited_states"]["d_ct"][0] = out[0][2][0]
                    jf[i]["results"]["excited_states"]["q_ct"][0] = out[0][2][1]
                    jf[i]["results"]["excited_states"]["mu_ct"][0] = out[0][2][2]
                    jf[i]["results"]["excited_states"]["e-_barycenter"][0] = out[0][2][3].tolist()
                    jf[i]["results"]["excited_states"]["hole_barycenter"][0] = out[0][2][4].tolist()
                    print("EDD done")
    
                if emi_osc == 0.0 :
                    print("The oscillator strength is 0.0. The emission spectra won't be plotted.")
                else :
                    # Spectrum width should min cm-1 - 9000 (min 5000), max cm-1 + 9000 (max 100000 cm-1)  
                    td_start = round(int(emi_energy) - 9000, -3)
                    td_end = round(int(emi_energy) + 9000, -3)
                    if td_start < 5000:
                        td_start = 5000
                    if td_end > 150000 :
                        td_end = 150000
                    numpts = int((td_end - td_start)/20)
                    # Calculating the Gaussian broadening based on wavenumbers
                    heights = [[x*2.174e8/FWHM for x in et_oscs]]
                    xvalues, spectrum = TD2UVvis.GaussianSpectrum(td_start,td_end,numpts,[emi_energy],heights,FWHM)
                    # If Rotational strength in calculation treat Circular Dichroism       
                    if emi_rotat  == 0.0 :
                        print("The rotational strength is 0.0. The emission circular dichroism spectra won't be plotted.")
                    else : 
                        heights = TD2UVvis.CDheights([emi_energy], [emi_rotat], FWHM)
                        xvalues, CDspectrum = TD2UVvis.GaussianSpectrum(td_start,td_end,numpts,[emi_energy],[heights], FWHM)

                    # Output calculated spectra in text files and figures
                    if not ((len(spectrum)) == 0 ):
                        if not ((len(CDspectrum)) == 0 ):
                            # Write both Absorption and Circular dichroism spectra in text file to compare with experimental data
                            visu_txt.textfile_emiCD(xvalues, spectrum, CDspectrum)                                 
                            # IF CD plot CD spectrum 
                            visu_plots.emiCD([emi_energy], [emi_rotat], xvalues, CDspectrum)
                            print("Circular dichroism spectrum done.")
                        else :
                            # If no CD write only Absorption spectrum in text file to compare with experimental data
                            visu_txt.textfile_emiUV(xvalues, spectrum)                                                                                                   
                        # Plotting UV spectrum 
                        visu_plots.emiUV([emi_energy], [emi_osc], xvalues, spectrum)
                        print("Absorption spectrum done.")

            
              
    # Autocrop generated images
    for file in os.listdir("."):
        if file.endswith(".png"):
            print("Autocrop image:", file)
            autocrop(file)
            
