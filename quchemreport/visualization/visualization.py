#-*- coding: utf-8 -*-

### DISCRETIZATION and VISUALIZATION jobs
# Depending on the job_types (OPT, FREQ, TD, etc.) different pictures will be generated

import os
import sys
import math
import numpy as np
import resource
from PIL import Image, ImageOps 

from quchemreport.processing import calc_orb, TD2UVvis
from quchemreport.visualization import visu_mayavi, visu_plots, visu_txt
from quchemreport.utils.parameters import FWHM, obk_step, obk_extand

extand = obk_extand
step = obk_step

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

def jobs(args, jf, data):
    
    # 3 report types are considered. 
    # Full report have an original layout. Pretty but inappropiate for pdf2word conversion.
    # Full report consider all data and needs discretization process: Full log verbosity. Orbkit calculations. 
    # SI mode will not calculate the Fukui values and pictures, Analysis populations, Excited states dipoles. 
    # text mode reports are the most simple ones and discard all discretization pictures and calculations
    # Based on the SI mode, we remove the FMO, EDD pictures, CDFT and Charge transfer data.
    report_type = args['mode']
    
    # Parameters list of logfiles and JSON data, list of job-types and key data obtained by the conformity tests.
    job_types = data['job_types']
    nres_noES = data['nres_noES']
    charges = data['charges']
    charge_ref = data['charge_ref']
    discret_proc = data['discret_proc']
    mo_viz_done = data['mo_viz_done']
    data_for_discretization = data['data_for_discretization']
    nproc = args['nproc']
    restart = args['restart']
    doMEP = args['MEP']
    maxMem = args['mem'] * 1e+9
    verbose = args['verbose']
    # MO list initialization
    MO_list = []
    # electronic transitions   
    Elec_T = []  
    #resource.setrlimit(resource.RLIMIT_AS, (maxMem*nproc, maxMem*nproc))
    
    if (discret_proc is False) and (report_type != "text"):
        print("If discretization is not possible, only a text report can be generated. Report mode changed to text")
        report_type = 'text'
        
    for i, jt in enumerate(job_types):	  
        if job_types[i] == ['OPT'] or job_types[i] == ['FREQ', 'OPT'] or job_types[i] == ['FREQ', 'OPT', 'TD'] : 
            # Test if Ground state by geometry, if atomic positions present: generate picture
            if jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"] == nres_noES[0] \
                              and jf[i]["results"]["geometry"]["elements_3D_coords_converged"] != 'N/A': 
                visu_mayavi.topo(jf[i],file_name="img")
                print("Picture generated for the topology of the ground state")
            # Test if Ground state by geometry and by charge and if MO coeffs presents 
            if report_type != 'text':
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
                        #MO_list = ['homo-7', 'homo-6', 'homo-5', 'homo-4', 'homo-3' ,'homo-2','homo-1','homo', 'lumo', 'lumo+1', 'lumo+2']
                        MO_list = ['homo-1', 'homo', 'lumo', 'lumo+1']
                        #MO_list = [HO_ind[0]-7, HO_ind[0]-6, HO_ind[0]-5, HO_ind[0]-4, HO_ind[0]-3, HO_ind[0]-2, HO_ind[0]-1, HO_ind[0], HO_ind[0]+1, HO_ind[0]+2, HO_ind[0]+3, HO_ind[0]+4]
                        MO_labels = MO_list
                        if (restart == 1) and (os.path.isfile("img-MO-homo-1.png")) and (os.path.isfile("img-MO-homo.png")) \
                                          and (os.path.isfile("img-MO-lumo.png")) and (os.path.isfile("img-MO-lumo+1.png")):
                            print("Molecular orbitals pictures already done!")
                        else:                                   
                            ## Calculations of MO
                            out, X, Y, Z = calc_orb.MO(data_for_discretization, MO_list, spin="none", grid_step=step, nproc=nproc)
                            ## Visulation of the MO
                            visu_mayavi.viz_MO(out, X, Y, Z, data_for_discretization, file_name="img", labels=MO_labels)
                            mo_viz_done = True

                    # Since the electrostatic potential is a very long process. Check if the png file exist. Therefore 
                    if (restart == 1) and (os.path.isfile("img-MEP.png")) :
                        print("Electrostatic potential map already done!")
                    elif doMEP: 
                        print("Starting calculations of Molecular Electrostatic Potential Map...")
                        ## Calculations of Molecular Electrostatic Potential Map
                        try:
                            dens, pot, X, Y, Z = calc_orb.Potential(data_for_discretization, grid_step=step, nproc=nproc)
                            ## Visulation of the MO
                            visu_mayavi.viz_Potential(dens, pot, X, Y, Z, data_for_discretization, file_name="img")
                        except MemoryError:
                            sys.stderr.write('\n\nERROR: Memory Exception during calculations of Molecular Electrostatic Potential\n')
    
        #if  job_types[i] == ['FREQ']  or job_types[i] == ['FREQ', 'OPT'] :
        # For now, no spectrum of IR. 
        
        if job_types[i] == ['TD'] or job_types[i] == ['FREQ', 'OPT', 'TD'] :
            if jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"] == nres_noES[0] \
                and int(jf[i]["molecule"]["charge"]) == charge_ref :                
                print("Electronic transitions detected. Calculating the UV absorption spectrum.")                     
                # For ORCA
                if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in jf[i]["results"]["excited_states"]["et_energies"] and \
                'ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS' in jf[i]["results"]["excited_states"]["et_energies"]:
                    jf[i]["results"]["excited_states"]["et_energies_vel"] = jf[i]["results"]["excited_states"]["et_energies"]['ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS'][0]
                    jf[i]["results"]["excited_states"]["et_energies"] = jf[i]["results"]["excited_states"]["et_energies"]['ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS'][0]

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
                if report_type != 'text':
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
                        # Dipolar moment of Ground state norm in x, y, z. 
                        gs_eldip = (np.array(jf[i]["results"]["wavefunction"]["moments"][1]), np.array(jf[i]["results"]["wavefunction"]["moments"][0])) 
                        #print(gs_eldip)
                        
                                   
                        ## Discretization of all MO used in the transitions
                        try:
                            out, X, Y, Z = calc_orb.TD(args, data_for_discretization, et_transitions, grid_step=step, nproc=nproc)
                        except MemoryError :
                            sys.stderr.write('\n\nERROR: Memory Exception during discretization of MO used in the transitions\n')
                            et_transitions = []
                        if (len(et_transitions)> 0):
                            for k, transitions in enumerate(et_transitions):
                                if restart == 0 :
                                    print("EDD visualization in progress for the transition:", k+1)
                                    visu_mayavi.viz_EDD([out[k][0]], X, Y, Z, data_for_discretization, et_sym[k], 
                                                                         file_name="img", labels=[k+1]) 
                                    if (et_oscs[k] > 0.1) or ((len(et_rotats) > 0) and (abs(et_rotats[k]) > 10.)): 


                                        if (len(et_rotats) > 0) and (abs(et_rotats[k]) > 10.):
                                            #calculate only the elect and magnetic dipole for chiral compounds
                                            print("Generating overlap image for the selected transition:", k+1)
                                            visu_mayavi.viz_Oif([out[k][3]], X, Y, Z, data_for_discretization, et_sym[k], 
                                                                                 file_name="img", labels=[k+1])                                     
                                            if "et_magdips" in jf[i]["results"]["excited_states"]:
                                                et_magdips = (np.array(jf[i]["results"]["excited_states"]["et_magdips"][k]), np.array([0,0,0]))

                                        # Get overlap transition dipoles
                                        O_dip = (np.array(out[k][4][1]),  np.array(out[k][4][0]))

                                        print("Generating transition dipole images for selected transition:", k+1),
                                        ct_dip = out[k][2][3:]
                                        data_dip = {"GSDIP" : gs_eldip, "OVDIP" : O_dip, "CTDIP" : ct_dip}
                                        # In Gaussian Transition dispoles have for origin the center of masses. 
                                        if "et_eldips" in jf[i]["results"]["excited_states"]:
                                            et_eldips = (np.array(jf[i]["results"]["excited_states"]["et_eldips"][k]), np.array([0,0,0]))
                                            data_dip["ELDIP"] = et_eldips
                                        if "et_veldips" in jf[i]["results"]["excited_states"]:
                                            et_veldips = (np.array(jf[i]["results"]["excited_states"]["et_veldips"][k]), np.array([0,0,0]))
                                            #data_dip["VELDIP"] = et_veldips # Disabled for now
                                        visu_mayavi.viz_dip(data_dip, data_for_discretization, et_sym[k], file_name="img", labels=[k+1])
  
                                
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
                if report_type != 'text':
                    if discret_proc is True :
                        print("Calculating the emission electronic density difference.")
                        ## Discretization of all MO used in the transition
                        out, X, Y, Z = calc_orb.TD(args, data_for_discretization, [emi_transition], grid_step=step, nproc=nproc)
                        print("EDD visualization in progress for the transition :", emi_state)
                        visu_mayavi.viz_EDD([out[0][0]], X, Y, Z, data_for_discretization, emi_sym, 
                                                                 file_name="img-emi", labels=[emi_state]) 
                        ct_dip = out[0][2][3:]
                        data_dip = {"CTDIP" : ct_dip}
                        if (emi_rotat != 0.0):
                            visu_mayavi.viz_dip(data_dip, data_for_discretization, emi_sym, file_name="img-emi", labels=[emi_state])
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
                    heights = [[emi_osc*2.174e8/FWHM]]
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

            
    ## FUKUI functions processing, only in Full mode  
    if report_type == 'full':     
        # Determine Fukui oxydation states
        charge_SPm = charge_ref +1
        charge_SPp = charge_ref -1
        
        # Test on charges and Treatment of SP_plus
        #the number of SP+ is determined by the count of charges.
        #the index of the last in the list of charges is taken to be reused as the index of the jsonfile in the jsonfile list.
        
        if charges.count(charge_SPp) > 0:
            print(charges.count(charge_SPp),"Fukui reduced state detected. Only the last one will be considered.")
            SPp_index = charges.index(charge_SPp)
            #The presence of Mulliken partial charges is tested in order to process the calculation of CDFT indices    
            try : Mpc_p = jf[SPp_index]["results"]["wavefunction"]["Mulliken_partial_charges"]
            except KeyError :
                Mpc_p = []
            if len(Mpc_p) > 0 :
                A, fplus_lambda_mulliken, fplus_lambda_hirshfeld = calc_orb.CDFT_plus_Indices(data_for_discretization,jf[SPp_index])
                data_for_discretization["results"]["wavefunction"]["A"] = A
                data_for_discretization["results"]["wavefunction"]["fplus_lambda_mulliken"] = fplus_lambda_mulliken
                data_for_discretization["results"]["wavefunction"]["fplus_lambda_hirshfeld"] = fplus_lambda_hirshfeld
            # Proceed with Fukui dicretization
            delta_rho_SPp, X, Y, Z = calc_orb.Fukui(data_for_discretization,jf[SPp_index], label=None, grid_step=step, nproc=nproc)
            visu_mayavi.viz_Fukui(delta_rho_SPp, X, Y, Z, data_for_discretization, file_name="img", labels="SP_plus")
        
        # Test on charges and Treatment of SP_minus
        #the number of SP- is determined by the count of charges.
        #the index of the last in the list of charges is taken to be reused as the index of the jsonfile in the jsonfile list.
        if charges.count(charge_SPm) > 0 :
            print(charges.count(charge_SPm),"Fukui oxydized state detected. Only the last one will be considered.")
            SPm_index = charges.index(charge_SPm)
            delta_rho_SPm, X, Y, Z = calc_orb.Fukui(data_for_discretization,jf[SPm_index], label=None, grid_step=step, nproc=nproc)   
            visu_mayavi.viz_Fukui(delta_rho_SPm, X, Y, Z, data_for_discretization, file_name="img", labels="SP_minus")
        #The presence of Mulliken partial charges is tested in order to process the calculation of CDFT indices    
            try : Mpc_m = jf[SPm_index]["results"]["wavefunction"]["Mulliken_partial_charges"]
            except KeyError :
                Mpc_m = []
            if len(Mpc_m) > 0 :       
                I, fminus_lambda_mulliken, fminus_lambda_hirshfeld = calc_orb.CDFT_minus_Indices(data_for_discretization, jf[SPm_index])
                data_for_discretization["results"]["wavefunction"]["I"] = I
                data_for_discretization["results"]["wavefunction"]["fminus_lambda_mulliken"] = fminus_lambda_mulliken
                data_for_discretization["results"]["wavefunction"]["fminus_lambda_hirshfeld"] = fminus_lambda_hirshfeld
       
        # If both SP_plus and SP_minus are present. Treatment of Dual Descriptor           
        if (charges.count(charge_SPp) > 0) and (charges.count(charge_SPm) > 0) :
            delta_rho_dual, X, Y, Z = calc_orb.Fdual(data_for_discretization, delta_rho_SPp, delta_rho_SPm, grid_step=step)
            visu_mayavi.viz_Fdual(delta_rho_dual, X, Y, Z, data_for_discretization, file_name="img")
            if (len(Mpc_p) > 0) and (len(Mpc_m) > 0)  :
                A, I, Khi, Eta, Omega, DeltaN, fplus_lambda_mulliken, fminus_lambda_mulliken, fdual_lambda_mulliken, fplus_lambda_hirshfeld, fminus_lambda_hirshfeld, fdual_lambda_hirshfeld = calc_orb.CDFT_Indices(data_for_discretization, jf[SPp_index],jf[SPm_index])
                data_for_discretization["results"]["wavefunction"]["A"] = A
                data_for_discretization["results"]["wavefunction"]["I"] = I
                data_for_discretization["results"]["wavefunction"]["Khi"] = Khi
                data_for_discretization["results"]["wavefunction"]["Eta"] = Eta
                data_for_discretization["results"]["wavefunction"]["Omega"] = Omega
                data_for_discretization["results"]["wavefunction"]["DeltaN"] = DeltaN
                data_for_discretization["results"]["wavefunction"]["fplus_lambda_mulliken"] = fplus_lambda_mulliken
                data_for_discretization["results"]["wavefunction"]["fminus_lambda_mulliken"] = fminus_lambda_mulliken
                data_for_discretization["results"]["wavefunction"]["fdual_lambda_mulliken"] = fdual_lambda_mulliken
                data_for_discretization["results"]["wavefunction"]["fplus_lambda_hirshfeld"] = fplus_lambda_hirshfeld
                data_for_discretization["results"]["wavefunction"]["fminus_lambda_hirshfeld"] = fminus_lambda_hirshfeld
                data_for_discretization["results"]["wavefunction"]["fdual_lambda_hirshfeld"] = fdual_lambda_hirshfeld

            
    # Autocrop generated images
    ImagesToCrop = [file for file in os.listdir(".") if file.endswith(".png")]
    for i, file in enumerate(ImagesToCrop):
        if False:
            print("Autocrop image:", file)
        else:
            end = '\r' if i+1 != len(ImagesToCrop) else '\n'
            print("Cropping images... (%d/%d)" % (i+1, len(ImagesToCrop)), end=end)
        autocrop(file)
