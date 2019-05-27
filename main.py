#!/usr/bin/python -W ignore
#-*- coding: utf-8 -*-

# QuchemReport reads a list of quantum chemistry computational results (log files)
# After a quality check, it converts the extracted information in JSON files
# QuchemReport can read Single Points, Optimization, Frequency and Time Dependent calculations
# The properties are automatically processed and put in a common report 

# Required informations :
# Ground state should be provided by at least one OPT or Freq 
# Geometry (atomic positions) has to be present for 3D molecular views (Topo)    
# Molecular coefficients (MO coeff) and basis sets have to be present in a log for discretization process
# TD, Freq, OPT_ES, Freq_ES job types means a spectrum (UV or IR) process
# TD and MO or OPT_ES/Freq_ES and MO launch an electron density differences (EDD) discretization process
# Different charges between oxydation states launch a Fukui process and EDD

# Conformity check :
# All log files must correspond to the same molecular entity / conformer
# That means the same formula and geometry for ground states and oxydation states
# Optimized parent excited states like S1, T1 are also accepted.

import os
import sys
import tempfile
import numpy as np
import math
import matplotlib.pyplot as plt
#  import modules
from quchemreport.parser import process_logfile_list
from quchemreport.processing import calc_orb
from quchemreport.externals import  GaussSum_IRspectrum, GaussSum_UVspectrum
from quchemreport.visualization import visu_mayavi, json2latex2pdf

## Usage
usage = ("Usage: {0}  list_of_computational_output_log_files").format(sys.argv[0])

# Create a temporary directory to store intermediate log file, etc.
tmpdirname = tempfile.mkdtemp()    
print('Creating temporary directory:', tmpdirname)

# Read a list of log_files, parse them with quality check and return a list of JSON 
# lf is the list of parsed log files in temp directory
# jf is the list of created json  
lf, jf = process_logfile_list(sys.argv[1:], log_storage_path=tmpdirname, verbose=False)

# checking that jf is not empty before any processing.
if (len(jf)) == 0:
     print("There is no valid log file. Program will exit.")
     sys.exit()
print(len(jf), " valid log file detected.")

### CONFORMITY/UNIFORMITY TESTS
# Conformity and Uniformity between log files depend on chosen JSON keys
# All log files must correpond to the same formula and theory (functional and basis set)     
# Same formulas, same theory and same NRE (nuclear repulsion energy) == same conformer
# If same conformer but different charges = Fukui process
# If same conformer but different NRE and job type = OPT_ES or FREQ_ES, related excited state

formulas =[]
functionals = []
ao_names = []
basis_sets = []
MO_coeff  = []
job_types = []
nres = []
charges = []
multiplicities = []
# NRE without the Excited states because they will have a different NRE. 
nres_noES = []
# Discard SP charge and multiplicity because oxydation states will differ from neutral state.
charges_noSP = []
multiplicities_noSP = []
charges_SP = []

# Create list of key values for conformity tests
for i, jsonfiles in enumerate(jf):
    formulas.append(jf[i]["molecule"]["formula"])
    nres.append(jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"])
    if ((jf[i]["comp_details"]["general"]["job_type"] == ['OPT_ES']) or
        (jf[i]["comp_details"]["general"]["job_type"] == ['FREQ_ES'])) :
        print ('Detected optimized excited state in :', os.path.basename(lf[i]))
    else :
        nres_noES.append(jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"])
    functionals.append(jf[i]["comp_details"]["general"]["functional"])
    # basis sets, ao_names and MO coeffs are not mandatory but necessary for the discretization process
    try : 
         basis_sets.append(jf[i]["comp_details"]["general"]["basis_set"])
    except KeyError :
        pass 
    try : 
         ao_names.append(jf[i]["comp_details"]["general"]["ao_names"])
    except KeyError :
        pass 
    try : 
         MO_coeff.append(jf[i]["results"]["wavefunction"]["MO_coefs"]) 
    except KeyError : 
        MO_coeff.append([])   # Problem we would like to add N/A is no MO to keep the same index hax jf[i] but test on len after   
    charges.append(jf[i]["molecule"]["charge"])
    if (jf[i]["comp_details"]["general"]["job_type"] != ['SP'])  :
        charges_noSP.append(jf[i]["molecule"]["charge"])
    if (jf[i]["comp_details"]["general"]["job_type"] == ['SP'])  :
        charges_SP.append(jf[i]["molecule"]["charge"])    
    job_types.append(jf[i]["comp_details"]["general"]["job_type"]) 
    multiplicities.append(jf[i]["molecule"]["multiplicity"])
    if (jf[i]["comp_details"]["general"]["job_type"] != ['SP'])  :
         multiplicities_noSP.append(jf[i]["molecule"]["multiplicity"])

# Problem in the formula: the charges appear and direct comparison fails
# example : [u'C4H7NO', u'C4H7NO', u'C4H7NO+', u'C4H7NO', u'C4H7NO-']
# Clean the charges in the formula to compare them later
formulas_nocharges = []
# If no charge, get formula, 
# If monocationic remove last caracter (+) because 1 is implicit
# If charge > 1, remove length of charge and the plus sign
# If charge is negative, length of charge include the minus sign

for i, charge in enumerate(charges) :
    if (int(charges[i]) == 0) :
        formulas_nocharges.append(formulas[i])
    elif (int(charges[i]) == 1) or (int(charges[i]) < -1) :
        formulas_nocharges.append(formulas[i][:-len(str(charges[i]))])
    elif (int(charges[i]) > 1) :
        formulas_nocharges.append(formulas[i][:-(len(str(charges[i]))+1)])
    elif (int(charges[i]) == -1) :
        formulas_nocharges.append(formulas[i][:-(len(str(charges[i]))-1)])
        
# Checking if jf has more than 1 file in order to process parenting and uniformity test
if (len(jf)) > 1 :
    # Test on cleaned formulas. There should be only one value
    if (len(set(formulas_nocharges))) > 1:
        print("You need to provide a set of logfiles for a unique conformer : Too many different formulas")
        # print log files names and corresponding formulas cleaned of charges
        # print(os.path.basename(lf)) Print of list of log names does not work !
        print(formulas_nocharges)
        print("Program will exit.")
        sys.exit()    
    elif (len(set(formulas_nocharges))) == 1: 
        print("All logfiles correspond to the same formula (cleaned of charges). Test OK")
        print("Detected formula:", set(formulas_nocharges))
    elif (len(set(formulas_nocharges))) == 0:
        print("There is no correct formulas in the list. Program will exit. Contact admin")
        sys.exit()
    
    # Tests on theory (same functional and same basis set)
    #print("Detected number of basis sets : ", len(set(basis_sets)))
    #print("Detected number of functionals : ", len(set(functionals)))
    
    if (len(set(functionals)) > 1) or (len(set(basis_sets)) > 1): 
        print("You need to provide a set of logfiles with the same level of theory")
        print ("That means a unique functional and basis set.")
        # print log files names and corresponding theories
        # print(os.path.basename(lf)) Print of list of log names does not work !
        print("Detected number of basis sets: ", len(set(basis_sets)))
        print("Detected functionals: ",functionals)
        print("Program will exit.")
        sys.exit()  
    else :
        print("All logfiles present the same functional and basis set. Test OK")
        print("Detected functional : ", set(functionals))
        print("Detected number of basis sets : ", len(set(basis_sets)))
    
    # Test on NRE to ensure that there is only one conformer
    if (len(set(nres_noES))) > 1:
        print("You need to provide a set of logfiles for a unique conformer without counting the optimized excited states.")
        # print log files names and corresponding NRE
        print("Detected number of geometries (no optimized excited state): ", len(set(nres_noES)))
        # print(os.path.basename(lf)) Print of list of log names does not work !
        print("Detected NRE (with optimized excited state): ",nres)
        print("Program will exit.")
        sys.exit()
    else :
        print("All logfiles correspond to the same conformer. Test OK")
    
    # Test on charge to ensure that there is only one oxydation state whitout counting the single points
    if (len(set(charges_noSP))) > 1:
        print("You need to provide a set of logfiles with the same charge without counting the single points.")
        # print log files names and corresponding charges
        print("Detected number of charges (No single points): ", len(set(charges_noSP)))
        # print(os.path.basename(lf)) Print of list of log names does not work ! 
        print("Detected charges (with single points): ",charges)
        print("Program will exit.")
        sys.exit()
    elif  (len (set(charges_noSP))) == 0 :
        print("There is no correct charges in the list. Program will exit. Contact admin")
        sys.exit()
    else :
        print("All logfiles present the same charge. Test OK")
    
    # Test on multiplicity to ensure that there is only one oxydation state whitout counting the single points
    if (len(set(multiplicities_noSP))) > 1:
        print("You need to provide a set of logfiles with the same electronic state.")
        # print log files names and corresponding multiplicity
        print("Detected number of multiplicities (No single points):", len(set(multiplicities_noSP)))
        # print(os.path.basename(lf)) Print of list of log names does not work !
        print("Detected multiplicities (with single points): ",multiplicities)
        print("Program will exit.")
        sys.exit()
    elif  (len(set(multiplicities_noSP))) == 0 :
        print("There is no correct multiplicity in the list. Program will exit. Contact admin")
        sys.exit()
    else :
        print("All logfiles present the same multiplicity. Test OK")
 
    # Test on Ground State information in case of multiples log files
    if job_types.count(['OPT']) == 0 and job_types.count(['FREQ']) == 0 and job_types.count(['FREQ', 'OPT']) == 0:
        print("In a list of logfiles, one OPT or FREQ of a Ground State should be provided. Program will exit.")
        sys.exit()
   
    
# End of uniformity / conformity tests. 
# At this point : we have a set of log files of a unique conformer
# Or related excited state or related Single point
        
# PARENTING AND INFORMATION TESTS
### DISCRETIZATION TESTS

#Setting a json file as reference for discretization informations like MO coefficients and basis set.      
charge_ref = charges_noSP[0]
multiplicity_ref = multiplicities_noSP[0]
MO_coeff_ref = [] 

if (len(basis_sets)) == 0:
    basis_sets_ref = []
else:
    basis_sets_ref = basis_sets[0]
if (len(ao_names)) == 0:
    ao_names_ref = []
else:
    ao_names_ref = ao_names[0]

for i, charge in enumerate(charges):
    if jf[i]["comp_details"]["general"]["is_closed_shell"] is 'False':
        if (charges[i] == charge_ref and  multiplicities[i] == multiplicity_ref and
            nres[i] == nres_noES[0] and (len(MO_coeff[i])) == 2):  
            jf_ref = jf[i]
            MO_coeff_ref = MO_coeff[i] 
            data_for_discretization = jf[i]
            print("A reference json file has been selected for the discretization process:", lf[i])
    else :
        if charges[i] == charge_ref and  multiplicities[i] == multiplicity_ref  \
        and nres[i] == nres_noES[0] and (len(MO_coeff[i])) == 1 :  
            jf_ref = jf[i]
            MO_coeff_ref = MO_coeff[i] 
            data_for_discretization = jf[i]
            print("A reference json file has been selected for the discretization process:", lf[i])

# Discretization depends on Orbkit and needs basis sets and MO coefficients.
# Orbkit can not work if the file have both sphericals and cartesians basis set functions.
discret_proc = False
cartesian = False
spherical = False
        
# Test on Basis set and MO coefficients to ensure that the discretization process is doable
# Than tests if basis set is cartesian or spherical or mixed
if (len(set(basis_sets_ref))) == 0 or (len(MO_coeff_ref)) == 0 :
    print("There is no basis set or MO coefficients in your logfiles. MO and EDD pictures cannot be generated.")
    discret_proc = False
elif (len(set(basis_sets_ref))) != 0 and (len(MO_coeff_ref)) != 0:
    print("MO coefficients and basis sets detected." )
    # Basis set names for D and F orbitals
    ao = []
    for i, orbitals in enumerate(ao_names_ref) :
        ao.append(ao_names_ref[i].split('_')[-1])
    ao_DF = []
    for i, orbitals in enumerate(ao) :
        # Isolate D and F atomic orbitals names
        # If cartesian, D may not appear in the name. We keep  the X coordinates for test
        if ("D" in ao[i]) or ("XX" in ao[i]) or ("F" in ao[i]) :
            ao_DF.append(ao[i])    
    for i, orbitals in enumerate(ao_DF) :        
        if ("XX" in ao_DF[i]) :
            cartesian = True
        elif ("+" in ao_DF[i]) :
            spherical = True
    if (cartesian is True) and (spherical is True) :
        discret_proc = False      
    elif (cartesian is True) and (spherical is False) :
        discret_proc = True
        print("Cartesian basis set detected")
    elif (cartesian is False) and (spherical is True) :
        discret_proc = True
        print("Spherical basis set detected")
    else :
        print("The basis set is neither cartesian nor spherical nor mixed. Error contact admin")
        sys.exit()
    if discret_proc is True :
        print("All necessary information for orbkit are present. Discretization process can be done.")
    else :
        print("The basis set is a mixture of cartesian and spherical atomic orbitals. Discretization cannot be done.")
                     
       
# SETUP for ORBKIT
## Get grid parameters and initialize grid
## Oversizing the grid aroung the molecule in Bohr radii
## ORBKIT uses 5 by default, tune this as required
over_s = 7
## Spacing/Number of points
par = -6
MO_list = []

Elec_T = []  
# electronic transitions

for i, jt in enumerate(job_types):	  
    if job_types[i] == ['OPT']  or job_types[i] == ['FREQ', 'OPT'] :
	        # Test if Ground state by geometry, if atomic positions present: generate picture
	        if jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"] == nres_noES[0] \
	        and jf[i]["results"]["geometry"]["elements_3D_coords_converged"] != 'N/A': 
	            visu_mayavi.topo(jf[i],file_name="img")
	            print("Picture generated for the topology of the ground state")
	        # Test if Ground state by geometry and by charge and if MO coeffs presents 
	        # So generate HOMO and LUMO pictures for the ground state at the reference oxydation state
	        if jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"] == nres_noES[0] \
	        and int(jf[i]["molecule"]["charge"]) == charge_ref \
	        and (discret_proc is True) :
	            # SETUP for the HOMO and LUMO (HOMO+1) index list
                #Orbkit starts MO_list counter from one
                 #MO_list = [jf[i]["results"]["wavefunction"]["homo_indexes"]+1] + [jf[i]["results"]["wavefunction"]["homo_indexes"]+2] 
#                if len(jf[i]["results"]["wavefunction"]["homo_indexes"]) == 2:
#                    MO_list_alpha = [str(jf[i]["results"]["wavefunction"]["homo_indexes"][0]+1)+"_a"] + [str(jf[i]["results"]["wavefunction"]["homo_indexes"][0]+2)+"_a"] 
#        	            MO_list_beta = jf[i]["results"]["wavefunction"]["homo_indexes"][1]+1 + jf[i]["results"]["wavefunction"]["homo_indexes"][1]+2 
                 MO_list = jf[i]["results"]["wavefunction"]["homo_indexes"] + [x+1 for x in jf[i]["results"]["wavefunction"]["homo_indexes"]] 
    	           # MO_list = ['homo', 'lumo']
                    ## Calculations of MO
                 out, X, Y, Z = calc_orb.MO(jf[i], MO_list, grid_par=par)
                 for series in out:
	                ## The length product works because all voxels of the ORBKIT grid have the same dimensions
    	                print ("MO norm by numpy :", np.sum(np.square(series))*(X[1,0,0] - X[0,0,0])*(Y[0,1,0] - Y[0,0,0])*(Z[0,0,1] - Z[0,0,0]))
	            ## Visulation of the MO
                 visu_mayavi.viz_MO(out, X, Y, Z, jf[i], file_name="img", labels=MO_list)
               
    if  job_types[i] == ['FREQ']  or job_types[i] == ['FREQ', 'OPT'] :
        print("Frequency calculation detected. Calculating the IR Spectrum.")
        # Calculating of IR spectrum
        GaussSum_IRspectrum.Vibfreq(jf[i],os.path.basename(lf[i]),600,4000,1700,4.0,"Gen",1,0,jf[i]["comp_details"]["freq"]["temperature"])
       # Plotting IR spectrum 
       #extracts the datas from the resulting txt file to generate a png file of the spectrum
        x,y = np.loadtxt("IRSpectrum.txt",skiprows=2,usecols=(0,1), unpack=True)
        plt.xlabel('Wavenumber (cm-1)')
        plt.ylabel('IR activity')
        plt.title(os.path.basename(lf[i]) + "_IR_spectrum")
        plt.plot(x,y, 'o', color='black',linestyle="dashdot",markersize= 2 )
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        plt.savefig(os.path.basename(lf[i]) + "_IR_spectrum.png")
        plt.show()
    if job_types[i] == ['TD'] :
         if jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"] == nres_noES[0] \
	        and int(jf[i]["molecule"]["charge"]) == charge_ref :                
            print("Electronic transitions detected. Calculating the UV absorption spectrum.")
             #end wave number and start wave number are in nm. The externals GaussSum_UVspectrum use the function convertor of the parser cclib to convert them in cm-1      
             # Get enegyrange in nm, convert in cm-1.
             # Spectrum width should min cm-1 - 9000 (min 5000), max cm-1 + 9000 (max 100000 cm-1)
            start =  min(jf[i]["results"]["excited_states"]["et_energies"]) - 9000
            end = max(jf[i]["results"]["excited_states"]["et_energies"]) + 9000
            if start < 5000:
                start = 5000
            if end > 100000 :
                end = 100000
            numpts = math.floor((end - start)/20)
            if ((len(jf[i]["results"]["excited_states"]["et_energies"])) > 0 ) and  ((len(jf[i]["results"]["excited_states"]["et_oscs"])) > 0 ):
                    try:
                        et_rotats=jf[i]["results"]["excited_states"]["et_rot"]
                    except KeyError:
                        et_rotats = []
                    GaussSum_UVspectrum.ET(os.path.basename(lf[i]),jf[i]["results"]["excited_states"]["et_energies"],\
                                           jf[i]["results"]["excited_states"]["et_oscs"], \
                                           et_rotats,start,end,numpts,3000)
            print("Absorption spectrum done.")
            # Plotting UV spectrum 
            #extracts the datas from the resulting txt file to generate a png file of the spectrum
            x,y = np.loadtxt("UVSpectrum.txt",skiprows=1,usecols=(1,2), unpack=True)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Abs')
            plt.title(os.path.basename(lf[i]) + "_UV spectrum")
            plt.plot(x,y, 'o', color='black', markersize= 2)
            plt.savefig(os.path.basename(lf[i]) + "_UV spectrum.png")
            plt.show()
            if discret_proc is True :
                print("Calculating the electronic density differences.")
                # Set the transition list to process = T_list. By default only the first transtion is taken into account.
                T_list = [0,1,2]
                for k, transitions in enumerate(jf[i]["results"]["excited_states"]["et_transitions"]):
                    transitions = [jf[i]["results"]["excited_states"]["et_transitions"][k] for k in T_list]
                print("EDD in progress for the transition(s):", T_list)
                    ## Calculations
                out, X, Y, Z = calc_orb.TD(data_for_discretization, transitions, grid_par=par) 
                visu_mayavi.viz_EDD([e[0] for e in out], X, Y, Z, data_for_discretization, file_name="img", labels=[x+1 for x in T_list])
                print("EDD done")
                ## Returns the calculated values of the tozer_lambda, d_CT, Q_CT, Mu_CT and e- barycenter and hole barycenter to the json
                jf[i]["results"]["excited_states"]["Tozer_lambda"] = []
                jf[i]["results"]["excited_states"]["d_ct"] = []
                jf[i]["results"]["excited_states"]["q_ct"] = []
                jf[i]["results"]["excited_states"]["mu_ct"] = []
                jf[i]["results"]["excited_states"]["e-_barycenter"] = []
                jf[i]["results"]["excited_states"]["hole_barycenter"] = []
                for e in out:
                    jf[i]["results"]["excited_states"]["Tozer_lambda"].append(e[1])
                    jf[i]["results"]["excited_states"]["d_ct"].append(e[2])
                    jf[i]["results"]["excited_states"]["q_ct"].append(e[2][1])
                    jf[i]["results"]["excited_states"]["mu_ct"].append(e[2][2])
                    jf[i]["results"]["excited_states"]["e-_barycenter"].append(e[2][3].tolist())
                    jf[i]["results"]["excited_states"]["hole_barycenter"].append(e[2][4].tolist())
            
            # for l in [(e[1], e[2]) charges.count(charge_SPm)for e in out]:
	    #     print l
            #with open("./test_bdm.json", "w+") as f:
         #   with open(argv[4], "w+") as f:
          #      dump(TD_data, f)
                
    if  jf[i]["comp_details"]["general"]["job_type"] == ['OPT_ES'] :
        print("Optimization of an excitated state detected. Calculating the UV emission spectrum")
        # TODO : get optimized excited state number and isolate the et_energy and osc_strgth       
        #Ploting of UV emission spectrum         
        print("electronic transition detected. Calculating the UV absorption Spectrum")
        GaussSum_UVspectrum.GaussianSpectrum(200,900,1300,
                             (jf[i]["results"]["excited_states"]["et_energies"] ,[[x*2.174e8/3000 for x in jf[i]["results"]["excited_states"]["et_oscs"] ]] ),
                             3000)
        print("UV emission spectrum done")
        if discret_proc is True : 
                T_list = [0]
                for k, transitions in enumerate(jf[i]["results"]["excited_states"]["et_transitions"]):
                    transitions = [jf[i]["results"]["excited_states"]["et_transitions"][k] for k in T_list]
                print("EDD in progress")
                    ## Calculations
                out, X, Y, Z = calc_orb.TD(data_for_discretization, transitions, grid_par=par) 
                visu_mayavi.viz_EDD([e[0] for e in out], X, Y, Z, data_for_discretization,file_name="img",  labels=T_list)
                print("EDD done")
                     
## FUKUI functions processing
# Determine Fukui oxydation states
charge_SPm = charge_ref +1
charge_SPp = charge_ref -1

# Test on charges and Treatment of SP_plus
if charges.count(charge_SPp) > 0:
    print(charges.count(charge_SPp),"Fukui reduced state detected. Only the last one will be considered.")
    SPp_index = charges.index(charge_SPp)
    delta_rho_SPp, X, Y, Z = calc_orb.Fukui(data_for_discretization,jf[SPp_index], label=None, grid_par=par)
    visu_mayavi.viz_Fukui(delta_rho_SPp, X, Y, Z, data_for_discretization, file_name="img", labels="SP_plus")
    try : Mpc_p = jf[SPp_index]["results"]["wavefunction"]["Mulliken_partial_charges"]
    except KeyError :
       Mpc_p = []
    if len(Mpc_p) > 0 :
       calc_orb.CDFT_plus_Indices(data_for_discretization,jf[SPp_index])

# Test on charges and Treatment of SP_minus
if charges.count(charge_SPm) > 0 :
    print(charges.count(charge_SPm),"Fukui oxydized state detected. Only the last one will be considered.")
    SPm_index = charges.index(charge_SPm)
    delta_rho_SPm, X, Y, Z = calc_orb.Fukui(data_for_discretization,jf[SPm_index], label=None, grid_par=par)   
    visu_mayavi.viz_Fukui(delta_rho_SPm, X, Y, Z, data_for_discretization, file_name="img", labels="SP_minus")
    try : Mpc_m = jf[SPm_index]["results"]["wavefunction"]["Mulliken_partial_charges"]
    except KeyError :
        Mpc_m = []
    if len(Mpc_m) > 0 :       
       if charges.count(charge_SPp) > 0  :
           calc_orb.CDFT_minus_Indices(data_for_discretization, jf[SPm_index])

# If both SP_plus and SP_minus are present. Treatment of Dual Descriptor           
if (charges.count(charge_SPp) > 0) and (charges.count(charge_SPm) > 0) :
    delta_rho_dual, X, Y, Z = calc_orb.Fdual(data_for_discretization, jf[SPp_index],jf[SPm_index] , grid_par=-6)
    visu_mayavi.viz_Fdual(delta_rho_dual, X, Y, Z, data_for_discretization, file_name="img")
    if (len(Mpc_p) > 0) and (len(Mpc_m) > 0)  :
        calc_orb.CDFT_Indices(data_for_discretization, jf[SPp_index],jf[SPm_index])



##Report Generation

for i,jsonfile in enumerate(jf):
    json2latex2pdf.json2latex2pdf(jf[i],mode="clean")  
  
    
         

 

        
    
                
                




