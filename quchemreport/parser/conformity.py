#-*- coding: utf-8 -*-

### CONFORMITY/UNIFORMITY TESTS
# Conformity and Uniformity between log files depend on chosen JSON keys
# All log files must correpond to the same formula and theory (functional and basis set)
# Same formulas, same theory and same NRE (nuclear repulsion energy) == same conformer
# If same conformer but different charges = Fukui process
# If same conformer but different NRE and job type = OPT_ES or FREQ_ES, related excited state

import sys

def tests(args, jf):
# Parameters list of logfiles and JSON data

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

    verbose = args['verbose']

    # Create list of key values for conformity tests
    for i, jsonfiles in enumerate(jf):
        package =  jf[i]["comp_details"]["general"]["package"]
        formulas.append(jf[i]["molecule"]["formula"])
        nres.append(jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"])
        if ((jf[i]["comp_details"]["general"]["job_type"] == ['OPT_ES']) or
            (jf[i]["comp_details"]["general"]["job_type"] == ['FREQ_ES'])) :
            print ('Detected optimized excited state in :', jf[i]['metadata']['log_file'])
        else :
            nres_noES.append(jf[i]["results"]["geometry"]["nuclear_repulsion_energy_from_xyz"])
        try :
            functionals.append(jf[i]["comp_details"]["general"]["functional"])
        except KeyError :
            functionals.append("ERROR reading functionals")
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
            print("ERROR: You need to provide a set of logfiles for a unique conformer : Too many different formulas")
            print(formulas_nocharges)
            print("Program will exit.")
            sys.exit()    
        elif (len(set(formulas_nocharges))) == 1: 
            print("\tDetected formula:                                ... ", set(formulas_nocharges))
            print("\tSame formula for all logfiles:                   ...  Test OK")
        elif (len(set(formulas_nocharges))) == 0:
            print("ERROR: There is no correct formulas in the list. Program will exit. Contact admin")
            sys.exit()

    # Tests on theory (same functional and same basis set)
    #print("Detected number of basis sets : ", len(set(basis_sets)))
    #print("Detected number of functionals : ", len(set(functionals)))

    if (len(set(functionals)) > 1) or (len(set(basis_sets)) > 1): 
        print("\tDetected number of basis sets:                   ... ", len(set(basis_sets)))
        print("\tDetected functionals:                            ... ", functionals)
        print("ERROR: You need to provide a set of logfiles with the same level of theory")
        print("That means a unique functional and basis set.")
        print("Program will exit.")
        sys.exit()  
    else:
        print("\tDetected number of basis sets:                   ... ", len(set(basis_sets)))
        print("\tDetected functionals:                            ... ", functionals)
        print("\tSame functional and basis set for all logfiles:  ...  Test OK")
    # Test on NRE to ensure that there is only one conformer
    if (len(set(nres_noES))) > 1:
        print("\tSame geometries (no optimized excited state):    ... ", len(set(nres_noES)))
        print("\tDetected NRE (with optimized excited state):       ... ",nres)
        print("ERROR: You need to provide a set of logfiles for a unique conformer without counting the optimized excited states.")
        print("Program will exit.")
        sys.exit()
    else:
        print("\tSame conformer for all logfiles:                 ...  Test OK")

    # Test on charge to ensure that there is only one oxydation state whitout counting the single points
    if (len(set(charges_noSP))) > 1:
        print("\tDetected charges (No single points):             ... ", len(set(charges_noSP)))
        print("\tDetected charges (with single points):           ... ",charges)
        print("ERROR: You need to provide a set of logfiles with the same charge without counting the single points.")
        print("Program will exit.")
        sys.exit()
    elif  (len (set(charges_noSP))) == 0 :
        print("ERROR: There is no correct charges in the list. Program will exit. Contact admin")
        sys.exit()
    else:
        print("\tSame charge for all logfiles:                    ...  Test OK")

    # Test on multiplicity to ensure that there is only one oxydation state whitout counting the single points
    if (len(set(multiplicities_noSP))) > 1:
        print("\tDetected multiplicities (No single points):      ... ", len(set(multiplicities_noSP)))
        print("\tDetected multiplicities (with single points):    ... ",multiplicities)

        print("ERROR: You need to provide a set of logfiles with the same electronic state.")
        print("Program will exit.")
        sys.exit()
    elif  (len(set(multiplicities_noSP))) == 0 :
        print("ERROR: There is no correct multiplicity in the list. Program will exit. Contact admin")
        sys.exit()
    else:
        print("\tSame multiplicity:                               ...  Test OK")

    # Test on Ground State information in case of multiples log files
    if job_types.count(['OPT']) == 0 and job_types.count(['FREQ']) == 0 and \
       job_types.count(['FREQ', 'OPT']) == 0 and job_types.count(['FREQ', 'OPT', 'TD']) == 0:
        print("ERROR: In a list of logfiles, one OPT or FREQ of a Ground State should be provided. Program will exit.")
        sys.exit()

    
    # End of uniformity / conformity tests. 
    # At this point : we have a set of log files of a unique conformer
    # Or related excited state or related Single point
    print("All conformity tests OK.")

    
    # PARENTING AND INFORMATION TESTS
    ### DISCRETIZATION TESTS
    print("\nStarting discretization tests.")
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
        # Previously we tested on jf[i]["comp_details"]["general"]["is_closed_shell"] : Problem some closed shell with both alpha and beta
        if (len(MO_coeff[i]) == 1) or (len(MO_coeff[i]) == 2) : # Restricted or unrestricted calculation
            if charges[i] == charge_ref and  multiplicities[i] == multiplicity_ref  \
                                        and nres[i] == nres_noES[0]  :
                jf_ref = jf[i]
                MO_coeff_ref = MO_coeff[i]
                data_for_discretization = jf[i]
                if verbose:
                    print("A reference json file has been selected for the discretization process:", jf[i]['metadata']['log_file'])
                break
        else :
            data_for_discretization = []
            if verbose:
                print("No reference json file has been found for the discretization process due to MO_coeff length")
    
    # Discretization depends on Orbkit and needs basis sets and MO coefficients.
    # Orbkit can not work if the file have both sphericals and cartesians basis set functions.
    discret_proc = False
    cartesian = False
    spherical = False
    # Define booelan to not repeat MO discretization process
    mo_viz_done = False

    # Test on Basis set and MO coefficients to ensure that the discretization process is doable
    # Than tests if basis set is cartesian or spherical or mixed
    if (len(set(basis_sets_ref))) == 0 or (len(MO_coeff_ref)) == 0 :
        print("There is no basis set or MO coefficients in your logfiles. MO and EDD pictures cannot be generated.")
        discret_proc = False
    elif (len(set(basis_sets_ref))) != 0 and (len(MO_coeff_ref)) != 0:
        if verbose:
            print("MO coefficients and basis sets detected. MO and EDD pictures can be generated." )
        # Basis set names for D and F orbitals
        ao = []
        for i, orbitals in enumerate(ao_names_ref) :
            ao.append(ao_names_ref[i].split('_')[-1])
        ao_DF = []
        for i, orbitals in enumerate(ao) :
            # Isolate D and F atomic orbitals names
            # If cartesian, D may not appear in the name. We keep  the X coordinates for test
            if ("D" in ao[i]) or ("XX" in ao[i]) or ("F" in ao[i]) or ("XY" in ao[i]):
                ao_DF.append(ao[i])    
        for i, orbitals in enumerate(ao_DF) :        
            if ("XX" in ao_DF[i]):
                cartesian = True
            elif ("+" in ao_DF[i]) or package == 'ORCA':
                spherical = True

        # Test if there is no D or F functions
        if (len(ao_DF)) == 0 :
            for i, orbitals in enumerate(ao) :
            # Enumerate again to test if spherical or cartesian based on the p orbitals
                if ("PX" in ao[i]) :
                    cartesian = True
            # See which test could isolate a spherical p orbital
            #    else : 
            #        spherical = True
        if (cartesian is True) and (spherical is True) :
            discret_proc = False      
        elif (cartesian is True) and (spherical is False) :
            discret_proc = True
            if verbose:
                print("Cartesian basis set detected")
        elif (cartesian is False) and (spherical is True) :
            discret_proc = True
            if verbose:
                print("Spherical basis set detected")
        else :
            print("ERROR: The basis set is neither cartesian nor spherical nor mixed. Contact admin.")
            sys.exit()
        if discret_proc is True :
            print("All discretization tests OK.")
        else :
            print("The basis set is a mixture of cartesian and spherical atomic orbitals. Discretization cannot be done.")

    data = {'job_types' : job_types,
            'nres_noES' : nres_noES,
            'charges' : charges,
            'charge_ref' : charge_ref,
            'discret_proc' : discret_proc,
            'mo_viz_done' : mo_viz_done,
            'data_for_discretization' : data_for_discretization}
    return (data, job_types, nres_noES, charges, charge_ref, discret_proc, mo_viz_done, data_for_discretization)



