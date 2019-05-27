# -*- coding: iso-8859-1 -*-
'''Module for reading the output files of quantum chemical software. 
'''
'''
Modification of the cclib part of  the orbkit read.py
Gunter Hermann, Vincent Pohl, and Axel Schild

Institut fuer Chemie und Biochemie, Freie Universitaet Berlin, 14195 Berlin, Germany

This file is part of orbkit.

orbkit is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or any later version.

orbkit is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with orbkit.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
#import copy
import numpy

from orbkit.core import l_deg, lquant, create_mo_coeff
from orbkit.display import display
from orbkit.qcinfo import QCinfo, get_atom_symbol

def convert_json(jData, all_mo=False, spin=None):
  '''Converts a scanlog JSON data instance to an instance of
  orbkit's QCinfo class.

  **Parameters:**
  
    jData : class
      Contains the input JSON data.
    all_mo : bool, optional
      If True, all molecular orbitals are returned.
    spin : {None, 'alpha', or 'beta'}, optional
      If not None, returns exclusively 'alpha' or 'beta' molecular orbitals.


  **Returns:**

    qc (class QCinfo) with attributes geo_spec, geo_info, ao_spec, mo_spec, etot :
          See :ref:`Central Variables` for details.
  '''
  aa_to_au = 1/0.52917720859
  # Initialize the variables 
  qc = QCinfo()
  
  # Converting all information concerning atoms and geometry
  qc.geo_spec = numpy.array(jData['results']['geometry']['elements_3D_coords_converged']).reshape((-1, 3)) * aa_to_au
  for ii in range(jData["molecule"]['nb_atoms']):
    symbol = get_atom_symbol(atom=jData["molecule"]['atoms_Z'][ii])
    qc.geo_info.append([symbol,str(ii+1),str(jData["molecule"]['atoms_Z'][ii])])
  
  # Convert geo_info and geo_spec to numpy.ndarrays
  qc.format_geo()
  
  # Converting all information about atomic basis set
  from pickle import loads
  gbasis = loads(bytes(jData['comp_details']['general']['basis_set'], 'utf-8'))
  for ii in range(jData["molecule"]['nb_atoms']):
    for jj in range(len(gbasis[ii])):
      pnum = len(gbasis[ii][jj][1])
      qc.ao_spec.append({'atom': ii,
                  'type': str(gbasis[ii][jj][0]).lower(),
                  'pnum':  pnum,
                  'coeffs': numpy.zeros((pnum, 2))
                  })
      for kk in range(pnum):
        qc.ao_spec[-1]['coeffs'][kk][0] = gbasis[ii][jj][1][kk][0]
        qc.ao_spec[-1]['coeffs'][kk][1] = gbasis[ii][jj][1][kk][1]

  if "ao_names" in jData['comp_details']['general']:
    # Reconstruct exponents list for ao_spec
    aonames = jData['comp_details']['general']['ao_names']
    cartesian_basis = True
    for i in aonames:
      if '+' in i or '-' in i:
        cartesian_basis = False
# There is a problem here with the 6D 7F basis sets, that are a mixture of cartesian and spherical basis sets.
    if not cartesian_basis:
        qc.ao_spherical = []
    
    count = 0
    for i,ao in enumerate(qc.ao_spec):
      l = l_deg(lquant[ao['type']],cartesian_basis=cartesian_basis)
      if cartesian_basis:
        ao['exp_list'] = []
        
      for ll in range(l):
        if cartesian_basis:
          ao['exp_list'].append((aonames[count].lower().count('x'),
                                 aonames[count].lower().count('y'),
                                 aonames[count].lower().count('z')))
        else:
          m = aonames[count].lower().split('_')[-1]
          m = m.replace('+',' +').replace('-',' -').replace('s','s 0').split(' ') 
          p = 'yzx'.find(m[0][-1])
          if p != -1:
            m = p - 1
          else:
            m = int(m[-1])
          qc.ao_spherical.append([i,(lquant[ao['type']],m)])
        count += 1
  
  # Converting all information about molecular orbitals
  ele_num = numpy.sum(jData["molecule"]['atoms_Z']) - numpy.sum(jData['comp_details']['general']['core_electrons_per_atoms']) - jData['molecule']['charge']
  ue = (jData['molecule']['multiplicity']-1)
  
  # Check for natural orbitals and occupation numbers
  is_natorb = False
  #if hasattr(ccData,'nocoeffs'):
  #  if not hasattr(ccData,'nooccnos'):
  #    raise IOError('There are natural orbital coefficients (`nocoeffs`) in the cclib' + 
  #                  ' ccData, but no natural occupation numbers (`nooccnos`)!')
  #  is_natorb = True
  
  restricted = (len(jData['results']['wavefunction']['MO_energies']) == 1)
  if spin is not None:
    if spin != 'alpha' and spin != 'beta':
      raise IOError('`spin=%s` is not a valid option' % spin)
    elif restricted:
      raise IOError('The keyword `spin` is only supported for unrestricted calculations.')
    else:
      display('Converting only molecular orbitals of spin %s.' % spin)
  
  import scipy.sparse
  sym = {}
  shape = (jData['results']['wavefunction']['MO_number_kept'], jData['comp_details']['general']['basis_set_size'])
  pre_mocoeffs = jData['results']["wavefunction"]["MO_coefs"]
  if restricted:
    add = ['']
    orb_sym = [None]
    mocoeffs = [numpy.asarray(scipy.sparse.csr_matrix(tuple([numpy.asarray(d) for d in pre_mocoeffs[0]]), shape=shape).todense())]
  else:
    add = ['_a','_b']      
    orb_sym = ['alpha','beta']
    mocoeffs = [numpy.asarray(scipy.sparse.csr_matrix(tuple([numpy.asarray(d) for d in pre_mocoeffs[0]]), shape=shape).todense()),
                numpy.asarray(scipy.sparse.csr_matrix(tuple([numpy.asarray(d) for d in pre_mocoeffs[1]]), shape=shape).todense())]

  nmo = jData['results']['wavefunction']['MO_number'] if "nmo" in jData['results']['wavefunction'] else len(mocoeffs[0])
  for ii in range(nmo):
    for i,j in enumerate(add):
      a = '%s%s' % (jData['results']['wavefunction']['MO_sym'][i][ii],j)
      if a not in sym.keys(): sym[a] = 1
      else: sym[a] += 1
      #if is_natorb:
      #  occ_num = ccData.nooccnos[ii]
      if not restricted:
        occ_num = 1.0 if ii <= jData['results']['wavefunction']['homo_indexes'][i] else 0.0
      elif ele_num > ue:
        occ_num = 2.0
        ele_num -= 2.0
      elif ele_num > 0.0 and ele_num <= ue: 
        occ_num = 1.0
        ele_num -= 1.0
        ue -= 1.0
      else:
        occ_num = 0.0
      qc.mo_spec.append({'coeffs': mocoeffs[i][ii],
              'energy': jData['results']['wavefunction']['MO_energies'][i][ii],
              'occ_num': occ_num,
              'sym': '%d.%s' %(sym[a],a)
              })
      if orb_sym[i] is not None:
        qc.mo_spec[-1]['spin'] = orb_sym[i]
        if spin is not None and spin != orb_sym[i]:
          del qc.mo_spec[-1]
  
  # Use default order for atomic basis functions if aonames is not present
  if 'ao_names' not in jData['comp_details']['general']:
    display('The attribute `aonames` is not present in the parsed data.')
    display('Using the default order of basis functions.')
    
    # Check which basis functions have been used
    c_cart = sum([l_deg(l=ao['type'], cartesian_basis=True) for ao in qc.ao_spec])
    c_sph = sum([l_deg(l=ao['type'], cartesian_basis=False) for ao in qc.ao_spec])
    
    c = create_mo_coeff(qc.mo_spec,'').shape[-1]
    if c != c_cart and c == c_sph: # Spherical basis
      qc.ao_spherical = get_ao_spherical(qc.ao_spec,p=[0,1])
    elif c != c_cart:
      display('Warning: The basis set type does not match with pure spherical ' +
              'or pure Cartesian basis!') 
      display('Please specify qc.mo_spec["exp_list"] and/or qc.ao_spherical by your self.')
  
  # Are all MOs requested for the calculation? 
  if not all_mo:
    for i in range(len(qc.mo_spec))[::-1]:
      if qc.mo_spec[i]['occ_num'] < 0.0000001:
        del qc.mo_spec[i]

  return qc

## End convert_json

def get_ao_spherical(ao_spec,p=[1,0]):
  ao_spherical = []
  for i,ao in enumerate(ao_spec):
    ii = ao['type']
    l = lquant[ii]
    for m in (range(0,l+1) if l != 1 else p):
      ao_spherical.append([i,(l,m)])
      if m != 0:
        ao_spherical.append([i,(l,-m)])
    #for m in (range(1,l+1) if l != 1 else p):
      #if m != 0:
        #ao_spherical.append([i,(l,-m)])
  return ao_spherical

def mo_select(mo_spec, fid_mo_list, strict=False):
  '''Selects molecular orbitals from an external file or a list of molecular 
  orbital labels.

  **Parameters:**
   
    mo_spec :        
      See :ref:`Central Variables` for details.
    strict : bool, optional
      If True, orbkit will follow strictly the fid_mo_list, i.e., the order of 
      the molecular orbitals will be kept and multiple occurrences of items 
      will evoke multiple calculations of the respective molecular orbitals. 
    fid_mo_list : str, `'all_mo'`, or list
      | If fid_mo_list is a str, specifies the filename of the molecular orbitals list.
      | If fid_mo_list is 'all_mo', creates a list containing all molecular orbitals.
      | If fid_mo_list is a list, provides a list (or a list of lists) of molecular 
        orbital labels.

  **Supported Formats:**
  
    Integer List (Counting from **ONE**!)::
    
      1       2       3
      5       4
      homo    lumo+2:lumo+4
    
    List with Symmetry Labels::
    
      1.1     2.1     1.3
      1.1     4.1
      4.1     2.3     2.1
  
  **Returns:**
  
    Dictionary with following Members:
      :mo: - List of molecular orbital labels.
      :mo_ii: - List of molecular orbital indices.
      :mo_spec: - Selected elements of mo_spec. See :ref:`Central Variables` for details.
      :mo_in_file: - List of molecular orbital labels within the fid_mo_list file.
      :sym_select: - If True, symmetry labels have been used. 
  
  ..attention:
    
    For **unrestricted** calculations, orbkit adds `_a` (alpha) or `_b` (beta) to
    the symmetry labels, e.g., `1.1_a`. 
    If you have specified the option `spin=alpha` or `spin=beta`, only the 
    alpha or the beta orbitals are taken into account for the counting 
    within the Integer List.
  '''
  import re
  display('\nProcessing molecular orbital list...')
  
  mo_in_file = []
  selected_mo = []
  sym_select = False
  
  def assign_selected_mo(selected_mo,mo_spec,strict=False,  
                          what=lambda x,y: y[x]['sym']):
    selected_mo_spec = []
    selected_mo_ii = [] 
    for i in selected_mo:
      is_present = False
      for k in range(len(mo_spec)):
        if (what(k,mo_spec) == i):
          is_present = True
          if strict or (i not in selected_mo_ii):
            selected_mo_spec.append(mo_spec[k])
            selected_mo_ii.append(what(k,mo_spec))
      if not is_present:
        raise IOError('Cannot find %s in mo_spec' % i)
    selected_mo_ii = numpy.array(selected_mo_ii)
    return selected_mo_spec,selected_mo_ii
  
  def expr2int(expr):
    if isinstance(expr,int):
      return expr
    x = 0
    for i in re.findall(r'\d+|[+-]\d+',expr):
      x += int(i)
    return x
  
  def get_selection(selected_mo):
    mo_occup = numpy.array([i['occ_num'] for i in mo_spec])
    homo = (mo_occup>0.).nonzero()[0][-1]   + 1 # molden numbering
    lumo = (mo_occup>0.).nonzero()[0][-1]+1 + 1 # molden numbering
    mo_energy = numpy.array([i['energy'] for i in mo_spec])
    last_bound = sum(mo_energy<=0.0)            # molden numbering
    sel = []
    for i in selected_mo:
      i = i.lower().replace('homo',str(homo)).replace('lumo',str(lumo))
      i = i.replace('last_bound',str(last_bound))
      if ':' in i:
        k = [1,len(mo_spec)+1,1]
        i = i.split(':')
        for ik,j in enumerate(i):
          if j != '': k[ik] = j
        i = list(range(*[expr2int(j) for j in k]))
        sel.extend(i)
      else:
        sel.append(expr2int(i))
    return sel
  
  if isinstance(fid_mo_list,str) and fid_mo_list.lower() == 'all_mo':
    selected_mo = numpy.array(numpy.arange(len(mo_spec))+1, dtype=numpy.str)
    mo_in_file = [selected_mo]
    selected_mo_spec = mo_spec
    selected_mo_ii = numpy.array([i['sym'] for i in selected_mo_spec])
  else:
    if isinstance(fid_mo_list,str) and not os.path.exists(fid_mo_list):
      if ',' in fid_mo_list:
        fid_mo_list = fid_mo_list.split(',')
      else:
        fid_mo_list = [fid_mo_list]
    if isinstance(fid_mo_list, list):
      for i in fid_mo_list:
        if not isinstance(i, list):
          i = i.split(',') if isinstance(i,str) else [i]
        selected_mo.extend(list(map(str,i)))
        mo_in_file.append(list(map(str,i)))
    else:
      try:
        fid=open(fid_mo_list,'r')
        flines = fid.readlines()
        fid.close()
        for line in flines:
          integer = line.replace(',',' ').split()
          mo_in_file.append(integer)
          selected_mo.extend(integer)
      except:
        raise IOError('The selected mo-list (%(m)s) is not valid!' % 
                      {'m': fid_mo_list} + '\ne.g.\n\t1\t3\n\t2\t7\t9\n')
    
    # Print some information
    for i,j in enumerate(mo_in_file):
      display('\tLine %d: %s' % (i+1,', '.join(j)))
    
    # Check if the molecular orbitals are specified by symmetry 
    # (e.g. 1.1 in MOLPRO nomenclature) or 
    # by the number in the input file (e.g. 1)
    
    try: # Try to convert selections into integer
      for i in selected_mo:
        if isinstance(i,int):
          continue
        i = i.replace('homo','1').replace('lumo','2').replace('last_bound','3')
        for r in ['-','+',':']:
          i = i.replace(r,'')
        int(i)
    except ValueError:
      sym_select = True
      errors = []
      for i in range(len(selected_mo)):
        if not '.' in selected_mo[i]:
          errors.append(i)      
      if errors:
        err = [selected_mo[i] for i in errors]
        raise IOError('`%s` are no valid labels according '% ', '.join(err) +
                      'to the MOLPRO nomenclature, e.g., `5.1` or `5.A1`.' +
                      '\n\tHint: You cannot mix integer numbering and MOLPRO\'s ' +
                      'symmetry labels')
    
    if sym_select:
      what = lambda x,y: y[x]['sym']
      selected_mo_spec,selected_mo_ii = assign_selected_mo(selected_mo,
                                                           mo_spec,
                                                           strict=strict,
                                                           what=what)
    else:
      selected_mo = get_selection(selected_mo)
      
      if not strict:
        selected_mo = list(map(int, selected_mo))            
        selected_mo.sort()
      selected_mo = list(map(str, selected_mo))
      what = lambda x,y: str(x+1)
      selected_mo_spec,selected_mo_ii = assign_selected_mo(selected_mo,
                                                           mo_spec,
                                                           strict=strict,
                                                           what=what)
      selected_mo = selected_mo_ii
      for i in range(len(mo_in_file)):
        mo_in_file[i] = list(map(str, get_selection(mo_in_file[i])))
    
    # Print some information
    display('\nThe following orbitals will be considered...')
    for i,j in enumerate(mo_in_file):
      display('\tLine %d: %s' % (i+1,', '.join(j)))
  
  display('')
  return {'mo': selected_mo, 'mo_ii': selected_mo_ii,
          'mo_spec': selected_mo_spec, 
          'mo_in_file': mo_in_file, 'sym_select': sym_select}
