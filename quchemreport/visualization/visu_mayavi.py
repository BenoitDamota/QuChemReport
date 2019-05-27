#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-




import numpy as np
from sys import exit
import os

#from utils import A_to_a0, V_to_Kcal_mol
from quchemreport.utils import units

try:
    from enthought.mayavi import mlab
except ImportError:
    from sys import stderr
    stderr.write("import enthought.mayavi failed -- trying mayavi\n")
    try:
        from mayavi import mlab
    except ImportError:
        stderr.write("import mayavi failed\n")
        exit(1)

## No screen
mlab.options.offscreen = True


## Initialize visualization details common to all jobs
with open('%s/../utils/Atoms_properties.csv' % os.path.dirname(__file__), "r") as f:
    tab = [line.split() for line in f]

## Elevation angle
angle = 10

def _init_scene(j_data):
    u"""Initializes the MayaVi scene.

    ** Parameters **
      j_data : dict
    Data on the molecule, as deserialized from the scanlog format.

    ** Returns **
      figure : mayavi.core.scene.Scene
    The MayaVi scene with the atoms plotted.
      normal : numpy.ndarray
    The normal vector of the viewing plane (used mostly for the topological view).
    """
    figure = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))
    figure.scene.disable_render= True
    geom = np.array(j_data["results"]["geometry"]["elements_3D_coords_converged"]).reshape((-1,3))/units.A_to_a0
    if len(geom) > 1:
        ## Eliminate hydrogens
        #mod_geom = geom[np.array(map(int, [p[1] for p in qc.geo_info if p[2] != '1.0']))]
        mod_geom = geom
        ## Calculate best fitting plane via PCA
        eival, eivec = np.linalg.eig(np.cov((mod_geom - np.mean(mod_geom, axis=0)).T))
        sort = eival.argsort()
        eival, eivec = eival[sort], eivec[:,sort]

        normal = eivec[:,0]
        ## Grab point from best fitting plane (NOT the view) to use as focal point
        #point = np.mean(geom, axis=0)

        from math import sqrt, acos, atan2
        ## Calculate viewing distance r
        r = sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)
        ## Calculate azimuth a and elevation e
        ## Python and Numpy use radians, but MayaVi uses degrees
        a, e = np.rad2deg(atan2(normal[1], normal[0])), np.rad2deg(acos(normal[2]/r))
        mlab.view(azimuth=a, elevation=e+angle, figure=figure)

        ## DEBUG: show normal and view vectors
        #print a, e
        #print mlab.view()
        #mlab.quiver3d(point[0], point[1], point[2], normal[0], normal[1], normal[2])
        #mlab.quiver3d([0]*3, [0]*3, [0]*3, [1,0,0], [0,1,0], [0,0,1], color=(0,0,0))
        #mlab.quiver3d(*np.concatenate((np.zeros((3,2)), eivec[1:])), color=(0,0,1))

    else:
        normal = np.array([0,0,1])
        mlab.view(azimuth=0, elevation=0, figure=figure)

    conn = j_data["molecule"]["connectivity"]["atom_pairs"]
    atom_nums = j_data["molecule"]["atoms_Z"]

    ## Draw atoms and bonds
    for i, atom in enumerate(atom_nums):
        p, color = geom[i], tuple(float(x)/255.0 for x in tab[atom][3:6])

        ## Requires >=MayaVi-4.6.0
        mlab.points3d([p[0]], [p[1]], [p[2]],
                      figure=figure, mode='sphere', color=color, resolution=15, scale_factor=0.5)

    for pair in conn:
        att1 = tab[atom_nums[pair[0]]]
        p1, p2 = geom[pair[0]], geom[pair[1]]
        color = tuple(float(x)/255.0 for x in att1[3:6])

        mlab.quiver3d([p1[0]],         [p1[1]],         [p1[2]], 
                      [p2[0] - p1[0]], [p2[1] - p1[1]], [p2[2] - p1[2]],
                      figure=figure, mode='cylinder', color=color, resolution=15, scale_factor=0.5)

    #mlab.axes(figure=figure)

    return figure, normal

## Set the Iso contour value for mayavi from a percent
# Choose the % (between 0 to 100) of the positive values to show in picture.
#IsoContourPercent=30

def CalcCutOffP(data,IsoContourPercent=30,i=0):
	#data are the function values for each voxels. Transforma as a list, sort and select the positive values
    datar = data.ravel()
    datas = np.sort(datar[datar>0])#.ravel()) AttributeError: 'list' object has no attribute 'ravel'
    #cumulative sum of function values, normalize
    cumdata1DsortP=np.cumsum(datas)
    cumdata1DsortPnorm=cumdata1DsortP/cumdata1DsortP[-1]*100.
        #np.save("data%d.npy" % i, cumdata1DsortPnorm)
        #return the value of the voxel that is more intense than the IsoContourPercent. that should be the CutOff value.
    return datas[cumdata1DsortPnorm>=(100.-IsoContourPercent)][0]
def CalcCutOffm(data,IsoContourPercent=30,i=0):
	#data are the function values for each voxels. Transforma as a list, sort and select the negative values
    datar = data.ravel()
    datas = np.sort(datar[datar<0])#.ravel()) AttributeError: 'list' object has no attribute 'ravel'
    #cumulative sum of function values, normalize
    cumdata1DsortM=np.cumsum(datas)
    cumdata1DsortMnorm=cumdata1DsortM/cumdata1DsortM[-1]*100.
        #np.save("data%d.npy" % i, cumdata1DsortPnorm)
        #return the value of the voxel that is more intense than the IsoContourPercent. that should be the CutOff value.
    return datas[cumdata1DsortMnorm>=(100.-IsoContourPercent)][0]


## Visualize
def topo(j_data, file_name=None, size=(600,600)):
    u"""Creates the topological view of the molecule.

    ** Parameters **
      j_data : dict
    Data on the molecule, as deserialized from the scanlog format.
      file_name : str, optional
    Base name of the file in which to save the image.
      size : tuple(int, int), optional
    The size of the image to save.

    ** Returns **
      figure : mayavi.core.scene.Scene
    The MayaVi scene containing the visualization.
    """

    figure, normal = _init_scene(j_data)
    geom = np.array(j_data["results"]["geometry"]["elements_3D_coords_converged"]).reshape((-1,3))/units.A_to_a0

    ## Show labels and numbers ( = indices + 1 )
    for i, atom in enumerate(j_data["molecule"]["atoms_Z"]):
        P, label = geom[i], tab[atom][1]
        mlab.text3d(P[0] - normal[0], P[1] - normal[1], P[2] - normal[2], label + str(i + 1), color=(0,0,0), scale=0.5, figure=figure)

    if file_name is not None:
        mlab.savefig("{}-TOPOLOGY.png".format(file_name), figure=figure, size=size)

    return figure

def viz_MO(data, X, Y, Z, j_data, file_name=None, labels=None, size=(600,600)):
    u"""Visualizes the molecular orbitals of the molecule.

    ** Parameters **
      data : list(numpy.ndarray)
    List of series of voxels containing the scalar values of the molecular orbitals to plot.
      X, Y, Z
    Meshgrids as generated by numpy.mgrid, for positioning the voxels.
      j_data : dict
    Data on the molecule, as deserialized from the scanlog format.
      file_name : str, optional
    Base name of the files in which to save the images.
      labels : list(str), optional
    Labels to append to `file_name` for each series in `data`. If None, the position of the series is appended.
      size : tuple(int, int), optional
    The size of the image to save.

    ** Returns **
      figure : mayavi.core.scene.Scene
    The MayaVi scene containing the visualization.
    """

    figure, normal = _init_scene(j_data)
    for i, series in enumerate(data):

        MO_data = mlab.pipeline.scalar_field(X, Y, Z, series, figure=figure)
        Cutoffp = CalcCutOffP(series,i=i)
        Cutoffm = CalcCutOffm(series,i=i)
        print("positive Cutoff :", Cutoffp)
        print("negative Cutoff :", Cutoffm)
        MOp = mlab.pipeline.iso_surface(MO_data, figure=figure, contours=[ Cutoffp ], color=(0.4, 0, 0.235))
        MOn = mlab.pipeline.iso_surface(MO_data, figure=figure, contours=[ Cutoffm ], color=(0.95, 0.90, 0.93))

        if file_name is not None:
                        mlab.savefig("{}-MO-{}.png".format(file_name, labels[i] if labels is not None else i), figure=figure, size=size)

        MOp.remove()
        MOn.remove()

    return figure

def viz_EDD(data, X, Y, Z, j_data, file_name=None, labels=None, size=(600,600)):
    u"""Visualizes the electron density differences for the transitions of the molecule.

    ** Parameters **
      data : list(numpy.ndarray)
    Voxels containing the scalar values of the electron density differences to plot.
      X, Y, Z : numpy.ndarray
    Meshgrids as generated by numpy.mgrid, for positioning the voxels.
      j_data : dict
    Data on the molecule, as deserialized from the scanlog format.
      file_name : str, optional
    Base name of the files in which to save the images.
      labels : list(str), optional
    Labels to append to `file_name` for each series in `data`. If None, the position of the series is appended.
      size : tuple(int, int), optional
    The size of the image to save.

    ** Returns **
      figure : mayavi.core.scene.Scene
    The MayaVi scene containing the visualization.
    """

    figure, normal = _init_scene(j_data)
    for i, series in enumerate(data):
        D_data = mlab.pipeline.scalar_field(X, Y, Z, series, figure=figure)
        Cutoff = CalcCutOffP(series,i=i)
        print("Cutoff :", Cutoff)
        Dp = mlab.pipeline.iso_surface(D_data, figure=figure, contours=[ Cutoff ], color=(0.0, 0.5, 0.5))
        Dn = mlab.pipeline.iso_surface(D_data, figure=figure, contours=[-Cutoff ], color=(0.95, 0.95, 0.95))
        #Dn.actor.property.representation = 'wireframe'
        #Dn.actor.property.line_width = 0.5

        if file_name is not None:
            mlab.savefig("{}-EDD-{}.png".format(file_name, labels[i] if labels is not None else i), figure=figure, size=size)

        Dp.remove()
        Dn.remove()

    return figure

def viz_BARY(data, j_data, file_name=None, labels=None, size=(600,600)):
    u"""Visualizes the barycenters of the electron density difference (for visualizing dipole moments).

    ** Parameters **
      data : tuple(numpy.ndarray((3,N)), numpy.ndarray((3,N)))
    Pair of column-major matrices containing the coordinates of the positive and negative barycenters, in that order.
      j_data : dict
    Data on the molecule, as deserialized from the scanlog format.
      file_name : str, optional
    Base name of the files in which to save the images.
      labels : list(str), optional
    Labels to append to `file_name` for each datum in `data`. If None, the position of the datum is appended.
      size : tuple(int, int), optional
    The size of the image to save.

    ** Returns **
      figure : mayavi.core.scene.Scene
    The MayaVi scene containing the visualization.
    """

    figure, normal = _init_scene(j_data)

    for i, D in enumerate(data):
        #Pp = mlab.points3d(D[0][0], D[0][1], D[0][2], figure=figure, mode='axes', scale_factor=0.3, color=(0.0, 0.5, 0.5))
        #Pm = mlab.points3d(D[1][0], D[1][1], D[1][2], figure=figure, mode='axes', scale_factor=0.3, color=(0.95, 0.95, 0.95))

        ## Chemistry convention (from negative to positive)
        Mu = mlab.quiver3d(D[1][0], D[1][1], D[1][2], D[0][0] - D[1][0], D[0][1] - D[1][1], D[0][2] - D[1][2], figure=figure, mode='arrow', scale_factor=1.0, color=(0.0, 0.5, 0.5))

        if file_name is not None:
            mlab.savefig("{}-BARY-{}.png".format(file_name, labels[i] if labels is not None else i), figure=figure, size=size)

        #Pp.remove()
        #Pm.remove()

        Mu.remove()

    return figure

def viz_Potential(r_data, V_data, X, Y, Z, j_data, file_name=None, size=(600,600)):
    u"""Visualizes the electrostatic potential difference of the molecule.

    ** Parameters **
      r_data, V_data : numpy.ndarray
    Voxels of the electron density and the potential difference of the molecule, respectively.
      X, Y, Z : numpy.ndarray
    Meshgrids as generated by numpy.mgrid, for positioning the voxels.
      j_data : dict
    Data on the molecule, as deserialized from the scanlog format.
      file_name : str, optional
    Base name of the file in which to save the image.
      size : tuple(int, int), optional
    The size of the image to save.

    ** Returns **
      figure : mayavi.core.scene.Scene
    The MayaVi scene containing the visualization.
    """

    figure, normal = _init_scene(j_data)

    src = mlab.pipeline.scalar_field(X, Y, Z, r_data, figure=figure)

    ## Add potential as additional array
    src.image_data.point_data.add_array(V_data.T.ravel()/units.V_to_Kcal_mol)

    ## Name it
    src.image_data.point_data.get_array(1).name = "potential"

    ## Update object
    src.update()

    ## Select scalar attribute
    srcp = mlab.pipeline.set_active_attribute(src, figure=figure, point_scalars="scalar")

    ## Plot it
    cont = mlab.pipeline.contour(srcp, figure=figure)
    cont.filter.contours=[0.001]

    ## Select potential
    cont_V = mlab.pipeline.set_active_attribute(cont, figure=figure, point_scalars="potential")
    #contp = mlab.pipeline.threshold(cont_V, figure=figure, up=)
    #contn = mlab.pipeline.threshold(cont_V, figure=figure, low=V_data.min()*0.95)

    ## And finally plot that
    mlab.pipeline.surface(cont_V, figure=figure, opacity=0.7)

    ## Continue with this until problems with potential calculation are fixed
    #V_data[np.isinf(V_data)] = np.nan

    #src = mlab.pipeline.scalar_field(X, Y, Z, V_data, figure=figure)
    #srcp = mlab.pipeline.iso_surface(src, figure=figure, contours=[ 0.4], color=(0.0, 0.5, 0.5))
    #srcn = mlab.pipeline.iso_surface(src, figure=figure, contours=[-0.05], color=(0.95, 0.95, 0.95))

    mlab.colorbar(title='Kcal/mol', orientation='vertical', nb_labels=3)

    #mlab.show()

    if file_name is not None:
        mlab.savefig("{}-Potential.png".format(file_name), figure=figure, size=size)

    return figure

def viz_Fukui(data, X, Y, Z, j_data, file_name=None, labels=None, size=(600,600)):
        u"""Visualizes the fukui density differences for the molecule.

        ** Parameters **
          data : list(numpy.ndarray)
        Voxels containing the scalar values of the electron density difference to plot.
          X, Y, Z : numpy.ndarray
        Meshgrids as generated by numpy.mgrid, for positioning the voxels.
          j_data : dict
        Data on the optimized state of the molecule, as deserialized from the scanlog format.
          file_name : str, optional
        Base name of the files in which to save the images.
          label : text
        Label to append to `file_name` for each series in `data`. "plus" for sp_plus - opt; "minus" for opt - sp_minus.
          size : tuple(int, int), optional
        The size of the image to save.

        ** Returns **
          figure : mayavi.core.scene.Scene
        The MayaVi scene containing the visualization.
        """

        figure, normal = _init_scene(j_data)
        F_data = mlab.pipeline.scalar_field(X, Y, Z, data, figure=figure)
        Cutoff = CalcCutOffP(data)
        print("Cutoff :", Cutoff)
        Fp = mlab.pipeline.iso_surface(F_data, figure=figure, contours=[ Cutoff ], color=(0.0, 0.5, 0.5))
        Fn = mlab.pipeline.iso_surface(F_data, figure=figure, contours=[-Cutoff ], color=(0.95, 0.95, 0.95))
        #Fn.actor.property.representation = 'wireframe'
        #Fn.actor.property.line_width = 0.5

        if file_name is not None:
            mlab.savefig("{}-fukui-{}.png".format(file_name, labels), figure=figure, size=size)
            Fp.remove()
            Fn.remove()

        return figure

def viz_Fdual(data, X, Y, Z, j_data, file_name=None, size=(600,600)):
        u"""Visualizes the fukui density differences for the molecule.

        ** Parameters **
          data : list(numpy.ndarray)
        Voxels containing the scalar values of the electron density difference to plot.
          X, Y, Z : numpy.ndarray
        Meshgrids as generated by numpy.mgrid, for positioning the voxels.
          j_data : dict
        Data on the optimized state of the molecule, as deserialized from the scanlog format.
          file_name : str, optional
        Base name of the files in which to save the images.
          size : tuple(int, int), optional
        The size of the image to save.

        ** Returns **
          figure : mayavi.core.scene.Scene
        The MayaVi scene containing the visualization.
        """

        figure, normal = _init_scene(j_data)
        F_data = mlab.pipeline.scalar_field(X, Y, Z, data, figure=figure)
        #D_data = mlab.pipeline.scalar_field(X, Y, Z, series, figure=figure)
        Cutoff = CalcCutOffP(data)
        print("Cutoff :", Cutoff)
        Fp = mlab.pipeline.iso_surface(F_data, figure=figure, contours=[ Cutoff ], color=(0.0, 0.5, 0.5))
        Fn = mlab.pipeline.iso_surface(F_data, figure=figure, contours=[-Cutoff ], color=(0.95, 0.95, 0.95))
        #Dn.actor.property.representation = 'wireframe'
        #Dn.actor.property.line_width = 0.5

        if file_name is not None:
            mlab.savefig("{}-Fdual.png".format(file_name), figure=figure, size=size)
            Fp.remove()
            Fn.remove()

        return figure


        pass
