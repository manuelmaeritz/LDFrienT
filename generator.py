#!/usr/bin/env python3

"""This module supplies several functions for generating objects of the class
`LdftModel` (or its inherited classes), picard iterate them and save them.
Those functions can either be loaded by another module or alternatively one
can directly call this script in the shell with some parameters determining
what it is supposed to do.

For calling the script in the shell several modes are supplied. The `mode` is
always the first parameter handed over and it defines what other parameters are
required. Depending on the chosen mode the script will do different things.
Following modes are supported:

`mode`='single':
    In the 'single'-mode, the script will create a single canonical (fixed
    average density) system, picard-iterate it and save it. It takes the
    following further parameters:
    `dftModel` :
        The DFT-model. Either '2d-highl', '3d-highl', '2d-mf' or '3d-mf'.
    `size` :
        The size of the (quadratic) system. E.g. '64' for a two-dimensional
        system of size 64x64 or a three dimensional system of size 64x64x64.
    `epsi` :
        The attraction strength of the systems particles (times inverse
        temperature).
    `dens` :
        The average density of the system.
    `init` :
        The initial state (density profile), the iteration will start with.
        This might either be 'hom' for an homogeneous system, 'sph' for a
        spherical nucleus, 'sl' for a planar nucleus or 'cyl' for a cylindrical
        nucleus (in tow dimensions 'sl' and 'cyl' is the same). One, however,
        could also pass the filename of an already existing system (which is of
        the same dft-model, size and attraction strength). Then, the latest
        state of this system is taken as the initial configuration. This
        filename has to start with the letters 'dens'.
    `alpha` :
        The iteration parameter, determining how fast the iteration goes. It is
        a value in the range (0, 1).
    `it_steps` :
        Number of iteration steps which shall be done.
    `chp` :
        A parameter determining in which frequency snapshots of the system are
        added to the history of the system. This might be an integer (e.g.
        '10000' for snapshots in regular iteration intervals) or 'dec#', where
        # stands for any integer. Then the intervals are close in the beginning
        and become more distant later. See the documentation of the class
        `LdftModel` for more detail.
`mode`='continue' :
    In the 'continue'-mode, the script loads an already existing system, and
    continues its iteration and saves the result. Note, that the loaded system
    generally will not be overwritten by the further iterated system. This is
    just the case if the save parent directory matches the source parent
    directory (which is not recommended). (See below for more details regarding
    the save and source path.) In difference to the 'single'-mode with the
    parameter `init` referring to an already existing system, here not only the
    density profile is taken over as initial state, but the hole system
    including its `history`-parameter, average density and so on. The further
    parameters in this mode are:
    `dftModel` :
        The DFT-model of the system which one wants to load. Either '2d-highl',
        '3d-highl', '2d-mf' or '3d-mf'.
    `size` :
        The size of the (quadratic) system, which one wants to load. E.g. '64'
        for a two-dimensional system of size 64x64 or a three dimensional system
        of size 64x64x64.
    `epsi` :
        The attraction strength of the particles (times inverse temperature) in
        the system which one wants to load.
    `alpha` :
        The iteration parameter, determining how fast the iteration goes. It is
        a value in the range (0, 1).
    `it_steps` :
        Number of iteration steps which shall be done.
    `chp` :
        A parameter determining in which frequency snapshots of the system are
        added to the history of the system. This might be an integer (e.g.
        '10000' for snapshots in regular iteration intervals) or 'dec#', where
        # stands for any integer. Then the intervals are close in the beginning
        and become more distant later. See the documentation of the class
        `LdftModel` for more detail.
    `syst` :
        The filename of the system one wants to load.
`mode`='series' :
    In the 'series'-mode, a series of canonical (fixed average density) systems
    is generated in a range of density [`startDens`, `maxDens`) . They are
    created, iterated and saved consecutively. It takes the following further
    parameters:
    `dftModel` :
        The DFT-model of the series. Either '2d-highl', '3d-highl', '2d-mf' or
        '3d-mf'.
    `size` :
        The size of the (quadratic) systems. E.g. '64' for two-dimensional
        system of size 64x64 or three dimensional system of size 64x64x64.
    `epsi` :
        The attraction strength of the systems particles (times inverse
        temperature).
    `startDens` :
        It defines one limit of the range of densities. This might also be the
        upper limit. It is the average density of the first system which is
        generated.
    `maxDens` :
        It defines the other limit of the range of densities (beside
        `startDens`). This value is not assumed by any system (it is an open
        limit).
    `stpWidth` :
        This parameter determines the step width of the series of average
        densities.
    `init` :
        The initial state (density profile), the iteration will start with.
        This might be a certain shape of nuclei ('hom', 'sph', 'sl', 'cyl') or
        the profile of an already existing system which shall be taken over.
        See the description of the `init`-parameter of the 'single'-mode for
        more detail. Weather this parameter holds only for the first system
        which is generated or the hole series is determined by the
        `cons`-parameter below.
    `alpha` :
        The iteration parameter, determining how fast the iteration goes. It is
        a value in the range (0, 1) and holds for the hole series.
    `it_steps` :
        Number of iteration steps which shall be done. It holds for the hole
        series.
    `chp` :
        A parameter determining in which frequency snapshots of the system are
        added to the history of the systems. See the description of the
        parameter `chp` of the 'single'-mode for more detail.
    `cons` :
        This parameter determines weather the `init`-parameter holds for the
        hole series (`cons`='0') or just the first system of the series
        (`cons`=1). In the second case the final density profile of the
        previous system is taken as the initial profile for the current system.
`mode`='searchTrans' :
    In the 'searchTrans'-mode the script will search for a phase transition
    within a certain range of densities. If there is a phase transition, it
    will echo a density interval in which the transition occurs. All systems
    which are created in this process are saved. If there is no transition it
    will also inform you.
    The possible interval in which the transition occurs is iteratively more
    and more restricted. This is done by picard iterating two systems in the
    middle of the current interval. Both systems have different initial density
    profiles, such that the final state results in different phases for both
    cases. The free energy of the two systems then is compared and depending on
    which system has the lower free energy the transition density is either
    left or right to the current density of the systems. In this way the
    Interval is halved in each step. The following further arguments are
    required: 
    `dftModel` :
        The DFT-model. Either '2d-highl', '3d-highl', '2d-mf' or
        '3d-mf'.
    `size` :
        The size of the (quadratic) system. E.g. '64' for two-dimensional
        system of size 64x64 or three-dimensional system of size 64x64x64.
    `epsi` :
        The attraction strength of the systems particles (times inverse
        temperature).
    `minDens` :
        It defines the lower bound of the initial interval the transition is
        suspected in.
    `maxDens` :
        It defines the upper bound of the initial interval the transition is
        suspected in.
    `alpha` :
        The iteration parameter, determining how fast the iteration goes. It is
        a value in the range (0, 1) and holds for the hole series.
    `it_steps` :
        Number of iteration steps which shall be done. It holds for all systems
        created in the searching process.
    `chp` :
        A parameter determining in which frequency snapshots of the system are
        added to the history of the systems. See the description of the
        parameter `chp` of the 'single'-mode for more detail. It holds for all
        systems created in the searching process.
    `accuracy` :
        This parameter determines the accuracy of the output interval. The
        smaller it is, the longer the searching process will take and the more
        precise (smaller) the output interval is.
    `initLeft` :
        The initial profile for the less dense phase. This might be a certain
        shape of nuclei ('hom', 'sph', 'sl', 'cyl') or the profile of an
        already existing system which shall be taken over. See the description
        of the `init`-parameter of the 'single'-mode for more detail.
    `initRight` :
        Same as for `initLeft` but for the more dense phase.

The generated systems are saved in a path:
**saveParentDirecotry/model/size/temperature**. Here, the folder names 'model',
'size' and 'temperature' are determined by the parameters of the corresponding
system. The save parent directory may by chosen by yourself in the function
``savePathDecorator''. The filename is of the form 'dens=0,3.pkl', where '0,3'
is the corresponding density of the system.
The modes, which require an already existing system, look for them in the
directory **sourceParentDirecotry/model/size/temperature**. Again the folder
names 'model', 'size', 'temperature' are determined by the parameters of the
system. The source parent directory may be chosen by yourself in the function
``sourcePathDecorator''.
"""


from ldft_classes_v2.lg_2d_highl import LG2dAOHighl
from ldft_classes_v2.lg_2d_mf import LG2dMf
from ldft_classes_v2.lg_3d_highl import LG3dAOHighl
from ldft_classes_v2.lg_3d_mf import LG3dMf
import numpy as np
import os
from functools import reduce
import sys

def create_path_decorator(motherPath):
    """Creates a decorator for the function `create_daughter_path`. The
    decorator appends the parameter ``motherPath`` to the ``daughterPath``
    created by `create_daughter_path`. Furthermore, if ``motherPath`` is
    a relative path, it turns the resulting path (i.e.
    motherPath/daughterPath) to an absolute path.

    Parameters
    ----------
    motherPath : `String`
        relative or absolute path, where the results generated by this
        script should be stored.

    Returns
    -------
    decorator : `function`
        Decorator for the `create_daughter_path`-function, which adds
        `motherPath` to the path created by `create_daugther_path`.
    """
    def pathDecorator(create_daughter_path):
        def path_wrapper(dftModel, size, epsi):
            path = os.path.join(motherPath,\
                    create_daughter_path(dftModel, size, epsi))
            path = os.path.abspath(path)
            return path
        return path_wrapper
    return pathDecorator


def create_daughter_path(Model, size, epsi):
    """Creates a file structure in which the data created by this script
    should be stored or data which needs to be loaded suppose to be
    found in. The structure is as following: 'dftModel/size/epsi' (e.g.
    '2d-highl/size=64/epsi=1,2').

    Parameters
    ----------
    Model : `Class`
        The model under which the data is generated.
    size : `int`
        The system size.
    epsi : `float`
        The attraction strength of the system which is generated (times
        the inverse temperature)

    Returns
    -------
    daughterPath : `String`
        The file structure under which the generated data should be
        stored or the data which needs to be loaded is assumed to be
        found in.
    """
    if Model == LG2dAOHighl: modelFolder = '2d-highl'
    elif Model == LG3dAOHighl: modelFolder = '3d-highl'
    elif Model == LG2dMf: modelFolder = '2d-mf'
    elif Model == LG3dMf: modelFolder = '3d-mf'
    sizeFolder = 'size='+str(size)
    epsiFolder = 'epsi='+str(np.round(epsi, 3)).replace('.', ',')
    path = os.path.join(modelFolder, sizeFolder, epsiFolder)
    return path

sourcePathDecorator = create_path_decorator('../samples')
"""Decorator for the function `create_daughter_path`. It adds the path
'../samples' to the path returned by `create_daughter_path`. This is the mother
path where the script can look for data.

Parameters
----------
create_daughter_path : `function`
    A function which creates the file structer under which the data is
    supposed to be found.

Returns
-------
path_wrapper : `function`
    A new function which returns the absolute path under which the
    script can look for specific data.
"""

savePathDecorator = create_path_decorator('./generatorOutput')
"""Decorator for the function `create_daughter_path`. It adds the path
'./generatorOutput' to the path returned by `create_daughter_path`. This is
the mother path where the script can dump data.

Parameters
----------
create_daughter_path : `function`
    A function which creates the file structer under which the generated
    data of this script should be stored.

Returns
-------
path_wrapper : `function`
    A new function which returns the absolute path under which the
    generated data of this script should be stored.
"""

create_source_path = sourcePathDecorator(create_daughter_path)
"""Creates a path under which the script can look for data.

Parameters
----------
dftModel : `String`
    The model under which the desired data was created.
size : `int`
    The system size of the desired data.
epsi : `float`
    The attraction strength of the desired system (times
    the inverse temperature).

Returns
-------
sourcePath : `String`
    The path under which the script can look for data.
"""

create_save_path = savePathDecorator(create_daughter_path)
"""Creates a path under which the script can dump the generated data.

Parameters
----------
dftModel : `String`
    The model under which the data was generated.
size : `int`
    The system size.
epsi : `float`
    The attraction strength of the System (times the inverse
    temperature).

Returns
-------
sourcePath : `String`
    The path under which the script can look for data.
"""

def create_sys_name(dens, spec=None):
    """Creates the name of the file under which a system generated by
    this script should be stored. E.g. 'dens=0,32(sl).pkl'

    Parameters
    ----------
    dens : `float`
        Average density of the system.
    spec : `String`; optional: default = '-'
        Any specification. E.g. 'sl' or 'sph' if the initial profile was
        slab or a point nucleus. Or 'ser' if the system was produced in
        a serial of average densities and the initial profile was the
        result of the last generated density.

    Returns
    -------
    sysName : `string`
        The name of the system under which it should be dumped.
    """
    dens = np.round(dens, 4)
    if spec != None: spec = '('+spec+')'
    else: spec = ''
    sysName = 'dens='+str(dens).replace('.', ',')+spec+'.pkl'
    return sysName

def create_sys(Model, size, epsi, dens, init, bd_cond='periodic'):
    """Creates a instance of the LdftModel `Model` under the parameters `size`,
    `epsi`, `dens`. The `init`-parameter determines the initial density profile
    of the created system.
    
    Parameters
    ----------
    Model : `class`
        The model, the system should be an instance of. Either ``LGAO2dHighl``,
        ``LGAO3dHighl``, ``LGAO2dMf`` or ``LGAO3dMf``.
    size : `Tuple` of `int`
        Size of the system (for each dimension).
    epsi : `float`
        Attraction strength (times inverse temperature).
    dens : `float`
        Average density
    init : `string` or `Model`
        Specifies the initial density profile. It is either a system providing
        its current density profile or a String. The following strings are
        supported: 'hom', 'sph', 'sl', 'cyl'. Note, that 'sl' and 'cyl' is the
        same in two dimensions. The density profile then either is homogeneous
        (`init`='hom'), a point-like nucleus (`init`='sph'), a slab
        (`init`='sl') or a cylinder (`init`='cyl'). The cylinder point in the
        0'th axis. The normal of the slab points in the last axis (2 in three
        dim. and 1 in two dim.). For bd_cond '11_if', '110_if' or '111_if'
        therefore the last entry of the `size`-tuple is supposed to be double
        of the others.
    bd_cond : `string`
        The boundary condition of the system. It supports the following values:
        'periodic', '11_if', '110_if', '111_if'. For bd_cond '11_if', '110_if'
        or '111_if' therefore the last entry of the `size`-tuple is supposed to
        be double of the others.

    Returns
    -------
    syst : `Model`
        System of the type `Model`, with initial density profile.
    """
    #Create the system in case of no profile of another system shall be
    #taken over as initial profile
    if type(init)==str:
        nucLen = max(size[0]//20, 1)
        #Determine the shape of the nucleus
        if init == 'hom':
            nuc = tuple([0 for s in size])
        elif init == 'sph':
            nuc = tuple([nucLen for s in size])
        elif init == 'cyl':
            nuc = [nucLen for s in size]
            nuc[0] = size[0]
            nuc = tuple(nuc)
        elif init == 'sl':
            nuc = [s for s in size]
            nuc[-1] = nucLen
            nuc = tuple(nuc)
        #Create the initial density profile
        if Model == LG2dAOHighl:
            r_pc = Model._cal_bulk_r_pc(dens, epsi)
            syst = Model(size, epsi=epsi, dens_c=dens,
                    bound_cond=bd_cond)
            syst.create_init_profile(dens=[dens, r_pc, r_pc],
                        shape=[nuc, (0, 0), (0, 0)])
        if Model == LG3dAOHighl:
            r_pc = Model._cal_bulk_r_pc(dens, epsi)
            syst = Model(size, epsi=epsi, dens_c=dens,
                    bound_cond=bd_cond)
            syst.create_init_profile(dens=[dens, r_pc, r_pc, r_pc],
                        shape=[nuc, (0, 0, 0), (0, 0, 0), (0, 0, 0)])
        if Model == LG2dMf or Model == LG3dMf:
            syst = Model(size, epsi=epsi, dens=dens, bound_cond=bd_cond)
            syst.create_init_profile(dens=[dens], shape=[nuc])
    #Create the system in case when a density profile of another system
    #shall be taken over as initial profile
    if type(init)==Model:
        initProfile = init.r
        if Model == LG2dAOHighl or Model == LG3dAOHighl:
            syst = Model(size, epsi=epsi, dens_c=dens, r=initProfile,
                    bound_cond=bd_cond)
        elif Model == LG2dMf or Model == LG3dMf:
            syst = Model(size, epsi=epsi, dens=dens, r=initProfile,
                    bound_cond=bd_cond)
    return syst

def iterate_and_save(sys, alpha, it_steps, chp, saveSpec=None,
        min_err=10**(-20)):
    """Picard iterates a system and saves the system on the hard drive disk.
    The save path is determined by the savePathDecorator (around line 290).

    Parameters
    ----------
    sys : `LdftModel`
        The system you want to find the equilibrium profile of.
    alpha : `Float`
        The iteration parameter determining how fast the iteration is. Value
        between (0, 1).
    it_steps : `integer`
        The number of iteration steps which shall be done.
    chp : `integer` or `string`
        A parameter determining in which frequency snapshots of the system are
        added to the history of the system. This might be an integer (e.g.
        '10000' for snapshots in regular iteration intervals) or 'dec#', where
        # stands for any integer. Then the intervals are close in the beginning
        and become more distant later. See the documentation of the class
        `LdftModel` for more detail.
    saveSpec : `string`; default value: None
        A string which is added to the default filename under which the
        resulting system is saved.
    min_err : `Float`; default value: 10**(-20)
        A minimal picard error. When the picard error fall below this value,
        the iteration is aborted early.
    """
    #Create save path
    size = sys.size
    epsi = sys.epsi
    dens = sys._dens[0]
    if reduce(lambda a, b: (a+b)/2, size)==size[0]:
        savePath = create_save_path(type(sys), size[0], epsi)
    else:
        savePath = create_save_path(type(sys), size, epsi)
    #iterate
    sys.make_picard_iteration(alpha, it_steps, chp, min_err=min_err)
    sys.save_syst(savePath, create_sys_name(dens, saveSpec))

def generate_series(Model, size, epsi, densRange, stepWidth, init,
        alpha, it_step, chp, saveSpec=None, consecutive=False):
    """
    A series of canonical (fixed average density) systems is generated in a
    rage of densities. They are created, iterated and saved consecutively.

    Parameters
    ----------

    Model : `LdftModel`
        The model the resulting systems of the series should be an instance of.
    size : `Tuple`
        The size of the systems in each dimension.
    epsi : `Float`
        The attraction strength of the systems particles (times inverse
        temperature).
    densRange : `Tuple` or `List`
        It defines of the density interval in which the series should be
        produced. The first value defines the starting density and the second
        value defines the end of the series (which will not be assumed). The
        starting value might also be the upper limit of the interval. In this
        case, however, the parameter `stepWidth` has to be chosen negatively.
    stepWidth : `float`
        This parameter determines the step width in which the intervall
        `densRange` is walked through.
    init : `string` or `LdftModel`
        The initial state (density profile), the iteration will start with.
        This might either be 'hom' for an homogeneous system, 'sph' for a
        spherical nucleus, 'sl' for a planar nucleus or 'cyl' for a cylindrical
        nucleus (in tow dimensions 'sl' and 'cyl' is the same). One, however,
        could also pass another, already existing system. Then, the latest
        state of this system is taken as the initial configuration.
    alpha : `float`
        The iteration parameter, determining how fast the iteration goes. It is
        a value in the range (0, 1) and holds for the hole series.
    it_steps : `int`
        Number of iteration steps which shall be done. It holds for the hole
        series.
    chp : `string` or `int`
        A parameter determining in which frequency snapshots of the system are
        added to the history of the system. This might be an integer (e.g.
        '10000' for snapshots in regular iteration intervals) or 'dec#', where
        # stands for any integer. Then the intervals are close in the beginning
        and become more distant later. See the documentation of the class
        `LdftModel` for more detail.
    saveSpec : `string`; default=`None`
        A sting which is appended to the filenames of the saved systems.
    cons : `bool`; defautl=`False`
        This parameter determines weather the `init`-parameter holds for the
        hole series (`cons`='0') or just the first system of the series
        (`cons`=1). In the second case the final density profile of the
        previous system is taken as the initial profile for the current system.
    """
    dens = np.arange(densRange[0], densRange[1], stepWidth)
    for rho in dens:
        sys = create_sys(Model, size, epsi, rho, init)
        iterate_and_save(sys, alpha, it_step, chp, saveSpec)
        if consecutive: init = sys

def search_trans(Model, size, epsi, densRange, alpha, it_step, chp,
        acuracy, initLeft, initRight, specLeft='left', specRight='right'):
    """
    The function searches for a phase transition within a certain range of
    densities. If there is a phase transition, it will echo a density interval
    in which the transition occurs. All systems which are created in this
    process are saved. If there is no transition it will also inform you.
    The possible interval in which the transition occurs is iteratively more
    and more restricted. This is done by picard iterating two systems in the
    middle of the current interval. Both systems have different initial density
    profiles, such that the final state results in different phases for both
    cases. The free energy of the two systems then is compared and depending on
    which system has the lower free energy the transition density is either
    left or right to the current density of the systems. In this way the
    Interval is halved in each step.

    Parameters
    ----------
    Model : `LdftModel`
        The DFT-model, the phase transition of which shall be searched.
    size : `Tuple`
        The system size, in which the transition shall be searched. This tuple
        has an entry for each dimension.
    epsi : `float`
        The attraction strength of the systems particles (times inverse
        temperature).
    densRange : `Tuple`
        The density interval in which the transition is suspected.
    alpha : `float`
        The iteration parameter, determining how fast the iteration goes. It is
        a value in the range (0, 1) and holds for the hole series.
    it_steps : `int`
        Number of iteration steps which shall be done. It holds for all systems
        created in the searching process.
    chp : `int`, or `string`
        A parameter determining in which frequency snapshots of the system are
        added to the history of the system. This might be an integer (e.g.
        '10000' for snapshots in regular iteration intervals). Alternatively,
        one can hand over a string of the form 'dec#', where # stands for any
        integer. Then the intervals are close in the beginning and become more
        distant later. See the documentation of the class `LdftModel` for more
        detail.
    accuracy : `float`
        This parameter determines the accuracy of the output interval. The
        smaller it is, the longer the searching process will take and the more
        precise (smaller) the output interval is.
    initLeft : `string` or `LdftModel`
        The initial profile for the less dense phase. This might be a certain
        shape of nuclei ('hom', 'sph', 'sl', 'cyl') or the profile of an
        already existing system which shall be taken over. See the description
        of the `init`-parameter of the 'single'-mode for more detail.
    initRight : `string` or `LdftModel`
        Same as for `initLeft` but for the more dense phase.
    specLeft : `string`
        A sting appended to the filenames of the saved files, where the left
        phase has won the competition.
    specRight : `string`
        A sting appended to the filenames of the saved files, where the right
        phase has won the competition.
    """
    if reduce(lambda a, b: (a+b)/2, size)==size[0]:
        print(Model)
        savePath = create_save_path(Model, size[0], epsi)
    else:
        savePath = create_save_path(Model, size, epsi)
    if densRange[1]-densRange[0] > acuracy:
        dens = (densRange[1]+densRange[0])/2
        print('current density: '+str(dens))
        left = create_sys(Model, size, epsi, dens, initLeft)
        left.make_picard_iteration(alpha, it_step, chp, min_err=10**(-20))
        right = create_sys(Model, size, epsi, dens, initRight)
        right.make_picard_iteration(alpha, it_step, chp, min_err=10**(-20))
        leftF = left.cal_semi_Om()
        rightF = right.cal_semi_Om()
        print('leftF='+str(leftF))
        print('rightF='+str(rightF))
        if abs(leftF-rightF)<10**(-3):
            print('no Transition here')
            left.save_syst(savePath, create_sys_name(dens,
                spec='noTrans'))
        elif leftF<rightF:
            print('leftF<rightF')
            left.save_syst(savePath, create_sys_name(dens,
                spec=specLeft))
            densRange = [dens, densRange[1]]
            search_trans(Model, size, epsi, densRange, alpha, it_step,
                    chp, acuracy, left, initRight, specLeft, specRight)
        else:
            print('rightF<leftF')
            right.save_syst(savePath, create_sys_name(dens,
                spec=specRight))
            densRange = [densRange[0], dens]
            search_trans(Model, size, epsi, densRange, alpha, it_step,
                    chp, acuracy, initLeft, right, specLeft, specRight)
    else:
        print('Finish!!')
        print('The transition is somewehre between '\
            +str(densRange[0])+' and '+str(densRange[1]))

if __name__=='__main__':
    """Here stands what happens when the script is called by the command line.
    The description is at the beginning of the file.
    """
    #Evaluate the parameters which were handed over
    mode = str(sys.argv[1])
    
    if mode == 'continue':
        dftModel = str(sys.argv[2])
        size = int(sys.argv[3])
        epsi = float(sys.argv[4])
        alpha = float(sys.argv[5])
        it_steps = int(sys.argv[6])
        chp = str(sys.argv[7])
        syst = str(sys.argv[8])
    
        if dftModel == '2d-highl':
            Model = LG2dAOHighl
        elif dftModel == '2d-mf':
            Model = LG2dMf
        elif dftModel == '3d-highl':
            Model = LG3dAOHighl
        elif dftModel == '3d-mf':
            Model = LG3dMf

        sourcePath = create_source_path(Model, size, epsi)
        syst = Model.load_syst(sourcePath, syst)
        iterate_and_save(syst, alpha, it_steps, chp)
    
    if mode == 'single':
        dftModel = str(sys.argv[2])
        size = int(sys.argv[3])
        epsi = float(sys.argv[4])
        dens = float(sys.argv[5])
        init = str(sys.argv[6])
        alpha = float(sys.argv[7])
        it_steps = int(sys.argv[8])
        chp = str(sys.argv[9])
    
        if dftModel == '2d-highl':
            Model = LG2dAOHighl
            size = (size, size)
        elif dftModel == '2d-mf':
            Model = LG2dMf
            size = (size, size)
        elif dftModel == '3d-highl':
            Model = LG3dAOHighl
            size = (size, size, size)
        elif dftModel == '3d-mf':
            Model = LG3dMf
            size = (size, size, size)

        if init[0:4] == 'dens':
            sourcePath = create_source_path(Model, size[0], epsi)
            init = Model.load_syst(sourcePath, init)
            
        syst = create_sys(Model, size, epsi, dens, init)
        iterate_and_save(syst, alpha, it_steps, chp)

    if mode == 'series':
        dftModel = str(sys.argv[2])
        size = int(sys.argv[3])
        epsi = float(sys.argv[4])
        startDens = float(sys.argv[5])
        maxDens = float(sys.argv[6])
        stpWidth = float(sys.argv[7])
        init = str(sys.argv[8])
        alpha = float(sys.argv[9])
        it_steps = int(sys.argv[10])
        chp = str(sys.argv[11])
        cons = bool(int(sys.argv[12]))
    
        if dftModel == '2d-highl':
            Model = LG2dAOHighl
            size = (size, size)
        elif dftModel == '2d-mf':
            Model = LG2dMf
            size = (size, size)
        elif dftModel == '3d-highl':
            Model = LG3dAOHighl
            size = (size, size, size)
        elif dftModel == '3d-mf':
            Model = LG3dMf
            size = (size, size, size)

        densRange = [startDens, maxDens]

        if init[0:4] == 'dens':
            sourcePath = create_source_path(Model, size[0], epsi)
            init = Model.load_syst(sourcePath, init)

        generate_series(Model, size, epsi, densRange, stpWidth, init,
                alpha, it_steps, chp, saveSpec='ser', consecutive=cons)
    
    if mode == 'searchTrans':
        dftModel = str(sys.argv[2])
        size = int(sys.argv[3])
        epsi = float(sys.argv[4])
        minDens = float(sys.argv[5])
        maxDens = float(sys.argv[6])
        alpha = float(sys.argv[7])
        it_steps = int(sys.argv[8])
        chp = str(sys.argv[9])
        acuracy = float(sys.argv[10])
        initLeft = str(sys.argv[11])
        initRight = str(sys.argv[12])
    
        if dftModel == '2d-highl':
            Model = LG2dAOHighl
            size = (size, size)
        elif dftModel == '2d-mf':
            Model = LG2dMf
            size = (size, size)
        elif dftModel == '3d-highl':
            Model = LG3dAOHighl
            size = (size, size, size)
        elif dftModel == '3d-mf':
            Model = LG3dMf
            size = (size, size, size)

        densRange = [minDens, maxDens]
        sourcePath = create_source_path(Model, size[0], epsi)

        if initRight[0:4] == 'dens':
            initRight = Model.load_syst(sourcePath, initRight)
            specRight = 'inhProf'
        elif initRight in ['hom', 'sph', 'cyl', 'sl']:
            specRight=initRight
        if initLeft[0:4] == 'dens':
            initLeft = Model.load_syst(sourcePath, initLeft)
            specLeft = 'inhProf'
        elif initLeft in ['hom', 'sph', 'cyl', 'sl']:
            specLeft=initLeft

        search_trans(Model, size, epsi, densRange, alpha, it_steps,
                chp, acuracy, initLeft, initRight, specLeft=specLeft,
                specRight=specRight)

