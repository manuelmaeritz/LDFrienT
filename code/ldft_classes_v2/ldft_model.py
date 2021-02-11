import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from functools import reduce
import abc

class LdftModel(abc.ABC):
    """LdftModel is an abstract base class for lattice density functional
    theory (ldft) models. It provides an interface for every ldft model.
    This means functions specific to a certain model are predefined and
    must be implemented in a class inheriting "LdftModel". Functions
    applicable to any model (eg. picard updates) are completely
    implemented. This safes a lot of work and redundancy.
    Further it supports a structure to administer a lattice system
    through various instance variables and functions. If you want to
    define your one model, just create a class inheriting "LdftModel"
    and overwrite the abstract methods.
    
    Parameters
    ----------
    size : `Tuple`
        It defines the amount of lattice sites for each dimension.
    mu : `List`, optional: default = `None`
        Contains the chemical potential for every species. Supports also
        type `None` as entry (eg. [0.3, None, None]), if just the
        chemical potential of certain species should be given. See also
        ``mu_fix``. Choose `None` if you do not want to set any entry.
    dens : `List`; Optional: default = `None`
        Contains the systems average density for every species. Supports
        also `None` as entry (eg. [None, 0.4, 0.4]), if just the
        density of a certain species should be given. See also
        ``mu_fix``. Choose `None` if you do not want to set any entry.
    mu_fix : `List`
        Contains a boolean value for each species deciding whether the
        chemical potential should be kept fixed during picard iteration
        or not. `True` value treads the corresponding species as
        grand-canonical, `False`-species are treated canonical. The
        ``True``-values require a number value in the corresponding
        ``mu``-argument whereas the ``False``-values require a number
        value in the ``density``-argument at a picard-iteration.
    v_ext : `List`; Optional: default = `None`
        List of numpy.array of same shape as the ``size``-parameter.
        Contains an external potential for every species. Choose `None`,
        for zero external potential for every species.
    r : `List`; Optional: default = `None`
        List of numpy.array of same shape as the ``size``-parameter. Each
        array corresponds to the density profile of a species. Choose
        `None` in case you hand over the ``r_hist``-parameter or in case
        you do not want to set the variable yet.
    r_hist : `List`; Optional: default = `None`
        Picard-history of a density profile. It contains the density
        profiles for certain picard-steps of a system which has already
        been evolved through picard iteration. Every entry is of the
        format of the ``r``-parameter. Use `None` if the system has no
        history yet.
    err_hist : `List`;  Optional: default = `None`
        Contains the error at the picard-steps corresponding to the
        entries of `r_hist`. The entries are lists containing an error
        for every species. Use `None` if no history available.
    it_hist : `List`; Optional: default = `None`
        List of the picard steps corresponding to the density profiles at
        the ``r_hist``-parameter. Use `None` if no history available.
        Note: if ``r_hist`` is given then also this argument should be
        assigned with an appropriate list.
    bound_cond : `String`; Optional: default value 'periodic'
        Determines the boundary condition (bc). Values available:
        'periodic' for periodic bc, '11_if' for 2d systems with 45°
        tilted bc (to create 45° slab interfaces (11-interfaces)),
        '110_if' for 3d systems with a 45° tilted bc with respect to one
        axis (for 110-interfaces), '111_if' for 3d systems with a 45°
        tilted bc with respect to two axis (for 111-interfaces). For
        '11_if', '110_if' and '111_if' the ``size``-argument should be
        chosen in that way, the first two axis are of equal length and
        one and the last one is of twice that length (cuboid with square
        front face and long edge twice the short edges). If one wants to
        make use of the ``bound_cond``-argument one needs to use the
        ``_boundary_roll``-method for rolling the density profile
        instead of numpy.roll in the class implementing the specific
        model. The variable can be easily extended to further accepted
        values by adapting the ``self._boundary_roll``-method.
   """

    _size = None
    """Here the ``size``-parameter is stored. See it's description for
    further detail. (`Tuple`)
    """

    _mu = None
    """Here the ``mu``-parameter is stored if not `None`. If `None`, a
    list of `None` of the same length as the ``mu_fix``-parameter is
    assigned. This attribute is being updated after every picard update.
    See the description of the function `_make_picard_update` and the
    parameter ``mu`` for further detail. (`List`)
    """

    _dens = None
    """Here the ``dens``-parameter is stored if not `None`. If `None`, a
    list of `None` of the same length as the ``mu_fix``-parameter is
    assigned. This attribute is being updated after every picard update.
    See the description of the function `_make_picard_update` and the
    parameter ``dens`` for further detail. (`List`)
    """
    
    _mu_fix = None
    """Here ``mu_fix``-argument is stored. See it's description for
    further detail. (`List`)
    """

    _v_ext = None
    """Here ``v_ext``-parameter is stored, if not `None`. If `None`,
    then a list of the same length as the ``mu_fix``-parameter is
    assigned. The entries of which are zero-arrays of the shape of the
    ``size``-parameter. See the description of the parameter ``v_ext``
    for further details. (`List`)
    """

    _r = None
    """Here the ``r``-parameter is stored if not `None`. If `None`, then
    the last entry of the parameter ``r_hist`` is taken, when
    ``r_hist``!=`None`. If both parameters are `None` then
    ``_r``=`None`. This argument gets updated after every picard-update
    (see description of ``_make_picard_update``). For more details see
    description of the parameter `r`. (`List`)
    """

    _r_hist = None
    """Here the ``r_hist``-parameter is stored if not `None`. Otherwise
    an empty list is assigned. This argument gets updated after certain
    picard-steps (see description of ``make_picard_iteration``). For
    more detail see the description of the parameter ``r_hist``. (`List`)
    """

    _err_hist = None
    """Here the ``err_hist``-parameter is stored if not `None`. If
    `None`, then an empty list is assigned. This argument gets updated
    after certain picard-steps (see description of
    ``make_picard_iteration``). For more detail see the description of
    the parameter ``err_hist``. (`List`)
    """

    _it_hist = None
    """Here the ``it_hist``-parameter is stored if not `None`. If
    `None`, then an empty list is assigned. This argument gets updated
    after certain picard-steps (see description of
    ``make_picard_iteration``). For more detail see the description of
    the parameter ``it_hist``. (`List`)
    """

    _bound_cond = None
    """Here the ``bound_condition``-parameter is stored. See its
    description for further information. (`String`)
    """

    _it_counter = None
    """Counts the number of picard-updates the system has gone through.
    If the parameter ``it_hist`` is set, its last entry is taken as its
    initial value. Otherwise it is initialised with `0`.  It is updated
    after every picard-update (see description of
    ``_make_picard_update``). Every time the ``set_r`` function is
    called, ``_it_counter`` is being reset to `0`.
    """

    _dim = None
    """Dimension of the system. Evaluates the length of the
    ``size``-parameter.
    """

    def __init__(self, size, mu_fix, mu=None, dens=None, v_ext=None,
            r=None, r_hist=None, err_hist=None, it_hist=None,
            bound_cond='periodic'):
        self._size=size
        self._dim=len(size)
        self._bound_cond=bound_cond
        self._it_hist = [] if it_hist ==None else it_hist
        self._r_hist = [] if r_hist==None else r_hist
        self._it_counter = 0 if it_hist == None else it_hist[-1]
        if self._r_hist != []:
            self._r=self._r_hist[-1]
        else:
            if r==None: self._r=None
            else: self.set_r(r)
        self._err_hist = [] if err_hist==None else err_hist
        self._size = tuple(int(L_ax) for L_ax in size)
        self._mu = mu if mu!=None else [None for i in mu_fix]
        self._dens = dens if dens!=None else [None for i in mu_fix]
        self._mu_fix = mu_fix
        self._v_ext = v_ext if v_ext!=None else\
                [np.zeros(size) for i in range(len(self._mu_fix))]
    
    def __str__(self):
        descrStr = 'This is a LdftModel with the following properties:'
        Str0 ='{0:<40s}: {1}\n'.format('System size', self._size)
        Str1 ='{0:<40s}: {1}\n'.format('mu_fix', self._mu_fix)
        Str2 ='{0:<40s}: {1}\n'.format('Chem. pot. \'mu\'', self._mu)
        Str3 ='{0:<40s}: {1}\n'.format('Density', self._dens)
        Str4 ='{0:<40s}: {1}\n'.format('External potential \'V_ext\'',\
                'on' if np.any(self._v_ext) else 'off')
        Str5 ='{0:<40s}: {1}\n'.format('Current dens prof',type(self._r))
        Str6 ='{0:<40s}: len={1}\n'.format('History',len(self._it_hist))
        Str7 ='{0:<40s}: last entry={1}\n'.format('',\
                self._it_hist[-1] if len(self._it_hist)>0 else '---')
        Str8 ='{0:<40s}: len={1}\n'.format('Error history',\
                len(self._err_hist))
        Str9 ='{0:<40s}: last value={1}\n'.format('',\
                self._err_hist[-1] if len(self._err_hist)>0 else '---')
        Str10='{0:<40s}: {1}\n'.format('Boundary Condition',\
                self._bound_cond)
        return descrStr+'\n\n'+Str0+Str1+Str2+Str3+Str4+Str5+Str6\
                +Str7+Str8+Str9+Str10

    ####################################################################
    #It follows public properties for access of instance variables to
    #the end-user. Some of them are read only.
    ####################################################################

    @property
    def size(self):
        """See description of ``_size`` (`Tuple`, read-only)
        """
        return self._size

    @property
    def mu(self):
        """See description of ``_mu`` (`List`)
        """
        return self._mu

    @mu.setter
    def mu(self, mu):
        self._mu = mu

    @property
    def dens(self):
        """See description of ``_dens`` (`List`)
        """
        return self._dens

    @dens.setter
    def dens(self, dens):
        self._dens = dens

    @property
    def mu_fix(self):
        """See description of ``_mu_fix`` (`List`)
        """
        return self._mu_fix

    @mu_fix.setter
    def mu_fix(self, mu_fix):
        self._mu_fix = mu_fix

    @property
    def v_ext(self):
        """See description of ``_v_ext`` (`List`)
        """
        return self._v_ext

    @v_ext.setter
    def v_ext(self, v_ext):
        self._v_ext = v_ext

    @property
    def r(self):
        """See description of ``_r`` (`List`)
        The setter method calls the function ``set_r``
        """
        return self._r

    @r.setter
    def r(self, r):
        self.set_r(r)

    @property
    def r_hist(self):
        """See description of ``_r_hist`` (`List`, read-only)
        """
        return self._r_hist

    @property
    def err_hist(self):
        """See description of ``_err_hist`` (`str`, read-only)
        """
        return self._err_hist

    @property
    def it_hist(self):
        """See description of ``_it_hist`` (`List`, read-only)
        """
        return self._it_hist

    @property
    def bound_cond(self):
        """See description of ``_boundary_cond`` (`str`, read-only)
        """
        return self._bound_cond

    @property
    def it_counter(self):
        """See description of ``_it_counter`` (`int`, read-only)
        """
        return self._it_counter

    @property
    def dim(self):
        """See description of ``_dim`` (`int`, read-only)
        """
        return self._dim

    ####################################################################
    #It follows the functionals defining the model.
    ####################################################################

    @abc.abstractmethod
    def cal_F(self):
        """Calculates the free energy of the models current density
        profile (meaning every species treated canonical, as if
        ``_mu_fix`` is ``False`` for every species)

        Returns
        -------
         The free energy : `Float`
        """
        pass

    def cal_Om(self):
        """Calculates the grand potential of the models current density
        profile (meaning every species treated grand canonical, as if
        ``_mu_fix`` is ``True`` for every species).

        Returns
        -------
        The grand potential : `Float`
        """
        Om = self.cal_F()
        for i in range(len(self._mu_fix)):
            Om -= self._mu[i]*np.sum(self._r[i])
        return Om

    def cal_semi_Om(self):
        """Calculates the semi grand potential of the models current
        density profile (meaning every species with ``_mu_fix==True``
        is treated grand canonically and every other canonical).

        Returns
        -------
        The semi-grand potential : `Float`
        """
        semi_Om = self.cal_F()
        for i in range(len(self._mu_fix)):
            semi_Om -= self._mu[i]*np.sum(self._r[i])\
                    if self._mu_fix[i] else 0
        return semi_Om

    @abc.abstractmethod
    def cal_mu_ex(self):
        """Calculates the excess chemical potential of the models current
        density profile

        Returns
        -------
        The excess chemical potential : `List`
        """
        pass

    ####################################################################
    #It follows help-functions for implementing the model specific
    #functionals
    ####################################################################

    @classmethod
    def _tilted_roll_3d(cls, array, steps, roll_axis, shift, shift_axis):
        """Rolls a 3d numpy array in the manner of numpy.roll in
        direction of ``roll_axis``, but with different boundary
        conditions. The padding happens after the opposite surface, but
        shifted. The shift corresponds to another rolling in direction of
        a ``shift_axis`` unequal the ``shift_axis``.

        Parameters
        ----------
        array : `numpy.array`
            A 3d array which should be rolled.
        steps : `int`
            Number of steps of the rolling. Negative numbers for rolling
            in negative direction.
        roll_axis : `int`
            Axis in which direction should be rolled. Possible values:
            1, 2 and 3.
        shift : `int`
            Shift of the padding area with respect to the opposite
            surface of the array.
        shift_axis : `int`
            Axis in which the shift should be done. Possible values: 1,
            2 and 3 but not the same value as in ``roll_axis``.

        returns
        -------
        Rolled array : `numpy.array`
        """
        array=np.copy(array)
        l0=array.shape[0]
        l1=array.shape[1]
        l2=array.shape[2]
        if roll_axis==0:
            if steps>0:
                array[-steps:l0,:,:]=np.roll(array[-steps:l0,:,:],
                        shift, shift_axis)
            else:
                array[0:-steps,:,:]=np.roll(array[0:-steps,:,:], shift,
                        shift_axis)
        if roll_axis==1:
            if steps>0:
                array[:,-steps:l1,:]=np.roll(array[:,-steps:l1,:],
                        shift, shift_axis)
            else:
                array[:,0:-steps,:]=np.roll(array[:,0:-steps,:], shift,
                        shift_axis)
        if roll_axis==2:
            if steps>0:
                array[:,:,-steps:l2]=np.roll(array[:,:,-steps:l2],
                        shift, shift_axis)
            else:
                array[:,:,0:-steps]=np.roll(array[:,:,0:-steps], shift,
                        shift_axis)
        array = np.roll(array, steps, roll_axis)
        return array

    @classmethod
    def _tilted_roll(cls, array, steps, roll_axis, shift, shift_axis):
        """See the description of ``_tilted_roll_3d``. This function
        makes the same but independent of the dimension of the array
        which should be rolled.

        Parameters
        ----------
        array : `numpy.array`
        A 2d or 3d array which should be rolled.
        steps : `int`
            Number of steps of rolling. Negative numbers for rolling in
            negative direction.
        roll_axis : `int`
            Axis in which direction should be rolled.
        shift : `int`
            Shift of the padding area with respect to the opposite
            surface.
        shift_axis : `int`
            Axis in which the shift should be done.

        Returns
        -------
        The rolled array : `numpy.array`
        """
        if len(array.shape)==2:
            array = np.array([array])
            array = cls._tilted_roll_3d(array, steps, roll_axis+1, shift,
                    shift_axis+1)
            array = array[0, :,:]
            return array
        elif len(array.shape)==3:
            return cls._tilted_roll_3d(array, steps, roll_axis, shift,
                    shift_axis)

    def _boundary_roll(self, r, steps, axis=0):
        """Performs the rolling of a density profile under consideration
        of the boundary condition in the class variable ``_bound_cond``.
        If the boundary condition is not 'periodic', then the function
        ``_tilted_roll`` is applied in an appropriate way to satisfy the
        given boundary condition while rolling.

        Parameters
        ----------
        r : `numpy.array`
            The density profile which should be rolled.
        steps : `int`
            Number of steps of rolling. Negative numbers for rolling in
            negative direction.
        axis : `int`
            Axis in which direction should be rolled.

        Returns
        -------
        The rolled array : `numpy.array`
        """
        if self._bound_cond=='periodic':
            return np.roll(r, steps, axis)
        elif self._bound_cond=='111_if' or self._bound_cond=='11_if':
            shift_axis = self._size.index(max(self._size))
            if shift_axis == axis:
                return np.roll(r, steps, axis)
            else:
                shift = max(self._size)//2
                return self._tilted_roll(r, steps, axis, shift,
                        shift_axis)
        elif self._bound_cond=='110_if':
            shift_axis = self._size.index(max(self._size))
            if shift_axis == axis:
                return np.roll(r, steps, axis)
            elif (shift_axis+1)%3 == axis:
                return np.roll(r, steps, axis)
            else:
                shift = max(self._size)//2
                return self._tilted_roll(r, steps, axis, shift,
                        shift_axis)

    def _cal_Phi_id(self):
        """Calculates the ideal gas part of the free energy density.

        Returns
        -------
        Result : `numpy.ndarray`
        """
        Phi_id_func = lambda r: r*(np.log(r)-1)
        Phi_id_spec = list(map(Phi_id_func, self._r))
        Phi_id = reduce(lambda a, b: a+b, Phi_id_spec)
        return Phi_id

    @staticmethod
    def _cal_Phi_0(x):
        """Calculates the free energy density of a 0d-cavity depending
        on the packing fraction.

        Parameters
        ----------
        x : `float`
            The packing fraction at which the 0d-cavity is evaluated

        Returns
        -------
        Result : `float`
            The free energy density (Result is multiplied with the
            inverse temperature to make its dimension 1).
        """
        return x+(1-x)*np.log(1-x)

    @staticmethod
    def _cal_dPhi_0(x):
        """Calculates the derivative of the free energy density of a
        0d-cavity with respect of the packing fraction.

        Parameters
        ----------
        x : `float`
            The packing fraction

        Returns
        -------
        Result : `float`
            Derivative of the free energy density (Result is multiplied
            with the inverse temperature to make its dimension 1).
        """
        return -np.log(1-x)

    @staticmethod
    def _cal_d2Phi_0(x):
        """Calculates the second derivative of the free energy density
        of a 0d-cavity with respect of the packing fraction.

        Parameters
        ----------
        x : `float`
            The packing fraction

        Returns
        -------
        Result : `float`
            Second derivative of the free energy density (Result is
            multiplied with the inverse temperature to make its
            dimension 1).
        """
        return 1/(1-x)

    #####################################################################
    #It follows functions concerning the picard-iteration
    #####################################################################
    def _make_picard_update(self, alpha):
        """Runs one Picard-Iteration. The instance variable ``_mu_fix``
        decides whether the density or the chemical potential is to be
        kept fixed during the iteration. When ``_mu_fix[i]``==`False` for
        one species ``i``, the density is kept fix for this species and
        the ``_mu``-attribute for the same is updated. In case of `True`,
        the chemical potential ``_mu[i]`` is kept constant and the
        density `_dens[i]` is going to be updated. The variable `_r` is
        being updated, where the updated `r` is a superposition of the
        old ``_r`` and the iterated ``r``. The `alpha`-parameter steers
        the contribution of the iterated ``r`` to that superposition.
        Finally 'self._it_counter'-Variable is increased by one.

        Parameters
        ----------
        alpha : `Float`
            Value between 0 and 1. Determines how 'fast' the iteration
            is done (The higher, the faster). In case of to high
            ``alpha`` the danger of divergence arises.

        Returns
        -------
        r : `List`
            The iterated density profile.
        error : `List`
            The error for each species.
            In case of divergence prints 'divergent!!!' and returns nothing.
        """
        mu_ex = self.cal_mu_ex()
        temp = [np.exp(mu_ex[i]-self._v_ext[i])\
                for i in range(len(self._mu_fix))]
        r_new = []
        V= reduce(lambda a, b: a*b, self._size)
        for i,temp_i in enumerate(temp):
            if np.any(np.isnan(temp_i)):
                print('divergent!!!')
                return
            if self._mu_fix[i]:
                r=temp_i*np.exp(self._mu[i])
                r_new.append(r)
                self._dens[i]=np.mean(r)
            else:
                temp_dens = np.sum(temp_i)/V
                z = self._dens[i]/temp_dens
                r_new.append(temp_i*z)
                self._mu[i] = np.log(z)
        error = [np.sum((r_new[i]-self._r[i])**2)\
                for i in range(len(r_new))]
        r =[alpha*r_new[i]+(1-alpha)*self._r[i]\
                for i in range(len(r_new))]
        self._r=r
        self._it_counter +=1
        return r, error

    def make_picard_iteration(self, alpha, it_steps, checkp_method,\
            min_err=None):
        """Calls ``it_steps`` times the method ``_make_picard_update``
        with the update parameter ``alpha``. The iteration can be
        prematurely aborted when the iteration error fall below a minimal
        error ``min_err``. When ``self._it_counter`` reaches certain
        values (checkpoints) the current profile is appended to the
        ``self._r_hist``-attribute by calling ``_append_hist``. The next
        checkpoint is calculated by ``_set_new_checkp`` according to the
        parameter ``checkp_method``. Before exiting the function the last
        profile is also appended to ``_err_hist`` with ``_append_hist``.

        Parameters
        ----------
        alpha : `Float`
            Value between 0 and 1. Determines how 'fast' the iteration is
            done (The higher, the faster). In case of to high ``alpha``
            the danger of divergence arises.
        it_steps : `Int`
            Number of iteration steps
         checkp_method : `String`
            Determines in which intervals the profile should be
            appended to the ``_r_hist``-attribute. Possible values:
            integer number, 'exp#', 'dec#' where # needs to be replaced
            by a number. See description of ``_set_new_checkp``.
        min_err : `Float`
            Determines at which error the iteration can be aborted
            prematurely.
        """
        checkp = self._set_new_checkp(checkp_method)
        for i in np.arange(it_steps):
            r, error = self._make_picard_update(alpha)
            self._err_hist.append(error)
            if self._it_counter - checkp ==0:
                print('checkpoint at: {0:>10}, Error: {1:>20}'.\
                        format(self._it_counter, self._err_hist[-1][0]))
                self._append_hist()
                checkp = self._set_new_checkp(checkp_method)
            if min_err!=None and all(np.array(error)<min_err): break
        if self._it_counter != self._it_hist[-1]:
            self._append_hist()

    def _set_new_checkp(self, checkp_method):
        """Calculates the next 'checkpoint' meaning an iteration number
        at which the current density profile ``_r`` should be appended to
        ``self._r_hist``. The next checkpoint is determined by the current
        value of ``_it_counter`` and the method defined by the
        parameter ``checkp_method``. 
        
        Parameters
        ----------
        checkp_method : `String` or `Int`
            Determines how the next checkpoint is calculated. Recommended
            value: 'dec2'. It can take the following values:
            integer value (for equidistant checkpoints with interval of
            the integer); 'exp#' where # is to be replaced by a float
            value (next checkpoint is last checkp to the power of float);
            'dec#' where # is replaced by an integer (if e.g. #==3, the
            checkpoints goes like this: 30, 60, 90, 100, 300, 600, 900,
            1000, 3000, ...).

        Returns
        -------
        checkp : `Int`
            The calculated next checkpoint
        """
        if type(checkp_method)== int:
            checkp = self._it_hist[-1]+checkp_method
        else:
            method = str(checkp_method)[0:3]
            parm = float(checkp_method[3:len(checkp_method)])
            if method == 'exp':
                checkp = int(self._it_hist[-1]**parm)
                if checkp == self._it_hist[-1]: checkp += 1
            elif method =='dec':
                exp = int(np.floor(np.log10(self._it_hist[-1])))\
                        if self._it_hist[-1] != 0 else 1
                steps = parm*10**exp
                checkp = self._it_hist[-1]+steps
        return checkp

    ####################################################################
    #It follow functions for constructing initial density profiles
    ####################################################################

    def create_init_profile(self, dens=None, shape=None):
        """Creates an initial density profile for each species the
        picard iteration can start with. A list of average density of
        each species is handed over via the ``dens``-parameter.
        Additionally a nucleus can be placed in the density profile of
        each species, the shape of which determined by the
        ``shape``-parameter. Calls the function ``self.set_r`` to set
        the density profile to the variable ``_r``. The Nucleus further
        satisfies the boundary condition ``_bound_cond``

        Parameters
        ----------
        dens : `List`
            Determines the average density of each species.
        shape : `List` of `Tuples`
            The tuples determines the shape of the nucleus for each
            species. E.g. (3, 4) for a 2d-system with a nucleus of
            expand 3x4.
        """
        r = []
        for i in range(len(self._mu_fix)):
            r.append(self.return_nuc_densProfile(dens[i], shape[i]))
        self.set_r(r)

    def return_hom_densProfile(self, dens):
        """Returns a homogeneous one species density profile with
        density according to the parameter ``dens``. The shape of which
        is determined by the `_size`-instance variable.

        Parameters
        ----------
        dens : `Float`
            Density of the homogeneous profile.

        Returns
        -------
        Profile : `np.array`
            The resulting density profile.
        """
        grid = np.ones(self._size)
        profile = grid*dens
        return profile
    
    def return_nuc_densProfile(self, dens, shape):
        """Returns a one species density profile with average density
        according to the ``dens``-parameter and a nucleus of shape
        determined by the ``shape``-parameter. The nucleus further
        satisfies the boundary condition ``_bound_cond``.

        Parameters
        ----------
        dens : `Float`
            Average density of the profile.
        shape : `Tuple`
            Determines the shape of the nucleus. E.g. (3, 4) for a
            2d-system with a nucleus of expand 3x4.

        Returns
        -------
        Profile :`np.array`
            The density resulting profile.
        """
        profile = self.return_hom_densProfile(dens)
        nuc_dens = dens+0.05 if dens <0.5 else dens-0.05
        if self._dim==2:
            profile[(self._size[0]-shape[0])//2:\
                    (self._size[0]+shape[0])//2,\
                    (self._size[1]-shape[1])//2:\
                    (self._size[1]+shape[1])//2] = nuc_dens
            if self._bound_cond=='11_if':
                sheer = profile.shape[0]
                for i in range(sheer):
                    profile[i,:]=np.roll(profile[i,:], sheer-i)
        if self._dim==3:
            profile[(self._size[0]-shape[0])//2:\
                    (self._size[0]+shape[0])//2,\
                    (self._size[1]-shape[1])//2:\
                    (self._size[1]+shape[1])//2,\
                    (self._size[2]-shape[2])//2:\
                    (self._size[2]+shape[2])//2] = nuc_dens
            if self._bound_cond in ['111_if', '110_if']:
                sheer = profile.shape[0]
                for i in range(sheer):
                    profile[:,i,:]=np.roll(profile[:,i,:], sheer-i,\
                            axis=1)
            if self._bound_cond=='111_if':
                for i in range(sheer):
                    profile[i,:,:]=np.roll(profile[i,:,:], sheer-i,\
                            axis=1)
        profile = profile*dens/np.mean(profile)
        return profile
   
   #####################################################################
   #It follows functions to administer the instance variables. 
   #####################################################################

    def set_r(self, r):
        """This function is used for assigning a new initial profile
        ``r`` to the instance variable ``_r``. Therefor the
        ``_it_counter`` is being reset to '0' and the history
        attributes ``_r_hist``, ``_it_hist``, ``_err_hist`` are updated.
        
        Parameters
        ----------
        r : `List` of `numpy.array`
            New initial density profile for each species.

        """
        self._r=r
        self._it_counter=0
        self._r_hist =[r]
        self._it_hist =[0]
        self._err_hist =[]

    def set_hist(self, r_hist, it_hist, err_hist):
        """This function is to manually set the internal history
        variables ``_r_hist``, ``_it_hist`` and ``_err_hist``. The last
        entry of the ``r_hist``-parameter is assigned to the instance
        variable ``_r``, which is the current density profile.

        Parameters
        ----------
        r_hist : `list` of `list` of `numpy.ndarray`
            Iteration history of the density profile. This parameter
            should be of the following format [profile_0, profile_1,...]
            where ``profile_i`` is the profile of the i'th iteration
            step and has the format [r_1, r_2, ...], where the entries
            are the profile of the corresponding species.
        it_hist : `list` of `int`
            This parameter lists the corresponding iteration steps of
            the ``r_hist`` parameter.
        err_hist : `list` of `list` of `float`
            History of the picard error. It is of the following format:
            [err_0, err_1,...] where err_i is the error of the i'th
            iteration step and is a list itself, with an error entry for
            every species.
        """
        self._it_hist = it_hist
        self._r_hist =r_hist
        self._err_hist = err_hist
        self._r = self._r_hist[-1]
    
    def _append_hist(self):
        """Updates the history variables ``_r_hist``,``_it_hist``, by
        appending the current density profile ``_r`` to ``_r_hist``
        and appending ``_it_counter`` to ``_it_hist``.
        """
        self._r_hist.append(self._r)
        self._it_hist.append(self._it_counter)
    
    def save_syst(self, path, filename):
        """Uses ``pickle.dump`` to save the instance variables of a
        system.

        Parameters
        ----------
        path : `String`
            Directory in which the system should be stored (needs to be
            a absolute path)
        filename : `String`
            The filename under which the system should be stored.
        """
        if not os.path.isdir(path): os.makedirs(path, exist_ok=True)
        f = open(os.path.join(path, filename), 'bw')
        pickle.dump(self, f)
        f.close()
    
    @classmethod
    def load_syst(cls, path, filename):
        """Uses ``pickle.load`` to load a system. It is strongly
        recommended to override this method in the inherited classes,
        as the returned system might be of an outdated type! A typecast
        should be implemented!

        Parameters
        ----------
        path : `String`
            Directory in which the system is stored which one want's to
            load (needs to be a absolute path)
        filename: `String`
            The filename under which the system of interest is stored.
        
        Returns
        -------
        Model : `LdftModel`
            The returned model probably has the type of an inherited
            class. It might also be the class of an outdated type.
        """
        path=os.path.join(path, filename)
        f = open(path, 'rb')
        syst = pickle.load(f)
        return syst

    ####################################################################
    #It follow functions supporting some predefined plots
    ####################################################################

    def print_error(self):
        """Returns a figure where the error history ``_err_hist`` is
        plotted.

        Returns
        -------
        Figure : `matplotlib.pyplot.figure`
            Plotted error history.
        """
        err = np.array(self._err_hist)
        height = 5
        width = height*1.62
        fig = plt.figure(figsize=(width, height))
        ax = fig.add_subplot(111)
        ax.set_yscale('log')
        ax.set_ylabel('err')
        ax.set_xlabel('it_step')
        for i in range(len(self._mu_fix)):
            err_i = err[:,i]
            ax.plot(err_i, label='Species '+str(i))
        ax.legend()
        return fig

    def print_2d_profile(self):
        """Creates a figure where the current profile is plotted. This
        function is just for 2d-systems.

        Returns
        -------
        Figure : `matplotlib.pyplot.figure`
            Plotted profile
        """
        rows = len(self._mu_fix)
        x_fix = int(self._size[0]//2)
        y_fix = int(self._size[1]//2)
        cols = 3
        fig = plt.figure(figsize=(15, 3*rows))
        for i in range(rows):
            r = self._r[i]
            ax1 = fig.add_subplot(rows, cols, cols*i+1,\
                    ylabel='Species '+str(i))
            ax1.imshow(r)
            #ax1.set_axis_off()
            ax2 = fig.add_subplot(rows, cols, cols*i+2,\
                    ylabel='Density')
            ax2.plot(r[x_fix,:])
            #ax2.set_ylim((0,1))
            ax3 = fig.add_subplot(rows, cols, cols*i+3,sharey=ax2)
            ax3.plot(r[:,y_fix])
        fig.axes[0].set_title('2D-Plot')
        fig.axes[1].set_title('x='+str(x_fix))
        fig.axes[2].set_title('y='+str(y_fix))
        fig.axes[-2].set_xlabel('y-coord')
        fig.axes[-1].set_xlabel('x-coord')
        return fig
    
    def print_2d_profile2(self):
        """Creates a figure where the current profile is plotted. This
        function is just for 2d-systems.

        Returns
        -------
        Figure : `matplotlib.pyplot.figure`
            Plotted profile
        """
        rows = len(self._mu_fix)
        cols = 2
        x_fix = int(self._size[0]//2)
        y_fix = int(self._size[1]//2)
        fig = plt.figure(figsize=(10, 3*rows))
        for i in range(rows):
            r = self._r[i]
            ax1 = fig.add_subplot(rows, cols, cols*i+1,\
                    ylabel='Species '+str(i))
            ax1.imshow(r)
            #ax1.set_axis_off()
            ax2 = fig.add_subplot(rows, cols, cols*i+2, ylabel='Density')
            ax2.plot(r[x_fix,:], label='x='+str(x_fix))
            ax2.plot(r[:,y_fix], label='y='+str(y_fix))
            ax2.legend()
        fig.axes[0].set_title('2D-Plot')
        fig.axes[1].set_title('x/y fixed')
        fig.axes[-1].set_xlabel('x/y-coord')
        return fig
        
    def print_2d_hist(self, species=0, rows=10, idx_list=None):
        """Creates a figure where the history ``_r_hist`` is plotted.
        Just one species can be plotted at the same time. Not the total
        history is plotted but certain iteration steps.

        Parameters
        ----------
        species : `int`; optional: default = 0
            The species, the iteration-history of which shall be
            plotted.
        rows : `int`; optional: default = 10
            Number of iteration-steps which shall be plotted. This
            parameter is just be considered when the parameter
            ``idx_list`` is `None`.
        idx_list : `List`; optional: default = None
            If `None`, the iteration steps which are plotted are chosen
            equidistant in the ``_it_hist``-list. Alternatively one can
            choose ones own list. This list, however, does not contain
            the iteration-steps which shall be plotted, but the indices
            of those.

        Returns
        -------
        Figure : `matplotlib.pyplot.figure`
            Plotted history
        """
        if idx_list == None:
           idx_list = np.linspace(0, len(self._r_hist)-1, rows)
           idx_list =idx_list.astype(int)
           idx_list=list(idx_list)
        else: rows=len(idx_list)
        fig, ax=plt.subplots(nrows=rows, ncols=3)
        fig.set_figheight(3*rows)
        fig.set_figwidth(15)
        x_fix = int(self._size[0]//2)
        y_fix = int(self._size[1]//2)
        for row, idx in enumerate(idx_list):
            r = self._r_hist[idx][species]
            ax[row, 0].imshow(self._r_hist[idx][species])
            ax[row, 0].set_ylabel('it step: '+str(self._it_hist[idx]))
            ax1 = ax[row, 1]
            ax1.plot(r[x_fix,:])
            ax1.set_ylabel('Density')
            ax1.set_ylim((0,1))
            ax2=ax[row, 2]
            ax2.plot(r[:,y_fix])
            ax2.set_ylim((0,1))
        ax[0, 0].set_title('2D-Plot')
        ax[0, 1].set_title('x='+str(x_fix))
        ax[0, 2].set_title('y='+str(y_fix))
        ax[-1, -2].set_xlabel('y-coord')
        ax[-1, -1].set_xlabel('x-coord')
        return fig
    
    def print_2d_hist2(self, species=0, rows=10, idx_list=None):
        """Creates a figure where the history ``_r_hist`` is plotted.
        Just one species can be plotted at the same time. Not the total
        history is plotted but certain iteration steps.

        Parameters
        ----------
        species : `int`; optional: default = 0
            The species, the iteration-history of which shall be
            plotted.
        rows : `int`; optional: default = 10
            Number of iteration-steps which shall be plotted. This
            parameter is just be considered when the parameter
            ``idx_list`` is `None`.
        idx_list : `List`; optional: default = None
            If `None`, the iteration steps which are plotted are chosen
            equidistant in the ``_it_hist``-list. Alternatively one can
            choose ones own list. This list, however, does not contain
            the iteration-steps which shall be plotted, but the indices
            of those.

        Returns
        -------
        Figure : `matplotlib.pyplot.figure`
            Plotted history
        """
        if idx_list == None:
           idx_list = np.linspace(0, len(self._r_hist)-1, rows)
           idx_list =idx_list.astype(int)
           idx_list=list(idx_list)
        else: rows=len(idx_list)
        fig, ax=plt.subplots(nrows=rows, ncols=2)
        fig.set_figheight(3*rows)
        fig.set_figwidth(10)
        x_fix = int(self._size[0]//2)
        y_fix = int(self._size[1]//2)
        for row, idx in enumerate(idx_list):
            r = self._r_hist[idx][species]
            ax[row, 0].imshow(self._r_hist[idx][species])
            ax[row, 0].set_ylabel('it_step: '+str(self._it_hist[idx]))
            ax1 = ax[row, 1]
            ax1.plot(r[x_fix,:], label='x='+str(x_fix))
            ax1.plot(r[:,y_fix], label='y='+str(y_fix))
            ax1.set_ylabel('Density')
            ax1.set_ylim((0,1))
            ax1.legend()
        ax[0, 0].set_title('2D-Plot')
        ax[0, 1].set_title('x/y='+str(x_fix)+'/'+str(y_fix))
        ax[-1, -1].set_xlabel('x/y-coord')
        return fig

    ####################################################################
    #It follows functions concerning the surface properties of the
    #system.
    ####################################################################
    
    @abc.abstractmethod
    def _cal_p(self, dens):
        """Calculates the pressure for a bulk system with given densities
        for each species. The other parameters (temperature, attraction
        strength, etc.) are taken from the current instance ``self``.

        Parameters
        ----------
        dens : `List`
            The density for each species.

        Returns
        -------
        The pressure : `Float` 
        """
        pass

    @abc.abstractmethod
    def _cal_coex_dens(self):
        """Calculates the coexisting densities of bulk system for each
        species under the parameters of the current instance ``self``.

        Returns
        -------
        Coexisting densities : `List` of `Tuple`
            The coexisting densities arranged in a List of Tuples. Each
            species corresponds to a Tuple of the form:
            (vapour_dens, liquid_dens)
        """
        pass

    def cal_p_vap(self):
        """Calculates the coexisting pressures under the current
        parameters of the system (``_mu``, ``_dens``) and returns the
        vapour pressure.

        Returns
        -------
        Vapour pressure : `Float`
            The vapour pressure of the current system
        """
        r_coex = self._cal_coex_dens()
        r_v = [r[0] for r in r_coex]
        p = self._cal_p(r_v)
        return p

    def cal_p_liq(self):
        """Calculates the coexisting pressures under the current
        parameters of the system (``_mu``, ``_dens``) and returns the
        liquid pressure.

        Returns
        -------
        Liquid pressure : `Float`
            The vapour pressure of the current system
        """
        r_coex = self._cal_coex_dens()
        r_v = [r[1] for r in r_coex]
        p = self._cal_p(r_v)
        return p

    def det_intface_shape(self):
        """Determines the shape of the interface of the current
        configuration. It requires the inhomogeneities to be centered in
        the system.

        Returns
        -------
        Shape : `String`
            The shape of the interface: 'Droplet', 'Cylinder', 'Slab',
            'Homogeneous'
        """
        if self._dim==2:
            devX = np.mean(self._r[0][0,:]
                    -self._r[0][self._size[0]//2,:])
            devY = np.mean(self._r[0][:,0]
                    -self._r[0][:,self._size[1]//2])
            devZ=1
        elif self._dim==3:
            devX = np.mean(self._r[0][0,:,:]
                    -self._r[0][self._size[0]//2,:,:])
            devY = np.mean(self._r[0][:,0,:]
                    -self._r[0][:,self._size[1]//2,:])
            devZ = np.mean(self._r[0][:,:,0]
                    -self._r[0][:,:,self._size[1]//2])
        devArray = abs(np.array([devX, devY, devZ]))
        eqAxes = np.count_nonzero(devArray<10**-6)
        if eqAxes == 0:
            return 'Droplet'
        elif eqAxes == 1:
            return 'Cylinder' if self._dim==3 else 'Slab'
        elif eqAxes == 2:
            return 'Slab' if self._dim==3 else 'Homogeneous'
        elif eqAxes == 3:
            return 'Homogeneous'

    def cal_del_Om(self):
        """Calculates the delta between the current grand potential and
        the one by a homogeneous system of (oversaturated) vapor with the
        same chemical potential as the current system.

        Returns
        -------
        delta Omega : `Float`
            Delta of the grand potential
        """
        V= reduce(lambda a, b: a*b, self._size)
        Om = self.cal_Om()
        p_v = self.cal_p_vap()
        del_Om = Om + p_v*V
        return del_Om
    
    def cal_R_s(self):
        """Calculates the radius of surface of tension. In case of a
        Cylinder configuration in three dimensions, the cylinder has to
        point in the 0th axis of the density profile ``self._r``.

        Returns
        -------
        Radius of s.o.t. : `Float`
            Radius of surface of tension
        """
        del_Om = self.cal_del_Om()
        p_v = self.cal_p_vap()
        p_l = self.cal_p_liq()
        del_p = p_l-p_v
        ifShape = self.det_intface_shape()
        if self._dim==3 and ifShape=='Droplet':
            R = (3*del_Om/2/del_p/np.pi)**(1/3)
        elif self._dim==2 and ifShape=='Droplet':
            R = (del_Om/del_p/np.pi)**(1/2)
        elif self._dim==3 and ifShape=='Cylinder':
            h = self._size[0] 
            R = (del_Om/del_p/np.pi/h)**(1/2)
        return R

    def cal_R_em(self, em_species=0):
        """Calculates the equimolar radius for the species given by
        ``em_species``. In case of cylinder configurations in three
        dimensions the cylinder has to point in the 0th axis of the
        density profile ``self._r``. This function does only work
        properly, if a droplet/cylinder is embedded in a supersaturated
        vapour. For configurations of bubbles or vapour cylinders
        embedded in liquid, the result will be wrong.

        Parameters
        ----------
        em_species : `Int`; Optional: default=0
            Decides for which species the equimolar radius should
            be calculated
        
        Returns
        -------
        equimolar radius : `Float`
            Radius or the equimolar surface of a specific species.
        """
        r_mean = np.mean(self._r[em_species])
        r_liq = np.max(self._r[em_species])
        r_vap = np.min(self._r[em_species])
        x = (r_mean-r_vap)/(r_liq-r_vap)
        V= reduce(lambda a, b: a*b, self._size)
        ifShape = self.det_intface_shape()
        if self._dim == 3 and ifShape=='Droplet':
            R_em = (x*V*3/4/np.pi)**(1/3)
        elif self._dim == 2 and ifShape=='Droplet':
            R_em = np.sqrt(x*V/np.pi)
        elif self._dim == 3 and ifShape=='Cylinder':
            h = self._size[0]
            R_em = np.sqrt(x*V/np.pi/h)
        return R_em

    def cal_gamma_R(self, R):
        """Calculates the surface tension for spheres/circles in 3d/2d of
        radius ``R``. In case of cylinder configurations in three
        dimensions the cylinder has to point in the 0th axis of the
        density profile ``self._r``.

        Parameters
        ----------
        R : `Float`
            Radius at which the surface tension should be calculated
        
        Returns
        -------
        surface tension : `Float`
            surface tension for radius R.
        """
        del_p = self.cal_p_liq()-self.cal_p_vap()
        del_Om = self.cal_del_Om()
        ifShape = self.det_intface_shape()
        if self._dim == 3 and ifShape=='Droplet':
            A = 4*np.pi*R**2
            V = 4*np.pi*R**3/3
        elif self._dim ==2 and ifShape=='Droplet':
            A = 2*np.pi*R
            V = np.pi*R**2
        elif self._dim ==3 and ifShape=='Cylinder':
            h = self._size[0]
            A = 2*np.pi*R*h
            V = np.pi*R**2*h
        gamma = del_Om/A+del_p*V/A
        return gamma

    def cal_gamma_s(self):
        """Calculates the surface tension for spheres/circles in 3d/2d at
        the surface of tensions. In case of cylinder configurations in
        three dimensions the cylinder has to point in the 0th axis of the
        density profile ``self._r``.

        Returns
        -------
        surface tension : `Float`
            surface tension at the surface of tension.
        """
        del_Om = self.cal_del_Om()
        del_p = self.cal_p_liq()-self.cal_p_vap()
        ifShape = self.det_intface_shape()
        if self._dim==3 and ifShape=='Droplet':
            gamma_s=(del_Om*del_p**2*3/np.pi/16)**(1/3)
        elif self._dim==2 and ifShape=='Droplet':
            gamma_s=np.sqrt(del_Om*del_p/np.pi)
        elif self._dim==3 and ifShape=='Cylinder':
            h=self._size[0]
            gamma_s=np.sqrt(del_Om*del_p/np.pi/h)
        return gamma_s

    def cal_gamma_em(self, species=0):
        """Calculates the surface tension for spheres/circles in 3d/2d at
        the equimolar surface of a given species. In case of cylinder
        configurations in three dimensions the cylinder has to point in
        the 0th axis of the density profile ``self._r``. This function
        does only work properly, if a droplet/cylinder is embedded in a
        supersaturated vapour. For configurations of bubbles or vapour
        cylinders embedded in liquid, the result will be wrong.
        
        Parameters
        ----------
        species : `Int`; Optional: default=0
            species for the equimolar surface
        
        Returns
        -------
        Surface tension : `Float`
            Surface tension at the equimolar surface.
        """
        R_em = self.cal_R_em(em_species=species)
        gamma_em=self.cal_gamma_R(R_em)
        return gamma_em

    def cal_adsorptionAtSurfOfTens(self, species=0):
        """Calculates the adsorption for spheres/circles in 3d/2d at the
        surface of tension for a given species. This function
        does only work properly, if a droplet/cylinder is embedded in a
        supersaturated vapour. For configurations of bubbles or vapour
        cylinders embedded in liquid, the result will be wrong.

        Parameters
        ----------
        species : `Int`; Optional: default =0
            Species for which the adsorption should be calculated
        
        Returns
        -------
        Adsorption : `Tuple` of `Float`
            First entry: Adsorbed particle number; Second entry:
            adsorption.
        """
        r_mean = np.mean(self._r[species])
        r_liq = np.max(self._r[species])
        r_vap = np.min(self._r[species])
        R_s = self.cal_R_s()
        ifShape = self.det_intface_shape()
        if self._dim==3 and ifShape=='Droplet':
            V_dr = 4*np.pi*R_s**3/3
            A = 4*np.pi*R_s**2
        elif self._dim==2 and ifShape=='Droplet':
            V_dr = np.pi*R_s**2
            A = 2*np.pi*R_s
        elif self._dim ==3 and ifShape=='Cylinder':
            h = self._size[0]
            A = 2*np.pi*R_s*h
            V_dr = np.pi*R_s**2*h
        V= reduce(lambda a, b: a*b, self._size)
        V_sur= V-V_dr
        N_x = V*r_mean-V_dr*r_liq-V_sur*r_vap
        adsorption = N_x/A
        return N_x, adsorption

    def cal_gamma_inf(self, area):
        """Calculates the surface tension of a flat interface. This
        function can not determine the area of the surface itself.
        Therefore it has to be passed as parameter.

        Parameters
        ----------
        area : `float`
            Area of the surface. There are always two surfaces
            separating the liquid and the vapour. Meant is the area of
            one of those
        
        **returns:**
            ``Tuple`` of ``Float``: First entry: Adsorbed particle
            number; Second entry: Adsorption.
        """
        del_Om = self.cal_del_Om()
        gamma_inf=del_Om/2/area
        return gamma_inf
