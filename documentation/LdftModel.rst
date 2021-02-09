Documentation on LdftModel
==========================

*class* ``ldft_classes_v2.ldft_model.LdftModel(size, mu_fix, mu=None, dens=None, v_ext=None, r=None, r_hist=None, err_hist=None, it_hist=None, bound_cond='periodic')``

Base: ``abc.ABC``

LdftModel is an abstract base class for lattice density functional
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
    
Attributes
----------
_size : `Tuple`
    Here the ``size``-parameter is stored. See it's description for
    further detail.

_mu :  `List`
    Here the ``mu``-parameter is stored if not `None`. If `None`, a
    list of `None` of the same length as the ``mu_fix``-parameter is
    assigned. This attribute is being updated after every picard update.
    See the description of the function `_make_picard_update` and the
    parameter ``mu`` for further detail.

_dens : `List`
    Here the ``dens``-parameter is stored if not `None`. If `None`, a
    list of `None` of the same length as the ``mu_fix``-parameter is
    assigned. This attribute is being updated after every picard update.
    See the description of the function `_make_picard_update` and the
    parameter ``dens`` for further detail.

_mu_fix : `List`
    Here ``mu_fix``-argument is stored. See it's description for
    further detail.

_v_ext : `List`
    Here ``v_ext``-parameter is stored, if not `None`. If `None`,
    then a list of the same length as the ``mu_fix``-parameter is
    assigned. The entries of which are zero-arrays of the shape of the
    ``size``-parameter. See the description of the parameter ``v_ext``
    for further details.

_r : `List`
    Here the ``r``-parameter is stored if not `None`. If `None`, then
    the last entry of the parameter ``r_hist`` is taken, when
    ``r_hist``!=`None`. If both parameters are `None` then
    ``_r``=`None`. This argument gets updated after every picard-update
    (see description of ``_make_picard_update``). For more details see
    description of the parameter `r`.

_r_hist : `List`
    Here the ``r_hist``-parameter is stored if not `None`. Otherwise
    an empty list is assigned. This argument gets updated after certain
    picard-steps (see description of ``make_picard_iteration``). For
    more detail see the description of the parameter ``r_hist``.

_err_hist : `List`
    Here the ``err_hist``-parameter is stored if not `None`. If
    `None`, then an empty list is assigned. This argument gets updated
    after certain picard-steps (see description of
    ``make_picard_iteration``). For more detail see the description of
    the parameter ``err_hist``.

_it_hist : `List`
    Here the ``it_hist``-parameter is stored if not `None`. If
    `None`, then an empty list is assigned. This argument gets updated
    after certain picard-steps (see description of
    ``make_picard_iteration``). For more detail see the description of
    the parameter ``it_hist``.

_bound_cond : `String`
    Here the ``bound_condition``-parameter is stored. See its
    description for further information.

_it_counter : `integer`
    Counts the number of picard-updates the system has gone through.
    If the parameter ``it_hist`` is set, its last entry is taken as its
    initial value. Otherwise it is initialised with `0`.  It is updated
    after every picard-update (see description of
    ``_make_picard_update``). Every time the ``set_r`` function is
    called, ``_it_counter`` is being reset to `0`.

_dim : `integer`
    Dimension of the system. Evaluates the length of the
    ``size``-parameter.

Properties
----------

size : `Tuple`, read-only
    Accesses the ``_size``-attribute 

mu : `List`, read and write
    Accesses the ``_mu``-attribute

dens : `List`, read and write
    Accesses the ``_dens``-attribute

mu_fix : `List`, read and write
    Accesses the ``_mu_fix``-attribute

v_ext : `List`, read and write
    Accesses the ``_v_ext``-attribute

r : `List`, read and write
    Read accesses the ``_r``-attribute
    The setter method calls the function ``set_r``

r_hist : `List`, read-only
    Accesses the ``_r_hist``-attribute 

err_hist : `string`, read-only
    Accesses the ``_err_hist``-attribute

it_hist : `List`, read-only
    Accesses the ``_it_hist``-attribute

bound_cond : `string`, read-only
    Accesses the ``_boundary_cond``-attribute

it_counter : `int`, read-only
    Accesses the ``_it_counter``-attribute

dim : `int`, read-only
    Accesses the ``_dim``-attribute

Methodes
--------

``__init__(self, size, mu_fix, mu=None, dens=None, v_ext=None, r=None, r_hist=None, err_hist=None, it_hist=None, bound_cond='periodic')``
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
**Parameters**

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


``__str__(self)``
'''''''''''''''''


*abstractmethod* ``cal_F(self)``
''''''''''''''''''''''''''''''''
Calculates the free energy of the models curent density
profile (meaning every species treated canonical, as if
``_mu_fix`` is ``False`` for every species)

**Returns**

The free energy : `Float`


``cal_Om(self)``
''''''''''''''''
Calculates the grand potential of the models curent density
profile (meaning every species treated grand canonicaly, as if
``_mu_fix`` is ``True`` for every species).

**Returns**

The grand potential : `Float`

``cal_semi_Om(self)``
'''''''''''''''''''''
Calculates the semi grand potential of the models current
density profile (meaning every species with ``_mu_fix==True``
is treated grand canonically and every other canonical).

**Returns**

The semi-grand potential : `Float`

*abstractmethod* cal_mu_ex(self)
''''''''''''''''''''''''''''''''
Calculates the excess chemical potential of the models current
density profile

**Returns**

The excess chemical potential : `List`


*classmethod* ``_tilted_roll_3d(cls, array, steps, roll_axis, shift, shift_axis)``
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Rolls a 3d numpy array in the manner of numpy.roll in
direction of ``roll_axis``, but with different boundary
conditions. The padding happens after the opposite surface, but
shifted. The shift corresponds to another rolling in direction of
a ``shift_axis`` unequal the ``shift_axis``.

**Parameters**

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

**returns**

Rolled array : `numpy.array`

*classmethod* ``_tilted_roll(cls, array, steps, roll_axis, shift, shift_axis)``
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
See the description of ``_tilted_roll_3d``. This function
makes the same but independent of the dimension of the array
which should be rolled.

**Parameters**

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

**Returns**

The rolled array : `numpy.array`

``_boundary_roll(self, r, steps, axis=0)``
''''''''''''''''''''''''''''''''''''''''''
Performs the rolling of a density profile under consideration
of the boundary condition in the class variable ``_bound_cond``.
If the boundary condition is not 'periodic', then the function
``_tilted_roll`` is applied in an appropriate way to satisfy the
given boundary condition while rolling.

**Parameters**

r : `numpy.array`
    The density profile which should be rolled.
steps : `int`
    Number of steps of rolling. Negative numbers for rolling in
    negative direction.
axis : `int`
    Axis in which direction should be rolled.

**Returns**

The rolled array : `numpy.array`


``_cal_Phi_id(self)``
'''''''''''''''''''''
Calculates the ideal gas part of the free energy density.

**Returns**

Result : `numpy.ndarray`

*Staticmethod* ``_cal_Phi_0(x)``
''''''''''''''''''''''''''''''''
Calculates the free energy density of a 0d-cavity depending
on the packing fraction.

**Parameters**

x : `float`
    The packing fraction at which the 0d-cavity is evaluated

**Returns**

Result : `float`
    The free energy density (Result is multiplied with the
    inverse temperature to make its dimension 1).

*staticmethod* ``_cal_dPhi_0(x)``
'''''''''''''''''''''''''''''''''''''
Calculates the derivative of the free energy density of a
0d-cavity with respect of the packing fraction.

**Parameters**

x : `float`
    The packing fraction

**Returns**

Result : `float`
    Derivative of the free energy density (Result is multiplied
    with the inverse temperature to make its dimension 1).

*staticmethod* ``_cal_d2Phi_0(x)``
''''''''''''''''''''''''''''''''''
Calculates the second derivative of the free energy density
of a 0d-cavity with respect of the packing fraction.

**Parameters**

x : `float`
    The packing fraction

**Returns**

Result : `float`
    Second derivative of the free energy density (Result is
    multiplied with the inverse temperature to make its
    dimension 1).


``_make_picard_update(self, alpha)``
''''''''''''''''''''''''''''''''''''
Runs one Picard-Iteration. The instance variable ``_mu_fix``
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

**Parameters**

alpha : `Float`
    Value between 0 and 1. Determines how 'fast' the iteration
    is done (The higher, the faster). In case of to high
    ``alpha`` the danger of divergence arises.

**Returns**

r : `List`
    The iterated density profile.
error : `List`
    The error for each species.
    In case of divergence prints 'divergent!!!' and returns nothing.

``make_picard_iteration(self, alpha, it_steps, checkp_method, min_err=None)``
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Calls ``it_steps`` times the method ``_make_picard_update``
with the update parameter ``alpha``. The iteration can be
prematurely aborted when the iteration error fall below a minimal
error ``min_err``. When ``self._it_counter`` reaches certain
values (checkpoints) the current profile is appended to the
``self._r_hist``-attribute by calling ``_append_hist``. The next
checkpoint is calculated by ``_set_new_checkp`` according to the
parameter ``checkp_method``. Before exiting the function the last
profile is also appended to ``_err_hist`` with ``_append_hist``.

**Parameters**

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

``_set_new_checkp(self, checkp_method)``
''''''''''''''''''''''''''''''''''''''''
Calculates the next 'checkpoint' meaning an iteration number
at which the current density profile ``_r`` should be appended to
``self._r_hist``. The next checkpoint is determined by the current
value of ``_it_counter`` and the method defined by the
parameter ``checkp_method``. 

**Parameters**

checkp_method : `String` or `Int`
    Determines how the next checkpoint is calculated. Recommended
    value: 'dec2'. It can take the following values:
    integer value (for equidistant checkpoints with interval of
    the integer); 'exp#' where # is to be replaced by a float
    value (next checkpoint is last checkp to the power of float);
    'dec#' where # is replaced by an integer (if e.g. #==3, the
    checkpoints goes like this: 30, 60, 90, 100, 300, 600, 900,
    1000, 3000, ...).

**Returns**

checkp : `Int`
    The calculated next checkpoint


``create_init_profile(self, dens=None, shape=None)``
''''''''''''''''''''''''''''''''''''''''''''''''''''
Creates an initial density profile for each species the
picard iteration can start with. A list of average density of
each species is handed over via the ``dens``-parameter.
Additionally a nucleus can be placed in the density profile of
each species, the shape of which determined by the
``shape``-parameter. Calls the function ``self.set_r`` to set
the density profile to the variable ``_r``. The Nucleus further
satisfies the boundary condition ``_bound_cond``

**Parameters**

dens : `List`
    Determines the average density of each species.
shape : `List` of `Tuples`
    The tuples determines the shape of the nucleus for each
    species. E.g. (3, 4) for a 2d-system with a nucleus of
    expand 3x4.

``return_hom_densProfile(self, dens)``
''''''''''''''''''''''''''''''''''''''
Returns a homogeneous one species density profile with
density according to the parameter ``dens``. The shape of which
is determined by the `_size`-instance variable.

**Parameters**

dens : `Float`
    Density of the homogeneous profile.

**Returns**

Profile : `np.array`
    The resulting density profile.

``return_nuc_densProfile(self, dens, shape)``
'''''''''''''''''''''''''''''''''''''''''''''
Returns a one species density profile with average density
according to the ``dens``-parameter and a nucleus of shape
determined by the ``shape``-parameter. The nucleus further
satisfies the boundary condition ``_bound_cond``.

**Parameters**

dens : `Float`
    Average density of the profile.
shape : `Tuple`
    Determines the shape of the nucleus. E.g. (3, 4) for a
    2d-system with a nucleus of expand 3x4.

**Returns**

Profile :`np.array`
    The density resulting profile.

``set_r(self, r)``
''''''''''''''''''
This function is used for assigning a new initial profile
``r`` to the instance variable ``_r``. Therefor the
``_it_counter`` is being reset to '0' and the history
attributes ``_r_hist``, ``_it_hist``, ``_err_hist`` are updated.

**Parameters**

r : `List` of `numpy.array`
    New initial density profile for each species.

``set_hist(self, r_hist, it_hist, err_hist)``
'''''''''''''''''''''''''''''''''''''''''''''
This function is to manually set the internal history
variables ``_r_hist``, ``_it_hist`` and ``_err_hist``. The last
entry of the ``r_hist``-parameter is assigned to the instance
variable ``_r``, which is the current density profile.

**Parameters**

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

``_append_hist(self)``
''''''''''''''''''''''
Updates the history variables ``_r_hist``,``_it_hist``, by
appending the current density profile ``_r`` to ``_r_hist``
and appending ``_it_counter`` to ``_it_hist``.

``save_syst(self, path, filename)``
'''''''''''''''''''''''''''''''''''
Uses ``pickle.dump`` to save the instance variables of a
system.

**Parameters**

path : `String`
    Directory in which the system should be stored (needs to be
    a absolute path)
filename : `String`
    The filename under which the system should be stored.

*classmethod* ``load_syst(cls, path, filename)``
''''''''''''''''''''''''''''''''''''''''''''''''
Uses ``pickle.load`` to load a system. It is strongly
recommended to override this method in the inherited classes,
as the returned system might be of an outdated type! A typecast
should be implemented!

**Parameters**

path : `String`
    Directory in which the system is stored which one want's to
    load (needs to be a absolute path)
filename: `String`
    The filename under which the system of interest is stored.

**Returns**

Model : `LdftModel`
    The returned model probably has the type of an inherited
    class. It might also be the class of an outdated type.


``print_error(self)``
'''''''''''''''''''''
Returns a figure where the error history ``_err_hist`` is
plotted.

**Returns**

Figure : `matplotlib.pyplot.figure`
    Plotted error history.

``print_2d_profile(self)``
''''''''''''''''''''''''''
Creates a figure where the current profile is plotted. This
function is just for 2d-systems.

**Returns**

Figure : `matplotlib.pyplot.figure`
    Plotted profile

``print_2d_profile2(self)``
'''''''''''''''''''''''''''
Creates a figure where the current profile is plotted. This
function is just for 2d-systems.

**Returns**

Figure : `matplotlib.pyplot.figure`
    Plotted profile

``print_2d_hist(self, species=0, rows=10, idx_list=None)``
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Creates a figure where the history ``_r_hist`` is plotted.
Just one species can be plotted at the same time. Not the total
history is plotted but certain iteration steps.

**Parameters**

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

**Returns**

Figure : `matplotlib.pyplot.figure`
    Plotted history

``print_2d_hist2(self, species=0, rows=10, idx_list=None)``
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Creates a figure where the history ``_r_hist`` is plotted.
Just one species can be plotted at the same time. Not the total
history is plotted but certain iteration steps.

**Parameters**

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

**Returns**

Figure : `matplotlib.pyplot.figure`
    Plotted history


*abstractmethod* ``_cal_p(self, dens)``
'''''''''''''''''''''''''''''''''''''''
Calculates the pressure for a bulk system with given densities
for each species. The other parameters (temperature, attraction
strength, etc.) are taken from the current instance ``self``.

**Parameters**

dens : `List`
    The density for each species.

**Returns**

The pressure : `Float` 

*abstractmethod* ``_cal_coex_dens(self)``
'''''''''''''''''''''''''''''''''''''''''
Calculates the coexisting densities of bulk system for each
species under the parameters of the current instance ``self``.

**Returns**

Coexisting densities : `List` of `Tuple`
    The coexisting densities arranged in a List of Tuples. Each
    species corresponds to a Tuple of the form:
    (vapour_dens, liquid_dens)

``cal_p_vap(self)``
'''''''''''''''''''
Calculates the coexisting pressures under the current
parameters of the system (``_mu``, ``_dens``) and returns the
vapour pressure.

**Returns**

Vapour pressure : `Float`
    The vapour pressure of the current system

``cal_p_liq(self)``
'''''''''''''''''''
Calculates the coexisting pressures under the current
parameters of the system (``_mu``, ``_dens``) and returns the
liquid pressure.

**Returns**

Liquid pressure : `Float`
    The vapour pressure of the current system

``det_intface_shape(self)``
'''''''''''''''''''''''''''
Determines the shape of the interface of the current
configuration. It requires the inhomogeneities to be centered in
the system.

**Returns**

Shape : `String`
    The shape of the interface: 'Droplet', 'Cylinder', 'Slab',
    'Homogeneous'

``cal_del_Om(self)``
''''''''''''''''''''
Calculates the delta between the current grand potential and
the one by a homogeneous system of (oversaturated) vapor with the
same chemical potential as the current system.

**Returns**

delta Omega : `Float`
    Delta of the grand potential

``cal_R_s(self)``
'''''''''''''''''
Calculates the radius of surface of tension. In case of a
Cylinder configuration in three dimensions, the cylinder has to
point in the 0th axis of the density profile ``self._r``.

**Returns**

Radius of s.o.t. : `Float`
    Radius of surface of tension

``cal_R_em(self, em_species=0)``
''''''''''''''''''''''''''''''''
Calculates the equimolar radius for the species given by
``em_species``. In case of cylinder configurations in three
dimensions the cylinder has to point in the 0th axis of the
density profile ``self._r``. This function does only work
properly, if a droplet/cylinder is embedded in a supersaturated
vapour. For configurations of bubbles or vapour cylinders
embedded in liquid, the result will be wrong.

**Parameters**

em_species : `Int`; Optional: default=0
    Decides for which species the equimolar radius should
    be calculated

**Returns**

equimolar radius : `Float`
    Radius or the equimolar surface of a specific species.

``cal_gamma_R(self, R)``
''''''''''''''''''''''''
Calculates the surface tension for spheres/circles in 3d/2d of
radius ``R``. In case of cylinder configurations in three
dimensions the cylinder has to point in the 0th axis of the
density profile ``self._r``.

**Parameters**

R : `Float`
    Radius at which the surface tension should be calculated

**Returns**

surface tension : `Float`
    surface tension for radius R.

``cal_gamma_s(self)``
'''''''''''''''''''''
Calculates the surface tension for spheres/circles in 3d/2d at
the surface of tensions. In case of cylinder configurations in
three dimensions the cylinder has to point in the 0th axis of the
density profile ``self._r``.

**Returns**

surface tension : `Float`
    surface tension at the surface of tension.

``cal_gamma_em(self, species=0)``
'''''''''''''''''''''''''''''''''
Calculates the surface tension for spheres/circles in 3d/2d at
the equimolar surface of a given species. In case of cylinder
configurations in three dimensions the cylinder has to point in
the 0th axis of the density profile ``self._r``. This function
does only work properly, if a droplet/cylinder is embedded in a
supersaturated vapour. For configurations of bubbles or vapour
cylinders embedded in liquid, the result will be wrong.

**Parameters**

species : `Int`; Optional: default=0
    species for the equimolar surface

**Returns**

Surface tension : `Float`
    Surface tension at the equimolar surface.

``cal_adsorptionAtSurfOfTens(self, species=0)``
'''''''''''''''''''''''''''''''''''''''''''''''
Calculates the adsorption for spheres/circles in 3d/2d at the
surface of tension for a given species. This function
does only work properly, if a droplet/cylinder is embedded in a
supersaturated vapour. For configurations of bubbles or vapour
cylinders embedded in liquid, the result will be wrong.

**Parameters**

species : `Int`; Optional: default =0
    Species for which the adsorption should be calculated

**Returns**

Adsorption : `Tuple` of `Float`
    First entry: Adsorbed particle number; Second entry:
    adsorption.

``cal_gamma_inf(self, area)``
'''''''''''''''''''''''''''''
Calculates the surface tension of a flat interface. This
function can not determine the area of the surface itself.
Therefore it has to be passed as parameter.

**Parameters**

area : `float`
    Area of the surface. There are always two surfaces
    separating the liquid and the vapour. Meant is the area of
    one of those

**returns:**
    ``Tuple`` of ``Float``: First entry: Adsorbed particle
    number; Second entry: Adsorption.
