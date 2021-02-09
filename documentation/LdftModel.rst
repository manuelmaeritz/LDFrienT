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
