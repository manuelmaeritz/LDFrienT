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
_size : `None`
    Here the ``size``-parameter is stored. See it's description for
    further detail. (`Tuple`)
    
Methodes
--------
``__init__(self, size, mu_fix, mu=None, dens=None, v_ext=None, r=None, r_hist=None, err_hist=None, it_hist=None, bound_cond='periodic')``
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

----

``cal_F(self)``
Calculates the free energy of the models curent density
profile (meaning every species treated canonical, as if
``_mu_fix`` is ``False`` for every species)

**Returns**

The free energy : `Float`

----

``cal_Om(self)``
Calculates the grand potential of the models curent density
profile (meaning every species treated grand canonicaly, as if
``_mu_fix`` is ``True`` for every species).

**Returns**

The grand potential : `Float`
