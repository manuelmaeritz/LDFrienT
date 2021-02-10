Documentation on LG2dAOHigh
============================
*class* ``ldft_classes_v2.lg_2d_highl.LG2dAOHighl(self, size, epsi, mu_fix_c=False, mu_c=None, dens_c=None, v_ext_c=None, bound_cond='periodic', r=None, r_hist=None, err_hist=None, it_hist=None)`` 
Base: ``ldft_model.LdftModel`` 

This class describes a single component lattice gas in 2d with
sticky next neighbour attractions on a simple cubic lattice. The
description is done within the framework of lattice density
functional theory (ldft). The free energy functional was constructed
by translating the model to the Asakura-Oosawa (AO) model and then
setting up the functional of the resulting colloid-polymer
dispersion by the Highlander version of dft. Therefor this class
works with three species instead of one, namely the species of the
extended AO-model (colloid, polymer clusters species accounting for
attraction in x-direction and polymer for the attraction in
y-direction). The free energy functional is the one for the three
species. It differs from the free energy functional of the AO-model
just by a correction term accounting for the zero- and one-body
interaction of the polymers. If one wants the free energy of the
lattice gas, one would have to calculate the semi-grand potential of
the previous free energy, where the polymer clusters are treated
grand-canonically and the colloids canonically. In this class extra
functions are supported for this. The colloids correspond to the
species in the lattice gas.

Parameters
----------
size : `tuple` of `int`
    Shape of the systems simulation box. Expects a `Tuple` of two
    integers, each for one dimensional axis.
epsi : `float`
    Attraction strength of the lattice gas particles (multiplied
    with the inverse temperature to make it's dimension 1). From
    this the value of the chemical potential of the polymer clusters
    is calculated.
mu_fix_c : `bool`, optional: default = False
    Determines whether or not the system is treated canonical or
    grand canonical. Meant is the lattice gas system. This parameter
    therefore only steers the colloid-species. The others are set
    `True` by default. `False` for canonical.
mu_c : `float`, optional: default = `None`
    The chemical potential for the colloid species (multiplied with
    the inverse temperature to make it's dimension 1). Just required
    when ``mu_fix==True``. The chemical potential of the polymer
    clusters is determined by the value of ``epsi``.
dens_c : `float`, optional: default = `None`
    The average density of the colloids. Just required when
    ``mu_fix==False``. The average density of the polymer clusters
    is not required, as for those ``mu_fix`` is set `True`.
v_ext_c : `numpy.ndarray`, optional: default=`None`
    An external potential for the colloids. Shape must be of the
    same shape as chosen in ``size``. This class does not consider
    the possibility of sticky walls. Therefore the external
    potential of polymers is set zero by default.
bound_cond : `string`, optional: default='periodic'
    The boundary condition. Supports 'periodic' for periodic
    boundary conditions and '11_if' for a 45° tilted system with
    respect to the lattice. The latter is for creating slab
    interface with (11) orientation. If '11_if' is chosen then one
    dimension has to be chosen twice as the other dimension in the
    ``size`` parameter e.g. (64, 128). Default value is 'periodic'.
r : `List` of `np.array`; Optional: default = `None`
    Density profile for all three species arranged in a `List`. Choose
    `None` in case you hand over the ``r_hist``-parameter or in case
    you do not want to set the variable yet.
r_hist : `List` of `List` of `np.array`; Optional: default = `None`
    Picard-history of a density profile. It contains the density
    profiles for certain picard-steps of a system which has already
    been evolved through picard iteration. Caution! Every entry is
    of the format of the ``_r``-instance variable, which is a list
    itself containing the profile for each species. Therefore in our
    case the list is of length one. Use `None` if the system has no
    history yet.
err_hist : `List` of `Tuple` of `Float`; Optional: default = `None`
    Contains the error at the picard-steps corresponding to the
    entries of `r_hist`. The entries are tuples containing an error
    for every species. Use `None` if no history available.
it_hist : `List`; Optional: default = `None`
    List of the picard-steps corresponding to the density profiles at
    the ``r_hist``-parameter. Use `None` if no history available.
    Note: if ``r_hist`` is given then also this argument should be
    assigned with an appropriate list.

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
Properties for external use
^^^^^^^^^^^^^^^^^^^^^^^^^^^
size : `Tuple`, read-only
    Accesses the ``_size``-attribute 

mu : `List`, read-only
    Accesses the ``_mu``-attribute

dens : `List`, read-only
    Accesses the ``_dens``-attribute

mu_fix : `List`, read-only
    Accesses the ``_mu_fix``-attribute

v_ext : `List`, read-only
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

epsi(self) : `Float`, read-only
    The attraction strength between the lattice-particles of the
    lattice gas.

mu_c(self) : `Float`, read and write
    The chemical potential of the colloids (times the inverse
    temperature to make its dimension 1). It is equals the chemical
    potential of the particles of the lattice gas.

mu_pc1(self) : `Float`, read-only
    The chemical potential of the polymer-cluster in x-direction
    (times the inverse temperature to make its dimension 1).

mu_pc2(self) : `Float`, read-only
    The chemical potential of the polymer-cluster in y-direction
    (times the inverse temperature to make its dimension 1).

dens_c(self) : `float`, read-only
    """The average density of the colloids. It is equals the average
    density in the lattice gas. 

dens_pc1(self) : `floag`, read-only
    The average density of the polymer clusters in x-direction.

dens_pc2(self) : `floag`, read-only
    The average density of the polymer clusters in y-direction.

mu_fix_c(self) : `Bool`, read and write
    Flag which determines Wether the colloids (a.k. the particles
    of the lattice gas) are treated canonical (`False`) or grand
    canonical (`True`). (`Bool`)

v_ext_c(self) : `np.array`, read and write
    External potential acting on the colloids (a.k. the particles
    of the lattice gas).

r_c(self) : `np.array`, read-only
    The density profile of the colloids (a.k. the particles of
    the lattice gas).

r_pc1(self) : `np.array`, read-only
    The density profile of the polymer clusters in x-direction.

r_pc2(self) : `np.array`, read-only
    The density profile of the polymer clusters in y-direction.

r_c_hist(self) : `List`, read-only
    Iteration history of the density profile of the colloids
    (a.k. the particles of the lattice gas).

err_c_hist(self) : `List`, read-only
    Iteration history of the picard-error at the colloidal
    density profile.

Properties for internal use
^^^^^^^^^^^^^^^^^^^^^^^^^^^

_mu_c(self) : `float`, read-only
    The chemical potential of the colloid species (times the
    inverse temperature to make its dimension 1)

_mu_pc1(self) : `float`, read-only
    The chemical potential of the polymer species in x-direction
    (times the inverse temperature to make its dimension 1).

_mu_pc2(self) : `float`, read-only
    The chemical potential of the polymer species in y-direction
    (times the inverse temperature to make its dimension 1).

_dens_c(self) : `float`, read-only
    The average density of the colloid species 

_dens_pc1(self) : `float`, read-only
    The average density of the polymer species in x-direction

_dens_pc2(self) : `folat`, read-only
    The average density of the polymer species in x-direction

_v_ext_c(self) : `np.array`, read-only
    The external potential acting on the colloids

_v_ext_pc1(self) : `np.array`, read-only
    The external potential acting on the polymer clusters in
    x-direction.

_v_ext_pc2(self) : `np.array`, read-only
    The external potential acting on the polymer clusters in
    y-direction.

_r_c(self) : `numpy.ndarray`, read-only
    The density profile of the colloid species.

_r_pc1(self) : `numpy.ndarray`, read-only
    The density profile of the polymer species in x-direction.

_r_pc2(self) : `numpy.ndarray`, read-only
    The density profile of the polymer species in y-direction.

Methods
--------
``__init__(self, size, epsi, mu_fix_c=False, mu_c=None, dens_c=None, v_ext_c=None, bound_cond='periodic', r=None, r_hist=None, err_hist=None, it_hist=None)``
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

**Parameters**

size : `tuple` of `int`
    Shape of the systems simulation box. Expects a `Tuple` of two
    integers, each for one dimensional axis.
epsi : `float`
    Attraction strength of the lattice gas particles (multiplied
    with the inverse temperature to make it's dimension 1). From
    this the value of the chemical potential of the polymer clusters
    is calculated.
mu_fix_c : `bool`, optional: default = False
    Determines whether or not the system is treated canonical or
    grand canonical. Meant is the lattice gas system. This parameter
    therefore only steers the colloid-species. The others are set
    `True` by default. `False` for canonical.
mu_c : `float`, optional: default = `None`
    The chemical potential for the colloid species (multiplied with
    the inverse temperature to make it's dimension 1). Just required
    when ``mu_fix==True``. The chemical potential of the polymer
    clusters is determined by the value of ``epsi``.
dens_c : `float`, optional: default = `None`
    The average density of the colloids. Just required when
    ``mu_fix==False``. The average density of the polymer clusters
    is not required, as for those ``mu_fix`` is set `True`.
v_ext_c : `numpy.ndarray`, optional: default=`None`
    An external potential for the colloids. Shape must be of the
    same shape as chosen in ``size``. This class does not consider
    the possibility of sticky walls. Therefore the external
    potential of polymers is set zero by default.
bound_cond : `string`, optional: default='periodic'
    The boundary condition. Supports 'periodic' for periodic
    boundary conditions and '11_if' for a 45° tilted system with
    respect to the lattice. The latter is for creating slab
    interface with (11) orientation. If '11_if' is chosen then one
    dimension has to be chosen twice as the other dimension in the
    ``size`` parameter e.g. (64, 128). Default value is 'periodic'.
r : `List` of `np.array`; Optional: default = `None`
    Density profile for all three species arranged in a `List`. Choose
    `None` in case you hand over the ``r_hist``-parameter or in case
    you do not want to set the variable yet.
r_hist : `List` of `List` of `np.array`; Optional: default = `None`
    Picard-history of a density profile. It contains the density
    profiles for certain picard-steps of a system which has already
    been evolved through picard iteration. Caution! Every entry is
    of the format of the ``_r``-instance variable, which is a list
    itself containing the profile for each species. Therefore in our
    case the list is of length one. Use `None` if the system has no
    history yet.
err_hist : `List` of `Tuple` of `Float`; Optional: default = `None`
    Contains the error at the picard-steps corresponding to the
    entries of `r_hist`. The entries are tuples containing an error
    for every species. Use `None` if no history available.
it_hist : `List`; Optional: default = `None`
    List of the picard-steps corresponding to the density profiles at
    the ``r_hist``-parameter. Use `None` if no history available.
    Note: if ``r_hist`` is given then also this argument should be
    assigned with an appropriate list.

``__str__(self)``
'''''''''''''''''

Translating the lattice gas to the AO-model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*staticmethod* ``translate_epsi_to_mu_pc(epsi)``
'''''''''''''''''''''''''''''''''
Maps the attraction strength of the lattice gas ``epsi`` to
the corresponding polymer cluster chemical potential.

**Parameters**

epsi : `float`
    The attraction strength (multiplied with the inverse
    temperature to make the quantity dimensionless).

**Returnsi**

mu_pc : The chemical potential (multiplied with the inverse
    temperature to make the quantity dimensionless). (`float`)

*staticmethod* ``translate_mu_pc_to_epsi(mu_pc)``
'''''''''''''''''''''''''''''''''''''''''''''''''
Maps the polymer cluster chemical potential to the attraction
strength of the lattice gas ``epsi``.

**Parameters**

mu_pc : `float`
    The polymer chemical potential (multiplied with the inverse
    temperature to make the quantity dimensionless).

**Returns**

epsi : The attraction strength (multiplied with the inverse
    temperature to make the quantity dimensionless). (`float`)

Methods concerning the functional
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``_cal_n(self)``
''''''''''''''''
Calculates the weighted densities necessary for the
calculation of the free energy and the excess chemical
potential.

**Returns**

Result : `tuple` of `numpy.ndaray`

``_cal_Phi_ex_AO(self)``
''''''''''''''''''''''''
Calculates the excess free energy of the AO-model.

**Returns**

Result : `np.array`
    Free energy density of the AO-model.

``cal_F(self)``
'''''''''''''''
Calculates the free energy of the three component system. It
differs from the free energy functional of the AO-model just by
a correction term accounting for the zero- and one-body
interaction of the polymers (see description of the class). For
getting the free energy of the lattice gas use ``cal_F_lg``,
which is the semi-grand potential, where the polymer clusters are
treated grand canonically and the colloids canonically.

**Returns**

Result : `float`
    Free energy of the three component system (times the inverse
    temperature to make the results dimension 1).

``cal_F_lg(self)``
''''''''''''''''''
Calculates the free energy of the lattice gas. If
``self.mu_fix==False`` this should give the same result as the
``cal_semi_Om``-function.

**Returns**

Result : `float`
    Free energy of the lattice gas.

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

``cal_mu_ex(self)``
'''''''''''''''''''
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


Methods concerning the bulk properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*classmethod* ``_cal_bulk_r_pc(cls, r_c, epsi)``
''''''''''''''''''''''''''''''''''''''''''''''''
Calculates the bulk polymer cluster density in dependence of
the colloid density and the chosen attraction strength

**Parameters**

r_c : `float` or `np.ndarray`
    The colloid density.
epsi : `float`
    Attraction strength (times inverse temperature).

**Returns**

r_pc : `float` or `np.ndarray`
    The polymer cluster density.

*classmethod* ``_cal_bulk_dr_pc(cls, r_c, epsi)``
'''''''''''''''''''''''''''''''''''''''''''''''''
Calculates the derivative of the bulk polymer cluster density
with respect to the colloidal density in dependence of
the colloid density and the chosen attraction strength

**Parameters**

r_c : `float` or `np.ndarray`
    The colloid density.
epsi : `float`
    Attraction strength (times inverse temperature).

Returns
-------
dr_pc : `float` or `np.ndarray`
    The derivative of the polymer cluster density.

*classmethod* ``cal_bulk_mu_lg(cls, r_c, epsi)``
''''''''''''''''''''''''''''''''''''''''''''''''
Calculates the chemical potential for a bulk lattice gas.

**Parameters**

r_c : `Float` or `np.ndarray`
    The colloidal density.
epsi : `Float`
    Attraction strength

**Returns**

mu : `Float` or `np.ndarray`
    The chemical potential for the lattice gas. 

*classmethod* ``cal_bulk_dmu_lg(cls, r_c, epsi)``
'''''''''''''''''''''''''''''''''''''''''''''''''
Calculates the derivative of the chemical potential from the
bulk lattice gas with respect to the colloidal density.

**Parameters**

r_c : `Float` or `np.ndarray`
    The colloidal density.
epsi : `Float`
    Attraction strength

**Returns**

dmu : `Float` or `np.ndarray`
    The derivative of the chemical potential from the lattice
    gas.

*classmethod* ``_cal_bulk_f_AO_id(cls, r_c, r_pc)``
'''''''''''''''''''''''''''''''''''''''''''''''''''
Calculates the ideal gas part of the free energy density of
a bulk AO-system under given colloid and polymer cluster
density.

**Parameters**

r_c : `float`
    Colloid density
r_pc : `float`
    Polymer cluster density

**Returns**

f_id : `float`
    The idea gas part of the free energy density.

*classmethod* ``_cal_bulk_f_AO_ex(cls, r_c, r_pc)``
'''''''''''''''''''''''''''''''''''''''''''''''''''
Calculates the excess part of the free energy density of a
bulk AO-system under given colloid and polymer cluster density.

**Parameters**

r_c : `float`
    Colloid density
r_pc : `float`
    Polymer cluster density

**Returns**

f_ex : `float`
    The excess part of the free energy density.

*classmethod* ``cal_bulk_f_lg(cls, r_c, epsi)``
'''''''''''''''''''''''''''''''''''''''''''''''
Calculates the free energy density of the bulk lattice gas
under given density. (The function is the same as in
``cal_F_lg`` but simplified for bulk systems.)

**Parameters**

r_c: `float` or `np.ndarray`
    Density
epsi: `float`
    Attraction strength (times inverse temperature)

**Returns**

f : `float` or `np.ndarray`
    The free energy density of a bulk lattice gas.

*classmethod* ``cal_bulk_om_lg(cls, r, epsi)``
''''''''''''''''''''''''''''''''''''''''''''''
Calculates the grand potential density for a bulk lattice gas
under given densities.

**Parameters**

r : `float` or `np.ndarray`
    The density.
epsi : `float`
    The attraction strength (times inverse temperature).

**Returns**

om : `Float` 
    The grand potential density

*classmethod* ``cal_bulk_p(cls, r, epsi)``
''''''''''''''''''''''''''''''''''''''''''
Calculates the pressure of a bulk lattice gas under given
density.

**Parameters**

r : `float` or `np.ndarray`
    The density.
epsi : `float`
    The attraction strength (times inverse temperature).

**Returns**

The pressure : `Float`

*classmethod* ``_cal_difMu(cls, r_c, *args)``
'''''''''''''''''''''''''''''''''''''''''''''
Calculates the difference between a certain chemical
potential of the lattice gas and the chemical potential
belonging to a certain density. This is a help-function for the
function ``cal_bulk_coex_dens``.

**Parameters**

r_c : `float`
    The colloid density of the system
*args:
    First argument: Attraction strength (times inverse
    temperature). (`float`)
    Second argument: The reference chemical potential which the
    chemical potential for at density r_c should be compared to.
    (`float`)

**Returns**

difMu : `float`
    The difference between the two colloidal chem. pot.

*classmethod* ``cal_bulk_coex_dens(cls, mu, epsi, init_min=0.01, init_max=0.99)``
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Calculates the coexisting densities of a bulk system lattice
gas system under given chemical potential.

**Parameters**

mu : `Float`
    The chemical potential of the lattice gas.
epsi : `Float`
    The attraction strength (times inverse temperature).

**Returns**

r_coex : `Tuple`
    The coexisting densities arranged in a tuple of the shape 
    (vapour_dens, liquid_dens)
