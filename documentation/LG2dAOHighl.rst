Documentation on LG2dAOHighl
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

Methodes
--------

