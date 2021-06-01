from .ldft_model import LdftModel
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
from functools import reduce

class LG3dMf(LdftModel):
    """This class describes a single component lattice gas in 3d with
    sticky next neighbour attractions on a simple cubic lattice. The
    description is done within the framework of mean field theory. The
    class inherits from the interface-class ``LdftModel``.
    
    Parameters
    ----------
    size : `tuple` of `int`
        Shape of the systems simulation box. Expects a `Tuple` of three
        integers, each for one dimensional axis.
    epsi : `float`
        Attraction strength (multiplied with the inverse temperature to
        make it's dimension 1)
    mu_fix : `bool`, optional: default = False
        Determines whether or not the system is treated canonical or
        grand canonical. ``False`` for canonical.
    mu : `float`, optional: default = `None`
        The chemical potential (multiplied with the inverse temperature
        to make it's dimension 1). Just required when ``mu_fix==True``.
    dens : `float`, optional: default = `None`
        The average density of the system. Just required when
        ``mu_fix``==`False`.
    v_ext : `numpy.ndarray`, optional: default=`None`
        An external potential. Shape must be of the same shape as
        chosen in ``size``.
    bound_cond : `string`, optional: default='periodic'
        The boundary condition. Supports 'periodic' for periodic
        boundary conditions and '11_if' for a 45° tilted system with
        respect to the lattice. The latter is for creating slab
        interface with (11) orientation. If '11_if' is chosen then one
        dimension has to be chosen twice as the other dimension in the
        ``size`` parameter e.g. (64, 128). Default value is 'periodic'.
    r : `np.array`; Optional: default = `None`
        Density profile. choose `None` in case you hand over the
        ``r_hist``-parameter or in case you do not want to set the
        variable yet.
    r_hist : `List`; Optional: default = `None`
        Picard-history of a density profile. It contains the density
        profiles for certain picard-steps of a system which has already
        been evolved through picard iteration. Caution! Every entry is
        of the format of the ``_r``-instance variable, which is a list
        itself containing the profile for each species. Therefore in our
        case the list is of length one. Use `None` if the system has no
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
    """

    _epsi = None
    """Here the ``epsi``-parameter is stored. See its description for
    further detail.
    """

    def __init__(self, size, epsi, mu_fix=False, mu=None, dens=None,
            v_ext=None, bound_cond='periodic', r=None, r_hist=None,
            err_hist=None, it_hist=None):
        self._epsi=epsi
        mu_fix = [mu_fix]
        mu = [mu] if mu!=None else None
        dens = [dens] if dens!=None else None
        v_ext= [v_ext] if type(v_ext)==type(np.array([])) else None
        r = [r] if type(r)==type(np.array([])) else None
        super().__init__(size, mu_fix=mu_fix, mu=mu,
                dens=dens, v_ext=v_ext, bound_cond=bound_cond, r=r,\
                        r_hist=r_hist, err_hist=err_hist,\
                        it_hist=it_hist)

    def __str__(self):
        descrLG2MF = 'This is a Lattice gas described with mean field'\
                +'DFT. It is an object of the Type \'LG3dMF\' and has'\
                +'the following properties:'
        epsiStr='{0:<40s}: {1}\n'.format('Attr. strength \'epsi\'',\
                self._epsi)
        descrLdftModel = 'It inherits from \'LdftModel\', with the'\
                +'following properties:'
        restStr=super().__str__()
        return descrLG2MF+'\n\n'+epsiStr+'\n'+descrLdftModel+\
                '\n'+restStr

    ####################################################################
    #In the following some properties are overwritten. As there is just
    #one species around in this model the lists summarizing a certain
    #property for different species are not necessary here. Therefore it
    #is more convenient to work without the lists of species. The
    #following properties account for this. But one must not overwrite
    #the instance variables of the LdftModel class. Internally the code
    #still works with the lists. Overwriting of instance variables might
    #cause errors.
    ####################################################################

    @property
    def epsi(self):
        """See description of ``_epsi`` (`Tuple`, read-only)
        """
        return self._epsi

    @LdftModel.mu.getter
    def mu(self):
        """Chemical potential. It is the (only) entry of the list
        ``_mu``. See its description for further details. (`Float`)
        """
        return self._mu[0]

    @mu.setter
    def mu(self, mu):
        self._mu[0] = mu

    @LdftModel.dens.getter
    def dens(self):
        """Average density of the System. It is the (only) entry of the
        list ``_dens``. See its description for further details. (`Float`)
        """
        return self._dens[0]

    @dens.setter
    def dens(self, dens):
        self._dens[0] = dens

    @LdftModel.mu_fix.getter
    def mu_fix(self):
        """Flag indicating whether the system is treated canonical or
        grand canonical. It is the (only) entry of the list ``_mu_fix``.
        See its description for further details. (`Bool`)
        """
        return self._mu_fix[0]

    @mu_fix.setter
    def mu_fix(self, mu_fix):
        self._mu_fix[0] = mu_fix

    @LdftModel.v_ext.getter
    def v_ext(self):
        """External potential. It is the (only) entry of the list
        ``_v_ext``. See its description for further details.
        (`numpy.ndarray`)
        """
        return self._v_ext[0]

    @v_ext.setter
    def v_ext(self, v_ext):
        self._v_ext[0] = v_ext

    @LdftModel.r.getter
    def r(self):
        """Current density profile. It is the (only) entry of the list
        ``_r``. See its description for further details. (`np.array`)
        """
        return self._r[0]

    @r.setter
    def r(self, r):
        self.set_r([r])

    @LdftModel.r_hist.getter
    def r_hist(self):
        """Iteration history of the density profile. It is the same list
        as ``_r_hist`` but with the density profiles not wrapped in
        another list. (`List`, read-only)
        """
        r_hist = [r[0] for r in self._r_hist]
        return r_hist

    @LdftModel.err_hist.getter
    def err_hist(self):
        """Iteration history of the picard-error. It is the same list as
        ``_err_hist`` but with the density profiles not wrapped in
        another list. (`List`, read-only)
        """
        err_hist =[err[0] for err in self._err_hist]
        return err_hist

    ####################################################################
    #In the following we define all the models necessary for overwriting
    #the abstract method ``cal_F`` and ``cal_mu_ex`` and eventually
    #overwrite them.
    ####################################################################

    def cal_F_id(self):
        """Calculates the ideal part of the free energy on the current
        density profile ``_r``.

        Returns
        -------
        Ideal free energy : `Float`
        """
        r = self._r[0]
        F_id = np.sum(r*(np.log(r)-1))
        return F_id

    def cal_F_hr(self):
        """Calculates the hard rod part of the excess free energy on the
        current density profile ``_r``.

        Returns
        -------
        Hard rod part of the excess free energy : `Float`
        """
        r = self._r[0]
        F_hr = np.sum(r+(1-r)*np.log(1-r))
        return F_hr

    def cal_F_sa(self):
        """Calculates the sticky attraction part of the excess free
        energy on the current density profile ``_r``.

        Returns
        -------
        Sticky attraction part of the excess free energy : `Float`
        """
        r = self._r[0]
        epsi = self._epsi
        r1 = self._boundary_roll(r, -1, axis=1)
        r2 = self._boundary_roll(r, 1, axis=1)
        r3 = self._boundary_roll(r, -1, axis=0)
        r4 = self._boundary_roll(r, 1, axis=0)
        r5 = self._boundary_roll(r, -1, axis=2)
        r6 = self._boundary_roll(r, 1, axis=2)
        nn = r1+r2+r3+r4+r5+r6
        F_sa = -0.5*epsi*np.sum(r*nn)
        return F_sa

    def cal_F(self):
        F = self.cal_F_id()+self.cal_F_hr()+self.cal_F_sa()
        return F

    @LdftModel._RespectBoundaryCondition()
    def cal_mu_ex(self):
        r = self._r[0]
        epsi = self._epsi
        r1 = self._boundary_roll(r, -1, axis=1)
        r2 = self._boundary_roll(r, 1, axis=1)
        r3 = self._boundary_roll(r, -1, axis=0)
        r4 = self._boundary_roll(r, 1, axis=0)
        r5 = self._boundary_roll(r, -1, axis=2)
        r6 = self._boundary_roll(r, 1, axis=2)
        mu_ex = np.log(1-r)+epsi*(r1+r2+r3+r4+r5+r6)
        return [mu_ex]

    ####################################################################
    #In the following we will provide some methods for bulk properties
    ####################################################################

    @classmethod
    def cal_bulk_f(cls, dens, epsi):
        '''Calculates the free energy density of a bulk system under
        given density.

        Parameters
        ----------
        dens : `float`
            Density
        epsi : `float`
            Attraction strength (times inverse temperature)

        Returns
        -------
        f : `float`
            The free energy density
        '''
        rho = dens[0] if type(dens)==list else dens
        f_id = rho*(np.log(rho)-1)
        f_hr = rho + (1-rho)*np.log(1-rho)
        f_sa = -3*epsi*rho**2
        f = f_id+f_hr+f_sa
        return f

    @classmethod
    def cal_bulk_mu(cls, dens, epsi):
        '''Calculates the chemical potential of a bulk system under
        given density.

        Parameters
        ----------
        dens : `float`
            Density
        epsi : `float`
            Attraction strength (times inverse temperature)

        Returns
        -------
        mu : `float`
            The chemical potential.

        '''
        rho = dens[0] if type(dens)==list else dens
        mu = np.log(rho/(1-rho))-epsi*6*rho
        return mu

    @classmethod
    def cal_bulk_om(cls, r, epsi):
        """Calculates the grand potential density for a bulk lattice gas
        under given densities.

        Parameters
        ----------
        r : `float` or `np.ndarray`
            The density.
        epsi : `float`
            The attraction strength (times inverse temperature).

        Returns
        -------
        om : `Float`
            The grand potential density
        """
        f = cls.cal_bulk_f(r, epsi)
        mu = cls.cal_bulk_mu(r, epsi)
        om = f-mu*r
        return om

    @classmethod
    def cal_bulk_p(cls, r, epsi):
        """Calculates the pressure of a bulk lattice gas under given
        density.

        Parameters
        ----------
        r : `float` or `np.ndarray`
            The density.
        epsi : `float`
            The attraction strength (times inverse temperature).

        Returns
        -------
        The pressure : `Float`
        """
        p = -cls.cal_bulk_om(r, epsi)
        return p

    @classmethod
    def _cal_difMu(cls, rho, *args):
        """Calculates the difference between a certain chemical
        potential and the chemical potential belonging to a certain
        density. This is a help-function for the function
        ``cal_bulk_coex_dens``.

        Parameters
        ----------
        rho : `float`
            Density for which the chemical potential should be compared
            with. Number between 0 and 1.
        *args:
            First argument: Attraction strength (times inverse
            temperature). (`float`)
            Second argument: The reference chemical potential which the
            chemical potential for at density r_c should be compared to.
            (`float`)

        Returns
        -------
        difMu : `float`
            The difference between the two chemical potentials
        """
        epsi = args[0]
        mu = args[1]
        mu_rho = cls.cal_bulk_mu(rho, epsi)
        return mu_rho-mu

    @classmethod
    def cal_bulk_coex_dens(cls, mu, epsi):
        """Calculates the coexisting densities of a bulk system lattice
        gas system under given chemical potential.

        Parameters
        ----------
        mu : `Float`
            The chemical potential of the lattice gas.
        epsi : `Float`
            The attraction strength (times inverse temperature).

        Returns
        -------
        r_coex : `Tuple`
            The coexisting densities arranged in a tuple of the shape
            (vapour_dens, liquid_dens)
        """
        r_coex = op.fsolve(cls._cal_difMu,\
                np.array([0.01, 0.5, 0.99]), args=(epsi, mu))
        r_coex = (r_coex[0], r_coex[2])
        return r_coex

    ####################################################################
    #In the following section the abc-methods concerning the surface
    #properties of the mother class are overridden
    ####################################################################

    def _cal_p(self, dens):
        epsi = self._epsi
        r = dens[0]
        p = self.cal_bulk_p(r, epsi)
        return p

    def _cal_coex_dens(self):
        mu = self._mu[0]
        epsi = self._epsi
        r_c_coex = self.cal_bulk_coex_dens(mu, epsi)
        return [r_c_coex]







