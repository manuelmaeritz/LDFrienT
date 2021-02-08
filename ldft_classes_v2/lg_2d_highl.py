from .ldft_model import LdftModel
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
from functools import reduce

class LG2dAOHighl(LdftModel):
    """This class describes a single component lattice gas in 2d with
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
        therefore only steres the colloid-species. The others are set
        `True` by defalut. `False` for canonical.
    mu_c : `float`, optional: default = `None`
        The chemical potential for the colloid species (multiplied with
        the inverse temperature to make it's dimension 1). Just required
        when ``mu_fix==True``. The chemical potential of the polymer
        clusters is determined by the value of ``epsi``.
    dens_c : `float`, optional: default = `None`
        The average density of the colloids. Just required when
        ``mu_fix``==`False`. The average density of the polymer clusters
        is not required, as for those ``mu_fix`` is set `True`.
    v_ext_c : `numpy.ndarray`, optional: default=`None`
        An external potential for the colloids. Shape must be of the
        same shape as choosen in ``size``. This class does not consider
        the possibility of sticky walls. Therefore the external
        potential of polymers is set zero by default.
    bound_cond : `string`, optional: default='periodic'
        The boundary condition. Supports 'periodic' for periodic
        boundary conditions and '11_if' for a 45° tilted system with
        respect to the lattice. The latter is for creating slab
        interface with (11) orientation. If '11_if' is choosen then one
        dimension has to be choosen twice as the other dimension in the
        ``size`` parameter e.g. (64, 128). Default value is 'periodic'.
    r : `List` of `np.array`; Optional: default = `None`
        Desity profile for all thre species arranged in a `List`. Choose
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
        for every species. Use `None` if no history availabe.
    it_hist : `List`; Optional: default = `None`
        List of the picardsteps corresponding to the density profiles at
        the ``r_hist``-parameter. Use `None` if no history availabe.
        Note: if ``r_hist`` is given then also this argument should be
        assigned with an appropriate list.
    """
    
    def __init__(self, size, epsi, mu_fix_c=False, mu_c=None,\
            dens_c=None, v_ext_c=None, bound_cond='periodic', r=None,\
            r_hist=None, err_hist=None, it_hist=None):
        mu_pc=self.translate_epsi_to_mu_pc(epsi)
        v_ext_pc = np.zeros(size)
        v_ext_c = v_ext_pc if type(v_ext_c)==type(None) else v_ext_c
        super().__init__(size=size, mu_fix=[mu_fix_c, True, True],
                mu=[mu_c, mu_pc, mu_pc], dens=[dens_c, None, None],
                v_ext=[v_ext_c, v_ext_pc, v_ext_pc], r=r, r_hist=r_hist,
                err_hist=err_hist, it_hist=it_hist,
                bound_cond=bound_cond)

    def __str__(self):
        descrLG2dHighl = 'This is a Lattice gas described with lattice'\
                +' DFT. It was translated to the AO-model and the'\
                +' functional was constructed by the Highlander method'\
                +' It is an object of the Type \'LG2dAOHighl\' and has'\
                +' the following properties:'
        epsiStr='{0:<40s}: {1}\n'.format('Attr. strength \'epsi\'',\
                self.epsi)
        motherClass = 'It inherits from \'LdftModel\', with the'\
                +' following properties:'
        descrLdftModel=super().__str__()
        return descrLG2dHighl+'\n\n'+epsiStr+'\n'+motherClass+\
                '\n\n'+descrLdftModel

    ####################################################################
    #Protected descriptors for internal use. These are for a more
    #convenient adresseing of the species specific instance variables.
    #Important to notice: do not override the protected variables of the
    #super class LdftModel. Otherwise the functionality of the instance
    #methodes in LdftModel can not be secured.
    ####################################################################

    @property
    def _mu_c(self):
        """The chemical potential of the coloid species (times the
        inverse temperatue to make its dimension 1)
        (`float`, read-only).
        """
        return self._mu[0]

    @property
    def _mu_pc1(self):
        """The chemical potential of the polymer species in x-direction
        (times the inverse temperatue to make its dimension 1).
        (`float`, read-only)
        """
        return self._mu[1]

    @property
    def _mu_pc2(self):
        """The chemical potential of the polymer species in y-direction
        (times the inverse temperatue to make its dimension 1).
        (`float`, read-only)
        """
        return self._mu[2]

    @property
    def _dens_c(self):
        """The average density of the coloid species (`float`,
        read-only).
        """
        return self._dens[0]

    @property
    def _dens_pc1(self):
        """The average density of the polymer species in x-direction
        (`float`, read-only).
        """
        return self._dens[1]

    @property
    def _dens_pc2(self):
        """The average density of the polymer species in x-direction
        (`float`, read-only).
        """
        return self._dens[2]
    
    @property
    def _v_ext_c(self):
        """The external potential acting on the colloids (`np.array`,
        read-only)
        """
        return self._v_ext[0]

    @property
    def _v_ext_pc1(self):
        """The external potential acting on the polymer clusters in
        x-direction. (`np.array`, read-only)
        """
        return self._v_ext[1]
    
    @property
    def _v_ext_pc2(self):
        """The external potential acting on the polymer clusters in
        y-direction. (`np.array`, read-only)
        """
        return self._v_ext[2]

    @property
    def _r_c(self):
        """The density profile of the coloid species. (`numpy.ndarray`,
        read-only)
        """
        return self._r[0]

    @property
    def _r_pc1(self):
        """The density profile of the polymer species in x-direction.
        (`numpy.ndarray`, read-only)
        """
        return self._r[1]

    @property
    def _r_pc2(self):
        """The density profile of the polymer species in y-direction.
        (`numpy.ndarray`, read-only)
        """
        return self._r[2]


    ####################################################################
    #Public descriptors. These are for the user to access the variables
    #of intrest. Some are already defined in the super class. Some of
    #them are reused, but others are overwritten.
    ####################################################################

    @property
    def epsi(self):
        """The attraction strengs between the lattice-particles of the
        lattice gas. (`Float`, read-only)
        """
        return self.translate_mu_pc_to_epsi(self._mu_pc1)
    
    @property
    def mu_c(self):
        """The chemical potential of the colloids (times the inverse
        temperatue to make its dimension 1). It is equals the chemical
        potential of the particles of the lattice gas. (`float`)
        """
        return self._mu[0]

    @mu_c.setter
    def mu_c(self, mu_c):
        self._mu[0]=mu_c

    mu_pc1=_mu_pc1
    """The chemical potential of the polymer-cluster in x-direction
    (times the inverse temperature to make its dimension 1).
    (`float`, read-only)
    """

    mu_pc2=_mu_pc2
    """The chemical potential of the polymer-cluster in y-direction
    (times the inverse temperature to make its dimension 1).
    (`float`, read-only)
    """

    @LdftModel.mu.setter
    def mu(self, mu):
        print('This setter has been deactivated in favour for \`mu_c\`')

    @property
    def dens_c(self):
        """The average density of the colloids. It is equals the average
        density in the lattice gas. (`float`)
        """
        return self._dens[0]
    
    dens_pc1=_dens_pc1
    """The average density of the polymer clusters in x-direction.
    (`float`, read-only)
    """
    
    dens_pc2=_dens_pc2
    """The average density of the polymer clusters in x-direction.
    (`float`, read-only)
    """

    @LdftModel.dens.setter
    def dens(self, dens):
        print('This setter has been deactivated in favour for \
                \`dens_c\`')
    
    @property
    def mu_fix_c(self):
        """Flag which determines wether the colloids (a.k. the particles
        of the lattice gas) are treated canonical (`False`) or grand
        canonical (`True`). (`Bool`)
        """
        return self._mu_fix[0]

    @mu_fix_c.setter
    def mu_fix_c(self, mu_fix_c):
        self._mu_fix[0]=mu_fix_c
    
    @LdftModel.mu_fix.setter
    def mu_fix(self, mu_fix):
        print('This setter has been deactivated in favour for \
                \`mu_fix_c\`')

    @property
    def v_ext_c(self):
        """External potential acting on the colloids (a.k. the particles
        of the lattice gas). (`np.array`)
        """
        return self._v_ext[0]

    @v_ext_c.setter
    def v_ext_c(self, v_ext_c):
        self._v_ext[0]=v_ext_c

    @LdftModel.v_ext.setter
    def v_ext(self, v_ext):
        print('This setter has been deactivated in favour for \
                \`v_ext_c\`')

    @property
    def r_c(self):
        """The density profile of the colloids (a.k. the particles of
        the lattice gas). (`np.array`, read-only)
        """
        return self._r[0]

    r_pc1=_r_pc1
    """The density profile of the polymer clusters in x-direction.
    (`np.array`, read-only)
    """

    r_pc2=_r_pc2
    """The density profile of the polymer clusters in y-direction.
    (`np.array`, read-only)
    """
    
    @property
    def r_c_hist(self):
        """Iteration history of the density profile of the colloids
        (a.k. the particles of the lattice gas). (`List`, read-only)
        """
        r_c_hist = [r[0] for r in self._r_hist]
        return r_c_hist

    @property
    def err_c_hist(self):
        """Iteration history of the picard-error at the colloidal
        density profile. (`List`, read-only)
        """
        err_hist =[err[0] for err in self._err_hist]
        return err_hist

    ####################################################################
    #Map the lattice gas to the AO-model:
    ####################################################################
    
    @staticmethod
    def translate_epsi_to_mu_pc(epsi):
        """Maps the attraction strength of the lattice gas ``epsi`` to
        the corresponding polymercluster chemical potential.

        Parameters
        ----------
        epsi : `float`
            The attraction strength (multiplied with the inverse
            temperature to make the quantitiy dimensionless).

        Returns
        -------
        mu_pc : The chemical potential (multiplied with the inverse
            temperature to make the quantitiy dimensionless). (`float`)

        """
        mu_pc=np.log(np.exp(epsi)-1)
        return mu_pc

    @staticmethod
    def translate_mu_pc_to_epsi(mu_pc):
        """Maps the polymercluster chemical potential to the attraction
        strength of the lattice gas ``epsi``.

        Parameters
        ----------
        mu_pc : `float`
            The polymer chemical potential (multiplied with the inverse
            temperature to make the quantitiy dimensionless).

        Returns
        -------
        epsi : The attraction strength (multiplied with the inverse
            temperature to make the quantitiy dimensionless). (`float`)

        """
        epsi=np.log(np.exp(mu_pc)+1)
        return epsi

    ####################################################################
    #The inhomogenious functional:
    #In this section all the functions concerning the model specific
    #free energy functional are defined.
    ####################################################################

    def _cal_n(self):
        """Calculates the weighted densities necessary for the
        calculation of the free energy and the excess chemical
        potential.

        Returns
        -------
        Result : `tuple` of `numpy.ndaray`
        """
        n1 = self._r_c + self._r_pc1
        n2 = self._boundary_roll(self._r_c, -1, axis=1) + self._r_pc1
        n3 = self._r_c + self._r_pc2
        n4 = self._boundary_roll(self._r_c, -1, axis=0) + self._r_pc2
        n5 = self._r_pc1
        n6 = self._r_pc2
        n7 = self._r_c
        return n1, n2, n3, n4, n5, n6, n7

    def _cal_Phi_ex_AO(self):
        """Calculates the excess free energy of the AO-model.

        Returns
        -------
        Result : `np.array`
            Free energy density of the AO-model.
        """
        n=self._cal_n()
        n1=n[0]
        n2=n[1]
        n3=n[2]
        n4=n[3]
        n5=n[4]
        n6=n[5]
        n7=n[6]
        Phi0=self._cal_Phi_0
        Phi_ex = Phi0(n1)+Phi0(n2)+Phi0(n3)+Phi0(n4)-Phi0(n5)-Phi0(n6)\
                -3*Phi0(n7)
        return Phi_ex
        
    def cal_F(self):
        """Calculates the free energy of the three component system. It
        differs from the free energy functional of the AO-model just by
        a correction term accounting for the zero- and one-body
        interaction of the polymers (see description of the class). For
        getting the free energy of the lattice gas use ``cal_F_lg``,
        which is the semigrand potential, where the polymer clusters are
        treated grand canonically and the colloids canonically.

        Returns
        -------
        Result : `float`
            Free energy of the three component system (times the inverse
            temperature to make the results dimension 1).
        """
        z_pc1 = np.exp(self._mu_pc1)
        z_pc2 = np.exp(self._mu_pc2)
        r_c = self._r_c
        r_pc1 = self._r_pc1
        r_pc2 = self._r_pc2
        Phi_id = self._cal_Phi_id()
        Phi_ex = self._cal_Phi_ex_AO()
        F_id = np.sum(Phi_id)
        F_ex_AO = np.sum(Phi_ex)
        F = (F_id + F_ex_AO
                - np.log(z_pc1+1)
                *np.sum(-1+r_c+self._boundary_roll(r_c, -1, axis=1))
                - np.log(z_pc2+1)
                *np.sum(-1+r_c+self._boundary_roll(r_c, -1, axis=0)))
        return F

    def cal_F_lg(self):
        """Calculates the free energy of the lattice gas. If
        ``self.mu_fix==False`` this schould give the same result as the
        ``cal_semi_Om``-function.

        Returns
        -------
        Result : `float`
            Free energy of the lattice gas.
        """
        F_lg = self.cal_F()
        mu_pc1 = self._mu_pc1
        mu_pc2 = self._mu_pc2
        r_pc1 = self._r_pc1
        r_pc2 = self._r_pc2
        F_lg -= (mu_pc1*np.sum(r_pc1)+mu_pc2*np.sum(r_pc2))
        return F_lg

    def cal_mu_ex(self):
        n = self._cal_n()
        n1=n[0]
        n2=n[1]
        n3=n[2]
        n4=n[3]
        n5=n[4]
        n6=n[5]
        n7=n[6]
        z_pc = np.exp(self._mu_pc1)
        mu_c_ex = np.log((1-n1)*(1-self._boundary_roll(n2, 1, axis=1))\
                *(1-n3)*(1-self._boundary_roll(n4, 1, axis=0))\
                /(1-n7)**3) + 4*np.log(z_pc+1)
        mu_pc1_ex = np.log((1-n1)*(1-n2)/(1-n5))
        mu_pc2_ex = np.log((1-n3)*(1-n4)/(1-n6))
        return mu_c_ex, mu_pc1_ex, mu_pc2_ex

    ####################################################################
    #The homogenious methodes:
    #The following section contains all the methodes concerning the bulk
    #properties of the system.
    ####################################################################

    @classmethod
    def _cal_bulk_r_pc(cls, r_c, epsi):
        """Calculates the bulk polymer cluster density in dependence of
        the the coloid density and the choosen attraction strength

        Parameters
        ----------
        r_c : `float` or `np.ndarray`
            The coloid density.
        epsi : `float`
            Attraction strengs (times inverse temperature).

        Returns
        -------
        r_pc : `float` or `np.ndarray`
            The polymer cluster density.
        """
        mu_pc = cls.translate_epsi_to_mu_pc(epsi)
        z_pc = np.exp(mu_pc)
        r_pc = ((1+2*z_pc*(1-r_c))/(2*(z_pc+1))
                - 1/(2*(z_pc+1))*np.sqrt((1+2*z_pc*(1-r_c))**2 -
                 4*z_pc*(z_pc+1)*(1-r_c)**2))
        return r_pc

    @classmethod
    def _cal_bulk_dr_pc(cls, r_c, epsi):
        """Calculates the derivative of the bulk polymer cluster density
        with respect to the coloidal density in dependence of
        the the coloid density and the chosen attraction strength

        Parameters
        ----------
        r_c : `float` or `np.ndarray`
            The coloid density.
        epsi : `float`
            Attraction strength (times inverse temperature).

        Returns
        -------
        dr_pc : `float` or `np.ndarray`
            The derivative of the polymer cluster density.
        """
        mu_pc = cls.translate_epsi_to_mu_pc(epsi)
        z_pc = np.exp(mu_pc)
        dr_pc = -z_pc/(z_pc+1)\
                *(1+(1-2*r_c)/np.sqrt(4*z_pc*(1-r_c)*r_c+1))
        return dr_pc

    @classmethod
    def cal_bulk_mu_lg(cls, r_c, epsi):
        """Calculates the chemical potential for a bulk lattice gas.

        Parameters
        ----------
        r_c : `Float` or `np.ndarray`
            The coloidal density.
        epsi : `Float`
            Attraction strength

        Returns
        -------
        mu : `Float` or `np.ndarray`
            The chemical potential for the lattice gas. 
        """
        r_pc = cls._cal_bulk_r_pc(r_c, epsi)
        mu_pc = cls.translate_epsi_to_mu_pc(epsi)
        z_pc = np.exp(mu_pc)
        mu_c = (np.log(r_c) +4*cls._cal_dPhi_0(r_c+r_pc)
                -3*cls._cal_dPhi_0(r_c)-4*np.log(z_pc+1))
        return mu_c

    @classmethod
    def cal_bulk_dmu_lg(cls, r_c, epsi):
        """Calculates the derivative of the chemical potential from the
        bulk lattice gas with respect to the colloidal density.

        Parameters
        ----------
        r_c : `Float` or `np.ndarray`
            The coloidal density.
        epsi : `Float`
            Attraction strength

        Returns
        -------
        dmu : `Float` or `np.ndarray`
            The derivative of the chemical potential from the lattice
            gas.
        """
        r_pc = cls._cal_bulk_r_pc(r_c, epsi)
        dr_pc = cls._cal_bulk_dr_pc(r_c, epsi)
        mu_pc = cls.translate_epsi_to_mu_pc(epsi)
        z_pc = np.exp(mu_pc)
        dmu = 1/r_c + 4*cls._cal_d2Phi_0(r_c+r_pc)*(1+dr_pc)\
                -3*cls._cal_d2Phi_0(r_c)
        return dmu

    @classmethod
    def _cal_bulk_f_AO_id(cls, r_c, r_pc):
        """Calculates the ideal gas part of the free energy density of
        a bulk AO-system under given colloid and polymer cluster
        density.

        Parameters
        ----------
        r_c : `float`
            Coloid density
        r_pc : `float`
            Polymer cluster density

        Returns
        -------
        f_id : `float`
            The idea gas part of the free energy density.
        """
        f_id = r_c*(np.log(r_c)-1) +2*r_pc*(np.log(r_pc)-1)
        return f_id

    @classmethod
    def _cal_bulk_f_AO_ex(cls, r_c, r_pc):
        """Calculates the excess part of the free energy density of a
        bulk AO-system under given colloid and polymer cluster density.

        Parameters
        ----------
        r_c : `float`
            Coloid density
        r_pc : `float`
            Polymer cluster density

        Returns
        -------
        f_ex : `float`
            The excess part of the free energy density.
        """
        n1 = n2 = n3 = n4= r_c+r_pc
        n5 = n6 = r_pc
        n7 = r_c
        f_ex = (cls._cal_Phi_0(n1)+cls._cal_Phi_0(n2)+cls._cal_Phi_0(n3)
                +cls._cal_Phi_0(n4)-3*cls._cal_Phi_0(n7)
                -cls._cal_Phi_0(n5)-cls._cal_Phi_0(n6))
        return f_ex

    @classmethod
    def cal_bulk_f_lg(cls, r_c, epsi):
        """Calculates the free energy density of the bulk lattice gas
        under given density. (The function is the same as in
        ``cal_F_lg`` but simplified for bulk sytems.)

        Parameters
        ----------
        r_c: `float` or `np.ndarray`
            Density
        epsi: `float`
            Attraction strength (times inverse temperature)

        Returns
        -------
        f : `float` or `np.ndarray`
            The free energy density of a bulk lattice gas.
        """
        r_pc = cls._cal_bulk_r_pc(r_c, epsi)
        f_AO_id = cls._cal_bulk_f_AO_id(r_c, r_pc)
        f_AO_ex = cls._cal_bulk_f_AO_ex(r_c, r_pc)
        mu_pc = cls.translate_epsi_to_mu_pc(epsi)
        z_pc = np.exp(mu_pc)
        f_tilde = f_AO_id+f_AO_ex-2*np.log(z_pc+1)*(2*r_c-1)
        f_eff = f_tilde-2*r_pc*np.log(z_pc)
        return f_eff
    
    @classmethod
    def cal_bulk_om_lg(cls, r, epsi):
        """Calculates the grand potential density for a bulk lattice gas
        under given densitie.

        Parameters
        ----------
        r : `float` or `np.ndarray`
            The denstity.
        epsi : `float`
            The attraction strength (times inverse temperature).

        Returns
        -------
        om : `Float` 
            The grand potential density
        """
        f = cls.cal_bulk_f_lg(r, epsi)
        mu = cls.cal_bulk_mu_lg(r, epsi)
        om = f-mu*r
        return om

    @classmethod
    def cal_bulk_p(cls, r, epsi):
        """Calculates the pressure of a bulk lattice gas under given
        density.

        Parameters
        ----------
        r : `float` or `np.ndarray`
            The denstity.
        epsi : `float`
            The attraction strength (times inverse temperature).

        Returns
        -------
        The pressure : `Float`
        """
        p = -cls.cal_bulk_om_lg(r, epsi)
        return p
    
    @classmethod
    def _cal_difMu(cls, r_c, *args):
        """Calculates the difference between a certain chemical
        potential of the lattice gas and the chemical potential
        belonging to a certain density. This is a help-function for the
        function ``cal_bulk_coex_dens``.

        Parameters
        ----------
        r_c : `float`
            The coloid density of the system
        *args:
            First argument: Attraction strength (times inverse
            temperature). (`float`)
            Second argument: The reference chemical potential which the
            chemical potential for at density r_c should be compared to.
            (`float`)

        Returns
        -------
        difMu : `float`
            The difference between the two colloidal chem. pot.
        """
        epsi = args[0]
        mu_c = args[1]
        mu = cls.cal_bulk_mu_lg(r_c, epsi)
        return mu-mu_c

    @classmethod
    def cal_bulk_coex_dens(cls, mu, epsi, init_min=0.01, init_max=0.99):
        """Calculates the coexisting densities of a bulk system lattice
        gas system under given chemical potential.

        Parameters
        ----------
        mu : `Float`
            The chemical potential of the lattice gas.
        epsi : `Float`
            The attraction strength (times inverse temperatue).

        Returns
        -------
        r_coex : `Tuple`
            The coexisting densities arrangend in a tuple of the shape 
            (vapour_dens, liquid_dens)
        """
        def dmu(rc, *args):
            epsi = args[0]
            mu_c = args[1]
            return np.diag(cls.cal_bulk_dmu_lg(rc, epsi))
        r_coex = op.fsolve(cls._cal_difMu,
                           np.array([init_min, init_max]),
                           args=(epsi, mu), fprime=dmu)
        while r_coex[0]==init_min or r_coex[1]==init_max:
            init_min = init_min/2
            init_max = (init_max+1)/2
            r_coex = op.fsolve(cls._cal_difMu,
                               np.array([init_min, init_max]),
                               args=(epsi, mu), fprime=dmu)
        r_coex = tuple(r_coex)
        return r_coex


    ####################################################################
    #In the following section the abc-methodes concerning the surface
    #prperties of the mother class are overridden.
    ####################################################################

    def _cal_p(self, dens):
        epsi = self.epsi
        r = dens[0]
        p = self.cal_bulk_p(r, epsi)
        return p

    def _cal_coex_dens(self):
        mu = self._mu[0]
        epsi = self.epsi
        init_min = np.min(self._r[0])
        init_max = np.max(self._r[0])
        r_c_coex = self.cal_bulk_coex_dens(mu, epsi, init_min=init_min,
                                           init_max=init_max)
        r_pc_coex = self._cal_bulk_r_pc(np.array(r_c_coex), epsi)
        r_pc_coex = tuple(r_pc_coex)
        return [r_c_coex, r_pc_coex, r_pc_coex]



