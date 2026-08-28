"""
Cohn-Kulsrud loss-cone flux, and the collision-capped equilibrium (condition C3)
of the TDE-engine paper.

FLUX (Cohn & Kulsrud 1978; Merritt 2013 ch.6; Vasiliev & Merritt 2013 eqs 44-47):

    F_lc(E) = 4 pi^2 J_c^2(E) f(E) [P(E)/t_r(E)] / ln(1/R_0(E))
    R_lc(E) = J_lc^2/J_c^2 = 4 E r_t /(G M_h)
    q(E)    = [P(E)/t_r(E)] / R_lc(E)
    alpha(q)= (q^4+q^2)^(1/4)   -> sqrt(q) for q<<1, -> q for q>>1
    R_0     = R_lc exp(-alpha)
    Gamma   = (1/m_*) int F_lc dE

CUSP normalisation is separable from the engine model so the flux code can be
validated on real nuclei before the paper's prescription is imposed.
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import beta as Beta, gamma as gammafn
from scipy.optimize import brentq

G, msun, rsun = 6.67428e-8, 1.9889225e33, 6.955e10
yr, pc, km = 3.1556926e7, 3.0856776e18, 1.0e5
GAMMA, MSTAR, RSTAR = 1.75, 0.5*msun, 0.4626*rsun
VESC2 = 2*G*MSTAR/RSTAR


class Cusp:
    def __init__(self, mh, rho_ref, r_ref, r_out, gamma=GAMMA,
                 mstar=MSTAR, rstar=RSTAR):
        self.mh, self.gamma, self.mstar, self.rstar = mh, gamma, mstar, rstar
        self.rho_ref, self.r_ref, self.r_out = rho_ref, r_ref, r_out
        self.rt  = (mh/mstar)**(1/3)*rstar
        self.lnL = np.log(0.4*mh/mstar)
        self.f0  = (rho_ref*r_ref**gamma
                    / (4*np.pi*np.sqrt(2)*Beta(gamma-0.5, 1.5)*(G*mh)**gamma))

    @classmethod
    def engine(cls, m6, sigma, **kw):
        """paper's prescription: a_h = 2GM/sigma^2, N_* = 2 M_h/m_* inside a_h"""
        mh = 1e6*m6*msun
        ah = 2*G*mh/sigma**2
        rho0 = MSTAR*(3-GAMMA)*(2*mh/MSTAR)/(4*np.pi*ah**3)
        c = cls(mh, rho0, ah, ah, **kw)
        c.ah, c.sigma, c.m6 = ah, sigma, m6
        return c

    def a_of_E(self, E): return G*self.mh/(2*E)
    def P(self, E):      return 2*np.pi*G*self.mh/(2*E)**1.5
    def Jc2(self, E):    return (G*self.mh)**2/(2*E)
    def Rlc(self, E):    return 4*E*self.rt/(G*self.mh)
    def f(self, E):      return self.f0*E**(self.gamma-1.5)
    def rho(self, r):    return self.rho_ref*(r/self.r_ref)**(-self.gamma)
    def sig_r(self, r):  return np.sqrt(G*self.mh/((1+self.gamma)*r))

    def t_r(self, E):
        """local relaxation time (retained for comparison only)"""
        a = self.a_of_E(E)
        return 0.34*self.sig_r(a)**3/(G**2*self.mstar*self.rho(a)*self.lnL)

    def Qtilde(self):
        """Stone & Metzger (2016) Eq.(A8): Q~ = 3Q~_{1/2} - Q~_{3/2} + 2Q~_0,
        accurate to 10% for 0.5 <= gamma < 2.95."""
        x = self.gamma - 0.5
        return 1.99/x - 0.0657 + 0.597*x - 0.192*x**2

    def muP(self, E):
        """mu_bar(E) P(E), the per-orbit angular-momentum diffusion, from the
        closed-form flux of Stone & Metzger (2016) App.A (their Eq. A4), which
        derives mu_bar from <(Delta v_perp)^2> rather than from t_r.

        F_A4 = 4 pi^2 J_c^2 f mu_bar P / ln(1/R_lc)  ->  invert for mu_bar P.
        """
        g, mh, m = self.gamma, self.mh, self.mstar
        K = g*(g-0.5)*gammafn(g)/gammafn(g+0.5)
        pref = (32*np.pi/(3*np.sqrt(2)))*G**5*mh**3*self.rho_ref**2 \
               * self.lnL*(G*mh/self.r_ref)**(-2*g)*K**2*self.Qtilde()
        # F_A4 = pref * E^(2g-11/2) / ln(1/Rlc)   [number flux per unit energy]
        # divide by 4 pi^2 J_c^2 f / ln(1/Rlc)  ->  mu_bar P
        return pref*E**(2*g-5.5)/(4*np.pi**2*self.Jc2(E)*self.f(E)/m)

    def q(self, E): return self.muP(E)/self.Rlc(E)

    def flux(self, E, mode='ck'):
        Rlc = self.Rlc(E)
        pre = 4*np.pi**2*self.Jc2(E)*self.f(E)
        if mode == 'pinhole':
            return pre*Rlc
        D = self.muP(E)
        if mode == 'diffusion':
            return pre*D/np.log(1/Rlc)
        q = D/Rlc
        return pre*D/(np.log(1/Rlc) + (q**4 + q**2)**0.25)

    def _erange(self):
        amin, amax = 4*self.rt, self.r_out
        return G*self.mh/(2*amax), G*self.mh/(2*amin)

    def rate(self, mode='ck'):
        Emin, Emax = self._erange()
        g = lambda lE: self.flux(np.exp(lE), mode)*np.exp(lE)
        return quad(g, np.log(Emin), np.log(Emax), limit=200)[0]/self.mstar

    def rate_coll(self, LamC=3.0):
        n0 = self.rho_ref/self.mstar
        return (8*np.pi**2*LamC*self.rstar**2*VESC2*n0**2
                * self.r_ref**3.5/np.sqrt(G*self.mh))

    def mass_check(self):
        Emin, Emax = self._erange()
        dos = lambda lE: (4*np.pi**2*self.P(np.exp(lE))*self.Jc2(np.exp(lE))
                          * self.f(np.exp(lE)))*np.exp(lE)
        return quad(dos, np.log(Emin), np.log(Emax), limit=200)[0]


def sigma_eq(m6, LamC=3.0, mode='ck', lo=2e6, hi=None):
    """solve C3: Gamma_coll(sigma) = Gamma_TDE(sigma).

    The upper bracket is capped where the cusp would shrink to the loss-cone
    floor (a_h = 40 r_t); beyond that the energy integral is not defined.
    """
    mh = 1e6*m6*msun
    rt = (mh/MSTAR)**(1/3)*RSTAR
    hi_phys = np.sqrt(2*G*mh/(40*rt))
    hi = hi_phys if hi is None else min(hi, hi_phys)
    def resid(ls):
        c = Cusp.engine(m6, np.exp(ls))
        return np.log(c.rate_coll(LamC)/c.rate(mode))
    if resid(np.log(lo))*resid(np.log(hi)) > 0:
        return np.nan
    return np.exp(brentq(resid, np.log(lo), np.log(hi)))


def main():
    print("="*72); print("1. LIMIT CHECKS"); print("="*72)
    E0 = None
    for s, lab, tgt in [(1e6, 'q >> 1', 'pinhole'), (1e-6, 'q << 1', 'diffusion')]:
        cc = Cusp.engine(1.0, 300*km); E = G*cc.mh/(2*cc.ah)
        o = cc.muP; cc.muP = lambda E_, o=o, s=s: o(E_)*s
        print(f"   q = {cc.q(E):8.1e}  ({lab}):  F_ck / F_{tgt:9s} = "
              f"{cc.flux(E,'ck')/cc.flux(E,tgt):.4f}")

    print(); print("="*72)
    print("2. VALIDATION on observed nuclei (paper's prescription NOT used)")
    print("="*72)
    rho1 = 1e6*msun*(3-GAMMA)/(4*np.pi*pc**3)
    mw = Cusp(4.3e6*msun, rho1, pc, 3*pc)
    print("   Milky Way: M_h=4.3e6 Msun, M(<1pc)=1e6 Msun, gamma=1.75, r_out=3pc")
    print(f"      Gamma_CK  = {mw.rate()*yr:.2e} /yr   [literature 1e-5 - 1e-4]")
    print(f"      q(1 pc)   = {mw.q(G*mw.mh/(2*pc)):.2f}")
    print()
    print("   M-sigma nuclei, M(<r_infl)=2M_h, gamma=1.75:")
    print(f"      {'M_h':>9} {'sig':>6} {'r_i(pc)':>9} {'Gamma(/yr)':>12} {'q(r_i)':>8}")
    Ms, Gs = [], []
    for lM in [6.0, 6.5, 7.0, 7.5, 8.0]:
        mh = 10**lM*msun; sig = 200*km*(10**lM/1e8)**(1/4.4)
        ri = G*mh/sig**2
        cc = Cusp(mh, 2*mh*(3-GAMMA)/(4*np.pi*ri**3), ri, ri)
        r = cc.rate()*yr; Ms.append(lM); Gs.append(np.log10(r))
        print(f"      {10**lM:>9.1e} {sig/km:>6.0f} {ri/pc:>9.2f} {r:>12.2e} "
              f"{cc.q(G*mh/(2*ri)):>8.2f}")
    print(f"      -> Gamma ~ M_h^{np.polyfit(Ms,Gs,1)[0]:+.2f}"
          "   [Stone & Metzger 2016: ~M_h^-0.4]")


if __name__ == '__main__':
    main()
