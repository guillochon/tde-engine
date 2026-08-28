"""NOTE: this script implements the SUPERSEDED monomial loss-cone rate
(the old Eq. 25, Gamma ~ m6^{4/3} n0 sigma_h^{-3}). That form is degenerate
with the collision rate in sigma and has been replaced in the manuscript by
the bridged Cohn-Kulsrud flux; see scripts/losscone.py. Retained only to
reproduce the belt (C1, C2, virial) solve and the C4 criterion, which are
unchanged. Equilibrium rates printed here are NOT the published values.
"""
import sympy as sp

# constants (cgs)
G   = sp.Float('6.67428e-8',30); msun=sp.Float('1.9889225e33',30); rsun=sp.Float('6.955e10',30)
yr  = sp.Float('3.1556926e7',30); pc=sp.Float('3.0856776e18',30); kb=sp.Float('1.3806503e-16',30)
mp  = sp.Float('1.67262158e-24',30); cc=sp.Float('2.99792458e10',30); sigT=sp.Float('6.6524e-25',30)
km  = sp.Float('1e5',30)

# stellar / model parameters
gam   = sp.Rational(7,4)
mstar = msun/2
M = sp.Float('0.5',30)
rstar = rsun*(sp.Float('1.715359',30)*M**sp.Rational(5,2)+sp.Float('6.597788',30)*M**sp.Rational(13,2)
        +sp.Float('10.08855',30)*M**11+sp.Float('1.012495',30)*M**19+sp.Float('0.07490166',30)*M**sp.Rational(39,2)) \
        /(sp.Float('0.01077422',30)+sp.Float('3.082234',30)*M**2+sp.Float('17.84778',30)*M**sp.Rational(17,2)
        +M**sp.Rational(37,2)+sp.Float('0.00022582',30)*M**sp.Rational(39,2))
vesc2 = 2*G*mstar/rstar
Tcl   = 10; muH = sp.Float('1.4',30); mucl = sp.Float('2.33',30)
cs    = sp.sqrt(kb*Tcl/(mucl*mp))               # isothermal sound speed, 10 K molecular gas
Lam0  = sp.Float('1.3e-27',30)                  # CO cooling coefficient  L = Lam0 n^2 per cm^3
LamC  = 3                                        # Coulomb-like log for collisions, ln(a_h/r_f)+1
epsacc= sp.Rational(1,20)

m6, sig, Mach, rc, rb = sp.symbols('m6 sigma Mach R_MC r_b', positive=True)
mh  = sp.Float('1e6',30)*m6*msun
Nst = 2*mh/mstar
ah  = 2*G*mh/sig**2
n0  = (3-gam)*Nst/(4*sp.pi*ah**3)
sigh= sp.sqrt(G*mh/ah)                           # = sigma/sqrt(2)

# ---- Step 1: cluster equilibrium sigma from  Gamma_coll = Gamma_full ----
Gfull = sp.Float('4e-5',30)/yr*m6**sp.Rational(4,3)*(n0/(sp.Float('1e5',30)*pc**-3)) \
        *(100*km/sigh)**3*(sigh/sig)**(2*gam-3)
Gcoll = 8*sp.pi**2*rstar**2*vesc2*n0**2*ah**sp.Rational(7,2)*(G*mh)**sp.Rational(-1,2)*LamC

lm = sp.Symbol('lm6'); ls = sp.Symbol('ls')
def loglin(expr, vars_):
    L = sp.expand_log(sp.log(sp.powsimp(expr.subs({v: sp.exp(s) for v,s in vars_.items()}), force=True)), force=True)
    return L
L = loglin(Gcoll/Gfull, {m6: lm, sig: ls})
ls_sol = sp.solve(sp.Eq(L,0), ls)[0]
sig_pref = sp.exp(ls_sol.subs(lm,0)); sig_slope = sp.nsimplify(sp.diff(ls_sol,lm), rational=True)
print(f"sigma_eq = {sp.N(sig_pref/km,4)} km/s * m6^{sig_slope} (={sp.N(sig_slope,4)})")
sig_eq = sig_pref*m6**sig_slope

def powlaw(name, expr, unit=1, unitname=''):
    e = sp.powsimp(expr.subs(sig, sig_eq), force=True)
    Lx = sp.expand_log(sp.log(sp.powsimp((e/unit).subs(m6, sp.exp(lm)), force=True)), force=True)
    slope = sp.nsimplify(sp.diff(Lx, lm), rational=True)
    pref  = sp.exp(Lx.subs(lm,0))
    print(f"{name:26s} = {sp.N(pref,4)} m6^{sp.N(slope,4)}  ({slope}) {unitname}")
    return pref, slope, e

Geq = powlaw('Gamma_eq [yr^-1]', Gfull, 1/yr, 'yr^-1')[2]
powlaw('Gamma_coll [yr^-1]', Gcoll, 1/yr, 'yr^-1')
powlaw('a_h [pc]', ah, pc, 'pc')
powlaw('n_0 [pc^-3]', n0, pc**-3, 'pc^-3')
powlaw('sigma_h [km/s]', sigh, km, 'km/s')
Ph = 2*sp.pi*sp.sqrt(ah**3/(G*mh))
powlaw('P_h=2pi(ah^3/GMh)^.5 [yr]', Ph, yr, 'yr')
trelaxst = sp.Float('0.34',30)*sig**3/(G**2*mstar**2*n0*10)
powlaw('t_relax,* [yr]', trelaxst, yr, 'yr')

# ---- Step 2: cloud belt ----
rhoMC = 9*mh/(4*sp.pi*rb**3)          # tidally limited at r_b
nMC   = rhoMC/(muH*mp)
sigcl = Mach*cs
tc    = 2*rc/sigcl
E     = sp.Float('1.8e50',30)*m6**sp.Rational(1,3)
trad  = sp.Float('1.7e4',30)*yr*(E/sp.Float('1e50',30))**sp.Rational(4,17)*nMC**sp.Rational(-9,17)
Rudr  = sp.Float('1.15',30)*(E/rhoMC)**sp.Rational(1,5)*trad**sp.Rational(2,5)
fUDR  = Rudr**2/(4*rb**2)
Nmc   = 4*rb**2/rc**2
Vmc   = sp.Rational(4,3)*sp.pi*rc**3
Mmc   = rhoMC*Vmc

eqC = Mach**2*cs**2 - sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2          # virial support
eqD = Geq*fUDR*tc - 1                                               # one hit per crossing
eqE = Geq*E - Lam0*nMC**2*Nmc*Vmc                                   # heating = cooling

lM, lr, lb = sp.symbols('lM lr lb')
vars_ = {m6: lm, Mach: lM, rc: lr, rb: lb}
A = sp.zeros(3,3); Bv = sp.zeros(3,1)
for i, ratio in enumerate([Mach**2*cs**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2),
                           Geq*fUDR*tc,
                           Geq*E/(Lam0*nMC**2*Nmc*Vmc)]):
    Lx = loglin(ratio, vars_)
    for j,v in enumerate([lM, lr, lb]):
        A[i,j] = sp.diff(Lx, v)
    Bv[i] = -(Lx.subs({lM:0,lr:0,lb:0,lm:0}) + sp.diff(Lx,lm)*lm)
sol = A.solve(Bv)
subs2 = {}
for v, name, unit, un in [(Mach,'Mach',1,''),(rc,'R_MC',pc,'pc'),(rb,'r_b',pc,'pc')]:
    idx = [Mach,rc,rb].index(v)
    e = sp.expand(sol[idx])
    slope = sp.nsimplify(sp.diff(e,lm), rational=True)
    pref  = sp.exp(e.subs(lm,0))
    subs2[v] = pref*m6**slope
    print(f"{name:26s} = {sp.N(pref/unit,4)} m6^{sp.N(slope,4)}  ({slope}) {un}")

def powlaw2(name, expr, unit=1, unitname=''):
    e = sp.powsimp(expr.subs(sig, sig_eq).subs(subs2), force=True)
    Lx = sp.expand_log(sp.log(sp.powsimp((e/unit).subs(m6, sp.exp(lm)), force=True)), force=True)
    slope = sp.nsimplify(sp.diff(Lx, lm), rational=True)
    pref  = sp.exp(Lx.subs(lm,0))
    print(f"{name:26s} = {sp.N(pref,4)} m6^{sp.N(slope,4)}  ({slope}) {unitname}")
    return pref, slope

print('--- derived belt/cloud ---')
powlaw2('rho_MC [g/cm3]', rhoMC)
powlaw2('n_MC [cm^-3]', nMC, 1, 'cm^-3')
powlaw2('sigma_cl [km/s]', sigcl, km, 'km/s')
powlaw2('tau_c [yr]', tc, yr, 'yr')
powlaw2('M_MC [Msun]', Mmc, msun, 'Msun')
powlaw2('N_MC', Nmc)
powlaw2('M_gas,belt [Msun]', Mmc*Nmc, msun, 'Msun')
powlaw2('R_UDR [pc]', Rudr, pc, 'pc')
powlaw2('tau_rad [yr]', trad, yr, 'yr')
powlaw2('f_UDR', fUDR)
powlaw2('L_UDR [erg/s]', Geq*E)
vc2 = G*mh/rb
powlaw2('v_c(r_b) [km/s]', sp.sqrt(vc2), km, 'km/s')
rhobg = rhoMC*sigcl**2/vc2
powlaw2('n_bg [cm^-3]', rhobg/(muH*mp), 1, 'cm^-3')
mdot = 4*sp.pi*rb**2*rhobg*sp.sqrt(vc2)*sp.Float('1e-3',30)
LEdd = 4*sp.pi*G*mh*mp*cc/sigT
powlaw2('L_AGN/L_Edd (fB=1e-3)', epsacc*mdot*cc**2/LEdd)
powlaw2('A_V (bg, r_b)', (rhobg/(muH*mp))*rb/sp.Float('2.2e21',30))
powlaw2('Mdot_TDE/Mdot_amb (fB=1e-3)', (mstar/2)*Geq/mdot)
powlaw2('t_relax,MC [yr]', sp.Float('0.34',30)*sig**3/(G**2*Mmc*rhoMC*10), yr, 'yr')
# loss-cone fullness from MC perturbers: q = (P_orb / t_J), t_J=(2 rt/a) t_rlx
rt = (mh/mstar)**sp.Rational(1,3)*rstar
powlaw2('q_MC (>>1 => full LC)', Ph*ah/(2*rt*sp.Float('0.34',30)*sig**3/(G**2*Mmc*rhoMC*10)))
powlaw2('r_t [pc]', rt, pc, 'pc')
powlaw2('a_h/r_b', ah/rb)
powlaw2('R_MC/r_b', rc/rb)

print('--- self-consistency ---')
# loss-cone fullness of compressed cusp via star-star relaxation
qstar = Ph*ah/(2*rt*trelaxst)
powlaw2('q_* (cusp)', qstar)
powlaw2('m_* Gamma_eq [Msun/yr]', mstar*Geq*yr/msun/yr*yr, 1, '')
tgas = (Mmc*Nmc)/(mstar*Geq)
powlaw2('t_gas = M_belt/(m* G_eq) [yr]', tgas, yr, 'yr')
powlaw2('N_* consumed in 1e8 yr / N_*', Geq*sp.Float('1e8',30)*yr/Nst)
powlaw2('Gamma_eq * tau_c (hits per crossing)', Geq*tc)
powlaw2('E_UDR [erg]', E)
powlaw2('Sigma_gas belt [Msun/pc^2]', Mmc*Nmc/(sp.pi*rb**2), msun/pc**2, 'Msun/pc^2')
powlaw2('t_engine consume cusp [yr]', Nst/Geq, yr, 'yr')
powlaw2('Jeans/cloud mass i.e. sigcl^2 R/(G M)', sigcl**2*rc/(G*Mmc))
powlaw2('L_Edd [erg/s]', LEdd)
powlaw2('Mdot_amb (fB=1e-3) [Msun/yr]', mdot, msun/yr, 'Msun/yr')
powlaw2('sig_eq/sig_Msigma(200*(m6/1e2.6)^0.23?)', sig_eq/km)

# ---- Step 3: can the clouds actually deliver this compression? (condition C4) ----
# C1-C3 fix sigma from the collision cap alone; nothing in them checks that the
# clouds are capable of compressing the cusp to that sigma in the first place.
# Impose instead that the clouds dominate the relaxation of the cusp,
#     t_relax,MC = t_relax,*,
# and re-solve the belt (virial, C1, C2, C4) for sigma.  This is the largest
# sigma the clouds can sustain.  The engine sits at the SMALLER of the two:
# if sigma_comp > sigma_coll the cap binds and the solution above stands;
# if sigma_comp < sigma_coll the clouds are the bottleneck and the rate is lower.
print('--- C4: compression criterion ---')
trMC_s = sp.Float('0.34',30)*sig**3/(G**2*Mmc*rhoMC*10)
A4 = sp.zeros(4,4); B4 = sp.zeros(4,1)
lsy = sp.Symbol('ls4')
vars4 = {m6: lm, sig: lsy, Mach: lM, rc: lr, rb: lb}
for i, ratio in enumerate([Mach**2*cs**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2),
                           Gfull*fUDR*tc,
                           Gfull*E/(Lam0*nMC**2*Nmc*Vmc),
                           trMC_s/trelaxst]):
    Lx = loglin(ratio, vars4)
    for j, v in enumerate([lsy, lM, lr, lb]):
        A4[i,j] = sp.diff(Lx, v)
    B4[i] = -(Lx.subs({lsy:0, lM:0, lr:0, lb:0, lm:0}) + sp.diff(Lx, lm)*lm)
s4 = sp.expand(A4.solve(B4)[0])
sc_pref = sp.exp(s4.subs(lm,0)); sc_slope = sp.N(sp.diff(s4, lm), 4)
print(f"sigma_comp = {sp.N(sc_pref/km,4)} km/s * m6^{sc_slope}   (clouds' capability)")
print(f"sigma_coll = {sp.N(sig_pref/km,4)} km/s * m6^{sig_slope}   (collision cap)")
mcrit = sp.N((sp.N(sig_pref/sc_pref))**(1/(sc_slope - sp.N(sig_slope))), 4)
print(f"crossover at m6 = {mcrit}")
print(f"  m6 > {mcrit}: collision-capped, solution above stands")
print(f"  m6 < {mcrit}: compression-limited, sigma and Gamma fall below the values above")
