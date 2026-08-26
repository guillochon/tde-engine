"""Disk-geometry variant of the TDE engine equilibrium.

Geometry parameters
-------------------
f_Om : aspect ratio H/r of the cloud disk = fraction of 4pi it subtends.
b    : fraction of UDRs that intercept the disk.
         b = f_Om  -> isotropic cusp (UDRs isotropic)
         b = 1     -> disk-fed cusp (parabolic orbits return debris to the plane)
f_st : fraction of 4pi occupied by the stellar cusp (flattening).
         f_st = 1     -> isotropic cusp
         f_st = f_Om  -> cusp as flat as the gas disk

Conditions C1 and C2 depend on geometry only through beta = b/f_Om.
  beta = 1        -> identical to the isotropic-shell solution (invariance)
  beta = 1/f_Om   -> beamed driving
Condition C3 gains a factor 1/f_st on Gamma_coll (a local-density process),
while the pinhole Gamma_TDE is set by N and P and is not directly enhanced.
"""
import sys
import sympy as sp

G   = sp.Float('6.67428e-8',30); msun=sp.Float('1.9889225e33',30); rsun=sp.Float('6.955e10',30)
yr  = sp.Float('3.1556926e7',30); pc=sp.Float('3.0856776e18',30); kb=sp.Float('1.3806503e-16',30)
mp  = sp.Float('1.67262158e-24',30); cc=sp.Float('2.99792458e10',30); sigT=sp.Float('6.6524e-25',30)
km  = sp.Float('1e5',30)

gam   = sp.Rational(7,4)
mstar = msun/2
M = sp.Float('0.5',30)
rstar = rsun*(sp.Float('1.715359',30)*M**sp.Rational(5,2)+sp.Float('6.597788',30)*M**sp.Rational(13,2)
        +sp.Float('10.08855',30)*M**11+sp.Float('1.012495',30)*M**19+sp.Float('0.07490166',30)*M**sp.Rational(39,2)) \
        /(sp.Float('0.01077422',30)+sp.Float('3.082234',30)*M**2+sp.Float('17.84778',30)*M**sp.Rational(17,2)
        +M**sp.Rational(37,2)+sp.Float('0.00022582',30)*M**sp.Rational(39,2))
vesc2 = 2*G*mstar/rstar
Tcl   = 10; muH = sp.Float('1.4',30); mucl = sp.Float('2.33',30)
cs    = sp.sqrt(kb*Tcl/(mucl*mp))
Lam0  = sp.Float('1.3e-27',30)
LamC  = 3
epsacc= sp.Rational(1,20)

# ---------------- geometry ----------------
F_OM = sp.Rational(1,5)      # H/r of the disk; fiducial 0.2
BEAM = 1                     # b: all UDRs return to the plane (disk-fed cusp)
F_ST = sp.Rational(1,5)      # cusp flattening; = F_OM for a disk-fed cusp
if len(sys.argv) > 1:
    F_OM = sp.Rational(sys.argv[1])
if len(sys.argv) > 2:
    BEAM = sp.Rational(sys.argv[2])
if len(sys.argv) > 3:
    F_ST = sp.Rational(sys.argv[3])
beta = BEAM/F_OM
print(f"# geometry: f_Om={F_OM}  b={BEAM}  f_*={F_ST}   -> beta=b/f_Om={beta}")

m6, sig, Mach, rc, rb = sp.symbols('m6 sigma Mach R_MC r_b', positive=True)
mh  = sp.Float('1e6',30)*m6*msun
Nst = 2*mh/mstar
ah  = 2*G*mh/sig**2
n0  = (3-gam)*Nst/(4*sp.pi*ah**3)      # spherical-equivalent normalization
sigh= sp.sqrt(G*mh/ah)

# ---- C3: collision cap.  Collisions are a local-density process: flattening
#      the cusp into a fraction f_* of 4pi raises n by 1/f_* and lowers the
#      volume by f_*, so Gamma_coll ~ n^2 V scales as 1/f_*.
#      The pinhole (full loss cone) rate is set by N and P, not by local
#      density, so it is left unenhanced here.
Gfull = sp.Float('4e-5',30)/yr*m6**sp.Rational(4,3)*(n0/(sp.Float('1e5',30)*pc**-3)) \
        *(100*km/sigh)**3*(sigh/sig)**(2*gam-3)
Gcoll = 8*sp.pi**2*rstar**2*vesc2*n0**2*ah**sp.Rational(7,2)*(G*mh)**sp.Rational(-1,2)*LamC/F_ST

lm = sp.Symbol('lm6'); ls = sp.Symbol('ls')
def loglin(expr, vars_):
    return sp.expand_log(sp.log(sp.powsimp(expr.subs({v: sp.exp(s) for v,s in vars_.items()}), force=True)), force=True)

L = loglin(Gcoll/Gfull, {m6: lm, sig: ls})
ls_sol = sp.solve(sp.Eq(L,0), ls)[0]
sig_pref = sp.exp(ls_sol.subs(lm,0)); sig_slope = sp.nsimplify(sp.diff(ls_sol,lm), rational=True)
print(f"sigma_eq = {sp.N(sig_pref/km,4)} km/s * m6^{sig_slope}")
sig_eq = sig_pref*m6**sig_slope

def powlaw(name, expr, unit=1, un=''):
    e = sp.powsimp(expr.subs(sig, sig_eq), force=True)
    Lx = sp.expand_log(sp.log(sp.powsimp((e/unit).subs(m6, sp.exp(lm)), force=True)), force=True)
    slope = sp.nsimplify(sp.diff(Lx, lm), rational=True); pref = sp.exp(Lx.subs(lm,0))
    print(f"{name:28s} = {sp.N(pref,4)} m6^{sp.N(slope,4)} {un}")
    return pref, slope, e

Geq = powlaw('Gamma_eq [yr^-1]', Gfull, 1/yr, 'yr^-1')[2]
powlaw('a_h [pc]', ah, pc, 'pc')
powlaw('n_0 [pc^-3]', n0, pc**-3, 'pc^-3')
Ph = 2*sp.pi*sp.sqrt(ah**3/(G*mh))
powlaw('P_h [yr]', Ph, yr, 'yr')
trelaxst = sp.Float('0.34',30)*sig**3/(G**2*mstar**2*n0*10)
powlaw('t_relax,* [yr]', trelaxst, yr, 'yr')

# ---- belt / disk ----
rhoMC = 9*mh/(4*sp.pi*rb**3)
nMC   = rhoMC/(muH*mp)
sigcl = Mach*cs
tc    = 2*rc/sigcl
E     = sp.Float('1.8e50',30)*m6**sp.Rational(1,3)
trad  = sp.Float('1.7e4',30)*yr*(E/sp.Float('1e50',30))**sp.Rational(4,17)*nMC**sp.Rational(-9,17)
Rudr  = sp.Float('1.15',30)*(E/rhoMC)**sp.Rational(1,5)*trad**sp.Rational(2,5)
fUDR  = Rudr**2/(4*rb**2)
Nmc   = F_OM*4*rb**2/rc**2          # clouds needed to pave a disk of aspect ratio f_Om
Vmc   = sp.Rational(4,3)*sp.pi*rc**3
Mmc   = rhoMC*Vmc

lM, lr, lb = sp.symbols('lM lr lb')
vars_ = {m6: lm, Mach: lM, rc: lr, rb: lb}
A = sp.zeros(3,3); Bv = sp.zeros(3,1)
for i, ratio in enumerate([Mach**2*cs**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2),   # virial
                           Geq*fUDR*tc*beta,                                        # C1
                           Geq*E*BEAM/(Lam0*nMC**2*Nmc*Vmc)]):                      # C2
    Lx = loglin(ratio, vars_)
    for j,v in enumerate([lM, lr, lb]):
        A[i,j] = sp.diff(Lx, v)
    Bv[i] = -(Lx.subs({lM:0,lr:0,lb:0,lm:0}) + sp.diff(Lx,lm)*lm)
sol = A.solve(Bv)
subs2 = {}
for v, name, unit, un in [(Mach,'Mach',1,''),(rc,'R_MC',pc,'pc'),(rb,'r_b',pc,'pc')]:
    e = sp.expand(sol[[Mach,rc,rb].index(v)])
    slope = sp.nsimplify(sp.diff(e,lm), rational=True); pref = sp.exp(e.subs(lm,0))
    subs2[v] = pref*m6**slope
    print(f"{name:28s} = {sp.N(pref/unit,4)} m6^{sp.N(slope,4)} {un}")

def powlaw2(name, expr, unit=1, un=''):
    e = sp.powsimp(expr.subs(sig, sig_eq).subs(subs2), force=True)
    Lx = sp.expand_log(sp.log(sp.powsimp((e/unit).subs(m6, sp.exp(lm)), force=True)), force=True)
    slope = sp.nsimplify(sp.diff(Lx, lm), rational=True); pref = sp.exp(Lx.subs(lm,0))
    print(f"{name:28s} = {sp.N(pref,4)} m6^{sp.N(slope,4)} {un}")
    return pref, slope

print('--- derived ---')
powlaw2('n_MC [cm^-3]', nMC, 1, 'cm^-3')
powlaw2('sigma_cl [km/s]', sigcl, km, 'km/s')
powlaw2('tau_c [yr]', tc, yr, 'yr')
powlaw2('M_MC [Msun]', Mmc, msun, 'Msun')
powlaw2('N_MC', Nmc)
powlaw2('M_gas,disk [Msun]', Mmc*Nmc, msun, 'Msun')
powlaw2('R_UDR [pc]', Rudr, pc, 'pc')
powlaw2('tau_rad [yr]', trad, yr, 'yr')
powlaw2('f_UDR', fUDR)
powlaw2('L_UDR [erg/s]', Geq*E*BEAM)
vc2 = G*mh/rb
rhobg = rhoMC*sigcl**2/vc2
powlaw2('n_bg [cm^-3]', rhobg/(muH*mp), 1, 'cm^-3')
mdot = 4*sp.pi*rb**2*rhobg*sp.sqrt(vc2)*sp.Float('1e-3',30)
LEdd = 4*sp.pi*G*mh*mp*cc/sigT
powlaw2('L_AGN/L_Edd (fB=1e-3)', epsacc*mdot*cc**2/LEdd)
powlaw2('A_V (bg, r_b)', (rhobg/(muH*mp))*rb/sp.Float('2.2e21',30))
powlaw2('Mdot_TDE/Mdot_amb', (mstar/2)*Geq/mdot)
rt = (mh/mstar)**sp.Rational(1,3)*rstar
powlaw2('a_h/r_b', ah/rb)
powlaw2('R_MC/r_b', rc/rb)
print('--- lifetimes ---')
powlaw2('t_cusp = N*/Geq [yr]', Nst/Geq, yr, 'yr')
powlaw2('t_gas [yr]', (Mmc*Nmc)/(mstar*Geq), yr, 'yr')
powlaw2('m_* Geq [Msun/yr]', mstar*Geq/msun*yr)
powlaw2('q_* (cusp)', Ph*ah/(2*rt*trelaxst))
powlaw2('Jeans check', sigcl**2*rc/(G*Mmc))
powlaw2('hits per crossing', Geq*tc)
