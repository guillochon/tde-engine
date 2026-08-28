"""Flattening scan with the Stone & Metzger normalisation."""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, sympy as sp, losscone as L
from scipy.optimize import brentq
exec(open('solve_equilibrium_4x4.py').read().split('lm,ls,lM,lr,lb')[0])
lm,ls,lM,lr,lb=sp.symbols('lm ls lM lr lb')
cf=np.load('/tmp/ckfit2.npy')
G_ck=sp.Float(str(np.exp(cf[0])),30)/yr*(sig/(300*km))**sp.Float(str(cf[1]),30)*m6**sp.Float(str(cf[2]),30)

def _solve(eqs,unk):
    sub={m6:sp.exp(lm),sig:sp.exp(ls),Mach:sp.exp(lM),rc:sp.exp(lr),rb:sp.exp(lb)}
    n=len(unk);A=sp.zeros(n,n);B=sp.zeros(n,1)
    for i,e in enumerate(eqs):
        Lx=sp.expand_log(sp.log(sp.powsimp(e.subs(sub),force=True)),force=True)
        for j,v in enumerate(unk): A[i,j]=sp.diff(Lx,v)
        B[i]=-(Lx.subs({v:0 for v in unk+[lm]})+sp.diff(Lx,lm)*lm)
    return A.solve(B)

def sigma_comp(F_OM):
    beta=1/F_OM; Nmc_=F_OM*4*rb**2/rc**2
    trMC=sp.Float('0.34',30)*sig**3/(G**2*(rhoMC*Vmc)*rhoMC*10)
    eqs=[Mach**2*cs**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2),G_ck*fUDR*tc*beta,
         G_ck*E/(Lam0*nMC**2*Nmc_*Vmc),trMC/trelax_st]
    x=sp.expand(_solve(eqs,[ls,lM,lr,lb])[0])
    return float(sp.exp(x.subs(lm,0))/km),float(sp.diff(x,lm))

def belt(Gp,Gs,F_OM):
    Geq=sp.Float(str(Gp),30)/yr*m6**sp.Float(str(Gs),30)
    beta=1/F_OM; Nmc_=F_OM*4*rb**2/rc**2
    eqs=[Mach**2*cs**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2),Geq*fUDR*tc*beta,
         Geq*E/(Lam0*nMC**2*Nmc_*Vmc)]
    s=_solve(eqs,[lM,lr,lb]);d={}
    for v,x in zip([Mach,rc,rb],s):
        x=sp.expand(x);d[v]=sp.exp(x.subs(lm,0))*m6**sp.nsimplify(sp.diff(x,lm),rational=True)
    vc2=G*mh/rb;rhobg=rhoMC*sigcl**2/vc2;nbg=rhobg/(muH*mp)
    cc=sp.Float('2.99792458e10',30);sigT=sp.Float('6.6524e-25',30)
    mdot=4*sp.pi*rb**2*rhobg*sp.sqrt(vc2)*sp.Float('1e-3',30)
    LEdd=4*sp.pi*G*mh*mp*cc/sigT
    ev=lambda e,u=1: float(sp.N(sp.powsimp(sp.together((e/u).subs(d)),force=True).subs(m6,1)))
    pw=lambda e,u=1: float(sp.N(sp.diff(sp.expand_log(sp.log(sp.powsimp((e/u).subs(d),force=True).subs(m6,sp.exp(lm))),force=True),lm)))
    return dict(rb=ev(rb,pc),RMC=ev(rc,pc),Mach=ev(Mach),nMC=ev(nMC),
                Mgas=ev(rhoMC*Vmc*Nmc_,msun),AV=ev(nbg*rb/sp.Float('2.2e21',30)),
                nbg=ev(nbg),LEdd=ev(sp.Rational(1,20)*mdot*cc**2/LEdd),
                LEdd_p=pw(sp.Rational(1,20)*mdot*cc**2/LEdd),NMC=ev(Nmc_))

def sig_coll(m6v,fst,LamC=3.0):
    mh_=1e6*m6v*L.msun;rt=(mh_/L.MSTAR)**(1/3)*L.RSTAR
    hi=np.sqrt(2*L.G*mh_/(40*rt))
    f=lambda x:(lambda c: np.log(c.rate_coll(LamC)/fst/c.rate()))(L.Cusp.engine(m6v,np.exp(x)))
    return np.exp(brentq(f,np.log(2e6),np.log(hi)))
def fitpow(fst,xs=(0.03,0.3,3,30,300)):
    s=np.array([sig_coll(m,fst) for m in xs]);g=np.array([L.Cusp.engine(m,sv).rate()*L.yr for m,sv in zip(xs,s)])
    a=np.log10(xs);ps=np.polyfit(a,np.log10(s/L.km),1);pg=np.polyfit(a,np.log10(g),1)
    return (10**ps[1],ps[0]),(10**pg[1],pg[0])
def ptde(Gv): return (1e33+3e31+7e32)*(Gv/1e-2)
def pagn(f): return 1e35*min(f,1)*f

fs=[0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30]
print(f"{'f_*':>5} {'sigma':>7} {'Gamma_eq':>10} {'pTDE/pAGN':>10} {'H/R_MC':>8} "
      f"{'M_floor':>9} {'A_V':>6} {'rho_0':>9}")
out=[]
for f in fs:
    (sp_,ss),(gp,gs)=fitpow(f); fr=sp.Rational(int(round(f*100)),100)
    sc,scs=sigma_comp(fr); b=belt(gp,gs,fr)
    ah=2*6.67428e-8*1e6*1.9889225e33/(sp_*1e5)**2/3.0856776e18
    rho0=0.5*1.25*(2e6/0.5)/(4*np.pi*ah**3)
    H=f*b['rb']; floor=((sp_/sc)**(1/(scs-ss)))*1e6
    out.append(dict(f=f,sig=sp_,sigs=ss,G=gp,Gs=gs,ratio=ptde(gp)/pagn(b['LEdd']),
                    HR=H/b['RMC'],floor=floor,rho0=rho0,**b))
    print(f"{f:>5} {sp_:>7.0f} {gp:>10.2e} {out[-1]['ratio']:>10.2f} {H/b['RMC']:>8.2f} "
          f"{floor:>9.1e} {b['AV']:>6.1f} {rho0:>9.1e}")
import pickle; pickle.dump(out,open('/tmp/scan2.pkl','wb'))
r=np.array([o['ratio'] for o in out]); h=np.array([o['HR'] for o in out]); fl=np.array([o['floor'] for o in out])
fa=np.array(fs)
print(f"\n  AGN dominance   p_TDE=p_AGN  at f_* = {np.interp(0,np.log(r),fa):.2f}")
print(f"  clouds fit disk H=R_MC       at f_* = {np.interp(0,np.log(h),fa):.2f}")
print(f"  engine reaches 1e6 Msun      at f_* = {np.interp(0,np.log(fl/1e6),fa):.2f}")
