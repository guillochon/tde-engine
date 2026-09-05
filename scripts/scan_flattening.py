#!/usr/bin/env python3
"""Flattening scan for the TDE-engine paper, and generator for window.pdf.

For each covering factor f_* = f_Omega (disk-fed, b = 1):
  C3  collision cap from the Cohn-Kulsrud flux (losscone.py), fitted
      as a power law over m6 = 0.03-300, the same grid as the published rates
  belt  virial + C1 + C2 at that Gamma(m6)
  C4  largest sigma the clouds can sustain, using a (sigma, m6) fit of the
      same CK flux

Three f_*-dependent ratios must exceed unity:
  p_TDE / p_AGN,  H / R_MC,  1e8 Msun pc^-3 / rho_0
They close the window 0.18 <= f_* <= 0.25 quoted in the manuscript.

Run from the repo root or from scripts/:
  python scripts/scan_flattening.py
"""
import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import sympy as sp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import losscone as L

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

OUTDIR = os.environ.get('TDE_FIG_OUT',
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'paper'))
os.makedirs(OUTDIR, exist_ok=True)

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Liberation Serif', 'Times New Roman', 'DejaVu Serif']
rcParams['mathtext.fontset'] = 'stix'
rcParams['axes.linewidth'] = 1.0
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = True
rcParams['ytick.right'] = True
rcParams['pdf.fonttype'] = 42

# ---------------- belt physics (same closures as solve_equilibrium_disk.py) ----
G   = sp.Float('6.67428e-8', 30)
msun = sp.Float('1.9889225e33', 30)
yr  = sp.Float('3.1556926e7', 30)
pc  = sp.Float('3.0856776e18', 30)
kb  = sp.Float('1.3806503e-16', 30)
mp  = sp.Float('1.67262158e-24', 30)
cc  = sp.Float('2.99792458e10', 30)
sigT = sp.Float('6.6524e-25', 30)
km  = sp.Float('1e5', 30)

gam = sp.Rational(7, 4)
mstar = msun / 2
Tcl = 10
muH = sp.Float('1.4', 30)
mucl = sp.Float('2.33', 30)
cs = sp.sqrt(kb * Tcl / (mucl * mp))
Lam0 = sp.Float('1.3e-27', 30)
epsacc = sp.Rational(1, 20)

m6, sig, Mach, rc, rb = sp.symbols('m6 sigma Mach R_MC r_b', positive=True)
mh = sp.Float('1e6', 30) * m6 * msun
Nst = 2 * mh / mstar
ah = 2 * G * mh / sig**2
n0 = (3 - gam) * Nst / (4 * sp.pi * ah**3)
rhoMC = 9 * mh / (4 * sp.pi * rb**3)
nMC = rhoMC / (muH * mp)
sigcl = Mach * cs
tc = 2 * rc / sigcl
E = sp.Float('1.8e50', 30) * m6**sp.Rational(1, 3)
trad = (sp.Float('1.7e4', 30) * yr
        * (E / sp.Float('1e50', 30))**sp.Rational(4, 17)
        * nMC**sp.Rational(-9, 17))
Rudr = sp.Float('1.15', 30) * (E / rhoMC)**sp.Rational(1, 5) * trad**sp.Rational(2, 5)
fUDR = Rudr**2 / (4 * rb**2)
Vmc = sp.Rational(4, 3) * sp.pi * rc**3
trelax_st = sp.Float('0.34', 30) * sig**3 / (G**2 * mstar**2 * n0 * 10)

lm, ls, lM, lr, lb = sp.symbols('lm ls lM lr lb')


def _solve(eqs, unk):
    sub = {m6: sp.exp(lm), sig: sp.exp(ls), Mach: sp.exp(lM),
           rc: sp.exp(lr), rb: sp.exp(lb)}
    n = len(unk)
    A = sp.zeros(n, n)
    B = sp.zeros(n, 1)
    for i, e in enumerate(eqs):
        Lx = sp.expand_log(sp.log(sp.powsimp(e.subs(sub), force=True)), force=True)
        for j, v in enumerate(unk):
            A[i, j] = sp.diff(Lx, v)
        B[i] = -(Lx.subs({v: 0 for v in unk + [lm]}) + sp.diff(Lx, lm) * lm)
    return A.solve(B)


def fit_ck_sigma_m6():
    """G_CK ~ C (sigma/300 km s^-1)^p m6^q, for the C4 monomial solve."""
    sigs = np.array([150, 200, 250, 300, 400, 500, 650]) * L.km
    m6s = np.array([0.3, 1.0, 3.0, 10.0])
    rows = []
    for mv in m6s:
        for sv in sigs:
            try:
                g = L.Cusp.engine(mv, sv).rate() * L.yr
            except Exception:
                continue
            if np.isfinite(g) and g > 0:
                rows.append((np.log(mv), np.log(sv / (300 * L.km)), np.log(g)))
    rows = np.array(rows)
    A = np.column_stack([np.ones(len(rows)), rows[:, 1], rows[:, 0]])
    coeff, _, _, _ = np.linalg.lstsq(A, rows[:, 2], rcond=None)
    pref, p_sig, p_m6 = np.exp(coeff[0]), coeff[1], coeff[2]
    print(f"# CK flux fit: Gamma = {pref:.4e} (sig/300km)^{p_sig:.3f} m6^{p_m6:.3f} /yr")
    return (sp.Float(str(pref), 30) / yr
            * (sig / (300 * km))**sp.Float(str(p_sig), 30)
            * m6**sp.Float(str(p_m6), 30))


G_ck = None  # filled in main() before C4 is called


def sigma_comp(F_OM):
    beta = 1 / F_OM
    Nmc_ = F_OM * 4 * rb**2 / rc**2
    trMC = sp.Float('0.34', 30) * sig**3 / (G**2 * (rhoMC * Vmc) * rhoMC * 10)
    eqs = [Mach**2 * cs**2 / (sp.Rational(4, 5) * sp.pi * G * rhoMC * rc**2),
           G_ck * fUDR * tc * beta,
           G_ck * E / (Lam0 * nMC**2 * Nmc_ * Vmc),
           trMC / trelax_st]
    x = sp.expand(_solve(eqs, [ls, lM, lr, lb])[0])
    return float(sp.exp(x.subs(lm, 0)) / km), float(sp.diff(x, lm))


def belt(Gp, Gs, F_OM):
    Geq = sp.Float(str(Gp), 30) / yr * m6**sp.Float(str(Gs), 30)
    beta = 1 / F_OM
    Nmc_ = F_OM * 4 * rb**2 / rc**2
    eqs = [Mach**2 * cs**2 / (sp.Rational(4, 5) * sp.pi * G * rhoMC * rc**2),
           Geq * fUDR * tc * beta,
           Geq * E / (Lam0 * nMC**2 * Nmc_ * Vmc)]
    s = _solve(eqs, [lM, lr, lb])
    d = {}
    for v, x in zip([Mach, rc, rb], s):
        x = sp.expand(x)
        d[v] = (sp.exp(x.subs(lm, 0))
                * m6**sp.nsimplify(sp.diff(x, lm), rational=True))
    vc2 = G * mh / rb
    rhobg = rhoMC * sigcl**2 / vc2
    nbg = rhobg / (muH * mp)
    mdot = 4 * sp.pi * rb**2 * rhobg * sp.sqrt(vc2) * sp.Float('1e-3', 30)
    LEdd = 4 * sp.pi * G * mh * mp * cc / sigT
    ev = lambda e, u=1: float(sp.N(
        sp.powsimp(sp.together((e / u).subs(d)), force=True).subs(m6, 1)))
    pw = lambda e, u=1: float(sp.N(sp.diff(sp.expand_log(sp.log(sp.powsimp(
        (e / u).subs(d), force=True).subs(m6, sp.exp(lm))), force=True), lm)))
    fEdd = sp.Rational(1, 20) * mdot * cc**2 / LEdd
    return dict(rb=ev(rb, pc), RMC=ev(rc, pc), Mach=ev(Mach), nMC=ev(nMC),
                Mgas=ev(rhoMC * Vmc * Nmc_, msun),
                AV=ev(nbg * rb / sp.Float('2.2e21', 30)),
                nbg=ev(nbg), LEdd=ev(fEdd), LEdd_p=pw(fEdd), NMC=ev(Nmc_))


def sig_coll(m6v, fst, LamC=3.0):
    mh_ = 1e6 * m6v * L.msun
    rt = (mh_ / L.MSTAR)**(1 / 3) * L.RSTAR
    hi = np.sqrt(2 * L.G * mh_ / (40 * rt))
    def f(x):
        c = L.Cusp.engine(m6v, np.exp(x))
        return np.log(c.rate_coll(LamC) / fst / c.rate())
    return np.exp(brentq(f, np.log(2e6), np.log(hi)))


def fitpow(fst, xs=(0.03, 0.3, 3, 30, 300)):
    s = np.array([sig_coll(m, fst) for m in xs])
    g = np.array([L.Cusp.engine(m, sv).rate() * L.yr for m, sv in zip(xs, s)])
    a = np.log10(np.asarray(xs))
    ps = np.polyfit(a, np.log10(s / L.km), 1)
    pg = np.polyfit(a, np.log10(g), 1)
    return (10**ps[1], ps[0]), (10**pg[1], pg[0])


def ptde(Gv):
    return (1e33 + 3e31 + 7e32) * (Gv / 1e-2)


def pagn(f):
    return 1e35 * min(f, 1) * f


def _crossing(x, y, target=1.0):
    """x at which y crosses target, interpolating in log y.

    np.interp needs an increasing xp array, so decreasing ratios (rho_0) are
    reversed first. If the root lies off the grid, linearly extrapolate from
    the two nearest points rather than clipping to an endpoint.
    """
    x = np.asarray(x, dtype=float)
    ly = np.log(np.asarray(y, dtype=float) / target)
    if ly.min() > 0 or ly.max() < 0:
        i = int(np.argmin(np.abs(ly)))
        j = 1 if i == 0 else i - 1
        k = (ly[j] - ly[i]) / (x[j] - x[i])
        return float(x[i] - ly[i] / k)
    order = np.argsort(ly)
    return float(np.interp(0.0, ly[order], x[order]))


def scan(fs):
    rows = []
    print(f"{'f_*':>5} {'sigma':>7} {'Gamma_eq':>10} {'pTDE/pAGN':>10} {'H/R_MC':>8} "
          f"{'M_floor':>9} {'A_V':>6} {'rho_0':>9}")
    for f in fs:
        (sp_, ss), (gp, gs) = fitpow(f)
        fr = sp.Rational(int(round(f * 100)), 100)
        sc, scs = sigma_comp(fr)
        b = belt(gp, gs, fr)
        ah_pc = (2 * 6.67428e-8 * 1e6 * 1.9889225e33
                 / (sp_ * 1e5)**2 / 3.0856776e18)
        rho0 = 0.5 * 1.25 * (2e6 / 0.5) / (4 * np.pi * ah_pc**3)
        H = f * b['rb']
        floor = ((sp_ / sc)**(1 / (scs - ss))) * 1e6
        rec = dict(f=f, sig=sp_, sigs=ss, G=gp, Gs=gs,
                   ratio=ptde(gp) / pagn(b['LEdd']),
                   HR=H / b['RMC'], floor=floor, rho0=rho0, **b)
        rows.append(rec)
        print(f"{f:>5.2f} {sp_:>7.0f} {gp:>10.2e} {rec['ratio']:>10.2f} "
              f"{rec['HR']:>8.2f} {floor:>9.1e} {b['AV']:>6.1f} {rho0:>9.1e}")
    return rows


def fig_window(rows, f_lo, f_hi):
    fa = np.array([r['f'] for r in rows])
    ratio = np.array([r['ratio'] for r in rows])
    hr = np.array([r['HR'] for r in rows])
    rho_r = 1e8 / np.array([r['rho0'] for r in rows])
    sig = np.array([r['sig'] for r in rows])
    gam = np.array([r['G'] for r in rows])

    BLUE = (0.15, 0.40, 0.80)
    RED = (0.80, 0.15, 0.15)
    GOLD = (0.85, 0.62, 0.10)
    GREEN = (0.55, 0.82, 0.55)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.6, 3.8),
                                   gridspec_kw={'wspace': 0.32})

    axL.axvspan(f_lo, f_hi, color=GREEN, alpha=0.45, lw=0, zorder=0)
    axL.axhline(1.0, color='k', lw=1.0, zorder=1)
    axL.plot(fa, ratio, 'o-', color=BLUE, ms=5.5, lw=1.6, label=r'$\dot{p}_{\rm TDE}/\dot{p}_{\rm AGN}$')
    axL.plot(fa, hr, 's-', color=RED, ms=5.5, lw=1.6, label=r'$H/R_{\rm MC}$')
    axL.plot(fa, rho_r, '^-', color=GOLD, ms=6.0, lw=1.6,
             label=r'$10^{8}\,M_{\odot}\,{\rm pc}^{-3}/\rho_{0}$')
    axL.set_yscale('log')
    axL.set_ylim(0.28, 22)
    axL.set_xlim(fa[0] - 0.005, fa[-1] + 0.005)
    axL.set_xlabel(r'$f_{\ast}$', fontsize=13)
    axL.set_ylabel(r'constraint ratio (need $> 1$)', fontsize=12)
    axL.tick_params(labelsize=11)
    axL.legend(loc='upper right', fontsize=9, frameon=True, fancybox=False,
               edgecolor='0.7')
    axL.set_title(rf'allowed window $f_{{\ast}} = {f_lo:.2f}$–${f_hi:.2f}$',
                  fontsize=12, pad=6)
    axL.text(0.5 * (f_lo + f_hi), 0.42, 'allowed', color=(0.15, 0.50, 0.15),
             fontsize=11, fontweight='bold', ha='center', va='center')

    axR.axvspan(f_lo, f_hi, color=GREEN, alpha=0.45, lw=0, zorder=0)
    axR.plot(fa, sig, 'o-', color=BLUE, ms=5.5, lw=1.6)
    axR.set_ylim(100, 320)
    axR.set_xlim(fa[0] - 0.005, fa[-1] + 0.005)
    axR.set_xlabel(r'$f_{\ast}$', fontsize=13)
    axR.set_ylabel(r'$\sigma_{\rm eq}$ (km s$^{-1}$)', fontsize=12, color=BLUE)
    axR.tick_params(labelsize=11, axis='y', colors=BLUE)
    axR.tick_params(labelsize=11, axis='x')
    axR.set_title('equilibrium across the window', fontsize=12, pad=6)

    axG = axR.twinx()
    axG.plot(fa, gam, 's--', color=RED, ms=5.5, lw=1.6)
    axG.set_yscale('log')
    axG.set_ylabel(r'$\Gamma_{\rm eq}$ (yr$^{-1}$)', fontsize=12, color=RED)
    axG.tick_params(labelsize=11, colors=RED)
    ymin, ymax = gam.min() / 1.4, gam.max() * 1.4
    axG.set_ylim(ymin, ymax)

    fig.tight_layout(pad=0.4)
    out = os.path.join(OUTDIR, 'window.pdf')
    fig.savefig(out)
    plt.close(fig)
    print(f'wrote {out}')


def main():
    global G_ck
    G_ck = fit_ck_sigma_m6()
    # Start below the p_TDE = p_AGN root (~0.12) so the left-panel blue
    # curve crosses unity on the plotted range.
    fs = np.round(np.arange(0.08, 0.341, 0.02), 2)
    rows = scan(fs)
    fa = np.array([r['f'] for r in rows])
    ratio = np.array([r['ratio'] for r in rows])
    hr = np.array([r['HR'] for r in rows])
    rho_r = 1e8 / np.array([r['rho0'] for r in rows])
    f_agn = _crossing(fa, ratio)
    f_geo = _crossing(fa, hr)
    f_rho = _crossing(fa, rho_r)
    print(f"\n  AGN dominance   p_TDE=p_AGN  at f_* = {f_agn:.3f}")
    print(f"  clouds fit disk H=R_MC       at f_* = {f_geo:.2f}")
    print(f"  rho0 = 1e8 Msun/pc^3         at f_* = {f_rho:.2f}")
    # Binding edges in the paper: H = R_MC and rho_0 = 1e8.
    # p_TDE = p_AGN (UDR+jets+outflows) is weaker; it fails only below ~0.12.
    f_lo, f_hi = f_geo, f_rho
    print(f"  allowed window               {f_lo:.2f} <= f_* <= {f_hi:.2f}")
    fig_window(rows, f_lo, f_hi)


if __name__ == '__main__':
    main()
