#!/usr/bin/env python3
"""Regenerate mom-rates.pdf, dominance.pdf, domratio.pdf for the TDE engine paper.

Momentum injection rates (dyn = g cm s^-2), consistent with Section 2 of the paper:
  p_UDR = 1e33 (Gamma/1e-2) m6^{1/3}
  p_TJ  = 3e31 (Gamma/1e-2)          (TDE jets)
  p_TO  = 7e32 (Gamma/1e-2)          (TDE super-Eddington outflows)
  p_TDE = p_UDR + p_TJ + p_TO
  p_SNe = 2e33 (Gamma_SNe/1e-3)
  p_AGN = 1e35 * min(f,1) * f * m6, with f = f_Edd coupled to the engine:
          f(q,Gamma,fB) = 0.06 (fB/1e-3) m6^{-0.2114} * (Gamma/Gamma_eq(m6))
          Gamma_eq(m6) = 0.037 m6^{-1/6}
  => f = 1.6 (fB/1e-3) * Gamma * m6^{-0.0447}
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import to_rgba
import os

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

# ---------------- physics ----------------
GAM_EQ_C, GAM_EQ_P = 0.037, -1.0/6.0          # Gamma_eq = 0.037 m6^{-1/6} yr^-1
# f_Edd = FEDD_C * fB3 * Gamma * m6^FEDD_P, normalised at the equilibrium.
# Recomputed from the belt solve driven by the Cohn-Kulsrud loss-cone rate.
# NOTE the disk-fed normalisation is very different (C=47.6, p=+0.786); the
# shaded regions below use the spherical one, as in the published figure.
# Set FEDD_SOURCE='published' to hold the AGN boundary at its published
# calibration, which isolates how much of the change in these figures comes
# from the new equilibrium curves alone rather than from renormalising f_Edd.
FEDD_SOURCE = os.environ.get('TDE_FEDD', 'ck')
if FEDD_SOURCE == 'published':
    FEDD_C, FEDD_P = 1.62, -0.0447          # old, Eq.(25) equilibrium
else:
    # normalised at the f_* = 0.22 equilibrium: L/L_Edd = 0.0511 at
    # Gamma_eq = 1.419e-2, with the mass slope from the belt solve
    FEDD_C, FEDD_P = 3.60, +0.826

def m6_of(lq):
    return 10.0**(np.asarray(lq) - 6.0)

def p_tde(lq, G):
    m6 = m6_of(lq)
    return (1e33 * m6**(1.0/3.0) + 3e31 + 7e32) * (G / 1e-2)

def p_sne(Gsne):
    return 2e33 * (Gsne / 1e-3)

def p_agn(lq, G, fB3):
    m6 = m6_of(lq)
    f = FEDD_C * fB3 * G * m6**FEDD_P
    fc = np.minimum(f, 1.0)
    return 1e35 * fc * f * m6

# Compression-limited branch (condition C4): the largest rate the clouds can
# sustain, from t_relax,MC = t_relax,*.  Nearly geometry-independent.
GAM_CM_C, GAM_CM_P = 2.169e-2, 0.466           # f_* = 0.22
GAM_CM_C_D = 2.158e-2                          # f_* = 0.18

def curve_angle(ax, x0, fn, dx=0.05):
    """Angle in degrees of log10(fn) at x0, measured in DISPLAY space so the
    label lies along the plotted line regardless of the axes aspect ratio."""
    import numpy as _np
    p0 = ax.transData.transform((x0 - dx, _np.log10(fn(x0 - dx))))
    p1 = ax.transData.transform((x0 + dx, _np.log10(fn(x0 + dx))))
    return _np.degrees(_np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))


def gamma_comp(lq, c=GAM_CM_C):
    return c * m6_of(lq)**GAM_CM_P

# Collision-capped branch, from C3 with the Cohn-Kulsrud bridged loss-cone flux
# (losscone.py).  Gamma_eq propto f_*^{3/2} no longer holds, so the disk-fed
# case carries its own fitted coefficients.
GAM_CK_C,   GAM_CK_P   = 1.419e-2, -0.844      # f_* = 0.22 (window midpoint)
GAM_CK_C_D, GAM_CK_P_D = 8.63e-3,  -0.853      # f_* = 0.18 (window lower edge)
GAM_CK_HI,  GAM_CK_HI_P = 1.933e-2, -0.860     # f_* = 0.25 (window upper edge)
GAM_CK_LO,  GAM_CK_LO_P = 8.63e-3, -0.853

def gamma_eq(lq):
    """Window midpoint f_* = 0.22, capped by C4 compression."""
    return np.minimum(GAM_CK_C * m6_of(lq)**GAM_CK_P, gamma_comp(lq, GAM_CM_C))

F_DISK = 0.18                                  # lower edge of the allowed window
def gamma_eq_disk(lq):
    return np.minimum(GAM_CK_C_D * m6_of(lq)**GAM_CK_P_D, gamma_comp(lq, GAM_CM_C_D))

def gamma_eq_band(lq):
    """f_* = 0.18 - 0.25 envelope: the range allowed by Sec 5.1."""
    lo = np.minimum(GAM_CK_LO * m6_of(lq)**GAM_CK_LO_P, gamma_comp(lq, GAM_CM_C_D))
    hi = np.minimum(GAM_CK_HI * m6_of(lq)**GAM_CK_HI_P, gamma_comp(lq, GAM_CM_C))
    return lo, hi

LQ_CRIT_S = 5.863                              # m6 = 0.730, C4 crossover, f_*=0.22
LQ_CRIT_D = 5.700                              # m6 = 0.501, C4 crossover, f_*=0.18

# colors matched to the Mathematica originals
ORANGE = (1.0, 0.5, 0.0)
MAGENTA = (0.66, 0.0, 0.66)      # Darker@Magenta
DGREEN = (0.0, 0.52, 0.0)        # Darker@Green
DPURPLE = (0.33, 0.0, 0.33)      # Darker@Purple
GOLD = (1.0, 0.80, 0.10)         # Blend[Yellow,Orange]
BLUE = (0.0, 0.45, 1.0)          # Blend[Cyan,Blue]

LQ_SWALLOW = 8.5

# =====================================================================
# Figure 1: mom-rates.pdf
# =====================================================================
def fig_mom_rates():
    lq = np.linspace(3.0, 10.3, 900)
    lG = np.linspace(-5.6, 1.9, 900)
    LQ, LG = np.meshgrid(lq, lG)
    Q, G = LQ, 10.0**LG

    PT = p_tde(Q, G)
    lfBs = [-4, -3, -2, -1, 0]      # log10 f_B
    lGsnes = [-5, -4, -3, -2, -1]   # log10 Gamma_SNe

    fig, ax = plt.subplots(figsize=(5.4, 5.4))

    # purple stacked regions: TDE beats SNe (and beats AGN at weakest coupling f_B=1e-4)
    agn_weak = p_agn(Q, G, 10.0**(lfBs[0] + 3))
    for lGs in lGsnes[:-1]:
        mask = (PT > p_sne(10.0**lGs)) & (PT > agn_weak) & (LQ <= LQ_SWALLOW)
        ax.contourf(LQ, LG, mask.astype(float), levels=[0.5, 1.5],
                    colors=[to_rgba(MAGENTA, 1.0/len(lGsnes))])
    # orange stacked regions: TDE beats AGN (and beats SNe at lowest SN rate)
    sne_weak = p_sne(10.0**lGsnes[0])
    for lfB in lfBs:
        mask = (PT > p_agn(Q, G, 10.0**(lfB + 3))) & (PT > sne_weak) & (LQ <= LQ_SWALLOW)
        ax.contourf(LQ, LG, mask.astype(float), levels=[0.5, 1.5],
                    colors=[to_rgba(ORANGE, 1.0/len(lfBs))])

    # boundary contour lines
    for i, lGs in enumerate(lGsnes):
        Z = np.where((LQ <= LQ_SWALLOW), PT - p_sne(10.0**lGs), np.nan)
        ax.contour(LQ, LG, Z, levels=[0.0],
                   colors=[MAGENTA], linewidths=0.75 + 0.25*i)
    for i, lfB in enumerate(lfBs):
        Z = np.where((LQ <= LQ_SWALLOW) & (PT > sne_weak),
                     PT - p_agn(Q, G, 10.0**(lfB + 3)), np.nan)
        ax.contour(LQ, LG, Z, levels=[0.0],
                   colors=[ORANGE], linewidths=0.75 + 0.25*i)

    ax.axvline(LQ_SWALLOW, color='k', ls=(0, (5, 3)), lw=1.5)
    ax.text(8.62, -3.6, r'$r_{t} \lesssim r_{s}$', rotation=-90, fontsize=15,
            fontweight='bold', ha='left', va='center')

    # line labels (rotated, at right edge like the original)
    # SNe boundaries slope downward: Gamma = 1e-2 * p_sne / (1e33 m6^{1/3} + 7.3e32)
    def lG_sne_bound(lqv, lGs):
        return np.log10(1e-2 * p_sne(10.0**lGs) / (1e33 * m6_of(lqv)**(1/3.) + 7.3e32))
    ax.text(7.0, lG_sne_bound(7.0, -1) + 0.16, r'$\Gamma_{\rm SNe} = 10^{-1}\ {\rm yr}^{-1}$',
            color=MAGENTA, fontsize=11, rotation=-13, ha='center')
    ax.text(7.0, lG_sne_bound(7.0, -5) + 0.16, r'$\Gamma_{\rm SNe} = 10^{-5}\ {\rm yr}^{-1}$',
            color=MAGENTA, fontsize=11, rotation=-13, ha='center')
    # f_B bands: the AGN boundary turns vertical once f_Edd saturates, so the
    # family is annotated in the empty strip beyond the swallowing limit.
    ax.text(8.32, 1.52, r'$f_{\rm B}: 10^{-4} \rightarrow 10^{0}$', color=ORANGE,
            fontsize=12, ha='right', va='center')
    ax.text(8.32, 1.25, '(left to right)', color=ORANGE, fontsize=10,
            ha='right', va='center')
    ax.text(8.32, -4.95, r'$\Gamma_{\rm SNe}: 10^{-5} \rightarrow 10^{-1}\ {\rm yr}^{-1}$',
            color=MAGENTA, fontsize=11, ha='right', va='center')
    ax.text(8.32, -5.24, '(bottom to top)', color=MAGENTA, fontsize=10,
            ha='right', va='center')

    ax.set_xlim(3, 10.3)
    ax.set_ylim(-5.6, 1.9)
    ax.set_xlabel(r'${\rm Log}_{10}[M_{h}/M_{\odot}]$', fontsize=16)
    ax.set_ylabel(r'${\rm Log}_{10}\,\Gamma_{\rm TDE}\ ({\rm yr}^{-1})$', fontsize=16)
    ax.tick_params(labelsize=13)
    fig.tight_layout(pad=0.4)
    fig.savefig(OUTDIR+'/mom-rates.pdf')
    plt.close(fig)

# =====================================================================
# Figure 2: dominance.pdf   (f_B = 1e-3, Gamma_SNe = 1e-4)
# =====================================================================
def fig_dominance():
    lq = np.linspace(5.0, 9.0, 1000)
    lG = np.linspace(-4.0, 1.0, 1000)
    LQ, LG = np.meshgrid(lq, lG)
    Q, G = LQ, 10.0**LG
    fB3, Gsne = 1.0, 1e-4

    PT, PA, PS = p_tde(Q, G), p_agn(Q, G, fB3), np.full_like(Q, p_sne(Gsne))

    tde = (PT > PA) & (PT > PS) & (LQ <= LQ_SWALLOW)
    agn = ((PA >= PT) & (PA > PS) & (LQ <= LQ_SWALLOW)) | ((LQ > LQ_SWALLOW) & (PA > PS))
    sne = ~(tde | agn)

    fig, ax = plt.subplots(figsize=(5.4, 5.4))
    for mask, col in [(tde, DGREEN), (agn, GOLD), (sne, DPURPLE)]:
        ax.contourf(LQ, LG, mask.astype(float), levels=[0.5, 1.5],
                    colors=[to_rgba(col, 0.55)])

    ax.axvline(LQ_SWALLOW, color='k', ls=(0, (5, 3)), lw=1.5)

    # equilibrium line, clipped to TDE-dominant region
    lqe = np.linspace(5.0, LQ_SWALLOW, 400)
    lGe = np.log10(gamma_eq(lqe))
    # shaded band = the flattening window of Sec 5.1, f_* = 0.18 - 0.25
    _lo, _hi = gamma_eq_band(lqe)
    ax.fill_between(lqe, np.log10(_lo), np.log10(_hi), color=BLUE, alpha=0.16, lw=0)
    for _y in (_lo, _hi):
        ax.plot(lqe, np.log10(_y), color=BLUE, lw=1.1, ls=(0, (1, 2)), alpha=0.9,
                zorder=4)
    on = p_tde(lqe, 10**lGe) > p_agn(lqe, 10**lGe, fB3)
    ax.plot(lqe[on], lGe[on], color=BLUE, lw=2.6, solid_capstyle='round')
    ax.plot([LQ_CRIT_S], [np.log10(gamma_eq(LQ_CRIT_S))], 'o', color=BLUE, ms=5,
            mec='w', mew=0.8, zorder=6)
    ax.plot([LQ_CRIT_D], [np.log10(gamma_eq_disk(LQ_CRIT_D))], 'o', color=BLUE, ms=5,
            mec='w', mew=0.8, zorder=6)
    _hiF = lambda q: gamma_eq_band(q)[1]
    ax.text(7.25, np.log10(_hiF(7.25)) + 0.20, r'$f_{\ast}=0.25$', color=BLUE,
            fontsize=10, fontweight='bold', rotation=curve_angle(ax, 7.25, _hiF),
            rotation_mode='anchor', ha='center')
    ax.text(6.15, np.log10(gamma_eq(6.15)) + 0.24,
            r'$\Gamma_{\rm eq}$, $f_{\ast}=0.22$', color=BLUE, fontsize=10.5,
            fontweight='bold', rotation=curve_angle(ax, 6.15, gamma_eq),
            rotation_mode='anchor', ha='center')
    ax.text(7.55, np.log10(gamma_eq_disk(7.55)) - 0.30,
            r'$f_{\ast}=0.18$', color=BLUE, fontsize=10,
            fontweight='bold', rotation=curve_angle(ax, 7.55, gamma_eq_disk),
            rotation_mode='anchor', ha='center')

    ax.text(5.60, -2.58, 'TDEs Dominate\nFeedback', color=DGREEN, fontsize=12,
            fontweight='bold', ha='center', va='center', linespacing=1.0)
    ax.text(7.15, 0.72, 'AGN Dominates Feedback', color=(0.72, 0.50, 0.0),
            fontsize=17, fontweight='bold', ha='center', va='center')
    ax.text(8.75, -1.5, 'MS Stars Swallowed Whole', rotation=-90, fontsize=18,
            fontweight='bold', ha='center', va='center')
    ax.text(6.0, -3.72, 'SNe Dominate Feedback', color=DPURPLE, fontsize=14,
            fontweight='bold', ha='center', va='center')

    ax.set_xlim(5, 9)
    ax.set_ylim(-4, 1)
    ax.set_xlabel(r'${\rm Log}_{10}[M_{h}/M_{\odot}]$', fontsize=16)
    ax.set_ylabel(r'${\rm Log}_{10}\,\Gamma_{\rm TDE}\ ({\rm yr}^{-1})$', fontsize=16)
    ax.tick_params(labelsize=13)
    fig.tight_layout(pad=0.4)
    fig.savefig(OUTDIR+'/dominance.pdf')
    plt.close(fig)

# =====================================================================
# Figure 3: domratio.pdf   ratio of TDE momentum to all others
# =====================================================================
def fig_domratio():
    lq = np.linspace(5.0, 9.0, 1000)
    lG = np.linspace(-4.0, 1.0, 1000)
    LQ, LG = np.meshgrid(lq, lG)
    Q, G = LQ, 10.0**LG
    fB3, Gsne = 1.0, 1e-4

    ratio = p_tde(Q, G) / (p_agn(Q, G, fB3) + p_sne(Gsne))
    lratio = np.log10(ratio)
    lratio_m = np.where(LQ <= LQ_SWALLOW, lratio, np.nan)

    fig, ax = plt.subplots(figsize=(5.4, 5.4))
    levels = np.arange(-1.5, 1.51, 0.5)
    # yellow (low) -> green (high)
    n = len(levels) + 1
    cols = [tuple(np.array((1.0, 0.90, 0.42)) * (1 - t) + np.array((0.15, 0.60, 0.15)) * t)
            for t in np.linspace(0, 1, n)]
    ax.contourf(LQ, LG, lratio_m, levels=[-10, *levels, 10], colors=cols)
    cs = ax.contour(LQ, LG, lratio_m, levels=levels, colors='k', linewidths=0.55)

    def fmt(v):
        if abs(v) < 1e-9:
            return '1'
        if abs(v - round(v)) < 1e-9:
            e = int(round(v))
            return r'$10^{%d}$' % e if e > 0 else r'$10^{%d}$' % e
        return r'$10^{%.1f}$' % v
    ax.clabel(cs, fmt=fmt, fontsize=9.0, inline=True, inline_spacing=6)
    # blank out the swallowed-whole zone
    ax.fill_betweenx([-4, 1], LQ_SWALLOW, 9.0, color='white', zorder=2)
    ax.axvline(LQ_SWALLOW, color='k', ls=(0, (5, 3)), lw=1.5, zorder=3)

    lqe = np.linspace(5.0, LQ_SWALLOW, 400)
    lGe = np.log10(gamma_eq(lqe))
    _lo, _hi = gamma_eq_band(lqe)
    ax.fill_between(lqe, np.log10(_lo), np.log10(_hi), color=BLUE, alpha=0.16,
                    lw=0, zorder=4)
    for _y in (_lo, _hi):
        ax.plot(lqe, np.log10(_y), color=BLUE, lw=1.1, ls=(0, (1, 2)), alpha=0.9,
                zorder=4)
    ax.plot(lqe, lGe, color=BLUE, lw=2.6, solid_capstyle='round', zorder=4)
    ax.text(6.2, np.log10(gamma_eq(6.2)) + 0.26, r'$\Gamma_{\rm eq}$, $f_{\ast}=0.22$',
            color=BLUE, fontsize=11.5, fontweight='bold',
            rotation=curve_angle(ax, 6.2, gamma_eq), rotation_mode='anchor', ha='center',
            zorder=6, bbox=dict(fc='white', ec='none', alpha=0.8, pad=1.5))
    ax.text(7.2, np.log10(gamma_eq_disk(7.2)) - 0.28, r'$f_{\ast}=0.18$',
            color=BLUE, fontsize=10, fontweight='bold',
            rotation=curve_angle(ax, 7.2, gamma_eq_disk), rotation_mode='anchor', ha='center',
            zorder=6, bbox=dict(fc='white', ec='none', alpha=0.8, pad=1.5))

    ax.text(8.75, -1.5, 'MS Stars Swallowed Whole', rotation=-90, fontsize=18,
            fontweight='bold', ha='center', va='center', zorder=5)
    ax.text(0.035, 0.035, r'contours: $\dot{p}_{\rm TDE}/\dot{p}_{\rm other}$',
            transform=ax.transAxes, fontsize=12, ha='left', va='bottom', zorder=5,
            bbox=dict(fc='white', ec='0.6', lw=0.6, pad=2.5))

    ax.set_xlim(5, 9)
    ax.set_ylim(-4, 1)
    ax.set_xlabel(r'${\rm Log}_{10}[M_{h}/M_{\odot}]$', fontsize=16)
    ax.set_ylabel(r'${\rm Log}_{10}\,\Gamma_{\rm TDE}\ ({\rm yr}^{-1})$', fontsize=16)
    ax.tick_params(labelsize=13)
    fig.tight_layout(pad=0.4)
    fig.savefig(OUTDIR+'/domratio.pdf')
    plt.close(fig)

if __name__ == '__main__':
    fig_mom_rates()
    fig_dominance()
    fig_domratio()
    print('done')
