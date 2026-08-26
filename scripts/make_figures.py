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
FEDD_C = 1.62                                  # f_Edd = FEDD_C * fB3 * Gamma * m6^{-0.0447}
FEDD_P = -0.0447

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

def gamma_eq(lq):
    return GAM_EQ_C * m6_of(lq)**GAM_EQ_P

F_DISK = 0.2                                   # fiducial disk aspect ratio
def gamma_eq_disk(lq):
    """Disk-fed equilibrium: Gamma_eq propto f_*^{3/2}."""
    return GAM_EQ_C * (F_DISK**1.5) * m6_of(lq)**GAM_EQ_P

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
    ax.text(9.42, 1.30, r'$f_{\rm B}: 10^{-4} \rightarrow 10^{0}$', color=ORANGE,
            fontsize=12, ha='center', va='center')
    ax.text(9.42, 1.03, '(left to right)', color=ORANGE, fontsize=10,
            ha='center', va='center')
    ax.text(9.42, -4.85, r'$\Gamma_{\rm SNe}: 10^{-5} \rightarrow 10^{-1}\ {\rm yr}^{-1}$',
            color=MAGENTA, fontsize=11, ha='center', va='center')
    ax.text(9.42, -5.12, '(bottom to top)', color=MAGENTA, fontsize=10,
            ha='center', va='center')

    ax.set_xlim(3, 10.3)
    ax.set_ylim(-5.6, 1.9)
    ax.set_xlabel(r'${\rm Log}_{10}[M_{h}/M_{\odot}]$', fontsize=16)
    ax.set_ylabel(r'${\rm Log}_{10}\,\Gamma_{\rm TDE}\ ({\rm yr}^{-1})$', fontsize=16)
    ax.tick_params(labelsize=13)
    fig.tight_layout(pad=0.4)
    fig.savefig('/home/claude/work/figs/mom-rates.pdf')
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
    lGd = np.log10(gamma_eq_disk(lqe))
    ax.fill_between(lqe, lGd, lGe, color=BLUE, alpha=0.16, lw=0)
    on = p_tde(lqe, 10**lGe) > p_agn(lqe, 10**lGe, fB3)
    ax.plot(lqe[on], lGe[on], color=BLUE, lw=2.6, solid_capstyle='round')
    ax.plot(lqe, lGd, color=BLUE, lw=2.2, ls=(0, (6, 3)), solid_capstyle='round')

    ax.text(6.05, np.log10(gamma_eq(6.05)) + 0.22,
            r'$\Gamma_{\rm eq}$: spherical, $f_{\ast}=1$', color=BLUE, fontsize=11,
            fontweight='bold', rotation=-6.5, ha='center')
    ax.text(6.15, np.log10(gamma_eq_disk(6.15)) - 0.36,
            r'$\Gamma_{\rm eq}$: disk-fed, $f_{\ast}=0.2$', color=BLUE, fontsize=11,
            fontweight='bold', rotation=-6.5, ha='center')

    ax.text(7.55, -3.05, 'TDEs Dominate\nFeedback', color=DGREEN, fontsize=16,
            fontweight='bold', ha='center', va='center', linespacing=1.0)
    ax.text(6.9, 0.62, 'AGN Dominates Feedback', color=(0.72, 0.50, 0.0),
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
    fig.savefig('/home/claude/work/figs/dominance.pdf')
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
    lGd = np.log10(gamma_eq_disk(lqe))
    ax.fill_between(lqe, lGd, lGe, color=BLUE, alpha=0.16, lw=0, zorder=4)
    ax.plot(lqe, lGe, color=BLUE, lw=2.6, solid_capstyle='round', zorder=4)
    ax.plot(lqe, lGd, color=BLUE, lw=2.2, ls=(0, (6, 3)), zorder=4)
    ax.text(6.5, np.log10(gamma_eq(6.5)) + 0.26, r'$\Gamma_{\rm eq}$ (spherical)',
            color=BLUE, fontsize=11.5, fontweight='bold', rotation=-6.5, ha='center',
            zorder=6, bbox=dict(fc='white', ec='none', alpha=0.8, pad=1.5))
    ax.text(6.6, np.log10(gamma_eq_disk(6.6)) - 0.30, r'$\Gamma_{\rm eq}$ (disk-fed)',
            color=BLUE, fontsize=11.5, fontweight='bold', rotation=-6.5, ha='center',
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
    fig.savefig('/home/claude/work/figs/domratio.pdf')
    plt.close(fig)

if __name__ == '__main__':
    fig_mom_rates()
    fig_dominance()
    fig_domratio()
    print('done')
