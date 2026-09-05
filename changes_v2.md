# Paper revision log (v2)

Tracking manuscript changes from commit `5640c66` (`acknowledgements change`) until the revised version is submitted. Newest entries first.

Figures are named when they change. If an entry does not mention a figure, none was regenerated.

---

## 2026-09-05 — Consistency with the Cohn–Kulsrud equilibrium

Replace leftover monomial-rate coefficients (the superseded \(\Gamma \propto m_6^{4/3} n_0 \sigma_h^{-3}\) solve) with the Cohn–Kulsrud power-law fit already used in the engine-state tables, the flattening window, and `scripts/make_figures.py`. The published \(\Gamma_{\rm eq}\) intercepts are unchanged: they remain the four-decade fit over \(m_6 = 0.03\)–\(300\), not the live integrator at \(10^6 M_\odot\).

**Figures:** none regenerated. `mom-rates.pdf`, `dominance.pdf`, and `domratio.pdf` already use the window fit quoted in the text. `window.pdf` and `engine-diagram.pdf` are unchanged.

### Section 2 (tidal disruptions)

- \(L_{\rm UDR}\): \(3\times10^{40}\Gamma_{-2}M_{\rm h,6}^{1/3}\) → \(6\times10^{40}\Gamma_{-2}M_{\rm h,6}^{1/3}\), so that \(L_{\rm UDR} = \Gamma E_{\rm UDR}\) with \(E_{\rm UDR} = 1.8\times10^{50} M_{\rm h,6}^{1/3}\) from the engine section. \(\dot{p}_{\rm UDR}\) left at \(10^{33}\) (snowplow calibration, not \(E/v\)).

### Section 4.2 (compression criterion)

- Dropped \(\sigma_{\rm coll} \propto M_{\rm h,6}^{1/6}\) from the C4 preamble (that slope is the old monomial C3). The following paragraph already had the CK result: \(\sigma_{\rm coll} = 6.5\times10^{2}\,f_{\ast}^{0.66}\,{\rm km\,s}^{-1}\), mass-independent to within \(b = 0.001\)–\(0.044\).
- Unified the isotropic mass floor, previously quoted as \(6\times10^{6}\), \(8\times10^{6}\), and \(8.9\times10^{6}\,M_\odot\), to \(9\times10^{6}\,M_\odot\), matching \(M_{\rm h,6} \gtrsim 8.9\,f_{\ast}^{1.64}\) at \(f_{\ast}=1\). Disk-fed floor at \(f_{\ast}=0.2\) (\(6\times10^{5}\,M_\odot\)) and window floor at \(f_{\ast}=0.22\) (\(7.3\times10^{5}\,M_\odot\)) unchanged.

### Section 4.3 (isotropic engine state)

- Stated that C3 is an energy integral and that the displayed power laws follow after replacing \(\Gamma_{\rm TDE}(\sigma,M_{\rm h})\) with a fit over \(3\times10^{4}\)–\(3\times10^{8}\,M_\odot\). Dropped “exact” from “exact power laws”.
- Validity paragraph: isotropic mass floor \(8\times10^{6}\,M_\odot\) → \(9\times10^{6}\,M_\odot\).

### Section 4.4 (disk geometry)

- \(M_{\rm disk} = f_{\Omega} M_{\rm belt}\): \(4.8\times10^{5}(f_{\Omega}/0.2)M_{\rm h,6}^{1.14}\) → \(1.3\times10^{5}(f_{\Omega}/0.2)M_{\rm h,6}^{1.57}\), i.e. \(0.2\times\) the current isotropic belt (\(6.7\times10^{5} M_{\rm h,6}^{1.57}\)) rather than the old monomial belt.
- Rewrote the following sentence. The old “\(2.4\times10^{6}\,M_\odot\) is uncomfortably large” argument applied to the superseded spherical mass; the CK isotropic belt is already \(6.7\times10^{5}\,M_\odot\). Flattening is now a further reduction to \(\sim10^{5}\,M_\odot\), below the Garcia-Burillo median.
- Flattened \(\sigma(f_{\ast})\): \(400\,f_{\ast}^{1/2}M_{\rm h,6}^{1/6}\) (old C3) → \(6.5\times10^{2}\,f_{\ast}^{0.66}\,{\rm km\,s}^{-1}\) (CK, as in C4). At \(f_{\ast}=0.2\) this is \(\sim230\,{\rm km\,s}^{-1}\), matching the disk-fed block.
- Restricted \(\Gamma_{\rm eq} \propto (f_{\ast}/0.22)^{2.4}\) to \(0.18\lesssim f_{\ast}\lesssim 0.25\). The isotropic rate remains Equation (eq:eqrate), not this extrapolation (which would give \(0.53\,{\rm yr}^{-1}\) at \(f_{\ast}=1\) instead of \(0.42\)).
- Disk-fed C4 comparison: isotropic floor \(8\times10^{6}\,M_\odot\) → \(9\times10^{6}\,M_\odot\).
- Disk-fed \(\tau_{\rm relax,\ast}\): \(2.1\times10^{9} M_{\rm h,6}^{3/2}\,{\rm yr}\) (old monomial \(\sigma\propto M^{1/6}\)) → \(1.1\times10^{9} M_{\rm h,6}^{1.9}\,{\rm yr}\) (CK \(\sigma\simeq230 M_{\rm h,6}^{0.03}\), \(n_0\propto M_{\rm h,6}^{-1.80}\)).

### Section 6 (discussion and implications)

- CMZ analogy: \(\mathcal{M} \simeq 40\) → \(\mathcal{M} \simeq 30\) (CK disk/window Mach is \(32\)–\(33\)).
- AGN level: \(L_{\rm AGN}/L_{\rm Edd} \sim 0.06 f_{\rm B,-3}\) → \(\sim 0.05 f_{\rm B,-3}\) (window \(0.051\), disk-fed \(0.053\)).
- k+A disk mass: \(\sim7\times10^{5} M_{\rm h,6}^{1.14}\,M_\odot\) (old disk-fed \(M_{\rm gas}\)) → \(\sim4\times10^{5} M_{\rm h,6}^{1.5}\,M_\odot\) (window \(3.9\times10^{5}\), disk-fed \(3.8\times10^{5} M_{\rm h,6}^{1.48}\)).
- LRD comparison: \(a_{\rm h} \sim 0.05\)–\(0.3\,{\rm pc}\) (old monomial range) → \(\sim 0.02\)–\(0.2\,{\rm pc}\) (CK isotropic \(0.021\,{\rm pc}\), disk-fed \(0.17\,{\rm pc}\)).
- LRD rate window: \(8.6\times10^{-3}\)–\(2.1\times10^{-2}\,{\rm yr}^{-1}\) → \(8.6\times10^{-3}\)–\(1.9\times10^{-2}\,{\rm yr}^{-1}\) (\(f_{\ast}=0.25\) fit is \(1.933\times10^{-2}\)).

### Left unchanged on purpose

- Abstract \(\Gamma = 1.4\times10^{-2} M^{-0.84}\), window midpoint, disk-fed table, isotropic table, \(\dot{M}\) ratios, \(A_V\simeq50\), \(\mathcal{F}_{\rm Edd}\) / \(\langle L_{\rm TE}\rangle\), \(\dot{p}_{\rm TDE}/\dot{p}_{\rm AGN}=9.4\).
- \(\sigma_{\rm comp} = 2.8\times10^{2} M_{\rm h,6}^{0.40}\) and the prose values \(275\) / \(272\) / \(271\).
- Caveat \(q_{\ast}\) and \(r_{\rm crit}\) (not leftover monomial coefficients; not re-derived in this pass).

---
