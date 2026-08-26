# A Self-Sustaining Black Hole Engine Powered by the Tidal Disruptions of Stars

James Guillochon & Abraham Loeb

Manuscript source, figures, and solver scripts for a paper proposing that tidal disruption
events can sustain their own elevated rate. The unbound debris of each disruption drives a
supernova-like remnant into the molecular clouds surrounding the nucleus. That momentum keeps
the clouds turbulent and dense, the dense clouds shorten the relaxation time of the nuclear
star cluster, and the compressed cluster disrupts stars faster — closing a feedback loop. The
loop is not runaway: once the cluster is dense enough, star-star collisions begin removing
stars from the loss cone as fast as they can be fed into it, which caps the rate.

The first draft dates from around 2018. The steady-state system was never solved at the time;
it is solved here in closed form.

## Repository layout

```
paper/     LaTeX source, bibliography, figures, compiled PDF
scripts/   equilibrium solvers and figure generation
```

| File | What it is |
|---|---|
| `paper/engine.tex` | Manuscript (AASTeX, `preprint2`) |
| `paper/engineNotes.bib` | Bibliography |
| `CHANGES.md` | What changed relative to the original draft, and why |
| `REFERENCES_AUDIT.md` | Provenance and verification status of every bibliography entry |

## Building the paper

```sh
cd paper
pdflatex engine && bibtex engine && pdflatex engine && pdflatex engine
```

It additionally needs `ulem` and `soul`, which live in TeX Live's `texlive-plain-generic`.

## The code

Every equilibrium condition is monomial in the unknowns, so the whole system reduces to a
linear solve in log space. The scripts therefore produce **exact closed-form power laws**
rather than numerical fits — there is no root-finding and nothing to converge.

```sh
pip install sympy matplotlib

python scripts/solve_equilibrium.py        # spherical belt of clouds
python scripts/solve_equilibrium_disk.py   # flattened disk geometry
python scripts/make_figures.py             # regenerates the three figures
```

`solve_equilibrium_disk.py` accepts three optional rational arguments: the disk aspect ratio
`f_Omega`, the fraction `b` of unbound remnants that intercept the disk, and the cusp
flattening `f_*`.

```sh
python scripts/solve_equilibrium_disk.py 1/5 1/5 1     # isotropic cusp, disk gas
python scripts/solve_equilibrium_disk.py 1/5 1   1/5   # disk-fed and beamed (fiducial)
```

The first of these reproduces the spherical solution exactly. That is not a coincidence: with
an isotropic star cluster, flattening the gas into a disk cancels out of both equilibrium
conditions, because the smaller number of clouds is compensated exactly by the smaller
fraction of remnants that strike them. Only the total gas mass changes.

Both scripts print self-consistency checks alongside the solution: heating against cooling,
remnant hits per cloud crossing time, and the Jeans ratio of the clouds.

## Main result

| | spherical belt | disk-fed, f = 0.2 |
|---|---|---|
| Γ<sub>eq</sub> | 3.7 × 10⁻² m₆<sup>−1/6</sup> yr⁻¹ | 3.3 × 10⁻³ m₆<sup>−1/6</sup> yr⁻¹ |
| σ | 396 m₆<sup>1/6</sup> km/s | 177 m₆<sup>1/6</sup> km/s |
| a<sub>h</sub> | 0.055 m₆<sup>2/3</sup> pc | 0.27 m₆<sup>2/3</sup> pc |
| cloud density | 8.0 × 10⁴ cm⁻³ | 2.4 × 10⁴ cm⁻³ |
| total gas mass | 2.4 × 10⁶ M<sub>⊙</sub> | 7.2 × 10⁵ M<sub>⊙</sub> |
| A<sub>V</sub> | ~51 | ~52 |

with m₆ = M<sub>h</sub>/10⁶ M<sub>⊙</sub>. More generally σ ∝ f<sup>1/2</sup> and
Γ<sub>eq</sub> ∝ f<sup>3/2</sup>, where f is the flattening of the stellar cusp.

The equilibrium rate lands tens to hundreds of times above the canonical 10⁻⁴ yr⁻¹, and the
clouds that emerge are close analogues of those in the Milky Way's Central Molecular Zone:
parsec-scale, ~10⁴–10⁵ cm⁻³, turbulently supported, and forming almost no stars. The total gas
mass is comparable to the circumnuclear molecular disks now resolved around nearby active
nuclei. The nucleus is dusty enough that sightlines through the gas are heavily obscured, so
the model predicts that the obscured fraction of disruptions should equal the covering factor
of the disk.

## A note on preparation

The 2026 revision was prepared with assistance from Anthropic's Claude models. All references
and edits were hand-checked by the authors.
