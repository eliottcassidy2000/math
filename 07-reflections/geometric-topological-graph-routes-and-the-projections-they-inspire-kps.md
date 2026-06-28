# The geometry of LRC(14): five lenses, the routes they open, and the proof technologies they project onto

*kind-pasteur-2026-06-28-S31ax. The owner asked to think about the underlying topology/geometry/graphs, how
they can be leveraged for visual or structural proof routes, and how those inspire other types of routes.
This maps the geometry of LRC(14) onto five lenses, reads a proof route off each, and shows that the named
proof technologies (equidistribution, moment-LP, Borsuk–Ulam, flatness, spectral) are **projections of one
geometric object** onto different machinery. The visual (the AP safety function equioscillating at the six
units mod 14) is that object made manifest.*

## The five geometric lenses (same object, different category)
| lens | the object | LRC(14) becomes | category |
|---|---|---|---|
| **1. Torus + geodesic** | closed geodesic `γ_v(t)=(v₁t,…,v₁₃t)` on `T¹³`, class `v∈H₁=ℤ¹³`; safe box `B=[1/14,13/14]¹³` | does `γ_v` enter `B`? | differential topology |
| **2. Circle + covering** | danger arcs `D_i={‖v_i t‖<1/14}` on `S¹`, total measure `13/7≈1.857>1` | do the arcs miss a point? | measure / combinatorics |
| **3. Polyhedron + lattice** | the lonely-runner polyhedron `P_v` | does `P_v` contain an integer point? | geometry of numbers |
| **4. Heptagon + D₇** | 7 sectors = heptagon, `14=|D₇|`; witnesses `Φ₁₄` = 3 antipodal pairs | odd-degree free `ℤ/2` (HYP-3241) | discrete / algebraic topology |
| **5. Graph** | binding graph at a witness; danger-arc overlap graph | is there an independent (uncovered) point? | graph theory |

The visual is lens 2/4 fused: `S(t)=min_i‖v_i t‖` (the depth of `γ_v` into `B`), touching `1/14` at the six
units `Φ₁₄`, in three antipodal pairs mirror-symmetric about `t=1/2`.

## A route off each lens — and the technology it projects onto
**1. Geodesic → equidistribution.** The depth `S(t)` is how far `γ_v` reaches into the box. The geodesic is a
closed 1-subtorus; its visiting measure on `T¹³` is the pushforward of Lebesgue on the `t`-circle. "Does it
reach depth `1/14`?" projects onto **Weyl/Erdős–Turán equidistribution** (mac-mini's bulk, my Fejér HYP-3218).
The geometry says *why* the analytic tool is the right one: the box is a positive-measure target and the
geodesic equidistributes — measure is the natural certificate **off the core** (Vitali wall, HYP-3237).

**2. Covering → moment-LP / circle method.** Arcs of total mass `1.857` must still miss a point, so the
**overlap** structure is everything. Inclusion–exclusion of the overlaps is exactly the **moment-LP**
(codex's `q_t`, THM-534/537): `meas(gap)=Σ(-1)^j S_j`. The geometry says the obstruction is the alternating
(even/odd) overlap sum — the positive/negative duality (HYP-3220) is the IE sign.

**3. Polyhedron → flatness / geometry of numbers.** No integer point ⟹ `P_v` is **flat** (thin) in some
lattice direction (Khinchin). The AP is the flatness-extremal `v`; the **de Moivre/Φ₇ direction is the flat
direction**. Projects onto the Beck–Hosten–Schymura / covering-radius route (HYP-2764, HYP-3215).

**4. Heptagon/D₇ → Borsuk–Ulam (topology).** `14=|D₇|`; the complement is the order-2 reflection, a free
`ℤ/2` for `p≡3 mod 4`; the witnesses are `(p−1)/2=3` antipodal pairs. Projects onto the **degree/index**
certificate (HYP-3239/3241): the odd index forces the witness — the *construction* analogue of the measure
bound.

**5. Graph → Chebyshev extremality (the new inspired route).** At each witness the **binders form a perfect
matching of antipodal pairs** `{1,13},{3,11},{5,9}` — a 3-edge matching on `Φ₁₄`. The safety function
**equioscillates** at `2·3=6` points (the visual). By the **Chebyshev equioscillation theorem**, a function
that equioscillates between its extremes at the maximal number of points is the **minimax** — so the AP, which
equioscillates at `φ(14)=6` points all pinned to `1/14`, is the candidate **M-minimizer**, and `M(any) ≥
M(AP)=1/14` would follow from a uniqueness/sufficiency argument. This is the **de Moivre = Chebyshev**
identity (HYP-3212) read structurally: *equioscillation count `= φ(2p)` is the signature of the extremal
config*, and the proof target is "equioscillation ⟹ global minimizer."

## The cross-inspiration (the dependency web)
```
                    γ_v depth S(t)  ── the visual (equioscillation at Φ₁₄) ──
                   /            |               |                 \
   equidistribution   moment-LP/IE      Borsuk–Ulam index     Chebyshev minimax
   (bulk, measure)    (overlap sign)    (core, degree)        (AP = M-minimizer)
        |                  |                  |                      |
   Erdős–Turán/Fejér   even/odd parity   free ℤ/2, (p-1)/2     de Moivre = Chebyshev
   HYP-3218            HYP-3220          HYP-3241              HYP-3212
```
Reading across: the **bulk** wants measure (geodesic equidistributes), the **core** wants degree (the
antipodal `ℤ/2`), and the **boundary between them** is the equioscillation — where the moment-LP is tight and
the Chebyshev minimax is attained. Vitali wall = the seam (HYP-3237). The visual sits exactly on the seam:
the six peaks are the moment-LP's tight constraints, the Borsuk–Ulam's antipodal pairs, and the Chebyshev's
alternation points — **all the same six points `Φ₁₄`**.

## Why this is useful (the inspired non-geometric routes)
- **Approximation theory (NEW here):** treat `M(v)` as a minimax over `t`; the Chebyshev/Markov machinery
  (equioscillation ⟹ extremal, Markov brothers' inequality for the derivative at the binders) becomes
  available. The AP's six-fold equioscillation is the entry point. *Action:* prove "any `v` whose `S(t)` does
  not reach `1/14` cannot equioscillate at 6 points," contradicting a degree bound on `S`.
- **Spectral/graph:** the binding matching (3 antipodal edges) is the adjacency of a `6`-cycle / the `D₇`
  representation; its spectrum is the de Moivre angles. A spectral lower bound on `M` via the binder graph's
  eigenvalues. *Action:* relate the smallest binder-graph eigenvalue to `M` (Perron/Carathéodory, HYP-3201).
- **Geometry of numbers:** the flatness route gives a *finite* search (the flat direction is cyclotomic), a
  concrete alternative to the moment-LP for the core.

## Net
- **FIVE LENSES** (torus-geodesic, circle-covering, polyhedron-lattice, heptagon-`D₇`, graph) are one object;
  each projects LRC(14) onto a different proof technology (equidistribution / moment-LP / flatness /
  Borsuk–Ulam / Chebyshev).
- **THE VISUAL = THE SEAM:** the six `Φ₁₄` equioscillation peaks are simultaneously the moment-LP tight
  constraints, the Borsuk–Ulam antipodal pairs, and the Chebyshev alternation points.
- **NEW INSPIRED ROUTE:** approximation theory — the AP equioscillates at `φ(2p)` points, the signature of
  the minimax; target "equioscillation ⟹ global `M`-minimizer." Plus a spectral binder-graph lower bound.

→ HYP-3242 (this), HYP-3241 (Borsuk–Ulam index), HYP-3237 (Vitali wall), HYP-3220 (parity), HYP-3218
(Fejér/equidist), HYP-3212 (de Moivre=Chebyshev), HYP-3201 (Carathéodory), HYP-2764 (zonotope/covering
radius), OPEN-Q-108.
