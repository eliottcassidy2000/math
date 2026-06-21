# The LRC sector-discrepancy is a sum of two cycle-graph resistances

*kind-pasteur, 2026-06-21 (THREAD A, the owner's modular-forms lead)*

## The convergence

Three independent framings of the L7 obstruction collapse onto one object.

The LRC(14) sector-route bottleneck (L7) reduced to a single sharp cell-discrepancy
of the `(q,p)` torus geodesic against the 7-sector grid (HYP-2739). That discrepancy is
*residue-only mod the apex prime*: `D_{p,q} = G_P(p,q)/(P p q)`. The question the owner
posed was modular: is `G_P` a value of a level-`P` modular form, a zonal/spherical
function on `P^1(F_P)`, a dessin d'enfant?

The answer is cleaner than any of those, and it is the same object three ways:

```
   G_P(p,q) = [ 2 A B (P-A)(P-B) + 2 C (P-C) ] / P,
   A = ||p||_P,  B = ||q||_P,  C = ||pq||_P,   ||x|| = min(x mod P, P - x mod P).
```

Each leg `g(t) = 2 t (P-t)` is simultaneously:

1. **A Bernoulli value.** `g(t) = -2P^2 B_2(t/P) + P^2/3`. The second Bernoulli
   polynomial is the seed of the weight-2 Eisenstein series `E_2`. There is no
   *holomorphic* weight-2 form — and that is exactly right: the discrepancy is an
   `L1` (absolute) quantity, and the holomorphic/signed object is its degree-1
   Fourier shadow (the classical Dedekind sum), the very split that HYP-2743 and
   MISTAKE-082 already flagged as why the cover needs the absolute discrepancy.

2. **A cycle resistance.** `t(P-t)/P = R_eff(0,t)` on the cycle graph `C_P` — the
   effective electrical resistance between two vertices `t` apart, equivalently the
   discrete Green's function of `C_P` (Chung–Yau). Its Fourier coefficients are
   `1/(2 - 2cos(2 pi k/P))`, the inverse Laplacian eigenvalues. The parabola that
   governs a number-theoretic discrepancy *is the same parabola* that governs current
   flow around a ring of `P` resistors.

3. **A discrepancy.** Literally `sum_j |P c_{0j} - pq|`, the `L1` defect of a sawtooth.

That a torus-geodesic discrepancy, a Bernoulli/Eisenstein value, and a graph Green's
function are the *same finite parabola* is the reflection. The cycle `C_P` is the apex
prime made into a graph; the LRC sectors are its vertices; the geodesic's mismatch with
them is the resistance you'd measure across the ring.

## What survives and what washes out

The symmetry is just as clean and now *proved for every prime*, not conjectured for 7:
the stabilizer of `G_P` inside `PGL_2(F_P)` is exactly the Klein four-group
`<z -> -z, z -> 1/z>` (order 4) — the `±` hyperelliptic involution times modular
inversion. Doubling `z -> 2z`, the quadratic-residue / cubic order-3 structure of the
apex prime, is **not** a symmetry at any `P >= 5`. The order-2 part of the prime
survives into the magnitude of the obstruction; the order-3 part does not. This is the
uniform-in-`P` form of "QR washes out" (HYP-2657/2692).

## The honest correction it forced

The earlier framing said `G` is "a function on `P^1(F_7)`." It is not. Every interior
slope `p/q` is `G`-multivalued, because the third leg `C = ||pq||` is a **product**
coordinate — it lives on the multiplicative group, invisible to the slope. `G_P` is a
class function on the *pair* torus `(Z/P)^2` modulo `<±, swap>`, finer than the
projective line. The Möbius `<±, inv>` group captured the symmetry but mis-stated the
domain. The geometry was right one level up: the data is the *multiplication table* of
`Z/P` folded by `±1`, which is exactly the cusp/elliptic-point bookkeeping of `X(P)` —
for `P = 7`, the Klein quartic. `G_7` is the hyperelliptic (order-2) shadow of `X(7)`,
not a full level-7 form. The order-3 PSL(2,7) structure of the Klein quartic is precisely
the part that does not act.

## Why this matters beyond LRC

The same parabola is the building block of: Dedekind sums (modular), cycle resistances
(spectral graph theory), the variance of a discrete uniform (probability), and the L1
discrepancy of `{k alpha}` (equidistribution). The LRC apex obstruction sits at their
intersection. When a single degree-2 object keeps reappearing across number theory,
electrical networks, and discrepancy, the lesson is that the "hard" L7 constant was never
analytic — it was the resistance of a ring, computable by Kirchhoff, bounded by `B_2`.
The triangle foundation of this project says everything is a side of the staircase; here
the side is literally a resistor.

-> HYP-2745, HYP-2743, HYP-2742, HYP-2739, MISTAKE-082, OPEN-Q-108.
