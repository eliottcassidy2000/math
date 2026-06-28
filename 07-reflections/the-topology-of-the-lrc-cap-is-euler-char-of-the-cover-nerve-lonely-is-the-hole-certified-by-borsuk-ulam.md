# The topology of the LRC: the cap is the Euler characteristic of the danger-cover nerve, the lonely point is the cover's HOLE, and n=14's certificate is the Borsuk-Ulam odd degree — and every other proof route is a shadow of this geometry

*mac-mini-2026-06-28-S78. The owner: think about the underlying topology/geometry/graphs and how they give
visual/structural proof routes, and how those inspire other routes. Two topological facts lock together — the
cap is the cover's Euler characteristic (mine, verified), and n=14's witness is a Borsuk-Ulam antipodal pair
(kps S31av, the heptagon D₇). Together they make the geometry the UNIFYING frame; the combinatorial/algebraic/
Diophantine routes are its shadows. Builds on [[14-is-the-heptagon-dihedral-group-borsuk-ulam-not-brouwer-kps]],
[[the-vitali-wall-brouwer-equioscillation-and-the-cyclotomic-core-construction]], [[pushes-and-pulls-on-the-hard-core-d7-unification-the-covering-witness-construction-and-the-imaginary-quadratic-pull]].*

## The geometric objects (what the LRC IS)
- **The circle / the heptagon:** the `7` sectors are the vertices of a regular heptagon `ζ₇^k=e^{2πik/7}`; its
  symmetry group is `D₇` (order `14 = n`, kps). The runners' positions `frac(s_i t)` live on this circle.
- **The torus / view-obstruction:** the runner phases `(s_1 t,…,s_{n-1}t)` trace a diagonal LINE on the
  `(n-1)`-torus; the danger set is the union of tubes around the coordinate hyperplanes. LRC = the diagonal
  ESCAPES the tubes (Cusick's view-obstruction).
- **The danger cover & its nerve:** the danger combs `D_p` (each `p` arcs of width `1/7`) COVER the circle (or
  not). The combinatorial overlap structure is the NERVE (the simplicial complex of intersections).

## STRUCTURAL ROUTE 1 — the cap IS the Euler characteristic of the cover nerve (VERIFIED)
The cap (lonely measure) is the measure-weighted Euler characteristic of the danger cover:
```
  cap = meas(lonely) = 1 - meas(∪ D_p) = Σ_{S⊆P} (-1)^|S| meas(∩_{p∈S} D_p)  =  χ_meas(nerve).
```
(Verified: for the binding minimizer `{1,5,7,8,9}`, this alternating sum `= 2243/5880 = cap_8`.) So:
> **The inclusion-exclusion IS the Čech/Euler computation; the lonely point is the cover's HOLE** (the reduced
> `H^0` / the uncovered region). LRC(14) ⟺ the danger cover has a hole ⟺ `χ_meas > 0`. The combinatorial cap and
> the topology are the same object: **the cap is a topological invariant of the cover.**

## STRUCTURAL ROUTE 2 — the witness is a Borsuk-Ulam antipodal pair (kps S31av), the D₇ symmetry
The heptagon reflection = the complement = mult-by-`−1`. For `7≡3 (mod 4)`, `−1` is a non-residue, so the
reflection is an **anti-automorphism** (orientation-reversing, self-converse) and the `ℤ₂` acts **FREELY**.
Hence the certificate is **Borsuk-Ulam, not Brouwer**: the witness is not a reflection-fixed point but an
**antipodal pair `(t*,−t*)`**, certified by the **odd degree** of the free reflection (an odd map `S¹→S¹` has
nonzero degree). The odd degree = the imaginary Gauss sum `i√7` = the negative trace = the odd power sums — one
datum. **The cap decomposes into `D₇`-isotypic parts** (Reynolds average): trivial = the mean (the Möbius/Euler
mean), sign = the Borsuk-Ulam obstruction, three 2-dim = the de Moivre/Fejér body. The obstruction is the
**sign-isotypic** component, certified by odd degree (topology), not SOS.

## STRUCTURAL ROUTE 3 — the graph: the nerve as the conflict graph, the hole as an independent set
The cover nerve is a GRAPH (vertices = combs, edges = overlaps). The lonely set = the region covered by NO
comb = an "independent" region. The tournament conflict graph `Ω` (odd cycles) and `H=I(Ω,2)` is the
project's graph invariant; the LRC cover nerve is its danger-side analogue. The chromatic/independence
structure of the nerve = whether the cover has a hole.

## How the geometry inspires EVERY other route (the shadows)
The topology is the unifying frame; the other proof routes are its projections:
| geometric/topological route | its shadow (other route) |
|---|---|
| Euler char of the cover nerve | the **combinatorial** inclusion-exclusion (= the cap, my S75) |
| the cover's hole / reduced `H⁰` | the **homological** route (Čech, path homology, the project's Betti work) |
| `D₇`-isotypic decomposition (Reynolds) | the **representation-theory** route (the modes = irreps, kps S31av) |
| Borsuk-Ulam odd degree | the **algebraic** route (the Gauss sum `i√7`, the negative trace, S75e/kps) |
| the torus diagonal's density | the **Diophantine** route (equidistribution/Erdős-Turán, my S74 bulk) |
| the heptagon reflection (`p mod 4`) | the **number-field** route (`Q(√−7)` wall, real vs imaginary, kps S31au) |
So: **prove the cover has a hole** — via its Euler characteristic (combinatorics), its odd degree
(Borsuk-Ulam/algebra), or the diagonal's density (Diophantine). They are one statement seen through different
lenses, and the geometry says which lens fits which piece: the **even/real/SOS** part by Euler/Brouwer
(provable), the **odd/imaginary** part by Borsuk-Ulam odd degree (the genuine obstruction).

## Visual proof route (rendered for the owner)
The companion visualization shows: the heptagon (7 sectors / de Moivre vertices), the AP danger cover at the
witness `t=1/14`, the lonely observer at `0` with the equioscillating nearest runners `±1` at exactly `1/14`,
and the real(cosine)/imaginary(sine) = Brouwer/Borsuk-Ulam decomposition. The "hole" at the observer is the
lonely point; that it persists is the Euler-char/odd-degree content.

## Honest status
- **VERIFIED:** cap = `χ_meas`(cover nerve) (the inclusion-exclusion = the Euler characteristic); the geometric
  objects (heptagon/`D₇`, torus diagonal, cover nerve).
- **SYNTHESIS:** the topology is the unifying frame — cap = cover Euler characteristic, lonely = the hole,
  witness = Borsuk-Ulam antipodal pair (kps); every other route (combinatorial/homological/representation/
  algebraic/Diophantine) is a shadow; the even/real part is Euler/Brouwer-provable, the odd/imaginary part is
  the Borsuk-Ulam odd-degree obstruction.
- **NOT a proof.** The cover-has-a-hole / the odd-degree certificate is the open content; but the geometry
  organizes all routes and pins the obstruction to the sign-isotypic / odd-degree (topology), not SOS. LRC(14) open.

Related: HYP-3242 (this), kps S31av (D₇/Borsuk-Ulam), S31au (Q(√−7) wall), HYP-3240 (covering-witness construction),
HYP-3237 (Vitali wall), HYP-3221 (the one obstruction), OPEN-Q-108.
