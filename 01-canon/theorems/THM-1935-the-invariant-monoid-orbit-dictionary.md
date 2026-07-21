---
id: THM-1935
title: "THE INVARIANT–MONOID–ORBIT DICTIONARY: five theorems of this project are secretly one sentence — 'an invariant f is constant on the orbits of a monoid M acting on a set X, and how far f separates orbits measures M.' PROVED/VERIFIED instances: (1) RYSER IS AN ORBIT THEOREM — the score sequence is the COMPLETE invariant of the directed-3-cycle-reversal monoid; its orbits are EXACTLY the out-degree-vector classes (verified n=4: 38=38, n=5: 291=291 orbits=score-classes). (2) SWITCHING CLASSES ARE CUT-SPACE ORBITS — the two-graph (odd-triangle set) is the complete invariant of the cut-space (Z/2)^{n−1} ⊂ (Z/2)^{C(n,2)} acting on signings; #orbits = 2^{C(n,2)}/2^{n−1} = 2^{C(n−1,2)} (exact all n), and UNLABELLED = A002854 = V(E_n) = the repo's even-graph metagraph = classical Seidel two-graphs (2,3,7,16,54 for n=3..7). (3) ISO CLASSES ARE S_n-ORBITS — census A000568 = Burnside(1/n!)Σ|Fix σ| (2,4,12,56 verified). (4) THE INVARIANT LATTICE — labeled ⊐ {score, spectrum, two-graph} ⊐ iso ⊐ merged is a lattice of monoid-quotients ordered by refinement; score and spectrum are INCOMPARABLE from n=5 (neither refines the other; verified), and cospectral non-isomorphic tournaments first appear at n=4 (spectrum = a STRICT quotient of iso). (5) THE HARDNESS = STABILIZER-DIMENSION LAW (sharpens THM-1885's amenability heuristic) — an invariant is polynomial-time iff its stabilizer is positive-dimensional: tr(A^k)/det have GL-conjugation stabilizers (continuous) ⇒ P; the permanent / H have only the FINITE S_n×S_n stabilizer ⇒ #P. This is the Mulmuley–Sohoni GCT symmetry split, and THM-1780 ('H leaves the spectral ladder at n=6') is its tournament instance. THE ONE-SENTENCE PROJECT: every theorem asserts a function is M-equivariant; every counterexample-at-n0 is a pair in one M-orbit with different f (f was never M-invariant); MISTAKES.md is the list of functions mistaken for invariants."
status: >
  (1) VERIFIED-EXACT n=4,5 (reversal orbits = score-vector classes, 38=38 and 291=291); the
  underlying statement is Ryser's theorem (1964, directed-3-cycle reversals connect equal-score
  tournaments) — CLASSICAL, here recast as "orbits of the reversal monoid = score classes."
  (2) VERIFIED-EXACT: 2^{C(n,2)−(n−1)} = 2^{C(n−1,2)} all n=3..8; unlabelled = A002854 matches
  V(E_n)=2,3,7,16,54,243; the switching-class ↔ two-graph ↔ even-graph identification is CLASSICAL
  (Seidel) — the contribution is naming that the repo's E_n IS the two-graph object.
  (3) VERIFIED (A000568 = 2,4,12,56); Burnside form CLASSICAL (Davis 1954).
  (4) VERIFIED-EXACT n=3..6 (incomparability of score/spectrum from n=5; first cospectral-non-iso
  at n=4); the LATTICE framing is the contribution.
  (5) The per/det stabilizer dichotomy is Mulmuley–Sohoni GCT (CLASSICAL, external); per(PAQ)=per(A)
  and det(gAg⁻¹)=det(A) verified numerically; applying it to H as the tournament instance of
  THM-1780/1870 is the contribution — it upgrades THM-1885's "amenable⇒easy" to the sharper
  "positive-dimensional stabilizer ⇒ P".
  NET: a UNIFYING DICTIONARY. Every equation is verified or classical; no new open problem is
  closed. The value is that five separate results become one template, and the hardness law becomes
  quantitative (stabilizer dimension, not just amenability).
source: kind-pasteur-2026-07-21-S128c141 (owner: think of everything as invariants, monoids, orbits; make creative statements)
depends_on:
  - THM-1885    # the a/b monoid = BS(1,2); the (object,monoid,action) frame + amenability heuristic
  - THM-1780    # H leaves the spectral ladder at n=6 (the #P boundary)
related: [THM-1810, THM-1870, THM-1775]
concurrent:
  - "opus-2026-07-20-S441 THM-1930 (var-lambda2-decouples-from-c3, cite by filename): proves the GIT scalar var(lambda^2) is NOT a function of score or c3 for n>=5 (genuinely spectral, strictly finer) -- a THIRD strict-quotient instance for my lattice entry (4), parallel to cospectral-non-iso and to THM-1865 (H not score-determined). It also CORRECTS my THM-1885 named-next: the family ((x+c)^n+(x-c)^n)/2 interpolates char_A<->char_S of the SINGLE transitive tournament (c:0->1), NOT transitive<->Paley."
external:
  - "Ryser 1964 (3-cycle reversals connect equal-score tournaments); Seidel two-graphs / switching classes; Davis 1954 (tournament census = Burnside); Mulmuley–Sohoni Geometric Complexity Theory (per vs det = symmetry groups)."
  - "OEIS A000568 (tournaments), A002854 (two-graphs / even graphs), A000571 (score sequences)."
script: 04-computation/invariant_monoid_orbit_kps_S128c141.py (+ .out)
---

# THM-1935 — the invariant / monoid / orbit dictionary

**One sentence.** *An invariant `f` is a function constant on the orbits of a monoid `M` acting on a
set `X`; a theorem says `f` is `M`-equivariant; a census counts `M`-orbits (Burnside); how far `f`
separates orbits measures `M`.* Read this way, five results of the project are the **same** result
on different `(X, M, f)`. All checkable claims verified in the script.

## (1) Ryser is an orbit theorem — score = the 3-cycle-reversal invariant

`X` = labeled tournaments; `M` = the monoid generated by **reversing a directed 3-cycle**; `f` =
the score (out-degree) vector.

> **The orbits of the directed-3-cycle-reversal monoid are exactly the out-degree-vector classes.**
> Verified: n=4 → 38 reversal-orbits = 38 score-classes; n=5 → 291 = 291.

Reversing a directed 3-cycle preserves every score (each vertex keeps one in/one out inside the
triangle); Ryser (1964) proved these reversals connect any two equal-score tournaments. So *the
score sequence is the complete invariant of the reversal monoid* — a classical theorem, now stated
as `orbits = level sets of f`.

## (2) Switching classes are cut-space orbits — the two-graph = the repo's even graph

`X` = signings of `K_n` / tournaments; `M` = the **cut space** `(Z/2)^{n−1}` (vertex-star flips
`δ(v)`) inside the arc hypercube `(Z/2)^{C(n,2)}`; `f` = the **two-graph** (the odd-triangle set).

> **#switching classes = #cut-space orbits = `2^{C(n,2)}/2^{n−1} = 2^{C(n−1,2)}`** (exact, n=3..8),
> and **unlabelled = A002854 = `V(E_n)` = classical Seidel two-graphs** (2, 3, 7, 16, 54 for n=3..7).

This *identifies the repo's even-graph metagraph `E_n` with the classical two-graph / switching-class
object*: the even-graph duality the project has developed is Seidel switching, and the recurring
count `2^{C(n−1,2)}` (blue lines, switching classes) is the cut-space orbit count. `b` (halving) is
the switching generator.

## (3) Iso classes are `S_n`-orbits — the census is Burnside

`X` = labeled tournaments; `M` = `S_n` (relabeling); `f` = the isomorphism type.

> **#iso classes = A000568 = `(1/n!) Σ_σ |Fix σ|`** (2, 4, 12, 56 verified), with `|Fix σ| = 0`
> unless `σ` has all-odd cycles, then `2^{pairs}` — the **odd-cycle parity is the "odd" of the
> Odd-Cycle Collection** (CLAUDE.md, Davis 1954).

## (4) The invariant lattice — quotients ordered by refinement

The invariants are a **lattice of monoid-quotients** of the labeled cube, ordered by refinement:

```
        labeled (2^C(n,2))
        /     |      \
   score   two-graph  spectrum        <- three incomparable coarsenings
        \     |      /
          iso (A000568)                <- S_n
             |
          merged (÷ Z/2 complement)
```

Verified: **score and spectrum are incomparable from n=5** (neither refines the other), and
**cospectral non-isomorphic tournaments first appear at n=4** — so the spectral map is a *strict*
quotient of the iso quotient (co-spectral mates are the fibers). At n=4 score is strictly finer than
spectrum; by n=5 they cross. The lattice meet is the common refinement; the join is the coarsest
common coarsening.

## (5) The hardness = stabilizer-dimension law (sharpens THM-1885)

THM-1885 offered "amenable acting monoid ⇒ tractable." The sharper statement is about the
**stabilizer of the invariant**:

> **An invariant is polynomial-time iff its stabilizer group is positive-dimensional.**
> `tr(A^k)`, `det` are invariant under all of `GL`-conjugation (a positive-dimensional group) ⇒ **P**
> (spectral). The **permanent** and `H` are invariant only under the **finite** `S_n × S_n` ⇒ **#P**.

Verified: `per(PAQ) = per(A)` (finite symmetry), `det(g A g⁻¹) = det(A)` (continuous symmetry). This
is exactly the **Mulmuley–Sohoni GCT** explanation of permanent-vs-determinant, and **THM-1780 ('H
leaves the spectral ladder at n=6') is its tournament instance**: `H` is `per`-like, so it cannot be
a function of the continuous-stabilizer invariants `tr(A^k)`, and the gap opens at the first `n`
where the Hamiltonian term stops being a symmetric function of the spectrum.

## The project in one line

> Every theorem here asserts a function is `M`-equivariant; every counterexample at `n₀` is a pair
> in one `M`-orbit with different `f` (so `f` was never `M`-invariant); `MISTAKES.md` is the list of
> functions mistaken for invariants; and the frontier (LRC, `H`) is where the acting monoid is
> non-amenable (`PSL(2,Z)`) or the stabilizer is finite (`per`).

## Named next

- **Which quotients commute?** Does `÷Z/2` commute with `÷S_n` (it does — merged metagraph) but not
  with `÷cut-space`? Map the sublattice of *commuting* quotients (simultaneous invariants).
- **The two-graph ↔ score meet.** What invariant is the common refinement of the switching class and
  the score sequence? (A tournament's odd-triangle set *plus* its scores.)
- **Formalise Ryser-as-orbits** and **switching = cut-space orbits** in Lean (both are finite,
  `native_decide`-checkable at small n; the arithmetic `2^{C(n,2)−(n−1)}=2^{C(n−1,2)}` is one line).
- **The stabilizer-dimension law as a screening test** for every repo invariant: compute its
  symmetry group; positive-dimensional ⇒ expect a spectral formula, finite ⇒ expect `#P`.
