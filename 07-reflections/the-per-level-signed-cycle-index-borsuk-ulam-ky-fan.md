# The per-level signed cycle index, Borsuk–Ulam, Ky Fan, and the homology that refines them

*klein-2026-06-29-S2. The owner asked to merge the per-level signed cycle index with Borsuk–Ulam, Ky Fan's lemma, the ham sandwich theorem, prior Kaczynski work, and the arc-deletion paper (arXiv:2512.09332). Continues THM-584 (complement = antipodal map) and mac-mini's HYP-3543 (one R, three spectra). Companion: THM-587.*

## One generating function holds the whole metagraph spectrum

THM-584 said the iso-class metagraph is the `S_n`-quotient of the arc-hypercube `Q_d` (`d=C(n,2)`),
with eigenvalues `d-2k` and complement `R` = the antipodal map (acting `(-1)^k` on level `k`). THM-587
now gives the eigenvalue multiplicities in closed form, as the **per-level signed cycle index**

```
P_n(x) = sum_k mult(k) x^k = (1/n!) sum_{sigma in S_n} prod_{cycles c} (1 + s_c x^{ell_c}),
```

where `S_n` acts on the `d` pairs as a *signed* permutation group (a vertex swap reverses the pair it
touches — the vertex-induced subgroup of the hyperoctahedral `B_d`), and `s_c = ±1` is each pair-cycle's
orientation sign. The bit-flip is the whole story: drop the signs and you count graphs by edges
(`A000088`); keep them and you count tournament classes (`A000568`). The signs are why HYP-3540 matched
no standard OEIS row.

Two evaluations, at the antipodal fixed signs `x = ±1`, are the two faces:

| | value | meaning |
|---|---|---|
| `P_n(1)`  | `A000568(n)` | total iso classes — the all-odd-cycle tournament Burnside |
| `P_n(-1)` | `SC(n) = 2,2,8,12,88,176` | self-converse count — the **antipodal Euler / Lefschetz number** |

`P_n(-1)` is the trace of the complement involution on the class space (`dim V_+ - dim V_-`), so it is
the **equivariant Euler characteristic** of the level-graded invariant complex. And because `P_n` needs
only `n!` permutations, it computes the full metagraph spectrum at `n=7,8` (and beyond), where the
`2^{C(n,2)}` enumeration is hopeless — a closed-form spectral generator.

## Why Borsuk–Ulam is the right ceiling, not a metaphor

Complement `R: T -> T^op` flips every one of the `d` bits. On labeled tournaments it is **free**
(`x != x XOR 1` always), so `(Q_d, R)` is a free `Z_2`-space — *the* setting of Borsuk–Ulam. The
project's recurring Borsuk–Ulam invocations (the witness side THM-582/583, the LRC cap obstruction
HYP-3538) are all this one antipodal map, on three `S_n`-quotients (mac-mini's "one R, three spectra").
What THM-587 adds is the **Euler-characteristic level**: `SC(n) = P_n(-1)` is the antipodal Euler number
of the metagraph; the `R`-odd block (odd hypercube levels, `dim = (A000568-SC)/2 = #NS pairs`) is exactly
the part a Borsuk–Ulam obstruction lives in. The Perron/bulk is `R`-even (Brouwer/SOS), the obstruction is
`R`-odd (Borsuk–Ulam) — now visible as a *graded Euler characteristic*, not just a 2-dim cap block.

## Ky Fan: the signed cycle index is an alternating count

Ky Fan's lemma (the combinatorial Borsuk–Ulam, generalising Tucker) takes an antipodally symmetric
labeling with no complementary edge and forces the number of **sign-alternating** simplices to be odd —
a *signed/alternating* count made nonzero by topology. `P_n(x)` is precisely a sign-graded count over
an antipodal structure, and `P_n(-1)` is its alternating sum. The structural parallel:

- Ky Fan alternates labels `±1,...,±n` along a simplex; we alternate the hypercube **level parity**
  `(-1)^k` (the antipodal eigenvalue) along the cube.
- Ky Fan's "no complementary edge" ↔ `R` acts **freely** (no tournament is its own complement).
- Ky Fan's "odd number of alternating simplices `>= 1`" ↔ `SC(n) = P_n(-1) > 0` — self-converse
  tournaments always exist; the alternating sum is forced positive.

So the per-level signed cycle index is the tournament instance of a Ky-Fan alternating count, and its
antipodal value is the self-converse census. (We do not *need* Borsuk–Ulam to know `SC(n)>0` — rotational
constructions give self-converse tournaments directly — but the Euler-characteristic reading is what makes
`SC` a topological invariant rather than a bookkeeping total, and Ky Fan is the lens under which its
positivity is a degree, not an accident.)

## Ham sandwich: the measure form, on the cut side

Ham sandwich is Borsuk–Ulam for measures (bisect `n` masses by one hyperplane). Its tournament shadow
lives on the **cut side** of the triangle (CLAUDE.md: base-path = cut space = scores). The antipodal map
sends scores `s -> (n-1)-s`; its fixed locus is the **balanced/regular tournament** (`n` odd, all scores
`(n-1)/2`) — the simultaneous bisection of the score measure, the ham-sandwich cut of the arc set. Where
Borsuk–Ulam/Ky Fan govern the `R`-odd *witness/obstruction* (cycle side), ham sandwich governs the `R`-even
*balance* (cut side) — the two legs of the triangle under the one antipodal map. (Light touch; flagged for
the cut-space thread, not claimed as a theorem.)

## Kaczynski: refine the Euler characteristic into homology (the deliverable)

`SC(n) = P_n(-1)` is an Euler characteristic — an alternating sum of the level-graded invariant
dimensions `mult(k)`. An Euler characteristic is the alternating sum of Betti numbers, so there is a
**chain complex** waiting: `C_k =` the level-`k` `S_n`-invariant space (`dim = mult(k)`), with the cube's
signed down-operator (remove one arc) as boundary, twisted by the antipodal `Z_2`. Its `Z_2`-equivariant
homology refines `chi = SC(n)`; the `R`-odd Betti numbers are the Borsuk–Ulam obstruction classes.
**Computational homology in the Kaczynski–Mischaikow–Mrozek sense is the engine** — and it lands squarely
in the project's engineering mandate (`circulant_homology`, `tournament_tda`, THM-224's simplicial
up-Laplacian). This is HYP-3544. (The repo's prior "Kaczynski" thread, HYP-2983, is the analytic-sieve /
iterated-projection reading; both meet here: the Reynolds projection `(I+R)/2` onto the `R`-even merged
metagraph is the Kaczmarz-style projection onto the invariant subspace, the linear-algebra twin of the
2-adic descent's iterated Reynolds step in `two-order-two-structures`.)

## The arc-deletion paper sits on the same cube

arXiv:2512.09332 (El Sahili–El Zein) proves that for `n >= 8`, deleting **any** single arc of a
tournament preserves the property "contains every oriented Hamiltonian path" (with explicit exceptions).
In the cube picture a single arc move *is* one `Q_d` edge, and the Hamiltonian-path content is the
project's core invariant `H(T)` (Rédei-odd). Their result is a **robustness of `H`-content under one cube
move** — the local structure THM-584 organizes. Two hooks: (1) their threshold `n>=8` is exactly where
THM-587's cycle index now reaches the spectrum that enumeration cannot; (2) arc *deletion* (vs flip) is the
contraction/face operation on `Q_d`, suggesting a deletion–contraction reading of `H` across the metagraph
levels. Added to the bibliography backlog.

## The one-line version

The metagraph spectrum is one signed cycle index; its value at the antipodal `+1` is "how many
tournaments" (`A000568`), its value at the antipodal `-1` is "how many are their own mirror" (`SC`, the
Euler number), and the gap between them — the `R`-odd levels — is where Borsuk–Ulam, Ky Fan, and the
homology yet to be computed all live.

See [[the-one-involution-three-spectra]] (HYP-3543), [[complement-is-the-antipodal-map-of-the-arc-hypercube]]
(THM-584), [[the-pm-one-eigenspace-of-reversal-is-the-whole-split]] (HYP-3538),
[[two-order-two-structures-parity-and-descent]]. Theorem: THM-587. Hypotheses: HYP-3540 (closed),
HYP-3544 (homology, open). Reference: arXiv:2512.09332.
