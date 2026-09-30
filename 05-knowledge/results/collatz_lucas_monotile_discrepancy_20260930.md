# Unit-minus-one clocks: the owner's Lucas identities are Cayley–Hamilton for the Fibonacci automorphism, the monotile paper's torsion orders and the Collatz cycle clocks are the same kind of number (fixed-point counts `|N(u^j - 1)|`), a torsion-point census of the Collatz shapes to `A = 22`, the prime-occurrence lattice of the clocks, and the discrepancy paper's "rigid object plus sparse random perturbation" as a mirror of the near-critical band

**Session:** opus, `collatz-posets-zeta5-20260927` (S15, eleventh note), 2026-09-30.
**Owner's directive:** "consider deeply how the following equation and two
attached papers relate to our continued investigation and merge in new ideas
you explore and discover on the way towards a collatz proof: `φ^8 + 1 = 7φ^4`,
`φ^10 = 1 + 11φ^5`."
**The two papers.** (1) Victor Reis and Zhao Song, *Three Standard Deviations
Suffice While One Does Not*, arXiv:2609.34471v1 (28 Sep 2026): every
`A ∈ [-1,1]^{n×n}` has a signing `x ∈ {±1}^n` with
`||Ax||_∞ ≤ sqrt(3 arsinh(10)) sqrt(n) + 4 < 2.9992 sqrt(n) + 4` (by minimising
a potential built from `φ(t) = 5/sinh(t^2/3)` over the cube and rounding the
minimiser to the nearest sign vector), and for every power of two `n ≥ 2^50`
there is a `±1` matrix with `||Ax||_∞ > (1 + 2^-22) sqrt(n)` for *every*
signing: a Hadamard matrix with a `2^-22` fraction of its columns replaced by
independent random signs (proving Bandeira's Conjecture 12,
`limsup D_n^± > 1`, the first improvement on the Olson–Spencer bound
`disc(H) ≥ sqrt(n)` that Parseval gives for every Hadamard `H`). (2) The
attached `fibonacci-monotile.pdf`, *A chiral aperiodic polygon with Fibonacci
monodromy* (author line "Pingyou Ltd"; no identifier in the file): a 27-gon
that is a strictly chiral aperiodic monotile, with supertile expansion
`λ = 4 + sqrt(15)` and boundary growth `θ = φ^3 = 2 + sqrt 5`; one level of its
boundary subdivision induces, on the boundary's homology, the Fibonacci
automorphism `Q = [[0,1],[1,1]]` (class `-Q^3`); the mapping tori `T_n` are
cyclic covers of the Gieseking manifold (of the figure-eight knot complement
at even `n`) with torsion homology `Z[φ]/((-φ)^{3n} - 1)`, of order the Lucas
number `L_{3n}` (`n` odd) or `L_{3n} - 2` (`n` even); the rings
`O_j = Z[φ]/(φ^j - 1)` have order `|L_j - 1 - (-1)^j|`, are fields exactly for
`j ∈ {3, 4}` or `j = p ≥ 5` prime with `L_p` prime, and the boundary's
Hausdorff dimension is the torsion growth against the expansion,
`δ = 6 log φ / log(4 + sqrt 15) = 3 log M(t^2 - 3t + 1)/log λ` (Mahler measure
of the figure-eight knot's Alexander polynomial).

**Status: PROVED (elementary) for Propositions 1–3; FINITE-EXACT for the
census, the coverage table, the discrepancies and the spectra; CITED for the
two papers, Gersonides, Zsigmondy, Böhm–Sontacchi, Olson–Spencer, Schmidt and
the repo results named by path; ANALOGY and NUMEROLOGY typed where used;
Collatz OPEN. Independent audit OWED (weekly subagent limit; self-checks
recorded: two counters agree on every shape with `A ≤ 12`, Cayley–Hamilton
verified in integers, spectra verified numerically).**
Scripts `04-computation/experiments/collatz_lucas_monotile_discrepancy_20260930.py`
(parts 1–4) and `collatz_lucas_monotile_20260930_census.py` (parts 5–9), outputs
alongside (`.out`).

**Reconciliation with the frontier that moved while this session wrote.**
Since the tenth note, sessions S16–S23 (`collatz-poset-dag-20260927`) entered
STICKY into the barrier atlas as control s7
([`collatz_sticky_20260927_size_coupled_persistence.md`](collatz_sticky_20260927_size_coupled_persistence.md)),
digested Mazur's positive-density theorem
([`mazur_positive_density_20260928.md`](mazur_positive_density_20260928.md)),
built the Fourier profile `M(h)` of the 3-adic Syracuse law and its renewal
structure ([`collatz_five_mirrors_20260929.md`](collatz_five_mirrors_20260929.md)),
and the character spectrum of the law
([`collatz_three_mirrors_20260929.md`](collatz_three_mirrors_20260929.md));
the mac-mini session proved THM-4515–4517 on cycle necklaces, perfect-power
clocks and root-uniform basin bounds
([`collatz_necklace_20260929_fair_splits_power_clocks_basins.md`](collatz_necklace_20260929_fair_splits_power_clocks_basins.md)).
A parallel session's note of the same day,
[`collatz_circulant_20260930_circulants_lucas_cubic_monotile.md`](collatz_circulant_20260930_circulants_lucas_cubic_monotile.md)
(THM-4520; the Pillai gaps `2^K - 3^X` that are Lucas or Fibonacci numbers or
Fibonacci-group orders, all with `K ≤ 8`), answers the same directive from
the circulant side and was found after this note was pushed; the two are
complementary (that note asks which clocks *are* Lucas numbers, this one what
kind of number all three families are). Nothing below duplicates them: this note is about *which numbers the clocks
are* (fixed-point counts), which primes divide them (a lattice), and the
exact torsion-point census; the necklace note is about the necklace
combinatorics and the perfect-power clocks, the mirrors notes about the
3-adic Fourier side.

## 1. The owner's two identities, read exactly (EXACT; the numerology typed)

`φ^{2n} + (-1)^n = L_n φ^n` for every `n`, because `φ^n` and `ψ^n = (-1/φ)^n`
are the roots of `x^2 - L_n x + (-1)^n` (trace `L_n`, norm `(-1)^n`). The
owner's equations are `n = 4` (`L_4 = 7`) and `n = 5` (`L_5 = 11`). In the
basis `(1, φ)` multiplication by `φ` is `Q = [[0,1],[1,1]]`, `Q^n =
[[F_{n-1}, F_n],[F_n, F_{n+1}]]`, and the identity is the Cayley–Hamilton
theorem for `Q^n`: `Q^8 - 7Q^4 + I = 0` with `Q^4 = [[2,3],[3,5]]`, and
`Q^10 - 11Q^5 - I = 0` with `Q^5 = [[3,5],[5,8]]` (script part 9, in
integers). The monotile paper's torsion rings are `O_j = coker(Q^j - I)`, of
order `|det(Q^j - I)| = |1 - L_j + (-1)^j|`: at the owner's two indices the
traces are `7` and `11` and the torsion orders are `|O_4| = 5` and
`|O_5| = 11`; at `j = 10`, `|O_10| = L_10 - 2 = 121 = 11^2 = 5! + 1`, and
`O_10 ≅ (Z/11)^2` is not a field (Proposition 7.3 and Corollary 7.8 of the
paper: `L_5^2 - 5F_5^2 = -4` makes `5` a square modulo `11`). So the two
torsion orders `5, 11` are the two smallest Brocard roots
(`4! + 1 = 5^2`, `5! + 1 = 11^2`) of the [eighth note](collatz_doubling_tower_brocard_20260927.md),
and `71` is not a torsion order (no Lucas number is `71` or `73`). Typed
NUMEROLOGY: exact statements, no mechanism.

The one exact thing the pair `(7, 11)` does in the Collatz thread is the
`-17` cycle, `(A, p) = (11, 7)` halvings and odd steps, `3^7 - 2^11 = 139`.
Section 5 says why the Lucas numbers `L_4, L_5` are there.

## 2. Three clocks, one formula (PROVED, elementary)

Three families of integers in the thread are values of the same expression,
the absolute norm of `u^j - 1` for a unit `u` of a ring of `S`-integers, and
each counts fixed points of a compact-group endomorphism:

| family | ring, unit | `|N(u^j - 1)|` | as a resultant | counts fixed points of |
|---|---|---|---|---|
| Mersenne `2^n - 1` (the tower orders of THM-447 / HYP-9162 are `2^k - 1`) | `Z[1/2]`, `u = 2` | `2^n - 1` | `|Res(t^n - 1, t - 2)|` | `x ↦ 2^n x` on `R/Z` |
| associated Mersenne `|L_j - 1 - (-1)^j|` (OEIS A001350: `0, 1, 1, 4, 5, 11, 16, 29, 45, 76, 121, 199, 320, ...`) | `Z[φ]`, `u = φ` | `|O_j|` | `|Res(t^j - 1, t^2 - t - 1)|` | `Q^j` on `R^2/Z^2` (Proposition 7.2(b) of the paper) |
| Collatz clock `|2^A - 3^p|` | `Z[1/6]`, `u = 3^p/2^A` | `|2^A - 3^p|` | on the ray `A = vp`: `|Res(t^p - 1, 2^v t - 3)|` | `x ↦ (3^p/2^A) x` on the 2-adic solenoid |

**Proposition 1 (the torsion-point reading of the cycle equation).** Let
`Σ_2` be the 2-adic solenoid, the Pontryagin dual of the discrete group
`Z[1/2]`. For a shape `(A, p)` the endomorphism `û` of `Σ_2` dual to
multiplication by `u = 3^p/2^A` on `Z[1/2]` has exactly `|2^A - 3^p|` fixed
points, forming the group `Fix(û) ≅ Z/(2^A - 3^p)`. For a word
`w = (v_1, ..., v_p)` of that shape (`v_i ≥ 1`, `Σ v_i = A`) with carry
`S_w = Σ_{i<p} 3^{p-1-i} 2^{d_i}` (`d_i = v_1 + ... + v_i`, `d_0 = 0`), the
affine cycle map `x ↦ (3^p x + S_w)/2^A` has the unique fixed point
`x_w = S_w/(2^A - 3^p)`, and its class `S_w mod (2^A - 3^p)` is a point of
`Fix(û)`. An integer Collatz cycle of shape `(A, p)` exists iff some word of
the shape has `S_w = 0` in `Fix(û)`; then every rotation of `w` does, and the
`p` rotations are the `p` points of the cycle.

*Proof.* `Fix(û)` is the dual of `coker(u - 1)` on `Z[1/2]`;
`(u - 1)Z[1/2] = (3^p - 2^A)Z[1/2]` since `2^A` is a unit there, and
`Z[1/2]/(3^p - 2^A) ≅ Z/(3^p - 2^A)` because `3^p - 2^A` is odd, so `2` is
invertible modulo it. The fixed-point equation is `2^A x = 3^p x + S_w`; the
integrality condition `S_w ≡ 0 mod (2^A - 3^p)` is the classical cycle
condition (Böhm–Sontacchi; in the repo, the cycle gate of
[`collatz_mod6_20260921_pillai_convergents_cycle_gates.md`](collatz_mod6_20260921_pillai_convergents_cycle_gates.md)
and THM-4484). Rotation invariance: the rotated word's fixed point is the
image of `x_w` under the corresponding partial composition, an integer iff
`x_w` is. ∎ (The content is the reading, not the condition: the cycle points
of a shape are the torsion points of the linear part's fixed-point group
that the words of the shape happen to hit, exactly as `O_j` is the group of
fixed points of `Q^j` on the torus. The resultant column is checked in
script part 9: `Res(t^7 - 1, t - 2) = -127`, `Res(t^5 - 1, t^2 - t - 1) =
-11`, `Res(t^3 - 1, 4t - 3) = 37 = 2^6 - 3^3`.)

**The dictionary this fixes.** The paper's zero rings `O_1 = O_2 = 0`
(`φ - 1 = φ^{-1}` and `φ^2 - 1 = φ` are units) are the Collatz shapes with
unit clock, Gersonides' `2^A - 3^p = ±1` at `(1,1), (2,1), (3,2)`: the *free*
cycles `-1, +1, -5` of THM-4484, where the fixed-point group is trivial and
every word of the shape closes a cycle (script part 5: all `1, 1, 2` words
hit). Both are instances of the finiteness of the `S`-unit equation
`u^j - 1 = unit`. The paper's field cases (`O_j` a field iff `j ∈ {3,4}` or
`j = p` with `L_p` prime) correspond to prime clocks: the sporadic `-17`
cycle lives in the field case, `Z/139`, and both `2` and `3` are primitive
roots modulo `139` (orders `138`, script part 7). Typed ANALOGY: the same
formula in three rings, no theorem transferred.

**The doubling chain.** Along `j = 2^k` the torsion is `|O_{2^k}| = L_{2^k}
- 2 = 5 F_{2^{k-1}}^2` (`5, 45, 2205, 4870845`; script part 9), whose traces
`L_{2^k} = 3, 7, 47, 2207` are the historian index's "self-observation chain"
(`L_{2n} = L_n^2 - 2` for even `n`). The Mersenne tower of HYP-9162 and this
Lucas tower share `3, 7` and then part (`15, 31` against `47, 2207`): the
two towers are the norms of `2^{2^k} - 1` and of `φ^{2^k} - 1`.

## 3. The torsion-point census (FINITE-EXACT)

Script part 5 counts, for every shape with `A ≤ 22` by direct enumeration of
the `C(A-1, p-1)` words (cut sets of `{1, ..., A-1}`) and for the five
near-critical shapes `(23,14), (23,15), (24,15), (25,16), (27,17)` by a
residue DP (each odd step permutes `Z/D` by `r ↦ 3r + 2^d`; the two counters
agree on every shape with `A ≤ 12`), the words with `S_w = 0` in
`Z/(2^A - 3^p)`; 258 shapes, 119 s.

| primitive hit | clock | zero words / all words | no-descent words |
|---|---|---|---|
| `(1,1)` | `-1` | `1 / 1` | `1` |
| `(2,1)` | `+1` | `1 / 1` | `1` |
| `(3,2)` | `-1` | `2 / 2` | `1` |
| `(11,7)` | `-139` | `7 / 210` | `30` |

Every other hit is a repeat `(ka, kb)` of one of these four (the `k`-fold
traversal of the same cycle: `(22,14)` has `7` zero words of `203490`, the
rotations of the doubled `-17` word), and no other primitive shape hits in
range. The shared convergent `8/5` of `φ` and `log_2 3` carries no cycle
(`35` words, none zero modulo `13`), nor does the golden shape `(13,8)`
(`792` words modulo `1631 = 7 · 233`; that `233 = F_13` is NUMEROLOGY). This
census is implied by the cycle-length theorems (Eliahou 1993; Hercher 2023;
Barina's verification) and is recorded as a consistency check of the
reading, not as evidence.

**Coverage.** For small clocks the words do not even cover the torsion group:
`(5,3)` has `6` words hitting `4` of the `5` residues, `(6,4)` `10` words
hitting `9` of `17`, `(7,4)` `20` words hitting `18` of `47` (part 5, table).
Collisions `S_w ≡ S_{w'}` are the rule, and `0` is missed; Direction D30.

**The uniform-residue heuristic and the dimension exponent.** If the
`C(A-1,p-1)` carries were uniform in `Z/D`, a shape would carry
`C(A-1,p-1)/(p|D|)` cycles; summed over the examined non-Gersonides shapes
this is `4.46` (`3.90` with the no-descent refinement, since the word read
from a cycle's least element never descends below it) against the one actual
sporadic cycle (whose own term is `210/973 = 0.216`). The heuristic
overcounts because carries are not uniform for small shapes; its tail is the
real content. With `W_nd(A,p) ≈ 2^{hA}·poly` no-descent words
(`h = h(log_3 2) = 0.9499555272`, the Hausdorff dimension of the never-descending
2-adic set, Theorem 1a of
[`collatz_procgen_20260922_sibling_dimension_ladder.md`](collatz_procgen_20260922_sibling_dimension_ladder.md);
the value is recomputed in part 6 to ten digits) against a clock
`|2^A - 3^p| = 2^A |1 - 3^p/2^A|`, the expected number of cycles of shape
`(A, p)` is `2^{-(1-h)A}/|1 - 3^p/2^A|` up to polynomial factors: the
per-halving decay exponent is the *codimension* `1 - h = 0.0500445` of the
no-descent set, and the denominators concentrate the mass on the convergents
of `log_2 3` (tail over the three nearest-critical shapes for `23 ≤ A ≤ 60`:
`0.358`, led by `(27,17)`, `(38,24)`, `(46,29)`, `(30,19)`). This is
Lagarias's classical cycle heuristic with the exponent identified as the
dimension deficit; the identification is the new sentence, the heuristic is
not.

**A Bowen-type reading (ANALOGY, both sides exact).** Theorem E(v) of the
monotile paper is "dimension = torsion growth / expansion":
`δ = lim log L_{3n} / (n log sqrt λ) = 1.399`. The thread's exact analogue is
Theorem 1a above: the never-descending set has `2^{hN}` admissible words of
`T`-length `N` against the expansion `2` per `T`-step, so `dim_H = h = 0.95`;
read archimedeanly, the residues below `2^N` that have not descended after
`N` steps number `(2^N)^h`, a set of "dimension `0.95`" against the
Krasikov–Lagarias basin bound `X^0.84` and the conjectured `X^1`. On the
constant-valuation rays the Collatz "Alexander polynomial" is `2^v t - 3`,
of Mahler measure `max(2^v, 3)`, and the clock `|2^{vp} - 3^p|` grows at
that rate; the critical slope `v = log_2 3` is where the two terms of the
Mahler measure tie, and it is not an integer: Gersonides again, no integer
ray is critical.

## 4. Which primes divide the clocks (PROVED, elementary), and Zsigmondy

**Proposition 2 (prime occurrence).** For a prime `ℓ ≥ 5`, the set of shapes
`(A, p) ∈ Z^2` with `ℓ | 2^A - 3^p` is the kernel of the homomorphism
`Z^2 → F_ℓ^×`, `(A, p) ↦ 2^A 3^{-p}`, a lattice of index `|<2, 3>|`, the order
of the subgroup generated by `2` and `3`. In particular every prime `ℓ ≥ 5`
divides some clock, and `ℓ` divides `2^A - 3^p` for a positive proportion
`1/|<2,3>|` of all shapes. *Proof.* Immediate. ∎ This is the Collatz form of
Proposition 7.4 of the monotile paper (the levels with `p`-torsion are the
multiples of the order of `-φ^3` modulo the primes above `p`): there the
lattice is a subgroup of `Z`, here of `Z^2`. Script part 7: among the primes
`5 ≤ ℓ < 2000`, `2` and `3` together generate `F_ℓ^×` for `222/301 = 0.7375`
(the two-generator Artin density; the GRH-conditional constant is CITED from
memory as Matthews 1976, unverified here); the first primes of index `2` are
`23` and `47`, where `2` and `3` are both squares and generate exactly the
squares. For `ℓ = 139` the lattice has index `138` (both `2` and `3` are
primitive roots) and is generated by `(11, 7)` and `(7, 17)`
(`11·17 - 7·7 = 138`); its points with `A, p < 30` are
`(0,0), (3,27), (7,17), (11,7), (18,24), (22,14), (26,4)`.

**Proposition 3 (repeats).** The `k`-fold repeat of a word has
`S_{w^k} = S_w (2^{kA} - 3^{kp})/(2^A - 3^p)`, so the same fixed point; the
cofactor `((2^A)^k - (3^p)^k)/(2^A - 3^p)` is the `k`-th "cyclotomic" factor
of the ray, and by Zsigmondy's theorem (CITED; here checked by factoring for
`k ≤ 12` on the rays of `(1,1), (2,1), (3,2)` and `k ≤ 5` on the ray of
`(11,7)`, script part 9) it carries a primitive prime divisor for every
`k ≥ 2`: the torsion group of a repeated shape always has new primes, exactly
as the monotile tower acquires new `p`-torsion at the ranks of apparition.
*Proof of the identity.* The cycle map of `w^k` is the `k`-th iterate of
`f(x) = (3^p x + S_w)/2^A`, namely
`f^k(x) = (3^{kp} x + S_w Σ_{i<k} 3^{pi} 2^{A(k-1-i)})/2^{kA}`, and the
geometric sum is the stated cofactor (checked on the `-1` cycle: carries
`1, 5, 19` at `k = 1, 2, 3`). ∎ (NUMEROLOGY, one line: on the ray of
`(2,1)`, `2^10 - 3^5 = 781 = 11 · 71`, and `3^7 - 2^7 = 2059 = 29 · 71`: the
Brocard roots `11, 71` as primitive divisors.)

## 5. Why the Lucas numbers `7, 11` are the `-17` cycle's shape (EXACT; the numerology dissolved)

`log_2 3 = [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, ...]` and `φ = [1; 1, 1, 1, ...]`
share the prefix `[1; 1, 1]`, so they share the convergents `1/1, 2/1, 3/2`,
and `8/5 = F_6/F_5` is a convergent of both (though `log_2 3` skips `5/3`).
The three free cycles sit on the shared convergents, `(A, p) = (1,1) =
F_2/F_1`, `(2,1) = F_3/F_2`, `(3,2) = F_4/F_3`; the sporadic cycle's `11/7` is
the semiconvergent of `log_2 3` between `3/2` and `19/12` — the mediant of
`3/2` and `8/5`, as the [coincidence atlas](collatz_coincidence_atlas_20260927.md)
§4 recorded — and the mediant of two consecutive Fibonacci ratios is a Lucas
ratio: `(F_{n+1} + F_{n+3})/(F_n + F_{n+2}) = L_{n+2}/L_{n+1}`, here
`(3 + 8)/(2 + 5) = 11/7 = L_5/L_4`. That is the whole mechanism behind the
owner's `7` and `11`: the golden ratio enters the Collatz cycle data only
through the two shared partial quotients, and the Lucas numbers only as the
mediant one step past the shared segment. The `-17` word
`(1, 1, 1, 2, 1, 1, 4)` has nothing golden in it, and `11/7` is a best
approximation of `log_2 3` of the first kind (`|log_2 3 - 11/7| = 0.0135 <
0.0150 = |log_2 3 - 8/5|`) but not of the second (`|7 log_2 3 - 11| = 0.095 >
0.075 = |5 log_2 3 - 8|`).

**The zeta function of `T` on `Z` (EXACT reformulation of the cycle half).**
The brute-force census of the shortcut map on `|x| ≤ 2·10^5` (script part 8)
finds the cycles through `0, -1, 1, -5, -17` with `T`-periods `1, 1, 2, 3,
11` (the `T`-period is the halving count `A`). The cycle half of the Collatz
conjecture over `Z` is the statement that the Artin–Mazur zeta function of
`T` on `Z` is exactly
`ζ_T(z) = 1/((1 - z)^2 (1 - z^2)(1 - z^3)(1 - z^11))`.
The periods `1, 1, 2, 3` are `F_1, ..., F_4` and `11 = L_5`, for the reason
just given. (The repo's earlier dynamical zeta,
[`collatz_procgen_20260923_trunk_rh.md`](collatz_procgen_20260923_trunk_rh.md)
§4.3, is of the escape map `Ψ`, a different object.)

## 6. The discrepancy paper as a mirror (ANALOGY; two FINITE-EXACT facts)

**What it proves.** `disc(A) = min_x ||Ax||_∞` over signings; for a Hadamard
matrix Parseval gives `||Hx||_∞ ≥ ||Hx||_2/sqrt n = sqrt n` for every `x` at
once (Olson–Spencer 1978); Schmidt's theorem gives `disc(H_{2^m})/sqrt(2^m) →
1` for the Sylvester matrices, with equality `disc = sqrt n` at even `m`
(from `H_4 y = 2y`, `y = (1,1,1,-1)`); Bandeira conjectured that Hadamard
matrices are not extremal, `limsup D_n^± > 1`, and Reis–Song prove it by
replacing a `2^-22` fraction of the columns by random signs; the upper bound
`sqrt(3 arsinh(10)) sqrt n + 4` comes from a potential
`Φ(z) = Σ_i φ((Δ - <a_i, z>)/sqrt(s(z)))`, `φ(t) = 5/sinh(t^2/3)`, minimised on
the cube with a slack `R(z) = n - ||z||^2 ≥ 4`, whose minimiser rounds to a
sign vector losing only `4` (Talagrand's inequality controls the `ℓ_1` mass
of `Hu` in the lower bound).

**Fact 1 (the repo's own Hadamard family).** The skew-Sylvester tower of
THM-447 / HYP-9162 produces skew-Hadamard matrices `S = M + I` with
`S S^T = nI`, so `disc(S) ≥ sqrt n` by the same line. Exactly (part 3):
`disc(S) = 2, 2, 4, 6` at `n = 2, 4, 8, 16` against Sylvester's `2, 2, 4, 4`;
the skew tower does *not* attain `sqrt n` at `16`. Reason (part 9): `M^2 =
-(n-1)I` gives `S` the spectrum `1 ± i sqrt(n-1)`, non-real, so there is no
real eigenvector `Sy = sqrt(n) y` to round, whereas Sylvester's spectrum is
`±sqrt n` with the attaining eigenvector. PROVED for the spectrum,
FINITE-EXACT for the discrepancies (local search gives `8` for `H_32`; the
paper cites `disc(H_512) ≤ 28 < 32`).

**Fact 2 (the structural mirror, typed ANALOGY).** The paper's shape is: an
exact `L^2` identity valid for all `2^n` signings at once forces a universal
`L^∞` bound, and the extremal object is not the most symmetric one but a
sparse random perturbation of it. The thread has both halves in its own
coordinates. The universal identity is Terras's bijection: every word occurs
exactly once modulo `2^N`, so no residue is forced to descend and no
low-discrepancy word is forbidden — the Collatz counterpart of "every matrix
has a signing with `||Ax||_∞ < 3 sqrt n`" is that the signing is *not free*:
the residue class fixes the word, and a bounded-discrepancy (Sturmian) word
to depth `N` is paid in size `2^N` (the universal size price of the
[sixth note](collatz_generic_price_topology_20260927.md); the
[guards note](collatz_guards_20260921_discrepancy.md) proves no integer orbit
keeps a bounded strip). The perturbation half mirrors S21's finding
([`collatz_five_mirrors_20260929.md`](collatz_five_mirrors_20260929.md)
§2c) that the maximal no-descent rate `0.947–0.950` sits slightly *above* the
pure critical value `3^{h*-1} = 0.9465`: the extremiser is a near-critical
band of width `≈ h/8`, not the pure Sturmian word — "Hadamard plus `2^-22`
random columns" against "critical word plus a thin band". No theorem
transfers; the mirror says where to look (D33).

## 7. Verdicts

| claim | status |
|---|---|
| `φ^{2n} + (-1)^n = L_n φ^n` is Cayley–Hamilton for `Q^n`; the owner's `7, 11` are the traces at `j = 4, 5`, the torsion orders there are `5, 11`, and `|O_10| = 121 = 5! + 1` | EXACT; the Brocard link NUMEROLOGY |
| Mersenne, associated Mersenne (Lucas) and Collatz clocks are `|N(u^j - 1)|` and count fixed points (Proposition 1; resultants checked) | PROVED, elementary |
| zero rings ↔ Gersonides' free cycles; field cases ↔ prime clocks (`139`) | ANALOGY (same formula, no transfer) |
| torsion-point census: hits only at `(1,1), (2,1), (3,2), (11,7)` and repeats, all shapes `A ≤ 22` plus five to `A = 27`; coverage deficits | FINITE-EXACT (implied by the cycle-length theorems; consistency check) |
| heuristic decay exponent per halving = `1 - h(log_3 2) = 0.0500445`, the codimension of the no-descent set; mass on the convergents | EXACT identification, heuristic itself classical |
| Proposition 2 (kernel lattice of index `|<2,3>|`; `139` generated by `(11,7), (7,17)`), Proposition 3 (repeat cofactor, Zsigmondy) | PROVED / CITED / FINITE-EXACT |
| `11/7 = L_5/L_4` = mediant of `3/2` and `8/5` = the first semiconvergent after the shared golden prefix; zeta `1/((1-z)^2(1-z^2)(1-z^3)(1-z^11))` | EXACT; the golden numerology dissolved |
| skew tower `disc = 2, 2, 4, 6` vs Sylvester `2, 2, 4, 4`; spectra `1 ± i sqrt(n-1)` vs `±sqrt n` | FINITE-EXACT / PROVED |
| discrepancy paper ↔ Terras universality, size price, near-critical band | ANALOGY |
| Collatz | OPEN |

## 8. Directions

* **D30 (coverage of the torsion group).** For which shapes do the
  `C(A-1,p-1)` carries fail to cover `Z/(2^A - 3^p)`, and is `0` missed for a
  structural reason? The zero fibre is a union of rotation classes (the
  rotation action on carries, `c_{rot w} = (3^{w_0} c_w + w_0 D)/2`, is the
  necklace note's sidecar), so the question is about the orbit structure of
  that action on `Z/D`; by Fourier inversion the zero count is
  `(1/|D|) Σ_{k mod D} Σ_w e(k S_w/D)`, and "no cycle of this shape" says the
  non-trivial characters cancel the main term `C(A-1,p-1)/|D|` exactly — a
  character sum over words modulo the clock, coprime to `6`, beside the
  3-adic character spectrum of S23.
* **D31 (index of `<2,3>`).** The primes of index `> 1` (`23, 47, ...`) are
  those where `2` and `3` lie in a common proper subgroup; the density of
  index-one primes (`0.7375` to `2000`) is a two-generator Artin constant
  (CITED-UNVERIFIED). Clocks divisible only by index-one primes versus the
  rest: does the cycle gate see the index?
* **D32 (Lefschetz for cycles).** Proposition 1 makes the cycle count of a
  shape a count of torsion points hit by a word map into a fixed-point group;
  is there a trace formula for `#{w : S_w = 0}` (a Lefschetz number of the
  word action on `Z/D`), so that the cycle half of Collatz becomes a
  vanishing statement for a trace?
* **D33 (the band as the extremiser).** Reis–Song's lower bound is a
  perturbation computation around an exactly solvable object; S21's
  near-critical band is the same phenomenon in the no-descent rate. Compute
  the rate as a function of the band width (the Collatz `δ`) and compare with
  their `(1 + δ) sqrt n` gain: is the excess over `3^{h*-1}` first order in
  the band width?
* **D34 (skew vs Sylvester at higher order).** Exact `disc` of the skew tower
  at `32, 64` (branch and bound), and whether Schmidt's limit `1` holds for the
  skew family too; if not, the skew tower is a natural candidate for
  `limsup D_n^± > 1` without randomness.
