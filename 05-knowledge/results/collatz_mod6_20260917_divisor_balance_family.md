# Divisor balance family: density, linear cells, k-free refinement, Omega blindness

**Status:** PROVED (scoped) classifications of the linear family `F = alpha S + beta U` (finitely many exponent profiles in every cell except `(1,0)`, `(2,0)`, which are described exactly) and of the `k`-free refinement `F = S_k + U` for all `k`; PROVED lattice reading of the inherited defect minimum; PROVED conditional leading asymptotic for `#{p^2qr <= x}` on three CITED classical inputs (Landau; Chebyshev plus Mertens; PNT), with no novelty claim for that leading term; FINITE-EXACT densities, crossover bookkeeping and sandwich census to `10^7`; one HEURISTIC scale remark; one REFUTED typed analogy and one REFUTED clause of the wave-one draft (recorded in 2.4); SCOPE statements for the named canon theorems (no map found). OPEN: the secondary term of `#{p^2qr <= x}`, an effective certificate that no sign change or tie occurs beyond `147924`, and explicit profile lists for `alpha, beta >= 4`. No claim here touches twin primes, Collatz, or Goldbach. Research result note, not a canon theorem ID.

**Session:** collatz-mod6-20260917 (mac-mini), lane `divisor_balance_family`; wave-one explorer, two independent audits (recompute, proof-audit), finalized 2026-09-21 from the recovered script.

**Convention.** `N >= 2` throughout. `N = 1` has `F = S = U = 0` and would trivially solve every cell of section 2 and every `k` of section 3; it is excluded once here, as in the inherited note.

## Inheritance and concept board

The closest proved mechanism is the exponent-box mechanism of section DB1 of [arithmetic_braids_20260917_divisors.md](arithmetic_braids_20260917_divisors.md) (itself resting on the multiplicative divisor fibre of THM-2422 cited there): divisors of `N = prod p_i^{a_i}` form the box `prod {0..a_i}`, squarefree divisors its Boolean subcube `{0,1}^r`, and `F - S` counts the box points outside the cube after removing `N`; equality `F = S + U` says this excess equals `r`. The inherited equations are, with the inherited tags kept verbatim,

```text
F = prod(a_i+1) - 2,   S = 2^r - 1 - [N squarefree],   U = r - [N prime]        (DB1)
prod(a_i+1) = 2^r + r + 1        (nonsquarefree case of F = S + U)               (DB2)
D = F - S - U >= 2^(r-1) - r - 1, equality exactly at profile (2,1,...,1)          (DB3)
```

with the section-level facts: `F = S + U iff N in {p, p^3, p^2qr}` (inherited section DB1), the defect table with `p^3q` at `D = 1` and `p^2q^2` at `D = 2` (inherited section DB2), and the almost-prime remark that an `Omega` class forgets which prime is repeated and that neither `Omega` nor `omega` alone determines `D` (inherited section DB3). Those are cited, not re-derived. The canonical hostile is `p^2` (and `p^2q`): `D = -1`, explained in 3.1 as the apex of the divisor lattice colliding with an atom image. The corrected near miss is the wave-one draft's clause "`(1,1)` is the unique prime-cube cell whose solution profiles have pairwise distinct `Omega`", REFUTED by the infinite cell `(2,0)` (2.4). The least-used sidecar is the set of *unreached* non-squarefree divisors `p^2 d`, `2 <= omega(d) <= r-2`, whose count is exactly the (DB3) minimum (3.1); a second, cheaper sidecar is the distinction between a *tie* and a *sign change* in the density race (1.3).

| Live concept | Retained coordinate | Predicate preserved | Hostile or decisive test |
|---|---|---|---|
| Linear cell `F = alpha S + beta U` | support `r`, multiplicative partitions of `T(r)` | finiteness of the profile list | `(1,0)`, `(2,0)` infinite; `(1,3)` at `r = 5` is tight |
| Divisor lattice `D(N)` | order filter `Q` of non-squarefree divisors | `F - S - U = |Q| - (r+1)` | `p^2`: apex = atom image |
| `k`-free refinement | `s = #{a_i >= k}` | `B(A - k^s) = r + 1` | `k = 1` is not covered by the `s >= 2` inequality |
| Density race | difference `#p^2qr - pi` | sign, not the predicate `> 0` | six of the 13 predicate flips are ties |
| Sandwich endpoints | side `6k -+ 1` and `p mod 6` | `p^3 = p (mod 6)` | `8`, `27` are never endpoints |

Other inheritance and SCOPE. The sandwich census SW3 of the inherited note supplies the positive control `37915` twin-prime centers below `6*10^6`, and its class-parity law (SW4) in section SW2 (`b(6k-1)` odd, `b(6k+1)` even) is consistent with the cube side law of 4.2. [arithmetic_braids2_20260917_squarefree_symmetry.md](arithmetic_braids2_20260917_squarefree_symmetry.md) already proves that the seven proper nontrivial squarefree divisors of `p^2qr` are the seven Fano points (`S = 7` there) and that `p^3qr` has the same seven with defect four, and [arithmetic_braids2_20260917_floor_reciprocity.md](arithmetic_braids2_20260917_floor_reciprocity.md) already computes `n/rad(n)` at the three shapes; both are cited where relevant and not re-derived. The Collatz and summand companions ([collatz](arithmetic_braids_20260917_collatz.md), [summand](arithmetic_braids_20260917_summand.md)) and the [blueprint audit](collatz_blueprint_20260921_synthesis.md) are unrelated to divisor counting and are not used. SCOPE: THM-3341 (U-spine square hypotenuses, Pell orbits), THM-4139 and THM-4146 (rational three-cycle of `x^2 - 29/16`; their "divisor" is the algebraic-geometric divisor of a fibre), THM-3333 and THM-1745 were checked by title and status header only; no map from any of them into divisor counting was found, and none is claimed.

Companion: `04-computation/experiments/collatz_mod6_20260917_divisor_balance_family.py`, regenerated stdout `05-knowledge/results/collatz_mod6_20260917_divisor_balance_family.out` (SHA-256 of source `8a82bcfb11451fffee92ba1bc3cc503898d7817af9f496f82af475e4fc1c20a1`, of output `559efb0487955e8871106a4853784394fccdca06547e73a01647fe60316bb5ec`). All checks are explicit `RuntimeError`s (active under `-O`); the `-O` run is identical modulo `[time]` lines. Section 0 of the script verifies (DB1) and the `S_k` formula of section 3 against direct divisor enumeration for every `2 <= N <= 200000`. The recovered draft of the script failed one check (it bounded the prime-cube list by the sieve limit `10^7` instead of the largest endpoint `6000001`); the true census values `21` and `19` were recomputed independently and the check now asserts them.

## 1. Density of the three shapes

### 1.1 Exact counts (FINITE-EXACT)

Sieve of `omega` and `Omega` to `10^7`; `p^2qr` is `(Omega,omega) = (4,3)`, `p^3` is `(3,1)`. An independent formula
`#{p^2qr <= x} = sum_p [ #{q<r: qr <= x/p^2} - (pi(x/p^3) - [p^3 <= x/p^2]) ]`, built only from the cumulative prime count, reproduces the sieve at every scale. Both audits reproduced every entry by direct triple generation.

| `x` | primes | `p^3` | `p^2qr` | `p^2qr/primes` | `P(2) x loglog x/log x` | asym/exact |
|---|---:|---:|---:|---:|---:|---:|
| `10^5` | 9592 | 14 | 9346 | 0.9744 | 9598 | 1.0270 |
| `10^6` | 78498 | 25 | 87338 | 1.1126 | 85955 | 0.9842 |
| `10^7` | 664579 | 47 | 804249 | 1.2102 | 780006 | 0.9699 |

Prime cubes equal `pi(floor(x^(1/3)))` at each scale (checked).

### 1.2 Asymptotic ordering (PROVED conditional on CITED classical inputs; no novelty claim)

**Theorem 1.** `#{p^2qr <= x} ~ P(2) * x loglog x / log x`, `P(2) = sum_p p^-2 = 0.4522474200...`. Consequently `#{p^2qr <= x} / pi(x) ~ P(2) loglog x -> infinity`, while `#{p^3 <= x} = pi(x^(1/3)) ~ 3 x^(1/3)/log x`; the eventual order is `p^3 << primes << p^2qr`. The leading term is a special case of the classical Landau family of counts of integers with a prescribed factorization pattern; only the elementary assembly below is this note's.

*Inputs (CITED, classical; not re-verified against the sources).*
(i) Landau (1909, Handbuch): `pi_2(y) := #{qr <= y, q <= r primes} ~ y loglog y / log y`.
(ii) An elementary uniform bound `pi_2(y) = O(y loglog y / log y)` for `y >= 3`: `pi_2(y) <= sum_{q <= sqrt y} pi(y/q)`; Chebyshev gives `pi(t) <= C t/log t`, and `y/q >= sqrt y` gives `log(y/q) >= (1/2) log y`, so the sum is at most `2C (y/log y) sum_{q <= sqrt y} 1/q = O(y loglog y/log y)` by Mertens' second theorem. (The wave-one draft cited Hardy–Ramanujan 1917 here; that inequality is stated for the `omega = k` count, which dominates `#{q<r: qr <= y}`, so it would also serve, but it is retained only as UNCITED-RECOLLECTION and is not needed.)
(iii) PNT, used both for `pi(x^(1/3))` and for `pi(x)` in the ratio statement.

*Proof.* Write `#{p^2qr <= x} = sum_p S_p(x)`, `S_p(x) = #{qr <= x/p^2 : q<r, q,r != p}`, and `M(x) = x loglog x/log x`.
Tail: for `p > x^(1/4)`, `S_p(x) <= x/p^2`, and `sum_{p > x^(1/4)} x/p^2 < x * x^(-1/4) = x^(3/4) = o(M(x))`.
Head: define `f_p(x) = S_p(x) [p <= x^(1/4)] / M(x)`. For `p <= x^(1/4)`, `y = x/p^2 >= x^(1/2)`, so `log y >= (1/2) log x` and `loglog y <= loglog x`; for `x >= e^(2e)` every head `y` has `loglog y >= 1`, so the `O(1)` in (ii) is absorbed and `f_p(x) <= C'/p^2` with `C'` independent of `x`: a summable dominator. For fixed `p`, `S_p(x) = pi_2(y) - pi(sqrt y) - (pi(y/p) - [p <= y/p])`, the two corrections are `O(y/log y) = o(M(x))`, and by (i) with `loglog y/log y ~ loglog x/log x` at fixed `p`, `f_p(x) -> 1/p^2`. Tannery's theorem (dominated convergence for series) gives `sum_p f_p(x) -> sum_p 1/p^2 = P(2)`; adding the tail, `#{p^2qr <= x}/M(x) -> P(2)`. The ratio and cube statements are (iii). ∎

The script's partial sum `sum_{p <= 10^7} p^-2 = 0.4522474142` with tail `< 1.0e-07` is consistent with the literature value `0.4522474200`.

### 1.3 Exact crossover (FINITE-EXACT) and a heuristic scale (HEURISTIC)

The first `n` with `#{p^2qr <= n} > #{primes <= n}` is

```text
n = 145119 = 3 * 13 * 61^2,   counts 13433 vs 13432.
```

Bookkeeping, corrected by the recompute audit and now asserted by the script. The *predicate* `#p^2qr(n) > pi(n)` flips 13 times on `[1,10^7]`, at
`145119, 145121, 145132, 145133, 145138, 145139, 145148, 147709, 147725, 147727, 147908, 147919, 147925`,
and the flip points alternate `p^2qr, prime, p^2qr, ...` (kinds `QPQPQPQPQPQPQ`): every second flip is a prime that drops the lead from `1` to a *tie*, not a sign change. The *difference* `#p^2qr - pi` changes sign only 3 times on `[1,10^7]`: `-` to `+` at `145119`, `+` to `-` at `147739`, `-` to `+` at `147908`. Primes are strictly ahead after `145119` only on the 155 integers of `[147739,147893]`, with minimum difference `-5`; there are 73 ties after `145119`, the last at `147924`; `p^2qr` leads strictly on `[147925,10^7]`, by `139670` at `10^7`.

HEURISTIC: the naive scale from `P(2) loglog x = 1` is `x = 9196`, and the actual crossover exceeds it by the factor `15.78`; the secondary term of `pi_2` is not negligible at this range (the asymptotic overshoots at `10^5` and undershoots at `10^7`, last column of the table in 1.1).

## 2. The linear family `F = alpha S + beta U`

### 2.1 Reduction (PROVED)

For `N >= 2` with profile `(a_1,...,a_r)`:
- prime: `0 = 0`, always a solution; hence the set of integers `N` solving any cell is infinite, and "finite" below always means *finitely many exponent profiles up to permutation*;
- squarefree composite (`r >= 2`): `(1-alpha)(2^r-2) = beta r`;
- nonsquarefree: `prod(a_i+1) = T(r) := 2 + alpha(2^r-1) + beta r`, with some factor `>= 3`.

**Lemma A.** If `b_1...b_r = T` with all `b_i >= 2` and `k = #{i : b_i >= 3}`, then `(3/2)^k <= T/2^r` and `2^(r-k) | T`. *Proof.* `prod(b_i/2) = T/2^r` with each factor `>= 1` and the `k` non-2 factors `>= 3/2`; the `r-k` factors equal to `2` divide `T`. ∎

**Support bounds (PROVED).** `alpha in {0,1}`: nonsquarefree forces `3*2^(r-1) <= T(r)`, giving `r <= 0,1,2,2` for `(0,beta)`, `beta = 0..3`, and `2^(r-1) <= 1 + beta r`, giving `r <= 1,3,4,5` for `(1,beta)`. `alpha in {2,3}`, `(alpha,beta) != (2,0)`: `T/2^r <= alpha + 3r/2^r <= 4.5`, so `k <= 3`; for `r >= 4`, `c_r := beta r + 2 - alpha` satisfies `0 < |c_r| <= 3r+1 < 2^r` (it vanishes only in cell `(2,0)` or at `r = 1`), so `v_2(T) = v_2(c_r) <= log2(3r+1)`; with `2^(r-k) | T` this gives `r <= 3 + log2(3r+1)`, i.e. `r <= 7`. Within these bounds the solutions are the multiplicative partitions of `T(r)` into `r` parts `>= 2` (finite, enumerated exactly, no exponent cap).

### 2.2 The table (PROVED; every row re-checked by a 646645-profile window search with support `<= 10`, exponents `<= 12`, and by direct integer control on `2 <= N <= 200000`; both audits re-enumerated all sixteen cells independently)

| cell | finitely many profiles | `r`-bound | nonsquarefree profiles | squarefree composites | #profiles | `#N <= 2e5` |
|---|---|---:|---|---|---|---:|
| (0,0) | yes | 0 | – | none | 1 (`p`) | 17984 |
| (0,1) | yes | 1 | `p^2` | `pq` | 3 | 63214 |
| (0,2) | yes | 2 | `p^3, p^2q` | `pqr` | 4 | 68867 |
| (0,3) | yes | 2 | `p^4, p^3q` | none | 3 | 22140 |
| (1,0) | **no** | 1 | `p^2` | all | inf | 121666 |
| (1,1) | yes | 3 | `p^3, p^2qr` | none | 3 | 36372 |
| (1,2) | yes | 4 | `p^4, p^2q^2` | none | 3 | 18122 |
| (1,3) | yes | 5 | `p^5, p^2q^2r, p^2qrst` | none | 4 | 21629 |
| (2,0) | **no** | inf | `p^3 q_1...q_(r-1)`, all `r >= 1` | none | inf | 32827 |
| (2,1) | yes | 7 | `p^4, p^4q, p^2q^2rs` | none | 4 | 21840 |
| (2,2) | yes | 7 | `p^5, p^3q^2, p^5q, p^4qrs` | none | 5 | 20334 |
| (2,3) | yes | 7 | `p^6, p^6q` | none | 3 | 18493 |
| (3,0) | yes | 7 | `p^4` | none | 2 | 17992 |
| (3,1) | yes | 7 | `p^5` | none | 2 | 17989 |
| (3,2) | yes | 7 | `p^6, p^4q^2` | none | 3 | 18042 |
| (3,3) | yes | 7 | `p^7, p^3q^3r, p^7qr` | none | 4 | 18565 |

Hand check of the `T(r)` factorizations, e.g. `(2,1)`: `T = 2^(r+1) + r = 5, 10 = 2*5, 19, 36 = 2*2*3*3, 69, 134, 263` giving `p^4, p^4q, p^2q^2rs`; `(3,3)`: `T = 3*2^r - 1 + 3r = 8, 17, 32 = 2*4*4 = 2*2*8, 59, 110, 209, 404` giving `p^7, p^3q^3r, p^7qr`; `(1,3)` at `r = 5`: `T = 48 = 3*2^4` is the tight bound, forcing `(2,1,1,1,1)`. The script also confirms that `alpha in {2,3}` has no nonsquarefree solution at supports `8..12`, consistent with the proved `r <= 7`.

### 2.3 The two infinite cells and the direct laws (PROVED)

- `(1,0)`: `F = S` iff every proper nontrivial divisor is squarefree iff `N` squarefree or `N = p^2` (if `p^2 | N`, `N != p^2`, then `p^2` is a proper non-squarefree divisor).
- `(2,0)`: nonsquarefree solutions satisfy `prod(a_i+1) = 2^(r+1)`, so every `a_i+1` is a power of two, exponents `e_i >= 1` sum to `r+1`, hence exactly one `e_i = 2`: profile `(3,1^(r-1))`. Squarefree composites give `-(2^r-2) = 0`, impossible. So `F = 2S` iff `N = p` or `N = p^3 m`, `m` squarefree, `gcd(m,p) = 1` (the case `m = 1` is `p^3`): exactly one profile per support size `r >= 2`, two at `r = 1` (`p`, `p^3`). The script enumerates this family to `r = 10`.
- `(0,1)`: `F = U` iff every proper nontrivial divisor is prime iff `N in {p, p^2, pq}`.
- Prime-power law: in every cell with `alpha + beta >= 1` there is exactly one nonsquarefree prime-power solution, `p^(alpha+beta+1)` (from `a + 1 = 2 + alpha + beta`); `p` itself is of course also a prime power.

**Theorem 2 (general finiteness of the profile list, PROVED).** For all integers `alpha, beta >= 0`, the set of exponent profiles (up to permutation) solving `F = alpha S + beta U` is finite unless `(alpha,beta) in {(1,0),(2,0)}`; the set of integers `N` is infinite in every cell, since every prime is a solution. *Proof.* Squarefree composites: `alpha = 1` forces `beta = 0`; `alpha = 0` gives `2^r - 2 = beta r`, finitely many `r`; `alpha >= 2` gives a negative left side. Nonsquarefree: `alpha = 0,1` are bounded by `3*2^(r-1) <= T(r)`, resp. `2^(r-1) <= 1 + beta r`. `alpha >= 2`: `T/2^r <= alpha + beta + 2 =: c`, so `k <= log_{3/2} c` by Lemma A; if `beta >= 1` then for large `r`, `0 < beta r + 2 - alpha < 2^r` and `v_2(T) = v_2(beta r + 2 - alpha) <= log2(beta r + 2)`, so `r <= k + log2(beta r + 2)` is bounded; if `beta = 0`, `alpha >= 3`, then `v_2(T) = v_2(alpha - 2)` for `r > v_2(alpha - 2)`, so `r <= k + v_2(alpha - 2)`. ∎ (The proof-audit probed 45 cells with `4 <= alpha <= 8`, `0 <= beta <= 8` and found every solution inside these bounds; the explicit lists for `alpha, beta >= 4` are not computed here: OPEN.)

### 2.4 What is and is not special about `(1,1)` (FINITE-EXACT reading of the proved table; one draft clause REFUTED)

- The prime cube `p^3` occurs exactly in the cells with `alpha + beta = 2`: `(0,2), (1,1), (2,0)`. Among these, only `(0,2)` has an `Omega` collision (`p^3, p^2q, pqr` all have `Omega = 3`, values `1,3,3,3`). Both `(1,1)` (values `1,3,4`) and `(2,0)` (values `1` and `r+2` for `r >= 1`, listed to `r = 10` as `1,3,4,...,12`) are `Omega`-injective, since `Omega(p^3 q_1...q_(r-1)) = r + 2` determines `r` and hence the profile. So the naming of shapes by almost-prime order (`A`, `A^3`, a 4-almost-prime) is injective in `(1,1)` *and* in `(2,0)`, not in `(0,2)`. **REFUTED (wave-one draft clause):** "`(1,1)` is the unique prime-cube cell whose solution profiles have pairwise distinct `Omega`"; minimal witness `(2,0)`. The recovered script printed `(2,0): collides` only because its flag conjoined injectivity with finiteness; the flag and the check are now separated. What survives: `(1,1)` is the unique *finite* `Omega`-injective prime-cube cell.
- The `Omega`-injective cells with finitely many profiles are `(0,0),(1,1),(2,1),(2,3),(3,0),(3,1)`; including the infinite cell, the full list is `(0,0),(1,1),(2,0),(2,1),(2,3),(3,0),(3,1)` (`(1,0)` collides: `p^2` and `pq` share `Omega = 2`). Those finite injective cells with exactly three profiles are `(1,1)` and `(2,3)`.
- The tight profiles `(2,1^(r-1))` of (DB3) occur as solutions exactly in `(0,1)` at `r = 1`, `(0,2)` at `r = 2`, `(1,0)` at `r = 1`, `(1,1)` at `r = 3`, `(1,3)` at `r = 5` (solve `3*2^(r-1) = T(r)`). The draft's further bullet singling out `(1,1)` by "largest-support solution tight with `r = 3`" was tautological (it built `r = 3` into the predicate) and is dropped.

**Typed analogy, REFUTED:** "almost like a cubic". Source: a degree-3 polynomial, at most three roots. Target: the three shapes of `(1,1)`. Map: none; the equation is (DB2), `prod(a_i+1) = 2^r + r + 1`, a factorization problem indexed by `r`, not a degree. Preserved predicate: the cardinality 3 only. Lost: everything else. Sidecar: none exists. Decisive test: `(1,3)`, `(2,1)`, `(3,3)` have four profiles, `(2,2)` five, `(3,0)`, `(3,1)` two, `(2,0)` infinitely many. The two objects are unrelated.

## 3. Lattice reading and the `k`-free refinement

### 3.1 Box = cube + directions + apex (identity inherited; tight-profile reading PROVED)

Let `N` be nonsquarefree with `r` distinct primes, `L = D(N)` its divisor lattice (the exponent box), `C = D(rad N)` the Boolean subcube `2^r`, and `Q = L \ C` the order filter of non-squarefree divisors. Then `F = |L| - 2`, `S = |C| - 1`, `U = r`, so

```text
F - S - U = |Q| - (r+1),      F = S + U  iff  |L| = |C| + r + 1.
```

This identity is the inherited DB1 mechanism (`F - S` = box points outside the cube minus one; equality iff the excess is `r`), restated in lattice words; it is cited, not new. New here is the reading at the tight profile `N = p^2 m` (`m` squarefree, `p` not dividing `m`, `omega(m) = r-1`): a divisor is non-squarefree iff `p^2` divides it, so `Q = p^2 * D(m)`, isomorphic to `B_(r-1)`, and `|Q| = 2^(r-1)`. The identity holds iff `2^(r-1) = r + 1` iff `r = 3`, i.e. iff `B_2 = {bottom} ∪ {two atoms} ∪ {top}` with the three parts disjoint, the only Boolean lattice with that property (`B_0`, `B_1` have coincidences). Concretely `12 = 8 + 3 + 1` at `N = 60 = 2^2*3*5`: `Q = {4,12,20,60}`, `Q/4 = {1,3,5,15} = D(15)`; the atoms map `l -> p*lcm(p,l)`: `2 -> 4, 3 -> 12, 5 -> 20`, apex `60`. For `r >= 3` this map plus the apex injects into `Q` (the `r+1` elements `p^2, p^2 q_j, p^2 m` are distinct exactly when `r - 1 >= 2`), and the leftover `{p^2 d : d | m, 2 <= omega(d) <= r-2}` has exactly `2^(r-1) - r - 1` elements: the (DB3) minimum read as *unreached non-squarefree divisors*. For `r = 1, 2` the apex coincides with an atom image (`p^2 = N`, `p^2 q = N`), which is precisely why `p^2` and `p^2q` have `D = -1`. Hostile: `d -> d/p` is not the bijection; it sends `Q` to the upper cube face `p*D(m) = {2,6,10,30}` (an earlier draft of the check asserted the wrong map and the run refuted it). The Fano reading of the seven squarefree proper divisors of `p^2qr` is inherited from the braids2 squarefree note and not repeated.

### 3.2 `k`-free divisors (PROVED)

Let `S_k = #{d | N : 1 < d < N, d k-free}`; `S_2 = S`. Then `S_k = prod min(a_i+1, k) - 1 - [N k-free]` (verified directly for `N <= 200000`, `k <= 6`).

**Theorem 3.** For `k >= 2`: `F = S_k + U` iff `N in {p, p^(k+1), p^k q r}`, or `N = p^k q^2` when `k >= 3`. For `k = 1`: `F = U` iff `N in {p, p^2, pq}`. Hence `k = 1` and `k = 2` have three shapes each and every `k >= 3` has four; among `k >= 2`, the squarefree case `k = 2` is the only three-shape case.

*Proof.* If `N` is `k`-free then `S_k = F`, so `U = 0` and `N` is prime. Otherwise let `I = {i : a_i >= k}`, `s = |I| >= 1`, `A = prod_I (a_i+1) >= (k+1)^s`, `B = prod_{i not in I} (a_i+1)`, each factor in `[2,k]`, so `prod min(a_i+1,k) = B k^s` and the equation reads `B(A - k^s) = r + 1` (with `U = r`, valid since `N` is composite).
If `s >= 2`: `B >= 2^(r-s) >= r-s+1` and, using `k >= 2`, `A - k^s >= (k+1)^s - k^s >= 1 + s k^(s-1) >= 2s+1`, so `B(A - k^s) >= (r-s+1)(2s+1) >= r+s+1 > r+1`; no solution. (At `k = 1` the inequality `(k+1)^s - k^s >= 2s+1` fails, which is why `k = 1` is handled by the direct law of 2.3.)
If `s = 1`, with `a` the large exponent: `B(a-k+1) = r+1` and `B >= 2^(r-1) > r+1` for `r >= 4`, so `r <= 3`. `r = 1`: `a-k+1 = 2`, `N = p^(k+1)`. `r = 2`: `B in [2,k]` divides `3`, so `B = 3`, which needs `k >= 3`, and `a = k`: `N = p^k q^2`. `r = 3`: `B >= 4 = r+1` forces `B = 4`, both other exponents `1`, `a = k`: `N = p^k q r`. ∎

| `k` | solution shapes | window `r <= 8`, exp `<= 12` | `#N <= 2e5` |
|---|---|---|---:|
| 1 | `p, p^2, pq` | 3 profiles, all predicted | 63214 |
| 2 | `p, p^3, p^2qr` | 3 | 36372 |
| 3 | `p, p^4, p^3q^2, p^3qr` | 4 | 25045 |
| 4 | `p, p^5, p^4q^2, p^4qr` | 4 | 21096 |
| 5 | `p, p^6, p^5q^2, p^5qr` | 4 | 19442 |
| 6 | `p, p^7, p^6q^2, p^6qr` | 4 | 18684 |

Boundary: at `k = 3`, `F - S_3 - U` is `0` at `p^4` (`F = 3, S_3 = 2, U = 1`), `p^3q^2` (`10, 8, 2`) and `p^3qr` (`14, 11, 3`), but `-1` at `p^3` (`2, 2, 1`) and `-3` at `p^2qr` (`10, 10, 3`): the `k = 2` solutions are not `k = 3` solutions.

## 4. `Omega` blindness and the sandwich endpoints

### 4.1 Profiles per `Omega` (FINITE-EXACT table; mechanism and witness inherited)

`F, S, U` depend only on the exponent profile; the profile-to-`Omega` map is `p(k)`-to-one.

| `Omega = k` | `#profiles = p(k)` | `#F = S+U` | which | defects `D` over the class |
|---:|---:|---:|---|---|
| 1 | 1 | 1 | `p` | `0` |
| 2 | 2 | 0 | – | `-2, -1` |
| 3 | 3 | 1 | `p^3` | `-3, -1, 0` |
| 4 | 5 | 1 | `p^2qr` | `-4, 0, 1, 2` |
| 5 | 7 | 0 | – | `-5, 2, 3, 4, 5, 6` |
| 6 | 11 | 0 | – | `-6, 3, 5, 8, 9, 10, 11, 12, 15` |
| 7 | 15 | 0 | – | `-7, 4, 7, 11, 12, 13, 18, 19, 20, 24, 25, 26, 27, 33, 34` |
| 8 | 22 | 0 | – | `-8, 5, 9, 14, 16, 17, 18, 24, 27, 28, 33, 36, 39, 42, 43, 51, 56, 57, 58, 60, 70, 73` |

The witness pair `p^2q^2` (`D = 2`) and `p^3q` (`D = 1`) is inherited: both defects sit in the inherited DB2 table and inherited DB3 already states that neither `Omega` nor `omega` alone determines `D`. The only addition is the joint phrasing: the pair `(Omega, omega) = (4,2)` does not determine `D` either. In the `{A,B,C}` notation the corrected description (inherited DB3) is `A`, `{a^3}`, `{a^2 b : b = qr squarefree, gcd(a,b) = 1}`; the `C` column of the sandwich matrix mixes `p^3` (`D = 0`) with `p^2q` (`D = -1`) and `pqr` (`D = -3`).

### 4.2 Endpoint census (FINITE-EXACT) and the cube side law (PROVED)

Over centers `6k`, `1 <= k <= 10^6` (largest endpoint `6000001`):

| shape | left `6k-1` | right `6k+1` |
|---|---:|---:|
| prime | 206502 | 206345 |
| `p^3` | 21 | 19 |
| `p^2qr` | 35061 | 34876 |

Side law: for every prime `p >= 5`, `p^2 = 1 (mod 6)` so `p^3 = p (mod 6)`; the boundary primes give `8 = 2` and `27 = 3 (mod 6)`, never endpoints. Hence over all primes, `p^3 = 6k-1` iff `p = 5 (mod 6)` and `p^3 = 6k+1` iff `p = 1 (mod 6)`; this is consistent with the inherited class-parity law (SW4) (`b(p^3) = 3` odd on the left). The script checks the law for every prime with `p^3 <= 6000001` (largest cube prime `181`) and that the `21` left and `19` right cube endpoints are exactly the `p = 5` and `p = 1 (mod 6)` cubes: first left cubes `125, 1331, 4913, 12167`; first right cubes `343, 2197, 6859, 29791`. Centers with both endpoints `F = S + U` solutions: `55613`, of which `37915` twin primes (matches the inherited SW3 cell) and `17698` mixed. Connection contract: source = solution shapes; target = sandwich endpoints; map = restriction to `6k -+ 1`; preserved = shape and side; lost = the pairing across the center; sidecar = `p mod 6` for cubes; decisive test = the cube side law. No asymptotic claim is made.

## Reproduction

```bash
python3 04-computation/experiments/collatz_mod6_20260917_divisor_balance_family.py > 05-knowledge/results/collatz_mod6_20260917_divisor_balance_family.out
python3 -O 04-computation/experiments/collatz_mod6_20260917_divisor_balance_family.py   # identical modulo [time] lines
```

Run from the repository root (the script uses no relative paths). About 36 s on one core (the `.out` carries its own `[time]` stamps) and under 750 MB peak RSS as measured with `/usr/bin/time -l` on the mac-mini (numpy sieve to `10^7`, pure-python direct audits to `2*10^5`). Independent audit companions from wave one: `04-computation/experiments/collatz_mod6_20260917_divisor_balance_family_audit_recompute.py` and `..._audit_proof-audit.py`.

## Stopping boundary / next question

The divisor-balance family is closed at three combinatorial levels: the `4x4` linear grid with general finiteness of the profile list, the `k`-free refinement for all `k`, and the lattice reading of the (DB3) minimum. What remains is analytic, not combinatorial: the secondary term `c` in `#{p^2qr <= x} = (x/log x)(P(2) loglog x + c + o(1))`, and an effective bound (explicit error terms for `pi_2` and `pi`) certifying that the three sign changes at `145119, 147739, 147908` and the last tie at `147924` are all of them. Repeating the census at a larger cutoff without such a term would not cross that boundary. A smaller combinatorial follow-up is the explicit profile lists for `alpha, beta >= 4` and the joint defect `D_k = F - S_k - U` as an inverse problem in the style of inherited DB2.
