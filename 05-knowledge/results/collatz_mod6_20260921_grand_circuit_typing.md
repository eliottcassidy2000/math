# Grand circuit typing: Giuga, Ankeny-Artin-Chowla, Littlewood and Sarnak against the mod-6 Collatz framework

**Status:** CITED: the four conjectures, stated exactly with sources (sections 1-4); the Bernstein-Lagarias 2-shift conjugacy; Einsiedler-Katok-Lindenstrauss; the AAC verification bound. PROVED (scoped, elementary): the Borwein-Borwein-Borwein-Girgensohn equivalence re-derived from the prime-power power-sum lemmas (S2); the user's criterion `p^2(p-1) | n-p` is EQUIVALENT to it (S4); the forced valuations `v_3(u)=0`, `v_2(u)<=1` of the AAC unit from the norm equation (S7b); the one-line "unbounded partial quotients imply Littlewood for every pair" (S11); the explicit Mobius-correlated 2-adic point (S15) and the periodicity argument for the finite-`j` correlations (S16); the box ring is the shift `x -> x+1` (S18); coprimality along fibres and odd orbits (S20-S21). FINITE-EXACT: the Giuga sum for every `n<=30000` and the criteria for `n<=10^5`; PARI `quadunit` for all 211 primes `p = 1 mod 4` below 3000; the first 60 partial quotients of `log_2 3`; the `N=10^6` Mobius correlations for `j<=6`; the repository greps. HEURISTIC: one float64 minimum (S12). REFUTED: `R(x)=4x+1` as a box-ring multiplication (S19). SCOPE / NO MAP: every proposed bridge from the four conjectures to the session's objects (S6, S9, S13, S17, S23, S24). OPEN: all four conjectures, and whether `log_2 3` has unbounded partial quotients. No novelty claim for any of the number theory; every proof here is textbook material re-derived so that the typing is self-contained. Results note of session `collatz-mod6-20260917` (machine `mac-mini`), wave 2026-09-21, lane `grand_circuit_typing`; not a reserved canon ID.

## Inheritance and concept board

The pasted "grand unified circuit" asserts that the session's Collatz framework proves, or is proved by, the Giuga, Ankeny-Artin-Chowla, Littlewood and Sarnak conjectures, through "`S_2 x S_3` projective lift of `S^6`", "`6/pi^2` square-free mask", "`T_{n-2}` tournament compression", "Wythoff Bragg peaks", "log-cosine Lyapunov energy", "23 mod 256 valve", "17-vertex tournament capacity" and "`B^3=-I`". Most of these tokens were already typed and are only cited here: the dyadic energy and the `S6` monodromy are REFUTED in the [blueprint synthesis](collatz_blueprint_20260921_synthesis.md) and its [energy companion](collatz_blueprint_20260921_energy.md); "`T_{n-2}` path edges carry zero bits" and the `4k+1` square peaks are REFUTED in the [scaffolding audit](collatz_mod6_20260917_scaffolding_audit.md) (its section 5 gives the true identity `sum_T h(T) = n! 2^{T(n-2)}` and section 8 the residue-stratified null); `6/pi^2` has no every-orbit content ([blueprint](collatz_blueprint_20260921_synthesis.md), [guards squarefree note](collatz_guards_20260921_squarefree.md)); the exact reset cylinder is `n = 23 mod 256`, not a "valve" ([guards valves](collatz_guards_20260921_valves.md)); `B^3=-I` is the order-six lift of the rational 3-cycle of `x^2-29/16` in [THM-4139](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md) and the wave-one [zsigmondy lane](collatz_mod6_20260917_zsigmondy_triad.md); "`F=S+U`" is the divisor identity of the [divisor_balance lane](collatz_mod6_20260917_divisor_balance_family.md) (`N in {p, p^3, p^2qr}`); the "4-vertex tournament matrix" is the score-profile object of the [cell_ordering lane](collatz_mod6_20260917_cell_ordering_scale.md); `C_2 x S_3` is the actual automorphism action on the fruit-curve triple in the [catalan_elliptic synthesis](catalan_elliptic_20260921_synthesis.md); the Pell hypotenuses are units of `Z[sqrt 2]` ([THM-3341](../../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md)); gcd-1 along orbit edges is the [odd_square edge encoding](odd_square_20260921_edges.md); the positive-entropy conjugacy of the greedy 3-adic map `G` is Theorem 1.3 of the [three_adic lane](collatz_mod6_20260917_three_adic_g_map.md); the continued fraction of `log_2 3` and the cycle gate `n_0 = bB/(2^K-3^L)` along its convergents are the object of the concurrent [pillai lane](collatz_mod6_20260921_pillai_convergents_cycle_gates.out) and of the [arithmetic braids sessions](arithmetic_braids2_20260917_synthesis.md). **Closest proved mechanism:** the Bernstein-Lagarias conjugacy of the 2-adic Collatz map with the full one-sided 2-shift (CITED), which decides the Sarnak question in one line (positive entropy, hypothesis void) and, by surjectivity of the parity map, even produces an explicit 2-adic point whose parity sequence is `[mu = 1]` (S15). **Canonical hostile:** the AAC norm equation `t^2 - p u^2 = -4` itself: reduced mod 3 and mod 16 it forbids `3 | u` and `4 | u` for every `p`, so the "large 2-adic or 3-adic valuation" the paste needs to link `u` to the session's 2/3-adic objects cannot occur for any prime (S7b, checked on all 211). **Corrected near miss:** the user's Giuga criterion `p^2(p-1) | n-p for every prime p | n` looks like a strengthening (it has a `p^2`) of the BBBG criterion, but `n-p = p(n/p-1)` turns it into exactly `p(p-1) | n/p-1`, which is BBBG with squarefreeness absorbed (S4); it is correct, and it is not new. **Least-used sidecar:** the parity of `T^j n` is periodic mod `2^(j+1)` (S16, checked `j<=6` at `N=10^6`), so every finite-`j` Mobius correlation of the parity is a sum over residue classes and tends to 0 by the prime number theorem in arithmetic progressions; that is the only "Mobius-Collatz" statement that is true, and it is a statement about the odometer on `Z/2^(j+1)`, not about Collatz. **No map found (SCOPE):** "`S^6` with `S_2 x S_3` projective lift", "17-vertex tournament capacity", "Wythoff Bragg peaks", "`5 pi/6` phase from `5 mod 6`", "silver ratio", "repunit prime breaks", "reversed Hamiltonian edge", "trapping loop" (S13, S17, S23, S24).

Typed analogy (the only one with a map), **Sarnak on `Z_2` -> Sarnak on the odometer**: source = the Collatz map `T` on `Z_2` with `f = parity`; target = the rotation `n -> n+1` on `Z/2^(j+1)` with `f_j = parity(T^j .)`; map = truncation of the parity vector to its first `j+1` coordinates; preserved = the finite-`N` correlation sums `c_j` exactly (S16); lost = the dynamics (the target has entropy 0, the source `log 2`); sidecar = `j`; test = `c_j -> 0` for every fixed `j` (PROVED via PNT in APs) versus disjointness for the whole system (FALSE, S15). Everything else in the paste: no map found.

## 1. Giuga

**Statement (Giuga 1950; CITED).** `n >= 2` is prime iff `sum_{i=1}^{n-1} i^(n-1) = -1 (mod n)`. The "if" is the conjecture; the "only if" is Fermat.

**Equivalence (Borwein, Borwein, Borwein, Girgensohn, Amer. Math. Monthly 103 (1996); CITED).** A composite `n` satisfies the congruence iff `n` is squarefree and, for every prime `p | n`, `(p-1) | (n/p-1)` (Korselt, so `n` is Carmichael) and `p | (n/p-1)` (Giuga).

**S1 (FINITE-EXACT).** The sum was computed directly, by vectorised binary exponentiation mod `n`, for every `n <= 30000`. Every prime gives `-1`; no composite does; the direct congruence and the BBBG criterion agree on every composite (0 disagreements).

**S2 (PROVED; re-derivation of BBBG).** Write `q = p^a || n` and `k = n-1 >= a`. Since `i^k mod q` depends only on `i mod q`, `sum_{i=1}^{n-1} i^k = (n/q) sum_{i=0}^{q-1} i^k (mod q)`, and the non-units contribute 0 because `k >= a`. For `a = 1` the unit sum is `-1` if `(p-1) | k` and `0` otherwise (primitive root). For `a >= 2` and `p` odd, the units are `C_{p-1} x (1+pZ/p^a)`; `k` is coprime to `p` (as `n-1` is), so `x -> x^k` permutes `1+pZ`, whose sum is `p^(a-1) (mod p^a)`, hence the unit sum is `0 mod p^(a-1)` and never `-1`. For `p = 2`, `a >= 2`, `k` is odd and permutes the units, whose sum is `2^(2a-2) = 0 mod 2^a`. So the congruence forces `n` squarefree, `(p-1) | (n-1)` and `n/p = 1 (mod p)` for every `p | n`; as `n-1 = n/p-1 (mod p-1)`, the first is `(p-1) | (n/p-1)`. QED. Both lemma tables were brute-forced (primes `p < 200` for all `k < p`; prime powers `q = p^a <= 3000`, `a >= 2`): 0 failures.

**S3 (FINITE-EXACT, `n <= 100000`).** Carmichael numbers: 16, namely 561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341, 41041, 46657, 52633, 62745, 63973, 75361. Giuga numbers (squarefree, `p | n/p-1` for all `p | n`): 30, 858, 1722, 66198, with `sum 1/p - 1/n = 1` in each case (`30 = 2*3*5`, `858 = 2*3*11*13`, `1722 = 2*3*7*41`, `66198 = 2*3*11*17*59`). Composites satisfying BBBG (Carmichael and Giuga): none. Composites satisfying the user's criterion: none. Composites where the two criteria disagree: 0. All primes below 2000 satisfy the user's criterion (vacuously, `n-p = 0`).

**S4 (PROVED; the user's criterion is EQUIVALENT to BBBG).** `n-p = p(n/p-1)`, so `p^2(p-1) | n-p` iff `p(p-1) | n/p-1` iff `p | n/p-1` and `(p-1) | n/p-1`, since `gcd(p, p-1) = 1`. Squarefreeness is implied: if `p^2 | n` then `p | n/p`, so `p` does not divide `n/p-1`. For prime `n` the criterion is vacuous. The paste's restatement is therefore correct, is the BBBG criterion in one line, and involves no Collatz object.

**S5 (FINITE-EXACT).** The trunk numbers `t_j = (4^j-1)/3 = R^j(1)` of the inherited inverse-fibre braid (the odd predecessors of the powers of two) are never Korselt and never Giuga for `3 <= j <= 40` (0 hits). Table for the small cases:

| `j` | `t_j` | squarefree | Korselt | Giuga | first failing prime |
|---|---|---|---|---|---|
| 3 | 21 | True | False | False | 7 |
| 4 | 85 | True | False | False | 5 |
| 5 | 341 | True | False | False | 11 |
| 6 | 1365 | True | False | False | 3 |
| 7 | 5461 | True | False | False | 43 |
| 8 | 21845 | True | False | False | 5 |
| 9 | 87381 | False | False | False | 3 |
| 10 | 349525 | False | False | False | 5 |
| 11 | 1398101 | True | False | False | 23 |
| 12 | 5592405 | True | False | False | 7 |
| 20 | 366503875925 | False | False | False | 5 |
| 30 | 18 digits | False | False | False | 5 |
| 40 | 24 digits | False | False | False | 5 |

(341 and 5461 are Cipolla's base-2 pseudoprimes, inherited from the [wild_typing lane](collatz_mod6_20260917_wild_typing.md); pseudoprime to base 2 is far weaker than Carmichael.)

**S6 (SCOPE / NO MAP).** `gcd(n, R(n)) = gcd(n, 4n+1) = 1` and `gcd(x_i, x_{i+1}) = 1` along odd orbits are facts about pairs; Giuga's criterion is a divisibility `p | n/p-1` inside one `n`. No session object carries it. The proposed "`T_{n-2}` bits versus Giuga factors" has no map: `T(n-2) = C(n-1, 2)` is an edge count of a tournament on `n` vertices (scaffolding audit, section 5) and is not attached to any integer `n` as a factorisation datum.

## 2. Ankeny-Artin-Chowla

**Statement (Ankeny, Artin, Chowla, Ann. of Math. 56 (1952); CITED).** For a prime `p = 1 (mod 4)` with fundamental unit `(t + u sqrt p)/2` of `Q(sqrt p)`, `p` does not divide `u`.

**S7 (FINITE-EXACT; PARI `quadunit`, ring of discriminant `p`).** 211 primes `p = 1 mod 4` in `[5, 3000)`. The norm `t^2 - p u^2` is `-4` for all 211 (`+4` for none). AAC violations (`u = 0 mod p`): none. The longest `u` is at `p = 2689` (38 digits), with `u mod p = 2646`.

| valuation | histogram over the 211 primes | "random integer" expectation |
|---|---|---|
| `v_3(u)` | `{0: 211}` | `v_3 >= 1`: 70.3, `v_3 >= 2`: 23.4 |
| `v_2(u)` | `{0: 79, 1: 132}` | `v_2 >= 1`: 105.5, `v_2 >= 2`: 52.8 |

**S7b (PROVED; why the histograms are degenerate).** The fundamental unit of `Q(sqrt p)`, `p = 1 mod 4` prime, has norm `-1` (CITED classical, Legendre/Dirichlet; re-seen on all 211), so `t^2 - p u^2 = -4`. Mod 3: `3 | u` would give `t^2 = 2 mod 3`, impossible, so `v_3(u) = 0` always. Mod 16: `4 | u` would give `t^2 = 12 mod 16`, impossible, so `v_2(u) <= 1` always. Mod 8: `t, u` both odd forces `1 - p = -4 mod 8`, i.e. `p = 5 mod 8`; so `p = 1 mod 8` forces `v_2(u) = 1`. Census: 101 primes `p = 1 mod 8`, all with `v_2(u) = 1`; 110 primes `p = 5 mod 8`, of which 79 have `v_2(u) = 0` and 31 have `v_2(u) = 1`. Hence "large 2-adic or 3-adic valuation of `u`" cannot occur for any `p`: the only conceivable bridge to the session's 2-adic and 3-adic objects is closed by the norm equation itself.

**S8 (CITED / UNCITED-RECOLLECTION).** AAC is verified for all `p < 2*10^11` (van der Poorten, te Riele, Williams, Math. Comp. 70 (2001), with a 2003 corrigendum; CITED). A 2024 preprint claiming a counterexample is an UNCITED-RECOLLECTION and is not asserted here; the conjecture is treated as OPEN.

**S9 (SCOPE / NO MAP).** No object of the session is a unit of a real quadratic field of prime discriminant `p = 1 mod 4`. The Pell hypotenuses of THM-3341 are units of `Z[sqrt 2]` (discriminant 8), and `p = 2` is excluded from AAC. "Silver ratio", "Fermat numbers `2^(2^r)+1`", "repunit prime breaks", "reversed Hamiltonian edge": no map found.

## 3. Littlewood

**Statement (Littlewood, c. 1930; CITED).** For all real `alpha, beta`, `liminf_n n ||n alpha|| ||n beta|| = 0`. Einsiedler, Katok, Lindenstrauss, Ann. of Math. 164 (2006) (CITED): the set of exceptional pairs has Hausdorff dimension 0.

**S10 (FINITE-EXACT; PARI at 400 digits).** First 60 partial quotients of `log_2 3`:
`[1, 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, 15, 1, 9, 2, 5, 7, 1, 1, 4, 8, 1, 11, 1, 20, 2, 1, 10, 1, 4, 1, 1, 1, 1, 1, 37, 4, 55, 1, 1, 49, 1, 1, 1, 4, 1, 3, 2, 3, 3, 1]`.
Record quotients `(index, value)`: `(0, 1), (3, 2), (5, 3), (7, 5), (9, 23), (14, 55)`; the running maximum is 55 from index 14 on (55 recurs at index 46). Sum 388, mean 6.467. Convergents `p/q` with `q <= 10^12` and the quality `q ||q alpha||`:

| `p/q` | `q ||q alpha||` | `p/q` | `q ||q alpha||` |
|---|---|---|---|
| 1/1 | 0.584963 | 301994/190537 | 0.017732 |
| 2/1 | 0.415037 | 16785921/10590737 | 0.799110 |
| 3/2 | 0.339850 | 17087915/10781274 | 0.189869 |
| 8/5 | 0.375937 | 85137581/53715833 | 0.269110 |
| 19/12 | 0.234600 | 272500658/171928773 | 0.443799 |
| 65/41 | 0.678036 | 357638239/225644606 | 0.547999 |
| 84/53 | 0.159665 | 630138897/397573379 | 0.060709 |
| 485/306 | 0.451282 | 9809721694/6189245291 | 0.854706 |
| 1054/665 | 0.041881 | 10439860591/6586818670 | 0.096198 |
| 24727/15601 | 0.409514 | 103768467013/65470613321 | 0.435645 |
| 50508/31867 | 0.334001 | 217976794617/137528045312 | 0.178303 |
| 125743/79335 | 0.419450 | 1193652440098/753110839881 | 0.129252 |
| 176251/111202 | 0.577584 | | |

**S11 (PROVED, one line).** If `alpha` has unbounded partial quotients then for every `beta`, `liminf n ||n alpha|| ||n beta|| <= liminf_i q_i ||q_i alpha|| * (1/2) <= liminf_i 1/(2 a_{i+1}) = 0`, using `||q_i alpha|| < 1/q_{i+1} <= 1/(a_{i+1} q_i)` and `||n beta|| <= 1/2`. Whether `log_2 3` has unbounded partial quotients is OPEN, so the pair `(log_2 3, phi)` is OPEN; and `phi = [1; 1, 1, ...]` has bounded quotients, so `phi` contributes nothing to any such argument.

**S12 (HEURISTIC, float64).** `min_{n <= 10^6} n ||n log_2 3|| ||n phi|| = 0.000209` at `n = 10946` (the Fibonacci number of index 21, as expected: Fibonacci `n` make `||n phi||` small); `min n ||n log_2 3|| = 0.017734`, `min n ||n phi|| = 0.381966` (`phi`'s Hurwitz constant is `1/sqrt 5 = 0.447214`). This is numerology, not evidence.

**S13 (SCOPE / NO MAP).** The Collatz discrepancy `K_j - j log_2 3` involves one irrational, an inhomogeneous one-number problem; Littlewood is a two-number simultaneous problem. Wythoff, Zeckendorf and the golden ratio occur in none of the audited notes (in the [pythagorean_semicircle lane](collatz_mod6_20260917_pythagorean_semicircle.md) `phi` is an angle). "Bragg peaks `q_{m,n} = (2 pi/phi^2)(m + n phi)`" is the diffraction module of the Fibonacci chain (Levine-Steinhardt 1984, CITED as a formula family; the paste's normalisation is not load-bearing); "`5 pi/6` phase from `5 mod 6`": no object.

## 4. Sarnak

**Statement (Sarnak 2009/2010; CITED).** For every topological dynamical system `(X, T)` of zero entropy, every `f in C(X)` and every `x in X`, `(1/N) sum_{n <= N} mu(n) f(T^n x) -> 0`.

**S14 (CITED).** The map `T(x) = x/2` (`x` even), `(3x+1)/2` (`x` odd) on `Z_2` is topologically conjugate, through the parity-vector map `Q`, to the one-sided full 2-shift (Lagarias 1985; Bernstein-Lagarias 1996), so `h_top = log 2 = 0.693147`. The greedy 3-adic map `G` of this session is conjugate to a 6-state SFT of entropy `log 3 = 1.098612` (three_adic lane, Theorem 1.3, PROVED there). Both have positive entropy: Sarnak's hypothesis does not apply, and no statement in either direction follows.

**S15 (PROVED + FINITE-EXACT; disjointness is FALSE on `Z_2`).** Since `Q` is a bijection, there is a 2-adic `x` with `parity(T^k x) = [mu(k+1) = 1]` for all `k >= 0`. The truncation `x = 890931253 mod 2^32` realises the first 32 values `10000100010001100000110001000000` of `[mu = 1]` (verified by direct iteration). For that `x` and `f = parity`, `(1/N) sum_{n <= N} mu(n) f(T^n x) = (1/N) #{n <= N : mu(n) = 1} -> 3/pi^2 = 0.303964`, not 0. So the Collatz system on `Z_2` is not Mobius-disjoint, exactly as any full shift is not; the point `x` is a 2-adic irrational.

**S16 (FINITE-EXACT, `N = 10^6`; then PROVED).** `c_j = (1/N) sum_{n <= N} mu(n) (-1)^{parity(T^j n)}`:

| `j` | `c_j` | parity mean | `parity(T^j n)` periodic mod `2^(j+1)` |
|---|---|---|---|
| 0 | -0.000068 | 0.500000 | True |
| 1 | +0.000344 | 0.500000 | True |
| 2 | +0.000306 | 0.500000 | True |
| 3 | +0.000086 | 0.500000 | True |
| 4 | -0.000054 | 0.500000 | True |
| 5 | +0.000600 | 0.500000 | True |
| 6 | -0.000760 | 0.499992 | True |

PROVED: `parity(T^j n)` depends only on `n mod 2^(j+1)` (the first `j+1` parities determine and are determined by `n mod 2^(j+1)`), so `(-1)^{parity(T^j n)}` is a finite combination of residue-class indicators, each of which is Mobius-orthogonal by the prime number theorem in arithmetic progressions (Landau; CITED). Hence `c_j -> 0` for every fixed `j`. This is the zero-entropy odometer on `Z/2^(j+1)`, not the Collatz dynamics, and it proves nothing about Collatz.

**S17 (SCOPE / NO MAP).** The proposed "Mobius-Collatz matrix" is not defined; the objects the paste names (`F = S + U`, the 4-vertex tournament matrix, `B^3 = -I`) are, respectively, a divisor identity, a score-profile object and a rational-3-cycle lift, and none is a dynamical system on which `mu` acts; no "trapping loop" is defined.

## 5. The box ring and the coprime fibres

**S18 (PROVED).** With `a boxtimes b := ab + a + b` one has `(a boxtimes b) + 1 = (a+1)(b+1)` (checked on `[-10, 10]^2`, 0 failures). With `a boxplus b := a + b + 1`, the map `x -> x + 1` is a ring isomorphism `(Z, boxplus, boxtimes) -> (Z, +, *)`: zero is `-1`, one is `0`, `a boxtimes a = (a+1)^2 - 1`. Transport preserves everything and adds nothing.

**S19 (REFUTED as a map).** `R(x) = 4x + 1` is not `c boxtimes x = (c+1)x + c` for any `c` (it would need `c+1 = 4` and `c = 1`); no solution in `[-50, 50]`.

**S20 (PROVED; inherited braid).** `x_j = 4^j x_0 + (4^j-1)/3 = R^j(x_0)` and `gcd(x_j, x_{j+1}) = gcd(x_j, 4x_j + 1) = 1` (max over `j < 30`, `x_0 = 7`: 1).

**S21 (PROVED).** Consecutive odd iterates `x, x' = (3x+1)/2^k` satisfy `gcd(x, x') | gcd(x, 3x+1) = 1` (checked on 326623 consecutive-odd-iterate pairs, max gcd 1). Both are odd, so the pair is a point of the odd-coprime chart of the pythagorean_semicircle lane and determines one primitive triangle after ordering; that is the odd_square lane's edge encoding. For Giuga it is content-free (S6).

## 6. Typing greps: "17-vertex tournament capacity" and "`S^6 / S_2 x S_3` projective lift"

**S22 (FINITE-EXACT; repository grep).** 182 canon theorems have "tournament" in their title; 10 of them have "17" in the title. In each the 17 is a coefficient of the 7-tournament spectrum `x(x^2+7)(x^4+14x^2+17)` (THM-1450, THM-1455, THM-1575), the Fermat prime 17 (THM-871), the fingerprint twin `(17, 13)` (THM-893), or a mod-16 / 2-adic statement (THM-1780, THM-1862, THM-1865, THM-1870, THM-1980). Tournament-titled theorems whose body mentions a 17-vertex, order-17 or `n = 17` tournament: 0. Results notes mentioning `S^6`, `S6 monodromy` or `S_2 x S_3`: 5, all audits of the same paste family.

**S23 (SCOPE).** The only 17s in tournament context in this session are the `h`-value 17 attained at `n = 6` (scaffolding audit, section 5; the `h`-spectrum is the odd numbers minus `{7, 21}`, [THM-1745](../../01-canon/theorems/THM-1745-leaf-graded-arborescence-filtration-721-shadow.md)) and the class 17 of the mod-30 wheel. No "capacity" of a 17-vertex tournament is defined anywhere: no object.

**S24 (SCOPE).** `S^6` (the 6-sphere) appears in no note. `S_2 x S_3 = C_2 x S_3` (order 12) appears once, as the actual automorphism action on the fruit-curve triple (catalan_elliptic synthesis, section 4). A "projective lift of `S^6`" is undefined; the blueprint audit already found no `S6` monodromy. No object.

## 7. The honest bridge: `2^K - 3^L` along the convergents of `log_2 3`

**S25 (FINITE-EXACT).** A convergent `p/q` of `log_2 3` gives `2^p ~ 3^q`, i.e. `(K, L) = (p, q)`; the inherited cycle gate `n_0 = bB/(2^K - 3^L)` needs `2^K - 3^L` small and positive relative to `B`.

| `K` | `L` | sign of `2^K - 3^L` | `|2^K - 3^L| / 3^L` | digits of `|2^K - 3^L|` |
|---|---|---|---|---|
| 1 | 1 | - | 3.333e-01 | 1 |
| 2 | 1 | + | 3.333e-01 | 1 |
| 3 | 2 | - | 1.111e-01 | 1 |
| 8 | 5 | + | 5.350e-02 | 2 |
| 19 | 12 | - | 1.346e-02 | 4 |
| 65 | 41 | + | 1.153e-02 | 18 |
| 84 | 53 | - | 2.086e-03 | 23 |
| 485 | 306 | + | 1.023e-03 | 144 |
| 1054 | 665 | - | 4.365e-05 | 313 |
| 24727 | 15601 | + | 1.819e-05 | 7439 |
| 50508 | 31867 | - | 7.265e-06 | 15200 |
| 125743 | 79335 | + | 3.665e-06 | 37847 |
| 176251 | 111202 | - | 3.600e-06 | 53052 |
| 301994 | 190537 | + | 6.451e-08 | 90903 |

CITED: Baker-type lower bounds `|2^K - 3^L| > 2^K / K^C` (Pillai's problem; Stroeker-Tijdeman 1982 for this pair) are the only formal input that bounds cycle lengths through the gate; that is the concurrent [pillai lane's](collatz_mod6_20260921_pillai_convergents_cycle_gates.out) object, not this lane's.

## 8. Verdict table

| Conjecture | Statement | Session-side computation | Proposed bridge | Typing |
|---|---|---|---|---|
| Giuga | CITED (Giuga 1950; BBBG 1996) | FINITE-EXACT `n <= 30000` direct, `n <= 10^5` by criterion; BBBG re-PROVED (S2); user's `p^2(p-1) \| n-p` PROVED equivalent (S4) | `T_{n-2}` bits / coprime fibres | NO MAP (S6) |
| AAC | CITED (1952); verified `p < 2*10^11` CITED | FINITE-EXACT 211 primes `p < 3000`, no violation; `v_3(u) = 0`, `v_2(u) <= 1` PROVED (S7b) | Pell / Fermat / repunits / silver ratio | NO MAP (`p = 2` excluded; S9) |
| Littlewood | CITED (c. 1930; EKL 2006) | FINITE-EXACT 60 quotients, max 55; one-liner PROVED (S11); `(log_2 3, phi)` OPEN | Wythoff / Bragg peaks / `(log_2 3, phi)` vectors | NO MAP (S13) |
| Sarnak | CITED (2009/2010) | entropy `log 2`, `log 3 > 0` CITED/PROVED (S14); disjointness FALSE on `Z_2` (S15); `c_j -> 0` trivially (S16) | Mobius-Collatz matrix / trapping loop | NO MAP (S17) |
| "17-vertex capacity" | no object | 0 body hits among 182 tournament theorems (S22) | tournament compression | SCOPE (S23) |
| "`S^6` / `S_2 x S_3` lift" | no object | 5 notes, all audits (S22) | projective lift | SCOPE (S24) |

**Answer to the closing question.** None of the three proposed formalisations (a Mobius-Collatz matrix; simultaneous approximation of the vector `(log_2 3, phi)`; `T_{n-2}` bits against Giuga factors) has a map from a defined session object to the conjecture it names. The honest, formalisable bridges are the Baker/Pillai lower bounds on `|2^K - 3^L|` through the cycle gate (S25; pillai lane) and the E-graph / greedy-`G` mirror of this session (three_adic and extended_collatz lanes), and neither touches any of the four conjectures.

## Reproduction

```
python3 04-computation/experiments/collatz_mod6_20260921_grand_circuit_typing.py > 05-knowledge/results/collatz_mod6_20260921_grand_circuit_typing.out
```

Dependencies: PARI/GP (`gp -q` on stdin, for `quadunit` and the 400-digit continued fraction), numpy (vectorised power sums and the Mobius sieve), sympy (`factorint`, `mobius`, `primerange`; every factorisation is re-multiplied). Every check is an explicit `raise`; the `python3 -O` run is identical. Runtime about two minutes, memory far below 1 GB. Every number quoted in this note appears in the `.out`.

## Stopping boundary / next question

Stopped at the typing boundary: the four conjectures are stated, their session-side content is exhausted (each is either decided against the paste's mechanism by a one-line argument, or has no object to attach to), and no computation beyond the stated ranges would change a verdict. Not attempted: any extension of the Giuga search beyond `10^5` (the literature bound is far larger and is not this lane's business), any statement about the alleged 2024 AAC counterexample, or any diophantine work on the partial quotients of `log_2 3`. The one live next question this lane exposes is the pillai lane's: the sign pattern of `2^K - 3^L` along the convergents alternates, so only every other convergent is a candidate `(K, L)` for the positive gate, and the Stroeker-Tijdeman bound on those is the quantity to compare with `B`.
