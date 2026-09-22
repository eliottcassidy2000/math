# Signed sandwich discrepancy: class drift versus pairing fluctuation

**Status:** PROVED (scoped, elementary): Lemma 2.1 (class-labelled exponent shape is determined by `(Omega, omega, b)` for `Omega<=3`, fails at `Omega=4`), Theorem 2.2 (semiprime class-bias identity, all `x>=1`), Lemma 2.3 (finite-modulus side symmetry, a restatement of inherited SW2). FINITE-EXACT: census on `10^7` centers with four shifted-pair controls, the nine mod-35 strata, the prime race, and a pooled-gap control with a block-level residual test. CITED: Chebyshev-bias literature (Rubinstein-Sarnak, Ford-Sneed, Meng, Bays-Hudson). HEURISTIC: the product-of-marginals drift explains the pooled data cleanly at `K=10^6` and only partially at `K=10^7`. REFUTED: the pointwise sign law `N_PS(K)>=N_SP(K)` (minimal witness `K=1453`). OPEN: the three start-class-1 residual deviations at `10^7`; the order of the pairing fluctuation; `S_1(x)>S_5(x)` for all large `x`. No novelty claim for the class-parity law (inherited SW4), for the product-of-marginals algebra, or for the Chebyshev-type biases (classical). No asymptotic theorem, no twin-prime or Chen-type statement is claimed. Results note of session collatz-mod6-20260917 (machine mac-mini), lane `sandwich_bias`, finalized 2026-09-21 from the recovered draft after two independent audits; not a reserved canon ID.

## Inheritance and concept board

The object and its exact local structure come from [the divisors note](arithmetic_braids_20260917_divisors.md): the ordered matrix `N_ij(K)` (SW1), the CRT sandwich identity `sum x^{h-} y^{h+} = prod (p-2+x+y)` with its negative covariance, the negation symmetry `k -> -k` that swaps sides (SW2), the `K=10^6` 4x4 table (SW3, reproduced here entry by entry as a gate), the class-parity law `b(6k-1)` odd / `b(6k+1)` even (SW4), and the Chen/Pintz literature audit. That note refuted exact equality `N_12(K)=N_21(K)` at `K=4` and left as the next step the signed discrepancy with the class parity retained; this note answers it. The [squarefree-symmetry note](arithmetic_braids2_20260917_squarefree_symmetry.md) is used only for typing: the sandwich matrix is an outcome table, not a tournament. **Closest proved mechanism:** SW4 (the sign of `n mod 6` is `(-1)^b`), which is what makes semiprimes on the right same-class-or-square and on the left mixed; everything class-driven here is bookkeeping on top of it. **Canonical hostile:** the class-reversed pair `(6k+1,6k+5)`, whose `N_PS-N_SP=+461` at `K=10^7` has the wrong sign against a predicted `-229.57`; together with the witness `K=1453` it fixes the scale of the pairing fluctuation. **Corrected near misses (found by the audits):** the draft's seven-stratum mod-35 partition omitted `35|left` and `35|right` and was off by the SP pair `(35,37)`; the headline "prime squares carry a third of the drift" is convention-dependent; the pooled-gap control is a consistency check, not decisive, because gap-pooled means are, block by block, the local product of marginals. **Least-used sidecar:** the class-labelled exponent shape (which residue class carries the repeated prime), invisible in `(Omega, b)` and the reason Lemma 2.1 stops at `Omega=3`. **Typed analogy (Chebyshev bias -> sandwich sign):** source = the prime race `pi_chi(x)` and the semiprime race `S_1-S_5`; target = `N_PS(K)-N_SP(K)`; map = product of marginals `(R_1C_2-R_2C_1)/K`; preserved = the sign of the expected drift at `K=10^5,10^6,10^7` and the block-level tracking (corr 0.99-1.00); lost = the pointwise sign (`K=1453`, `K=2*10^5`) and any pairing correlation; sidecar = the block residual `r_(g,b)`; test = the residual `z` of section 5(iv). No map found (SCOPE) from this lane to [THM-4139](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md), [THM-4146](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md), [THM-3341](../../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md) (Gaussian squaring, Pell hypotenuses), THM-3333, THM-1745, the [Collatz](arithmetic_braids_20260917_collatz.md) inverse braid, or the [blueprint audit](collatz_blueprint_20260921_synthesis.md): the sandwich question shares only the modulo-6 typing with them, and no connection is manufactured.

Notation. `Omega`, `omega`, `b(n)` (number of prime factors `5 mod 6` with multiplicity); `chi(n)=+1` for `n=1 mod 6`, `-1` for `n=5 mod 6`. `P,S,C,Q` = `Omega = 1,2,3,>=4`. `R_i` = number of left endpoints `6k-1`, `k<=K`, in class `i` (= number of `n=5 mod 6`, `n<=6K+1`, with `Omega=i`); `C_j` the same for right endpoints `6k+1`. `pi_c, S_c, C_c` = counts of primes, semiprimes, 3-almost-primes `<=6K+1` in residue class `c mod 6`. `pi_chi(y)=sum_{5<=p<=y} chi(p)`, `pi'(y)=#{5<=p<=y}`.

Hypothesis H1 under test: the sign `N_PS>N_SP` is systematic, driven by (1) the exact class structure of semiprimes on each side and (2) Chebyshev's bias `pi_5>pi_1`.

## 1. Census (FINITE-EXACT)

`Omega`, `omega`, `b` were sieved exactly on `[1, 6*10^7+301]` (uint8, prime-power slicing), audited against independent trial division on 33,000 integers (all `n<=30000` plus 3,000 pseudo-random `n`, as in the script), and SW4 was re-checked for all `k<=10^7`.

Ordered matrices, rows `Omega(6k-1)`, columns `Omega(6k+1)`; the `K=10^6` matrix equals the inherited SW3 table entry by entry (the script raises otherwise):

| K=10^5 | 1 | 2 | 3 | >=4 |
|---|---:|---:|---:|---:|
| 1 | 5,330 | 10,091 | 6,537 | 2,614 |
| 2 | 10,050 | 17,794 | 10,719 | 3,851 |
| 3 | 6,561 | 10,751 | 5,600 | 1,664 |
| >=4 | 2,583 | 3,844 | 1,666 | 345 |

| K=10^7 | 1 | 2 | 3 | >=4 |
|---|---:|---:|---:|---:|
| 1 | 280,557 | 632,905 | 539,140 | 328,636 |
| 2 | 632,766 | 1,380,659 | 1,118,019 | 632,028 |
| 3 | 539,061 | 1,118,744 | 837,851 | 416,682 |
| >=4 | 328,491 | 631,688 | 416,902 | 165,871 |

Signed discrepancies `N_ij-N_ji`:

| K | PS-SP | PC-CP | SC-CS |
|---|---:|---:|---:|
| 10^5 | 41 | -24 | -32 |
| 10^6 | 412 | -286 | 111 |
| 10^7 | 139 | 79 | -725 |

Running sign for `K'<=10^7`, where "sign changes" counts transitions between `+` and `-` in the nonzero-sign subsequence of the running difference (zeros skipped): `N_PS-N_SP` is `>0` for 7,716,244 values of `K'`, `=0` for 14,947, `<0` for 2,268,809; 992 sign changes, first at `K'=1453`, last at 8,714,332; last `K'` with difference `<=0` is 8,714,401; running range `[-436, 650]`. `N_PC-N_CP`: 639 changes, first negative at `K'=21`, range `[-517, 570]`. `N_SC-N_CS`: negative for 9,114,559 values, 470 changes, last at 1,880,886, range `[-1323, 319]`. Among the ten tracked scales `10^4, 2*10^4, ..., 10^7` the difference `N_PS-N_SP` is negative at `K=2*10^5` (value `-26`), so "the class mechanism predicts the sign at all three scales" refers to `10^5, 10^6, 10^7` only.

**REFUTED (pointwise form of H1).** `N_PS(K)>=N_SP(K)` for all `K` fails first at `K=1453`: `N_PS=245`, `N_SP=246` (re-verified in the script by independent trial division of all endpoints up to `k=1453`). The strongest survivor is the drift statement of sections 4-5.

## 2. Exact class structure (PROVED)

**Lemma 2.1 (class-labelled shape determination).** For `gcd(n,6)=1` with `1<=Omega(n)<=3`, the class-labelled exponent shape of `n` (the multiset of pairs `(class mod 6, exponent)` over its prime factors) is determined by `(Omega, omega, b)`, and `n = (-1)^b mod 6`. For every `n` coprime to `6` the bare class multiset with multiplicity is `{1^(Omega-b), 5^b}`, determined by `(Omega, b)` alone (this and `n=(-1)^b` are inherited SW4). The determination of the labelled shape fails for `Omega=4`.

*Proof.* Every prime `p>=5` is `+-1 mod 6`, so `n = prod p_i = (-1)^b mod 6` and exactly `b` of the multiplicity-weighted factors are `5 mod 6`. For `Omega<=3`, `omega` fixes the exponent shape (`p`; `p^2`, `pq`; `p^3`, `p^2q`, `pqr`), and in each shape the exponents are either all equal or there is one distinguished factor of exponent `2` or `3` and the rest have exponent `1`; `b` then determines which class carries the distinguished factor, because the two classes contribute different multiplicities to `b` (for `p^2q`: `b in {0,1,2,3}` corresponds to `p1^2q1, p1^2q5, p5^2q1, p5^2q5` respectively). Hence the 2+5+10 = 17 keys with `Omega<=3` each carry one shape. For `Omega=4` the argument breaks because two different classes can carry the repeated prime while `b` stays fixed: `(4,3,2)` is realised both by `p5^2 q1 r1` and by `p1^2 q5 r5`, which have the same class multiset `{1,1,5,5}` but a different squared class; likewise `(4,2,0)` by `p1^3 q1` and `p1^2 q1^2`, and `(4,2,4)` by `p5^3 q5` and `p5^2 q5^2`. QED.

The script verifies this by full factorization of every `n<=300000` coprime to `6` with `Omega<=4`: the 17 keys with `Omega<=3` each carry exactly one shape; for `Omega=4` there are 17 keys of which 3 collide, with minimal witnesses `(4,3,2)`: `2275=5^2*7*13` vs `2695=5*7^2*11`; `(4,2,0)`: `4459=7^3*13` vs `8281=7^2*13^2`; `(4,2,4)`: `1375=5^3*11` vs `3025=5^2*11^2`.

Consequently (verbatim inherited SW4, not new) semiprimes on the right are `(1,1)`, `(5,5)` or `p^2`; on the left they are mixed `(1,5)`. What is new here is only the census. At `x=6*10^7+1`: `(1,1)` distinct 1,732,047; `(5,5)` distinct 2,030,969; mixed 3,763,472; `p^2` with `p=1 mod 6`: 484; `p=5 mod 6`: 496. Hence `S_1=3,763,996 > S_5=3,763,472`. The ten 3-almost-prime patterns at `x=6*10^7+1`: right (`b` even) `p1^3` 36, `p1^2q1` 84,521, `p5^2q1` 138,521, `p1q1r1` 522,639, `p1q5r5` 2,166,195 (total `C_1=2,911,912`); left (`b` odd) `p5^3` 39, `p1^2q5` 85,178, `p5^2q5` 139,167, `p1q1r5` 1,909,337, `p5q5r5` 778,617 (total `C_5=2,912,338`).

**Theorem 2.2 (semiprime class-bias identity).** For all `x>=1`,

```text
S_1(x) - S_5(x) = (1/2) [ sum_{5<=p<=x/5} chi(p) pi_chi(x/p) + pi'(sqrt x) ],
```

where `S_c(x)=#{n<=x: gcd(n,6)=1, Omega(n)=2, n=c mod 6}`; both sides vanish for `x<25`.

*Proof.* For `p,q>=5`, `chi(pq)=chi(p)chi(q)` and `chi(p^2)=1`, so `S_1-S_5 = sum_{p<q, pq<=x} chi(p)chi(q) + pi'(sqrt x)`. The ordered sum `O = sum_{p,q>=5, pq<=x} chi(p)chi(q) = sum_p chi(p) pi_chi(x/p)` (the range `p<=x/5` is exact since `pi_chi(x/p)=0` for `p>x/5`) equals `2 sum_{p<q} + pi'(sqrt x)`. Substituting gives the claim; for `x<25` there is no `p<=x/5` with `p>=5` and `pi'(sqrt x)=0`. QED.

Exact checks (brute force in the script): `x = 1, 24, 25, 26, 35, 49, 100, 1000, 10007` (at `x=10007` both sides equal `7`); and from the prime list, `x=600001`: `66=(-3+135)/2`; `x=6000001`: `244=(127+361)/2`; `x=60000001`: `524=(68+980)/2`.

**Attribution is split-dependent; only the total is canonical.** In the Ford-Sneed form `O/2 + pi'/2` the squares contribute `67.5, 180.5, 490.0` of `66, 244, 524`. With squares counted in full the squares contribute `pi'(sqrt x) = 135, 361, 980` and the distinct-pair semiprimes contribute `(O-pi')/2 = -69, -117, -456` (= same-class distinct minus mixed, checked against the pattern census): mixed pairs outnumber same-class distinct pairs at all three scales, and the class-1 surplus of semiprimes is entirely the squares. The cross term `O` is small and of variable sign (`-3, 127, 68`; partial sums through `p=5` alone `31, 44, 242`), so `S_1(x)>S_5(x)` for all large `x` is HEURISTIC only (OPEN).

**Lemma 2.3 (finite-modulus side symmetry).** For any pair family `(6k+c, 6k+c+g)` and any modulus `M` coprime to `6`, the map `k -> k* = -k-(2c+g)*6^{-1} (mod M)` (with `6^{-1}` the inverse of `6` modulo `M`; e.g. `c=1, g=6` gives the shift `13 mod 35`) is an involution of `Z/MZ` sending the endpoint pair `(n, n+g)` to `(-(n+g), -n) mod M`. Divisibility patterns by the primes of `M` are therefore exactly swapped between the two sides on any complete residue system, and no finite-modulus sieve datum distinguishes `N_ij` from `N_ji` at the level of residue-system distributions.

*Proof.* `6k*+c = -6k-(2c+g)+c = -(6k+c+g)` and `6k*+c+g = -(6k+c) mod M`; `k** = k`; negation preserves divisibility. This is the finite-modulus form of inherited SW2. QED. (Checked in the script through the nine-stratum CRT sizes of section 5(ii).)

## 3. Chebyshev-type biases by Omega (FINITE-EXACT; literature CITED)

`chi`-sums (`#class 1 - #class 5`) up to `x=6K+1`:

| K | Omega=1 | 2 | 3 | >=4 |
|---|---:|---:|---:|---:|
| 10^5 | -48 | +66 | -54 | +36 |
| 10^6 | -157 | +244 | -194 | +107 |
| 10^7 | -363 | +524 | -426 | +265 |

The prime race `pi_chi` is `<=0` at every one of the 3,562,131 primes `>=5` up to `60,000,301` (zero at exactly 9 prefixes, after `7, 13, 19, 37, 43, 79, 163, 223, 229`; never positive; `-363` at `6*10^7+1`). CITED (all four verified against publisher pages by the session's proof-audit lens): Rubinstein-Sarnak, *Chebyshev's bias*, Experimental Math. 3 (1994), no. 3: under GRH and the linear-independence hypothesis the logarithmic density of `{x: pi(x;3,2)>pi(x;3,1)}` is about `0.9990`. Bays-Hudson, *Details of the first region of integers x with pi_{3,2}(x) < pi_{3,1}(x)*, Math. Comp. 32 (1978), no. 142, 571-576: the first crossing is at `x = 608,981,813,029`. Ford-Sneed, *Chebyshev's bias for products of two primes*, Experimental Math. 19 (2010), no. 4, 385-398: for products of two primes the bias is reversed. Meng, *Chebyshev's bias for products of k primes*, Algebra & Number Theory 12 (2018), no. 2, 305-341 (arXiv:1606.04877): `Omega(n)=k` odd prefers the nonresidue class, `k` even the residue class. The sign pattern above matches this alternation strictly for `k=1,2,3`; the `>=4` column is a union over `k>=4`. The cited results are logarithmic-density statements under GRH+LI, not statements about every `x`; no asymptotic statement is imported.

## 4. Independence prediction and its exact splits (FINITE-EXACT algebra; HEURISTIC interpretation)

With `R=(R_i)`, `C=(C_j)` the observed marginals (proved equal to the class counts; checked), the product-of-marginals value of the antisymmetric part is `(R_1C_2-R_2C_1)/K = 2(Sbar*dpi + pibar*dS)/K`, `pibar=(R_1+C_1)/2`, `dpi=(R_1-C_1)/2`, `Sbar=(R_2+C_2)/2`, `dS=(C_2-R_2)/2` (exact algebra). By Theorem 2.2, `dS=(O+pi'(sqrt x))/4`, so the semiprime part splits further, in either convention:

| K | actual | pred | prime part `2Sbar dpi/K` | semi part | FS-square `pibar pi'/2K` | FS-cross `pibar O/2K` | full-square `pibar pi'/K` | distinct `pibar(O-pi')/2K` | act/pred |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10^5 | 41 | 36.58 | 20.37 (55.7%) | 16.20 | 16.57 (45.3%) | -0.37 | 33.14 (90.6%) | -16.94 | 1.121 |
| 10^6 | 412 | 113.23 | 62.87 (55.5%) | 50.37 | 37.26 (32.9%) | 13.11 | 74.52 (65.8%) | -24.15 | 3.639 |
| 10^7 | 139 | 229.95 | 136.62 (59.4%) | 93.33 | 87.27 (38.0%) | 6.06 | 174.54 (75.9%) | -81.22 | 0.604 |

So the deterministic fact that prime squares are `1 mod 6` carries between a third and (nearly) the whole of the predicted drift depending on the split; the prime-Chebyshev part is `55-59%` in both. The class mechanism predicts the sign at `10^5, 10^6, 10^7` (not at the tracked `K=2*10^5`), and the residuals (`4.42, 298.77, -90.95`) are as large as the prediction, so a single `K` cannot test H1. Ratio form: `N_PS/N_SP = 1.004080, 1.005263, 1.000220` against `(pi_5 S_1)/(pi_1 S_5) = 1.003516, 1.001371, 1.000343`; the user's "equal probability" survives only as a density statement (HEURISTIC). For `(P,C)` the `Omega=1` and `Omega=3` biases have the same sign and cancel in the prediction (`-1.47, 2.86, 29.84`); for `(S,C)` they reinforce (`-39.12, -144.37, -312.93`), predicting `N_CS>N_SC`.

HEURISTIC order of the drift. The identity's own terms carry `pibar/K ~ 3/log x` and `Sbar/K ~ 3 loglog x/log x` multiplying biases of size about `sqrt(x)/log x`, which suggests drift of order `sqrt(x) loglog x/log^2 x` rather than `sqrt(x)/log x`. Observed: `pred/(sqrt x/log x) = 0.63, 0.72, 0.53` (drifting down) and `pred/(sqrt x loglog x/log^2 x) = 3.23, 4.10, 3.30`. Three scales cannot settle this; it is recorded as the working hypothesis and it strengthens the fluctuation-dominated conclusion below.

## 5. Hostile controls (FINITE-EXACT)

**(i) Single shifted pairs.** Class-reversed `(6k+1,6k+5)`: `N_PS-N_SP = -5, -43, +461` (pred `-36.15, -112.83, -229.57`); reversed sign as predicted at `10^5, 10^6`, opposite at `10^7` (running range `[-441, 679]`, 839 sign changes). Same-class `(6k-1,6k+5)`: `+39, -170, +129` (running range `[-708, 336]`); `(6k+1,6k+7)`: `+109, -187, +108` (running range `[-436, 389]`); predictions `0.38-0.67`. These give the empirical single-pair noise floor: of the order of the prediction. The `+461` at `10^7` for the reversed pair is the strongest single datum against H1.

**(ii) The nine mod-35 strata.** Each of `5` and `7` divides the left endpoint, the right endpoint, or neither (never both, since `gcd(6k-1,6k+1)=1`), giving nine strata; per block of 35 centers their sizes are none 15, cross `1+1`, one-sided `5,5,3,3`, `35|left` 1, `35|right` 1 (sum 35; the script checks the nine cover every center exactly once). The draft's seven-stratum version omitted `35|left` and `35|right`; the seven sum to `42, 413, 140`, the nine to the true totals `41, 412, 139`, the missing `-1` being the SP pair `(35,37)` at `k=6`, the only `P` or `S` datum in either `35|` stratum. Census of `PS-SP` (`act` / `pred` from the stratum's own marginals):

| stratum | K=10^5 | K=10^6 | K=10^7 |
|---|---:|---:|---:|
| none (no 5, no 7) | 21 / 5.60 | 182 / 139.37 | -334 / 255.22 |
| cross (5 left, 7 right) | 0 / 0.03 | 0 / 0.03 | 0 / 0.02 |
| cross (7 left, 5 right) | 0 / 0.00 | 0 / 0.00 | 0 / 0.00 |
| 5 left only | -1656 / -1684.48 | -11205 / -11657.91 | -82573 / -85296.36 |
| 5 right only | 1657 / 1696.16 | 11422 / 11683.72 | 83023 / 85340.87 |
| 7 left only | -1080 / -1113.73 | -7434 / -7685.91 | -54349 / -56070.81 |
| 7 right only | 1100 / 1123.00 | 7448 / 7680.06 | 54373 / 56033.04 |
| 35 left | -1 / -0.36 | -1 / -0.30 | -1 / -0.26 |
| 35 right | 0 / 0.00 | 0 / 0.00 | 0 / 0.00 |

The cross strata contain no primes other than the PP pair `(5,7)` at `k=1` (the `(5|left,7|right)` stratum has `PP=1`, all other `P/S` cells zero). The `none` stratum carries all PP pairs but one (`5329` of `5330`, `37914` of `37915`, `280556` of `280557`) yet only `62.4-62.5%` of the left-endpoint primes (`1,113,235` of `1,781,238` at `10^7`); the one-sided strata carry `PS` or `SP` counts of order `8*10^4` each whose near-cancellation is part of the total. Lemma 2.3 is why no stratum can create the asymmetry from the sieve alone; the symmetric strata combined give `21, 182, -334` against `5.63, 139.40, 255.24`.

**(iii) Prime race.** See section 3.

**(iv) Pooled-gap control with block-level residual test.** For each orientation `[class n | class n+g]` take the 50 gaps `g` in `2..300` with both endpoints coprime to `6` (`[5|1]`: `n=6k-1`, `g=2 mod 6`; `[1|5]`: `n=6k+1`, `g=4 mod 6`; `[5|5]`, `[1|1]`: `g=0 mod 6`), `D_g=N_ij-N_ji` over `k<=K`, `pred_g` from the pair's own marginals. H1 predicts the same drift for every gap of an orientation (it depends only on marginals). Split `K` into 100 blocks (size `10^4` at `K=10^6`, `10^5` at `K=10^7`), with block prediction `p_(g,b)` and block residual `r_(g,b)=D_(g,b)-p_(g,b)`; report the pooled mean of `D_g`, its across-gap sd (ddof 1), the number of positive `D_g`, the mean prediction, the correlation and slope of the block means `m_b=mean_g D_(g,b)` against `p_b=mean_g p_(g,b)`, the mean inter-gap correlation of `D_(g,b)` and of `r_(g,b)` over blocks, and the residual test `RES = sum_b mean_g r_(g,b)`, `SE = sd_b(mean_g r_(g,b)) sqrt(100)`, `z=RES/SE`. The block counts at the four control gaps reproduce the matrices of (i) (checked).

`PS-SP`:

| K | orient | mean D_g | sd | #D_g>0 | mean pred | corr(m_b,p_b) | slope | inter-gap corr D / r | RES | SE | z |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| 10^6 | [5\|1] | +100.7 | 197.4 | 32/50 | +117.5 | 1.00 | 1.00 | 0.50 / -0.01 | -26.1 | 21.1 | -1.24 |
| 10^6 | [1\|5] | -113.0 | 168.1 | 14/50 | -109.1 | 1.00 | 0.99 | 0.48 / -0.01 | 4.2 | 22.9 | 0.18 |
| 10^6 | [5\|5] | +9.7 | 215.1 | 25/50 | +4.1 | 0.48 | 0.94 | -0.01 / -0.01 | 4.8 | 22.5 | 0.21 |
| 10^6 | [1\|1] | +12.5 | 174.1 | 23/50 | +4.4 | 0.61 | 1.14 | -0.01 / -0.01 | 7.5 | 20.3 | 0.37 |
| 10^7 | [5\|1] | +195.1 | 541.2 | 30/50 | +234.2 | 0.99 | 1.03 | 0.39 / -0.01 | -59.2 | 64.9 | -0.91 |
| 10^7 | [1\|5] | -66.4 | 493.3 | 25/50 | -225.6 | 0.99 | 0.99 | 0.37 / -0.01 | 177.8 | 60.9 | 2.92 |
| 10^7 | [5\|5] | +96.7 | 585.9 | 28/50 | +4.3 | 0.17 | 1.21 | -0.01 / -0.01 | 91.6 | 85.1 | 1.08 |
| 10^7 | [1\|1] | -32.7 | 500.6 | 24/50 | +4.4 | 0.14 | 0.81 | -0.01 / -0.01 | -37.7 | 79.3 | -0.48 |

`SC-CS` pooled means (pred): `10^6`: `[5|1] -168.7 (-141.6)`, `[1|5] +171.5 (+146.2)`, `[5|5] +19.1 (+1.8)`, `[1|1] -39.1 (+2.9)`; `10^7`: `-369.0 (-310.6)`, `+498.2 (+316.2)`, `+116.8 (+3.2)`, `-66.4 (+2.5)`; residual `z`: `-0.98, 0.78, 0.48, -1.35` and `-0.59, 2.00, 1.03, -0.63`. `PC-CP` pooled means (pred): `10^6`: `27.1 (7.2)`, `-15.1 (0.9)`, `-0.6 (3.7)`, `19.2 (4.5)`, all `|z|<=1.05`; `10^7`: `74.8 (34.2)`, `-250.1 (-25.0)`, `-75.7 (4.9)`, `-1.8 (4.6)`, with `z=-3.16` for `[1|5]`.

What the control shows, in order of strength. (a) Block by block, the gap-pooled discrepancy is the block's own product of marginals: for the opposite-class orientations of `PS-SP` and `SC-CS`, `corr(m_b,p_b)=0.99-1.00` with slope `0.98-1.03` (for `PC-CP`, whose predictions are near zero, `0.93-0.94` with slope `0.98-1.05`), and the unconditional inter-gap correlation (`0.37-0.63` for opposite-class `PS-SP` and `SC-CS`, `0.07-0.09` for `PC-CP`, `about 0` for same-class orientations) is entirely this shared local class-bias fluctuation, since the residual inter-gap correlation is `-0.01` or `-0.00` in every cell. Hence the naive across-gap SE is valid for comparing the pooled mean with the local marginal prediction, but not against zero; and the agreement of a 50-shift average with a product of marginals is close to tautological (a 50-shift average reproduces the local product of marginals unless there is genuine pairing correlation). The control is therefore a consistency check that asymmetric short-range coupling, averaged over gaps `2..300`, is absent, which Lemma 2.3 already predicts at sieve level; it is not decisive and says nothing about `g=2` specifically. (b) Residual test: at `10^6` all twelve cells have `|z|<=1.35` (block SE `20-35`); at `10^7` the three start-class-1 families deviate from the local-marginal prediction: `PS-SP [1|5] +177.8 (SE 60.9, z=2.92)`, `PC-CP [1|5] -240.5 (SE 76.0, z=-3.16)`, `SC-CS [1|5] +180.6 (SE 90.5, z=2.00)`; all other cells `|z|<=1.08`. Three simultaneous `2-3 sigma` deviations, all for pairs starting at `6k+1`, are left OPEN (correlated-block artifact or genuine coupling beyond the marginals). (c) Pooled means versus prediction: at `10^6` the orientation signature is present with the predicted signs and magnitudes (`PS-SP` within `14%` and `4%`, `SC-CS` within `19%` and `17%`), and the same-class means are small. At `10^7` it is present for `[5|1]` (`17%` off) and for `SC-CS [5|1]` (`19%`), but `PS-SP [1|5]` is `71%` off and `SC-CS [1|5]` `58%` off, and the same-class `[5|5]` pooled mean `+96.7` exceeds the opposite-class `[1|5]` mean `-66.4` in magnitude: the draft's falsifier condition ("a same-class pooled mean as large as an opposite-class one") is met at `10^7` for that pair of orientations. (d) The sum of block predictions exceeds the global prediction at `g=2` (`122.02` vs `113.23` at `10^6`; `249.36` vs `229.95` at `10^7`): local class biases of primes and semiprimes co-vary positively across blocks. The `g=2` residual against the block predictions is `+290.0` (`1.47` per-gap residual sd) at `10^6` and `-110.4` (`-0.20` sd) at `10^7`; `4` and `25` of the other 49 `[5|1]` gaps have larger `D_g`, so the sandwich gap shows no special pairing correlation. (e) Noise versus drift: the across-gap sd (`168-215` at `10^6`, `493-586` at `10^7` for `PS-SP`) exceeds the predicted drift; as quotients of the table's numbers, `117.5/197.4 = 0.60` and `234.2/541.2 = 0.43` for `[5|1]`, decreasing with scale.

## 6. Verdict on H1 (HEURISTIC synthesis)

H1's two structural inputs are exact (Lemma 2.1, Theorem 2.2, the census of section 3), and a third deterministic input, prime squares being `1 mod 6`, carries between a third and most of the predicted drift depending on the split convention. The class mechanism produces a positive expected drift of `N_PS-N_SP` whose block-level form (local discrepancy = local product of class biases, corr about 1) is confirmed, and the pooled residual test is clean at `10^6`. At `10^7` the fit is partial: the `[5|1]` and `SC-CS` orientations fit, the three `[1|5]` families deviate at `z=2.0-3.2` from the local-marginal prediction (OPEN), and the draft's falsifier condition fires for `[5|5]` against `[1|5]`. "The sign is systematic" is false pointwise (`K=1453`, also `K=2*10^5` among the tracked scales), and the single-`K` sign is fluctuation-dominated at every scale examined: H1 is confirmed as a drift statement at `10^6` and partially at `10^7`, and refuted as a sign law.

## Reproduction

```bash
python3 04-computation/experiments/collatz_mod6_20260917_sandwich_bias.py > 05-knowledge/results/collatz_mod6_20260917_sandwich_bias.out
```

About 30 s and max RSS 0.87 GB on the finalizing mac-mini (the `[t=.. maxrss=..]` stamps are printed in the output); `python3 -O` gives identical output modulo those stamps; every check is an explicit `raise`, and the script raises on any disagreement with the inherited SW3 table, on the wrong minimal witness, or if the nine strata fail to partition the centers. SHA-256 at finalization: script `f65540f1e73ac0b4a97b4248d252b454a19aaa98462ea38b384143ffca2df778`, output `bc45bc1c3434c77d60861d8515078b7a69cdab9092b829bc1776fd36daf3f290` (the output hash changes with the timing stamps). The recovered 2026-09-17 draft lacked section 4(iv), the three-way split, and the two `35|` strata, and peaked at 1.7 GB when re-run during finalization; the two independent audit recomputations of 2026-09-17 (spf-sieve based) reproduced every number of the draft and identified the corrections applied here.

## Stopping boundary / next question

Repeating the census at larger `K` will not change the verdict: the class drift and the pairing fluctuation both grow, and the sign at a single `K` stays fluctuation-dominated. Three questions are substantive. (1) The variance: is the pairing fluctuation of `N_PS(K)-N_SP(K)` genuinely of order `sqrt(K)`, as the across-gap sd suggests, while the drift is of order `sqrt(x) loglog x/log^2 x`? A conditional second-moment analysis of the pair correlation, or a disjoint-block paired sign test with an honest standard error, would decide whether the class drift can ever dominate. (2) The start-class-1 anomaly at `10^7`: are the three `2-3 sigma` residual deviations for pairs starting at `6k+1` a block-correlation artifact, or a coupling between prime and 3-almost-prime endpoints beyond the marginals? A 500-block or disjoint-scale replication is the cheapest test. (3) Prove `S_1(x)>S_5(x)` for all large `x` from Theorem 2.2 by bounding the cross term `O(x)`, and extend the exact three-way split to the `(S,C)` and `(P,C)` pairs via an analogous identity for `sum_{Omega=3} chi(n)`.
