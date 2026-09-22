# Scaffolding audit: every peripheral claim of the geometric Collatz blueprint, checked exactly

**Status:** PROVED (scoped, elementary): the corrected core local identity and the unique solution of the two flow equations (section 2); the d-dimensional simplex convolution law containing the tetrahedral and pentatope laws (section 3); the unique symmetric reading of the g-operator (section 4); the double-counting identity `sum_T h(T) = n! 2^{T(n-2)}` and the resulting cost bound for the proposed tournament code (section 5); the residue law "every shortcut peak of an odd seed is 2 mod 6", the ladder `C^j(a 2^k - 1) = a 3^j 2^{k-j} - 1`, the exclusion "no shortcut peak is 6 mod 8" (Lemma 8.3) and the multiplicity lemma "the odd-preimage chain of a peak is its ladder, of length `v_3(p+1)`; the seed is unique iff `9 ∤ p+1`" (Lemma 8.4) (section 8). FINITE-EXACT: all grids, censuses and tables below (`[0,11]^2` simplex checks through `d=8`; labelled tournaments `n<=6`; Bang/Zsigmondy tables `n<=40`; Cipolla through `p=23`; 49,999 odd seeds `<=10^5`, their 20,341 distinct peaks, and the residue-stratified null of section 8 under which both square windows are at chance). CITED: Bang 1886 / Zsigmondy 1892, Mihailescu 2004 (Catalan), Cipolla 1904, Redei 1934, Shannon 1948, Poonen 1998, Bott 1959. REFUTED (minimal witnesses in the table): the stated core identity, the g-operator readings, "path edges carry 0 bits", the 341 bottleneck, the peak endpoints 1/26/80, the 4k+1 square tracking, "196 = sum of the first 12 primes", and "growth bound by the area of squares". SCOPE (no map found): every "physical insight", the horizons, the unified cycle, Bott/tetration, the quasicrystalline Fourier transform, the Wythoff diffraction. OPEN: exactly what the [blueprint audit](collatz_blueprint_20260921_synthesis.md) left open, a positive-margin word for every fixed `n>1`. No novelty claim for any of the surviving identities (Chu-Vandermonde, Redei, Bang, Cipolla, the inherited descent certificate; section 6 is entirely inherited from the wave-one [wild_typing lane](collatz_mod6_20260917_wild_typing.md)); the two small exact observations are Lemmas 8.3-8.4 (elementary; not found in the inherited notes; whether they are in the literature is UNCITED). The first version's "ladder explanation of the square-minus-one peaks" was a tautology and is withdrawn: adversarially audited 2026-09-21 ([audit script](../../04-computation/experiments/collatz_mod6_20260917_scaffolding_audit_audit.py), [audit output](collatz_mod6_20260917_scaffolding_audit_audit.out)). Results note of session collatz-mod6-20260917 (machine mac-mini), lane `scaffolding_audit`, written 2026-09-21; not a reserved canon ID.

## Inheritance and concept board

The pasted block is the scaffolding of the geometric blueprint whose proof claim was REFUTED in the [blueprint synthesis](collatz_blueprint_20260921_synthesis.md), with the exact affine/carry model in its [affine companion](collatz_blueprint_20260921_affine.md) and the [source document](../reference/COLLATZ-BLUEPRINT-2026-09-21-SOURCE.md). Everything Collatz-arithmetic here is inherited and cited, not re-derived: the inverse-fibre braid `R(n)=4n+1`, `(B1)-(B3)` of the [first Collatz note](arithmetic_braids_20260917_collatz.md); the descent certificate `2^K T^L(n) = 3^L n + B` of the blueprint synthesis, section 6; the Catalan reading of Bang's exception `2^6-1=7*9` from the wave-one [zsigmondy lane](collatz_mod6_20260917_zsigmondy_triad.out); Cipolla's trunk pseudoprimes, the exact criterion `6j | 4^j-4`, Zsigmondy without exception for `4^n-1` and the refutation of the 341 "bottleneck" from the wave-one [wild_typing lane](collatz_mod6_20260917_wild_typing.md), section 3; the rational 3-cycle of `x^2-29/16` from [THM-4139](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md) and [THM-4146](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md); the tournament `h`-spectrum holes `{7,21}` and `|Aut| | h` from [THM-1745](../../01-canon/theorems/THM-1745-leaf-graded-arborescence-filtration-721-shadow.md). **Closest proved mechanism:** the descent certificate itself, which the paste restates verbatim as `3^L n + B_L < 2^{K_L} n` and then asserts without argument. **Canonical hostile:** the ladder seed `2^k-1`, which climbs to `3^k-1` (section 8); it refutes "peaks bound by area" and, at the same time, is the whole explanation of the only real signal in the paste's square-tracking data. **Corrected near miss:** two. The core local identity is off by exactly one in both triangular arguments; with both shifted it is a true polynomial identity and the paste's own flow equations are its partial differences. And this note's own first version "explained" the 93 square-minus-one peaks by exhibiting, for each, a ladder seed `a 2^k - 1` whose peak is `a 3^k - 1 = p`; the audit found that this holds for every one of the 20,341 distinct peaks (`k=1` is the peak's own predecessor `(2p-1)/3`), so it was a tautology. The true comparison is the residue-stratified null of section 8, under which the square windows are at chance. **Least-used sidecar:** the peak rate by residue class mod 24 (`0.688, 0.889, 0, 0.685` on `2, 8, 14, 20 mod 24` below the seed box), which is what any window statistic on peaks has to be compared with; the paste's windows `s^2-1` and `s^2+1` sit in the classes `0 mod 8` and `2 mod 8`. **Typed analogy (ladder seeds -> odd-preimage chain):** source = the forward ladder `a 2^k - 1 -> a 3^k - 1` of Lemma 8.2; target = the odd-preimage chain of a peak `p` in the inverse tree; map = `k -> m^(k) = ((p+1)/3^k) 2^k - 1`; preserved = the values and the monotone climb to `p`; lost = nothing for `k <= v_3(p+1)`, and the map does not exist beyond it; sidecar = `v_3(p+1)` = chain length = ladder height; test = multiplicity 1 iff `9 ∤ p+1`, which passes on all 20,341 distinct peaks (Lemma 8.4). **No map found (SCOPE):** from "shear as swept volume", "Bott periodicity", "quasicrystalline Fourier transform", "Wythoff diffraction", "unified arithmetic cycle", "inverse square horizon" to any defined mathematical object.

Notation. `T(n)=n(n+1)/2`, `Te(n)=C(n+2,3)`, `Pt(n)=C(n+3,4)`, and generally `S_d(n)=C(n+d-1,d)` with `S_0=1`, `S_1(n)=n`, `S_2=T`. Shortcut Collatz map `C(n)=n/2` (even), `(3n+1)/2` (odd); full map `n/2`, `3n+1`. Accelerated Syracuse `T(n)=(3n+1)/2^k` on odd `n`. The peak of a seed is the maximum of its orbit down to 1.

## 1. Verdict table

| claim (paste numbering) | verdict | witness or proof pointer |
|---|---|---|
| 1 core identity `T((A+1)(B+1)-1)-T(AB-1)=T(A)+T(B)+Shear` | FALSE / TRUE-with-repair | fails on all 49 pairs of `[1,7]^2`; `A=B=1`: 6 vs 2; true after shifting both arguments by `-1` (section 2) |
| 1 horizontal and vertical flows of `S(A,B)` | TRUE-with-repair | compatible; unique solution `S=AB(A+B)+C(A,2)+C(B,2)` = the paste's right side, not its left side (section 2) |
| 1 "shear = swept volume", "conservation law", "multiplication = integrated accumulation" | SCOPE | the true shear is `AB(A+B+2)=2A T(B)+2B T(A)`; no volume or rotation is defined |
| 2 tetrahedral law | TRUE | Chu-Vandermonde, split `(2,1)` (section 3) |
| 2 pentatope law | TRUE | Chu-Vandermonde, split `(2,2)` (section 3) |
| 2 "prisms locking against a plate", "multiplication of two triangular surfaces" | SCOPE | a two-block set-partition count |
| 3 `AgB = AB*(AgB)` | FALSE | forces `(AB-1)(AgB)=0` |
| 3 `AgB=(AB)^{AB}` together with `1gB=B!` | FALSE | `B=2`: 4 vs 2 |
| 3 `AgA=(A!)^2` | TRUE-with-repair | holds for the unique symmetric reading `AgB=A!B!` (section 4) |
| 3 "squares the information capacity; Pell and square-triangular families" | SCOPE | no map found |
| 4 free edges `= T_{n-2}` | TRUE | `C(n,2)-(n-1)=C(n-1,2)` |
| 4 "Hamiltonian path edges carry 0 bits" | FALSE | code costs `C(n,2)+log2(n!/2^{n-1})` bits, excess `0.585` at `n=3` (section 5) |
| 4 "Wythoff diffraction peaks", "topological genus of randomness" | SCOPE | no object |
| 5 horizon lines `(Z,2Z),(Z,Z^2),(W,W^3),(W,W^4)` | SCOPE | definitions of diagonals |
| 5 63 bottleneck "short-circuits Bang" | TRUE-with-repair | `63=2^6-1` is Bang's unique exception `n>1`, via `2^3+1=3^2`; no "loop" is defined (section 6) |
| 5 341 bottleneck "stalls primitive primes" | FALSE | `4^5-1=3*11*31` has two primitive primes 11, 31 (section 6) |
| 5 341 `= 4^4+4^3+4^2+4+1 = 11*31` | TRUE | arithmetic; 341 is `R^4(1)` on the trunk of 1, a Cipolla pseudoprime |
| 5 rational 3-cycle at `x_0=-7/4`, `Q=-29/16` | TRUE (content = THM-4139/4146) | `-7/4 -> 5/4 -> -1/4 -> -7/4`; "horizon cancels translation": SCOPE (section 7) |
| 6 unified cycle `{f,+,*,g}`, Bott/tetration, quasicrystal Fourier of mod-30 primes | SCOPE | only the mod-30 wheel `{1,7,11,13,17,19,23,29}` is a defined object |
| II peaks of 5, 7, 23 are 1, 26, 80 | FALSE | shortcut peaks 8, 26, 80; full-map peaks 16, 52, 160 (section 8) |
| II peaks track squares of `4k+1` | FALSE | bases 1 mod 4 and 3 mod 4 equally represented (46/47 and 49/50 distinct peaks); both windows at chance under the residue-stratified null (section 8) |
| II 196 = sum of the first 12 primes | FALSE | the sum is 197 |
| II 160 `= 2^5*5` on the trunk of 5 | TRUE | 160 is the full-map peak of 15, 23, 35, 53; it is unrelated to 169 |
| II `3^L n + B_L < 2^{K_L} n` | TRUE (inherited) | the descent certificate; the obligation is OPEN (section 9) |
| II "peaks bound by the area of squares" proves descent | FALSE | `2^k-1` climbs to `3^k-1`; no proof supplied (sections 8, 9) |

## 2. Claim 1: the core local identity and the flows (PROVED)

**Lemma 2.1.** For all integers `A,B`,

```text
T((A+1)(B+1)-1) - T(AB-1) = sum_{k=AB}^{AB+A+B} k = (A+B+1)(2AB+A+B)/2
                          = T(A) + T(B) + AB(A+B+2) = T(A) + T(B) + 2A T(B) + 2B T(A).
```

*Proof.* `T(N)-T(M-1)` is the sum of the integers from `M` to `N`; with `M=AB` and `N=AB+A+B` there are `A+B+1` terms of mean `(2AB+A+B)/2`. Expanding, `(A+B+1)(2AB+A+B)/2 = A^2B+AB^2+2AB+(A^2+A)/2+(B^2+B)/2`, and `AB(A+B+2)=AB(B+1)+AB(A+1)=2A T(B)+2B T(A)`. The script checks these as exact polynomial identities in `Q[A,B]`. QED.

The paste's shear is `A(B^2-1)+B(A^2-1)=AB(A+B)-A-B`, so

```text
LHS - [T(A)+T(B)+A(B^2-1)+B(A^2-1)] = 2AB + A + B = A(B+1) + B(A+1),
```

never zero for positive `A,B`: the stated identity fails on all 49 pairs of `[1,7]^2` (witness `A=B=1`: `T(3)-T(0)=6` against `T(1)+T(1)+0=2`).

**Lemma 2.2 (no quadratic reading rescues the stated form).** For any `T(n)=a n^2+b n+c`, the stated form forces `a=1/2` from the `A^2B` coefficient, after which the left side keeps a `2AB` term that no right-side term produces. So no quadratic reinterpretation of `T` (and in particular no shift of its argument) makes the stated form an identity.

**Lemma 2.3 (the repair).** As polynomials,

```text
T((A+1)(B+1)-2) - T(AB) = T(A) + T(B) + A(B^2-1) + B(A^2-1).
```

The exhaustive search over `T(N+s_1)-T(AB+s_2)=T(A+t_1)+T(B+t_2)+Shear` with all shifts in `[-3,3]`, on `[0,7]^2`, returns exactly two solutions, `(s_1,s_2,t_1,t_2)=(-1,0,0,0)` and `(-3,-2,-2,-2)`; the second is the reflection of the first under `T(n-1)=T(-n)`. So the paste is off by one in both triangular arguments (witness of the repair: `A=B=1`, `T(2)-T(1)=2=T(1)+T(1)+0`).

**Lemma 2.4 (the flows).** The increments `g_h(A,B)=A+B^2+2AB+B` and `g_v(A,B)=B+A^2+2AB+A` are compatible (both mixed differences equal `2A+2B+2`), and the solutions of `S(A+1,B)-S(A,B)=g_h`, `S(A,B+1)-S(A,B)=g_v` on `Z_{>=0}^2` are exactly

```text
S(A,B) = AB(A+B) + C(A,2) + C(B,2) + S(0,0) = T(A)+T(B)+A(B^2-1)+B(A^2-1) + S(0,0) = T((A+1)(B+1)-2) - T(AB) + S(0,0).
```

*Proof.* `Delta_A[A^2B+AB^2]=2AB+B+B^2`, so the residual horizontal increment is `A=Delta_A C(A,2)`; symmetrically for `B`. Two solutions of both difference equations differ by a function constant in each variable, hence by a constant. The equality with the paste's right side is `T(A-1)+A=T(A)`. QED.

So the paste's flows are the partial differences of its right side, and the paste's left side is not their `S`: its horizontal flow is `g_h+(2B+1)`. Verdict TRUE-with-repair. The "physical insight" (a swept volume, a conservation law for discrete area, "multiplication is an integrated accumulation of line components") is SCOPE: nothing is defined beyond `AB=sum_{k=1}^{B} A`.

## 3. Claim 2: the simplex laws are Chu-Vandermonde (PROVED)

Both stated laws hold on `[0,11]^2` with 0 failures, as does the 2D member `T(A+B+1)=T(A)+T(B)+(A+1)(B+1)` (the correct "locking" identity that Claim 1 garbles).

**Theorem 3.1 (general d-simplex law).** For every `d>=1` and every split `a+b=d` with `a,b>=0`,

```text
S_d(A+B+1) = C(A+B+d, d) = sum_{i=0}^{d} C(A+a, i) C(B+b, d-i).
```

*Proof.* A `d`-subset of a set of `A+B+d=(A+a)+(B+b)` elements is chosen by how many elements it takes from the first block. QED. Checked exactly for all `d<=8`, all splits, on `[0,11]^2`.

The tetrahedral law is the split `(a,b)=(2,1)` after Pascal: `C(B+1,3)+(A+2)C(B+1,2)=Te(B)+(A+1)T(B)`, `C(A+2,2)(B+1)=(B+1)T(A)+(A+1)(B+1)`, `C(A+2,3)=Te(A)`. The pentatope law is the split `(2,2)`: `C(B+2,4)+(A+2)C(B+2,3)=Pt(B)+(A+1)Te(B)` (and symmetrically), with middle term `C(A+2,2)C(B+2,2)=T(A+1)T(B+1)`. A symmetric-looking companion, `S_d(A+B+1)=sum_{i+j=d} S_i(A+1) S_j(B)`, is the identity `sum_k C(x+k,k)C(y+d-k,d-k)=C(x+y+d+1,d)` at `x=A`, `y=B-1`, checked for `d<=8`.

What survives of the "physical insight": a convolution identity of two blocks. There is no shear, rotation, or plate; every term is `C(x,i)C(y,d-i)`.

## 4. Claim 3: the g-operator (REFUTED, with the consistent reading)

`AgB=AB*(AgB)` forces `(AB-1)(AgB)=0`, so it defines nothing unless `AB=1`. The two displayed values are mutually inconsistent: `(1*B)^{1*B}=B^B` while `1gB=B!`, and `B^B != B!` for `B>=2` (minimal witness `B=2`: 4 vs 2; then 27 vs 6, 256 vs 24).

**Lemma 4.1.** Any `f` with `f(A,B)=AB f(A-1,B-1)` (for `A,B>=1`) and `f(1,B)=B!` has `f(0,k)=k!` and `f(A,A)=(A!)^2`; the unique symmetric solution is `f(A,B)=A!B!`, which satisfies all three of the paste's statements at once. Symmetry is load-bearing: `f(k,0)` is free for `k>=2`, and `f(A,B)=A!B! 5^{A-B-1}` for `A>B` (else `A!B!`) also satisfies both conditions (checked in the audit). The alternative recursion `AgB=AB*(Ag(B-1))` gives `A^B B!`, which has `1gB=B!` but `AgA=A^A A!` (`A=2`: 8 vs 4); `(AB)^{AB}` satisfies neither recursion.

"Shifting from asymmetry to symmetry squares the information-carrying capacity, mapping the 2D Pell and square-triangular families": SCOPE, no map found. `(A!)^2` is a square by construction; square-triangular numbers are the Pell equation `(2m+1)^2-8k^2=1` and involve no factorials.

## 5. Claim 4: tournament compression (identity TRUE, "0 bits" REFUTED)

`C(n,2)-(n-1)=C(n-1,2)=T(n-2)` for all `n>=2` (checked `n<=12`). Exact census of all labelled tournaments:

| n | tournaments `2^{C(n,2)}` | `sum_T h(T)` | `n! 2^{T(n-2)}` | max h | h-values seen |
|---:|---:|---:|---:|---:|---|
| 2 | 2 | 2 | 2 | 1 | 1 |
| 3 | 8 | 12 | 12 | 3 | 1, 3 |
| 4 | 64 | 192 | 192 | 5 | 1, 3, 5 |
| 5 | 1,024 | 7,680 | 7,680 | 15 | 1, 3, 5, 9, 11, 13, 15 |
| 6 | 32,768 | 737,280 | 737,280 | 45 | 1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45 |

Every `h` is odd (Redei 1934) and 7, 21 never occur (THM-1745's holes, re-seen here for `n<=6`; at `n=6` the odd values 35 and 39 are also absent, which is not a hole, since the spectrum statement of THM-1745 is over all `n`).

**Theorem 5.1 (double counting).** Every ordering `pi` of `[n]` is a Hamiltonian path of exactly `2^{T(n-2)}` labelled tournaments (its `n-1` edges are forced, the other `T(n-2)` are free). Hence `sum_T h(T)=n! 2^{T(n-2)}` and the mean number of Hamiltonian paths is `n!/2^{n-1}` (`n=6`: `45/2`). QED.

**Corollary 5.2 (the proposed code does not compress).** The code "(a Hamiltonian path order, `T(n-2)` free bits)" has `n! 2^{T(n-2)}` codewords for `2^{C(n,2)}` tournaments; it is surjective and `n!/2^{n-1}`-to-one on average, so it costs `log2(n!)+T(n-2)=C(n,2)+log2(n!/2^{n-1})` bits: excess 0.585 (`n=3`), 1.585 (`n=4`), 2.907 (`n=5`), 4.492 (`n=6`), 8.299 (`n=8`), 12.791 (`n=10`). The path edges carry `log2(n!)>=n-1` bits, not 0. The uniform distribution on labelled tournaments has entropy exactly `C(n,2)` bits, and no lossless code beats the entropy (Shannon 1948). "Defects at the diffraction peaks of the quasicrystalline Wythoff wave" and "file sizes scale with the topological genus of randomness" are SCOPE: no object is defined.

## 6. Claim 5: the 63 and 341 "bottlenecks"

**63 (TRUE-with-repair).** `63=3*7*3=3^2*7=2^6-1=(2^3-1)(2^3+1)=7*9`. Among `2^n-1`, `1<=n<=40`, the indices with no primitive prime divisor are exactly `[1,6]`; by Bang/Zsigmondy `n=6` is the only exception with `n>1`. The mechanism is `2^3+1=9=3^2` (Catalan `3^2-2^3=1`, Mihailescu 2004), so the would-be primitive part `Phi_6(2)=3` already divides `2^2-1`; equivalently no prime has `ord_p(2)=6`, while `ord_9(2)=6`. This is the wave-one zsigmondy lane's Catalan reading, cited from its [regenerated output](collatz_mod6_20260917_zsigmondy_triad.out). The paste's "closed loop switching between a multiplication edge and an addition edge (3 x 7 x 3), recycling old prime factors" defines no loop; only the factorisation is true.

**341 (REFUTED).** `341=11*31=(4^5-1)/3=4^4+4^3+4^2+4+1`. `4^5-1=3*11*31` and both 11 and 31 have `ord_p(4)=5`: two primitive primes, so 341 does not "stall primitive prime generation". Zsigmondy's exceptions are `n=1` with `a-b=1`, `n=2` with `a+b` a power of 2, and `(2,1,6)`; since `4-1=3` and `4+1=5`, every `4^n-1` has a primitive prime (table for `n<=40` in the output; none missing). Everything in this paragraph, including the refutation, is already in the wave-one [wild_typing lane](collatz_mod6_20260917_wild_typing.md), section 3 (Cipolla with proof, the exact criterion `6j | 4^j-4`, Zsigmondy without exception for `4^n-1`, the tower closure); this section only re-verifies it. What is true about 341 is Cipolla (1904): for prime `p>=5`, `N_p=(4^p-1)/3` is a composite base-2 pseudoprime, `2^{N_p-1}=1 mod N_p`; verified for `p=5,7,11,13,17,19,23`, giving 341, 5461, 1398101, 22369621, 5726623061, 91625968981, 23456248059221 with their factorisations (proof sketch in the output: `N_p-1=4(4^{p-1}-1)/3` is divisible by `2p`, and `2^{2p}=3N_p+1=1 mod N_p`). Under `R(n)=4n+1` the orbit of 0 is `0,1,5,21,85,341,1365=(4^j-1)/3`, and `3(4^j-1)/3+1=4^j`: these are the odd predecessors of the powers of two, the trunk of 1 in the inherited inverse-fibre braid `(B1)-(B2)`. "A zero-friction state in `G_4`" is SCOPE.

**Horizons.** `(Z,2Z)`, `(Z,Z^2)`, `(W,W^3)`, `(W,W^4)` are the diagonals `X=Y` of the sum and product graphs: definitions, SCOPE.

## 7. Claim 5c: the rational 3-cycle (SCOPE beyond THM-4139/4146)

`x -> x^2-29/16` sends `-7/4 -> 5/4 -> -1/4 -> -7/4`; the three points form an arithmetic progression of difference `3/2`. What is PROVED in the repo and is not re-derived here: [THM-4139](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md) (the complete rational preperiodic graph is eight quarter-integers plus this 3-cycle; no rational 6-cycle; uniqueness among centred monic quadratics over `Q` with a nondegenerate AP-supported 3-cycle; a determinant-one lift `B` with `B^3=-I`) and [THM-4146](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md) (trace-one `SL_2` lift, integral hexagon on `X^2+2XY+13Y^2=48`, the signed Pythagorean template forcing the 3:4:5 class and `29=5^2+2^2`); the general classification of rational 3-cycles of `x^2+c` is Poonen 1998. "The vector of the inverse square horizon cancels the spatial translation after exactly three steps" is not a statement about this map: SCOPE, no map found.

## 8. Second part: Collatz peaks and squares (REFUTED, with the exact structure behind its data)

**Peaks of 5, 7, 23.** Shortcut orbits: `5, 8, 4, 2, 1`; `7, 11, 17, 26, 13, 20, 10, 5, 8, 4, 2, 1`; `23, 35, 53, 80, 40, 20, 10, 5, 8, 4, 2, 1`. Peaks 8, 26, 80 (full map: 16, 52, 160). "Sink 1 = 1^2" is not a peak. Incidentally `8=3^2-1`, `26=5^2+1`, `80=9^2-1` with bases 3, 5, 9.

**Lemma 8.1 (peak residue).** For an odd seed `n>1` the shortcut peak `p` is even (an odd `p>1` is followed by `(3p+1)/2>p`), so `p=(3m+1)/2` with `m` odd and `3m+1=0 mod 4`, i.e. `m=1 mod 4` and `p=2 mod 6`; under the full map the peak is `3m+1=4 mod 6`. Consequently a peak in the `s^2-1` window has `3|s`, and one in the `s^2+1` window has `3 ∤ s`; the middle window `s^2` is never hit (`s^2` is never `2 mod 6`). FINITE-EXACT: all 49,999 shortcut peaks of odd seeds in `[3,10^5]` are 2 mod 6 and all full-map peaks are 4 mod 6.

**Lemma 8.2 (ladder; the hostile).** `C^j(a 2^k-1)=a 3^j 2^{k-j}-1` for `j<=k` and odd `a`. So `2^k-1` climbs to `3^k-1` (checked `k=5,10,20`): peak/seed is unbounded, and the largest ratio among the census seeds is `10112.041` at seed 77671 (peak 785412368). "Growth is bound by the area of the squares" is REFUTED.

**The statistic.** The seed-weighted count reproduces the session lead's 248 of 49,999 shortcut peaks in `{(4k+1)^2-1,(4k+1)^2,(4k+1)^2+1}`; within `+-1` of any square: 535; by base residue mod 4: 0, 248, 0, 287; full-map peaks in the `4k+1` window: 0. The seed-weighted count is the wrong statistic because seeds share peaks: there are only 20,341 distinct peak values (top multiplicities 125252 with 448 seeds, 4616 with 408, 638468 with 347, 172772 with 161). Over distinct peaks, with the null "a distinct peak is uniform on the class 2 mod 6 of its dyadic range":

| window | base mod 4 | `3|s`? | distinct peaks | seeds | null expected (distinct) |
|---|---:|---|---:|---:|---:|
| `s^2-1` | 1 | yes | 46 | 152 | 20.4 |
| `s^2-1` | 3 | yes | 47 | 195 | 20.1 |
| `s^2+1` | 1 | no | 49 | 96 | 40.8 |
| `s^2+1` | 3 | no | 50 | 92 | 40.9 |

Under this crude null the `+1` window is at chance (99 distinct vs 81.8, Poisson `z=1.91`) and the `-1` window is not (93 vs 40.6); the split 46/47 between bases 1 and 3 mod 4 already shows that nothing selects `4k+1`. But the crude null is the wrong comparison, for two exact reasons.

**Lemma 8.3 (no shortcut peak is 6 mod 8; PROVED).** If `p = 6 mod 8` then `p/2 = 3 mod 4` is odd, `C(p/2) = (3p+2)/4` is odd (`3(p/2)+1 = 2 mod 4`), and `C((3p+2)/4) = (9p+10)/8 > p`; so the orbit of `p` exceeds `p` and `p` is not a peak. With Lemma 8.1, every shortcut peak is `2, 8` or `20 mod 24` (census of distinct peaks: 6190, 7967, 6184). For odd `s`, `s^2-1 = 0 mod 8` and `s^2+1 = 2 mod 8`: the two windows live in different residue classes, whose peak rates differ.

**Lemma 8.4 (multiplicity, and the ladder is the odd-preimage chain; PROVED).** The odd seeds with shortcut peak `p` are the odd vertices of the preimage tree of `p` kept below `p`. The `C`-preimages of `x` are `2x` and, when `x = 2 mod 3`, the odd `(2x-1)/3`; `2p` and `2m = (4p-2)/3 > p` (`m = (2p-1)/3`) are excluded, so below `p` the tree continues from `m` only through its odd preimage, which exists iff `m = 2 mod 3` iff `9 | p+1`. Iterating, the odd-preimage chain of `p` is `p <- m^(1) <- ... <- m^(v)` with `m^(k) = ((p+1)/3^k) 2^k - 1` and `v = v_3(p+1)`: these are exactly the ladder seeds of Lemma 8.2, each climbs monotonically to `p`, and each has peak `p` because the orbit of `p` never exceeds `p`. Consequences (all checked on the census): (a) a peak has a unique odd seed iff `9 ∤ p+1` (9,461 of the 20,341 distinct peaks; the converse holds whenever `m <= 10^5`); (b) "`p = a 3^k - 1` is the peak of the ladder seed `a 2^k - 1`" holds for every peak and every `k <= v_3(p+1)`, so it explains nothing about squares (the first version of this note offered it as the explanation of the `-1` window: `k=1` is the predecessor `m`, and the check was vacuous); (c) 28,483 of the 49,999 census seeds are ladder seeds and 21,516 branch off through an even vertex `2 m^(k) < p`, `k >= 2`; the ladder tops `6560=3^8-1` (21 seeds), `59048=3^10-1` (48), `164024=25*3^8-1=405^2-1` (17), `4782968=3^14-1` (17) have `v_3(p+1) = 8, 10, 8, 14`, and multiplicity grows with `v_3(p+1)`.

**The seed box (FINITE-EXACT).** Below the box the distinct peaks are exactly the 9,423 values `q = 2 mod 6`, `q < 10^5`, whose orbit never exceeds `q` (no selection by seeds is involved), with peak rate `2868/4166 = 0.688`, `3703/4167 = 0.889`, `0/4167 = 0.000`, `2852/4166 = 0.685` on the classes `2, 8, 14, 20 mod 24`. Above the box a peak must be reached from a seed `< 10^5`, and the ladder seed `((p+1)/3^k) 2^k - 1 ~ p (2/3)^k` is small when `v_3(p+1)` is large: the share of peaks with `9 | p+1` is 0.335 below the box (null `1/3`) and 0.707 above it.

**The stratified null.** Rate of being a distinct peak by (dyadic range, `q mod 24`, `min(v_3(q+1),4)`), with exact stratum sizes:

| window | observed (below box, above) | expected (below, above) | `z` | base 1 mod 4 obs / exp | base 3 mod 4 obs / exp |
|---|---|---|---:|---|---|
| `s^2-1` | 93 (44, 49) | 88.7 (45.5, 43.1) | 0.46 | 46 / 44.6 | 47 / 44.1 |
| `s^2+1` | 99 (68, 31) | 100.9 (71.2, 29.7) | -0.19 | 49 / 50.7 | 50 / 50.2 |

Both windows and both base classes are at chance (the audit also checks that below the box the result is unchanged under `q mod 1152` and `q mod 27648`, and that stratifying by `v_2` instead of `q mod 24` would wrongly report a `+1` excess, because `v_2 = 1` lumps `2 mod 8` with the rate-zero class `6 mod 8`). The crude `-1` excess is the product of `8 | s^2-1` with the peak rate 0.889 of the class `8 mod 24` (Lemma 8.3) and the box bias to `9 | p+1` (Lemma 8.4); neither fact is about squares. The shifted-window control `{(4k+1)^2+s-1,+s,+s+1}` gives 23, 0, 24, 115, 91, 152, 248, 96, 0, 150, 150, 21, 21 for `s=-12,-10,...,12`: shifts 6 and 8 each score 150 seeds, so the window at `s=0` is not distinguished. Verdict: "the accelerated Syracuse paths track the boundaries of perfect squares as they scale" is REFUTED; the surviving exact content is Lemmas 8.2-8.4, of which 8.2 is the hostile to the paste's own bound.

**196, 169, 160.** The first 12 primes `2,3,5,7,11,13,17,19,23,29,31,37` sum to 197, not 196 (REFUTED). `196=4 mod 6` is never a shortcut peak; its orbit `196, 98, 49, 74, 37, 56, 28, 14, 7` reaches 7 in 8 steps, so 196 is a predecessor of 7, not a "mirror". Odd seeds `<=10^5` with shortcut peak 170 are 75 and 113; none has peak 168 or 160. `160=2^5*5` is TRUE arithmetic and 160 is the full-map peak of 15, 23, 35, 53 (`53 -> 160`, shortcut image 80); `160=4 mod 6` can never be a shortcut peak, and no "compression deficit between 169 and 160" is defined (SCOPE).

## 9. The descent certificate and the open obligation (inherited)

With halving word `k_1,...,k_L`, `K_j=k_1+...+k_j`, `K=K_L`, `B=sum_{j=0}^{L-1} 3^{L-1-j} 2^{K_j}`, the inherited identity is `2^K T^L(n)=3^L n+B`, hence `T^L(n)<n iff 3^L n+B<2^K n` ([blueprint synthesis](collatz_blueprint_20260921_synthesis.md), section 6). The script re-verifies the identity and prints the margin `2^K n-(3^L n+B)`: `n=27, L=5`: `-5120`; `n=7, L=2`: `-40`; `n=23, L=3`: `2304`; `n=1, L=1`: `0`; `n=31, L=5`: `-5760`; `n=2^20-1, L=20`: `-3638566269747200`; `n=27, L=41`: positive (`K=70`, the word that reaches 1). For `n=27` the first `L` with positive margin is `L=37` (`K=59`, `T^37(27)=23`; `37+59=96` is the classical stopping time `sigma(27)`, the first full-map step below 27; the total stopping time, the first arrival at 1, is `41+70=111`; both re-verified independently in the script). The paste's inequality `3^L n+B_L<2^{K_L} n` is exactly this certificate; its claim that "peaks run out of volume" is not a statement about any fixed `n` and proves nothing. The obligation is unchanged: for every fixed odd `n>1`, exhibit a word of its own orbit with positive margin.

## 10. Claim 6 and the remaining phrases (SCOPE)

The mod-30 wheel has the 8 reduced classes `1, 7, 11, 13, 17, 19, 23, 29`; it is the only defined object in "the mod 30 primes framework". "Unified arithmetic cycle `{f,+,*,g}`", "modular stabilizer `q=x (mod AB)`", "macro-growth fractures back into local modular grids", "8-step Bott periodicity transforming tetration into pentation" (Bott periodicity is the 8-fold periodicity of the stable homotopy of `O`, Bott 1959, unconnected to hyperoperations), "quasicrystalline Fourier transform matrix": no definitions, no map found.

## 11. What survives

- The corrected core identity `T((A+1)(B+1)-2)-T(AB)=T(A)+T(B)+A(B^2-1)+B(A^2-1)`, its true unshifted form `T((A+1)(B+1)-1)-T(AB-1)=T(A)+T(B)+2A T(B)+2B T(A)`, and the flows' unique solution `S=AB(A+B)+C(A,2)+C(B,2)` (classical, no novelty claim).
- The simplex convolution laws, all instances of Chu-Vandermonde: `C(A+B+d,d)=sum_i C(A+a,i)C(B+b,d-i)`.
- The count `T(n-2)=C(n-1,2)` of non-path edges, and `sum_T h(T)=n! 2^{T(n-2)}`; Redei's parity and THM-1745's `{7,21}` holes; no compression.
- Bang's unique exception `63=7*9` through Catalan `9=2^3+1`; Zsigmondy without exception for `4^n-1`; Cipolla's trunk pseudoprimes `(4^p-1)/3` on the `R`-orbit of 1.
- The descent certificate `2^K T^L(n)=3^L n+B` and the exact inverse-fibre braid `R(n)=4n+1` (inherited).
- The ladder `C^j(a 2^k-1)=a 3^j 2^{k-j}-1`, which refutes bounded growth and is the odd-preimage chain of every peak (Lemma 8.4, multiplicity 1 iff `9 ∤ p+1`); the exclusion `p ≠ 6 mod 8` (Lemma 8.3); and the residue-stratified null under which the paste's square windows are at chance.

## 12. What this means for a Collatz proof

Nothing in the paste advances the descent problem. The two arithmetic facts it states about Collatz, the certificate `3^L n+B_L<2^{K_L} n` and the trunk `(4^j-1)/3`, are already in the repo and are exactly the objects the blueprint audit identified as the correct replacements for the geometric language. The quadratic-horizon picture is refuted at the level of its own data (bases 1 and 3 mod 4 indistinguishable; both square windows at chance once peaks are compared within their residue class mod 24 and their `v_3(p+1)` stratum; the only structures present are the mod-24 exclusion and the `3^k-1` ladder, which is a growth mechanism, not a bound), and the "area" argument is not an argument about any specific `n`. Consistent with the [blueprint audit](collatz_blueprint_20260921_synthesis.md), the open obligation is unchanged and precise: for each fixed odd `n>1`, a word of its own orbit with `2^K>3^L` and `n>B/(2^K-3^L)`. Residue coverage, finite-word completion, and any density statement supply witnesses for prefixes, not this statement for a given `n`.

## Reproduction

```text
python3 04-computation/experiments/collatz_mod6_20260917_scaffolding_audit.py > 05-knowledge/results/collatz_mod6_20260917_scaffolding_audit.out
```

Stdlib plus sympy (factorisation only; every factorisation is re-multiplied). Exact integer, Fraction and polynomial arithmetic; all checks are explicit `raise`s and the `python3 -O` run is byte-identical to the normal run. Runtime about 20 s, memory far below 1 GB. Independent audit: `python3 04-computation/experiments/collatz_mod6_20260917_scaffolding_audit_audit.py > 05-knowledge/results/collatz_mod6_20260917_scaffolding_audit_audit.out` (about 3 s; sympy symbolic identities, cyclotomic-value primitive-prime criterion, permutation-matching tournament census, independent peak census; also `-O`-identical). Corrections made after computing the truth: the crude seed-weighted null of section 8 (replaced by the distinct-peak cell table); a stale assertion that 27 first descends at `L=41` (the truth is `L=37`; 41 is its total number of odd steps); the first version's tautological "ladder explanation" of the `s^2-1` window (replaced by Lemmas 8.3-8.4 and the stratified null, under which the window is at chance); and the label "total stopping time" for 96 in section 9 (96 is the stopping time; 111 is the total). Section 6 was found to duplicate the wild_typing lane and is now cited as inherited.

## Stopping boundary / next question

Stopped once every sentence of the paste had a verdict with a witness or a proof pointer and the only signal in its data was reduced to residue classes and the seed box. The question the first version left ("for which ladder seeds is `a 3^k-1` the orbit maximum?") is answered by Lemma 8.4: for every `k <= v_3(p+1)` whenever `p` is a peak at all, so it was circular. The small real question that remains is the residue-class peak rate itself: below the box the values `q = 2 mod 6` whose orbit never exceeds `q` have densities `0.688, 0.889, 0.685` on `2, 8, 20 mod 24` (and 0 on `14 mod 24`, proved). Whether these densities exist as `X -> infinity` and are computable from the parity-word distribution (the probability that a 2-adic word never climbs above its start, restricted to a residue class) is OPEN; no claim is made that it helps the global obligation, which is unchanged.
