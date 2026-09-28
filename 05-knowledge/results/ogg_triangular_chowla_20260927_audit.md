# Audit of `ogg_triangular_chowla_20260927.md` (Ogg's fifteen primes, the nodal cubic against triangular numbers, the Chowla map)

**Verdict: SOUND WITH CORRECTIONS.** Every theorem-level claim (Proposition 1, Proposition 2, the genus / class-number / Dirichlet / Hasse-polynomial statements, the group-theoretic facts) is correct and every number in sections 2–6 reproduces, with these corrections: the eleven "`nP + mQ`" identities hold only up to sign (seven of eleven need `-(nP + mQ)`); the "no relation" claim is not tested by the note's script (it is true, and more: `P, Q` are independent, so `rank E(Q) >= 2` is now PROVED, not "consistent with"); the nodal cubic is the Legendre fibre at `lambda = 1` (not `lambda = 0`, whose `x -> -x` image is the `-1` twist) and is a cusp, not a split node, at `p = 2`; Cattaneo's theorem is stronger than Prop 1(c), not its "odd case"; `118` is the next `PSL(2,q)` Hurwitz genus, not the next Hurwitz genus (17 precedes it); the aliquot comparison figure `-5.34` is quoted from another range; the paper does contain a (minor) inaccurate statement. No unsupported PROVED label.

**Auditor:** independent subagent, 2026-09-27, worktree `math-wt-collatz-poset-dag-20260927`.
**Script:** `04-computation/experiments/ogg_triangular_chowla_20260927_audit.py` (written from scratch, exact arithmetic; 19 s), output `05-knowledge/results/ogg_triangular_chowla_20260927_audit.out` (232 lines; 1 flagged check, which is the provenance issue C10 below, not an error in the mathematics).

**Method.** The statements of Propositions 1 and 2 and the classical statements of sections 2 and 3 were re-derived before comparison (the derivations are recorded in the `.out` header of section A and in the items below); the Read tool returned the whole note, so "blind" means independently re-derived, not unseen. The paper text (`paper30.txt`, 4 pages) was read in full for section 1. The parallel notes were consulted only to locate the provenance of quoted numbers (STICKY Theorem 1, the `-5.34` figure) and to confirm that the cited IDs exist.

---

## 1. Itemized claims

### A. Proposition 1 (parity alternation of `s'(n) = sigma(n) - n - 1`)

| claim | status | evidence |
|---|---|---|
| `sigma(n)` odd iff `n` is a square or twice a square | CONFIRMED | re-derived (`sigma(2^a)` odd; `sigma(q^e) = e + 1 (mod 2)` for odd `q`); checked for all `n <= 2·10^6` |
| (a) `n = 2^a m` even: `s'(n)` odd iff `m` not a square | CONFIRMED | `s'(n) = sigma(n) + 1 (mod 2)`; checked `n <= 2·10^6` |
| (b) `n` odd: `s'(n)` even iff `n` not a square | CONFIRMED | `s'(n) = sigma(n) (mod 2)`; checked `n <= 2·10^6` |
| (c) fixed points are squares or twice squares | CONFIRMED | parity; no fixed point `<= 2·10^6` |
| (c) odd cycles and same-parity 2-cycles contain a square-type element | CONFIRMED | an odd number of parity flips cannot close a cycle; a same-parity 2-cycle has `s'(n) = n (mod 2)` at a member |
| (c) proportion with `s'(n) = n (mod 2)` is `O(N^(-1/2))` | CONFIRMED, constant sharpened | the set is EXACTLY the square-type numbers `2^a k^2` (`k` odd): count in `[M, 2M)` is `~ sqrt(M/2)`, share `(2M)^(-1/2)`; see C11 |
| "Cattaneo's theorem ... is the odd case" | CORRECTED (C5) | Cattaneo's theorem also excludes the even (twice-a-square) case; Prop 1(c) does not |
| attribution "Cattaneo 1951: quasiperfect numbers are odd squares" | UNVERIFIED (recollection, consistent) | P. Cattaneo, Boll. Un. Mat. Ital. 1951, is the standard reference for "quasiperfect => odd square"; not checked against the source |
| Hagis–Lord 1977 (quasi-amicable pairs) | UNVERIFIED (recollection, consistent) | Math. Comp. 31 (1977); the eighteen pairs below are all in the known list |

### B. Genus of `X_0(p)`, fixed points of `w_p`, Ogg's list

| claim | status | evidence |
|---|---|---|
| `g(X_0(p)) = 1 + (p+1)/12 - nu_2/4 - nu_3/3 - 1`, `nu_2 = 1 + (-4/p)`, `nu_3 = 1 + (-3/p)`, `nu_3 = 0` at `p = 2`, `nu_2 = 0` at `p = 3` | CONFIRMED | own Kronecker implementation gives `(nu_2, nu_3) = (1, 0)` at 2, `(0, 1)` at 3; integrality for all `p < 1000`; equals `floor((p+1)/12) - [p = 1 (12)]` for `p >= 5`; matches the standard table to `p = 113` |
| `f = h(-4p) + h(-p)·[p = 3 (4)]`, `p = 2` handled as `h(-8) + h(-4) = 2` | CONFIRMED, caveat (C8) | literature convention (Ogg 1974, Prop. 3, recollection): `h(-4p)` for `p = 1 (4)`, `h(-4p) + h(-p)` for `p = 3 (4)`, `p > 3`; `p = 2, 3` have 2 fixed points. Independent numerical confirmation: `f = 2·#{supersingular j in F_p}` (point counting, no class numbers) for all `5 <= p < 200`. The naive `h(-8) = 1` fails at `p = 2`; `p = 3` needs no patch (`h(-12) + h(-3) = 2`) |
| class numbers "counted as reduced primitive forms" | CONFIRMED | own reduced-form counter agrees with Dirichlet's analytic formula for all 305 fundamental `D in (-1000, 0)`, and with `h(-4p) = 3h(-p)` (`p = 3 (8)`), `h(-p)` (`p = 7 (8)`) for `3 < p < 1000` |
| `g(X_0(p)^+) = (2g + 2 - f)/4` | CONFIRMED | Riemann–Hurwitz for an involution; `2g + 2 - f = 0 (mod 4)` for all `p < 1000` |
| genus 0 exactly for the fifteen primes, `p < 120` | CONFIRMED and extended | exactly `{2, ..., 71}` for ALL `p < 1000`; `g^+ = 1` for `p < 200`: `37, 43, 53, 61, 79, 83, 89, 101, 131`; `g^+ = 2`: `67, 73, 103, 107, 167, 191` (agree with the standard lists, recollection) |
| table rows `11: 1,4,0; 37: 2,2,1; 71: 6,14,0; 73: 5,4,2; 97: 7,4,3` and the whole `p <= 100` row of the `.out` | CONFIRMED | identical |
| `h(-3), h(-7), h(-11), h(-19), h(-23), h(-31), h(-47), h(-59), h(-71) = 1,1,1,1,3,3,5,3,7` | CONFIRMED | |
| "Ogg's condition in one line: `2g + 2 = h(-4p) + h(-p)[p = 3 (4)]`" | CORRECTED (C8, minor) | needs `p >= 3` |

### C. Supersingular route

| claim | status | evidence |
|---|---|---|
| Hasse polynomial `H_p(lambda) = sum C(m,k)^2 lambda^k`, roots in `F_(p^2)`; `j = 256(lambda^2 - lambda + 1)^3/(lambda^2(lambda - 1)^2)` | CONFIRMED | `H_p` has `(p-1)/2` simple roots, all in `F_(p^2)` (sympy factorisation over `GF(p)`: degrees `<= 2`, multiplicity 1) for `5 <= p <= 97`; the `j`-sets obtained from the roots agree with an independent enumeration of `j in F_p` by point counting (`a_p = 0`) and by the Hasse invariant (coefficient of `x^(p-1)` in `(x^3 + ax + b)^((p-1)/2)` on an `F_p`-model), satisfy the Eichler–Deuring mass formula `sum w_j = (p-1)/12`, and have the count `floor(p/12) + eps(p)` |
| counts `1, 1, 2, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 5, 6, 7, 6, 7, 8, 8, 8` | CONFIRMED | |
| all in `F_p` exactly for `{5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71}`; `37: 3, 1`; `97: 8, 2` | CONFIRMED | full table in the `.out` (section C) |
| "Deuring" for the Hasse polynomial | UNVERIFIED (recollection, consistent) | Deuring 1941 / Igusa 1958 |

### D. Dirichlet's half-row formula and the Paley reading

| claim | status | evidence |
|---|---|---|
| `h(-p) = (sum_(0<a<p/2) (a/p))/(2 - (2/p))`, `p = 3 (4)`, `p > 3` | CONFIRMED | exact for all 154 primes `p = 3 (4)`, `3 < p < 2000`; equivalently the sum is `h(-p)` for `p = 7 (8)` and `3h(-p)` for `p = 3 (8)`; fails at `p = 3` (`1/3`), so `p > 3` is necessary |
| the sum is the out-minus-in degree of vertex 0 into `{1, ..., (p-1)/2}` in the Paley tournament | CONFIRMED | arc `a -> 0` iff `-a` is a residue iff `(a/p) = -1`; checked explicitly for `p = 7, 11, 19, 23, 31, 47, 59, 71` |
| typing "STRUCTURAL and classical ... not a mechanism" | CONFIRMED as typed | |

### E. Proposition 2 (the nodal cubic and `E: Y^2 = X^3 - 4X + 1`)

| claim | status | evidence |
|---|---|---|
| (a) points `(t^2 - 1, t^3 - t)`; group `G_m` via `(y + x)/(y - x) = (t + 1)/(t - 1)`; tangents `y = ±x`; integer points `y = 6 C(t+1, 3)` | CONFIRMED | own chord–tangent implementation on `y^2 = x^3 + x^2`: the map is multiplicative on 300 random rational pairs; derived parameter law `t(P1 + P2) = (1 + t1 t2)/(t1 + t2)`; integer points have `x + 1 = t^2`, `t` integer |
| (b) `x` triangular iff `(2m+1)^2 - 8t^2 = -7`; `t = 1, 2, 4, 11, 23, 64, 134, 373, 781, 2174, 4552, 12671, ...`; two orbits under `3 + sqrt 8`; `m = 0, 2, 5, 15, 32, 90, 189, 527` | CONFIRMED | list to `t <= 10^6` (17 terms); the two interleaved sequences are the orbits of `(1, 1)` and `(5, 2)` under `(u, t) -> (3u + 8t, u + 3t)`, each obeying `t_(k+2) = 6t_(k+1) - t_k` |
| (c) `y` triangular iff `8t^3 - 8t + 1` square iff `(2t, 2k+1)` on `Y^2 = X^3 - 4X + 1`; discriminant `2^4·229` | CONFIRMED | model minimal (`v_p(Delta) < 12`); conductor `916 = 2^2·229` by Tate's algorithm (type IV at 2: `b_6 = 4` not divisible by 8; `I_1` at 229) — own computation, no table lookup |
| `P = (0, 1)` of infinite order (`3P` has `x = -7/4`) | CONFIRMED | `2P = (4, 7)`, `3P = (-7/4, 13/8)`; torsion is in fact trivial (`#E(F_3) = 7`, `#E(F_5) = 9`) |
| integer points `-2 <= X <= 10^6`: the eleven listed | CONFIRMED and extended to `X <= 10^7` | no further point |
| "every one of them is `nP + mQ`: `P + Q, P - Q, P, Q, 2P + Q, 2P, 2P - Q, 2Q, 2P + 2Q, 4P + Q, 2P - 3Q`" | CORRECTED (C1) | true up to sign only; exact list in C1 |
| "no relation with `|n|, |m| <= 12`; consistent with rank two, not proved" | CORRECTED (C2) | the note's script does not test relations; audit: none with `|n|, |m| <= 30`; and `P, Q` are independent (mod-37 test), so `rank >= 2` is proved |
| `E(F_p)` structure "incompatible with `Q` a multiple of `P`" (task (e)) | REPORTED | `p = 37`: `#E = 44`, `ord P = ord Q = 22`, `<P, Q> = E(F_37)` NOT cyclic; `p = 53`: `#E = 48`, orders `12, 24`, `<P, Q>` of order 48, NOT cyclic; hence `Q` is not a multiple of `P` modulo torsion (torsion being trivial, not at all) |
| `t <= 10^6`: `t = 1, 2, 5, 6, 10, 57, 637`, `y = 0, 6, 120, 210, 990, 185136, 258474216 = T_0, T_3, T_15, T_20, T_44, T_608, T_22736` | CONFIRMED and extended to `t <= 5·10^6` | |
| "The even `X` give `t = 1, 2, 5, 6, 10, 57, 637`" | CORRECTED (C9, trivial) | also `X = -2, 0`, i.e. `t = -1, 0` (both `y = 0`) |
| (d) both triangular, `t <= 10^6`: only `(0, 0)` and `(3, 6)` | CONFIRMED to `t <= 5·10^6` | |
| "the perfect numbers are not on the Pell family" (`p = 2, 3` fail, "none below `p = 31`") | CONFIRMED for prime `p` (all `p <= 100`), addition (C13) | but `T_15 = 120 = 2^3(2^4 - 1)` (Euclid shape, `k = 4` not prime) IS on the family: `121 = 11^2`, the note's own `t = 11` point |
| cannonball: `n = 1, 24` only below `10^6` | CONFIRMED | |
| "the Legendre curve at `lambda = 0` (after `x -> -x`)", "split multiplicative reduction everywhere" | CORRECTED (C3, C4) | `lambda = 1` after `x -> x + 1` exactly; `lambda = 0` with `x -> -x` is the `-1` twist; cusp at `p = 2` |

### F. Section 5 dynamics (all `n <= 10^6`, own sieve to `10^7`, sympy check of the cycles)

| quantity | note / `.out` | audit |
|---|---|---|
| endings | `0: 982099`, `cycle2: 17900`, no fixed point, no longer cycle | identical |
| longest | 48 steps, `n = 948375` | identical |
| largest peak | `6.3 n` at `n = 980100` | `6148079/980100 = 6.273`, `n = 980100 = 990^2` |
| mean length | 9.32 | 9.321 |
| the eighteen 2-cycles and their basin counts | as listed | identical list AND identical counts (9724, 3485, 558, 1971, 963, 408, 441, 226, 72, 6, 7, 14, 5, 3, 12, 3, 1, 1) |
| each pair a genuine 2-cycle; all opposite parity | claimed | CONFIRMED with `sigma` from sympy; no member is square-type |
| drifts `-0.048` (even), `-2.82` (odd composite) | | `-0.048`, `-2.817` |
| aliquot comparison `-5.34` (odd) | hard-coded string in the script | CORRECTED (C10): `-5.26` on `n <= 10^6`; `-5.34` is the STICKY value on `n <= 2·10^6` |
| persistence `0.0254, 0.0079, 0.0025, 0.0008` | | `0.02543, 0.00792, 0.00245, 0.00076` (identical) — and equal to (square-type composites)/(composites) exactly (C11) |
| odd `n <= 10^5` with `s'(n) > n`: 210, first 945 | | identical; identical to the odd abundant numbers below `10^5` |
| "values beyond the sieve `4·10^6`: 3" | | consistent (peak `6148079`) |
| addition | | nine further betrothed pairs with smaller member in `(10^6, 2·10^6]` exist (`(1000824, 1902215)`, `(1081184, 1331967)`, ...), none reached from a start `<= 10^6`; all of opposite parity |

### G. Numbers of sections 2, 3, 6

CONFIRMED: `378 = T_27 = 2·3^3·7` = sum of the fifteen (`160 + 218`); Monster `sopfr = 637 = 7^2·13`; primes `< 72` missing: `37, 43, 53, 61, 67`; excess `1/6, 1/12, 1/30, 0, -1/42` for `q = 3..7`, orders `12, 24, 60`; `|PSL(2,7)|, |PSL(2,8)|, |PSL(2,13)| = 168, 504, 1092 = 84(g - 1)` for `g = 3, 7, 14`; factorisations `2^3·3·7`, `2^3·3^2·7`, `2^2·3·7·13`; `|Aut(Paley T_7)| = 21` (orders `1^1 3^14 7^6`, nonabelian, `F_21`); in `PSL(2,7)` (168 elements constructed) the normaliser of a Sylow 7-subgroup has order 21, is nonabelian, index 8, and there are 8 Sylow 7-subgroups; `120 = 5! = T_15`, `210 = 7# = T_20`, `24, 35` at `t = 5, 6`; `T_189 = 134^2 - 1 = 17955`; `189 = 3^3·7`; Mersenne primes in the list `3, 7, 31` (`127` not); `j = 0` supersingular at `p = 2 (3)` for `p >= 5`.
CORRECTED: `118` "next arithmetic Hurwitz genus" (C6); `j = 0` "exactly at `p = 2 (mod 3)`" also at `p = 3` (C12, trivial).
Cross-references: THM-640, THM-448, HYP-2220 exist; HYP-3805 is listed in the investigation backlog as REFUTED-BUT-CLOSE (corrected formula open) — the note cites it as "the Paley heptagon as the LRC extremal object" without that status (C14).

### H. Section 1: the typing of the paper

| claim | status |
|---|---|
| the paper distinguishes Ogg's list from Lang–Trotter primes and states a "regime-separation principle" and Proposition 1 | CONFIRMED (read in full) |
| Theorem 1 (Ogg) as stated in the paper is correct | CONFIRMED (both routes reproduced) |
| the paper never mentions the Monster | CONFIRMED (0 occurrences; "moonshine" 3) |
| "overclaim in framing": "a complete explanation ... without appeal to moonshine" explains the list, not the coincidence with the Monster | AGREED — fair reading; the paper's abstract/conclusion do use "moonshine-independent ... complete explanation" for a genus computation Ogg made in 1975 |
| reference [2] is Ogg, Bull. SMF 102 (1974) "Hyperelliptic modular curves" | CONFIRMED (paper's reference list) |
| the list / supersingular characterisation is in Ogg 1975 (Sém. Delange–Pisot–Poitou, "Automorphismes de courbes modulaires"), not the 1974 paper | UNVERIFIED (recollection, consistent with the auditor's): the 1975 exposé contains the fifteen-prime observation and the Jack Daniels remark; the 1974 paper contains the `w_N` fixed-point count (used here), so it is not irrelevant to the genus computation |
| "nothing in it is wrong as a statement" | CORRECTED (C7) |
| does the note over/understate the paper? | Neither materially. It could add that the paper's own "structurally complete" explanation is a two-line consequence of Deuring's count (`#ss j in F_p = f/2`) plus `#ss = g + 1`, which the paper does not spell out. |

### I. Labels

`PROVED` for Propositions 1 and 2 (a)–(c) as stated: justified. Section 2's "PROVED classical + FINITE-EXACT": the class-number reading is verified numerically (`p < 120`, `p < 200`), not re-derived; the two-route equivalence is checked, not derived (derivation: `#ss(p) = g + 1`, `#ss in F_p = f/2`, so all-in-`F_p` iff `f = 2g + 2` iff `g^+ = 0`). "Rank at least one" is correct but weaker than what the data prove (C2). No misattribution found beyond C5's mis-description of Cattaneo's theorem. NUMEROLOGY typing of section 6 is appropriate and honest.

---

## 2. Corrections (verbatim -> repaired)

**C1 (section 4(c), sign of the combinations).** Verbatim: "and every one of them is `nP + mQ` with `Q = (2, 1)`: `P + Q, P - Q, P, Q, 2P + Q, 2P, 2P - Q, 2Q, 2P + 2Q, 4P + Q, 2P - 3Q`". The script compared `|Y|` only. Exact: `(-2, 1) = -P - Q`, `(-1, 2) = -P + Q`, `(0, 1) = P`, `(2, 1) = Q`, `(3, 4) = -2P - Q`, `(4, 7) = 2P`, `(10, 31) = -2P + Q`, `(12, 41) = -2Q`, `(20, 89) = 2P + 2Q`, `(114, 1217) = -4P - Q`, `(1274, 45473) = -2P + 3Q` (e.g. `2Q = (12, -41)`, `P + Q = (-2, -1)`). Repaired statement: "every one of them is `±(nP + mQ)` with the listed `(n, m)`", or the list above.

**C2 (section 4(c), relations and rank).** Verbatim: "(no relation with `|n|, |m| <= 12`; the data are consistent with rank two, not proved)". The note's script never tests `nP + mQ = O` (it only records integral results). Audit: no relation with `|n|, |m| <= 30`; moreover `E(Q)_tors = 0` (`gcd(#E(F_3), #E(F_5)) = gcd(7, 9) = 1`) and at `p = 37` the cubic `X^3 - 4X + 1` splits, `E(F_37) = Z/2 x Z/22`, and `P, Q, P + Q` all lie outside `2E(F_37)`, so `P, Q` are independent in `E(F_37)/2E(F_37) = (Z/2)^2`; a relation `nP + mQ = O` with `gcd(n, m) = 1` (torsion trivial) would reduce to one there. Repaired: "`P` and `Q` are independent: `rank E(Q) >= 2` (PROVED by reduction mod 37; also mod 53, 173, 241, ...); `rank = 2` is not proved (needs a 2-descent over `Q(theta)`, `theta^3 - 4theta + 1 = 0`, or a table lookup at conductor 916)". Section 0's "of rank at least one" and T1556's "rank `>= 1`" can be upgraded to `>= 2`.

**C3 (section 2, the Legendre fibre).** Verbatim: "The owner's `y^2 = x^2(x+1) = x(x - 0)(x + 1)` is the Legendre curve at `lambda = 0` (after `x -> -x`)". `x -> -x` sends `y^2 = x^2(x - 1)` to `y^2 = -x^2(x + 1)`, the quadratic twist by `-1` (non-split node at `p = 3 (mod 4)`). Repaired: "is the Legendre fibre at `lambda = 1` after `x -> x + 1` (`y^2 = x(x - 1)^2 -> y^2 = (x + 1)x^2`), or equivalently the fibre at `lambda = 0` up to the twist by `-1`".

**C4 (section 2, reduction at 2).** Verbatim: "it has split multiplicative reduction everywhere (`a_p = 1`, `#E_ns(F_p) = p - 1`)". At `p = 2`, `y^2 = x^3 + x^2` has `(y - x)(y + x) = (y + x)^2`: a cusp, `E_ns = G_a`, `#E_ns(F_2) = 2 = p`. Repaired: "at every odd prime `p` (split node, `#E_ns(F_p) = p - 1`); at `p = 2` it is cuspidal".

**C5 (section 5, Prop 1(c), and section 0).** Verbatim: "(odd fixed points are odd squares: Cattaneo's theorem for quasiperfect numbers is the odd case)" and "Cattaneo's 'odd squares' is the theorem's odd case". Prop 1(c) gives only "a fixed point is a square or twice a square"; Cattaneo's theorem is strictly stronger: it also excludes `n = 2^a m^2` (`a` odd, `m` odd), since `sigma(n) = 2n + 1` forces `m^2 = -1 (mod 2^(a+1) - 1)` with `2^(a+1) - 1 = 3 (mod 4)`. Repaired: "fixed points are squares or twice squares (parity); Cattaneo (1951) proves more, that a quasiperfect number is an odd square, the twice-a-square case being excluded by a residue argument that parity alone does not give".

**C6 (section 6 table).** Verbatim: "| `118 = 2·59` | next arithmetic Hurwitz genus |". The Hurwitz genera are `3, 7, 14, 17, 118, 129, 146, ...`; `17` (group of order `1344 = 84·16`, recollection) precedes `118`, and every Hurwitz surface is arithmetic. Repaired: "`118` = genus of the `PSL(2,27)` Hurwitz surface (`|PSL(2,27)| = 9828 = 84·117`), the next `PSL(2,q)` Hurwitz genus after `14`".

**C7 (section 1).** Verbatim: "nothing in it is wrong as a statement". The paper's "Lang and Trotter conjectured that such primes are infinite but extremely sparse, with asymptotic growth of order `sqrt(X)/log X`" omits the non-CM hypothesis (for CM curves the supersingular primes have density `1/2`); and "The finiteness of the list follows from the classification of genus-zero modular curves" is a loose justification (finiteness is immediate from `g(X_0(p)^+) >= (p/12 - O(sqrt p log p))/2 -> infinity`). Repaired: "nothing in it is wrong beyond a missing non-CM hypothesis in the Lang–Trotter sentence and a loosely justified finiteness remark".

**C8 (section 0 and section 2, bold line).** Verbatim: "`2 g(X_0(p)) + 2 = h(-4p) + h(-p)·[p ≡ 3 (4)]`". Holds for `p >= 3`; at `p = 2` the right side is `h(-8) = 1`, not the 2 fixed points of an involution of `P^1` (the note's `h(-8) + h(-4) = 2` is a convention, stated in section 2 only). Repaired: add "`p >= 3`; `p = 2` has `f = 2`".

**C9 (section 4(c), trivial).** Verbatim: "The even `X` give `t = 1, 2, 5, 6, 10, 57, 637`". The even `X` are `-2, 0, 2, 4, 10, 12, 20, 114, 1274`, i.e. `t = -1, 0, 1, 2, 5, 6, 10, 57, 637`; `t = -1, 0, 1` give `y = 0 = T_0`. Repaired: "the even `X >= 2` give `t = 1, ...`".

**C10 (section 5, provenance of a comparison figure).** Verbatim: "(against `-5.34` for `s` on odd `n`, where primes give `s = 1`)". The number is a literal in the script's print statement (`(aliquot s(n): even -0.048, odd -5.34)`), not computed; it is the STICKY note's value on `n <= 2·10^6` (`-5.342`). On the range of the Chowla figures (`n <= 10^6`) the aliquot odd drift is `-5.26` (`-4.93` on `n <= 10^5`, `-5.53` on `n <= 10^7`; `-5.14, -5.43, -5.60` on `[10^5, 2·10^5)`, `[10^6, 2·10^6)`, `[5·10^6, 10^7)` — it drifts with the range because primes contribute `-log_2 p`). Repaired: "`-2.82` (Chowla, odd composite `n <= 10^6`) against `-5.26` (aliquot, odd `n <= 10^6`, primes included)", or cite the STICKY range.

**C11 (section 5 and the `.out`, the persistence constant).** Verbatim (`.out`): "squares' share ~ 1/sqrt(M)"; note: "against the square-type share `M^(-1/2)`". By Prop 1 the persisting `n` in `[M, 2M)` are exactly the square-type numbers `2^a k^2`, `k` odd, whose count is `~ sqrt(M/2)` (`22, 71, 224, 707` at `M = 10^3..10^6`), share `(2M)^(-1/2) = 0.707 M^(-1/2)` among all `n` and `(2M)^(-1/2)·M/#composites` among composites: `0.0259, 0.0079, 0.0024, 0.0008`, which is what the note "measured" (`0.0254, 0.0079, 0.0025, 0.0008`). Repaired: "persistence = (square-type composites)/(composites) exactly, `~ (2M)^(-1/2)`"; the `≍ N^(-1/2)` wording of section 0 is right, the constant `M^(-1/2)` is off by `sqrt 2`, and the numbers are a count forced by the theorem, not an independent statistic.

**C12 (section 3, trivial).** Verbatim: "the CM curve `y^2 = x^3 + 1` with `j = 0`, supersingular exactly at `p ≡ 2 (mod 3)`". True for `p >= 5`; `j = 0` is also supersingular at `p = 3` (and `p = 2`, which is `2 (mod 3)` anyway). Repaired: add "`p >= 5`".

**C13 (section 4, "Perfect numbers", addition rather than error).** Verbatim: "none below `p = 31`): the perfect numbers are not on the Pell family. No connection." Confirmed for all prime `p <= 100`. But the Euclid-shaped `T_(2^k - 1)` with `k = 4`, `T_15 = 120 = 2^3(2^4 - 1)`, IS on the Pell family (`121 = 11^2`; it is the `t = 11`, `m = 15` point listed in (b)). Suggested addition: "the only Euclid-shaped `2^(k-1)(2^k - 1)` with `k <= 100` on the family is `k = 4` (not prime)".

**C14 (cross-reference status).** HYP-3805 is cited as "the Paley heptagon as the LRC extremal object" (sections 0, 2, 3, 6); the repo's investigation backlog lists HYP-3805 as REFUTED-BUT-CLOSE ("corrected formula open"). The citation should carry that status.

---

## 3. Numbers recomputed: agreement summary

All numbers in sections 2, 3, 4, 6 agree (class numbers, genera, fixed points, supersingular counts, Pell list and indices, elliptic-curve integer points, `T_k` indices, `378`, `637`, group orders, factorisations, `24/120/210/189`). All numbers of section 5 agree to the printed precision (endings, longest, peak, mean, the eighteen pairs with basin counts, drifts, persistence, odd abundant count), with the single exception of the quoted `-5.34` (C10). Sign-sensitive identities in section 4(c): 7 of 11 disagree (C1).

## 4. Recollection flags (not checked against sources)

Ogg 1975 Sém. DPP as the home of the fifteen-prime observation (consistent); Ogg 1974 Prop. 3 as the `w_N` fixed-point count (consistent); Cattaneo 1951 (consistent); Hagis–Lord 1977 (consistent); the genus-17 Hurwitz group of order 1344; the standard `g^+ = 1, 2` lists (my computed lists agree with what I recall).

## 5. What remains open

* `rank E(Q)` exactly (proved `>= 2` here; `= 2` would need a 2-descent or LMFDB at conductor 916); the complete list of integer points (Siegel gives finiteness; searched to `X <= 10^7`, `t <= 5·10^6`).
* Termination of the Chowla map: no theorem, as the note says; the parity theorem alone does not exclude alternating growth.
* Whether the Chowla note's "three regimes" table is more than an organising analogy: not a mathematical claim, not audited.

## 6. Reproduction

```text
python 04-computation/experiments/ogg_triangular_chowla_20260927_audit.py > 05-knowledge/results/ogg_triangular_chowla_20260927_audit.out
```
numpy + sympy, about 19 s (sigma sieve to `10^7`: 11 s).
