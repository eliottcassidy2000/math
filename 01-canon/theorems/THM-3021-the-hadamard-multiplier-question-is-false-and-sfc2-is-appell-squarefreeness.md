---
id: THM-3021
title: "Settling the Hadamard multiplier question: it is FALSE as posed, and SFC(2) is squarefreeness of an Appell sequence"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. THM-3020 (K6)
  reduced SFC(2) at p = 0 to: can a positive-coefficient A = sum a_j lam^j
  share a root with its Hadamard product B = sum a_j w_j lam^j, w_j =
  prod_{i=1}^s (sj+i)? This file settles that question.
  (H1) EXPLICIT ROOT CRITERION. B = s!A + V(theta)A with V = W - W(0), so
  gcd(A,B) = gcd(A, V(theta)A), and at a root alpha of A,
  B(alpha) = sum_{i=1}^s c_i alpha^i A^{(i)}(alpha) with c_i = Delta^i W(0)/i!
  STRICTLY POSITIVE for i = 1..s. For s = 1 this is B(alpha) = alpha A'(alpha),
  giving gcd(A,B) = gcd(A,A'): SFC(2) at s = 1 IS squarefreeness.
  (H2) THE VERDICT: the Hadamard statement is FALSE as posed. For s = 2,
  A = (lam+5)(lam+9) = 45 + 14 lam + lam^2 has positive coefficients and
  B = 30 lam^2 + 168 lam + 90 = 6(lam+5)(5 lam+3) shares the root lam = -5.
  So positivity of the coefficients is NOT the operative hypothesis and no
  Polya-Schur / multiplier-sequence argument can prove SFC(2).
  (H3) WHAT IS operative: a_j = C(m,j) mu_{sj} with mu_n = n! a STIELTJES
  MOMENT sequence, i.e. A(lam) = int (1+lam u)^m dnu(u) for a positive measure
  nu. The counterexample's implied moments (45,7,1) have Hankel determinant
  -4 < 0; ours have +20 > 0. This certifies that the log-convexity hypothesis
  in opus's window-0 theorem is NOT removable.
  (H4) MAX-MODULUS THEOREM (PROVED, all m): for s <= 2 no root of A of
  MAXIMAL MODULUS is a root of B, via Re 1/(1-w) >= 1/2 for |w| <= 1, w != 1;
  general s needs only |arg zeta_k| < pi/(2(s-1)). Exact trace identity
  sum_k B(alpha_k)/(alpha_k A'(alpha_k)) = m(c_1 + c_2(m-1)) > 0 from the
  pairing 1/(1-r) + 1/(1-1/r) = 1.
  (H5) APPELL REFORMULATION. Phi_n(z) = int (z-u)^n dnu(u) =
  sum_j C(n,j)(-1)^j (sj)! z^{n-j} satisfies Phi_n' = n Phi_{n-1}, so a common
  root of (I_m, I_{m+1}) is exactly a MULTIPLE root of Phi_{m+1}:
  SFC(2) at window k (p = 0)  <=>  Phi_{k+2} is SQUAREFREE.
  Census: 250 cells, s <= 10, m <= 25 (windows to k = 24), zero failures.
  Scope: the Hadamard question is SETTLED (negatively, as posed) and replaced
  by Appell squarefreeness; SFC(2) for s >= 2 remains OPEN in that form.
source: death-star-2026-07-31-coinC2
depends_on:
  - THM-3020
  - THM-3019
related:
  - THM-2836
  - THM-3018
external:
  - "Edo, van den Essen: the Strong Factorial Conjecture."
  - "Polya, Schur: Uber zwei Arten von Faktorenfolgen (multiplier sequences)."
  - "Appell: sur une classe de polynomes."
script: 04-computation/sfc_hadamard_verdict_appell_thm3021.py
output: 05-knowledge/results/sfc_hadamard_verdict_appell_thm3021.out
script2: 04-computation/sfc_hadamard_census_thm3021.py
output2: 05-knowledge/results/sfc_hadamard_census_thm3021.out
---

# THM-3021 -- the Hadamard multiplier question, settled

THM-3020 (K6) reduced `SFC(2)` at `p = 0`, window `k`, to a resultant of two
polynomials of the **same** degree `m = k+1` differing only by a diagonal
multiplier:

```text
A(lam) = sum_j a_j lam^j,        a_j = C(m,j) (sj)!  > 0,
B(lam) = sum_j a_j w_j lam^j,    w_j = prod_{i=1}^s (sj+i) = (sj+s)!/(sj)!,
SFC(2) at window k   <=>   Res(A,B) != 0,   m = k+1.
```

and closed with the question: *can a positive-coefficient polynomial and its
Hadamard product with the increasing multiplier `w_j` share a root?* This file
answers it.

## 1. The Euclidean step and the root criterion (H1)

Write `W(X) = prod_{i=1}^{s}(sX+i)`, so `B = W(theta)A` with `theta = lam
d/dlam`, and `W(0) = s!`. Put `V = W - W(0)`, so `V(0) = 0` and `deg V = s`.
Then `B - s! A = V(theta)A`, and since `A(0) = 1 != 0` the constant `lam` is a
unit modulo `A`:

```text
gcd(A,B) = gcd(A, V(theta)A).                                          (H1a)
```

Expanding `theta^k = sum_i S(k,i) lam^i d^i` (Stirling, second kind) and
killing the `A` terms at a root, `V(X) = sum_i c_i X(X-1)...(X-i+1)` gives

```text
B(alpha) = sum_{i=1}^{s} c_i alpha^i A^{(i)}(alpha)   whenever A(alpha) = 0,
c_i = Delta^i W(0) / i!   (Newton forward differences).                (H1b)
```

**Every `c_i` is strictly positive.** `W` has positive coefficients, and
`Delta^i (X^k)(0) = i! S(k,i) >= 0` with `S(s,i) > 0` for `1 <= i <= s`; the
degree-`s` coefficient of `W` is `s^s > 0`, so `Delta^i W(0) > 0`. Computed
values: `s=1: (1)`, `s=2: (10,4)`, `s=3: (114,135,27)`, `s=4:
(1656,4272,2176,256)`, `s=5: (30120,150000,145000,40625,3125)`.

**`s = 1` is immediate.** There `B(alpha) = alpha A'(alpha)`, so (H1a) reads
`gcd(A,B) = gcd(A,A')`: SFC(2) at `s = 1` *is* squarefreeness of `A`, which
holds because `A` is a truncated exponential (`e_{m+1} = e_m + x^{m+1}/(m+1)!`).
This re-derives THM-3019 (S4) / THM-3020 (K5) in one line.

## 2. The verdict: the Hadamard statement is FALSE as posed (H2)

Take `s = 2`, so `w = (2,12,30)`, and

```text
A = 45 + 14 lam + lam^2 = (lam+5)(lam+9),        coefficients ALL POSITIVE,
B = 90 + 168 lam + 30 lam^2 = 6(lam+5)(5 lam+3),
gcd(A,B) = lam + 5.
```

A positive-coefficient polynomial **can** share a root with its `w`-Hadamard
product. Hence:

**No argument that uses only positivity of the coefficients -- in particular
no Polya-Schur / multiplier-sequence argument, and no appeal to `w_j` being
increasing or to `{w_j}` being a multiplier sequence -- can prove SFC(2).**
The question as posed in THM-3020 section 6 is settled, negatively.

(That `{w_j}` *is* a multiplier sequence is true and irrelevant: `W(X)` has
only real negative zeros `-i/s`, so by Polya-Schur `{W(j)}` maps real-rooted
polynomials to real-rooted polynomials. The counterexample above is real
rooted, and the transform is real rooted, and they share a root.)

## 3. What the operative hypothesis really is (H3)

The `a_j` are not merely positive: they are **moments**. Since
`(sj)! = int_0^inf t^{sj} e^{-t} dt`, pushing `e^{-t}dt` forward along
`u = t^s` gives a positive measure `dnu(u) = (1/s) u^{1/s-1} e^{-u^{1/s}} du`
on `(0,inf)` with

```text
a_j = C(m,j) int u^j dnu(u),    A(lam) = int (1+lam u)^m dnu(u),
                                B(lam) = int u (1+lam u)^m dnu(u).      (H3)
```

The counterexample fails exactly here: its implied moment vector `(45,7,1)`
has Hankel determinant `mu_0 mu_2 - mu_1^2 = 45 - 49 = -4 < 0`, so it is not a
Stieltjes moment sequence; ours has `1*24 - 2^2 = +20 > 0`.

This is the same mechanism as opus's window-`0` theorem of 2026-07-31
(`L(f)` and `L(f^2)` cannot both vanish, by coordinatewise log-convexity of
`Gamma` plus AM-GM). **Section 2 certifies that that hypothesis is not
removable**: drop Hankel positivity and keep only positivity of coefficients,
and the conclusion fails at the very next window.

## 4. Max-modulus theorem and a trace identity (H4)

For a root `alpha` of `A` (assumed squarefree; verified squarefree for
`s <= 5`, `m <= 8`), divide (H1b) by `alpha A'(alpha) != 0` and use
`A^{(i)}(alpha)/A'(alpha) = i! e_{i-1}(z_2,...,z_m)`, `z_k = 1/(alpha-alpha_k)`.
With `zeta_k = alpha z_k = 1/(1 - alpha_k/alpha)`,

```text
B(alpha) = 0   <=>   sum_{i=1}^{s} c_i i! e_{i-1}(zeta_2,...,zeta_m) = 0. (H4a)
```

**Lemma.** If `|w| <= 1` and `w != 1` then `Re 1/(1-w) >= 1/2`. *Proof.*
`Re 1/(1-w) = (1-Re w)/|1-w|^2` and `|1-w|^2 = 1 - 2 Re w + |w|^2 <=
2(1 - Re w)`. QED

**Theorem (max modulus).** For `s <= 2`, no root of `A` of maximal modulus is
a root of `B`. *Proof.* For `s = 1` (H4a) is `c_1 = 0`, false. For `s = 2` it
is `c_1 + 2c_2 sum_k zeta_k = 0`; at a maximal-modulus root every
`|alpha_k/alpha| <= 1` with `alpha_k != alpha`, so by the Lemma
`Re sum zeta_k >= (m-1)/2 >= 0`, whence the left side has real part
`>= c_1 > 0`. QED

More generally, if every `|arg zeta_k| < pi/(2(s-1))` then every
`e_{i-1}(zeta)` with `i <= s` has argument of modulus `< pi/2`, so the sum in
(H4a) has strictly positive real part and `B(alpha) != 0`. For `s = 2` that
condition is just `Re zeta_k > 0`, which maximal modulus supplies; for
`s >= 3` it is strictly stronger than the unit disk gives, which is exactly
where the argument stops.

**Trace identity.** For any `r != 0`, `1/(1-r) + 1/(1-1/r) = 1`: each
unordered pair of roots contributes exactly `1`, so
`sum_k sum_{j != k} 1/(1-alpha_j/alpha_k) = C(m,2)` and, for `s = 2`,

```text
sum_k B(alpha_k)/(alpha_k A'(alpha_k)) = m (c_1 + c_2 (m-1)) = 2m(2m+3) > 0.
```

Verified exactly for `m = 2..7` (28, 54, 88, 130, 180, 238). The criterion can
therefore never vanish *on average* -- only individually.

## 5. The right reformulation: Appell squarefreeness (H5)

Substituting `lam = -1/z` in (H3) and clearing `z^{-m}` gives the Appell
sequence of the measure `nu`:

```text
Phi_n(z) = int (z-u)^n dnu(u) = sum_{j} C(n,j) (-1)^j (sj)! z^{n-j},
Phi_n'(z) = n Phi_{n-1}(z),      Phi_n(0) = (-1)^n (sn)! != 0.          (H5a)
```

The nonzero roots of `I_m` correspond bijectively to the roots of `Phi_m` via
`lam = -1/z`. Now `Phi_{m+1}(z_0) = 0` **and** `Phi_{m+1}'(z_0) =
(m+1)Phi_m(z_0) = 0` is precisely the SFC(2) failure condition. Hence

```text
SFC(2) at window k, support {0,s}   <=>   Phi_{k+2} is SQUAREFREE.      (H5b)
```

Verified for `s = 1..4`, `n = 1..8`: `Phi_n' = n Phi_{n-1}`, every `Phi_n` is
squarefree, and squarefreeness agrees cell-by-cell with `gcd(I_m,I_{m+1}) = 1`.

This is the correct home for the problem: it removes the auxiliary polynomial
`B` entirely and makes the whole `p = 0` family one statement about one
sequence. It also explains the difficulty gradient. The Appell generating
function is

```text
sum_n Phi_n(z) xi^n / n! = e^{xi z} N(xi),   N(xi) = int_0^inf e^{-xi t^s - t} dt.
```

For `s = 1`, `N(xi) = 1/(1+xi)` is **meromorphic** with radius of convergence
`1`, and the classical entire-function machinery (Laguerre-Polya,
Hermite-Biehler, interlacing) applies -- which is why `s = 1` is elementary.
For `s >= 2` the Taylor coefficients of `N` are `(-1)^n (sn)!/n!`, so the
radius of convergence is **zero**: `N` is a Gevrey-`s` divergent series,
Borel-summable but not analytic at the origin, and none of that machinery is
available. **That is the precise structural reason `s >= 2` is hard**, and it
is invisible in the Hadamard formulation.

## 6. Census

`Res(A,B) != 0` certified modulo large primes (nonzero mod one prime is a
sound one-way certificate) for

```text
1 <= s <= 10,   1 <= m <= 25      -- 250 cells, 0 failures,
```

i.e. windows `k <= 24`, extending THM-3019's SFC(2) box (`s <= 8`, `k <= 8`)
in both directions.

## 7. Scope

Sections 1, 2, 4, 5 are proofs; section 6 is a finite census proving nothing
outside its box. The Hadamard multiplier question of THM-3020 section 6 is
**settled**: as posed it is false (section 2), so that route is closed, and
the surviving content is (H5b), squarefreeness of the Appell sequence
`Phi_n`. `SFC(2)` for `s >= 2` remains **OPEN** in that form, with the
max-modulus theorem of section 4 as the proved fragment and the Gevrey
obstruction of section 5 as the explanation of what a proof must overcome.
The `p >= 1` families are untouched here and remain open as in THM-3020.

Referee: `sfc_hadamard_verdict_appell_thm3021.py` (H1 criterion and `c_i > 0`
for `s <= 5`; H2 counterexample exactly; H4 lemma numerics and the trace
identity for `m <= 7`; H5 Appell identities, squarefreeness and cell-by-cell
agreement for `s <= 4`, `n <= 8`) and `sfc_hadamard_census_thm3021.py`
(250 cells).
