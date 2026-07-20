---
id: THM-1540
title: "GMC(2) REDUCED TO A NULLCONE STRUCTURE THEOREM — two of three pieces proved, and the remaining gap named exactly. NOT FINISHED. (0) THE REDUCTION: for one complex Gaussian, if the nullcone N = {P : E[P^m] = 0 ∀m} equals {P with all charges > 0} ∪ {all charges < 0} ∪ {0}, then GMC(2) follows in one line, since for one-sided P every monomial of QP^m has charge of modulus ≥ m|k| − deg Q, non-zero once m > deg(Q)/|k|. (I) L1 PROVED: a purely charge-0 P is g(|Z|²) with |Z|² ~ Exp(1), and E[P^m] = ∫g^m e^{−s}ds ~ c^m(md)! for deg g = d ≥ 1 — the top term grows factorially faster than every correction, so it cannot cancel; hence g is constant and then 0. (II) L2 PROVED (modulo the classical one-variable constant-term fact): writing P_top = y^d h(x/y), the charge-0 coefficient of (P_top)^m is [t^{md/2}]h(t)^m and its contribution to E[P^m] carries (md/2)!, which dominates every lower-degree term — so vanishing for all m forces h's support off d/2, i.e. THE TOP-DEGREE PART MUST BE ONE-SIDED. (III) THE GAP, stated plainly: the induction from 'P_top one-sided' to 'P one-sided' is NOT done. If P_top is one-sided while lower parts straddle, the leading asymptotics vanish and the argument must descend a degree, and that descent is unwritten. (IV) EVIDENCE: exhaustive search finds NO mixed-charge nullcone element — 16 at degree 2 and 160 at degree 3 with coefficients in {−1,0,1}, and 48 at degree 2 with coefficients in {−2..2}, every one of them one-sided"
status: >
  NOT FINISHED — GMC(2) is not proved here and must not be cited as such.
  (0) The reduction is PROVED and is one line.
  (I) PROVED, by factorial dominance; checked on sample g.
  (II) PROVED modulo the one-variable statement "CT(f^m) = 0 for all m implies supp(f)
  lies strictly on one side of 0" for Laurent f in one variable, which is classical
  (the n = 1 case of Duistermaat-van der Kallen) and is cited, not reproved.  The
  dominance step — that the (md/2)! term outgrows all lower-degree contributions — is
  an asymptotic argument and is stated as such.
  (III) OPEN.  This is the whole remaining content and it is not small: it is exactly
  where a mixed-charge counterexample would have to hide.
  (IV) SEARCH, with the detection floor stated: coefficients in {-1,0,1} over degree
  <= 3 and in {-2..2} over degree <= 2, depth m <= 9.  A search bounds a counterexample;
  it does not prove a theorem, and after the m <= 3 sieve episode of S128c117 the depth
  is stated rather than assumed adequate.
source: kind-pasteur-2026-07-20-S128c118 (owner: finish GMC(2) by finishing the stronger 2-dimensional nullcone conjecture)
depends_on:
  - THM-1525    # GMC(4) false; the mechanism that makes n = 4 work and why it needs 4
related: []
external:
  - "Duistermaat & van der Kallen (1998): constant terms of powers of Laurent polynomials."
script: 04-computation/gmc2_nullcone_kps_S128c118.py, gmc2_structure_kps_S128c118.py (+ .out)
---

# THM-1540 — GMC(2) reduced, two thirds proved, one third open

**Headline, so it is not mis-cited: GMC(2) is NOT proved here.** What is done is a
reduction to a clean structure theorem, a proof of two of its three pieces, and an exact
statement of what is left.

## 0. The reduction

Two real Gaussians = one complex Gaussian `Z`, with `E[Z^a conj(Z)^b] = δ_{ab} a!`. Call
`a − b` the **charge** of `Z^a conj(Z)^b`; it is the Fourier mode in the argument of `Z`,
and `E` kills every non-zero charge.

> **Structure conjecture.** `N := {P : E[P^m] = 0 ∀ m ≥ 1}` equals
> `{P : all charges > 0} ∪ {P : all charges < 0} ∪ {0}`.

**This gives GMC(2) immediately.** If every charge of `P` is `≥ k > 0`, then every
monomial of `Q P^m` has charge at least `mk − deg Q`, which is positive — so the
expectation vanishes — as soon as `m > deg(Q)/k`. Symmetrically for negative charges. So
the Mathieu property holds with an explicit threshold.

That is why the structure statement is the right target: it is stronger, and it settles
GMC(2) in one line.

## I. L1 — a purely charge-0 element is zero (PROVED)

Charge-0 monomials are `Z^a conj(Z)^a = s^a` with `s = |Z|² ~ Exp(1)`, so a charge-0 `P`
is `g(s)` and

> `E[P^m] = ∫₀^∞ g(s)^m e^{−s} ds`.

If `deg g = d ≥ 1` with leading coefficient `c`, the integrand is dominated by
`c^m s^{md}` and `E[P^m] ~ c^m (md)!`. Every correction is smaller by a factor that is
polynomial in `m` against a factorial, so **no cancellation is possible for large `m`**.
Hence `g` is constant, and then `E[P^m] = g^m = 0` forces `g = 0`. ∎

## II. L2 — the top-degree part must be one-sided (PROVED, modulo a classical fact)

Let `d = deg P` and `P_top` its degree-`d` homogeneous part. Writing `x = Z`,
`y = conj(Z)`, `P_top = y^d h(x/y)` for a polynomial `h` of degree `≤ d`. Then

> charge-0 coefficient of `(P_top)^m` = `[t^{md/2}] h(t)^m`,

and that monomial is `Z^A conj(Z)^A` with `A = md/2`, contributing with weight `(md/2)!`
— which **dominates every lower-degree contribution to `E[P^m]`**, since lower degree
means a strictly smaller factorial. So the leading asymptotics of `E[P^m]` is

> `(md/2)! · CT( (h(t)·t^{−d/2})^m )`.

For this to vanish for all `m`, the Laurent polynomial `h(t)t^{−d/2}` must have support
strictly on one side of `0` — the one-variable constant-term fact (the `n = 1` case of
Duistermaat–van der Kallen), cited here rather than reproved. Support strictly on one
side of `d/2` says exactly that **`P_top` is one-sided in charge**. ∎

Verified on samples: `h = 1+t`, `1−t`, `1+t²` (straddling) give central coefficients
`2, 6, 20, …` — never all zero; `h = t`, `1`, `t²` (one-sided) give identically zero.

## III. The gap, named exactly

**The induction from "`P_top` one-sided" to "`P` one-sided" is not done.** If `P_top` is
one-sided while a lower-degree part straddles, the leading term vanishes identically and
the argument has to descend to the next order in the degree filtration. That descent is
where a mixed-charge nullcone element would have to live, and it is unwritten.

This is not a formality. The GMC(4) counterexample (THM-1525) works precisely by having
*two* charge families conspire — `(1+Z)` contributing charges `{0,+1}` against
`(conj(Z) − Y)` contributing `{−1, 0}` — so mixed charge is exactly the mechanism that
does produce counterexamples one dimension up. The reason it cannot be imported to `n = 2`
is that it needs an independent `Exp(1)`, which costs the second complex variable
(THM-1525 §IV). But "this particular mechanism does not fit" is not "no mechanism fits".

## IV. Evidence, with its detection floor

Exhaustive search over `N`:

| box | nullcone elements | one-sided | **mixed** |
|---|---|---|---|
| coeffs `{−1,0,1}`, deg ≤ 2 | 16 | 16 | **0** |
| coeffs `{−1,0,1}`, deg ≤ 3 | 160 | 160 | **0** |
| coeffs `{−2,…,2}`, deg ≤ 2 | 48 | 48 | **0** |

all at depth `m ≤ 9`. **No mixed-charge element exists in any box searched.**

The floor is stated deliberately. In S128c117 a sieve at `m ≤ 3` produced three
"counterexamples" that dissolved on proper testing; the fix there was depth, and the
lesson is that a search of this shape earns belief only with its box and depth attached.

## Named next

- **Close the induction in III.** The natural route is a degree filtration: assume the
  top `j` graded pieces are one-sided and show the `(j+1)`-st must be, by the same
  factorial-dominance argument applied to the leading *surviving* order.
- Widen (IV) to degree 4 and coefficients `{−2..2}` — the box is large but the early
  filter `E[P] = E[P²] = 0` prunes hard.
- The one-variable constant-term fact is used as a citation in II; for a self-contained
  write-up it should be proved directly, which in one variable is short.
