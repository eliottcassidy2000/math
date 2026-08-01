---
id: HYP-9078
title: "The exponential-integral non-vanishing claim, its FC bridge, and what THM-3022 already settles there"
status: >
  PARTIAL + PROVENANCE NOTE. External claim considered: int_0^1 e^{Q(t)} dt != 0
  for every nonconstant Q in Qbar[t]. Established here: (i) the real-coefficient
  case is TRIVIAL by positivity, so all content is over C -- the same shape
  HYP-9076 sec 6 found for FC; (ii) DEGREE ONE follows from
  Lindemann-Weierstrass, since the integral is e^b(e^a-1)/a and vanishing needs
  a in 2 pi i Z, transcendental, so algebraicity of the coefficients is exactly
  the operative hypothesis; (iii) the precise relation to FC is that
  int_0^1 e^Q = sum_m M(Q^m)/m! with M(t^j) = 1/(j+1), so FC kills EVERY moment
  separately (weight j!) while this claim kills ONE weighted sum (weight
  1/(j+1)); (iv) THM-3022 applies directly to the [0,1] Lebesgue weight, which
  is strictly log-convex, so the FC-analogue for that measure is settled at TWO
  SLOTS. NOT established: the asserted implication "this conjecture implies
  FC(2)" -- no argument was supplied and none is reconstructed here.
  SEPARATELY: a pasted abstract claiming nonsofic groups have been constructed
  does NOT match its own cited arXiv id; see section 4.
source: opus-2026-07-31-amm12592-writeup
related:
  - THM-3022  # two-slot moment dichotomy for weight sequences
  - HYP-9076  # FC as moment determinacy; the real/complex split
  - THM-2173
external:
  - "Lindemann-Weierstrass theorem."
  - "Siegel-Shidlovskii theory of E-functions (the natural degree >= 2 tool)."
script: 04-computation/fc_exponential_integral_bridge.py
---

# HYP-9078 -- the exponential-integral claim and the FC bridge

## 1. The real case is trivial; all content is complex

If `Q` has REAL algebraic coefficients then `e^(Q(t)) > 0` on `[0,1]`, so the
integral is positive. **The entire content of the claim is for `Q` with
non-real coefficients.**

This is the same structure established for FC in HYP-9076 sec 6, where `L` is
a positive functional and so `L(f^2) > 0` for real `f != 0`, making every
real-coefficient census vacuous. Two independent problems, one shape: the
positivity of the underlying measure trivializes the real case and pushes all
difficulty into the complex one.

## 2. Degree one is Lindemann--Weierstrass

For `Q = a t + b` with `a != 0`,

```text
int_0^1 e^(at+b) dt = e^b (e^a - 1)/a,
```

which vanishes iff `e^a = 1`, i.e. `a in 2 pi i Z \ {0}`. Every such `a` is
transcendental, so no algebraic `a` works. **Algebraicity of the coefficients
is exactly the hypothesis doing the work**, and at degree one it is precisely
Lindemann--Weierstrass. At degree two the integral is an error-function value
and the natural machinery is Siegel--Shidlovskii E-function theory, which
matches the description supplied with the claim.

## 3. The bridge to FC, stated precisely

Expanding the exponential,

```text
int_0^1 e^(Q) dt  =  sum_(m>=0) M(Q^m)/m!,      M(t^j) = 1/(j+1).
```

So both problems are moment problems for a functional against a measure, and
they differ in WHAT they demand:

```text
FC:                L(f^m) = 0  for EVERY m,      weight w_j = j!
this claim:        sum_m M(Q^m)/m! = 0,          weight w_j = 1/(j+1)
```

FC kills every moment separately; the integral claim kills a single weighted
sum of them. That is the exact relationship. It makes an implication in either
direction a substantive theorem rather than a reformulation, and **the asserted
implication "the conjecture implies FC(2)" is not verified here** -- no
argument was supplied with the claim, and none is reconstructed.

## 4. What THM-3022 settles immediately

The `[0,1]` Lebesgue weight `w_j = 1/(j+1)` is strictly log-convex:
`w_j^2 = 1/(j+1)^2 < 1/(j(j+2)) = w_(j-1) w_(j+1)`, since
`(j+1)^2 > j^2+2j`. So THM-3022 applies verbatim and gives

```text
Q_w(a,b) = w_(2a) w_b^2 - 2 w_a w_b w_(a+b) + w_a^2 w_(2b) > 0
```

for all `a != b` -- computed exactly, e.g. `Q_w(0,1) = 1/12`,
`Q_w(1,3) = 11/1680`, `Q_w(2,5) = 7/3960`, positive throughout `0 <= a < b <= 9`.

**Hence the FC-analogue for the interval measure is settled at two slots by the
same log-convexity-plus-AM-GM argument** (and at THREE slots by THM-3028,
whose resultant invariant `R_w` is nonzero for `w_j = 1/(j+1)` on every tested
triple): `M(f) = M(f^2) = 0` forces `f = 0`.
The two-slot dichotomy is therefore not special to the factorial weight; it
covers the exponential-integral setting too, and THM-3022 sec 4's Fibonacci
counterexample shows the hypothesis it needs (log-convexity, or merely
`Q_w != 0`) is the real content.

## 5. Provenance note on the accompanying nonsofic claim

A second item supplied alongside was an abstract headed "NONSOFIC GROUPS
EXIST", asserting the construction of an infinite finitely presented nonsofic
group, cited to `arXiv:2604.19174`. **That id does not carry that paper.**
Fetched directly, `arXiv:2604.19174` is

```text
"On minimal non-sofic and omega-non-sofic groups", Kivanc Ersoy,
21 Apr 2026 (v1), rev. 8 Jun 2026 (v4).
Abstract opens: "We investigate structural properties of non-sofic groups,
                 ASSUMING THAT SUCH GROUPS EXIST."
```

It does not contain the phrase "nonsofic groups exist" as a claim, does not
construct such a group, and is explicitly conditional -- it proves that IF a
non-sofic group exists THEN the class of omega-non-sofic groups is non-empty.
So the cited source is evidence that the existence question was still open,
not that it was settled. The pasted abstract is recorded here as an
unverified external claim and is NOT used anywhere.

The repo has no sofic / property-(T) / expander thread to merge it into in any
case: grepping canon returns only incidental hits (THM-438 Paley/Catalan,
THM-139 chirality, THM-2228 Mahler, HYP-3832 PSL(2,7) cochain).


## 6. The claim strengthened to TRANSCENDENCE, and degree one PROVED

The claim as later supplied is stronger: `int_0^1 e^(Q(t)) dt` is
TRANSCENDENTAL for every nonconstant `Q in Qbar[t]`, not merely nonzero.

**Theorem (degree one).** For algebraic `a != 0` and algebraic `b`,
`int_0^1 e^(at+b) dt` is transcendental.

*Proof.* The integral is `(e^(a+b) - e^b)/a`. Suppose it equals an algebraic
`gamma`. Then

```text
1 * e^(a+b)  +  (-1) * e^b  +  (-a gamma) * e^0  =  0,
```

a `Qbar`-linear relation among `e^alpha` at the distinct algebraic exponents
`a+b, b, 0`. Lindemann--Weierstrass says those exponentials are linearly
independent over `Qbar`, so every coefficient vanishes -- but the first is
`1`. Contradiction. Edge cases: if `b = 0` use the exponents `a, 0`; if
`a+b = 0` the relation gives `a gamma = 1 - e^b` with `e^b` transcendental.
QED

So the transcendence form is true at degree one and, exactly as with the
non-vanishing form, **algebraicity of the coefficients is the whole
hypothesis**. Degree two is an error-function value, where Siegel--Shidlovskii
E-function theory is the natural tool, matching the description supplied.

## 6a. The log bridge -- and why it does NOT give the claimed implication

There is an exact change of variables linking the two functionals. With
`x = -log u`,

```text
L(g) = int_0^inf g(x) e^(-x) dx = int_0^1 g(-log u) du,
```

verified numerically: `int_0^1 (-log u)^j du = j!` for `j = 0..5`. **So the
factorial weight and the `[0,1]` weight are the same functional in different
coordinates.**

That is the obvious candidate route from the integral claim to FC, and it does
NOT work: under the substitution the exponent becomes a polynomial in
`log u`, not in `u`, so `L(e^(t f))` is `int_0^1 e^(t f(-log u)) du` with a
LOGARITHMIC exponent, whereas the claim concerns `int_0^1 e^(Q(u)) du` with
`Q` polynomial in `u`. These are different objects. **Whatever the asserted
implication "the conjecture implies FC(2)" is, it is not this change of
variables**, and it remains unreconstructed here.

## 7. Follow-up on the second item's source

**CORRECTION to the earlier entry.** A verifiable source does exist:
`github.com/openai/ten-proofs`, fetched directly, is a real repository of
Lean 4 formalizations whose README reads *"This repository contains Lean 4
formalizations of the results presented in Ten advances in mathematics and
theoretical computer science by OpenAI"*, with top-level files including
`NonSoficGroup.lean`, `SpherePacking.lean`, `GapCVP.lean`,
`ConnesRigidity.lean`, `Permanent.lean`, `MulticolorTriangleRamsey.lean`,
`EhrhartVolumeInequality.lean`, `QuantumParallelRepetition.lean`,
`MetricCodes.lean`, `CompactnessAndDegeneracy.lean`. My previous "no
verifiable source" was wrong and is withdrawn; the arXiv-id mismatch of
section 5 stands as a separate fact about that identifier only.

**What was checked.** Fetching `NonSoficGroup.lean` (~2900 lines) reports
**one occurrence of `sorry`**, in the proof of
`elementaryGroup_finitelyGenerated`. So as fetched the file is not
sorry-free. Caveats: this is a single automated read, not a compile; a
`sorry` in a finite-generation support lemma is not the main theorem; and
the companion page `openai.com/index/ten-advances-in-mathematics` returns
HTTP 403, so its framing (preliminary vs peer-reviewed) could not be read.
The honest status is therefore: **artifact exists, one admitted lemma
observed, correctness unverified.**

The repo still has no sofic / property-(T) / expander thread, and nothing
from this claim is used in any result here.
