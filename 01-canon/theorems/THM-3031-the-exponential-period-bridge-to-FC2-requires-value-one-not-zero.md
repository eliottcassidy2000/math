---
id: THM-3031
title: "The exponential-period bridge to FC(2): the operative constant is 1, not 0, and transcendence is what is needed"
status: >
  PROVED. The owner supplied an external claim -- int_0^1 e^{Q(t)} dt != 0 for every
  nonconstant Q in Qbar[t], and moreover that this integral is always TRANSCENDENTAL --
  together with the assertion that it implies FC(2). This file determines exactly what
  is true.
  (B1) A FC(2) counterexample forces the exponential period to the value ONE, not zero:
  by THM-3018, FC(2) is equivalent to "int_0^1 g^m du = 0 for all m >= 1 implies g = 0",
  and a counterexample g gives int_0^1 e^{g(u)} du = sum_m (1/m!) int_0^1 g^m du
  = 1 + 0 = 1. Also g is necessarily NONCONSTANT, since m = 1 kills constants.
  (B2) SPECIALISATION: the moment conditions are polynomial in the coefficients of g with
  RATIONAL coefficients (int_0^1 u^j du = 1/(j+1)), so the counterexample locus is a
  variety defined over Q; Qbar-points are Zariski dense, so a counterexample may be taken
  with ALGEBRAIC coefficients. Hence Q := g is a legitimate input to the external claim.
  (B3) THEREFORE the literal nonvanishing statement "int e^Q != 0" is INSUFFICIENT: the
  counterexample yields the value 1, which is != 0, so no contradiction arises and FC(2)
  does NOT follow. The bridge as stated does not close.
  (B4) THE CORRECT MINIMAL BRIDGE is
        int_0^1 e^{Q(t)} dt != 1   for every nonconstant Q in Qbar[t]   ==>   FC(2),
  a four-line derivation. This is strictly weaker than transcendence and strictly stronger
  than nonvanishing. In particular the TRANSCENDENCE form of the external claim DOES imply
  FC(2) (the value 1 is rational), which is why the transcendence statement -- not the
  nonvanishing statement -- is the operative one.
  Verified symbolically: FC(2) holds in degrees 1, 2, 3 (moment system has only g = 0), the
  moment conditions have all-rational coefficients, and the identity int_0^1 e^g = 1 is
  formal.
source: death-star-2026-08-01-coinC2
depends_on:
  - THM-3018
related:
  - THM-3019
  - THM-3021
  - kps-S166
external:
  - "external AI-generated claim: int_0^1 e^Q dt is transcendental for nonconstant Q in Qbar[t] (provenance unverified here)."
  - "Lindemann-Weierstrass; Siegel-Shidlovskii; Beukers, E-functions."
script: 04-computation/expq_fc2_bridge_thm3031.py
output: 05-knowledge/results/expq_fc2_bridge_thm3031.out
---

# THM-3031 -- the bridge needs `!= 1`, and transcendence delivers it

## 1. What was claimed

The owner relayed an external result: for every nonconstant `Q` with algebraic
coefficients,

```text
int_0^1 e^{Q(t)} dt != 0,        and indeed this integral is TRANSCENDENTAL,
```

with the remark that this implies `FC(2)`, and that the transcendence statement
is the more impressive half. **Both halves of that remark turn out to be
correct, but not for the reason the phrasing suggests: the nonvanishing form is
not merely less impressive, it is insufficient.**

## 2. A counterexample forces the value 1 (B1)

By THM-3018 the Factorial Conjecture in two variables is *exactly* the moment
problem on the 1-simplex:

```text
FC(2)   <=>   [ g in C[u],  int_0^1 g(u)^m du = 0 for all m >= 1   =>   g = 0 ].
```

Suppose `FC(2)` fails and `g != 0` witnesses it.

* `g` is **nonconstant**: for a constant `c` the `m = 1` condition reads
  `int_0^1 c du = c = 0`.
* The exponential period is then

```text
int_0^1 e^{g(u)} du  =  sum_{m >= 0} (1/m!) int_0^1 g(u)^m du  =  1 + 0  =  1.   (B1)
```

The series converges absolutely (`g` is bounded on `[0,1]`), so the interchange
is legitimate. **The period is forced to `1`, not to `0`.**

## 3. The counterexample may be taken algebraic (B2)

Fix a degree bound `d`. The conditions `int_0^1 g^m du = 0` are polynomial in the
coefficients `c_0,...,c_d` of `g`, with **rational** coefficients, because
`int_0^1 u^j du = 1/(j+1)`. (Verified: for `d = 3`, `m = 1` gives
`c_0 + c_1/2 + c_2/3 + c_3/4`, and all coefficients are rational at `m = 2` as
well.) By Noetherianity finitely many `m` cut out the same locus, so the
counterexample locus `V_d` is an affine variety **defined over `Q`**.

If `FC(2)` fails then some `V_d` has a `C`-point off the origin. Since `Qbar` is
algebraically closed and `V_d` is of finite type over `Q`, the `Qbar`-points are
Zariski dense in `(V_d)_C`; a component not contained in `{0}` therefore carries
`Qbar`-points `!= 0`. **So the counterexample may be taken with algebraic
coefficients**, and `Q := g` is an admissible input to the external claim.

## 4. Consequently the nonvanishing form does not close (B3)

Putting (B1) and (B2) together, a failure of `FC(2)` produces a nonconstant
`Q in Qbar[t]` with

```text
int_0^1 e^{Q(t)} dt = 1.
```

Against the hypothesis `int e^Q != 0` this is **consistent** -- `1 != 0` -- so no
contradiction arises and `FC(2)` does not follow. The literal nonvanishing
statement is too weak for the stated application.

## 5. The correct minimal bridge (B4)

```text
   If  int_0^1 e^{Q(t)} dt != 1  for every nonconstant Q in Qbar[t],  then FC(2).   (B4)
```

*Proof.* If `FC(2)` fails, (B2) supplies a nonconstant `Q in Qbar[t]` and (B1)
gives `int_0^1 e^Q = 1`, contradicting the hypothesis. QED

Three remarks.

* **Transcendence suffices.** `1` is rational, so the transcendence form of the
  external claim implies (B4) and hence `FC(2)`. This is the precise sense in
  which the transcendence statement is the operative one.
* **(B4) is much weaker than transcendence.** It only forbids a single value.
  Anyone wanting `FC(2)` from this route needs to exclude `1`, not to prove full
  transcendence -- a considerably smaller target, and a natural next objective.
* **The converse fails.** `int_0^1 e^Q = 1` is a *single* equation on `Q`, while a
  `FC(2)` counterexample needs *all* moments to vanish. So (B4) is strictly
  stronger than `FC(2)`; this is an implication, not an equivalence, and it must
  not be quoted as one.

## 6. Sanity checks

The moment system has only the zero solution in degrees `1, 2, 3` (exact
Groebner elimination), so `FC(2)` holds there and no low-degree counterexample
exists to test the bridge against -- consistent with THM-3018 R3.

## 7. Scope

(B1), (B2), (B4) are proofs. **Nothing here proves or disproves the external
claim in either form**, and its provenance is not verified in this repo; what is
established is which form of it is usable. Nothing here proves `FC(2)`. The
relation to kps-S166 is that the object identified there --
`Phi_f(t) = int e^{tf - |x|} dx` forced *constant* -- is the same phenomenon:
in both settings the counterexample pins an exponential period to a **rational
value**, and it is arithmetic rigidity of that value, not its nonvanishing, that
does the work.

Referee: `expq_fc2_bridge_thm3031.py`.
