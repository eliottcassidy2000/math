---
id: THM-3039
title: "The FC(n) exponential-period bridge: the forced level is the simplex volume 1/(n-1)!, and the FC(2) proof's crux step is one-dimensional"
status: >
  PROVED + VERIFIED-EXACT. Generalises THM-3031 from two variables to all n, and locates
  precisely where the (owner-supplied, unverified) FC(2) proof strategy fails for n = 3.
  (L1) FORCED LEVEL = SIMPLEX VOLUME. By THM-3018, for homogeneous f of degree d in n
  variables L(f^m) = (dm+n-1)! int_{Delta_{n-1}} g^m dV. A counterexample makes every
  m >= 1 moment vanish, so int_{Delta_{n-1}} e^g dV = Vol(Delta_{n-1}) = 1/(n-1)!.
  The value 1 of THM-3031 was NOT special: it is Vol(Delta_1) = 1. For n = 3 the level is
  Vol(Delta_2) = 1/2.
  (L2) THE BRIDGE, ALL n:  int_{Delta_{n-1}} e^{Q} dV != 1/(n-1)! for every nonconstant
  Q in Qbar[u_1..u_{n-1}]  ==>  FC(n) for homogeneous f. Since 1/(n-1)! is RATIONAL,
  transcendence of the period still suffices, and nonvanishing still does not.
  (L3) THE CRUX STEP OF THE FC(2) PROOF IS ONE-DIMENSIONAL. That proof excludes nonconstant
  projective leading forms by constellation monodromy (Pakovich-Muzychuk, positive-interval),
  a theory of the polynomial moment problem ON AN INTERVAL, and finishes with transcendence of
  a ONE-dimensional exponential integral int_A^B e^q dx (E-functions, Siegel-Shidlovskii,
  Beukers). For n = 3 the domain is the 2-simplex, a TRIANGLE, and the required period
  int_{Delta_2} e^Q dA is a 2-dimensional period. Siegel-Shidlovskii/Beukers is a theory of
  solutions of linear ODEs, i.e. one variable. So both crux steps are dimensionally specific,
  and n = 3 is not a matter of pushing the same argument further.
  (L4) THE CONCRETE ROUTE THAT REMAINS is iterated reduction: int_{Delta_2} e^Q dA =
  int_0^1 [ int_0^{1-u} e^{Q(u,v)} dv ] du, an outer integral of an inner exponential integral
  with u-DEPENDENT polynomial and ALGEBRAIC-in-u endpoints. Whether that lands in a known
  transcendence framework is the sharp open question for FC(3).
  Verified: Vol(Delta_2) = 1/2; the model h(u,v) = exp(2 pi i (2u - u^2)), which has all
  m >= 1 moments zero because T = 1-(1-u)^2 is uniform on [0,1] under the triangle's marginal
  2(1-u), gives int_Delta e^h dA = 1/2 to 6.8e-37. THM-3018 section 4's S_3 formula
  int_T g^{3k} = 1/((3k+1)(3k+2)) independently reconfirmed (m=3 -> 1/20, m=6 -> 1/56).
source: death-star-2026-08-01-coinC2
depends_on:
  - THM-3018
  - THM-3031
related:
  - THM-3021
  - kps-S166
external:
  - "octonion/mathematics/fc: draft proof of general FC(2) (NOT peer-reviewed, NOT verified here)."
  - "Pakovich, Muzychuk: solution of the polynomial moment problem."
  - "Siegel, Shidlovskii; Beukers: E-functions."
script: 04-computation/fc3_bridge_simplex_volume_thm3039.py
output: 05-knowledge/results/fc3_bridge_simplex_volume_thm3039.out
---

# THM-3039 -- the forced level is the simplex volume

## 1. The level (L1)

THM-3018 reduces the homogeneous Factorial Conjecture to a simplex moment problem:
for `f` homogeneous of degree `d` in `n` variables and `g = f|_Delta`,

```text
L(f^m) = (dm+n-1)! * int_{Delta_{n-1}} g^m dV,
```

with `Delta_{n-1} = {u_i >= 0, sum u_i <= 1}` and `dV` the coordinate volume, so
`Vol(Delta_{n-1}) = 1/(n-1)!`. If `g != 0` has all moments `m >= 1` vanishing then

```text
int_{Delta_{n-1}} e^{g} dV = sum_{m>=0} (1/m!) int g^m dV = Vol(Delta_{n-1}) + 0
                           = 1/(n-1)!.                                        (L1)
```

**The `1` in THM-3031 was not a special constant** -- it is `Vol(Delta_1) = 1`, the
length of `[0,1]`. For `n = 3` the level is `Vol(Delta_2) = 1/2`. (This is exactly the
normalisation trap flagged when the FC(3) analogue was first proposed: the simplex
volume is *not* `1` in general.)

*Verification.* `Vol(Delta_2) = 0.5`. For a model with all `m >= 1` moments vanishing,
note that under the uniform measure on `Delta_2` the marginal density of `u` is
`2(1-u)`, so `T = 1-(1-u)^2` is uniform on `[0,1]`; hence
`h(u,v) = exp(2 pi i (2u - u^2))` has `int_Delta h^m dA = Vol * int_0^1 e^{2 pi i m t}dt
= 0` for every `m >= 1`. Computed moments `m = 1..5` are `< 6.5e-35`, and

```text
int_{Delta_2} e^{h} dA = 0.5 + O(6.8e-37).
```

## 2. The bridge for every n (L2)

```text
If  int_{Delta_{n-1}} e^{Q} dV != 1/(n-1)!   for every nonconstant Q in Qbar[u_1..u_{n-1}],
then FC(n) holds for homogeneous f.                                            (L2)
```

*Proof.* As in THM-3031: a counterexample is nonconstant (the `m = 1` condition kills
constants), may be taken with algebraic coefficients (the moment conditions are
polynomial in the coefficients with rational coefficients, since
`int_{Delta} u^alpha dV` is rational, so the locus is a `Q`-variety with Zariski-dense
`Qbar`-points), and by (L1) forces the period to `1/(n-1)!`. QED

Since `1/(n-1)!` is **rational**, transcendence of the period suffices and mere
nonvanishing does not -- THM-3031 (B3) survives verbatim in every dimension.

## 3. Why n = 3 is not "more of the same" (L3)

The FC(2) draft's two crux steps are **dimensionally specific**:

* **Constellation monodromy** (Pakovich-Muzychuk, positive-interval) excludes nonconstant
  projective leading forms. It is a theory of the polynomial moment problem **on an
  interval** -- one-dimensional.
* **E-function transcendence** of `int_A^B e^{q(x)} dx` closes the flat-top residual, via
  Siegel-Shidlovskii and Beukers. That is a theory of **solutions of linear ODEs**, i.e.
  functions of one variable.

For `n = 3` the domain is the 2-simplex and the object in (L2) is a **2-dimensional
period** `int_{Delta_2} e^Q dA`, `Q in Qbar[u,v]`. Neither crux step applies as stated.
**So FC(3) is not blocked by needing to push the FC(2) argument further; it is blocked by
both of its essential tools being one-variable theories.** Naming that is the point of
this section: it redirects effort away from "generalise the constellation argument".

## 4. The route that does remain (L4)

Iterating the integral keeps one variable at a time:

```text
int_{Delta_2} e^{Q} dA = int_0^1 [ int_0^{1-u} e^{Q(u,v)} dv ] du.
```

The inner integral is an exponential integral in `v` whose polynomial `Q(u, .)` has
coefficients algebraic **in `u`** and whose endpoints `0, 1-u` are algebraic in `u`. So
it is exactly an `int_A^B e^q` of the FC(2) type, but **with algebraic function
coefficients rather than algebraic constants**. Whether transcendence survives that
relative setting -- and whether the outer `u`-integration can then be controlled -- is
the sharp open question for FC(3). This is the concrete handoff.

## 5. Scope

(L1), (L2) are proofs; the verification is numerical at `dps = 30` and the identity is
formal. (L3) is a statement about the **hypotheses** of the cited theories, not a proof
that no `n = 3` argument exists -- it says the *stated* tools are one-dimensional, and
should not be quoted as an impossibility. **Nothing here proves FC(3), and nothing here
verifies the external FC(2) draft**, whose provenance is unchecked.

Note also that (L2) is an implication only: `int e^Q != Vol` is strictly stronger than
`FC(n)`, since a counterexample needs *all* moments to vanish while the period condition
is a single equation.

Independently reconfirmed in passing: THM-3018 section 4's `S_3` mechanism, whose surviving
moments `int_T g^{3k} dA = 1/((3k+1)(3k+2))` were reproduced exactly (`m = 3 -> 1/20`,
`m = 6 -> 1/56`, and `k = 0` giving `1/2 = Vol`). That section is unaffected by the
audit flag on section 4b.
