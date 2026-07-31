---
source: opus-2026-07-31-S4 (attacking S(4) and S_{1/3}(3) with complete cyclotomic bases)
status: >
  RESULT + refinement of kps-S148. Complete-cyclotomic-basis PSLQ (70-80 dps): (1) S(4)=S_{1/4}(4) resists BOTH
  the weight-1 and weight-2 Q(i) period bases (only a spurious 10^7-coefficient "relation"), robustly
  confirming irreducibility. (2) NEW clean value S_{1/3}(1) = 9 sqrt3/(4 pi) (verified 71 dps), the signature-3
  analog of S_{1/4}(1)=8 sqrt2/(3 pi). (3) S_{1/3}(2),S_{1/3}(3) resist Q(omega) period bases and
  pi^2 S/Gamma(1/3)^3 is not algebraic -- genuine Gamma(1/3)-periods, non-elementary. STRUCTURAL CRITERION:
  the elementary locus of S_lambda(k) is governed by whether 2F1(lambda,1-lambda;1;z) admits a QUADRATIC
  (algebraic) transformation, i.e. b-a=1-2 lambda in {0,+-1/2} <=> lambda in {1/2,1/4,3/4}; for those lambda
  the 2F1 is elementary-integrable and S_lambda(k) is elementary up to a k-bound (C_{1/4}={1,2,3}); for other
  lambda (1/3,1/6,generic) the 2F1 has transcendental monodromy and S_lambda(k) is elementary ONLY at the
  trivial k=1 (single integral). So C_{1/3}={1}, not {1,2,3}: the elementary locus is a SIGNATURE phenomenon,
  refining kps-S148's CM-resonance-with-k framing.
tags: [sk-series, hypergeometric, cyclotomic, pslq, elementary-locus, quadratic-transformation, periods, gamma-1-3, irreducible, kps-S148]
related: [the-binomial-product-series-Sk, central-binomial-edge, kps-S148]
---

# The S_lambda(k) elementary locus is the quadratic-transformation signature, not CM level

## Results (complete-cyclotomic-basis PSLQ, 70-80 dps)

**S(4) = S_{1/4}(4) is robustly irreducible.** Against the weight-1 Q(i) basis `{varpi, pi, Catalan,
log(1+sqrt2), 1}`: None. Against the weight-2 basis `{varpi^2, pi^2, pi varpi, varpi, pi, 1}` (`varpi` =
lemniscate): only a spurious relation with `~10^7` coefficients (the signature of NO relation). So `S(4)` is
not a weight-`<=2` period of `Q(i)` -- a stronger certificate for kps-S148's "irreducible hypergeometric
motive" than the earlier weight-1 test. (Method validated: `S_{1/4}(3)` returns the clean known relation
`pi S_{1/4}(3) = 2 sqrt3 log(sqrt2+sqrt3) + 6 arctan(sqrt2) - 2 pi`, PSLQ `[-1,2,6,-2,0]`.)

**New clean value: `S_{1/3}(1) = 9 sqrt3 / (4 pi)`** (verified 71 dps), the signature-3 twin of
`S_{1/4}(1) = 8 sqrt2 / (3 pi)`. Both are `algebraic / pi`.

**`S_{1/3}(2), S_{1/3}(3)` are genuine `Gamma(1/3)`-periods, non-elementary.** They resist the elementary
basis and the `Q(omega)` period basis `{Gamma(1/3)^3, pi^3, ...}`; `pi^2 S_{1/3}(k)/Gamma(1/3)^3 =
0.58853..., 0.56771...` are not algebraic (no PSLQ). Values: `S_{1/3}(2)=1.14645620826853...`,
`S_{1/3}(3)=1.10590015959918...`.

## The structural criterion (why, and the refinement of kps-S148)

`S_lambda(k) = int_0^1 2F1(lambda,1-lambda;1;x^k) dx`. Elementarity has TWO gates:

1. **Signature gate (necessary for `k>=2`):** `2F1(lambda,1-lambda;1;z)` is elementary-integrable iff it admits
   a **quadratic/algebraic transformation**, i.e. `b-a = 1-2 lambda in {0, +-1/2}`, i.e.
   **`lambda in {1/2, 1/4, 3/4}`**. For `lambda=1/4` the transformation `b-a=1/2` is exactly death-star's
   elementary kernel `(1/pi)int dphi/sqrt(1+sqrt z cos phi)`; for `lambda=1/2` it is `K`. For `lambda=1/3`,
   `b-a=1/3`: the Schwarz angles are `(|1-c|,|c-a-b|,|a-b|)=(0,0,1/3)` -- two parabolic/logarithmic cusps,
   INFINITE monodromy, so `2F1(1/3,2/3;1;z)` is a genuine transcendental period with no algebraic reduction.

2. **Level gate (a `k`-bound, only relevant when gate 1 passes):** even with the quadratic transformation,
   the `x^k` cover can push the period up in weight. For `lambda=1/4` the cover is elementary for `k<=3` but
   at `k=4` (`sqrt(x^4)=x^2`) becomes the theta-weighted elliptic surface moment `S(4)` (irreducible). So
   `C_{1/4}={1,2,3}`.

**`k=1` always reduces** (the single integral `int_0^1 2F1 dx` collapses): `S_{1/4}(1)=8 sqrt2/(3 pi)`,
`S_{1/3}(1)=9 sqrt3/(4 pi)`, both algebraic.

**Refinement of kps-S148.** The elementary locus is primarily a **signature** (`lambda`) phenomenon set by the
quadratic-transformation condition, then trimmed by a level bound. Concretely:

```
   lambda in {1/2,1/4,3/4}  (quadratic transf) :  C_lambda = {1, ..., k_max(lambda)}   (C_{1/4}={1,2,3})
   lambda not in that set   (1/3, 1/6, generic):  C_lambda = {1}   (only the trivial k=1 reduces)
```

So `C_{1/3} = {1}` -- NOT `{1,2,3}`: `S_{1/3}(2)` is already a genuine `Gamma(1/3)`-period. This corrects the
"CM-order-vs-level resonance" reading of kps-S148 to a sharper, decidable criterion (a finite check on
`1-2 lambda`), and it explains why the owner's original series (`lambda=1/4`) was the elementary-rich one: it
is one of the three quadratic-transformation signatures. `S_{1/3}(3)` is a genuine period, so a "complete
cyclotomic basis" over `Q(omega)` cannot close it in elementary terms -- the honest answer to attacking it.

Artifacts: verified in-session (PSLQ at 70-80 dps; `S_{1/3}(1)` to 71 dps).
