# Session 127 independent peer review

**Verdict: PASS in the stated scopes. No mathematical proof gap found.**

Reviewed:

- `independent_norm/UNIFORM_NORM_TO_CF.md` and its checker;
- `inheritance_gsma/INHERITANCE_AND_KERNEL.md` and its checker.

## All-order positive norm-one implication

The proof holds for every integer `k>=3` and **positive** integers `p,q`.
Homogeneity forces primitivity. The exclusions `p<=q` and `p>2q` are
correct; the endpoint leaves precisely `(2,1)`, an actual convergent.
For the interior, the normalized polynomial
`h_k(x)=1-sum_{j=1}^k x^{-j}` has derivative strictly greater than `1/2`
on `[1,2]`. Therefore

`|alpha_k-p/q| < 2/p^k <= 2/(q+1)^3 < 1/(2q^2)`.

The last inequality is valid even at `q=1`, since its numerator is
`q^2(q-1)+3q+1>0`. Legendre's conclusion follows. The general residual
gate also correctly permits nonprimitive fractions: reducing their
denominator only weakens the required approximation threshold.

The signed tail identity passes an independent determinant/sign check:
for index `n`, `alpha-p_n/q_n` has sign `(-1)^n`, so `F_k` has sign
`(-1)^(n+1)`. At index four this gives the required negative cubic norm.
The proof correctly distinguishes the defining polynomial (all orders)
from the independently certified number-field norm interpretation
(`k=3,4,5`). Positivity is essential; the stated quartic negative-input
hostile prevents an invalid all-sign extension.

The exact Legendre declaration and interpretation were independently
checked in the [official mathlib documentation](https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/DiophantineApproximation/Basic.html#Real.exists_rat_eq_convergent).
This review makes no local Lean-build claim.

## GSMA denominator, kernel, and Thue–Mahler scope

For primitive `(p,q)` with `p!=0`, the modular gcd calculation proves
that `p*q^(k-1)/F_k(p,q)` has denominator exactly `|F_k(p,q)|`.
The case `q=0` is correctly retained and gives zero with denominator
one; irrationality excludes a zero denominator. Consequently both
integrality iff unit norm and fixed-`S` integrality iff `S`-supported
norm are genuine equivalences, not just necessary conditions.

The localized-unit statement is valid in the stated order, without
identifying a maximal order: the adjugate supplies an inverse when
the determinant is invertible in `Z[S^{-1}]`. The inverse-power update,
the universal exponent-`k` return, and the local gcd sidecar are correct.

The rank-two object is correctly an **additive kernel**, not a subring
or the entire unit group. The finite cyclic source census is not used
as a full Thue census. The cubic sign reversal and the exceptional
generating-function input outside its Taylor disk are both correct.
Quartic/quintic global classifications remain open in this packet.

The fixed-finite-`S`, nonzero-right-side, irreducible-degree-at-least-three
hypotheses match equation (2) and the following finiteness statement in
[Gherga–Siksek's primary paper](https://arxiv.org/html/2207.14492v2).
No effective bound or complete list is inferred from that citation.
The GSMA name is honestly presented as a proposed typed instantiation,
not as a recovered formalism.

## Independent execution

Both checkers were copied to separate temporary working directories
and run in normal and optimized Python modes. All four runs exited
zero, and each normal/optimized output pair was byte-identical.
`PEER_RUNS.json` records source and output hashes. This confirms the
finite certificates; the all-order conclusion rests on the proof above,
not on finite testing. No reviewed files were edited.

## Final-byte and synthesis review

The final uniform checker, SHA-256
`7ba382b7516bec990598927603d42cd19bdb419640cce956fa50abc0f264480d`,
was re-replayed after the certificate-only addition of the full 174-term
prefix through the next-digit-503 witness. Its output JSON, SHA-256
`c0b9526512bb6e4fa103eb5713cc8471c844fdf5e9cb848365c6b2ee0aa8eac8`,
is included in the final reviewed-source pins.

The root `SYNTHESIS.md` also passes the requested mathematical scope review.
Its cubic and Lucas closures explicitly depend on the cited cubic theorem;
the quartic/quintic claims remain bounded and do not imply global solution
lists. The binary-area formula and its degree-greater-than-two condition
are correct, and field degree is kept separate from binary rank. One
metric-convention clarification was requested and is now present: the
length bound greater than 180 is for the specified `sqrt(2)`-weighted
Minkowski embedding, not every Euclidean coordinate convention.
