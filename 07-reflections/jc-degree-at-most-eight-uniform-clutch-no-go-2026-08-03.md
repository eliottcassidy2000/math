# Degree-at-most-eight polynomial clutch no-go

**PROMOTED AS
[THM-3276](../01-canon/theorems/THM-3276-degree-at-most-eight-polynomial-clutch-critical-no-go.md),
VERIFIED-EXACT, AND INDEPENDENTLY HOSTILE-AUDITED.** The independent audit
rederived the gradient resultant, degree/boundary ledgers, exact triangular
jet slope, good-reduction coprimality and critical-point consequence through
separate engines.

## Inheritance

- Closest proved mechanism: THM-3225, whose affine `B` family leaves at
  least 50 units of critical-resultant multiplicity.
- Canonical hostile: `B` meeting the owner divisor `g=S*T`, which gives the
  explicit critical point `(alpha,-1/A'(alpha))`.
- Corrected near miss: THM-3257 shows one tuned degree-eight term can control
  a strict transform, but still leaves 52 Morse critical points; a resultant
  root is not an inverse-cover sheet or a branchwise Keller cofactor.
- Least-used sidecar: the complete eight-jet of `B` at the simple `S` root.

## Candidate result

Retain either THM-3212 accessory response pair

```text
deg V=16,             deg A=8,
2VA'-AV'=2V,          g=S*T,
```

and put

```text
P_B(x,z)=(V(x)z^2+B(x)z)^2+A(x)z+x.
```

For every polynomial `B` of degree at most eight over an algebraic closure
of the accessory field, `P_B` has a critical point. More precisely:

1. if `gcd(B,g)!=1`, the inherited finite-fibre formula gives an explicit
   critical point;
2. if `gcd(B,g)=1`, at least 43 units of the saturated critical-resultant
   multiplicity remain away from `g`.

Thus no such `P_B` is the first coordinate of a polynomial Keller pair.
This is a family exclusion, not `JC(2)`.

## Mechanism

The degree ledger closes uniformly. For `deg B<=8`, the covariant

```text
J_B=2VB'-BV'
```

has degree at most 22. In the twelve-term universal resultant `K_B`, the
unique degree-96 term is `64V^6`. The inherited monic boundary has degree
44, so its saturated quotient `H_B` has degree 52.

At a `T` root where `ord(V)=m in {3,4,5,6}`, the response identity gives
`A=2t/(2-m)+...`. For `b=B(alpha)!=0`, the first surviving coefficient is

```text
[t^(3m-1)]K_B = 16m(m-1)/(m-2) b^5 v^3 !=0.
```

Hence `H_B` has no `T` root in the coprime lane.

At the simple `S` root write `c=B(s)!=0` and expand `B` in its nine Taylor
coefficients. The `S^3` boundary is automatic. Successively killing the
coefficients of `K_B` in degrees 3 through 10 uniquely solves for the eight
remaining jets of `B`. The exact slope at jet `j` is
`16j*c^4*v_1^3`, hence nonzero. The first
two coefficients no longer controlled by a degree-eight jet are

```text
[t^11]K_B=F_11(c),       [t^12]K_B=F_12(c).
```

In the two exact cubic accessory fields their Laurent profiles are

```text
F_11: numerator degree 13, denominator c^8,
F_12: numerator degree 15, denominator c^10,
gcd(numerator(F_11),numerator(F_12))=1.
```

Therefore `ord_S(K_B)<=12` and `ord_S(H_B)<=9`. Since `deg H_B=52`, at
least `52-9=43` units remain away from the whole owner boundary.

The constant leading `y` coefficient of the first gradient equation turns
each such resultant root into an affine common zero of the localized
gradient equations. Distinctness is not claimed; existence is enough.

## Exact evidence

The companion

```text
04-computation/jc_degree_at_most_eight_uniform_clutch_no_go_wildcard_20260803.py
```

rebuilds both cubic fields and response pairs, verifies the response and
boundary identities, uses an exact Laurent-polynomial engine for the eight
successive jet eliminations, and computes the two gcds in the exact cubic
fields. The two normalized untunable-pair digests are

```text
(4,1,1,1): 6241cc936e35e6a2999f0ea2e1d251cd4b5f313d0a53331723aa1397e8a47a05
(3,2,1,1): bace7cc5c6c0123742a242cbbddc147595c87f307507b5f67d9dcec6e97951e4
```

Normal and optimized executions agree byte-for-byte with
`05-knowledge/results/jc_degree_at_most_eight_uniform_clutch_no_go_wildcard_20260803.out`.

## Failure anatomy and next boundary

The first failed implication is

```text
eight freely tuned B-jets -> absorption of all 52 residual units    FALSE.
```

The strongest survivor is a finite-dimensional contact theorem: eight
available jets buy at most nine saturated units at the simple owner root.
This closes every `B`-only clutch through degree eight in the current chart.
Degree nine is the first degree that changes the infinity face (THM-3237),
but the known strict-transform walls still leave 50 or more critical points.

The next lawful probes should change a different coordinate of the gradient
system: nonconstant `C_0`, nonlinear `E_0`, or a construction which supplies
an actual inverse-cover/cofactor sidecar. Another `B` retuning below degree
nine is now exhausted.
