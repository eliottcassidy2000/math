# The cyclic third coordinate has a two-root certificate in every grade

**Status: independent exact audit of the all-grade `z` gate in THM-3517.**
This route imports neither Sympy nor the existing coefficient-recurrence
companion.  It concerns only THM-3448's explicit cyclic weighted family.

## Inheritance and the disjoint route

THM-3517 proves that the numerator

```text
H=gamma[gamma(gamma-1+a)-aw]
```

has a nonconstant remainder modulo

```text
T_n(w)=w^n-w^(n-1)+Pw-Q.
```

Its primary proof computes the top remainder coefficient.  A disjoint test
is available at the rational target

```text
P=2^(-(n-1)),       Q=0,       C=1.                  (1)
```

The inverse polynomial then has the two visible roots `0` and `1/2` for
every `n>=3`.  Put `m=n-3` and `delta=1/(2n-4)`.  Direct substitution into
the cyclic seed gives

```text
gamma(0)=P,                  gamma(1/2)=-mP,
H(0)=P^2(P-2-delta),
H(1/2)=-mP[mP(mP+2+delta)+(1+delta)/2].               (2)
```

For `n=3`, the two values satisfy `H(0)<0=H(1/2)`.  For `n>=4`,

```text
|H(1/2)| > P/2,
|H(0)| <= (9/32)P,
```

because `m>=1`, `P<=1/8`, and `delta<=1/4`.  Thus the values are different
in every grade.  If `H mod T_n` were constant over `Q(P,Q)`, every
specialization of its nonconstant coefficients would vanish; (2) disproves
that.  Therefore `z` is not in the target field.  THM-3438's maximal
`S_(n-1)` point stabilizer then leaves no proper intermediate field, so `z`
is primitive.

## Boundary and forward replay

At `n=3`, the half-root has `gamma=0`.  It is therefore used only as a
polynomial-identity witness, not as a finite reconstructed source point.
This is the sharp boundary of the stronger geometric replay.

For every `4<=n<=256`, exact rational reconstruction gives two finite source
points above the same target `(A,B,C)=(0,2^(1-n),1)`, and direct evaluation
of the weighted map returns that target from both points while their actual
`z` coordinates differ.  The all-degree claim rests on the uniform
inequality, not on the finite range.

The certificate loses the other `n-2` sheets, all branch-divisor
multiplicities, and every Jelonek boundary statement.  It proves no claim
for arbitrary weighted seeds, the unstored historical `E_m`, arbitrary
Keller maps, `JC(2)`, LRC, ancestry, current, or `H^1`.

## Exact artifacts

```text
04-computation/jc_weighted_cyclic_z_primitivity_two_root_independent_audit_20260816.py
05-knowledge/results/jc_weighted_cyclic_z_primitivity_two_root_independent_audit_20260816.out
```

The script uses only `fractions.Fraction`, has no randomness or elapsed
fields, and has byte-identical normal/optimized transcripts.  Its semantic
digest is

```text
2c8d9f191ae5a8dddd6ef15feeb20d6564491bd0be6ea34c8e662428f74b729b.
```
