---
id: THM-3978
title: "Linear-seam submersion and exact response ideal"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Over an algebraically
  closed field of characteristic zero, for n>=2 and c nonzero, the height-n
  completion coordinate A=x+c(z-1) is globally submersive. Every rational
  solution of J(A,Q)=R(A) is classified explicitly. The plane response ideal
  is (A^(n-1)), whereas the completion response ideal is
  ([A(A+c)]^(n-1)); the added boundary contributes exactly the second factor.
  Thus A has an exact rational constant-Jacobian mate but no polynomial mate,
  and the first nonzero completion response has an explicit primitive. No
  Darboux pair or counterexample to JC(2) is constructed.
source: jc-cohn3709 + all-frontiers / post-THM-3973 critical-free linear-seam lane, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS AFTER REPAIR (audit-3978-candidate,
  2026-08-24). The first pass
  passed the theorem's submersion, rational classification, valuation, and
  response-ideal mathematics, but caught an ill-typed quotient: J(A,-) does
  not preserve B_n. The theorem now uses the typed image intersection
  J(A,B_n) intersection k[A]. It also repairs the rational-RHS reduction via
  k[x,t] intersection k(A)=k[A]. The former assertion-only scout has been
  replaced by an optimization-safe 53-gate companion that freezes the
  non-endomorphism pole and the c=0,n=1 hostiles. The re-audit checked the
  typed response ideal, rational classification, both valuation iff
  directions, general solution family, and all hostiles. Normal, optimized,
  and frozen outputs agree after LF normalization; all hashes pass.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
related:
  - THM-3348-linear-in-time-coordinate-response-ideal
  - THM-3975-danielewski-one-arm-modification-cubic-control-and-hyperelliptic-no-mate
script: 04-computation/jc2_linear_seam_response_ideal_thm3978.py
output: 05-knowledge/results/jc2_linear_seam_response_ideal_thm3978.out
script_sha256: 0600aa8bd898fc3f579cac377e3f3ad8b7aada4f315b95f539d652a02eefce73
output_sha256: 1c1623136f1c6b68e07225d76b8d609d8bd98b5e95df8f0e737f9d591806b253
semantic_sha256: 42ca0b98af0e696a99d7386caabc23cfec5eb26aad434075cd8d24535b693bf3
hash_basis: raw LF bytes
---

# THM-3978 -- the completion doubles the linear-seam response factor

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. For `n>=2`, inside
`k[x,t]` put

```text
z=1+x^n t,              p=zt,              y=x^(n-1)zt^2,
B_n=k[x,z,p,y],          X_n=Spec(B_n).                    (1)
```

Fix `c in k^*` and define the linear seam

```text
A=x+c(z-1)=x+c x^n t.                                   (2)
```

Let `D_A=J_(x,t)(A,-)` act on the rational function field. Then:

1. `A` is a submersion on the source plane and on `X_n`;
2. for arbitrary `R,H in k(A)`, every rational solution is

   ```text
   D_A(Q)=R(A)
   iff
   Q=R(A)/(c(n-1)x^(n-1))+H(A);                         (3)
   ```

3. the two polynomial response ideals are

   ```text
   D_A(k[x,t]) intersection k[A]=(A^(n-1)),
   D_A(B_n) intersection k[A]=([A(A+c)]^(n-1)).          (4)
   ```

The second line is an intersection inside `k(x,t)`: `D_A` does **not**
preserve `B_n`. This type distinction is load-bearing.

In particular

```text
Q_0=1/(c(n-1)x^(n-1)),             J(A,Q_0)=1,           (5)
```

but there is no constant-Jacobian mate in `k[x,t]`, hence none in `B_n`.
The completion's first nonzero response has an explicit primitive given in
Section 5.

## 1. Global submersivity

The source partials are

```text
A_x=1+cnx^(n-1)t,                 A_t=cx^n.              (6)
```

They satisfy the exact Bezout identity

```text
(1-cnx^(n-1)t)A_x+cn^2x^(n-2)t^2 A_t=1.                 (7)
```

Thus `A` has no source-plane critical point. THM-3973 identifies the only
added boundary as `D=V(x,z,p)`. In its smooth chart

```text
z(z-1)^2=x^(n+1)y,
```

the cotangent relation along `D` is `dz=0`; hence

```text
dA|D=dx.                                                 (8)
```

The coordinate `x` is a local parameter there, so `(8)` is nonzero. This
proves submersivity on all of `X_n`.

## 2. Complete rational solution classification

Since `c!=0`, equation `(2)` gives

```text
t=(A-x)/(cx^n),                    k(x,t)=k(A,x).         (9)
```

At fixed `A`, the Hamiltonian derivation is

```text
D_A=-cx^n partial_x.                                    (10)
```

Integrating `D_A(Q)=R(A)` in the rational function field gives `(3)`.
This is a classification, not merely one displayed solution: the kernel of
`partial_x` on `k(A,x)` is exactly `k(A)`.

Suppose first that `Q in k[x,t]` and `D_A(Q)=R(A)`. The bracket is a
polynomial, so

```text
k[x,t] intersection k(A)=k[A].                          (11)
```

forces `R in k[A]`. Equality `(11)` follows by applying a Bezout identity for
coprime numerator and denominator polynomials in `k[A]`. Now the first
summand of `(3)` is regular at every point with `x!=0`. Hence `H(A)` cannot
have a finite pole: every fibre `A=a` contains such a point, where that pole
would remain visible. Thus `H in k[A]`.

## 3. The source-plane response ideal

On the source plane,

```text
A=x(1+cx^(n-1)t),                    ord_(x=0)(A)=1.      (12)
```

For `F in k[A]`, the first summand of `(3)` is polynomial exactly when

```text
ord_(A=0)(F)>=n-1,
```

which is the first equality in `(4)`. The minimal plane primitive is

```text
H_0=A/x=1+cx^(n-1)t,
J(A,H_0^(n-1)/(c(n-1)))=A^(n-1).                        (13)
```

Taking `F=1` proves the polynomial constant-mate nonentry already in the
ambient plane.

## 4. The added boundary contributes a second target factor

On `X_n`, the zero divisor of `x` has two reduced components

```text
div_Xn(x)=D+L_1,                  L_1=V(x,z-1,y).         (14)
```

The target values and local orders are

```text
A|D=-c,                    A|L_1=0,
ord_D(A+c)=1,              ord_L1(A)=1,                  (15)
(A+c)|L_1=c,               A|D=-c.
```

The first summand of `(3)` has no possible pole away from `D union L_1`.
Since `X_n` is normal, it is regular exactly when its orders on those two
primes are nonnegative. Equations `(14)--(15)` give

```text
ord_0(F)>=n-1,                    ord_(-c)(F)>=n-1.       (16)
```

Because `0!=-c`, condition `(16)` is equivalent to

```text
[A(A+c)]^(n-1) divides F(A).                             (17)
```

This proves necessity in the second equality of `(4)`. More precisely, write
`F=[A(A+c)]^(n-1)G(A)`. Whenever `(17)` holds, every `B_n` solution is
`G(A)Q_min+H(A)`, where `Q_min` is the primitive below and `H in k[A]`.

For a rational right-hand side, no extra row is hidden. If `Q in B_n` and
`D_A(Q)=R(A)`, then `Q in k[x,t]` makes the bracket polynomial, so
`R(A) in k[x,t] intersection k(A)=k[A]`; the preceding argument applies.

## 5. Explicit minimal completion primitive

The two required factors cancel the two components of `div(x)` together:

```text
S=A(A+c)/x
 =x+c(2z-1)+c^2x^(n-1)p in B_n.                         (18)
```

Consequently

```text
Q_min=S^(n-1)/(c(n-1)) in B_n,
J(A,Q_min)=[A(A+c)]^(n-1).                              (19)
```

This proves sufficiency in `(4)` constructively. Geometrically, the fibres
`A=0` and `A=-c` each split into a distinguished component (`L_1` or `D`)
and an interior companion seam. A pole in `H(A)` introduced to cancel one
component would also appear on its companion, explaining why the two target
factors cannot be bypassed.

## 6. The derivation is not an endomorphism of the completion

Although the response ideal in `(4)` is intrinsic, the derivation itself does
not map every completion generator back to `B_n`. Directly,

```text
J(A,y)=(z-1)^2/x+2x^(n-1)p+c(n+1)x^(n-1)y.              (20)
```

Along `D`, the first numerator is a unit and `x` is a uniformizer, so `(20)`
has order `-1`. Hence expressions such as `B_n/D_A(B_n)` are ill-typed. The
correct object is the `k[A]`-ideal `D_A(B_n) intersection k[A]` inside the
function field, as used throughout.

## 7. Hostile endpoints and scope

The assumptions are sharp for this mechanism:

- if `c=0`, then `A=x` and `t` is a source-plane mate;
- if `n=1`, then `A=x+cxt` has the critical point `(0,-1/c)`;
- in positive characteristic, `n-1` can vanish and derivation kernels can
  enlarge.

THM-3348 already contains the source-plane factor `A^(n-1)` in its broader
linear-in-time grammar. The new content is the completion factor
`(A+c)^(n-1)`, its exact boundary mechanism, and the full comparison in
`(4)`. Unlike THM-3975's coordinate `p`, the generic `A`-fibre admits an exact
rational time form; here the obstruction arises only when extending it across
the two split target fibres.

The theorem closes this linear-seam coordinate, not arbitrary coordinates on
`B_n`. It produces a nonconstant response pair `(19)`, not a Darboux pair, a
finite map, or a counterexample to `JC(2)`.

## 8. Exact companion

The companion has 53 optimization-safe gates. At heights `2<=n<=9` it checks
the source Bezout identity, rational mate, both minimal responses, the split
factorization, and the non-endomorphism pole. It also freezes both boundary
values and the `c=0,n=1` hostile endpoints. The bounded height loop is a
transcription control; equations `(7)--(20)` prove the all-height statement.

Reproduce with

```bash
python3 04-computation/jc2_linear_seam_response_ideal_thm3978.py
python3 -O 04-computation/jc2_linear_seam_response_ideal_thm3978.py
```

Both runs print `CHECKS=53` and, after LF normalization, byte-match the frozen
output. **QED.**
