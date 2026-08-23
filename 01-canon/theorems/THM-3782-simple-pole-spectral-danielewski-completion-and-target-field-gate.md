---
id: THM-3782
title: "Simple-pole spectral Danielewski completion and target-field gate"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  A reduced
  planar rational Keller seed whose poles are cleared by U automatically
  has a finite constant component spectrum for its log-canonical blowdown
  V=UP, and admits a canonical polynomial completion
  E=Phi(V)/U.  Its complete relation and Poisson laws are those of a smooth
  exponent-one Danielewski surface.  One spectral value gives an actual
  polynomial mate.  Two or more values give the nonexact logarithmic
  symplectic obstruction.  If the induced etale map meets every target
  divisor, the Danielewski ring is the complete polynomial intersection
  with the rational target field and all target-field Keller pairs are
  excluded.  Codimension-one coverage is an explicit hypothesis, not an
  automatic conclusion.
source: root + jc_quartic_c3_construct / THM-3779 abstraction, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The exact companion checks the product
  Jacobian identity behind automatic component constancy, the universal
  relation, three Poisson laws, squarefree Bezout/unit-minor certificate for
  one through six spectral values, the synchronized polynomial-mate branch,
  a two-value failed-surjectivity/Darboux hostile, the nonreduced and
  horizontal-pole boundaries, and THM-3779's m=2 positive control.  Normal
  and optimized replay and independent audit of the tangent-constant-field,
  conditional DVR-descent, and THM-3600 steps remain due.
depends_on:
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
related:
  - THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate
  - THM-3779-three-component-tower-maximal-danielewski-polynomial-observable
  - THM-3777-three-component-tower-reciprocal-payment-horizontal-divisor-transfer
script: 04-computation/jc2_simple_pole_spectral_danielewski_gate_thm3782.py
output: 05-knowledge/results/jc2_simple_pole_spectral_danielewski_gate_thm3782.out
script_sha256: f166248449b5735e4f21b4ffc2fa5abcfa44673941646ad017d6cb55c5379836
output_sha256: 93900416967c17e500242e2fe3dab1c1144d3c9132bca272f65ce603a9262abc
semantic_sha256: 12e3dccda5c9f9413862e1bd9110e59b97302ce7201f656e40d5831ca40a73a0
hash_basis: raw LF bytes
---

# THM-3782 -- simple poles canonically complete to exponent one

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  This is
the abstract mechanism exposed by THM-3779.  It separates an unconditional
spectral completion from a conditional maximal-intersection theorem.  It
constructs no planar Jacobian counterexample.

Let `k` be an algebraically closed field of characteristic zero and put

```text
R=k[x,y],                         L=Frac R.             (1)
```

Suppose:

1. `U in R` is nonconstant and reduced, with factorization

   ```text
   U=gamma product_(i=1)^r D_i,       gamma in k^*,    (2)
   ```

   where the `D_i` are pairwise nonassociate irreducibles;
2. `P in L` satisfies the rational Keller equation

   ```text
   J(U,P)=1;                                              (3)
   ```

3. `P` has no horizontal poles and at most simple poles on `U=0`, in the
   exact algebraic sense that

   ```text
   V:=UP lies in R;                                      (4)
   ```

The constant component spectrum is then automatic.  Put `H_i=U/D_i`.
From `(3),(4)`,

```text
J(U,V)=U.                                              (4a)
```

Reducing

```text
J(D_i H_i,V)=D_i J(H_i,V)+H_i J(D_i,V)=D_i H_i        (4b)
```

modulo `D_i` gives `J(D_i,V)=0` in the domain `R/(D_i)`, because
reducedness makes the class of `H_i` nonzero.  On the irreducible geometric
curve `D_i=0`, the tangent derivation

```text
delta_i=(D_i)_y partial_x-(D_i)_x partial_y           (4c)
```

is nonzero on the function field `k(D_i)` and kills the class of `V`.
In characteristic zero the constant field of a nonzero derivation on a
one-variable function field has transcendence degree zero over `k`; since
`k` is algebraically closed, it equals `k`.  Hence the coordinate-ring
class itself is a scalar:

```text
V=lambda_i mod D_i,                 lambda_i in k.     (5)
```

Let `Lambda` be the set of distinct values among the `lambda_i`, put
`h=|Lambda|>=1`, and define the squarefree spectral polynomial

```text
Phi(z)=product_(lambda in Lambda)(z-lambda).            (6)
```

## 1. Unconditional spectral completion

For every `i`, equations `(5),(6)` give

```text
D_i divides Phi(V).                                     (7)
```

Reducedness in `(2)` and unique factorization in `R` therefore give

```text
U divides Phi(V).                                       (8)
```

Consequently

```text
E=Phi(V)/U lies in R.                                   (9)
```

This is the canonical polynomial observable missed by looking only at
`U,V`.  It satisfies the complete relation

```text
UE=Phi(V).                                             (10)
```

Equation `(4a)` makes `U,V` algebraically independent, and

```text
k(U,V)=k(U,P),                 trdeg_k k(U,V)=2.       (11)
```

Hence the map

```text
D_Phi:=k[u,v,e]/(ue-Phi(v)) -> R,
(u,v,e) |-> (U,V,E),                                  (12)
```

is injective: both sides of its image have transcendence degree two and the
prime relation `(10)` already has height one.

The remaining Poisson laws are equally rigid.  In `k(U,P)`, use `V=UP`
and `E=Phi(V)/U` to obtain

```text
{U,V}=U,
{U,E}=Phi'(V),
{V,E}=E.                                              (13)
```

Because `Phi` is squarefree, choose `A_0,B_0 in k[v]` with

```text
A_0(v)Phi(v)+B_0(v)Phi'(v)=1.                         (14)
```

On `D_Phi`, this becomes

```text
A_0(v)ue+B_0(v)Phi'(v)=1.                             (15)
```

Thus the three minors in `(13)` generate the unit ideal.  The surface

```text
Y_Phi=Spec D_Phi                                      (16)
```

is smooth for the same reason, and

```text
phi:A2 -> Y_Phi,                (x,y) |-> (U,V,E),     (17)
```

is etale: it is a morphism between smooth surfaces with everywhere
isomorphic differential.

This completes the unconditional part.  No image or maximality assertion
has been used.

## 2. One spectral value is exactly the equalizer branch

If `h=1`, write `Lambda={lambda}`.  Then

```text
Phi(V)=V-lambda,
E=(V-lambda)/U=P-lambda/U,
J(U,E)=1.                                             (18)
```

Thus `(U,E)` is a polynomial Keller pair.  The spectral surface itself is
an affine plane because `(10)` eliminates `V=UE+lambda`.

This is the exact synchronized-residue branch of THM-3770: a single value
does not merely avoid an obstruction; it explicitly supplies the polynomial
mate.  The formula is sharp and requires neither surjectivity nor an
intersection theorem.

## 3. Two or more values force the exponent-one class

Now let `k=C` and `h>=2`.  Put

```text
b=-v,                  c=u,                  Sigma(b)=Phi(-b). (19)
```

Equations `(10),(13)` become THM-3600's standard exponent-one Danielewski
surface and bracket:

```text
ce=Sigma(b),
{b,c}=c,              {c,e}=-Sigma'(b),       {b,e}=-e. (20)
```

The polynomial `Sigma` is squarefree and has `h>=2` roots.  THM-3600,
equations `(23),(24)`, proves that the symplectic form

```text
omega=db wedge dc/c                                   (21)
```

has nonzero algebraic de Rham class.  Equivalently, the arm-plane
Liouville potentials differ by logarithmic forms with nonzero root
residues.  If `F,G in D_Phi` had `{F,G}=kappa!=0`, then

```text
omega=kappa^(-1)d(F dG),                              (22)
```

a contradiction.  Therefore

```text
h>=2  =>  D_Phi contains no polynomial Darboux pair.  (23)
```

This conclusion concerns the spectral surface.  It does not yet identify
every source polynomial in the target field.

## 4. Conditional maximal-intersection theorem

Add the exact geometric hypothesis

```text
the etale image phi(A2) contains the generic point of
every prime divisor of Y_Phi.                          (24)
```

Since an etale morphism is open, `(24)` is equivalent to saying that its
complement has codimension at least two.  Put

```text
K=k(U,P)=k(U,V)=Frac D_Phi.                            (25)
```

Then

```text
R intersection K=D_Phi                 inside L.      (26)
```

The inclusion from right to left is `(9),(12)`.  Conversely, take
`H in R intersection K` and a height-one prime `q` of the normal ring
`D_Phi`.  Hypothesis `(24)` and etaleness supply a height-one prime `r` of
`R` over `q`.  The induced extension of DVRs is unramified, hence

```text
ord_r(H)=ord_q(H).                                    (27)
```

The left side is nonnegative because `H in R`; therefore `ord_q(H)>=0` for
every height-one `q`.  Normality and the Krull intersection of height-one
DVRs prove `(26)`.

Combining `(23),(26)` gives the all-target-field obstruction:

```text
h>=2, (24), F,G in k(U,P) intersection R
       => J(F,G) is not a nonzero constant.            (28)
```

Thus no rational symplectic target change, of any word length, can make both
outputs polynomial under the codimension-one coverage hypothesis.

## 5. Every hypothesis has a hostile boundary

### 5.1 Failed codimension-one coverage

Take

```text
U=x(1+xy),                  P=1/x.                    (29)
```

Then `J(U,P)=1`, and

```text
V=UP=1+xy,
Lambda={0,1},
E=V(V-1)/U=y.                                      (30)
```

All unconditional conclusions hold on

```text
Y: ue=v(v-1).                                        (31)
```

But the image misses the target divisor `v=e=0`: when `v=0`, the source
has `1+xy=0` and hence `u=0`, while the points `(u,0,0)` with `u!=0`
remain on `Y`.  The target-field function

```text
x=U/V                                                  (32)
```

extends polynomially on the source but does not lie in `D_Phi`.  The
enlarged intersection contains the Darboux pair `(x,E)=(x,y)`.  Therefore
`(23)` cannot be transferred to the source without `(24)`.

THM-3779's `m=1` inverse-discriminant observable is the finite-cover version
of this same failure.  Its `m>=2` members satisfy `(24)` and are the first
positive controls for `(26)`.

### 5.2 Nonreduced U

Take

```text
U=x^2,                         P=y/(2x).               (33)
```

Then `J(U,P)=1` and `V=UP=xy/2` is polynomial with value zero on the sole
reduced component.  But

```text
Phi(V)=V is divisible by x, not by x^2,
Phi(V)/U=y/(2x) is not polynomial.                    (34)
```

Thus the spectral product pays only `rad(U)`; reducedness in `(2)` is
load-bearing.

### 5.3 A horizontal pole

Take

```text
U=x,                         P=y+1/(x-1).              (35)
```

Again `J(U,P)=1`, but `UP=xy+x/(x-1)` is not polynomial.  Multiplication by
`U` cannot clear the pole on the distinct target fibre `U=1`.  Such a seed
must first pay all horizontal fibres; it is outside hypothesis `(4)`.

## 6. Design consequence and scope

The theorem gives a sharp dichotomy for reduced simple-pole seeds:

```text
one spectral value       -> explicit polynomial mate,
at least two values      -> exponent-one logarithmic obstruction,
missing target divisor   -> extra observables may evade that obstruction.
                                                               (36)
```

Consequently a counterexample design based on finitely many simple vertical
residues cannot stop at the log-canonical blowdown.  To escape `(36)` it must
produce genuine codimension-one noncoverage and control the resulting extra
observables, leave the field `k(U,P)`, retain horizontal poles for a later
operation, or change the pole/multiplicity grammar so the completion is not
exponent one.

The companion checks the universal algebra and unit-minor Bezout identities
for `h=1,...,6`, plus every hostile and positive control in Section 5.
Normal and optimized executions must byte-match the frozen transcript before
promotion.  The theorem proves `(5)` but does not assert that the
codimension-one coverage condition `(24)` is automatic, and it proves no
positive-characteristic analogue or planar Jacobian counterexample.  **QED,
conditional only on independent hostile audit.**
