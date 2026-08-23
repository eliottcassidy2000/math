---
id: THM-3833
title: "ABC-conditional cube radical growth and hyperbolic power finiteness"
status: >
  PROVED CONDITIONAL IMPLICATION + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED.  Assuming the standard ABC conjecture, primitive positive
  two-cube values have an explicit radical floor and only finitely many can
  satisfy rad(m)<=m^delta for fixed delta<1/3.  ABC also implies
  non-effective finiteness for every fixed primitive hyperbolic
  Fermat--Catalan signature.  ABC itself remains OPEN.  The published but
  disputed IUT-to-ABC implication is not a dependency.
source: root + abc_repo_transfer + iut_structure_audit / abc-IUT bridge session, 2026-08-23
audit: >
  PASS.  abc_repo_transfer independently rederived both conditional
  inequalities, the hyperbolic exponent, and the primitive-scale boundary.
  iut_structure_audit separately checked the literature/status firewall and
  the address-versus-evaluator boundary.  The assertion-free companion checks
  19,314 primitive cube ratios, 5,855 decoder-carrier states, the first 15
  shortlex exponent tasks, and scale, radical-fibre, Pythagorean,
  nonprimitive, Euclidean, and
  six-term hostiles with 154,651 active gates.  Normal and optimized raw LF
  streams match the frozen transcript.
related:
  - THM-3832-nonlinear-cubic-root-ratio-triangular-birational-chart
  - THM-3375-berggren-positive-two-cube-pell-ray
  - THM-3656-berggren-two-cube-frobenian-three-halves-sieve-bound
  - THM-3730-positive-distinct-two-cube-support-abscissa
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
  - THM-3756-odd-square-ordinal-berggren-affine-descent
  - THM-3793-inert-prime-sum-all-scale-two-cube-singleton
  - THM-3825-prime-colour-valuation-two-cube-decoder
script: 04-computation/abc_cube_hyperbolic_task_atlas_thm3833.py
output: 05-knowledge/results/abc_cube_hyperbolic_task_atlas_thm3833.out
script_sha256: 8a25f6cf0bf979adbfe81095fdbb0f873eed1856594c803e64d517433f72e660
output_sha256: d6b0aa913ef95ebd199ff535c3460b702ad34d9e2183eeea72dddaed57b5af61
hash_basis: raw LF bytes
---

# THM-3833 -- what ABC would actually buy on the current arithmetic frontier

**PROVED CONDITIONAL IMPLICATION + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  Every conclusion below has the ABC conjecture as an
explicit antecedent.  The theorem proves no case of ABC and imports no claimed
proof of ABC from inter-universal Teichmuller theory.

## 1. Truth gate and normalization

Write

```text
rad(n)=product_(prime ell|n) ell.
```

The sole conjectural hypothesis is the standard positive-integer form of
**ABC**:

> For every `epsilon>0` there is a constant `K_epsilon>0` such that every
> pairwise-coprime positive triple `A+B=C` satisfies
>
> ```text
> C <= K_epsilon rad(ABC)^(1+epsilon).                 (ABC)
> ```

ABC is **OPEN**.  Constants `K_epsilon` are not known here, so consequences
using them are non-effective.

Primitivity is structural, not cosmetic.  Before `(ABC)` is applied, a common
factor in `A,B,C` must be removed.  The normalized projective equation and an
ambient scaled integer state are different types.

A useful immediate address statement is the following.  If `S` is a fixed
finite prime set, `R_S=product_(ell in S)ell`, and a primitive positive
`A+B=C` has `rad(ABC)|R_S`, then conditionally

```text
C <= K_epsilon R_S^(1+epsilon).                       (1)
```

Thus each fixed radical support has a bounded primitive height fibre.  It may
still have infinitely many ambient scaled representatives if normalization is
forgotten.  Classical fixed-`S` unit-equation theory already gives primitive
finiteness; `(1)` records the simple ABC-conditional height shape, not a new
finiteness theorem.

## 2. Primitive two-cube radical floor

Let

```text
0<a<b,             gcd(a,b)=1,             m=a^3+b^3. (2)
```

### Theorem 2.1

Assume `(ABC)`.  For every `epsilon` with `0<epsilon<1/2`,

```text
rad(m) >= K_epsilon^(-1/(1+epsilon))
          m^(1/(1+epsilon))/rad(ab),                  (3)

rad(m) >  K_epsilon^(-1/(1+epsilon))
          b^((1-2epsilon)/(1+epsilon)),                (4)
```

and

```text
rad(m) >= 2^(2/3) K_epsilon^(-1/(1+epsilon))
          m^((1-2epsilon)/(3(1+epsilon))).             (5)
```

Consequently, for every real `delta<1/3`, only finitely many primitive
positive two-cube values satisfy

```text
rad(m)<=m^delta.                                      (6)
```

Equivalently, along every sequence of such primitive values with `m_n->infinity`,

```text
liminf log(rad(m))/log(m) >= 1/3,
limsup log(m)/log(rad(m)) <= 3.                       (7)
```

### Proof

Primitivity gives

```text
gcd(ab,m)=1,
rad(a^3 b^3 m)=rad(ab)rad(m).                         (8)
```

Apply `(ABC)` to `a^3+b^3=m` and rearrange.  This proves `(3)`.  Since
`rad(ab)<=ab<b^2` and `m>b^3`, equation `(4)` follows.

The arithmetic--geometric mean inequality gives

```text
m=a^3+b^3 >= 2(ab)^(3/2),
ab <= (m/2)^(2/3).                                    (9)
```

Substituting `(9)` into `(3)` proves `(5)`.  Its exponent

```text
theta_epsilon=(1-2epsilon)/(3(1+epsilon))             (10)
```

tends upward to `1/3` as `epsilon` tends to zero.  Given `delta<1/3`, choose
`epsilon>0` with `delta<theta_epsilon`.  Equations `(5)` and `(6)` then bound
`m^(theta_epsilon-delta)` by a fixed constant, hence bound `m`.  This proves
finiteness and `(7)`.  QED

## 3. What the radical sees, and what the decoder restores

On the primitive carrier of
[THM-3825](THM-3825-prime-colour-valuation-two-cube-decoder.md), write

```text
s=a+b,             q=a^2-ab+b^2,             m=sq.    (11)
```

The inert support of `s` and split support of `q` are disjoint.  Equation
`(5)` therefore becomes the conditional prime-colour budget

```text
log rad(s)+log rad(q)
  >= theta_epsilon(log s+log q)-O_epsilon(1).          (12)
```

This bounds aggregate valuation excess.  It does **not** decode a pair, count
the support, imply a singleton fibre, or retain coordinate labels.  THM-3825
is strictly richer on its restricted carrier: its coloured valuation decoder
and square test recover the primitive shell, Eisenstein cofactor, pair, and
inert scale, while a separate 3-adic tag carries label placement and
orientation.

The exact scale kernel is already infinite:

```text
M_k=(2^k)^3+(3*2^k)^3=28*8^k,
rad(M_k)=14                         for every k>=0.     (13)
```

In particular, `28=1^3+3^3` and `224=2^3+6^3` share the same value radical;
even `rad(xy(x^3+y^3))=42` is unchanged.  ABC first divides out the common
cube and sees `1+27=28`.  THM-3825 instead retains `k` through the inert-scale
valuation quotient.  Radical address and ambient scale address answer
different questions.

The taxicab identity

```text
1729=1^3+12^3=9^3+10^3                             (14)
```

is the collision hostile: a radical floor does not imply injectivity or a
collision-tax estimate.  The Frobenius-colour incidence in
[THM-3656](THM-3656-berggren-two-cube-frobenian-three-halves-sieve-bound.md)
is another lost sidecar; a bare radical cannot replace that sieve.

## 4. Fixed hyperbolic Fermat--Catalan signatures

Fix integers `p,q,r>=2` and put

```text
sigma=1/p+1/q+1/r,              kappa=1-sigma.         (15)
```

### Theorem 4.1

Assume `(ABC)` and `sigma<1`.  Then the positive pairwise-coprime solutions of

```text
x^p+y^q=z^r                                            (16)
```

are finite.  More precisely, for

```text
epsilon=kappa/(2sigma),                               (17)
```

every such solution satisfies

```text
(z^r)^(kappa/2)=z^(r kappa/2) <= K_epsilon.           (18)
```

### Proof

The three powered terms in `(16)` are pairwise coprime, so `(ABC)` applies.
Moreover

```text
rad(x^p y^q z^r)=rad(xyz)<=xyz<(z^r)^sigma.           (19)
```

The strict inequality uses `x^p<z^r`, `y^q<z^r`, and
`z=(z^r)^(1/r)`.  Combining `(ABC)`, `(17)`, and `(19)` leaves exponent

```text
1-(1+epsilon)sigma=kappa/2,                           (20)
```

which proves `(18)`.  It bounds `z`, and then `(16)` bounds `x,y`; hence only
finitely many solutions exist.  QED

The primitive qualifier is indispensable.  The first shortlex hyperbolic
signature already has the infinite nonprimitive family

```text
(9t^4)^3+(18t^4)^3=(9t^3)^4             (t>=1).       (21)
```

This family is removed by projective normalization, exactly as `(13)` is.

## 5. Natural-number task atlas: address is not curvature

Let

```text
H={(p,q,r):2<=p<=q<=r and 1/p+1/q+1/r<1}.             (22)
```

Order `H` first by `p+q+r`, then lexicographically.  Each weight shell is
finite, so its ordinal rank is a bijection `H<->N_(>0)`.  The beginning is

```text
1:(3,3,4),  2:(2,4,5),  3:(3,3,5),  4:(3,4,4),
5:(2,3,7).                                             (23)
```

This is a scheduler for unordered exponent patterns, not a numerical
substitute for `sigma` or `kappa`.  A concrete signed equation also needs an
output/sign slot and variable-order sidecar.

The unique least positive curvature is

```text
(p,q,r)=(2,3,7),              kappa=1/42.              (24)
```

This does not make it the first shortlex task.  To prove `(24)`, split on
`p`.  If `p=q=2`, no `r` is hyperbolic.  If `p>=3`, every positive gap is at
least the `(3,3,4)` gap `1/12`.  If `p=2,q=3`, hyperbolicity begins at `r=7` and the gap
`1/6-1/r` increases thereafter.  If `p=2,q=4`, the first gap is `1/20` at
`r=5`; if `q>=5`, the gap is at least `1/10`.  Thus `1/42` is uniquely
minimal.

The Euclidean boundary patterns are exactly

```text
(2,3,6),                 (2,4,4),                 (3,3,3). (25)
```

Indeed, if `p>=4` the reciprocal sum is below one.  For `p=3`, equality
forces `q=r=3`.  For `p=2`, equality forces `q=3,r=6` or `q=r=4`; if `q>=5`
the sum is again below one.

The Pythagorean signature `(2,2,2)` has `sigma=3/2>1`, so Theorem 4.1 does not
apply.  ABC would instead give a radical budget on every primitive ordered
Pythagorean triple `(A,B,C)`:

```text
rad(ABC) >= K_epsilon^(-1/(1+epsilon))
            C^(2/(1+epsilon)).                         (PYTH-ABC)
```

Indeed, apply `(ABC)` to `A^2+B^2=C^2`; primitivity makes the three terms
pairwise coprime and powers do not change their radical.  Along `C_n->infinity`
this would imply

```text
liminf log(rad(A_n B_n C_n))/log(C_n) >= 2.            (25a)
```

In the odd-square chart of
[THM-3756](THM-3756-odd-square-ordinal-berggren-affine-descent.md), put

```text
q=2r-1,       d=2s-1,
(A,B,C)=(qd,(q^2-d^2)/2,(q^2+d^2)/2).                 (25b)
```

The outer rank `r` alone does not determine either side of `(PYTH-ABC)`.
At `r=3`, the two fibre points

```text
(5,12,13),      (15,8,17)
```

share `B+C=25` but have `rad(ABC)=390,510`.  The complementary coordinate
`s` is the missing sidecar.  The budget gives no finiteness of the Berggren
tree and no obstruction to its odd-square task rank.

Similarly, the primitive collision

```text
3^6+19^6+22^6=10^6+15^6+23^6=160426514              (26)
```

is a six-term additive-energy object, not a three-term ABC input.  Reciprocal
exponent syntax alone does not create a transfer.

## 6. Exact finite atlas

The companion makes no use of `(ABC)`.  It checks unconditional interfaces
and hostile controls on the exact THM-3743 support-two universe

```text
U_356={(a,b):1<=a<b, gcd(a,b)=1, a+b<=356}.            (27)
```

There are `19,314` pairs, of which `5,855` lie in THM-3825's primitive
decoder carrier.  For each pair it checks `(8)`, factors `abm`, and records
the diagnostic ABC quality

```text
Q(a,b)=log(m)/log(rad(abm)).                           (28)
```

The largest-quality row in this finite universe is

```text
(a,b)=(1,80),      m=512001,      rad(abm)=9030,
Q=1.443306744511... .                                  (29)
```

This is a finite diagnostic, not evidence for or against ABC.  The executable
also checks the first 15 shortlex tasks, the unique least-curvature control in
`2<=p<=q<=r<=40`, all Euclidean controls, the infinite-family identities
`(13)` and `(21)` on finite ranges, the equal-radical fibre

```text
1+2=3,                 1+8=9,                 rad(ABC)=6, (30)
```

the Pythagorean shell/radical split, the six-term non-transfer, and all
factorization gates.

Run

```bash
python3 -B 04-computation/abc_cube_hyperbolic_task_atlas_thm3833.py
python3 -B -O 04-computation/abc_cube_hyperbolic_task_atlas_thm3833.py
```

Both raw LF streams equal
`05-knowledge/results/abc_cube_hyperbolic_task_atlas_thm3833.out`.

## 7. Explicit non-consequences

- **No IUT import.**  IUT I--IV are published works, and IUT IV claims an
  implication to Vojta/ABC/Szpiro.  That implication remains publicly
  disputed.  This theorem starts from `(ABC)` and does not adjudicate or use
  IUT III Corollary 3.12.
- **No LRC(14) consequence.**  THM-3743 already makes `(27)` finite, while
  `K_epsilon` is unknown.  A cube value loses the other eleven speeds, owner,
  phase, arrival, and loneliness; higher-support relations are multi-term.
- **No Jacobian consequence.**  ABC logarithmic arithmetic height is not the
  total polynomial degree/height in THM-3550.  THM-3832's triangular `(z,C)`
  chart is birational and retains a structural denominator/regularity
  passport; it is not an integer ABC packet.  No Keller-invariant map to an
  unbounded family of primitive ABC triples has been proved.
- **No support or collision asymptotic.**  The radical floor neither improves
  THM-3730's proved support abscissa nor supplies a two-cube support asymptotic.
- **No Pell-ray exclusion.**  THM-3375's infinite positive
  `x^3+y^3=A^2+2` rays are not equations of the form `(16)`; the constant
  summand is not a third powered variable.

The durable output is a typed budget: ABC would constrain primitive height
relative to prime support.  It does not replace exact valuation colours,
scale digits, operation words, or ambient semantic sidecars.
