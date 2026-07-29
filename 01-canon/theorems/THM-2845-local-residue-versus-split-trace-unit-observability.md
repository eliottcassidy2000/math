---
id: THM-2845
title: "Local residue versus split trace unit observability"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let A be a
  finite-dimensional unital
  associative K-algebra and let ell:A->K be linear.  If K is not F_2,
  then ell is nonzero on every unit exactly when ell is a nonzero scalar
  multiple of a K-algebra character A->K.  Over F_2 the complete extra
  family consists of the odd sums of distinct characters.  Equivalently,
  a detector exists exactly when A/J(A) has a K factor.  Every detector
  kills J(A).  A two-way scalar unit test, unit iff ell is nonzero,
  exists exactly when A is local with residue K, in which case ell is a
  scalar residue map.  Thus modular p-group augmentation is the unique
  local detector, while a split weighted trace detects all units only
  on one coordinate, except for the sharp odd-parity F_2 family.
source: root/local-versus-split-unit-observability-2026-07-28
depends_on: []
related:
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - THM-2839-prime-power-unit-mass-full-spectrum-and-q11-response-provenance
  - THM-2842-ordered-positive-cone-vandermonde-multiplier-observability
script: 04-computation/algebra_local_split_unit_observability_thm2845.py
output: 05-knowledge/results/algebra_local_split_unit_observability_thm2845.out
script_sha256: 296abc5461741f96d491defaa744fe2fe78fbcb14f4018cbeba675889dc78d44
output_sha256: 807a59ef8329419fb61d012b0bdf6400f7276e4d3e2808622c19a30b8de4799e
independent_script: 04-computation/local_unit_detector_classification_thm2845.py
independent_output: 05-knowledge/results/local_unit_detector_classification_thm2845.out
independent_script_sha256: 7575576d3c48c8bb3d92e6eed7fa6ce9f8df84a88403c6ec01dd8fad551c452a
independent_output_sha256: 6179e5606c3e3c293d4f1cdc04f052b556d4a9aaf00ccc266f6a9929f7aef216
hash_basis: LF-normalized bytes
---

# THM-2845 -- local residue versus split trace unit observability

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `A` be a finite-dimensional unital associative algebra over a field
`K`.  No commutativity is assumed.  Write

```text
X_K(A)=Hom_(unital K-algebras)(A,K)                     (1)
```

for its set of scalar characters, and call a linear functional
`ell:A->K` a **unit detector** when

```text
ell(u)!=0                         for every u in A^x.   (2)
```

Condition `(2)` is one-way: it asks the kernel hyperplane to contain no
unit.  It does not say that every element outside the hyperplane is a
unit.

## 1. Complete classification

If `K` is not the two-element field, the unit detectors are exactly

```text
{c chi : c in K^x, chi in X_K(A)}.                     (3)
```

If `K=F_2`, they are exactly

```text
{sum_(chi in S) chi :
       S subset X_K(A), |S| is odd}.                   (4)
```

Distinct odd subsets give distinct sums in `(4)`, because the characters
are the independent coordinate projections onto the scalar factors of
`A/J(A)`.  In particular, `(4)` is nonempty only when `X_K(A)` is
nonempty.

Equivalently, over every field:

```text
A has a linear unit detector
 iff A/J(A) has a simple direct factor isomorphic to K
 iff A has a K-algebra character.                       (5)
```

The first important correction to the reserved formulation is therefore:
locality is not necessary for the one-way condition `(2)`.  For example,
the first projection `K x K -> K` is a unit detector.

## 2. The radical is always invisible

Every unit detector satisfies

```text
ell(J(A))=0.                                           (6)
```

Indeed, `ell(1)!=0`.  If `j in J(A)` and `ell(j)!=0`, choose

```text
t=-ell(1)/ell(j).
```

But `1+tj` is a unit, whereas `ell(1+tj)=0`, contradicting `(2)`.
Thus `ell` descends to

```text
A/J(A),                                                (7)
```

and an element of `A` is a unit exactly when its image in `(7)` is a
unit.  This proves that nilpotent and noncommutative radical structure
creates no additional detector.

## 3. Infinite-field rigidity

Assume first that `K` is infinite.  Put `n=dim_K A` and define the left
regular norm polynomial

```text
N(a)=det_K(L_a),               L_a(x)=ax.               (8)
```

An element `a` is a unit exactly when `N(a)!=0`.  By `(2)`, every
`K`-point of the hyperplane `ker ell` lies in `N=0`.  A polynomial over
an infinite field which vanishes on a whole affine space vanishes
identically there.  Hence

```text
ell divides N in K[A].                                  (9)
```

Let `u` be a unit.  Since

```text
N(ua)=N(u)N(a),                                        (10)
```

the linear form `ell o L_u` is another linear factor of `N`.  We claim
that it is proportional to `ell`.

Consider the affine unit path

```text
u_t=(1-t)1+tu.                                         (11)
```

The polynomial `N(u_t)` is nonzero at `t=0`, so all but finitely many
`t in K` give units.  For every such `t`,

```text
ell o L_(u_t)=(1-t)ell+t(ell o L_u)                    (12)
```

is a linear factor of `N`.  The polynomial `N` has only finitely many
projective linear factors.  If `ell` and `ell o L_u` were independent,
the right side of `(12)` would give infinitely many distinct projective
linear factors.  Therefore

```text
ell(ua)=c(u)ell(a)                         for all a.   (13)
```

Normalize

```text
chi=ell/ell(1).                                        (14)
```

Putting `a=1` in `(13)` gives `c(u)=chi(u)`, and hence

```text
chi(ua)=chi(u)chi(a)                      for units u. (15)
```

For arbitrary `x in A`, the monic polynomial

```text
det(L_x+tI)
```

has only finitely many roots.  Choose `t in K` such that `x+t1` is a
unit, apply `(15)`, expand by linearity, and cancel the two `t chi(a)`
terms.  This gives

```text
chi(xa)=chi(x)chi(a)                      for all x,a. (16)
```

Thus `chi` is a unital `K`-algebra character, proving `(3)` for infinite
fields.  Conversely, every scalar multiple of a character detects
units because

```text
chi(u)chi(u^(-1))=1.                                  (17)
```

This proof also explains the mechanism: over an infinite field a
hyperplane contained in the regular-norm hypersurface must be a linear
norm factor, and the affine pencils containing infinitely many units
force that factor to be multiplicative.

## 4. Finite-field classification

Now let `K=F_q`.  By `(6)` it is enough to work in the semisimple
quotient.  Artin--Wedderburn and Wedderburn's little theorem give

```text
A/J(A) =
  product_i M_(n_i)(F_(q^(d_i))).                      (18)
```

Call a factor **scalar** when `(n_i,d_i)=(1,1)`, so that it is `K`.

### Lemma: every nonscalar simple component has full unit-value set

Let

```text
B=M_n(E),                    E=F_(q^d),                 (19)
```

with `(n,d)!=(1,1)`, and let `lambda:B->K` be nonzero and `K`-linear.
Then

```text
lambda(B^x)=K.                                         (20)
```

If `n=1` and `d>1`, every fiber of `lambda:E->K` has
`q^(d-1)` elements.  A nonzero fiber contains only nonzero elements, and
the zero fiber contains a nonzero element because `d>1`.  Since every
nonzero field element is a unit, `(20)` follows.

Suppose `n>=2`.  The finite-field trace pairing is nondegenerate, so
there is a nonzero `C in M_n(E)` with

```text
lambda(X)=Tr_(E/K)(tr(CX)).                            (21)
```

Left-right equivalence and cyclicity of matrix trace reduce `C` to

```text
D=diag(I_r,0),                         1<=r<=n.        (22)
```

For any `alpha in E`, use the invertible block

```text
B_beta = [ beta  1 ]
         [   1   0 ],                 det B_beta=-1.  (23)
```

If `r=1`, choose `beta=alpha`; if `r>=2`, choose
`beta=alpha-(r-2)`.  Then

```text
Y=B_beta direct-sum I_(n-2)
```

is invertible and `tr(DY)=alpha`.  The field trace `E->K` is surjective,
so `(20)` follows.

Consequently, if a descended detector were nonzero on any nonscalar
factor, that factor could contribute the negative of the fixed
contribution from chosen units in all the other factors.  Their product
would be a unit in the kernel.  Thus every detector vanishes on every
nonscalar simple factor.

Suppose there are `r` scalar factors.  On them write

```text
ell(x_1,...,x_r)=c_1x_1+...+c_rx_r.                    (24)
```

Let `m` be the number of nonzero `c_i`.  On units, each nonzero summand
`c_ix_i` can be any element of `K^x`.

If `q>2`, any `m>=2` nonzero field elements can be chosen to sum to
zero.  For even `m`, use pairs `1,-1`.  For odd `m>=3`, choose

```text
1, a, -1-a,                 a notin {0,-1},            (25)
```

and then pairs.  Hence `(2)` holds exactly when `m=1`, proving `(3)`.

If `q=2`, every nonzero summand in `(24)` equals one on every unit.
Therefore

```text
ell(u)=m mod 2,                                        (26)
```

and `(2)` holds exactly when `m` is odd.  Scalar-coordinate projections
are precisely the distinct characters, so this proves `(4)`.

Finally, every character kills `J(A)` and selects exactly one scalar
factor in `(18)`.  Conversely, projection to any scalar factor gives a
character.  This proves `(5)` and completes the classification.

For completeness, the equivalence `(5)` over an arbitrary field uses
the same Artin--Wedderburn argument.  A character descends to

```text
A/J(A)=product_i M_(n_i)(D_i).
```

The central factor idempotents map to zeros and a single one, so the
character selects one simple factor.  Its restriction to that factor is
a unital homomorphism to `K`; simplicity makes it injective.  Since it is
also `K`-linear, that factor must have dimension one and hence be `K`.
The converse is again coordinate projection.

## 5. When one scalar is an exact unit test

Strengthen `(2)` to the two-way condition

```text
a in A^x  iff  ell(a)!=0.                              (27)
```

Such an `ell` exists exactly when

```text
A is local and A/J(A)=K.                               (28)
```

In that case every such functional is a nonzero scalar multiple of the
residue map

```text
rho:A->A/J(A)=K,                                       (29)
```

and `(27)` is the standard unit-lifting criterion.

For the converse, apply the complete detector classification.  If
`K!=F_2`, a detector is supported on one scalar factor.  Any additional
semisimple factor can be set to zero while the observed scalar
coordinate remains nonzero, producing a nonunit outside the kernel.  If
`K=F_2`, the same argument applies to an odd character sum: a single
observed scalar coordinate already gives a nonunit with value one
unless it is the only semisimple coordinate.  Thus `(27)` forces
`A/J(A)=K`, which is `(28)`.

This is the precise role of locality.  A split quotient supplies a
one-way detector; a local residue `K` supplies an exact scalar
classification of units.

## 6. Two sharp applications

### Modular p-group augmentation

Let `G` be a finite `p`-group.  THM-2839 proves

```text
F_p[G]/J(F_p[G])=F_p
```

and identifies the residue map with augmentation `epsilon`.  The local
corollary therefore shows:

```text
every linear unit detector on F_p[G] is c epsilon.     (30)
```

This includes nonabelian `G`.  THM-2839 proves the stronger reverse
direction

```text
epsilon(a)!=0  =>  a is a unit,                        (31)
```

because the augmentation ideal is the whole Jacobson radical.
By contrast, the easy direction

```text
a is a unit  =>  epsilon(a)!=0                         (32)
```

uses only that augmentation is a character and holds for every group,
not just a `p`-group.  Thus THM-2839's `p`-group hypothesis is needed for
`(31)`, full-spectrum invertibility, and signed reconstruction, not for
the detector implication `(32)`.

### Split weighted traces and the Laguerre carrier

On a split algebra `K^r`, a weighted trace

```text
T_w(x)=sum_i w_ix_i                                   (33)
```

detects every unit exactly when:

```text
K!=F_2:  precisely one weight is nonzero;
K=F_2:   an odd number of weights is nonzero.          (34)
```

In particular, the full trace on `F_2^r` is a detector exactly for odd
`r`.  This is the sharp `F_2` exception; `UT_3(F_2)` in the exact
companion is a noncommutative radical-bearing realization.

After extending THM-2815's Laguerre quotient to its complex splitting
field, THM-2842's readout has the form

```text
Lambda_D(f)=sum_(i=1)^D w_i f(x_i),          w_i!=0.   (35)
```

For `D>=2`, `(34)` forces a trace-zero unit.  THM-2842 does more: it
constructs the rational canonical witness `ell_(D-1)` before splitting
and proves its sharp `D-1`-step Krylov invisibility.  The present theorem
also sees rational existence directly: factorial exactness gives

```text
Lambda_D(s)=1,                  Lambda_D(s^2)=2,
```

so `Lambda_D` is not multiplicative and therefore cannot be a detector
over the infinite field `Q`.

The selector statement also aligns exactly.  A multiplier `g` changes
the weights to `w_i g(x_i)`.  To detect **all** units, `(34)` requires
`g` to be a nonzero multiple of one cardinal idempotent.  THM-2842 asks
the weaker conditional question only on `ker Lambda_D`; adding
`c Lambda_D` is then invisible, so its complete family

```text
g=c+d lambda_i,                         d!=0,          (36)
```

is precisely the same coordinate selector modulo the original readout.

## 7. Exact hostile controls

The primary exact companion exhausts:

1. five truncated local algebras and every scalar functional on them;
2. all scalar functionals on `M_2(F_q)` for `q=2,3,5`;
3. every `F_2`-linear scalar on `F_4`;
4. `1,269` weighted split-product scalars; and
5. the finite-field zero-sum count

   ```text
   #{(y_1,...,y_D) in (F_q^*)^D:sum_i y_i=0}
    =((q-1)^D+(q-1)(-1)^D)/q
   ```

   in `24` direct universes through `q=7,D=6`.

The independent companion separately exhausts split products
`F_q^r` for `q=2,3,5` and `r<=5`, the dual numbers, `UT_2(F_3)`,
`UT_3(F_2)`, `M_2(F_2)`, and `F_4/F_2`.  The independent audit further
checked `M_2(F_3)`, `F_9/F_3`, and mixed scalar/matrix products directly.

Reproduce with

```text
python3 04-computation/algebra_local_split_unit_observability_thm2845.py
python3 -O 04-computation/algebra_local_split_unit_observability_thm2845.py
python3 04-computation/local_unit_detector_classification_thm2845.py
python3 -O 04-computation/local_unit_detector_classification_thm2845.py
```

Both companions byte-match their stored transcripts in normal and
optimized modes.

## 8. Scope

This theorem classifies scalar hyperplanes disjoint from the unit group.
It does not turn a one-way detector into a positive inverse, a
node-selector supplied by scalar moment nullity, a physical LRC
multiplier, or a Jacobian injectivity certificate.  In particular,
THM-2839's signed inverse and THM-2842's external multiplier-access
boundaries remain load-bearing.

**QED.**
