# Reserved three-direction candidate: analytic and H99 audit

**Proof-candidate status: analytic PASS + independent FINITE-EXACT H99 PASS.**
The reserved THM-4431 namespace remains unproved until explicit promotion.
This note proves the missing colored-lattice basis and projective
classification directly; it does not use the finite generic-body sampler as
an all-height argument.

Precisely, let `w=(a,b,c)` be sorted, primitive, distinct, positive, odd and
ternary-unit.  If its **complete** live raw-carrier support has exactly three
primitive unoriented directions, then every THM-4414 raw network projection
is strictly below `6/77`.  The proof inherits only the raw projection formula
and THM-4386's short-relation rigidity.  It supplies no chart entry,
synchronization, or LRC(14).

## 1. Native colored lattice

For an eligible primitive `w=(a,b,c)`, let

```text
L={C in Z^3 : C.w=0},
Gamma={C in L : aC_1=bC_2=cC_3 (mod 3)},
H=3L.
```

The common residue defines a homomorphism `tau:Gamma->F_3` with kernel
`H`.  Modulo `3L`, the lattice `L` is the two-dimensional kernel of `w`,
while `Gamma/3L` is its one-dimensional equal-coordinate line.  Hence

```text
[L:Gamma]=3,  Gamma/H is cyclic of order 3.
```

The live set is exactly `(Gamma\H) intersect K`, where `K` is the open,
bounded, centrally symmetric convex roof body in the real kernel plane.
Projection to the first two coordinates identifies `L` with

```text
{(x,y) in Z^2 : ax+by=0 (mod c)},
```

which has determinant `c` because `gcd(a,b,c)=1`.  Therefore `Gamma` has
first-two-coordinate determinant `3c`.

Primitive reduction preserves liveness.  If `C=g d` is live, then `3` does
not divide `g`; division by the unit `g mod 3` keeps `d` in `Gamma\H`, and
central convexity keeps `d` in `K`.  Thus every live direction has its
primitive representative in the complete live body.

## 2. Colored-basis descent

**Lemma.** If the complete live body has two nonparallel directions, two live
primitive representatives form a basis of `Gamma`.

Choose a nonparallel live primitive pair `u,v` minimizing its positive
`Gamma`-determinant `D`.  If `D>1`, choose a nonzero coset representative

```text
x=t u+s v,  |t|,|s|<=1/2,
```

in a centered fundamental parallelogram of `Zu+Zv`.  Neither coefficient is
zero: a rational point on the line of a primitive lattice vector is integral
only at an integral multiple.  Since `|t|+|s|<=1`, central convexity puts
`x` in `K`.

If `x` is live, primitive reduction gives another live direction and

```text
0<|det(u,x)|=|s|D<D,
```

contradicting minimality.  If `x` is dead, then `x in H`.  When
`|s|<=|t|`, put

```text
y=u-sign(t)x=(1-|t|)u-sign(t)s v.
```

The absolute coefficient sum is `1-|t|+|s|<=1`, so `y in K`; moreover
`y=u (mod H)`, so it is live, and

```text
0<|det(u,y)|=|s|D<D.
```

When `|t|<|s|`, the symmetric repair
`y=v-sign(s)x` has the same properties and determinant `|t|D`.
Primitive reduction again preserves the nonzero color and only decreases the
determinant.  Both cases contradict minimality, so `D=1` in `Gamma`.

In the original first-two coordinates this is precisely a pair with
determinant magnitude `3c`.  The subgroup character is essential: subtracting
a dead representative preserves a live coset.  An arbitrary deletion set
does not have this property.

## 3. Exactly three directions: complete classification

Assume there are exactly three unoriented primitive live directions.  Choose
the native basis above and orient its vectors `u,v` so that
`tau(u)=tau(v)=1`.  Write the third primitive direction as

```text
z=m u+n v,  gcd(m,n)=1,  m+n !=0 (mod 3).
```

At least one of `|m|,|n|` equals one.  Otherwise:

- If `mn>0`, orient `z` so `m,n>=2`.  The live primitive vector `u+v` lies
  in the convex hull because

  ```text
  u+v=[z+(n-1)u+(m-1)v]/(m+n-1).
  ```

  It is a fourth direction.

- If `mn<0`, orient `z=m u-r v` with `m,r>=2`.  Coprimality gives `m!=r`.
  If `m>=r+1`, then the fourth live direction `2u-v` lies in the hull:

  ```text
  2u-v=[2z+2(r-1)u+(m-r-1)(-v)]/(m+r-1).
  ```

  If `r>=m+1`, the symmetric hull identity puts `u-2v` inside:

  ```text
  u-2v=[2z+(r-m-1)u+2(m-1)(-v)]/(m+r-1).
  ```

  Each candidate has nonzero owner color and cannot be parallel to `z` under
  `m,r>=2` and coprimality.

Swap the basis if necessary and reverse `z` to write `z=m u+v`.  On the
segment from `v` to `z`, every `(k,1)` lattice coordinate between `0` and
`m` is in the complete convex body.  Liveness deletes exactly
`k=2 (mod 3)`.  Distinctness excludes `m=0`; `m=-1` and `m=2` are themselves
dead; `m>=3` forces the extra live direction `(1,1)`; and `m<=-3` forces
`(-2,1)`.  Therefore the only possibilities are

```text
z=u+v               or               z=-2u+v,
```

up to the declared swaps, signs, and unoriented rays.

## 4. Reciprocal tail

Every direction in a multi-direction body has `l1>=16` by THM-4386, hence
`M_d=max_i|d_i|>=7`.  In either normal form every pairwise first-two
determinant has magnitude at least `3c`.  Therefore, with `P=3c/2`,

```text
M_i M_j>=P for every pair.
```

Let `m=min_i M_i` and `R=sum_i 1/M_i`.  The other two maxima are at least
`max(m,P/m)`.

- If `P<=m^2`, then `R<=3/m<=3/sqrt(P)`.
- If `P>m^2` and `P>=14m`, then

  ```text
  R<=1/m+2m/P<=1/7+14/P,
  ```

  because the difference from the last expression is
  `(m-7)(14m-P)/(7mP)<=0`.
- If `m^2<P<14m`, put `rho=P/m^2`.  Since `m>=7`,
  `1<rho<14/m<=2`; and

  ```text
  R<=(1+2/rho)/m,  (rho+2)^2<=9rho,
  ```

  so again `R^2<=9/P`.

Thus always

```text
R<=1/7+14/P     or     R^2<=9/P=6/c.
```

The exact three residue-deleted ray counts give

```text
E_i < (12/49)R+12/(7c), every i.
```

In the linear branch this is at most `12/343+4/c`, whose value at `c=99`
is `2560/33957<6/77`.  In the radical branch it is at most
`(12/49)sqrt(6/c)+12/(7c)`.  At `c=99` the remaining positive comparison
squares exactly to

```text
96/26411 < 4/1089.
```

Both envelopes decrease with `c`, so `c>=99` is closed analytically.

## 5. Independent H99 replay

The no-import middle-coordinate verifier independently enumerates the whole
`c<99` universe and obtains

```text
primitive eligible rows:        5,409
multi-direction rows:            3,500
exactly-three-direction rows:    1,791
normal form -2:                    684
normal form  1:                  1,107
```

Every selected row satisfies all three projection bounds.  The maximum of
the three projections over the entire head is

```text
18/301 at w=(5,37,43),
E=(240/11137,18/301,2822/55685).
```

Wide controls independently reproduce the six-carrier `(19,23,29)` row,
the 16-carrier exactly-three-direction `(5,191,199)` row, and the
13-direction frontier `(7,611,613)`.

Artifacts:

```text
04-computation/lrc14_colored_basis_three_ray_h99_independent_referee.py
  SHA256 227576fef415844edb01dd2ec66873caa4de77102059322548f6f69ad46fedbb
05-knowledge/results/lrc14_colored_basis_three_ray_h99_independent_referee.out
  SHA256 7dc3cfe799cebbf3b4337a6114c008d60fb0fa566380dc3e16b4d87f6e5160eb
semantic rows
  SHA256 6ecbb5a2e079338af842a2f5d6c474ccdb8eb359941d9a8b192490d677115729
```

Normal, optimized, and frozen output are byte-identical; `65,789` explicit
checks pass.  This establishes a promotion-ready proof candidate for the
exactly-three-direction network closure.  Four or more directions, entry,
synchronization, and LRC(14) remain open.
