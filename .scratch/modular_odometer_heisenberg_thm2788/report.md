# Physical modular odometer versus the Heisenberg extension

**Status:** proof-complete scratch candidate + exact controls; not canon and
not independently audited.

## 1. Statement

Let `p` be prime and put `Omega=Z/p^2 Z`.  On `Omega` define

```text
X(n)=(1+p)n,             Y_phys(n)=n+1,             Z(n)=n+p.
```

Write `n=v+pw` with standard digits `v,w in F_p`, and define the
carry-suppressed low-digit map

```text
Y_0(v,w)=(v+1,w).
```

Then:

1. `M_p=<X,Y_phys>` has order `p^3` and presentation

   ```text
   Y_phys^(p^2)=X^p=1,
   Z=Y_phys^p central of order p,
   X Y_phys X^(-1)=Y_phys^(1+p)=Z Y_phys.
   ```

   Thus `M_p=C_(p^2) semidirect C_p`, the nonabelian modular group of
   order `p^3`.

2. `H_p=<X,Y_0>` has order `p^3` and the same central commutator relation

   ```text
   X Y_0 X^(-1)=Z Y_0,
   ```

   but `Y_0^p=1`.  For odd `p`, every nonidentity element of `H_p` has
   order `p`, whereas `M_p` has `p^2(p-1)` elements of order `p^2`.
   Hence the two groups, and therefore their permutation subgroups on
   `Omega`, are not isomorphic for odd `p`.

3. The two central extensions have the same alternating commutator form.
   With quotient coordinates `(a,b) in F_p^2` and the section
   `s(a,b)=Y^b X^a`, their normalized cocycles are

   ```text
   c_H((a,b),(a',b')) = a b',
   c_M((a,b),(a',b')) = a b' + floor((b+b')/p) mod p.
   ```

   Their difference is exactly the base-`p` carry cocycle, the Bockstein
   of the low digit.  Antisymmetrization kills that symmetric carry term:

   ```text
   c(u,v)-c(v,u)=a b'-a' b.
   ```

   For odd `p`, the `p`-power map on the quotient is zero for `H_p` and
   is `(a,b)->b Z` for `M_p`.

4. At `p=2`, both permutation groups are the same `D_8`.  Explicitly,

   ```text
   Y_phys=(X Y_0)^(-1).
   ```

   This is the sharp boundary: the binary carry is absorbed by a generator
   change, while for every odd prime the exponent obstruction survives.

5. Both groups have minimal faithful permutation degree `p^2`.  The displayed
   actions attain that degree.  Conversely, every subgroup of index `p` is
   normal and contains the commutator `Z`; if all orbits of a permutation
   action had size below `p^2`, `Z` would lie in the kernel of every orbit.

For `p=13`, `Y_phys=Y_0` at the `156` points off the low-digit wall and
`Y_phys=Z Y_0` at the `13` points with `v=12`.  Therefore any local scan that
does not close the low-digit cycle sees the Heisenberg law, while the full
physical successor necessarily sees the modular Bockstein.

## 2. Proof of the two group laws

Multiplication by `1+p` has order `p` modulo `p^2`, and conjugates unit
translation by

```text
X Y_phys X^(-1)(n)=n+(1+p)=Y_phys^(1+p)(n).
```

The translation subgroup has order `p^2`; its intersection with `<X>` is
trivial because every nonidentity translation moves zero while `X` fixes
zero.  Hence `|M_p|=p^3`, `Z=Y_phys^p`, and the displayed presentation
follows.

In digits,

```text
X(v,w)=(v,w+v),       Y_0(v,w)=(v+1,w),       Z(v,w)=(v,w+1).
```

Direct substitution gives `X Y_0=Z Y_0 X`.  Normal forms
`Z^c Y_0^b X^a` are distinct, so `|H_p|=p^3`.  For odd `p`, the standard
class-two power formula has commutator exponent `binom(p,2)`, a multiple of
`p`; hence every element has order dividing `p`.  The group is nontrivial,
so every nonidentity element has order `p`.

For `M_p`, write an affine element as

```text
n -> (1+p)^a n+b,      a in F_p, b in Z/p^2 Z.
```

For odd `p`, its `p`th power is translation by `pb`; it has order `p^2`
exactly when `b` is nonzero modulo `p`.  This gives

```text
order 1:       1 element,
order p:       p^2-1 elements,
order p^2:     p^2(p-1) elements.
```

The order spectrum separates `M_p` from `H_p` for every odd prime.

## 3. The carry is the missing extension coordinate

Use the section `s(a,b)=Y^bX^a`, with `a,b` represented in
`{0,...,p-1}`.  In either group the commutator contributes `ab'` when
moving `X^a` past `Y^(b')`.  In `M_p`, the exponent `b+b'` of the physical
translation also crosses `p` exactly when the two low digits carry, adding
one central `Z=Y_phys^p`.  This proves the two cocycle formulas.

The carry term is symmetric in `b,b'`, so both extensions have the same
alternating form.  This explains why a commutator-only, determinant-only,
or three-point central-direction check cannot distinguish them.

For odd `p`, raise a lift `Y^bX^a` to the `p`th power.  The class-two
commutator term again vanishes modulo `p`.  In `H_p` both generators have
order `p`, so the result is one.  In `M_p`, it is `Y_phys^(pb)=Z^b`.
Thus the power/Bockstein coordinate, not the symplectic commutator, is the
first distinguishing invariant.

When `p=2`, direct digit calculation gives
`Y_phys=(XY_0)^(-1)`, so the two generated subgroups coincide.  Their order
spectrum is `1^1 2^5 4^2`, the spectrum of `D_8`.

## 4. Minimal faithful carrier

Let `G` be either nonabelian order-`p^3` group above.  Its commutator
subgroup is `<Z>` of order `p`.  Every subgroup of index `p` is normal, and
the quotient has order `p`, so that subgroup contains `[G,G]=<Z>`.

Suppose a faithful permutation action had no orbit of size at least `p^2`.
Every point stabilizer would then have index one or `p`; its core would
contain `Z`.  The intersection of the orbit kernels would still contain
`Z`, contradicting faithfulness.  Hence the degree is at least `p^2`.
Both displayed actions on `Omega` are faithful, proving equality.

This agrees with THM-2779's root-degree obstruction but refines its typing:
the common `p^2` carrier size does not identify the two extension classes.

## 5. LRC(14) consequence and stopping boundary

THM-2779's endpoint action uses `X,Y_0,Z`; THM-2782's physical arms move in
the central `Z` direction.  Since `Z=Y_phys^p` is common to both extension
classes, THM-2782's three-point central segment is compatible with this
theorem and cannot detect the carry Bockstein.

A complete physical successor cycle is different.  In standard digits,

```text
Y_phys(v,w)=Y_0(v,w)             for v!=p-1,
Y_phys(v,w)=Z Y_0(v,w)           for v=p-1.
```

Thus for odd `p` no full-cycle current intertwiner can simply identify the
actual physical successor with THM-2779's carry-suppressed generator while
retaining the same group law.  A positive construction must instead retain
the carry/Bockstein cochain, work with the modular extension, or supply a
graded/enlarged carrier on which the power defect is represented.

This is a group/action obstruction.  It does not prove that a physical
current realizes `M_13`, construct a same-ancestry endpoint/root map, exclude
one of the `165` rows, or prove LRC(14).

## 6. Exact controls

The scratch companion `probe.py` checks primes `2,3,5,7,11,13`.  It
constructs both permutation groups, verifies their orders and relations,
checks every cocycle pair in `F_p^2 x F_p^2`, computes the complete order
spectra, verifies the `p=2` subgroup equality and generator relabel, checks
the `p=13` `156+13` carry-wall census, and contains no Python `assert` node.

Ordinary and optimized transcripts agree.  Current LF hashes:

```text
script bd64ea35ea093519c23b12484e1cf3857a74361d8633b541868916de62f3995b
output 2792a695832dc2c05227813bbf2e3182dce33c251856711146320f564f9ea35e
```
