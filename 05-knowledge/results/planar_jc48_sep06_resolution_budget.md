# A linear resolution test and the first surviving sextic boundary

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED. Nori's theorem is
CITED; explicit curve identities are FINITE-EXACT. No realization is
asserted at the equality boundary. JC(2) remains OPEN.** September 6, 2026.

**Forward status:** the equality curve constructed here is now
[excluded by marked infinity plumbing](planar_jc48_sep06_boundary_plumbing.md),
with an [independent exact braid proof](planar_jc48_sep06_boundary_braid.md).
The strict numerical theorem and frozen equality calculation are unchanged.
The [new three-node curve](planar_jc48_sep06_next_infinity.md) has margin-2
and passes the marked infinity test.

## 1. Recover the global obstruction before choosing another curve

The [audited global curve note](planar_jc48_sep06_global_curve.md) excludes
the former target `(t^4+t^2,t^6+t^5+t^2)`: its actual affine complement
has abelian fundamental group. The closest proved mechanism is Nori's
Proposition 3.27, applied after resolving the finite cusp **and the curve
together with the line at infinity**. The hostile is the
[alternating semilocal passport](planar_jc48_sep06_alternating.md), which
passes all its stated local and transitivity conditions but does not give
a representation of that complement. The corrected near miss is to
treat the abstract passport as a geometric survivor. The least-used
sidecar now is the canonical intersection number on the resolved target.

The concept board is: actual complement; resolution multiplicities;
adjunction; odd-cusp sheet degree; and the distinction between a strict
inequality and its equality boundary. This note rewrites the **classical
Nori criterion**, rather than claiming a new fundamental-group theorem.
It then constructs and verifies a different sextic exactly at the point
where that criterion stops deciding the problem.

## 2. The criterion in terms of the sum of multiplicities

Let `C subset A2_C` be an irreducible reduced curve. Write `e` for the
degree of its projective closure `Cbar`, and `g` for the genus of its
smooth projective normalization. Choose a finite sequence of point blowups

```text
sigma:X -> P2
```

supported over the line at infinity and finitely many points `T subset C`.
Let `D` be the strict transform of `Cbar` and let `E` be the reduced union
of the total inverse image of the line at infinity and **all** exceptional
curves over `T`. Require:

- `D` has only ordinary nodes, with exactly `N` of them;
- `D` meets `E` transversely at smooth points of `D` and `E`;
- all finite centres lie over `T`, and every point of `T` was blown up.

These are properties of the actual simultaneous resolution, not of a
Puiseux list with its line at infinity forgotten. At the `j`th centre let
`m_j` be the multiplicity of the current strict transform of `Cbar`; a
centre off that transform contributes zero. Set `B=sum_j m_j`.

Then the exact identity is

```text
D^2 - 2N = 3e - 2 + 2g - B.                              (1)
```

In particular,

```text
B < 3e-2+2g  =>  pi_1(A2_C\C) is abelian.                (2)
```

Such a curve cannot be the **whole** nonproperness support of a
nonautomorphic polynomial Keller map. For a rational support, a necessary
condition is therefore `B>=3e-2` for **every** resolution satisfying the
listed conditions. Existence of just one strict-inequality resolution
excludes it. The sum `B` depends on the chosen resolution: unnecessary
extra blowups may increase it. No intrinsic minimum is computed here.

**Proof of (1).** Initially `K_P2.Cbar=-3e`. At a blowup with multiplicity
`m`, the canonical divisor is `b^*K+F` and the strict transform is
`b^*D-mF`. Since `F^2=-1` and pullbacks have zero intersection with `F`,
the canonical intersection increases by `m`. Thus

```text
K_X.D=-3e+B.
```

Adjunction gives `D^2+K_X.D=2p_a(D)-2`. A nodal irreducible curve with
normalization genus `g` has `p_a(D)=g+N`. Substitution proves (1). This
also explains why the number of nodes cancels in the linear expression;
discarding them before adjunction would not justify that cancellation.

**Proof of (2) and its Keller consequence.** By construction,

```text
X\E isomorphic to A2\T,
X\(D union E) isomorphic to A2\C.
```

Complex affine two-space minus finitely many points is simply connected
(a contracting disk can be perturbed away from finitely many points in
real dimension four). If the right side of (1) is positive, Nori's
Proposition 3.27 makes the kernel of
`pi_1(X\(D union E))->pi_1(X\E)` abelian. Its target is trivial, proving
(2). The primary statement is in Madhav V. Nori,
[*Zariski's conjecture and related problems*](https://www.numdam.org/article/ASENS_1983_4_16_2_305_0.pdf),
*Ann. Sci. École Norm. Sup.* (4)16 (1983), Proposition3.27, journal p331;
the full statement, proof passage, and surface hypotheses were read.

For the Keller implication, use the actual-cover argument of
[global curve, §4](planar_jc48_sep06_global_curve.md). The monodromy group
of the connected covering over `A2\C` is transitive. A generic smooth
meridian fixes an actual affine sheet; that sheet exists because the
pullback of a defining irreducible polynomial of `C` is a nonconstant
polynomial, whose zero divisor cannot have zero-dimensional image under
a quasi-finite map. In an abelian transitive permutation group, an element
fixing one label fixes every label. Since these meridians normally
generate the complement group, its monodromy is trivial, forcing degree
one and hence an automorphism. This contradicts nonempty nonproperness.
The mapping degree used here is unrelated to the plane-curve degree `e`.

## 3. The original curve and a genuine strictness hostile

For the excluded five-node sextic the finite multiplicities are
`(2,2,1,1)` and the simultaneous infinity multiplicities are
`(2,2,2,1,1)`. Therefore

```text
e=6, g=0, B=14,   D^2-2N=18-2-14=2>0.
```

For a rational curve whose only nonnodal singularities are one finite
`(2,m)` cusp and an infinity `(2,n)` cusp, with odd `m,n>=3`, the standard
embedded cusp resolutions have respective multiplicity sums `m+1,n+1`.
If those sequences already resolve the pair with the line at infinity,
(1) becomes `3e-4-m-n`. The explicit simultaneous-resolution condition
is essential; an extra contact with the line can require more blowups.

Strict inequality in (2) cannot be replaced by equality. The ordinary
affine cusp `(t^2,t^3)` has projective degree three. Resolve its finite
cusp with `(2,1,1)` and its smooth infinity branch, tangent to the line
at infinity to order three, with `(1,1,1)`. Thus `B=7=3e-2`, `N=0`,
and `D^2=0`. Nevertheless its affine complement has the trefoil group
`<s,t | sts=tst>`, with the nonabelian quotient `S3` obtained by
`s=(12),t=(23)`. Weighted radial contraction identifies the affine
complement with the link complement times a contractible radial factor;
the marked `(2,3)` presentation is the specialization of Shimada's
[primary local topology notes, §6](https://www.math.sci.hiroshima-u.ac.jp/shimada/LectureNotes/LNZV.pdf)
already checked in the odd-cusp supplier. The source verifies the
noncommuting permutations and the braid relation. This is a hostile to
**abelianity at equality**, not an example of a Keller map.

## 4. A new explicit curve at equality

Consider

```text
U(t)=t^4+t^3+t^2,
V(t)=2t^6+3t^5+2t^3+2t^2.                                (3)
```

This is a finite birational normalization onto a rational sextic, with
one affine `(2,5)` cusp, **four** ordinary nodes, and an infinity `(2,9)`
cusp. Its simultaneous resolution has

```text
finite:   (2,2,1,1),
infinity: (2,2,2,2,1,1),
B=16,  D^2=8,  N=4,  D^2-2N=0.                          (4)
```

Thus the strict Nori test does not decide its complement. Applying the
audited [odd-cusp spectrum](planar_jc48_sep06_odd_cusp.md), any Keller map
having (3) as its whole nonproperness support would still have degree
six, three generic retained sheets, one actual cusp point, and every
node omitted. Applying the [alternating reduction](planar_jc48_sep06_alternating.md)
would give `A6` and a single boundary prime of index three and residue
degree one. None of this supplies a full complement representation or
an affine source. Both remain **OPEN** for this new target in this note.

Here is a complete exact singularity check. The derivative gcd is `t`,
and `gcd(U,V)=t^2`. Set `W=V-2U+2U^2`; then

```text
W=7t^5+8t^6+4t^7+2t^8.
```

So zero is the unique critical parameter, its image is unshared, and
its analytic branch is `(2,5)`. Monicity of `U` proves finiteness.
For the divided differences `N(s,t),M(s,t)` of `U,V`, put `p=s+t,q=st`.
They give the complete symmetric pair equations

```text
(2p+1)q=p(p^2+p+1),
M=-p^2 h(p)/(2p+1),
h(p)=p^4+3p^3+5p^2+p-7.                                  (5)
```

The first equation makes `2p+1` a unit: at `p=-1/2` the other side
is `-3/8`. The root `p=0` gives `q=0`, hence only the diagonal pair
at the cusp. For the four roots of `h`,

```text
Disc h=-80787,
Res(h,p(2p^2+3p+4)(2p+1))=610785,
(s-t)^2=-p(2p^2+3p+4)/(2p+1).
```

There are exactly eight off-diagonal ordered pairs. The pair scheme
is finite, so the normalization map is generically one-to-one.
An independent exact ideal for coincident tangent lines is

```text
(N,M,U'(s)V'(t)-V'(s)U'(t))=(s+t^3+t^2+t,t^4).
```

Thus all off-diagonal pairs are transverse. They have distinct node
images, rather than combining into a triple point: modulo `U-a`,

```text
V-b = 4t^3+(2a+5)t^2+a t-3a-b = R.
```

A triple fibre would make `R` divide `U-a`. Multiplying the remainder
by sixteen gives the coefficient equations

```text
r2=4a^2+8a+21=0,
r1=2a^2+13a+4b=0,
r0=-6a^2-2ab-19a-b=0.
```

Using `r1` to eliminate `b`, the last equation implies
`j=4a^3+4a^2-63a=0`. The exact identity

```text
(-76a^2-97a+1617) r2 + (76a+173) j = 33957
```

excludes every triple fibre. It also has an independent Groebner check.
The nonzero cubic already excludes four distinct parameters over one
image. This completes the affine singularity inventory.

For the projective extension, homogenize the coordinates in (3) to
degree six with last coordinate `S^6`. They have no common zero, and
at `S=0` the middle coordinate is `2T^6`. Therefore the image is a
sextic and the unique infinity point is `[0:1:0]`. Put `z=1/t` and

```text
X=U(1/z)/V(1/z),   Z=1/V(1/z).
X=z^2/2+O(z^3),
Z-4X^3+30X^4=(35/16)z^9+O(z^10).                         (6)
```

This gives the single Puiseux pair `(2,9)`. The following explicit
charts also keep the original line `Z=0`. Set
`rho=Z/X^3-4`, `tau=rho/X+30`. The successive pairs are

```text
(X,Z), (X,Z/X), (X,Z/X^2), (X,rho), (X,tau),
(X/tau,tau), (X/tau^2-1/2450,tau).
```

Their parameter orders are `(2,6),(2,4),(2,2),(2,2),(2,1),(1,1),(1,1)`.
Every transition is a point blowup followed by centering at the strict
curve's point. The line at infinity separates at the third blowup:
it is `rho=-4`, away from the branch at `rho=0`. At the fourth blowup
the strict curve is smooth and tangent to the new exceptional. The
last two blowups separate that tangency and the resulting intersection
of exceptional axes. In the final chart `tau=35z+O(z^2)`, so the curve
meets only the final exceptional, transversely.

At the finite cusp use `(x,y)=(U,W)` and the charts

```text
(x,y), (x,y/x), (x,y/x^2), (x^3/y,y/x^2),
(x^5/y^2-1/49,y/x^2).
```

Their orders are `(2,5),(2,3),(2,1),(1,1),(1,1)`. In the last chart
the second coordinate is `7t+O(t^2)`. The same two-axis separation
proves final transversality. These charts prove (4), including the
simultaneous divisor requirement. As an independent global check,
`delta_finite+delta_infinity+N=2+4+4=10`, the arithmetic genus of a sextic.

## 5. Connection contract and reproduction

Source: a resolved pair consisting of a rational target curve and its
actual line at infinity. Target: the integer `D^2-2N`, then the abelian
kernel in Nori's theorem. Map: canonical intersection under blowups and
adjunction, followed by the actual complement isomorphisms. Preserved:
the fundamental group of the affine complement and the retained sheets
of a hypothetical Keller cover. Lost by keeping only the integer:
the group at zero or negative margin, branch attachment, and source
realization. The restoring data are the full divisor and strict
transversality. The ordinary cusp is the decisive equality hostile.

The next test is the **actual global complement group of (3)**, or an
additional affine-boundary obstruction. Repeating the strict-inequality
argument at zero would not be progress. The earlier semilocal `A6`
fixture remains useful as a positive control, but is not that group.
No new META-PATTERNS card is needed: the existing controlled-forgetting
and retained-boundary cards already describe this move.

The standalone [source](../../04-computation/planar_jc48_sep06_resolution_budget.py)
and [frozen output](planar_jc48_sep06_resolution_budget.out) check the
complete pair algebra, critical/collision/triple controls, every displayed
blowup chart, the symbolic all-degree identity (including genus), the
positive inherited sextic, and the nonabelian equality hostile.

```sh
python3 04-computation/planar_jc48_sep06_resolution_budget.py
python3 -O 04-computation/planar_jc48_sep06_resolution_budget.py
```

The finite universe consists of these two sextics and the ordinary cusp,
not an enumeration of sextics or Keller maps. The analytic criterion is
proved by adjunction and the cited Nori theorem, not inferred from the
finite examples. Both modes pass **73 always-active exact gates** and
match the frozen output byte for byte. Pins:

```text
source SHA-256:
47f1e5775e8447ddd199fe6c4f7f358ed3091031cb54cd1f54dd9a0480e2a387
output SHA-256:
d371ee9ccd8e4de711e5f243817bf428da0521e97921c216d8d14fd118e7087b
```

The [independent audit](planar_jc48_sep06_resolution_budget_audit.md)
passes the full analytic argument, source and replay, and independently
reconstructs every chart using standard-library Fraction arithmetic.
Its SHA-256 is
`25ab4e737a49e67c955000854256262868ad0d3772dd979c533f954357a8ff97`.
