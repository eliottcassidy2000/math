# Three coordinate cubics, one sign class

**Status: proved synthesis from THM-2546, THM-3438, THM-3448, THM-3508,
THM-3517, THM-3519, THM-3531, THM-3532, THM-3533, THM-3535, THM-3537, and
THM-3538; no
arbitrary-map classification claim.**  The fixed sporadic tower now has full iterated
wreath monodromy, all-level constant-linear primitivity, an intrinsic
discriminant class, and a reduced newest-prime different at every level.  The
explicit THM-3448 cyclic weighted family is proved in every grade.  The
unstored historical family formerly attributed to THM-1605 remains
unidentified under MISTAKE-416.

**Closure update, 2026-08-16:** THM-3508 proves that all three source
coordinates of the fixed level-two composite are primitive with common class
`[H]`; THM-3519 proves the level-three class `[-2J]`.  THM-3517 proves that
all three actual coordinates are primitive in every explicit cyclic weighted
grade.  Its exact `m=3` quintic atlas has common class `[L5]` but Jelonek set
`V(C) union V(L5)`, so the sign quotient misses a genuine component with
even `C3` inertia.  THM-3531 then closes the basis-independent fixed-tower
class at every level, and THM-3533 upgrades its newest odd divisor from a
parity signal to exact normalization-discriminant multiplicity one.
THM-3535 finally proves that every nonzero rational constant linear form,
including each literal coordinate, is primitive simultaneously at all
levels.  THM-3538 separates primitivity from local maximality by an exact
internal-factor/resultant carrier and proves index zero for `y,z,1/x` at the
newest primes through levels two--four.

## 1. The three discriminants collapse to one square class

Let

```text
K = Q(a,b,c)
```

be the generic target function field of the fixed sporadic Keller map.  The
THM-2546 referee computation proves, for its three coordinate cubics,

```text
Delta_x = -4 A_x^2 L,
Delta_r = -4 A_r^2 L,
Delta_z = -4 A_z^2 L,                                (1)
```

with the three explicit nonzero polynomial square factors recorded there.
Since `4` and every `A_i^2` are squares in `K`, equation (1) gives

```text
[Delta_x]=[Delta_r]=[Delta_z]=[-L] in K*/K*2.          (2)
```

Equivalently, on the generic locus,

```text
K(sqrt(Delta_x)) = K(sqrt(Delta_r))
                 = K(sqrt(Delta_z))
                 = K(sqrt(-L)).                       (3)
```

Thus the three coordinate eliminants have one common quadratic sign
resolvent.  Their finite affine odd-valuation carrier is the same Jelonek
divisor `V(L)`.  A statement about the complete branch divisor would also
have to audit valuations at infinity.

The equality `-4=-(det JF)^2` explains why the constant disappears in square
class.  It does not turn (2) into a degree-four statement: `4` is a square
coefficient, while the surviving information is the binary class `[-L]`.

## 2. Exactly what the common class remembers

For a separable cubic, its discriminant square class is the sign character of
the permutation action on its three roots.  Consequently (2) synchronizes
the three sign quotients.  It does **not** identify:

- the three cubic root sets;
- their full splitting fields or labelled `S3` actions;
- the rational maps between coordinate roots;
- boundary-effective polynomial sections; or
- the multiplicities of the discriminant divisor.

Different cubic extensions can share the same quadratic resolvent.  The
lawful picture is therefore a common-source cospan,

```text
                         x-coordinate cubic
                       /
[-L] sign class  ----  r-coordinate cubic
                       \
                         z-coordinate cubic,          (4)
```

not three pairwise identifications.

This also identifies the correct tournament boundary.  There is no intrinsic
pairwise orientation among the three coordinate cubics: all three expose the
same binary observable.  A three-vertex tournament would manufacture order
that (2) does not contain.  The useful finite carrier is one binary hub with
three views, plus the sidecars lost by the sign quotient.

## 3. Why the `V4`/four-object analogy stops

A generic `V4` quartic-resolvent package has three nonzero binary characters.
The exact effectivity hostile in
`jc-torsor-killing-versus-boundary-effectivity-exact-hostiles-codex-20260814.md`
shows that even all three Kummer characters and their class relations can be
present while natural rational views fail to be polynomial.

Equation (2) is a different phenomenon: the three displayed cubic views land
on the **diagonal**

```text
([-L],[-L],[-L]) in (K*/K*2)^3.                       (5)
```

So neither the three nonzero characters of `V4` nor a tournament on four
vertices may be read off from the three discriminants.  Recovering such a
four-object carrier would require an explicit resolvent map and an oriented
boundary-valuation sidecar.  The numerical coincidence with `-4` supplies
neither.

## 4. Composition proves that the class is dynamical

THM-2576 proves for the fixed map

```text
S_(F o F)=V(LH),                                      (6)
```

with two distinct irreducible Jelonek components.  THM-2582 then proves that
the degree-nine **x-coordinate** eliminant of `F o F` has square class

```text
[Disc(E_2)]=[H],                                      (7)
```

not `[LH]`: the old `L` contribution cancels in parity between the norm
denominator and the odd-block alternating factor.  Hence the atom-level
class `[-L]` is not a static label that can simply be multiplied along the
Keller monoid.  At the first composite rung the sign quotient selects the
new image component instead.

THM-3508 closes the other two coordinate views: all three source coordinates
are primitive in the same degree-nine trace algebra, hence all three have
class `[H]`.  THM-3519 repeats the primitive-element gate at level three and
proves

```text
[Disc(m_x)]=[Disc(m_y)]=[Disc(m_z)]=[-2J]             (7a)
```

for the fixed degree-27 extension.  Thus the proved fixed-tower sequence is

```text
degree 3:   [-L],
degree 9:   [ H],
degree 27:  [-2J].                                    (7b)
```

At these three levels, primitivity turns the coordinate power bases into
bases of one trace space; their discriminants differ by basis-change squares.
The visible unit `[-2]` at level three is retained under MISTAKE-413
discipline.

THM-3531 now separates that coordinate issue from the intrinsic tower.  If
`K_n/K_0` is the degree-`3^n` extension and `P_j` is the raw image-prime
tower, its basis-independent trace-discriminant class is

```text
d(K_n/K_0)=[(-1)^n P_(n-1)]                           (7c)
```

for every `n>=1`.  Thus the newest prime is the sole odd horizontal divisor
at every depth.  In the named normalizations the next rows are

```text
n=4: [2G],       n=5: [-2R_5],       n=6: [2R_6].     (7d)
```

Equation `(7c)` is not a level-four theorem for literal `x,y,z`: a
nonprimitive coordinate can have zero eliminant discriminant even when the
trace form is nondegenerate.  It is also not formal multiplication of the
earlier classes; trace-discriminant transitivity, odd block degree, and the
odd raw clearing exponent produce the newest-prime rule.

THM-3533 supplies the positive-multiplicity refinement missing from square
classes.  In the normalization discriminant, `P_(n-1)` has coefficient
exactly one for every `n`.  For any integral primitive element `theta`,

```text
v_(P_(n-1))(Disc(theta))=1+2 length(B_n/A[theta]).      (7e)
```

Locally maximal primitive elements attain one, but primitivity alone does
not: already `u` and `Lu` at level one have valuations `1` and `7`.  Thus the
square class, normalization different, and coordinate order are three
strictly different layers of the passport.

THM-3535 closes the remaining primitivity quantifier without countable
avoidance.  The newest prime isolates one transposition on one bottom ternary
block; conjugation and leaf transitivity force the full iterated wreath group

```text
Mon(F^n)=S3 wr ... wr S3.                               (7f)
```

Its only blocks are ancestor cylinders.  Two exact three-point fibres then
show that no nonzero constant linear form descends one step, so every such
form generates `K_n/K_0` for every `n`.  Consequently the literal `x,y,z`
minimal polynomials have the intrinsic class `(7c)` at all levels.  This
still does not make their power orders locally maximal: THM-3533's index term
in `(7e)` remains a separate coordinate question.

THM-3537 begins the old-prime effective ledger at depth two.  The old factor
`L` has normalization coefficient four and geometric inertia

```text
(4)(2)(1)^3.                                           (7g)
```

The proof is a parity squeeze: transitivity gives `3+c`, the finite-sheet
unit removes one possible contribution, a tame cubic gives `0<=c<=2`, and
the absence of `L` from the odd class forces `c=1`.  On the canonical
transverse DVR, the literal `x_2` order has exponent eight and index length
two.  Thus the passport now contains one exact positive old multiplicity and
one prescribed-coordinate defect, but not their all-level recurrences.

THM-3538 identifies the missing prescribed-coordinate invariant at every
newest prime.  After the finite-etale predecessor cover is strictly split,
let `f_(theta,i)` be the last-step cubic blocks for `theta=y,z`, or for the
integral reciprocal coordinate `theta=u=1/x` when its chart gate holds.  Then

```text
C_(theta,n)
 =product_i h_theta(q_i)
  product_(i<j) Res(f_(theta,i),f_(theta,j)),
length(B_n/A[theta_n])=v(C_(theta,n)).                 (7h)
```

The first product detects collisions internal to one cubic; the second
detects equal residue values on distinct ancestry blocks.  Local maximality
is therefore equivalent to one forced dual-number shadow and otherwise
squarefree, pairwise-coprime reduced blocks.  Exact characteristic-zero/good-
reduction witnesses prove

```text
i_(y,n)=i_(z,n)=i_(u,n)=0          for n=2,3,4.         (7i)
```

This does not make literal `x_n` integral.  With `N=3^n` and reciprocal
index `I`, its monic discriminant exponent is `3-2N+2I`, while the integral
reciprocal and primitive raw-cleared reversal both have exponent `1+2I`.
Thus the three cubics share one field/sign passport, but a chosen coordinate
still needs its chart and resultant carrier.  Unitness of `(7h)` for every
`n>=5` remains open.

THM-3532 fixes the equivalence relation under which this tower passport is
covariant.  For a two-sided change `G=tau o F o sigma`, every one-step object
transports exactly, but iteration inserts the middle automorphism
`sigma o tau`.  The complete tower therefore transports only under honest
conjugacy `phi o F o phi^(-1)` unless a separate intertwiner is proved.
W1/W2 are target postcompositions, and their explicit second iterates are the
hostile counterexamples to naive tower covariance.  Under nonlinear
conjugacy, a primitive linear observation becomes `ell o phi^(-1)`, generally
polynomial; only affine conjugacy preserves the standard linear-observation
family.  Thus “move within a grade” is weaker than “preserve the iterated
passport.”

## 5. Relation to the explicit cyclic weighted family

MISTAKE-416 repairs the provenance boundary: the legacy THM-1605 record has
no recoverable literal formula for its historically reported `E_m`, and its
degree string cannot define a lawful test object.  THM-3517 instead uses the
explicit cyclic weighted family `E_m^cyc` from THM-3448.  It has

```text
det J(E_m^cyc)=1,
generic degree=2m-1,
global monodromy=S_(2m-1),
coordinate degrees=(10m-13,10m-14,4).                 (8)
```

The first-coordinate degrees begin `7,17,27,...`, not the historical
`7,13,26,...` string.  The common `m=2` numerical degree does not identify
the constructions.

THM-3517 proves a strong all-grade statement for the explicit family:

```text
K(x)=K(y)=K(z)=K(w)                                   (8a)
```

for every cyclic grade.  The `z` gate is the nonzero remainder coefficient

```text
((ell+1)^3+1)P^2
 +(4ell^2+ell-2)P/2+4(ell+1)Q.                        (8b)
```

Together with global `S_n` monodromy and the already proved `x,y` gates,
this makes all three coordinates primitive.  Trace-form congruence therefore
puts all three coordinate eliminants in one inverse discriminant square
class in every grade.  The common-class phenomenon is not special to cubic
degree.

At the first beyond-cubic member `m=3`, the inverse is

```text
T(w)=w^5-w^4+Pw-Q,                                   (8c)
```

with `S5` monodromy and irreducible discriminant

```text
D5=256P^5-27P^4-36P^3Q-50P^2Q^2
   -2500PQ^3+3125Q^4+256Q^3.                         (8d)
```

The `x,y,z` eliminants are irreducible quintics with `17,29,191` terms and
all have class `[D5]`.  After `P=BC,Q=AC^2`,

```text
D5(BC,AC^2)=C^4L5,
S_(E_3^cyc)=V(C) union V(L5),
[Disc_x]=[Disc_y]=[Disc_z]=[L5].                      (8e)
```

The missing `C=0` component has three escaping sheets and local inertia a
3-cycle, which is even as a permutation.  The `L5=0` component has
transposition inertia and remains visible.  This is the exact classification
boundary: the common sign class records parity of local monodromy, not the
full effective Jelonek set.  THM-3517 proves the same sign-blind `C`
component for every odd cyclic grade `m>=3`.

A serious classification record therefore needs at least

```text
degree/fibre grade
+ fixed/free involution orbit profile
+ full monodromy group and block systems
+ quadratic sign-resolvent class
+ Jelonek divisor and its image tower
+ oriented boundary valuations/effective cone.        (9)
```

The first and second entries can agree while the remaining four differ.
Conversely, equality of the sign class alone is far too coarse to identify a
map: it retains only one `C2` quotient of the root permutation data.  The
fixed compositional tower and the explicit weighted family now supply two
proved laboratories with the same trace-form mechanism but different degree,
monodromy, and effectivity ledgers.

## 6. What is exactly classified, and what is not

THM-3438 gives a complete classification of **degree values** in dimension
at least three:

```text
KDeg(m)={1} union {3,4,5,...}       for every m>=3.      (10)
```

Degree two is impossible, and the explicit weighted lift realizes every
degree `n>=3` with global monodromy `S_n`.  These maps are composition atoms:
an imprimitive nontrivial composition cannot have primitive `S_n` monodromy.
For every factorization by allowed degrees, compositions provide additional
maps in the same numerical grade with visible block systems.  At degree
nine, for example, an atomic weighted `S9` map and the fixed composite
`F o F` with full imprimitive monodromy `S3 wr S3` have the same degree but
cannot be equivalent under source/target automorphisms.

This settles the semigroup of attainable integers, not the moduli problem.
An exact classification of maps within one degree would still have to sort
at least the passport `(9)`.  The three cubic equations are one row of that
passport: they prove a common sign quotient for one extension, not that every
degree-three Keller map is conjugate to the sporadic example.  Stable
extension carries all of these three-dimensional examples to every higher
dimension, so it enlarges the ambient dimension without changing this
classification boundary.

## 7. Connection contract and decisive next tests

| field | exact content |
|---|---|
| source | three coordinate cubic discriminants of the fixed sporadic map |
| target | the single class `[-L]` in `K*/K*2` |
| map | take discriminant and quotient by squares |
| preserved | root-permutation sign and finite affine odd divisor `V(L)` |
| destroyed | labelled roots, full monodromy, multiplicities, poles at infinity, effectivity |
| required sidecar | splitting-field/block data plus oriented boundary valuations |
| cheapest hostile | two nonisomorphic cubic extensions with the same discriminant square class |

The square-class and normalized newest-prime tests are closed intrinsically
at every fixed-tower level.  The prescribed-coordinate index has an exact
all-level criterion and equality through levels two--four, while the explicit
weighted family is closed for all three actual coordinates.  The next exact
tests are:

1. prove or refute unitness of THM-3538's internal-factor/resultant carrier for
   `y_n,z_n,u_n` at every newest level `n>=5`; simultaneous primitivity alone
   cannot decide it;
2. extend THM-3537's exact `delta_2(L)=4` to old-prime multiplicities at
   depth three and beyond; the newest prime is already exactly reduced by
   THM-3533;
3. retain inertia cycle type component-by-component, so even cycles cannot
   disappear behind the sign quotient;
4. classify the explicit cyclic family within a fixed degree only after
   adding its Jelonek, monodromy, block, and boundary-effectivity data; and
5. keep the unstored historical THM-1605 family separate unless an explicit
   source formula is recovered.

The strongest current conclusion is now broader but still sharply bounded:
the fixed tower has full wreath monodromy, all-level literal-coordinate
primitivity, one intrinsic newest-prime sign character, a reduced newest
different at every depth, newest-coordinate index zero through levels two--
four for `y,z,1/x`, and the first exact old-prime inertia row
`(4)(2)(1)^3`.  Every literal coordinate shares the intrinsic character at
every level, as do the three coordinates in every explicit cyclic weighted
grade.  The degree-value spectrum is completely classified, but the maps are
not.  The `m=3` sign-blind component `C=0`, the newest/old index split, and
the atomic/composite degree-nine pair are three decisive witnesses to the
missing classification data.
