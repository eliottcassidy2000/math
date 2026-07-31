# The THM-2825 collar is a partial \(V_4\)-action, and its sharp carry base is Witt-nonsplit

Status: **PROVED abstract algebra + FINITE-EXACT scratch application.**
No tracked file was edited.  This note gives no row exclusion and no
LRC(14) conclusion.

The finite application is rebuilt from the hash-pinned physical constructor
behind
[`THM-2818`](../../01-canon/theorems/THM-2818-right-cofiber-positive-copy-stratification-and-alternating-half-step-chains.md).
The exact collar and physical boundaries are the proved
[`THM-2825`](../../01-canon/theorems/THM-2825-nearest-half-step-common-right-collar-and-semantic-parity-boundary.md).
The nonsplit carry comparison uses proved
[`THM-2851`](../../01-canon/theorems/THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar.md).

## 1. Verdict

The punctured square has a canonical mathematical meaning:

> The \(M_3\otimes I_{587}\) collar algebra is the corner
> \(p(M_4\otimes I_{587})p\) obtained by restricting the translation
> groupoid of a free \(V_4\)-torsor to three of its four objects.

As a completion of this specified partial translation action, it has a
unique minimal enveloping \(V_4\)-set up to equivariant isomorphism fixing
the three old objects: add one parity-shifted copy \(X=\Pi R\), of dimension
\(587\), at the missing degree \((1,0)\).  This is not uniqueness among
arbitrary algebraic or Morita completions.  The corresponding four-corner
Koszul square is exact and contractible.  Thus ordinary linear algebra,
mapping cones, Morita equivalence, and scalar cocycle twists can complete
the diagram formally.

They cannot complete it in the inherited physical translation category
while retaining carrier, native-factor, and endpoint predicates as
equivariant object data.  The obstruction occurs before holonomy:

```text
R:     source carrier = empty;   at least one native-factor hole;
M1,M2: source carrier = delta_0; all six native factors present.
```

The empty carrier mask and a delta mask lie in disjoint translation orbits.
A sign or phase multiplied by zero is still zero.  A quotient which
identifies them has discarded the physical presence predicate it was asked
to retain.

There is a second, independent nonsplitting.  The q/carry sidecar of
THM-2851 cannot be the split affine plane \(\mathbf F_{13}^2\), even though
that plane has the desired \(169\) points.  Every \(13\)-primary element of
\(\operatorname{AGL}_2(13)\) has order at most \(13\), whereas the natural
lift \(L_1\) has order \(169\) and satisfies

```text
L_1^13=T != 1.
```

The sharp abstract state object is instead

```text
C_169 = (W_2(F_13),+),
```

the nonsplit extension of the residue coordinate by the ancestry carry.
This is a structural distinction, not another state-count estimate.  It
does **not** say that the physical bank contains no carry quotient.  A
concurrent long-horn audit finds a literal \(T^8\) kernel step inside the E3
block; since \(8\) is a unit modulo \(13\), this is a faithful generator of
the same \(C_{13}\) carry kernel.  What is still missing is its gluing across
the E3/complement and q/address macro boundary, together with the fourth
\(V_4\) object.

## 2. Inheritance

- Closest positive mechanism: THM-2825's exact ladder
  \(R\xrightarrow dM_1\xrightarrow aM_2\), with \(ad=S\).
- Canonical hostile: all \(587\) roots cross
  `source empty -> source delta_0` and a native-factor boundary.
- Corrected near miss: a group/cocycle degree may be added only when the
  arrows are composable on one object; equal degree arithmetic is not a
  fibre product.
- Least-used sidecars: THM-2606's \(V_4\) partial-cube classification and
  THM-2851's nonsplit \(C_{169}\) realization groupoid.
- Incoming positive sidecar:
  [THM-2857](../../01-canon/theorems/THM-2857-endpoint-galois-carry-torsor-and-phase-alignment-sidecar.md)
  gives a free coefficient-field \(C_{13}\) torsor but proves that the
  present ancestry action is scalar-linear rather than semilinear.

The live concept board was:

```text
partial V4 action; strong grading; cocycle twist; partial cube;
Koszul filler; physical mask orbit; Witt carry.
```

## 3. The exact partial transformation groupoid

Write

```text
V4=F_2^2,
U={00,11,01}=V4\{10},

R =00,             M1=11,             M2=01.
```

Translation by \(g\in V_4\) is defined at \(u\in U\) exactly when
\(u+g\in U\).  For \(g\ne0\), its domain has two points.  Every ordered pair
of points of \(U\) has a unique difference, so this partial transformation
groupoid is the three-object pair groupoid.

On one root fibre let \(E_{ij}\) be the matrix unit from corner \(j\) to
corner \(i\), and grade it by

```text
deg(E_ij)=grade(i)+grade(j).
```

Then

```text
dim A_00=3,
dim A_11=dim A_10=dim A_01=2.                         (1)
```

Tensoring the arrows with the \(587\)-dimensional root space gives the
THM-2825 linking algebra.  The rootwise decorated groupoid has the same
support calculation in every root fibre.

This grading is not strong.  Its nine nontrivial ordered products each miss
one matrix unit:

```text
A_11 A_11 misses E_22,       A_10 A_10 misses E_00,
A_01 A_01 misses E_11;                                     (2)

each of the six A_g A_h, g,h distinct nonzero,
spans one orientation of A_(g+h) and misses the other.      (3)
```

Each missing unit has rank \(587\) after tensoring with the root identity.
The audit finds exactly

```text
7 ordered degree products with defect 0,
9 ordered degree products with defect 587.                 (4)
```

In particular, the observed forward path

```text
R --d,(11)--> M1 --a,(10)--> M2
```

supplies one order of the degree-\(01\) product.  The reverse ordering would
need the alternate path through the missing object.

## 4. Universal minimal globalization

Let the full object set be \(V_4\), let \(B=M_4(K)\) carry the translation
grading, and put

```text
p=e_00+e_11+e_01.
```

Then

```text
A=pBp=M_3(K),                                             (5)
```

and the absent diagonal idempotent is exactly \(e_{10}\).  Relative to the
partial algebra, the full globalization adds homogeneous dimensions

```text
degree 00: 1,
each nonzero degree: 2.                                  (6)
```

Thus it adds seven matrix units per root, \(7\cdot587=4109\) rootwise
decorated arrows in total.

This enveloping \(V_4\)-set completion is forced.  A \(V_4\)-set orbit has size
\([V_4:H]\in\{1,2,4\}\); no quotient orbit has size three.  Any completion
which keeps the three existing corners distinct therefore has a free
four-point orbit and one new object per root.  Equivalently, a strong
completion must repair the rank-\(587\) defect in `(2)` and hence cannot add
fewer than \(587\) missing object dimensions.

The corner `(5)` is full and hence Morita-equivalent to \(M_4\).  This is
why ordinary ungraded \(K\)-theory or Morita invariants are blind to the
problem.  The lost information is the domain of the partial action and its
physical object types.

## 5. Twisted torsors and signed Schur multipliers do not add support

A nonzero scalar two-cocycle replaces a supported product by a nonzero
multiple of the same product.  It changes neither homogeneous dimensions
nor the span of \(A_gA_h\).  Therefore all nine defects in `(4)` survive
every scalar cocycle twist.

There is an even sharper finite control for signs.  For the normalized
three-object pair groupoid, exact \(\mathbf F_2\) cochain elimination gives

```text
dim C1=6,      dim Z1=2,
dim C2=12,     dim Z2=dim B2=4,
H2=0.                                                       (7)
```

Thus every object-dependent sign two-cocycle is a coboundary.  If one
restricts the allowable gauge to a translation-invariant gauge, a group
cocycle may remain visible, as in a Pauli/projective representation.  That
does not weaken the support no-go: a projective phase still cannot create
the absent object or an empty carrier value.

The same statement holds for signed Schur multipliers.  Multiplying existing
matrix entries by signs preserves their zero pattern.  It can orient the
formal Koszul square below; it cannot manufacture its fourth physical
vertex.

There is a useful extension-stability form of the same no-go.  Let
\(\widetilde G\to V_4\) be any central extension and decorate every supported
partial translation by a lift in \(\widetilde G\).  After forgetting the
decoration, the degree-\(g\) support is still

```text
R_g={(u,u+g):u and u+g lie in U}.
```

A lifted product can exist only when its two projected arrows have an
intermediate object in \(U\).  Therefore every missing matrix unit in
`(2)--(3)` remains missing after *any* central carry, phase, or sign
extension.  The extension may distinguish two already-composable paths; it
cannot create the absent intermediate object.  In particular the
\(C_{169}\) carry sidecar in Section 10 and the missing \(V_4\) corner solve
orthogonal invoices.  They can be combined only after paying both.

## 6. Partial cube: the composite is a diagonal, not an edge

Choose the two \(V_4\) directions \(11\) and \(10\).  Their Cayley graph on
the full torsor is the \(C_4\) partial cube classified by THM-2606.  Removing
the vertex \(10\) leaves the isometric path

```text
R --d--> M1 --a--> M2.                                    (8)
```

Its two edges are the two Djoković--Winkler coordinates.  The even collar
\(S=ad\) joins the endpoints of this path and therefore has partial-cube
distance two.  Promoting \(S\) to another unit edge makes \(K_3\), which is
not bipartite and hence not a partial cube.

Consequently the unique minimal partial-cube repair is the fourth vertex
and the missing two-edge route.  It does not turn \(S\) into a primitive
translation.  This recovers, from graph geometry, THM-2825's warning that
the \(+2h\) collar factors rather than teleports.

## 7. The minimal Koszul/mapping-cone filler

Let \(H_R\) be the \(587\)-dimensional root space.  Introduce a parity-shifted
copy

```text
X=Pi(H_R)
```

in degree \(10\).  Let \(x:H_R\to X\) be the identity on root labels and
let \(y:X\to H_{M2}\) be the copy of \(S=ad\).  The square commutes:

```text
ad=yx=S.                                                  (9)
```

With the ordinary Koszul orientation, define

```text
D0(r)=(d r,x r),
D1(m,z)=a m-y z.                                         (10)
```

Then

```text
D1 D0=ad-yx=0,

0 -> H_R --D0--> H_M1 direct-sum X --D1--> H_M2 -> 0     (11)
```

is exact.  Indeed \(D_0\) is injective, \(D_1\) is surjective, and
\(\ker D_1=\{(dx,x):x\in H_R\}=\operatorname{im}D_0\).
The exact dimensions and ranks are

```text
587 -> 1174 -> 587,
rank D0=rank D1=587,
all homology zero.                                       (12)
```

This is the universal minimal filler.  Any alternate factorization
\(yx=S\), or \(yx=-S\) under the anticommuting sign convention, satisfies

```text
587=rank(S)<=dim(X).                                     (13)
```

Equality forces a full root-labelled copy.  Thus a mapping cone can close
the algebra, but the result is contractible: it kills the linear defect
rather than producing a physical current.

The underlying arrows in `(9)` may be positive rootwise partial
isometries; the minus sign in `(10)` is the external Koszul orientation.
It is not a new negative physical transport.

## 8. Semantic population gives a separate no-go

The independent bank reconstruction gives

```text
573 triples: (live,dead,live),
 14 triples: (dead,live,dead).                            (14)
```

Therefore no semantic-flipping permutation of the existing \(587\) right
roots exists.  At most the \(14\) dead roots can be paired with \(14\) live
roots:

```text
28 vertices covered,       559 live roots unmatched.     (15)
```

Across the current three corners there are

```text
1160 live, 601 dead.
```

The formal missing copy must contain

```text
14 live, 573 dead,
```

after which the four-corner object has \(1174\) of each semantic colour.
Thus the absolute semantic character defect is \(559\), and the missing
copy cancels it.

The relative \(V_4\) grading gives an even cleaner Fourier invoice.  The
three nontrivial character sums of the punctured object are

```text
+587, -587, +587,
```

and the missing corner contributes the opposite three values.  A free
\(V_4\) torsor has zero nontrivial character sums.

A signed Grothendieck class can record this parity-shifted copy.  It does
not turn a zero delayed value into a nonzero delayed value: semantic
zero/nonzero is an object label, not a coefficient sign.

## 9. Physical typing is an \(H^0\)/orbit obstruction

The all-root reconstruction gives the right-native hole census

```text
E3:       319,
E3+c2:    37,
c2:      217,
q1:       14.                                            (16)
```

Every \(M_1\) and \(M_2\) piece has all six native factors.  More decisively,
the thirteen-twist carrier masks are uniformly

```text
R:     (source,target)=(empty,delta_0),
M1,M2: (source,target)=(delta_0,delta_0).                 (17)
```

Translation fixes the empty mask and permutes the thirteen delta masks.
The two mask orbits are disjoint.  Hence no carrier twist, groupoid
1-cochain, signed multiplier, or projective phase gives even one typed
arrow \(R\to M_i\) while retaining source presence.

The canonical endpoint audit supplies a weaker but compatible obstruction:

```text
74/587 source endpoint pairs are translation-related, only by delta=0;
513/587 lie in distinct source endpoint-translation orbits;
587/587 target endpoint pairs agree with delta=0.          (18)
```

Thus path acyclicity is irrelevant on the \(513\) failures.  A tree makes
an already group-valued edge cochain exact, but it does not put two endpoint
masks into the same group orbit.  The obstruction is object-orbit
membership before \(H^1\), not loop holonomy.

A minimal explicit witness is

```text
R  =[142004190428100,142004216872980),
M1 =[142004591508780,142004617953660),
M2 =[142004992589460,142005019034340).
```

The steps are \(h\) and \(2h\), with \(h=401080680\).  The root already
fails native `E3`, has empty source carrier, and its source endpoint mask
has no present address; \(M_2\) has all six native factors, source
`delta_0`, and \(81\) present source endpoint addresses.  No translation
or sign changes those three facts.

The formal \(X=\Pi R\) preserves the root annotations by copying them.
Its new edge \(X\to M_2\) then inherits exactly the same
empty-to-delta/factor-hole boundary.  The formal completion relocates the
physical obstruction; it does not remove it.

## 10. Why the \(169\)-state sidecar is \(C_{169}\), not \(\mathbf F_{13}^2\)

Write a state of \(C_{169}\) in ordinary base-\(13\) digits as

```text
n=13a+q,       a,q in {0,...,12}.
```

The natural residue lift is

```text
L_h(a,q)
 =(a+floor((q+h)/13), q+h mod13).                        (19)
```

Exact exhaustion of all \(13^2\cdot13^2=28561\)
state/ordered-step instances gives

```text
L_k L_h=T^floor((h+k)/13) L_(h+k mod13),
L_1^13=T,
ord(L_1)=169,       ord(T)=13.                           (20)
```

Every one of the thirteen lifts \(1+13c\) of the residue generator in
\(C_{169}\) is prime to \(169\), hence has order \(169\).  Therefore

```text
0 -> C13 -> C169 -> C13 -> 0                             (21)
```

has no group section.

Now consider a putative affine realization on the split plane
\(\mathbf F_{13}^2\).  If \(g=(A,b)\in\operatorname{AGL}_2(13)\) had order
\(169\), its linear part would be a \(13\)-element.  Such an \(A\) is
unipotent with \((A-I)^2=0\).  Hence

```text
A^13=I,
I+A+...+A^12=13I+binom(13,2)(A-I)=0,
g^13=(I,0).                                               (22)
```

This is a contradiction.

The independent exhaustive control finds

```text
|GL_2(F13)|=26208;
169 linear 13-primary candidates
  = identity + 168 order-13 unipotents;
169*169=28561 affine candidates;
all 28561 are killed by the thirteenth power.             (23)
```

So the distinction is exact:

```text
F13^2 additive:       169 states, exponent 13, split;
W_2(F13) additive:    169 states, exponent 169, nonsplit.  (24)
```

The latter additive group is \(C_{169}\).  A general non-affine
\(169\)-point permutation action could also contain a \(169\)-cycle, but it
has then recreated the same cyclic torsor rather than used the split affine
plane.

## 11. The q3/q11/q7 horn meets the collar only as a formal degree

THM-2851's exact triangle is

```text
(a,3) --L8--> (a,11) --L9--> (a+1,7),
   \----------------L4--------------> (a,7),             (25)
```

so the two routes differ by \(T\).  The affine-address comparison is flat
and forgets this carry.

On the THM-2847 20-cell horn, the outer-word/E3 diagonal has formal
\(V_4\) degree \((1,1)\).  The semantic/source-carrier collar has formal
degree \((0,1)\).  These coordinates have different physical names; only
after a typed sidecar identifies them may one form the formal sum

```text
(1,1)+(0,1)=(1,0),                                      (26)
```

the missing corner degree.

This is a useful locator, not a construction.  The
[concurrent exact hinge scout](../lrc_q3q11_collar_bridge_20260728/cell_path_scout.out)
finds \(52\) collar roots across the twenty horn cells.  Exactly one
distinguished root in each cell—twenty labelled roots in all—uses the same
physical collar triple displayed in Section 9.  All twenty are root index
zero, are live/dead/live, fail only native `E3` at the root, and collapse to
one physical interval triple.  Their path lengths are

```text
14:4, 40:4, 118:4, 144:8.
```

The common pulled allocation atom contains q0 and excludes q3, q7, and q11:

```text
q0:true, q3:false, q7:false, q11:false.                   (27)
```

Thus \(M_2\) is the q0 atom, while the q3/q11/q7 horn atoms are outside
that collar triple.  The two arrows whose degrees are added in `(26)` do
not share the required source/target object.  The twenty labels collapse to
one physical hinge; they do not supply twenty carry fibres.

The companion
[endpoint/carry intertwiner audit](../lrc_q3q11_collar_bridge_20260728/endpoint_carry_intertwiner_audit.opt.out)
first reports that no *normalized* \((a,1)\) carry vertex occurs on a
distinguished rooted path.  Normalization is not an invariant obstruction:
the same audit finds that rooted steps \(2\to68\) have the unique endpoint
translation

```text
(0,8)=T^8.
```

Because \(8\in\mathbf F_{13}^{\times}\), \(T^8\) generates the full carry
kernel.  The strengthened incoming check verifies that the \(+66h\) move
from rooted step \(2\) to step \(68\) translates both source and target
endpoint masks by \((0,8)\).  The literal ancestry label sets agree, the
semantic step is even, and both pieces are common, all-factor, and
carrier-present in exactly the twelve cells

```text
clock=1,  sigma in {0,3,12},  target in {5,6,9,10}.
```

This is a genuine positive physical \(C_{13}\) carry quotient inside the E3
block, not merely a rank shadow.  The full endpoint-mask convolution has
rank \(169\), its \(C_{169}\) vertical derivative has rank \(156\), and
vertical projection has rank \(12\), consistently detecting that charged
direction.  The twelve cell labels collapse to one physical interval edge,
so this is one faithful physical generator with labelled rank \(12\), not
twelve independent carry fibres.

The generator is local, not a completed thirteen-cycle.  Along the long
path the source-frame \(T^8\) plateau is steps \(68,\ldots,90\), while the
target-frame plateau ends at \(89\); step \(90\) is an exact target boundary.
Semantic values alternate, with step \(88\) the last two-frame
semantic-preserving point.  This is the first concrete failure boundary for
iterating the carry edge.

The unit normalization is algebraically harmless.  The central presentation

```text
T^13=1,  [L,T]=1,  L^13=T^8
```

still defines \(C_{169}\): replace the kernel generator by \(T'=T^8\).
Equivalently, the nonzero class \(8\) and the normalized class \(1\) are in the
same automorphism orbit in \(H^2(C_{13},C_{13})\).  The missing datum is the
physical lifted \(L\) and its macro descent, not the choice between \(T\)
and \(T^8\).

THM-2857 writes its free Galois orbit as

```text
c_r=A-B omega^(3r).
```

An arbitrary abstract rechart \(r=5b\) would turn \(b\mapsto b+8\) into a
unit step, but no physical map identifies endpoint \(b\) with Galois \(r\).
Calling that rechart forced would manufacture precisely the missing clutch.
The exact inherited exponent comparison instead gives

```text
b -> b+8,       r=10b,       Delta r=2,
candidate centered phase=omega^(3*2)=omega^6.            (28)
```

The target endpoint current has this phase, but the source current has its
inverse:

```text
(source x-sweep phase, target endpoint-sum phase)
  =(omega^7,omega^6)=(omega^(-6),omega^6).               (29)
```

Their paired phase is one.  Therefore the \(+66h\) edge is a genuine
physical carry generator and a one-sided character match, but not the
THM-2857 semilinear clutch.  The proved ancestry action remains
\(K_0\)-linear, the q7 endpoint mask is empty, and the endpoint-address
arrow is not an allocation-q arrow.  The cheapest next test is a lawful
orientation-sensitive clutch that prevents the inverse source/target
phases from cancelling while preserving the already typed \(T^8\) edge.

This positive result does not create the missing \(V_4\) object and does not
yet glue the q3/q11/q7 states to the q0 collar across the
E3/complement--q/address macro boundary.  It replaces “carry is absent” by
the sharper statement “carry is present on one macro block but has not
descended across the required typed quotient.”

Categorically, the \(T^8\) edge becomes kernel isotropy over the
E3-common macro object.  The degree-\((1,0)\) defect instead asks for a new
base object.  This is the concrete reason a positive carry generator and
the punctured-\(V_4\) no-go coexist without contradiction.

Abstractly one could propagate this isotropy through the three-object pair
groupoid by conjugation, producing the gauge groupoid

```text
U x C13 x U
```

and the coefficient algebra \(M_3(K[C_{13}])\).  The physical conjugation
already fails on the even arrow \(S:R\to M_2\): its endpoint fibres have
cardinalities \(0\) and \(81\), so there is no bijective intertwiner on
which \(S^{-1}T^8S\) could act.  Hence the carry kernel is currently
localized to the common E3 block.  An actual solution must construct the
missing endpoint/current fibre, not merely declare the abstract isotropy
transport.

The first legitimate common base for the *full* q/carry lift in `(26)` is a
nonsplit \(C_{169}\) torsor.  The split \(\mathbf F_{13}^2\) plane can host
the order-\(13\) kernel generator \(T^8\), as the positive edge shows, but it
cannot host the lifted order-\(169\) generator by `(22)--(24)`.  The Witt
sidecar is therefore no longer merely formal: its kernel is physically
visible.  What remains is its descent to the q states and across Section
9's empty/delta carrier and native-factor boundary.

This explains the parallel with THM-2847's rank-one mapping cone.  Its
vector-space exact sequence splits over the cyclotomic field, just as
Section 7's Koszul completion is linearly exact.  The surviving obstruction
is in the realization category:

```text
descent of the local carry band + typed object-orbit + physical basepoint,
```

not in ordinary linear `Ext`, rank, support cardinality, or Fourier
nonsingularity.  THM-2852's full convolutional spectrum cannot repair an
absent object for the same reason.

## 12. Sharp next test

The next object should be a **Witt-decorated physical groupoid**, not another
flat \(13\times13\) array:

1. use \(C_{169}\) states \((a,q)\) with the natural lifts \(L_8,L_9,L_4\);
2. anchor \((a,0)\) at the distinguished q0 collar hinge;
3. attach the full source carrier, six native factors, source endpoint mask,
   and ancestry label to every proposed state;
4. test whether q3, q11, and q7 can inhabit the same physical orbit with
   \(L_9L_8=TL_4\) and \(T\ne1\);
5. require the degree-\((1,0)\) arrow from `(26)` to be an actual composable
   morphism, not merely an equality of grade names.

The distinguished q0 cofiber hinge still fails at step 3: its root has empty
source carrier and the common endpoint has \(81\) present addresses.  But
the long-horn \(2\to68\) edge now supplies a positive \(T^8\) kernel control
entirely inside the E3 block.  The next test is therefore to transport that
proved generator across the E3/complement and q/address macro gluing and ask
whether it extends from the kernel to the full nonsplit \(C_{169}\) lift
while also supplying the missing \(V_4\) object.  An enlarged cofiber state
or genuinely new physical/current carrier is required only at the boundary
where the existing typed edge ceases to extend.

Equivalently, this is now a descent problem rather than a state-count
problem: \(T^8\) supplies the local \(C_{13}\) band, the q3/q11/q7 triangle
supplies its nontrivial two-cocycle, and the missing endpoint/current
intertwiner is the absent descent datum across macro charts.  The fourth
\(V_4\) vertex remains a separate base-object completion.

## 13. Reproduction

Run

```bash
python3 .scratch/lrc_v4_obstruction_20260728/audit.py
python3 -O .scratch/lrc_v4_obstruction_20260728/audit.py

python3 .scratch/lrc_v4_obstruction_20260728/witt_hinge.py
python3 -O .scratch/lrc_v4_obstruction_20260728/witt_hinge.py
```

The scripts use explicit failure gates and contain no Python `assert`
statements.  The main audit independently reconstructs the bank, semantic
population, native factors, and carrier masks.  It reads the endpoint
`74/513` boundary from the hash-pinned canonical transcript rather than
repeating its multi-hour construction.  The Witt companion is
dependency-free.

Final normal and `-O` transcripts are byte-identical.  During validation an
earlier zero-byte run correctly stopped at an audit-only matrix-index bug:
enumerating the full \(V_4\) reordered an old corner before taking a
matrix-unit set difference.  The repair uses
`FULL_GRADES=GRADES+(MISSING,)`, preserving all three old indices; the exact
added dimensions are then `(1,2,2,2)`.  That failed run established no
mathematical claim and is not one of the frozen transcripts.

LF-byte SHA-256:

```text
audit.py          5a13e2cc0c5afc07d5a2df382f19aa55b245ac4bb103195a324abaf76c815e92
audit.out         4cd1e72050faa749cbd764deaa72874dfd9eaa61dfda5277c70af9eceae313b7
audit.opt.out     4cd1e72050faa749cbd764deaa72874dfd9eaa61dfda5277c70af9eceae313b7
witt_hinge.py     66f2e84d82e93fdbf88483984d45d5d731b4152986678151b711c6274703eef5
witt_hinge.out    d55d931fb66b38b62b371aee651ea15702d8b5537e2cb3a23fa87b163147eaac
witt_hinge.opt.out
                   d55d931fb66b38b62b371aee651ea15702d8b5537e2cb3a23fa87b163147eaac
```
