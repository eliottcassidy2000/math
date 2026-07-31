# Horn–collar q0 hinge, minimal `V_4` globalization, and Witt endpoint obstruction

**Status:** proof draft plus `FINITE-EXACT` primary computation, scratch only.
No row is excluded and LRC(14) remains open.  The exact hinge imports the
promoted THM-2825, THM-2847, and THM-2851 universes by pinned hashes.  The
universal groupoid and affine-exponent arguments are proved below.

## 1. Verdict

Three objects that previously looked merely analogous meet on one exact
twenty-cell locus.

1. Every cell in THM-2847's `E3`-only horn is a nonempty THM-2825
   common/right cell.  In each of them the selected q0 atom

   ```text
   I=(142004992589460,142005019034340)
   ```

   is exactly the even `+2h` collar image

   ```text
   R=I-2h  ->  M1=I-h  ->  M2=I.
   ```

   There are twenty labelled hinges but only one physical interval triple.

2. THM-2825's `M_3 tensor I_587` is exactly the compression of the regular
   `V_4` pair groupoid obtained by deleting its missing `(1,0)` object.
   Twisting signs or cocycles cannot restore a missing matrix-unit support.
   The universal minimal globalization adds one new 587-dimensional object
   and gives `M_4 tensor I_587`.

3. THM-2851's sharp joint target/carry carrier has 169 states, exactly the
   cardinality of the endpoint-address plane.  But the operations differ:
   the inherited endpoint plane is split `F_13^2`, while carry fidelity needs
   nonsplit `C_169`, equivalently the additive group of `W_2(F_13)`.
   No affine rechart of `F_13^2` has order 169.

Thus cardinality, grade, and cell labels are all already paid.  Moreover, a
nonlocal common-path move physically realizes the first carry character:
on twelve long horn cells, step 2 to step 68 shifts both endpoint masks by
`T^8`, where `8` generates `C_13`, while retaining factors, carriers,
semantics, and literal ancestry.  The missing operation is now narrower: a
typed **full Witt lift** coupling that carry reference to q and to the
`E3`/complement macro blocks, together with the new cofiber object required
by the `V_4` globalization.

The stopping boundary is exact.  Only the q0 allocation interval belongs to
the collar; all q1 through q12 intervals, in particular q3/q7/q11, are
outside.  The hinge root has empty source endpoint support while `M2=I` has
81 source endpoint addresses.  No permutation can identify those masks.

## 2. The exact twenty-cell hinge

THM-2847's `E3`-only horn is

```text
clock=1,
sigma in {0,3,8,9,12},
target in {5,6,9,10}.                                    (1)
```

The primary audit reconstructs the THM-2825 common/right bank on all twenty
cells in `(1)`.  Their aggregate census is

```text
20 cells,
52 labelled right roots,
4,076 labelled common pieces.                            (2)
```

Put `h=401080680`.  In every cell, direct lookup in the complete labelled
interval bank gives

```text
R  =(142004190428100,142004216872980)=I-2h,
M1 =(142004591508780,142004617953660)=I-h,
M2 =(142004992589460,142005019034340)=I.                 (3)
```

The arrows in `(3)` are literal THM-2825 half-step arrows.  Their semantic
pattern is uniformly

```text
(live,dead,live).                                         (4)
```

The rooted common-path lengths containing `M1,M2` are

```text
sigma 0:  144  (8 cells after including sigma 3),
sigma 3:  144,
sigma 8:   14,
sigma 9:   40,
sigma 12: 118,                                           (5)
```

more compactly

```text
144:8, 14:4, 40:4, 118:4.                                (6)
```

The twenty labels forget to one physical triple in `(3)`.  This is an exact
hinge, but not twenty independent physical copies.

In native-factor order

```text
(E3,clock,q1,q2,c2,c3),                                   (7)
```

the four source/target factor masks are uniformly

```text
R:  011111 / 111111 / 111111 / 011111,
M1: 111111 / 111111 / 111111 / 111111,
M2: 111111 / 111111 / 111111 / 111111.                   (8)
```

Across the thirteen carrier twists they are

```text
R:  source empty / target delta0,
M1: source delta0 / target delta0,
M2: source delta0 / target delta0.                        (9)
```

The selected endpoint masks inherited from independently audited THM-2818
have cardinalities

```text
R:  source 0,  target 81,
M2: source 81, target 81.                                (10)
```

The target masks in `(10)` agree by the identity translation.  The source
masks cannot agree under **any** permutation because their cardinalities
differ.

Let `I_q=I+q(T/13)` be the thirteen allocation intervals in the pullback
coordinate.  The complete membership table on every cell of `(1)` is

```text
q=0:       I_q=M2,
q=1,...,12:I_q is outside both common and right support.  (11)
```

Therefore the q3/q7/q11 horn legs have empty physical interval fibre product
with the collar.  The q0 endpoint row is a genuine shared hinge, but its
`QA/QAB` semantic columns are zero.  The semantic columns appear only after
moving to q11/q7, precisely where `(11)` leaves the collar.

### 2a. A physical `T^8` first-carry reference on the long paths

The failure of a normalized one-step endpoint translate is not the end of the
carry lane.  On the twelve cells

```text
clock=1,
sigma in {0,3,12},
target in {5,6,9,10},                                   (11a)
```

the distinguished path has length at least 118.  Let `K` be its step-68
common piece.  Since `M2` is step 2,

```text
K=M2+66h.                                                (11b)
```

Direct endpoint evaluation gives, uniquely on both sides,

```text
Endpoint_source(K)=T_(0,8) Endpoint_source(M2),
Endpoint_target(K)=T_(0,8) Endpoint_target(M2).          (11c)
```

Both masks have 81 points.  The move `(11b)` is even, so it preserves the
nonzero semantic value.  Both intervals are common pieces in every cell of
`(11a)`, hence retain all six native/pulled factors and source/target carrier
`delta0`.

The literal contributor computation gives at both endpoints

```text
|U|=966606, |V|=28534,
digest=15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd,
```

with exact equality of the two `U` sets and of the two `V` sets.  The
distinguished sheet remains supplied.  Thus `(11b)--(11c)` retain the complete
literal ancestry, not only its cardinality or digest.

Write the q0 endpoint mask as

```text
H=1_A tensor 1_B,
A={0,1,2,3,4,5,6,7,12},
B={0,1,2,3,4,5,8,9,10}.                                 (11d)
```

It is a `9 x 9` rectangle.  Its augmentation is
`81=3 mod13`, so THM-2839's p-group lemma makes `H` a unit in
`F_13[C_13^2]`.  The physical signed difference from `(11c)` is

```text
D_8=(Y^8-1)H.                                            (11e)
```

It has eighteen `+1`, eighteen `-1`, and 133 zero coefficients.  Since

```text
Y^8-1=(Y-1)(1+Y+...+Y^7)
```

and the second factor has augmentation `8`, it is a group-algebra unit.
Consequently

```text
rank(D_8 on F_13[C_13^2])=13*12=156.                    (11f)
```

Projecting along the first endpoint coordinate gives rank `12`, nonzero on
every nontrivial carry character.  Because `8` is invertible modulo thirteen,
`T^8` generates the same first carry quotient as `T`; demanding the literal
coordinate `T` would be an unnecessary normalization.

This is a genuine positive sidecar: a two-vertex, oriented, fully decorated
first-carry current inside the `E3` block.  It is not a `C_169` action, does
not reach q3/q7/q11, and supplies no transporter to the complementary `E3`
block.  Those are now the precise remaining debts.

## 3. The punctured `V_4` object is a compressed pair groupoid

Let

```text
G=V_4=F_2^2,
U=G minus {m},
H_U=direct_sum_(u in U) R_u,        dim R_u=r.            (12)
```

Grade the matrix unit `E_(v,u):R_u->R_v` by `v-u in G`.
The resulting algebra is

```text
End(H_U)=M_3 tensor I_r.                                  (13)
```

For each nonzero `g in G`, the homogeneous piece `A_g` contains exactly two
directed matrix units: the two orientations of the unique `g`-edge contained
in `U`.  This is exactly THM-2825's three collar degrees.

The partiality is visible in every product of nonzero degrees.

- If `g=h`, then `A_g A_g` contains two of the three diagonal units in
  `A_0`.
- If `g!=h`, then `A_g A_h` contains one of the two units in `A_(g+h)`.

Thus each of the nine ordered products `(g,h)` of nonzero degrees has a
one-object support defect.  Multiplying existing matrix units by signs or
nonzero cocycle values preserves support and cannot fill it.

There is no hidden scalar-twist obstruction: the pair groupoid on three
objects is contractible, so every normalized scalar two-cocycle is a
coboundary.  Equivalently, choosing a base object explicitly gauges every
arrow weight to one.  The defect is missing support, not phase.

Nor can a three-object set be a transitive `V_4`-set: orbit sizes divide
four and are `1,2,4`.  Therefore any global `V_4` action extending `(13)`
must add at least one object.

Adding

```text
X=R_m,             dim X=r,                              (14)
```

gives the regular four-object pair groupoid

```text
End(H_U direct_sum X)=M_4 tensor I_r.                    (15)
```

It is universal and dimension-minimal.  To see the rank invoice directly,
write the first collar maps as isomorphisms

```text
d:R->M1,       a:M1->M2,       S=ad:R->M2.              (16)
```

For any chosen identification `X=R`, the Koszul filler

```text
0 -> R --D0--> M1 direct_sum X --D1--> M2 -> 0,

D0(v)=(dv,v),
D1(m,x)=a(m)-S(x)                                        (17)
```

is exact and contractible.  Indeed `D1 D0=ad-S=0`; if `D1(m,x)=0`, then
`a(m)=ad(x)` and the isomorphism `a` gives `m=dx`; and `D1` is surjective
through its `M1` coordinate.

Conversely, any mapping-cone filler cancelling `S` must carry rank at least

```text
rank(S)=r.                                                (18)
```

For THM-2825, `r=587`.  Hence one 587-dimensional missing object, or seven
new matrix units per root, is sharp.

This is still coefficient algebra.  Physically, every right root has empty
source carrier support and at least one native-factor hole, whereas `M1,M2`
have source `delta0` and all six factors.  No sign twist changes a support
zero into a one.  The semantic population mismatch `573/14` independently
prevents filling `(14)` by permuting existing right roots.

## 4. The endpoint plane has the right size and the wrong extension class

THM-2851 proves that the natural residue/carry lift is the nonsplit extension

```text
0 -> C_13 -> C_169 -> C_13 -> 0,                         (19)
```

whose lifted generator `L_1` satisfies

```text
L_1^13=T,       order(L_1)=169.                          (20)
```

A faithful joint target/carry realization therefore needs at least 169
states, and equality is attained by the `C_169` torsor.

The endpoint address bank also has 169 states, indexed by `F_13^2`.  Its
inherited translations are split and have exponent 13.  Even allowing every
affine rechart does not change this:

> **Affine exponent lemma.**  For every odd prime `p`, every `p`-element of
> `AGL(2,F_p)` has order dividing `p`.  In particular
> `AGL(2,F_13)` contains no element of order `169`.

Embed an affine map in homogeneous coordinates:

```text
x |-> A x+b
      corresponds to
[[A,b],[0,1]] in GL(3,F_p).                              (21)
```

A `p`-element is unipotent, so it is `I+N` with `N^3=0`.  In characteristic
`p>=3`,

```text
(I+N)^p=I+N^p=I.                                         (22)
```

This proves the lemma.  The boundary is sharp: at `p=2`, `AGL(2,F_2)=S_4`
contains elements of order four.

The exact `p=13` enumeration finds

```text
169 linear p-primary candidates,
28,561 affine p-primary candidates,
orders: 1 once, 13 on the other 28,560,
order 169: none.                                         (23)
```

The same 169-element set can be equipped abstractly with the cyclic law:

```text
(W_2(F_13),+) is isomorphic to Z/169Z.                   (24)
```

But `(24)` is a nonlinear change of operation, not a relabelling inside the
inherited affine endpoint category.  It must be constructed and physically
typed.

THM-2852's Cayley-tournament nonsingularity does not alter this conclusion.
An invertible convolution on the split group algebra pays spectral rank, but
it commutes with the split translation action and cannot increase its
exponent.  Full tournament spectrum and carry holonomy are orthogonal
invariants.

## 5. Why the formal grade cancellation is real but insufficient

The THM-2851 horn attachment changes the abstract
`(outer-word,E3-macro-truth)` grade by

```text
(1,1).                                                    (25)
```

THM-2825's even collar changes `(semantic,source-carrier)` by

```text
(0,1).                                                    (26)
```

If their second coordinates could be physically identified, diagonal
composition would have degree

```text
(1,1)+(0,1)=(1,0),                                       (27)
```

exactly the missing corner of the punctured `V_4`.

Equation `(27)` is not numerology: the two constructions meet at the exact
q0 atom `M2=I` on all twenty horn cells.  But `(10)--(11)`, `(19)--(23)`,
and the absent q11-to-q7 physical action show the information destroyed by
the grade quotient:

```text
cell label,
grade,
cardinality
    do not retain
physical interval, source endpoint, or extension law.    (28)
```

The algebraic composition exists only after forgetting precisely the three
coordinates that the LRC realization needs.

## 6. Connection contract and next decisive test

```text
source:
  THM-2825 rooted collar, THM-2847 twenty-cell E3 horn,
  THM-2851 q3/q11/q7 natural carry triangle;

target:
  a four-object V4 globalization carrying a faithful C169 endpoint law;

maps:
  R --(+2h)--> M2=I at q0;
  formal horn attachment of degree (1,1);
  proposed Witt carry on the 169 endpoint states;

preserved now:
  all twenty cell labels, the literal q0 atom, semantic parity,
  target endpoint mask, target carrier delta0, and grade arithmetic;

destroyed now:
  q3/q7/q11 interval support, source endpoint support, native E3 at R,
  the missing fourth cofiber object, and nonsplit carry holonomy;

needed sidecar:
  a physical object X of rank 587; extend the proved T^8 first-carry
  reference to a lawful W_2(F13) q/carry action and transport it between
  E3 and its complement;

cheapest decisive tests:
  (a) prove every lawful endpoint rechart is affine, which makes the
      no-order-169 obstruction terminal for this bank; or
  (b) extend the physical D_8 reference from the q0/E3 block across the
      q11-to-q7 E3/complement attachment, and test the resulting Witt
      triangle on (sigma,target,clock)=(0,5,1).
```

The current result does not choose between `(a)` and `(b)`.  It replaces a
vague request for “another 169-state sidecar” with the exact missing
operation and the sharp minimal cofiber invoice.
