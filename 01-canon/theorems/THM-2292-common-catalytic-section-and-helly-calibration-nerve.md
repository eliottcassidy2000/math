---
id: THM-2292
title: "Common catalytic section and the Helly calibration nerve"
status: >
  PROVED. For every finite labelled packet in a commutative nonexpansive
  integer-valued metric monoid, one context simultaneously attains the
  catalytic length of every subset sum. Hence the complete Boolean defect
  table of the catalytic length is one honest diagonal continuation slice,
  not a coordinatewise splice; the persistent continuation floor is bounded
  above by this realized catalytic defect. In a finite-dimensional
  seminormed space, the zero-defect complex is exactly the nerve of the
  atoms' exposed calibration faces. Helly's theorem therefore bounds every
  minimal nonface by effective dimension plus one. Conversely, for every
  m>=2, m atoms in l_infinity^m have zero defect on every proper subpacket
  and defect one on the full packet. This boundary-simplex interaction is
  translation-invariant, context-rigid, catalytic, and stable, so no
  fixed-arity defect shadow, or tournament constructed only from such a
  shadow, classifies stable interaction in general. Applied to the rational
  stable knot group this gives a finite-rank Helly gate, but no higher-order
  knot example or new unknotting number.
source: codex-2026-07-25-catalytic-calibration-nerve
depends_on:
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2259-boolean-continuation-hasse-field-and-signed-interaction-dividends
  - THM-2272-persistent-interaction-packing-and-calibrated-catalytic-defect-spectrum
  - THM-2281-common-optimal-context-for-finite-catalytic-families
related:
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
---

# THM-2292 -- catalytic tables are sections; stable tables are nerves

There are two different compatibility questions behind a packet of
composable objects:

```text
context compatibility:
  can one continuation attain all requested catalytic minima?

calibrator compatibility:
  can one additive dual functional certify all requested stable equalities?
```

The first has finite Helly number one because optimal-context loci are
additive upper ideals. The second is a convex-intersection problem whose
effective dimension controls the first possible missing face. Confusing the
two is precisely how a separately minimized scalar table can be mistaken
for a compositional state.

## 1. One context realizes the whole catalytic Boolean table

Let `(M,+,0)` be a commutative monoid with an integer-valued metric `d`
satisfying joint nonexpansivity

```text
d(a+c,b+e)<=d(a,b)+d(c,e).                           (1)
```

Put

```text
ell(x)=d(x,0),
rho_z(x)=d(x+z,z),
ell_cat(x)=min_z rho_z(x).                           (2)
```

Fix a labelled packet `x_1,...,x_r`, and for `S subset [r]` write

```text
x_S=sum_(i in S)x_i.                                 (3)
```

Apply THM-2281 to the finite family of all `2^r` subset sums `x_S`.
There is one context `z_*` such that, simultaneously for every subset,

```text
rho_(z_*)(x_S)=ell_cat(x_S).                         (4)
```

Define the catalytic Boolean defect

```text
Beta_cat(S)
 =sum_(i in S)ell_cat(x_i)-ell_cat(x_S).             (5)
```

At the single context in (4), THM-2259's diagonal defect is

```text
Delta_(z_*)(S)
 =sum_(i in S)rho_(z_*)(x_i)-rho_(z_*)(x_S)
 =Beta_cat(S)                         for every S.   (6)
```

Thus the entire table `(Beta_cat(S))_(S subset [r])` is one actual stalk
of the continuation field. It is not assembled from incompatible
face-dependent contexts.

Several consequences are immediate.

1. `Beta_cat` is normalized, nonnegative, monotone, and superadditive on
   disjoint subsets, because the same is true of the one context table
   `Delta_(z_*)`.
2. Its zero sets form an abstract simplicial complex.
3. Its Hasse increments and signed Möbius dividends are exactly the Hasse
   field and dividends of one context:

   ```text
   c_cat(i|A)
    =ell_cat(x_i)+ell_cat(x_A)-ell_cat(x_(A union {i}))
    =c_(z_*)(i|A),                                  (7)

   h_cat(A)
    =sum_(B subset A)(-1)^(|A|-|B|)Beta_cat(B)
    =h_(z_*)(A).                                    (8)
   ```

4. If `underlineDelta` is THM-2272's persistent continuation floor, then

   ```text
   underlineDelta(S)
    =min_z Delta_z(S)
    <=Beta_cat(S).                                  (9)
   ```

Equation (9) is only an upper bound. The context which minimizes all
individual and merged diagonal responses need not minimize their
difference. No equality between persistent and catalytic defect is claimed.

For knots, (4) says that every finite packet has one knot `J_*` satisfying

```text
d_G(K_S#J_*,J_*)=u_cat(K_S)             for every S. (10)
```

Hence all catalytic subset defects, Hasse increments, and dividends of the
packet occur simultaneously at one genuine Gordian context.

## 2. Stable additivity is a calibration-face nerve

Let `(V,p)` be a finite-dimensional real seminormed vector space. Quotient
by `ker(p)` when necessary. Fix labelled vectors `x_1,...,x_r`, let

```text
Y=span_R{x_1,...,x_r},             d=dim(Y/ker(p|Y)), (11)
```

and let `B_*` be the unit ball of the dual of the resulting normed quotient.
For each atom define its exposed calibration face

```text
F_i={phi in B_*:phi(x_i)=p(x_i)}.                   (12)
```

Hahn--Banach makes every `F_i` nonempty; if `p(x_i)=0`, the condition is
automatic on the quotient.

For `S subset [r]`, put

```text
Delta_p(S)=sum_(i in S)p(x_i)-p(x_S).                (13)
```

> **Calibration-nerve theorem.**
>
> ```text
> Delta_p(S)=0
>   iff
> intersection_(i in S)F_i is nonempty.             (14)
> ```

If `phi` lies in the intersection, then

```text
p(x_S)>=phi(x_S)=sum_(i in S)p(x_i)>=p(x_S),        (15)
```

so equality holds. Conversely, suppose `Delta_p(S)=0`. Choose
`phi in B_*` calibrating `x_S`; it may be chosen of norm one when
`p(x_S)>0`, while `phi=0` is allowed when the calibrated sum vanishes. Then

```text
sum_i p(x_i)
 =p(x_S)
 =sum_i phi(x_i)
 <=sum_i p(x_i).                                    (16)
```

Every summand in the final inequality is termwise bounded, so equality of
the sums forces `phi(x_i)=p(x_i)` for every `i in S`. This proves (14).

Consequently the stable zero-defect complex is literally the nerve of the
finite family of compact convex faces `(F_i)`.

## 3. The finite-rank Helly gate

The faces in (12) live in the `d`-dimensional real vector space dual to the
effective span in (11). Helly's theorem and (14) give:

> If `Delta_p(A)>0`, then there is a subset
>
> ```text
> B subset A,              |B|<=d+1,                (17)
> ```
>
> with `Delta_p(B)>0`.

Equivalently, every minimal nonface of the stable zero-defect complex has
arity at most `d+1`. Thus subsets of arity at most `k` determine the
zero-versus-positive defect support whenever the effective stable span has
dimension at most `k-1`; no claim about recovering larger defect magnitudes
is made.

This is a rank-sensitive completeness statement, not a claim that pair
data are normally complete. Even in low finite dimension the nerve need not
be flag, and in unbounded dimension there is no uniform arity ceiling.

For the Gordian metric, THM-2191 puts the homogenized length `u_hash` on
the rational Grothendieck group and identifies its dual unit ball with
additive Gordian-1-Lipschitz calibrators. Restricting to the finite rational
span of a knot packet and extending scalars to the corresponding real
finite-dimensional quotient gives (14)--(17) for

```text
Delta_hash(S)
 =sum_(i in S)u_hash(K_i)-u_hash(K_S).               (18)
```

If every atom is individually stable-tight,
`u_hash(K_i)=u(K_i)`, then (18) is exactly THM-2272's common-calibrator
stable interaction rate on every subset.

## 4. Boundary-simplex defects occur in honest normed spaces

The Helly bound does not make binary data sufficient. Fix `m>=2`, take

```text
V=R^m,                   p(x)=||x||_infinity,

mathbf1=(1,...,1),       x_i=mathbf1-e_i,
                                      1<=i<=m.       (19)
```

Every atom has norm one. If `S` is proper and `s=|S|`, then

```text
sum_(i in S)x_i
 =s*mathbf1-sum_(i in S)e_i.                        (20)
```

The coordinates indexed by `S` equal `s-1`, while at least one coordinate
outside `S` equals `s`. Therefore

```text
p(x_S)=s=sum_(i in S)p(x_i)             if S proper. (21)
```

For the full packet,

```text
x_[m]=(m-1)mathbf1,

p(x_[m])=m-1,
Delta_p([m])=1.                                    (22)
```

Thus its zero-defect complex is exactly the boundary of the
`(m-1)`-simplex. Every proper subpacket is additive, while the full packet
has a one-unit saving.

The dual unit ball is the `l_1` ball. The calibration faces make the
mechanism visible:

```text
F_i=conv{e_j^*:j!=i},

intersection_(i in S)F_i
 =conv{e_j^*:j notin S}                 for S proper,

intersection_(i=1)^m F_i=empty,                     (23)
```

where the middle identity is for nonempty proper `S` (the intersection of
the empty family is the whole dual unit ball).

Give `V` its translation-invariant metric

```text
d(x,y)=p(x-y).                                      (24)
```

Then for every context `z`,

```text
rho_z(x)=p(x),

ell_cat(x)=ell_hash(x)=p(x).                        (25)
```

The missing full face in (22) is therefore simultaneously:

```text
a root interaction;
a persistent all-context interaction;
a catalytic interaction;
a homogenized stable interaction.                  (26)
```

There is a literal same-shadow comparator. Keep the labelled packet size
`m` but put

```text
y_i=e_1                         for every i.         (27)
```

Then every atom again has norm one and

```text
p(y_S)=|S|                      for every S,         (28)
```

so every defect of the `y` packet is zero. The `x` and `y` packets
therefore have identical defect tables on every proper subset, but their
full defects are respectively one and zero.

Choosing `m>k` proves that no fixed `k`-ary truncation determines stable
interaction in general. In particular, for `m>=3` the two packets have the
same empty pair-defect graph. No tournament or binary relation constructed
only from that pair shadow can recover which packet has the unique
`m`-body defect.

This strengthens the abstract stopping boundary in THM-2248 from an
integer word metric to a genuine normed and already homogenized object. It
does not realize the boundary simplex by knots.

## 5. Information ledger and scope

The three finite carriers now separate cleanly:

```text
catalytic contexts:
  additive upper ideals;
  every finite family has one common section;

stable calibrators:
  exposed convex faces;
  their nerve is exact and obeys a finite-rank Helly gate;

pair/tournament shadow:
  a useful low-arity certificate;
  incomplete as soon as the calibration nerve has a higher missing face.
                                                               (29)
```

The source-to-target maps are:

```text
all subset sums x_S
  -> one common optimal context z_*
  -> one realized catalytic Boolean table;

stable atoms x_i
  -> calibration faces F_i
  -> their intersection nerve
  -> zero stable-defect complex.                    (30)
```

The first map loses the size and identity of a smallest context. The second
loses the full off-diagonal continuation kernel and retains only stable
additivity. Neither operation computes a new unknotting number, produces a
positive knot catalyst, or exhibits a higher-order Gordian minimal
nonface. QED.
