# Atom-dependent four-corner positive witness (preserved before task switch)

**Status:** exact exploratory witness in two certified THM-2625 fields.
The calculation was run inline and still needs extraction into a standalone
normal/optimized companion before canonization.

The aligned objects are

```text
A = source one-sided carrier,
B = pullback(target one-sided carrier),
C = A intersect B.
```

Neither exclusive wing `A\B` nor `B\A` is used as the common atom.

Selected physical cell and ancestry copy:

```text
(s,t,e) = (0,4,1),
common weighted piece index = 1,
interval = [142004992589460,142005019034340),
multiplicity = 27581135604.
```

One of the positive multiplicity copies is retained with unit weight.  It
has the exact delayed semantic coefficient

```text
source: carry 12, root one, value 103478815440;
target: carry 6,  root one, value 103478815440.
```

The atom index was retained while constructing all `169*13` source and
target dual carrier cells.  Only after that atomwise construction were the
endpoint transforms taken.  In the first certified field the dual supports
are `522/2197` on both sides and the `k=l=0` endpoint supports are `169/169`.

Use

```text
L=R=(0,0),
source step s=(0,1),
target step t=(12,0),
det(s,t)=1.
```

The four target differences are

```text
(0,0), (0,1), (1,0), (1,1),
```

and the determinant corners are `(0,0,0,1)`, so their mixed face is `1`.

For prime `352341050142921841`, the four endpoint amplitudes are

```text
(345283638862540358,
 180096700909697779,
 167108689435423878,
 144423867160388279).
```

The four products are

```text
(347197813355232446,
 3438130309897238,
 320064096166488418,
 146419401152414063),
```

and their Hadamard coordinates are

```text
(112437340698188483,
 236493496489149044,
 165063327916487722,
 170114988031260853).
```

Every displayed field element is nonzero.

For prime `956354278959359281`, the corresponding endpoint amplitudes are

```text
(77115584218272234,
 620226426191367084,
 950936471790018101,
 445878683356068526).
```

The products are

```text
(788257626617883404,
 251268809443191686,
 161379668316310808,
 381695617330433150),
```

and the Hadamard coordinates are

```text
(626247442748459767,
 496451150414331132,
 316672868160569376,
 757304766188814060).
```

Again every entry is nonzero, and in both fields

```text
D_0 D_3 - D_1 D_2 = 0.
```

This is strong evidence for a literal same-atom coefficient square and a
nonzero THM-2620 determinant mixed face.  Remaining audits before promotion:

1. extract the inline calculation into an assertion-free pinned companion;
2. replay ordinary/optimized/stored output;
3. verify that choosing one unit copy from the positive multiplicity retains
   all hidden ancestry/semantic-owner labels rather than only their aggregate;
4. state precisely which exact `(r,k,l,L,R)` lift represents the displayed
   `k=l=0` section; and
5. retain the first-failure label `common_after_all_factors` and the
   representative gauge exhaustively.
