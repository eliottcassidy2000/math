# Exact square-factor carry across adjacent high products

This is a physical all-height carry in the two-block growing-spacer
seam of [RADIX4.md](RADIX4.md). It is **not** a dense low bank or a
solution of prescribed-target selection.

Fix r≥4, D=2^r, left width k=r+1, right width r, and spacer t=r−1; let K=2^t=D/2.
Choose any a=3+8j with j≥0 and a+1<D/2, and b=D−a. Define m=D+a.
The unshifted masks (a,b) have high product P=m²; shifting both up by
one gives factors (m−1,m+1) and product P−1. Shifting both down by one
gives the same product P−1 with its factors exchanged. All three sets of
masks are admissible. Define the middle remainder C=C(a,b).

**Exact result.** Counts of the downshifted, upshifted, and unshifted states,
respectively, are
```
(P−1)K+C,
(P−1)K+C+D·3/8+3j−[j odd],
PK+C.
```
Consequently the signed upward carry gap is D/8−3j+[j odd]. At fixed r
the possible signed gaps comprise two finite step-six progressions as j
ranges over its even and odd admissible values. These are gap progressions,
not progressions of absolute counts with a common baseline.

To get a **literal +2** carry for every r≥4 choose
```
a=(D−7)/3  if r is even,
a=(D+1)/3  if r is odd.
```
These masks have a≡3 (mod 8), and bit 3 of a (indexed from zero) is 0 on even r and 1
on odd r. Thus the three outputs become
```
(P−1)K+C,     (P−1)K+C+K−2,      PK+C.
```
The low phase drops by K−2 precisely as the high product advances by one.
The two exhibited low digits are the bottom and top of a radix-K phase,
not an assertion that the intervening low digits occur.

## Exact proof of the low reset

The upshift changes the low bits of a from 011 to 100 and those of b from
101 to 110. The reversal difference changes by D/8. For the internal pair
energy J, the change for a is −2−5(a−3)/8+[bit 3 of a], and for b it is
−(b−5)/4. To see this, in a the old pair (0,1) of energy 2 is removed.
Each higher set bit i≥3 has net change −5·2^(i−3), plus the new adjacency
at i=3 if present. In b the sole internal low pair keeps energy 2, and
each higher set bit has net change −2^(i−2). Thus
```
C(a+1,b+1)−C = (3D+3a−9)/8−[bit 3 of a]
               = 3D/8+3j−[j odd].
```
For the downshift, both masks a−1,b−1 are even; adding bit zero to each
increases the left-minus-right reversal by D/2. The added pair energy is
(a+b−2)/2+1=D/2 (one adjacency, since bit 1 of a−1 is set and that
of b−1 is not). These cancel, proving C(a−1,b−1)=C.

## When this is a canonical carry

For the explicit +2 choices,
```
C = (6r+1−D)/9     if r is even,
C = (6r+7−2D)/9    if r is odd.
```
For a proof, increasing r by two adds just a bit at position r to a and at
r+1 to b, with no adjacency to an old set bit. Each reversal R therefore
goes to 4R+2 and each J goes to J+R. Thus X=R_b−2R_a goes to 4X−2, while
C goes to C+2X. The given expressions follow by induction from r=4,5
and their respectively observed C=1,−3 and X=−2,−10.

At r=4 the canonical residues are 7→1 with K=8 (counts 2887→2889).
For all r≥5, −K+2<C<0, so the two adjacent totals actually share the
same quotient in the *unshifted* canonical radix. The coordinate P is
nevertheless physical and genuinely changes. Shifting the phase by the
explicit value W=C−1 gives digits K−1→1 and increments the normalized
quotient; equivalently the two low states then have digits 1,K−1.
A full physical low bank would have to supply the missing digits at that
same phase and baseline. A finite exact-divisor scan for r=4,…,110 finds
no adjacent third point at either of these two high products; that scan
does not exclude other factor-crossing families or multi-root repairs.

### Unbounded canonical radix along a subsequence

Although the full seam radix K is not a canonical rollover for r≥5,
there are infinitely many unshifted canonical carries at every *smaller*
dyadic radix. For r=2s≥6,
`v₂(C−1)=2+v₂(3s−2)`; for r=2s+1≥5,
`v₂(C−1)=2+v₂(3s+1)`. In each case the omitted dyadic term in the
numerator of C−1 has strictly higher valuation, since even its magnitude
exceeds the positive linear term. Given q≥0, the congruence
3s≡2 (mod 2^q) in the even branch (or 3s≡−1 in the odd branch) has
infinitely many solutions. Taking r sufficiently large, the two
adjacent physical counts are then respectively 2^(q+2)−1 and 1 modulo
2^(q+2), and their quotients differ by one, without a phase shift.
This makes the canonical radix unbounded across the family, but its
maximum for any fixed r is O(r), not K=2^(r−1).

## A third consecutive point from touching the spacer

Keep the same outer masks, and use the *high* state d=e=m. Temporarily
allow a mark at position j of the spacer, with 0≤j<t; require it not to
be adjacent to either outer block, so no adjacency bonus crosses the
boundary. Its change in the readout before interactions with other new
marks is
```
d·2^j − e·2^(t−1−j).
```
This follows by summing its interaction with the old marks on each
side as in RADIX4.md, and subtracting from its bare signed linear
term. Thus any reflection-symmetric set U of interior spacer marks
changes the readout by exactly `−J(U)`: the single-mark terms cancel
in reflected pairs, while all pair terms between new marks are
subtracted exactly once.

For odd r≥7 the spacer t is even and at least six. Choose the two
mirror positions j=t/2−2 and t/2+1; they are interior and are three apart,
so J(U)=4. Since the upshifted low state is already H−2, we have three
consecutive physical odd counts
```
H−4 (two central contacts at separation 3), H−2 (high product P−1), H (high product P).
```
The two mirror positions nearest the center instead have separation 1
and J=2, giving a second physical realization of H−2. This is a
larger contact architecture: the extra point does not belong to the
original unmarked-spacer factor bank. It works uniformly along odd
r, but the missing long run below the carry is still open.

### A low-digit ceiling for these fixed outer factors

For any set of spacer marks, including the boundary positions, the
change from the unmarked high state is
```
m S(U) − J(U) − [0 in U] − [t−1 in U],
S(U)=Σ_(j in U)(2^j−2^(t−1−j)).
```
Both outer-adjacent marks are present here because a,b are odd. The
full energy of all t positions is (2^t-2): each upper position i≥1
contributes (2^i-1) from the geometric series and one adjacency.
Hence if S≠0, the absolute change is at least m−2^t. A small offset
must have S=0, which by uniqueness of binary expansions means U is
reflection-symmetric. In even t≥6 an interior two-point symmetric
set has energies 2,4,16,…; larger nonempty symmetric sets have at
least four marks and energy at least 12, and a symmetric pair at the
boundary has energy at least 16. Therefore H−6 is impossible by
spacer-only marking with these fixed outer masks. In odd t≥5, H−4 is
impossible in the same way (two-point symmetric distances are even,
with energies 2,8,…, and three marks already cost at least 6).
In particular filling the low digit requires changing an outer
factor or a genuinely different architecture, not merely adding more
central contacts near the balanced state.

## Whole-spacer complement coalesces into the high state

The bottom low state (both outer masks downshifted) has outer factors
(d=m+1,e=m-1), count H−K, and no old marks adjacent to the spacer.
Marking **the entire spacer** changes the count by
```
(d−e)(K−1)−J(full) = 2(K−1)−(K−2) = K.
```
Thus this low outer factorization meets the unmarked high factorization
at exactly the same total H, at every height. This is different from
having a +2 path between the totals.

More generally, with spacer marks U, put V=full\U and
B(V)=[0 in V]+[t−1 in V]. The energy contribution of a single index j
when the spacer is full equals
`2^j+2^(t−1−j)−[j=0]−[j=t−1]`, including all pairs incident to j.
Subtract the missing indices' incident sums and add back J(V) once
to get J(U). It follows that the resulting count is
```
H + m S(U) − J(V) − B(V).
```
Under complementation S(U)=−S(V). Reflection changes the sign of S and
keeps J and B, so the entire palette of spacer marks in this low
outer state is identical to that of spacer marks in the high outer
state, via complement followed by reflection. This equality preserves
the count, not the marks or the factor pair. It holds for all masks,
including asymmetric masks and boundary contacts.

The top low state (both outer masks upshifted), by contrast, has readout
`H−2 + m S(U) − (F(U)+R(U)) − J(U)`, with F and R the forward and
reverse binary values of U. For t≥3, S(U) can never equal ±1: if a binary
integer x and its reversal differ by exactly one, the set of changed
bits under increment, a beginning interval of bit indices, must also
be reflection-invariant. Thus every bit would change, forcing
x=2^(t−1)−1 and rev(x)=2^t−2=x+1, possible only at t=2.
So any S>0 is at least two, and the shift is then at least
2m−3K+4; negative S makes a large negative shift. When S=0 the mask
is symmetric, and any nonempty symmetric mask at t≥5 contributes at
least eight to F+R. Consequently for r≥6 no extra top-low spacer mask
changes its total by ±2, ±4 or ±6. In particular these three outer
states, even allowing all spacer masks, cannot provide the missing
fourth point H−6 at odd r≥7.

## A test for changing the reader beyond the square family

For arbitrary outer masks, write their factors as d,e,
their remainder as C, and α,β for whether their two boundary bits
are marked (left then right). For an interior mask U including
possible spacer boundary positions, put F(U)=Σ_(j∈U)2^j,
R(U)=Σ_(j∈U)2^(t−1−j), and let u₀,u_last indicate its endpoint
contacts. Directly from the head formula its readout is

    Kde+C+d F(U)−e R(U)−J(U)−αu₀−βu_last.

Set W=reflection(full spacer \ U). Then
F(U)=K−1−R(W), R(U)=K−1−F(W), and
J(U)=K−2−F(W)−R(W)+[0∈W]+[t−1∈W]+J(W).
The last equality deletes the total pair incidence at missing
sites from the full spacer energy and adds back their internal
pairs. Consequently the entire readout becomes the same formula
in W with
`(d',e')=(e+1,d−1)`, `(α',β')=(1−β,1−α)`, and
`C'=C+e−d+2−α−β`.
This acts on physical reader masks when the new factors are
admissible and their actual flags match; if the new remainder
differs, the discrepancy is an exact fixed count translation.
In particular it gives a concrete non-scalar compatibility test
for a next factor pair. This algebra does not guarantee that
every sought pair of factors or remainder has a physical reader.

For the widths k=r+1 used here, the admissible partner is given by
`a'=D−1−b, b'=D−1−a`. Whenever admissible, this particular whole-palette
match is the tournament converse: reverse the order of the full
spine including the sentinel and complement every head bit.
The old sentinel moves to the first site and the old first site
moves to the sentinel. So coalescence of the low-bottom and high
palettes is explained by a graph bijection, not an independent
new source of additional low digits. The **+2** between the low-top
and high states is not a converse identity and is the actual carry.
