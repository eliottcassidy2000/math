# Final focused extension: all-r9-grid levels1 through6

**PROVED constructors and all-order fan classification. FINITE-EXACT source
search. OPEN:** levels s>=7 in this packet, a uniform s-family, and Chan–Pak.
This note is separate from the preceding checkpoint. It extends the lower
three levels, and records the useful simpler sources that the incoming
nonunit-offset strategy revealed.

## 1. Five actual offset sources

For every r>=1 and1<=s<=6, the number

    N(9,r,s)=1+2^s+9*2^(r+s)

has an explicit tournament constructor. Level1 is the earlier ordinal
product3(1+3*2^(r+1)). For levels2 through6 use these selected bases:

    s2: B=(2,12,9,1),               O={2}, D=5,  U=9;
    s3: B=(2,28,25,17,1),           O={0}, D=9,  U=9;
    s4: B=(2,4,57,3,11,27),         O={0}, D=17, U=18;
    s5: B=(2,12,25,49,35,7),        O={0}, D=33, U=36;
    s6: B=(6,12,24,33,107,7,47),    O={6}, D=65, U=288.

Every selected set is a singleton, so W=0 and D=H(B). Let U=9*2^j.
The corresponding j values are0,0,1,2,5. Attach the ordinary selected tail
of length

    t=r+s−j.

It is positive for every r>=1. The exact inherited tail identity gives

    D+2^t U=1+2^s+9*2^(r+s),

proving every r, not just finitely many outputs. The tail lengths are
r+2,r+3,r+3,r+3,r+1 respectively. All vertices and arcs are emitted by
`realize9_six_levels.py`.

The simpler s2 source means no finite head is needed if this offset
architecture is used. The earlier kernel-plus77 construction remains valid
and is one vertex smaller for every r, so it is retained as a distinct
compact construction. The s3 singleton-mark source is likewise a simpler
alternative to the earlier ordinal product. No claim of new numerical
spectrum values or unsearched historical priority is made.

## 2. What the complete small search establishes

`search_nine_offsets.py` generated all56 unlabelled six-vertex classes and
every one-vertex extension of them, then tested every selected subset. Its
3640 presentations comprise56 six-vertex bases and3584 seven-vertex extensions.
Every seven-vertex tournament is represented because deleting one vertex
leaves a member of the complete six-vertex classification. Smaller sources
can be padded by unselected source vertices without changing D or U, so
this also detects any source of order below six.

The search finds sources exactly for s=2,3,4,5,6 among the tested targets
D=1+2^s and oddpart(U)=9. In particular no such s7 source has order<=7.
This is an exact bounded statement, not a nonexistence result at larger order.
A further focused search tests all128 one-vertex extensions of the displayed
s6 source and all selected subsets; none has D129 and oddpart(U)=9. Its
receipt is `offset_extension_receipt.json`. It does not cover every
8-vertex tournament or exclude another carrier.

The positive result is five infinite all-r rows supplied by actual sources,
not a collection of sample N-values. The census is only their provenance and
a boundary on this particular search.

## 3. Exact reason the tempting Fermat-fan continuation stops

Let F_n=C_3[TT_n,1,1], n>=1. Write its chain in increasing order, followed
cyclically by y,z, so chain→y→z→chain. Then H(F_n)=2^n+1: paths using y→z
are indexed by assigning each chain vertex before y or after z, and there
is one additional path z→whole chain→y.

For one selected chain vertex x, let p be the number of preceding chain
vertices and q the number following it, so n=p+q+1. **For all p,q>=0,**

    U_x=2^(p+1)3^q+2^(q+1)−0^p,                    (1)

where0^0=1. The two exterior singleton responses are

    U_y=2^n+1,       U_z=2^(n+1)−1.                 (2)

Here W=0 automatically for every singleton mark.

Proof of (1): count covers(P,Q) with Q starting at x, allowing P empty.
When y,z occur together in order y→z, their hosting path can be P or Q.
Each choice gives2^p3^q assignments of the other chain vertices to its
before/after runs and to the other path. If z is in P and y in Q, the
prefix is forced into P and the suffix has2^q choices. If both y,z are in P
in the order z→chain→y, there are2^q−0^p choices: the intervening chain must
be nonempty. The remaining split y in P,z in Q is impossible, since Q starts
in the chain and cannot reach z without y. These disjoint cases total (1).
For (2), Q starting at y either uses y alone, or y→z followed by an arbitrary
chain subset, giving1+2^n. Q starting at z has2^n choices if y belongs to P
and2^n−1 if y belongs to Q, giving the other formula.

**All-order classification.** A singleton mark in F_n has oddpart(U)=9
exactly in the following three cases:

    n2: first chain vertex, U9;
    n3: exterior vertex y, U9;
    n4: last chain vertex, U18.

For a chain mark with p=0, (1) is odd and increases strictly with q, reaching9
only at q=1. If p>q, its odd part is2^(p−q)3^q+1, which equals9 only at
q=0,p=3. If1<=p<q, its odd part is3^q+2^(q−p)>9. If p=q>=1, U has odd part
coprime to3 because3^p+1 is1 modulo3. The exterior formulas in (2) give only
n=3. These exhaust every singleton position and every n.

In particular, cloning vertex0 of the s4 source as a growing chain and
marking its first new clone gives

    H_m=8*2^m+1,
    U_first(m)=16*3^(m−1)+2^m,       m>=1,

not9*2^m. Its first values are18,52,152. Marking the last clone gives
U_last=8*2^m+2. Thus the proposed uniform s>=4 continuation fails structurally,
not just at the first hostile. The s5 and s6 sources above genuinely change
the carrier; they are not hidden reparametrizations of this fan.

## 4. Precise remaining offset/head gate

For a level s, it suffices to construct any selected tournament with

    D=1+2^s,     U=9*2^j,

and j<=s. More generally, j>s covers r>=j−s+1 and leaves exactly the finite
heads1<=r<=j−s to be supplied by separate constructors. A uniform family of
such sources, or another direct count-level construction for all s>=7, would
extend this packet's result. Neither the bounded census nor the Fermat-fan
classification gives that missing family.

## Reproduction and status

    python3 realize9_six_levels.py
    python3 realize9_six_levels.py --r 20 --s 6

`realize9_six_receipt.json` contains the five exact source profiles,18 literal
graph controls for r=1,2,3 and s=1,...,6, and16 controls of the all-position
fan formula. The all-r and all-position statements are proved above;
these controls are independent checks rather than extrapolation.
