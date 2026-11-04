# Interior ternary repairs and source-preserving barriers

**PROVED scoped word identities; FINITE-EXACT complete attachment menus.**
The first balanced source of the preceding session has H=59. Two interior
insertions give a balanced H=179 source, and a further interior insertion
gives a balanced H=541 source. However no extension retaining the original
H=59 tournament as an induced subtournament has H=181. This is an
obstruction to source-preserving induction, not a forbidden Hamilton value.
The full prime and 7p realization problem remains OPEN.

## Inheritance and the operation being tested

The actual source and its physical names come from
[the previous moving-mark construction](../../orbit_spectrum_session_20261103/coverage/NOTE.md),
not an arbitrary primitive integer pair. The closest mechanism is the
two-birth reset NR, whose output rooted profile is (b,t,t,b). The canonical
hostile is the previous H15 all-marking suffix obstruction. The corrected
near miss is the
[all-root reset barrier, READY section 4](../../two_entry_funnels_20261103/READY.md):
re-marking any standard reset of a strong source of order at least four
cannot give a one-ready or two-ready root. The least-used data here are
the interior insertion cut and the complete incidence of a new vertex
with the protected root and its named child.

The live board is actual balanced source / interior cut / arbitrary
attachment mask / root-child incidence / source-changing return. Its
first test is whether both target values can be attained without changing
any old edge. The positive step and the complete obstruction below show
that these two branches behave differently on the same actual source.

The numerical family 45*2^m-1 already occurs as H(G(T_m)) in
[Session 164 section 4](../../session164_emergent_20261103/nonuniform_operator/NOTE.md).
The added result is a specified insertion program from the preceding
balanced source, its retained interface, and its exact primitive guard.
No new numerical spectrum family or finite coverage record is claimed.

## An infinite interior plus-two step with retained ports

Use the inherited two-apex convention: the apices have edge x->y; a chain
letter a has x->v iff a&1 and y->v iff a&2; chain vertices are in ascending
order. The preceding balanced C_m has an actual presentation

    W_(m,1) = 202 0^(m-1) 3 2,                  m>=2.

Its root is the first chain vertex, and its child is the second. Both are
old physical vertices. For every r>=1 put

    W_(m,r) = 202 0^(m-1) 3^r 2.

This is the two-birth NR reset, at the second apex, of the strong base
word 2 0^(m-1) 3^r 2. The root and child therefore have the balanced
profile (b,t,t,b), where

    b = 5(2^r+1)2^m-1,
    t = (2^r+1)(2^m+2)-1.                       (I1)

One exact proof of the two counts uses the inherited 36-state transfer.
For each prefix 2 and 202, the state after the first 0 satisfies

    (M_0^2-3M_0+2I)v = 0.

For the states with one and two zeroes, after the first 3 the residual
(M_3^2-3M_3+2I)v has zero output after terminal 2 and after each of the
next 35 powers of M_3. Cayley-Hamilton makes these zero outputs valid at
every further power. Thus both parameter directions have scalar
recurrence f_(k+2)=3f_(k+1)-2f_k. The four values m=2,3 and r=1,2
determine (I1). All state identities and all 144 scalar residuals are
listed in the certificate. This is an exact finite matrix proof, not
extrapolation from a numerical grid. The additional grid is only a check.

The restriction r>=1 matters: W_(2,0) has H=31, whereas substitution
of r=0 into (I1) would give 39.

Insert two new 3 letters immediately before the old 3. This retains the
old apices, every old arc, the named root and child, and their order. It
takes r=1 to r=3, hence

    (b,t) = (15*2^m-1,3*2^m+5)
       -> (45*2^m-1,9*2^m+17) = (3b+2,3t+2).    (I2)

The output pair is primitive exactly when m is not 13 modulo 14:

    gcd(45*2^m-1,9*2^m+17) = gcd(9*2^m+17,43),
    2^m = 22 (mod43) iff m = 13 (mod14).

Here 2^7=-1 mod43 and ord_43(2)=14. Thus the old mod13 input guard and
the new mod43 output guard are different actual attachment conditions.
Primality is tested only for explicitly supplied values; no prime
infinitude follows from (I2).

The plus-two counts 45*2^m-1 are never divisible by 7: their residues for
m=0,1,2 modulo 3 are 2,5,4. The competing plus-four targets 45*2^m+1
are divisible by 7 exactly when m=1 modulo 3. The source m=4 therefore
tests an actual 7p obligation: 3*239+4=721=7*103.

This particular repair does not iterate within one unsplit 3-run. To get
3b+2 at fixed m by increasing r to s, (I1) requires

    2^s = 3*2^r+2.

For r>=2 the right side has 2-adic valuation one and is greater than two,
so it is not a power of two. The only r>=1 solution is r=1,s=3. A split
run, another attachment relation, or a genuine source change is needed.

## A literal two-step path and an exit from the two-apex class

At m=2 the following two insertions are legal and retain both old ports:

    202032  -> 20203332 -> 202030332,
       59       179          541,
    (59,17)  (179,53)    (541,147).

The first step inserts 3 letters at zero-based positions 4 and 5 of the
new word. The second inserts 0 at position 5. Their counts are respectively
3H+2 and 3H+4. All three H values are prime. Full graph rows, fresh names,
old-to-new maps, and insertion cuts appear in check.json.

There are exactly two one-vertex ternary attachments to the displayed
H179 source. In the word's vertex order, with the two apices first, the
in-neighbour masks of the new vertex are:

| Target | In-neighbour mask | Retained balanced pair | Two-apex presentation |
|---:|---:|---|---|
|539=3*179+2|638|(539,157)|none|
|541=3*179+4|124|(541,147)|202030332|

The first is an actual exit from the two-apex word class while preserving
the named balanced root and child. The count 539=49*11 is outside the
prime and 7p lanes with p prime. The example separates a genuine
operation restriction from a requirement of the prime-lane induction.

For the resulting H541 source, all one- and two-vertex attachments miss
both 1625 and the prime 1627. There are 2048 first masks; exactly 36 have
H at most 1627, and all 147456 second masks extending those parents were
checked. This is only a two-birth exclusion: the strong-order bound still
allows three added vertices. The three-birth arbitrary extension problem
for this H541 source remains open in this packet.

For the original C_4 source at H239, the complete one/two/three-insertion menu over all
three initial two-apex markings has 33948 addresses. It has 30 plus-two
addresses and no plus-four address. The inherited strong-order bound at
721 permits at most 13 vertices, so three new vertices exhaust the
strong two-apex extension budget.

There is no hidden gain from choosing a fresh vertex as a final apex.
Every C_m contains the two vertex-disjoint directed triangles
(apex1,root,child) and (apex0,old3,last2), and therefore requires at least
two deletions to become transitive. If an extension had a two-apex
presentation with a fresh apex, deleting its other apex would already
make the old source transitive, a contradiction. Thus both final apices
must be old, and the complete old-pair superword menu applies. Trimming
leading 0 and trailing 3 padding outside the old strong source leaves a
strong superword. Consequently the m=4 source has no H721 extension
whose feedback vertex number is at most two, at any order. This does
not exclude extensions with feedback vertex number at least three.

## Complete attachment obstruction at the actual H59 source

Let T be the actual eight-vertex word 202032. For any new vertex z, let
I be its old in-neighbours and O the complementary set. Let e_i,s_j be
the Hamilton path end and start counts of T, and let W_ij count ordered
two nonempty spanning paths, first ending at i and second starting at j.
Deleting z from a Hamilton path gives the exact identity

    H(T+z) = sum_(i in I)e_i + sum_(j in O)s_j
             + sum_(i in I,j in O)W_ij.          (I3)

The three terms correspond to z last, first, and internal. The ordered
path-cover matrix is computed from independent subset end/start tables.
All 256 first masks are evaluated. Only 24 have H<=181; the other 232
cannot lead to either target because insertion never decreases H. Every
one of the 24*512=12288 second masks is evaluated. No one-birth target
occurs. Exactly six ordered two-birth pairs give 179:

    (63,63), (63,319), (63,383),
    (127,63), (127,127), (127,383).

All six retain the old balanced root with pair (179,53). None gives 181.
The unpruned ordered two-birth universe is 2^17=131072. No root choice,
transitive-complement assumption, strongness filter, or attachment-shape
restriction is imposed on these masks.

This yields an all-order obstruction for this source. In any tournament
containing strong T, the SCC containing T has H at least 59, and its H
divides the total H. If the total is 181, primality forces that SCC itself
to have H=181. The inherited strong-order floor gives at most ten vertices
at this H, so it contains at most two new vertices beyond T. The complete
two-birth certificate excludes it. Thus no induced extension of this
particular source has H=181, despite the inherited positive realization
of 181 by a different source in READY section 3.

## How much source change is necessary

All 28 single arc reversals of T change H. Among the 378 unordered pairs
of distinct arc reversals, exactly five retain H=59. Four retain the old
root-to-child arc and are strong; every one of these four still misses
181 under all one/two births, hence under all induced extensions by the
same SCC argument.

The fifth reverses arcs (2,3) and (5,7), destroys the protected root-child
edge, and is weak. Its two-birth exclusion must not be promoted to the
same all-order statement: its strong core has fewer old vertices. The
certificate records this boundary explicitly.

Consequently an endpoint H-preserving arc repair that retains the
protected root-to-child edge and enables an induced 181 extension must
change at least three old arcs. This is a lower bound on final edge
changes, not on the number of native moves or the H values of intermediate
states. The coordinator's constructive return supplies the matching
three-arc positive example.

## Reproduction and information retained

Run from the repository root:

```sh
python3 -B 05-knowledge/results/balanced_ternary_returns_20261103/interior/check.py
python3 -B -O 05-knowledge/results/balanced_ternary_returns_20261103/interior/check.py
```

Both outputs must equal check.json. The checker uses explicit exceptions,
and 56 small arbitrary births are independently compared with literal
subset path DP. The inherited transfer and prior source receipt are
hashed. The scalar 179 does not identify which old 3 lies among the three
new 3-run positions; the six named attachments retain that distinction.
The map to a balanced pair preserves the stated reset queries but loses
the whole next attachment menu. The H59 obstruction is the decisive
witness that the latter cannot be supplied by primitivity alone.
