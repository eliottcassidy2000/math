# Least-counterexample descent: decode defects zero, two, and four

**Proved in the canonical marked binary encoding:** every odd target
with defect 0, 2, or 4 occurs as a tournament Hamilton count except 7
and 21. The defect-two positive cases need a change of old arcs or
of source: no induced extension of the canonical source can read
its target. This is not a proof of Chan–Pak for larger defects.

## The code and the defect

Use the retained-root law in ANS_INVERSE.md. Start from 0→1 rooted at
0, with ordinary count h=1, end count e=0, two-cover marked start
count t=2 and deletion count d=1. Appending a sink (bit 0), or a sink
with only its arc to the root reversed (bit 1), leaves d=1. The
remaining state is

    0: (h,e,t) -> (h,0,2t).
    1: (h,e,t) -> (h+t−e,1,2t−1).

After k bits, t runs bijectively over [2^k+1, 2^(k+1)]:
t=2^(k+1) minus the ordinary integer of the chronological bits.
Thus every odd target N>=3 has an unambiguous physical source with
t=N, and its word ends in 1. Define the defect to be t−h. The
proof idea is to classify the location of the missing contribution
by its first carry rather than reason from the numerical target
alone. A zero is a carry that deposits t into the defect; a first
1 after a zero cancels exactly one, and consecutive 1s do not
change the defect.

## Exact classification at defect two

The first run 0^a 1 creates defect 2^(a+1)−2. Every later nonempty
zero run 0^b 1, starting at marked count u, adds u(2^b−1)−1 to
the defect. Here u>=3; the new contribution is at least 2, with
equality only if u=3 and b=1. The only way to be at u=3 is the
first letter 1. Subsequent trailing 1s preserve the defect.

Consequently the defect-zero odd targets come from all-one words and
have h=t directly. The defect-two odd targets are exactly these:

* words 01^r with r>=1: t=6·2^(r−1)+1, h=t−2;
* words 101^r with r>=1: t=10·2^(r−1)+1, h=t−2.

## Physical readers for the two families

**First family.** r=1 gives the forbidden 7. For r>=2, take the
five-vertex source of word 011, masks [6,28,24,17,1], with vertex 4
marked. The remainder of the word simply inflates that mark by
TT_m with m=r−1: later 1 vertices have the same outside contacts
and form a transitive block. Reverse the arc between unmarked
vertices 1 and 3 to get [6,20,24,19,1]. The marked visit profiles
are respectively (11,6,0) and (13,6,0), zeros thereafter.
Transitive inflation reads c_1+c_2(2^m−2), giving t−2 and t for
every m>=1. So one literal arc reversal on the canonical source
decodes this whole ray after its forbidden first value.

**Second family.** r=1 gives 11, witnessed by [24,1,3,6,14]; r=2
gives the forbidden 21. For r>=3, use the alternate seven-vertex
source [110,72,2,4,15,30,60], marking vertex 0. Its complete
visit profile is (41,20,0,0), zeros thereafter. Inflating the
mark by TT_(r−2) gives 20·2^(r−2)+1 = t. This is a replacement
of the source, not asserted to be a single arc reversal.

The marked coefficients count legal compressed path words; the
literal replay enumerates them. Substituting a transitive block
fills two visits in 2^m−2 ordered nonempty ways, and one visit
in one way, so the certificates apply at all widths.

## The next two defect layers

Defect four comes exactly from words 1101^r, r>=1. Indeed the
initial zero run would contribute 0, 2, or at least 6. After an
initial contribution of 2 any later zero run adds at least 6,
since t>=7. From initial contribution 0, the first later
zero run starts at u=2^s+1 after s ones; length b=1 contributes
2^s. A run of length at least two contributes at least 8, and
two later zero runs cost more than 4. Thus s=2,b=1.

The targets in this layer are t=18·2^(r−1)+1.
At r=1 use the count-19 witness [40,1,3,6,15,30].
For r>=2, use the seven-vertex source [30,104,2,4,14,29,61]
marked at 0, whose complete visit profile is (37,18,0,0).
Transitive inflation at length r−1 gives the target, all heights.
The source is inherited from RECURRING_READER_WINDOWS.md, where
its profile was recorded for a different purpose.

The same classification at defect six gives only 001^r, r>=1,
hence t=14·2^j+1, j=r−1>=0. An initial contribution of 2 cannot
continue below 8. Starting with no initial zeros, a single
unit-length later run costs a power of two, not 6; a longer run
costs at least 8. After even the cheapest first later run, t>=11,
so two later runs cannot stay within defect six. This leaves
exactly the initial two zeros.

For a *least* hypothetical counterexample in defect six, j
must be odd: even j forces 3 to divide t, and t/3 is neither
7 nor 21 (t is 1 mod 7), so the lower factor would be realized.
In fact j≡7 mod10 is also excluded, since then 11 divides t
(2^7=7 mod11). More generally any composite t in this ray is
excluded for a least counterexample: t is 1 mod7, so neither
proper factor in an odd factorization could be 7 or 21.
No assertion is made that all the prime targets here occur.

## Why an add-only descent cannot prove even this lane

For all these canonical sources h>=5. Induced addition is
nondecreasing in H by vertex insertion, and each increment is
even by Rédei oddness. If any induced supertournament had h+2,
along single-vertex additions exactly one step would add two,
starting at count at least five. UNIT_VERTEX_RETURN.md excludes
this: the only one-vertex +2 steps are 1→3 and 3→5. Thus all
defect-two sources fail as induced readers, even the many for
which the exact target exists by the constructions above.

The coincidence with 7 and 21 is exact in this encoding, but
does not imply that other defects introduce no holes. The
least hypothetical counterexample (or any counterexample)
outside the known exclusions must have even defect at least 6.
A possible induction is to classify all next nonzero defects and
find physical rewrites connecting them to already decoded types.

For each fixed defect D, the search for source *types* is finite:
at the last zero run of length b starting at u, we have
u(2^b−1)−1<=D and therefore u2^b<=2(D+1). After k bits
t>=2^k+1, so the last zero is within the first O(log(D+1))
positions, independent of the final address. The remaining
1s act by a transitive block, reducing exact continuation checks
to finitely many marked visit signatures for a given D.
This finite reduction does not supply those signatures or a
uniform bound over unbounded D.

Run python3 two_defect_replay.py.
