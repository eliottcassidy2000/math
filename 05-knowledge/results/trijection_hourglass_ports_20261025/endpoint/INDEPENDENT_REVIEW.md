# Independent structural review of the cyclic-prefix and tail repairs

Reviewer: separate trijection lane, 2026-10-25. Reviewed `THEOREM.md` and
`TAIL.md`; no producer source was edited. **PASS within their stated
domains.** The mathematical review below is distinct from the producer's
random or exhaustive controls and from the root lane's independent
Gaussian-elimination audit. No global connectivity conclusion follows.

## Cyclic-prefix theorem

The legality proof keeps the actual port h. Its new joins are exactly
`a->h`, `e->b`, `(h-1)->y`; the reflected guards are respectively the
incoming port, `3->c`, `(n-h)->x`. The XOR cancellation restores the
same labelled target. There is no hidden module or clone assumption.

The stronger cyclic-cut threshold n>=11 is valid here. On an even-label
orbit the only possible incoming neighbours of a are c and h; y and z
are odd labels. Transitivity on that orbit forces its induced regular
indegree to be at most two. Odd n instead makes the entire tournament
regular and forces n=9. This uses the actual pair y,c, not merely the
number of endpoint ports.

The Q2 image formula has exactly one exceptional pattern: `w=b` when
`h=(n-2)/2`. In all other cases w is neither b nor an incoming port.
At n12,h5 I checked the explicit positional reflection: the equations
`S(1,9)=S(0,8)` and `S(2,10)=S(11,7)` have exactly the stated pairs.
Their nonpath complements, P-blue relation and z-port exclusions give
the displayed contradiction. No unstated c-versus-x direction is used.

For Q3, the two moving images lie consecutively in the cyclic middle
order. The only middle ports are1,h, and they are adjacent in that order
only for h2 or h=e; those two cases have separate nonport witnesses.
The positional reflection product fixes a,z, so preserving a's exact
four remaining outports is necessary. This establishes the claimed
uniform blackness.

## Tail theorem

The new x->4 ledger is legal for every allowed middle completion, whether
the middle is open or closed. The new Q2 joins reverse `a->(n-5)` and
`e->a`. The omitted x->3 and3->4 joins are restored after two complements.
All labels and path-edge XOR terms agree with the stated target.

Independent reflection witnesses, written in the opposite orientation
to the producer's equivalent witnesses, are:

| Move | Compared arc bits | Values |
|---|---|---|
| Q1 at T | `T(a,c)` and `T(x,3)` | 0,1 |
| Q2 at S | `S(x,3)` and `S(b,z)` | 0,1 |

The final reflection product was derived independently as the composition
of +1 rotations on the first and last n-4 positions. Its pointwise formula
agrees with `TAIL.md`. For odd n it is transitive. For even n it has two
n/2-cycles, and the cycle containing a meets a's five outports only in b.
Thus the induced regular tournament would have outdegree one, forcing
an orbit of size three. This proves generic blackness at every n>=12.

At n11 the transitive cycle and the specified a-row determine **at most
one** labelled B, not merely one isomorphism class. Its a-row is compatible
with an orientation because its five cycle differences select exactly one
from each opposite pair. Direct reconstruction agrees with the published
source rows after applying K_P. Therefore every potential failure of this
particular generic final word is the single stated source.

The complementary 4->x branch needs no new source assumption. The old
ell1 ledger is legal under the tail ports; the sharper cyclic-cut argument
applies to its first reversal. The intermediate degree comparison is
5 against4. The final reflection product consists of two disjoint
four-cycles because n>=11. This proves that branch at every stated order.

## Independent order-eleven replay

I decoded the11 published rows directly into a Boolean adjacency matrix,
without importing either endpoint producer. The source is a tournament,
P is directed and blue, and x->4 holds. I independently reconstructed B
as the circulant whose cyclic order is
`(0,1,2,3,4,6,7,8,9,10,5)` and whose a-outneighbours are1,7,8,9,10.
It equals K_P of the stated rows. The generic final word is actually blue.

The alternative certificate has three directed black words. Independent
literal anti-reversal mismatches at their actual sources are:

| Word | Pair | Reflected reversed pair | Bits |
|---|---|---|---|
| first | `(0,6)` | `(3,10)` | 1,0 |
| second | `(0,2)` | `(10,8)` | 0,1 |
| third | `(0,3)` | `(10,8)` | 1,0 |

Complementing every nonpath unordered pair after each word yields exactly
K_P(T), with all11 labels unchanged. Thus the unique real failure of the
generic ledger is repaired, not discarded by a source filter.

## Consequence and limits

The whole tail family `I={n-4,n-3,n-2}` is proved repaired for n>=11,
including both middle states. The arbitrary-port theorem is a separate
sufficient repair. Neither proof asserts that every remaining h has one
of the requisite chords. Failure of the resulting bank still means
uncovered by the bank, not a black bridge.

The improved cyclic-cut lemma also permits a strengthening of the sibling
overlapping-prefix theorem to n>=max(11,2ell+3), provided its two new
small-overlap cases are explicitly handled: n13,L8 has a transitive
reflection product, and n14,L8 has a's orbit
`(0,3,6,9,12,2,5)`, containing at most the two outports c,h. This is an
audited integration: the sibling theorem's statement and normal/optimized
receipts were subsequently updated to include orders11..14.

## Final all-layer generalization audit

I subsequently reviewed `GENERAL_TAIL.md`. **PASS:** for every k>=4 and
n>=2k+3, the consecutive-tail source has a literal three-black detour.
This is an all-k structural theorem; it is not inferred from the finite
checks, and it does not cover arbitrary nonconsecutive port sets.

For the generic x->4 ledger, independently composing its final reversal
with r_P gives +1 rotations on `[0,e-1]` and `[k,n-1]`, with e=n-k.
At odd n their product is transitive; at even n a's orbit is exactly
`(0,1,...,k,k+2,...,e-2)` and has n/2 vertices. Its intersection with
`N^+_B(a)={b} union [e,n-1]` is just b. The regularity contradictions
are therefore independent of k, except for the one odd candidate
n=2k+3. Both initial blackness witnesses remain literal for arbitrary k.

In the complementary 4->x ledger, the cyclic-cut obstruction remains
strict: for even n the incoming a-neighbours on the even-label orbit
number exactly floor(k/2), forcing n=4floor(k/2)+2; for odd n regularity
forces n=2k+1. Both lie below n>=2k+3. The intermediate degrees k+1,k
and the two final four-cycles are also valid uniformly in k.

At the boundary n=2k+3 I independently checked the claimed cyclic order
and difference set. They give exactly one candidate B, and r_P acts in
those coordinates as t->2k+1-t. The nonzero differences are paired into
D and -D, with no overlap. All consecutive P arcs have difference1 or
k+2, so the candidate is an actual source of the stated type after K_P,
not a hypothetical symmetry assignment.

The exceptional repair's directedness is exact. Its first new source
arcs k->w,w->u,v->s have B differences2,k,2 outside D. Its second new
joins complement first(Z)->u and v->z, whose B differences are k-2,k,
also outside D. The last word restores k->u,v->w,w->s after two flips.
An explicit XOR certificate is

    E(C1) xor E(P) = E(C3) xor E(C2)
                  = {ku,vw,ws,kw,wu,vs}.

The three stated anti-reversal mismatches are correct at their actual
sources. In particular the last one uses B(z,s) of cyclic difference
k+5; this lies in D exactly at the required boundary k>=4. It is not
an unchecked extrapolation from the old eleven-vertex certificate.

As a separate finite control, I reconstructed B directly as Boolean
adjacency matrices and replayed the three exceptional words at
k=4..12,40,80,128, without importing the general-tail producer. Every
source has directed blue P and x->4; every intermediate word is directed
and black; every final matrix equals B. These controls supplement the
symbolic proof and do not set its unlimited quantifier.

No correction to `GENERAL_TAIL.md` was required. The documents and this
review are frozen after this audit.
