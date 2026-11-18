# Random-source synthesis: certify a continuation by its finite reader jet

**All-order elementary corollaries of physical substitution:** a sharp visit bound, finite exact interpolation certificate, bounded crossing count, and exact search interval. No universal Chan–Pak selector or historical novelty claim for the substitution kernel. This changes the useful question from “do these sources agree at this count?” to “can their entire continuation be certified by a finite collection of readings?”

## Random provenance

A reproducible draw used Python Random(938257) on top-level Markdown files under `05-knowledge/results`, excluding filenames containing tournament, chanpak, paley, metagraph or hamilton. The six draws were:

1. `continuing1_20260906_jets_audit.md`;
2. `mixed_triangle_refinement_20261004.md`;
3. `collatz_circulant_20260930_circulants_lucas_cubic_monotile.md`;
4. `collatz_golden_carriers_20261004.md`;
5. `planar_jc48_sep08_two_cusps_audit.md`;
6. `boundary_support_20261020_seven.md`.

Filename exclusion is not a claim that the contents have no tournament connections. The useful bridge was the first note's **complete jet rather than endpoint**: enough exact observations certify a whole source response. Mixed refinement independently warns that agreement of counts does not retain parent incidence. The circulant note warns that fixed-operator spectra do not justify changing-source dynamics. Those cautions are retained below.

Prior direct tournament theory: THM-3121's path-cover/walk-content substitution kernel and THM-3134's endpoint-jet/Newton transform already provide the general compositional framework. The explicit single-vertex readout and bounds here are a specialization and synthesis, not a claim to supersede them.

## Literal physical reader

Let T have n>=1 vertices and distinguished vertex v. Replace v by a transitive tournament TT_m of m>=1 vertices, all with the same external contacts v had. Write h_T,v(m) for its Hamiltonian count.

For r>=1 define c_r as follows. Order the n-1 outside vertices, choose r of the n slots before/between/after them, insert the symbol v in those slots, and count the words whose consecutive symbols are arcs of T. This is a finite nonnegative integer defined from the actual source orientation. It includes the complete marked contact information needed here, not just H(T).

Then, for R=floor((n+1)/2),

    h_T,v(m) = sum_(r=1)^R c_r r! S(m,r),       c_1=H(T),

where S is the Stirling number of the second kind.

**Proof.** Compress each maximal run through the inflated block to v. This gives exactly one of the counted words. Distributing the m distinct transitive-block vertices over its r ordered nonempty visits is a surjection onto r labels: each label's assigned vertices must be traversed in their increasing order. Conversely every such surjection yields a unique Hamiltonian path. There are r!S(m,r) of them. Between two visits there cannot be just one outsider u, since this would require both v->u and u->v. Consequently n-1>=2(r-1), giving r<=R. The one-visit words are precisely the Hamiltonian paths of T.

The visit bound is sharp in order: take a directed word v,u_1,w_1,v,...,u_(r-1),w_(r-1),v and orient all listed arcs accordingly, completing the remaining pairs arbitrarily. On2r-1 source vertices this gives c_r>0.

## Complete finite certificate

The values h(1),...,h(R) recover every c_r recursively:

    c_r = [h(r)-sum_(j<r)c_j j! S(r,j)]/r!.

The quotients must be nonnegative integers for a physical source. In particular, if two rooted sources with maximum visit bound R agree at these R inflation sizes, they agree for **every** m>=1. This certifies a repeatable scalar reader replacement along this particular transitive-inflation operation. It does not certify equality of endpoint counts, other port readings, or arbitrary external contexts.

More strongly, r!S(m,r)=sum_(j=1)^r (-1)^(r-j)binom(r,j) j^m. Any difference of two such reader sequences is therefore an exponential polynomial with bases1,...,R. Unless it is identically zero, it has at most R-1 distinct real zeros on m>=1. Proof: divide by its smallest exponential and differentiate; the constant vanishes and Rolle's theorem reduces the number of terms by one. Induct on the number of nonzero terms. Thus agreement at any R distinct positive integer sizes certifies all sizes. A claimed stable repair with too many exact returns must therefore be a full identity in this fixed-source model, not an endlessly intermittent coincidence.

## Neutrality and target-height consequences

c_2>0 exactly when v lies in a directed triangle. Necessity follows from the first outsider of an inter-visit path being beaten by v and its last beating v: the first transition along that path from the former side to the latter produces a triangle through v. For sufficiency begin with v_1,u,w,v_2 in the two-twin inflation of a triangle v->u->w->v, and insert each other vertex into this directed path. Such insertion is always possible in a tournament: insert before the first beaten path vertex after a prefix beating the newcomer, with the endpoint cases included. This never removes the outsiders separating the twins. Compressing gives a two-visit word.

If c_2=0, the same necessity argument excludes every multiple-visit word, so h(m)=H(T) for all m. Otherwise

    h(m) >= H(T)+2^m-2,
    h(m+1)-h(m) >= 2^m.

The second bound follows from c_2>=1 and the nondecrease of every other surjection count as m grows (append the largest new vertex to any fixed occupied visit to inject assignments). Hence fixed-source active inflation cannot give an unlimited fine-step repair window. A necessary condition for target N>=H(T) is

    m <= floor(log_2(N-H(T)+2)).

The reader is strictly increasing, so binary search within this exact interval decides whether this source hits N. Formula evaluation uses only the finite c-profile and integer powers. This is a fully physical, bounded target test for one source, not a construction of a suitable source for every target.

For a directed triangle with v marked, c=(3,1), so h(m)=2^m+1. This simple example already shows why complete source knowledge is not coverage: the exact decoder has a sparse image.

## Next useful experiment and its limit

Store sources by their complete c-vector rather than by H alone. A proposed source swap preserving an infinite inflation continuation can be certified by a finite reader jet, and unequal profiles cannot agree at more than R-1 heights. To close Chan–Pak one must still find physical profiles whose ranges jointly cover the targets, or change contacts/reader so this fixed-profile obstruction no longer applies. Neither arbitrary nonnegative c-vectors nor scalar equality are declared physically realizable.

## Replay

`reader_interpolation_replay.py` enumerates every labelled rooted source through order4, computes c by literal slot words, and independently evaluates inflated tournaments for m=1,...,5 by subset dynamic programming. It checks the sharp visit bound, triangle criterion and exact inverse interpolation. Normal/optimized JSON agree. These finite checks support the code; the all-order claims have the combinatorial proofs above.


## A literal infinite +2 lift

Use the five-vertex outgoing-mask rows [24,1,3,6,14] and
[18,8,3,5,14], marking vertex 1 (zero-indexed) in both.
Counting the valid slot words gives respectively (c1,c2,c3)=(11,6,0)
and (13,6,0). The third coordinate can be checked by its only possible
slot pattern v,a,b,v,c,d,v; the replay enumerates all patterns
independently. Consequently the common marked inflation yields

    (6*2^m-1, 6*2^m+1), for every m>=1.

Here both actual counts grow while their difference stays exactly 2.
Both sources are strong; replacing a marked vertex of a strong
tournament by a nonempty module preserves strong connectivity:
paths to or from the module use its uniform contacts, and different
module vertices can travel via an outside directed cycle through the
original marked vertex. The replay verifies the source connectivity
and direct counts through m=6. This family has only logarithmically
many centers below a bound; it solves the repeatability question for
these two sources and this reader, not universal placement.
