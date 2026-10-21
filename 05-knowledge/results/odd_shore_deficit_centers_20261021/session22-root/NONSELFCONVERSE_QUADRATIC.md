# Quadratic order bound for non-self-converse presentation atoms

Status: PROVED, independently reviewed. This is a scoped extension of the cyclic-shore argument, not a classification of all atoms attaining a bound.

Let T be strong and non-self-converse, with h=H/|Aut(T)|<=d, d a positive integer. Set r to the least power of two at least max(4,2d). Then

n<=d(r-1)<4d^2 for d>=2.

For d=1, nonexistence follows directly from 2<=q<=h<=1 in the proof below; no external rigid classification is needed.

Proof. Choose a Hamilton cycle, successor R, automorphic-rotation index q, subgroup order m=n/q odd. No cycle cut is blue, so q<=h<=d (unique cut-equivalence map is a rotation). If q=1, T is an odd circulant; negation of cyclic coordinates is an anti-automorphism, contradicting non-self-converseness. Thus q>=2.

Suppose m>=r+1. A rotation shore A is regular and strong and has a directed r-cycle S by Moon pancyclicity. It is a proper subset. The other intact strong shores form a chain under the marked-cycle perfect matchings; any remaining vertex of A closes it, and all remaining A vertices link in and out. Hence T-S is strong.

For any other shore B, each A vertex beats the same number k of B vertices. If every B vertex met S in exactly r/2 arcs in each direction, summing would give rk=(r/2)m, impossible for odd m. Thus some v has at least r/2+1 exits from S, or at least r/2+1 entries into S. In the exit case take the corresponding r-cycle cuts followed by a common Hamilton word of T-S starting at v. Distinct such words differ by a nonidentity power of an r-cycle, an even-order permutation; it cannot be a tournament automorphism. There are at least r/2+1>=d+1 distinct path orbits, contradicting h<=d. For the entry case apply the same construction to the opposite tournament and reverse the resulting words. Opposite preserves h and non-self-converseness; it requires no ambient anti-automorphism. Therefore odd m<=r-1 and n=qm<=d(r-1).

For d>=2 the power-of-two minimality gives r<4d, and d(r-1)<4d². Notice this orbit-distinctness argument needs no half-supported reversal condition: unlike the self-converse deficit argument it never needs to distinguish converse mates, because every orbit counts toward h.

Combination with the separately audited self-converse quadratic black-deficit bound supplies a quadratic order cutoff for every bounded-h strong atom. This improves the finite search scale underlying the fixed-presentation rational series, but does not evaluate their unknown atom polynomials. The h3 classification remains separately inherited and sharper.

Source/preserved/lost: marked Hamilton cycle -> regular rotation shores and literal matchings; retains enough adjacency for distinct Hamilton words and an order obstruction. Forgets isomorphism duplicates, full automorphism structure, blue data in general, and metagraph ports. A cyclic presentation is not a class census.
