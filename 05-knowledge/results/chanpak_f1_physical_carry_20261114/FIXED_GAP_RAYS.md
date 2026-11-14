# Fixed-gap theorem for interfaces with at least two apex contacts

Let H be the ordinary Hamiltonian count, v its apex, and U the
vertices beaten by v, with r=|U|≥2. In the [finite-support defect](ORDINAL_DEFECT.md),
R=R₊=h+P is the exposed response and d=R−H=E+S. Here E counts
complementary nonempty path pairs whose left end and right start
are both in U, and S counts full puncture paths ending in U.

**Theorem.** At each fixed positive d, the H attainable with at
least two contacts form a finite union of isolated values and
dyadic rays a+b·2^s (b>0, s above a specified integer minimum).
For odd d only the isolated values occur. Thus at most O_d(log X)
different H occur through X. This includes every r≥2 at once;
it is not a result about one-contact or unrestricted scalar counts.

## A core bound for one marked pair

Choose {u,w}⊆U. Temporarily ignore all other contact labels.
Delete u,w to make W. In the ordered list of strong components
of W peel the maximal prefix beating both u,w and maximal suffix
losing to both. They are ordinal wings of the full puncture, so the
remaining C_uw contains both marks and is an interval of full
puncture components. We claim |C_uw|≤2+E_uw, where E_uw is the
part of E using exactly the two chosen endpoints.

Write E_uw=A₀+O₀+J by the [AND/OR contraction formula](CONTACT_PAIR_DEFECT.md),
and h_W=H(W), with H(empty)=1.

1. The total size of all nonsingleton components of W is at most
   h_W: each component's order is at most its path count by
   Hamiltonian cycle rotations, and the sum of the nonsingleton
   counts is at most their product.
2. Each surviving singleton beating both contacts has an earlier
   vertex not beating both: it supplies a distinct marked triangle
   in AND. Each surviving singleton losing to both has a later
   vertex beating at least one: it supplies a distinct marked
   triangle in OR. In each case the counted singleton is the unique
   endpoint of its required pure type. The [odd-cycle collection formula](https://arxiv.org/src/2307.05569v1/anc/redeiberge-long.pdf#page=14)
   preserves the contributions of W and adds at least two per
   marked triangle. Thus at most (A₀+O₀−2h_W)/2 pure singletons
   survive.
3. Any two singleton components of the same mixed type contribute
   one distinct placement to J. End the left path at the later
   singleton; start the right path at the earlier one. Put all
   other components earlier than the later singleton on the
   left, and all later ones on the right. The two possible
   singleton contact choices are distinct. Hence if the mixed
   types occur k₁,k₂ times, then k₁+k₂≤J+2.

Summing gives at most (A₀+O₀)/2+J+2 unmarked vertices, which
is at most E_uw if A₀+O₀≥4. The remaining possibility is 2
(the counts are odd). Then h_W=1, so there are no nonsingleton
components or surviving pure singletons, and J+2=E_uw is still
a bound. This proves the claim including W empty.

## The overlap step: charge every contact only once

The C_uw are intervals in the **same** full puncture SCC order.
Whenever two marked pairs share a vertex their intervals
intersect. The graph of marked pairs joined by shared contacts is
connected, so their union C is an interval containing U. Its
complement consists of an ordinal prefix G and suffix D without
contacts. Consequently

    |C| ≤ r + Σ_{u<w}|C_uw − {u,w}| ≤ r + Σ_{u<w} E_uw = r + E ≤ r+d.

The pair budget r(r−1)≤d bounds r in terms of d. Equivalently,
remove all contacts at once, then peel the initial SCCs beating all
contacts and the final SCCs losing to all: a universal wing belongs
to every pair wing, and the intersection of the pair prefixes
(respectively suffixes) is itself universal.

## Exact absolute decode, with unbounded tail length

If the suffix D is empty, the literal common-apex join gives

    d=m(E_C+S_C),     R=m(h_C+P_C),       m=H(G).

If D is nonempty, all of it beats the apex. The [no-contact reader law](ORDINAL_DEFECT.md)
then gives

    d=m h_D E_C,      R=m(h_D h_C+f_D P_C).

Here E_C is positive and even, since every unordered contact
pair has positive even contribution. Thus odd d requires
an empty suffix. In either equation, each required wing count
divides d. Discard all singleton components in G. In D move
the singleton components to a transitive tail of total length s;
this preserves h_D and f_D since they multiply across ordinal
sums. The remaining nonsingleton wings have total order at most
their respective path counts, so range over a finite menu. For the
suffix's nonsingleton part B, h_D=h_B and f_D=f_B·2^s.
Require s≥1 for empty B, and s≥0 otherwise.

Enumerate the bounded cores and wings, verify the appropriate
equation for d, and compare the given H+d with the readout. In the
suffix branch, divide H+d−m h_B h_C by m f_B P_C and test for a
power of two with the required minimum width. Every accepted
choice is backed by the actual ordinal join. This proves exact
membership and logarithmic counting, not just residue coverage.
A finite menu of fixed gaps cannot cover the scalar odd targets
using this interface, even allowing arbitrarily large tournaments
and every contact set of size at least two.

## Gap four beyond the original five rays

For d=4 there can only be two contacts. Both wings must be
transitive (their odd counts divide 4), and the core has at most
six vertices. Enumeration with the first two vertices marked
gives the no-suffix values 5 and 15, plus the following (a,b)
for H=a+b·2^s, s≥1:

| a | b |
|---:|---|
| −3 | 9 |
| −1 | 9, 15 |
| 1 | 11, 12, 16 |
| 5 | 22, 25, 30, 31, 34 |
| 7 | 30, 33 |
| 9 | 33 |

The constant 15 also lies on the first ray. For example
(1,12) comes from core masks (2,12,9,1) and a nonempty
transitive suffix. Run **node replay_fixed_gap_rays.js** to
replay the entire gap-four menu and check wing peeling and core
bounds across 33,866 marked-pair cases and 16,712 additional
multi-contact marks (these are labelled counts).

The marked-pair estimate was found in a background pass. The
overlap argument, all-contact theorem and complete gap-four menu
are derived here.
