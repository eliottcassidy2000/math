# Independent odd-cycle proof of second-bit response faithfulness

Status: PROVED, using the cited Rédei modulo-four refinement plus a separate Boolean cancellation argument. This addendum preserves the earlier frozen exact and modular faithfulness files. It gives a second proof of the same second-bit XOR detector, not a new claim about arbitrary digraphs.

## Classical input: all odd cycles, not just triangles

For a tournament T, let c_odd(T) be the number of all simple directed cycles of odd length at least three. A cycle is counted up to cyclic rotation, with its directed orientation retained. It need not be spanning. The classical refinement used here is

H(T) ≡ 1+2 c_odd(T) (mod 4).

Hence the second response bit b(T)=((H(T)−1)/2) mod 2 is precisely c_odd(T) mod 2. A verified primary source is Grinberg–Stanley, “The Redei–Berge symmetric function of a directed graph,” arXiv:2307.05569v1, Theorem 7.1, with proof in Section 7: https://arxiv.org/html/2307.05569v1 . This addendum cites that input and does not claim to derive the entire modulo-four theorem independently.

## Setup

Use the fixed-endpoint experiment from `MODULAR_FAITHFULNESS.md`: Q has k≥2 labelled core vertices; u,v are distinct; the k−1 ordered base vertices are b_1,...,b_{k−1}; and

R=(v,b_1,w_1,b_2,...,w_{k−2},b_{k−1})

is the alternating path through every extension vertex except u. Its d=2k−3 underlying cross edges are varied independently. Fix b_{k−1}→u and freeze all other attachments. Let T(x) be the resulting tournament.

Then directly from odd-cycle enumeration,

XOR_{x∈{−1,+1}^d} c_odd(T(x)) mod 2 = Q_uv.

Together with the classical input, this independently proves the same identity with the Hamiltonian second bit b(T(x)). It does not use the Hamiltonian-word survivor argument.

## Boolean cancellation proof

Expand each c_odd(T(x)) as a sum of indicators over all possible labelled simple odd directed cycles, counted once up to rotation. Exchange this finite sum with the sum over orientation assignments, all in F_2.

If a candidate cycle omits an underlying edge of R, its indicator does not depend on the corresponding orientation bit. Pair assignments that differ only at that bit; its contributions cancel modulo two. Thus a surviving cycle must contain every underlying edge of R.

Such a simple cycle contains R as a contiguous path block: every internal vertex of R already has its two incident cycle edges used by R. The cycle must close this block between its endpoints v and b_{k−1}. Since R already visits all extension vertices except u, there are only two possibilities:

1. Close R directly by the edge between b_{k−1} and v, without u. This would have 2k−2 vertices and even length, so it is not among the cycles being counted. For k=2 this would be a prohibited two-cycle; it is excluded in any event.
2. Close R through u. This uses all 2k−1 vertices, so the resulting cycle has odd length. Because the extra edge is fixed b_{k−1}→u, its direction forces R forward, followed by b_{k−1}→u→v. This cycle exists in exactly one assignment of the R orientations if u→v is a core arc, and in none if v→u.

Therefore precisely Q_uv survives modulo two. No base/base edge appears in the surviving odd cycle, so its prescribed transitive orientation creates no additional condition. The only detected surviving cycle is the full odd cycle R,u, but the cancellation proof must begin with ALL odd cycles, not merely Hamiltonian cycles or triangles.

## High-degree Boolean coefficient interpretation

For fixed Q and frozen attachments, b(T(x)) is a Boolean-valued function of the d orientation coordinates. XOR over the entire d-dimensional Boolean cube is its top multilinear coefficient over F_2. The proof identifies this coefficient with Q_uv. Thus the second-bit faithfulness statement is a concrete high-order coefficient detector: low scalar parity is constant, while this marked contextual coefficient remembers each labelled arc.

This interpretation preserves the exact requirements: the ports are labelled, the frozen endpoint orientation matters, the queried family has k−1 base vertices, and the cube contains 2^(2k−3) points per arc. It does not claim noise robustness, a low-query algorithm, or practical tournament enumeration.

## Explicit triangle-only hostile

Use increasing unordered core pairs to encode a five-vertex tournament; a set mask bit means smaller label→larger label.

- Mask 41 has sorted scores (1,2,2,2,3), H=13, four directed triangles, and two directed five-cycles. Its second Hamiltonian bit is 0.
- Mask 12 has the same sorted scores (1,2,2,2,3), H=11, four directed triangles, and one directed five-cycle. Its second Hamiltonian bit is 1.

The triangle counts agree, but the relevant all-odd-cycle parities differ: 4+2 is even, while 4+1 is odd. Replacing all odd cycles by triangles would give the wrong second-bit interpretation. This is an exact small hostile, not a general characterization of same-score classes.

## Independent executable controls

`odd_cycle_faithfulness.py` explicitly counts directed simple cycles of each odd length using vertex subsets and permutations with the smallest cycle label fixed first. This removes rotation multiplicity without identifying opposite orientations. It checks the odd-cycle XOR and, separately, compares each extension's odd-cycle parity with its Hamiltonian second bit computed by endpoint subset DP.

The completed audit contains 206 fixed-endpoint experiments and 1,864 full extensions: every labelled core of orders two/three and all ordered endpoint pairs, four deterministic backgrounds at order three, and ten order-four sample experiments. It also enumerates both five-vertex triangle-only hostiles above. Normal and optimized outputs agree; freeze recorded in `ODD-CYCLE-FAITHFULNESS-SHA256SUMS`.
