# Checkpoint 2: complete cubic relabeling and explicit cube/octonion maps

**The 3-by-3 magic square of nine distinct integer squares remains OPEN.**
This checkpoint adds independently audited structure to the requested magic,
9-cycle and cube/octonion lanes; it does not identify the user's unresolved
“lambda cube” with a chosen familiar meaning. Session143 began
2026-11-02 at 23:55 UTC; this checkpoint was checked on November3 UTC.

## Complete degree-at-most-three classification

For real `A,B,C,D`, the nine indexed values

    f(k)=A k^3+B k^2+C k+D,  k=-4,...,4,

can be rearranged into a 3-by-3 magic square **if and only if**

    B=0 and (A=0 or C=-15A).

This includes repeated and negative values. In the nonlinear case the nine
values are distinct. A quadratic's nonzero third centered moment gives a
short analytic obstruction; the full cubic result classifies all 945
center/opposite-pair linear systems, then rejects the four nonodd candidate
rays and classifies all 384 signed odd-cubic matching branches exactly.
The survivor `k^3-15k` has an explicit magic grid. It is an ordinary magic
array, **not** a square-entry construction. In particular, squared roots
in any nonconstant real arithmetic progression cannot give a magic square,
regardless of rearrangement, signs or repetitions.

The cycle `4,8,3,7,2,6,10,5,9,4` is translation by4 modulo9. Its ternary
coordinate model needs a nonlinear carry: no relabeling makes it affine on
`F3^2`. Its odd cycle cannot embed edge-preservingly in a hypercube graph.
A graph-picture identification cannot remove that obstruction.

## Genuine cube/octonion correspondences, with the missing data retained

- Signed octonion multiplication is XOR **plus a sign cochain**. The
  multiplicative associator on basis labels is `(-1)^det_F2`; the ordinary
  additive associator is zero or `+/-2` times a basis unit, not `+/-1`.
- Binary-supported subalgebras form the 16-element subspace lattice of
  `F2^3`. Choosing three independent labels selects an eight-corner Boolean
  sublattice with dimensions1,2,4,8; there are28 such choices. **If** the
  intended lambda cube is Barendregt's, this gives only a chosen inclusion-
  lattice comparison, not transport of terms, typing or normalization.
- A multiquadratic field and the rational octonion algebra share an explicit
  eight-character diagonal action and matching fixed-subobject lattices.
  The linear identification is not an algebra isomorphism: its associator
  and multiplication coefficients differ.
- A chosen clockwise labeling of eight perimeter roots gives an explicit
  octonionic left-multiplication operator `J` of order4. After squaring the
  coordinates, `J` becomes the grid half-turn and preserves **all eight
  magic-line equations**. It preserves the squared-value multiset, so cannot
  construct a new magic square or prove existence. It is a linear isometry,
  not an octonion automorphism. Other XOR translations fail that predicate.

A concrete positive-integer-root hostile satisfies every opposite-pair
condition and the expected norm, yet its outer line sums disagree. Norm or
central symmetry alone is therefore insufficient. The nine-cycle cannot act
faithfully on eight basis labels; this does not exclude continuous linear
order-nine actions on an eight-dimensional vector space.

## Read and reproduce

- [Magic/cycle theorem](magic_cycle/THEOREM.md)
- [Cube/octonion proofs and losses](lambda_octonions/PROOF.md)
- [Independent full audit](independent/REVIEW.md)
- [Checkpoint 2 synthesis](SYNTHESIS_CHECKPOINT2.md)
- [Checkpoint 2 verification](VERIFICATION_CHECKPOINT2.md)
- [Checkpoint 1 tournament and C5 bridge](CHECKPOINT1.md)

Run `python3 /path/to/packet/replay_checkpoint2.py` from any working directory.
It replays the portable scientific packet in a fresh relocation, with normal
and optimized execution and explicit corruption controls. No external formal
build, GAP reconstruction, unspecified representation/Hecke map, novelty
priority, or general magic-square-of-squares closure is asserted.

## Explicit projective-point / puzzle / tournament construction

A further **PROVED / FINITE-EXACT** bridge starts with explicitly normalized
points of `PG(2,3)` and the legal Conway-puzzle rule: a hole move swaps its
endpoints and the other two points on their four-point line.

The closed word `[0,1,3,2,4,0]` induces an order-five permutation `g`.
Its selected orbit is `(3,1,4,2,12)`. The closed word
`[0,1,7,12,9,4,0]` induces an order-four permutation `h` with
`hgh^-1=g^2`. They generate an explicit order20 subgroup acting faithfully
on that orbit. Ordering it by `g` and taking forward steps1 and2 gives the
regular five-tournament; `h` rotates its four orientation choices.

This is a concrete projective-points → legal puzzle words → cyclic orbit →
tournament construction, with every word and permutation verified. It is
**not a projective collineation**: `g` sends one named line to a nonline.
It is not an action of all puzzle permutations on this fixed five-set, not
an individual graph's `5:4` automorphism group, and not a Baby Monster
intertwiner. The source choices, unused eight positions, line incidence and
word history are not hidden by the graph output.

See [the coordinate construction](puzzle_coordinates/report.md) and its
[separate independent audit](independent_puzzle/REVIEW.md).
