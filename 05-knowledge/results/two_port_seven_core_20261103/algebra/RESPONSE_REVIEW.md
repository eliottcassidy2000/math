# Independent review of the order-eleven full-response controls

**Status: VERIFIED / PASS.** This review covers the two named controls in
`../response/run.py` and `../response/check.cpp`, with the imported sidecar
engine `../self_converse/census.cpp`. It does not independently certify the
complete self-converse census or universal two-port reconstruction.

The reviewed old tournaments are transitive11 and Paley11, upper-triangle
codes36028797018963967 and25822532702336285. Each retains the entire named
S,E,M,C,U,V,D arrays; D has11^4=14,641 entries per control. In the producer's
notation V is the earlier proof's W tensor. The ordered indices agree with
equation(1) and the inherited full three-path formula.

The literal engine visits every old permutation. It cuts at one or two
positions, requires every uncut arc to exist, and counts exactly the
corresponding ordered nonempty path cover. Two consecutive cut positions
give a singleton middle path; these are included. The independent subset
engine counts all named start/end path polynomials by dynamic programming
and multiplies over disjoint complementary supports. Its three-path loop
visits every ordered triple of nonempty disjoint supports covering V.
Equality checks include every array entry, including zeros.

The response evaluator includes adjacent ports, both nonadjacent port
orders, full joint two-path endpoints, and the full four-index D term.
For each control every one of the67^2=4,489 mask pairs with each mask of
weight at most two is independently counted on the actual thirteen-vertex
extension. This battery determines the bidegree-(2,2) polynomial. Nine
explicit high-weight pairs additionally agree under coefficient
interpolation, tensor evaluation and direct Hamilton counting.

The automorphism check is independent of colour refinement: complete
partial-bijection backtracking imposes only S,E and ordered C equalities,
and every completed sidecar automorphism is checked against every actual
arc. Its counts1 and55 agree with the canonical engine. Each control also
has64 specified pseudorandom relabellings checking the full canonical
object. Those64 checks are a sample of relabellings, not all11! relabellings;
the automorphism search itself is exhaustive for each named control.

Integer ranges are safe for the stated universe. A fixed endpoint count
is at most10!, and the sum of all ordered two-path counts is at most
10*11!=399,168,000, below2^32. Thus the census engine's uint32 sidecars
cannot overflow here. Literal three-path tensors and actual13-vertex
Hamilton counts use uint64; their crude factorial bounds fit comfortably.

The wrapper pins all three source files before and after compiling and
running, and checks the expected H values1 and95,095 and automorphism
orders1 and55. An independent optimized replay (`python -O
../response/run.py`) completed successfully and its parsed JSON matched
the frozen `../response/normal.json` exactly. Each control records4,498
actual response checks. The C++ checks use explicit exceptions and remain
active regardless of Python optimization.

Reviewed and replayed source hashes:

- `response/check.cpp`: `4a1e47f72e9a365220f797039d6aeab87c09abdf4f148db0dcf49e70071fee8b`
- `response/run.py`: `481205f6bc93ce2969fdbd62fcc1281677448d924ae3422c8180663544f77e6e`
- `self_converse/census.cpp`: `37312242df9c03eff19d3c9587e8d7a36942bec109e9d6424da5003944e9109a`

No correction is required. The result is two fully checked order-eleven
controls with complete response information; it is not a census of all
order-eleven tournaments or a proof of injectivity at that order.
