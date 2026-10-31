# Mixed-radix square-sum rays: three compatible operators

**Status: PROVED, with finite universal certificates and independent replay.**
This packet constructs new explicit prescribed prefixes, not a new proof of
ordinary square-sum Hamiltonicity and not arbitrary prescribed-path embedding.
Tournament Hamiltonian-path count infinite completeness remains **OPEN**.

## Main result

The three operators are:

| Code | Radix | Offset | Fixed prefix | Added vertices |
|---|---:|---:|---|---:|
| A | 25 | 12 | (1,3) | 2 |
| B | 25 | 8 | (1,3,6,19) | 4 |
| C | 49 | 4 | (1,3,6,19) | 4 |

For an operator with radix Q, offset c, and prefix H, apply
`F(P)=H || (Q*p_j+(-1)^j*c)` in one-based order. For every nonempty finite
outer-to-inner code word w, the prefix `P_w=F_w(1,3)` extends **unscaled**
to a square-sum Hamiltonian cycle on every interval [1,N] with

    N >= (85 * product(radices in w) - 1) / 2.

Different codes produce a prefix-decodable family. Every infinite schedule
produces a distinct simple square-sum ray, each finite prefix of which has an
eventual all-interval cycle extension. These are not claimed Hamiltonian rays
on all positive integers. Compatible-code antichains also give explicit lower
bounds on distinct finite Hamiltonian cycles, not simultaneous embeddings.

See [the complete theorem and proof](TERNARY_SCHEDULE_THEOREM.md), including
uniform vertex/root formulas, exact coding, threshold induction and cycle-count
scope. B and C have different rays; the smaller radix does not improve the
threshold for the old C-ray itself.

## Why the operators compose

All consume and return the same weak B pair: interval paths of consecutive
orders, starting 1,3, ending at 8/10 according to parity, and agreeing on the
position parity of every common label. The new B operator does **not** require
a separated head to share its tail's shift inside a word. Their own colors
must agree separately between output words.

- [Split-head semantic relaxation](HEAD_TAIL_RELAXATION.md).
- [Exact piece-overlap classification](PIECE_OVERLAP_COMPATIBILITY.md).
- Independent family audits: [A](A/INDEPENDENT_AUDIT.md), [B](B/INDEPENDENT_AUDIT.md).

Earlier frozen certificates imposed a stronger sufficient condition and remain
valid. The relaxation removes only unnecessary within-word equalities, not the
cross-output compatibility required by subsequent recursion.

## Files and reproduction

[REPRODUCE.md](REPRODUCE.md) provides portable standard-library commands,
source-relative dependencies and independent audit interfaces.
`construct_cycle.py` checks a full streamed cycle for any certified code/order;
`verify_package.py` covers all mixed pairings and depth-three controls.
A/B's 200 new JSON certificates are included. The old large base bundle and
C's existing certificates are referenced by hash, not duplicated.

Independent audits have 51,027,430 checks for A and 51,028,631 for B, including
all universal states, every required base, hostile controls and cycles through
N=16,601,562. Normal and optimized receipts agree byte-for-byte. The universal
induction, not a large finite sample, proves the all-size theorem.

## Additional exact consequences

- `python schedule_frontier.py 1000000000` gives a certified lower bound of
  145 distinct cycles, using the exact maximum antichain inside this coded
  family. This is not an optimality claim about all square-sum cycles.
- [Independent piece-overlap lemma](PIECE_OVERLAP_LEMMA.md) and
  `test_overlap.py` classify precisely which piece shifts must agree.
- [Affine-ray arithmetic appendix](affine_rays/AFFINE_RAY_LEMMAS.md) gives
  general formulas and parameter classifications. Its broader arithmetic
  families are not claimed to have Hamiltonian extensions without catalogs.

Run `test_frontier.py`, `test_overlap.py`, and `affine_rays/test_identities.py`
for the corresponding exact finite controls. Each accompanies a proof; none
substitutes sampled arithmetic for a universal extension certificate.

The [independent ternary audit](ternary/AUDIT.md) passes 31,877,061 checks,
including 30 mixed cycles, closed-form decoding and exact frontier controls.
It also verifies that the maximum coded-prefix antichain count is
Θ(X^α), where `2*25^(-α)+49^(-α)=1` and `X=(2N+1)/85`.
The upper bound concerns this prefix-certificate method only, not all cycles.
