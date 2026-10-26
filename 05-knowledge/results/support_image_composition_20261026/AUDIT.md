# Proof and reproduction boundary

**PROVED:** the result-local statements S1–S4, I1–I5, E1–E4 and C1–C5
have producer proofs and independent cross-audits. The abstract composition
criteria have a separate proof review. Classical ancestry is CITED;
finite censuses are FINITE-EXACT. This is repository research review, not
external peer review or a literature-priority certification.

| Lane | Producer gates | Independent gates | Independent proof review |
|---|---:|---:|---|
| Four-marker sandwich | 3,878 | 36,688 | [REPORT](audit_sandwich/REPORT.md) |
| Hamilton-support image and black incidence | 76,618 | 1,682 | [AUDIT](audit_image/AUDIT.md) |
| Signed even-graph cards | 17,990 | 25,468 | [PROOF_AUDIT](audit_even/PROOF_AUDIT.md) |
| Arithmetic clocks and finite realization | 489,131 | 89,397 | [PROOF_AUDIT](audit_arithmetic/PROOF_AUDIT.md) |

That is **740,852 explicit gates across eight implementations**, each also
run with Python optimization enabled. Normal and optimized receipts are
byte-identical for each program. Optimized reruns and producer replays do
not count as additional independent evidence. Gate counts measure checks,
not distinct mathematical objects or theorem strength.

## Universes and hostile controls

- **Sandwich:** producer checks all joint rooted tournament types through
  core order five and the fixed-context collisions. The independent test
  additionally enumerates all 32,768 labelled tournaments at encoded
  order six: 800 are images, whereas 2,960 meet only the degree-cardinality
  condition. It independently checks all order-two/three rooted cores
  against the five shorter nonempty two-sided contexts. The all-order
  proof handles the sink/source triangle rotations and their exact
  automorphism kernel. An earlier degree-description typo was corrected
  during independent audit, before PASS.
- **Image:** producer exhausts all oriented partial graphs through order
  five and compares Hamilton-word enumeration with weighted subset DP.
  The independent implementation uses edge sets and permutations for all
  partial graphs through order four, then separately checks the first
  five-vertex false candidate, ordered path covers, the joint four-cycle
  hostile, ordinary-type transport sums at orders three/four and 16 padded
  hub paths. Its all-order proof audit uses the inherited promised-image
  reconstruction theorem explicitly. It does not claim the four-cycle
  hostile passes the local decoder.
- **Even graphs:** all 209 graph types and 77 tournament types through
  order six, with all 663 rooted attachment types to parents through
  order five. The independent verifier examines 6,169 marked card maps.
  Direct tournaments, graph sign characters and the classical partition
  formula supply separate routes. K2 root transport, P3's odd deletion
  parent, Eulerian K3 and unmarked deletion cards test the boundaries.
  Formula evaluation above order six is not a graph census there.
- **Clocks:** independent enumeration includes all 2,418 path cylinders
  for M≤12 and L≤4, every divisor pair M|N≤24, fixed-height counts,
  parity/native-time examples and signed boundary prefixes. The producer
  has the larger grid recorded in its receipt. The proof establishes
  arbitrary-modulus finite composition and **outgoing** path lifting;
  it makes no incoming-degree assertion. The all-zero/all-one paths
  separate finite positive realization from one infinite positive source.

The three cross-auditors wrote their checkers without importing producer
implementations. Root wrote the image producer and the integration notes;
its image proof and summary received another lane's review. See
[COMPOSITION_AUDIT](audit_integration/COMPOSITION_AUDIT.md) for the abstract
criteria and the integration manifest for final summary bindings.

## Source and byte integrity

Inherited results are pinned to research base
`4950589d056a64f8807774cf7da893dcafcd63a0`; historical verification uses
`git show` at that revision, not a subsequently edited truth surface.
The incoming THM-4808 retains its pre-audit status in the arithmetic
ancestry; the elementary vector-collision family is proved directly.
Every claimed source is named by slug/file as well as any theorem number.

Each proof candidate's hash is preserved through status/link-only
promotion. The cross-reviewer checked that reversing those replacements
recreates the audited candidate. Producer and audit manifests retain
final file bindings without mutual manifest hash cycles. The root
manifest closes the packet file universe, checks all lane/review bindings,
and records historical source pins and eight replay suites.

[MAINTAINED_PATCHES.json](MAINTAINED_PATCHES.json) records exact changes to
HYP-9416 and navigation. Only the all-core four-marker clause is closed;
other hypothesis clauses are retained. No new theorem ID is reserved.

Run from this directory:

```
python3 -I -B verify_packet.py --live
```

For a later historical replay, omit `--live`: maintained navigation may
legitimately have advanced. Use `--integrity-only` to omit computational
reruns. Byte integrity binds a review to an artifact; it does not prove
mathematical validity by itself.

## Open consequences retained

The image criterion is exact but exponential. The A000568 augmentation
uses graph states and marked transports, not earlier scalar counts alone.
The fixed-port hub and exact black-line incidence do not prove global
black-core connectivity. The clock/sieve construction gives finite-height
certificates, not twin-prime infinitude, Collatz convergence, abc estimates,
a new F1 equivalence, a machine halting classification or a reduction
between famous open problems.
