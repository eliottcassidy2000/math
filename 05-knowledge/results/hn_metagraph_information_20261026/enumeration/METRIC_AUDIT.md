# Independent audit of planar whole-frame transport

**Verdict: PASS for MT1--MT3 as stated in
[metric_transport/RESULT.md](../metric_transport/RESULT.md).** No producer
proof defect found. This audit does not promote a geometric endpoint section
to a proper six-colouring of an arbitrary unit-distance graph.

Reviewed independently: all three proof steps, the signed-permutation target
convention, the sufficiency direction for arbitrary multigraph networks, and
the producer checker. The producer checker was run successfully. The separate
checker here imports no producer functions.

## Proof and scope audit

1. Three distinct planar direction lines give rank two and a one-dimensional
   dependence with every coefficient nonzero. Its two sign choices give
   positive barycentric coordinates after normalization, hence exactly two
   balanced endpoint triangles. They are antipodal. No equal-length assumption
   is needed for this part if the direction lines remain distinct.
2. The centre is the mean of the whole six-point frame. A whole-frame isometry
   therefore preserves the centre and its intrinsic antipodal pairs. Its
   endpoint action is a signed coordinate permutation. Strict containment of
   the centre in an endpoint triangle is preserved, so the balanced pair is
   covariant without making arbitrary global choices of local markings.
3. Complementing all three bits changes parity. Each balanced pair thus has
   exactly one choice of each parity. Any solution of the parity seam equations
   uniquely lifts inside these pairs and satisfies every full seam equation.
   This works simultaneously on all edges, loops, and parallel edges. The
   proof does not assume a one-cycle graph or test loops separately.
4. A parity-zero map preserving an antipodal pair must fix both members:
   swapping them would reverse parity. The two displayed abstract holonomies
   have disjoint fixed sets, hence cannot both preserve the balanced pair of
   any one admissible planar frame. Realizing each at a different marking does
   not realize their simultaneous loop constraints at one vertex.

The stated exclusions are necessary. General relational gadgets may not be
whole-frame isometries; the theorem does not exclude those mechanisms. The
parent's qualification about arbitrary unit-distance graphs is retained.

## Independent exact controls

`check_metric_independent.py` uses exact rational Gaussian elimination for the
origin's barycentric coordinates. It tests metric admissibility directly from
all pairwise six-point distances, recovering the endpoint permutation from the
three coordinate derivatives of the eight-state map. It uses neither the
producer's Gram formula nor its cofactor function.

- 9,216 signed transport tests across all sign markings of three frames:
  a regular unit flower, different rational unit directions, and a frame with
  unequal spoke lengths. All accepted maps preserve the independently found
  balanced pair.
- 27,648 networks on two vertices, with two parallel edges and a loop. These
  have two independent cycles, so they test a common section rather than
  individual loop solvability. Full sections are enumerated among all 64
  assignments; parity sections are enumerated separately. The full-section
  histogram is `0:20736, 2:4992, 4:1728, 8:192`; exactly the nonzero cases pass
  parity, and every passing network retains its balanced section.
- A sharp dimension hostile: in three dimensions the unit frame
  `+/-e1,+/-e2,+/-e3` admits every signed permutation as a whole-frame isometry.
  In particular both parent's incompatible parity-zero maps are isometries.
  The two-loop network passes parity but has no common eight-state section.
  Planarity, and its unique dependence sign pair, is therefore essential even
  for actual unit-spoke isometries. This is outside the stated theorem.

Normal and optimized outputs are byte-identical. Their receipts pin the exact
producer proof, producer checker, and independent checker bytes.

```text
python 05-knowledge/results/hn_metagraph_information_20261026/enumeration/check_metric_independent.py
python -O 05-knowledge/results/hn_metagraph_information_20261026/enumeration/check_metric_independent.py --output 05-knowledge/results/hn_metagraph_information_20261026/enumeration/metric-independent-optimized.json
```
