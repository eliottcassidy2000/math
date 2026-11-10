# Independent audit: skew matching restoration and Hamilton-path transport

2026-11-10. **PASS within the stated finite universe and general action hypotheses.** No new global extremality, all-prime classification, or Chan–Pak closure claim follows.

## Independent finite results

`verify_independent.py` imports no producer files, uses only the Python standard library, and reconstructs each incidence structure from nonzero squares modulo p. Run from the repository root:

    python3 05-knowledge/results/session240_seven_eleven_symmetry_20261110/audit_design/verify_independent.py

The recorded three JSON objects are in `output.jsonl`. All verification conditions use an explicit `check` helper that raises `RuntimeError`; they remain active under Python optimization. Normal and `python3 -O` runs agree on all JSON fields except elapsed time. Both modes reject `--self-test-failure` with exit code 1 and the intended RuntimeError. See `output.optimized.jsonl` and `optimization_checks.json` for replay receipts. The script initially considered using a graph library, which was unavailable; the delivered implementation is dependency-free and uses a different group-enumeration route from the producer.

| p | compatible matchings | design group order | matching stabilizer | H | H/stabilizer |
|---|---:|---:|---:|---:|---:|
| 7 | 8 | 168 | 21 | 189 | 9 |
| 11 | 12 | 660 | 55 | 95095 | 1729 |
| 19 | 1 | 171 | 171 | 1172695746915 | 6857869865 |

All matchings form one orbit in each case. At 7 and 11, the induced matching action has kernel of size one, and its canonical-matching stabilizer has orbits [1,7] and [1,11]. This independently confirms faithfulness and 2-transitivity. At 19 the matching action is trivial; its kernel is the entire order-171 group (the JSON leaves the optional action-kernel computation null at 19).

### Why these are exhaustive, independent computations

Matching domains initially contain every block not incident with the assigned point. Arc consistency deletes a candidate only if another point has no candidate satisfying both distinct block labels and the skew orientation constraint. A branch assigns each remaining value in a smallest nonsingleton domain. Thus every valid bijection survives exactly one branch. This differs from the producer's sequential matching assignment and its incidence-partition group search.

For each completed matching m, the audit enumerates all tournament isomorphisms alpha from the canonical tournament to T_m, pruning only violated arcs between already assigned vertices. Each such alpha lifts uniquely to an incidence automorphism with beta(i)=m(alpha(i)). Conversely, any side-preserving incidence automorphism transports the canonical matching to exactly one m and occurs in precisely that enumeration. Every lift is independently checked on every point-block incidence. Positive lift count for each matching proves transitivity; the total lift count proves the design-group order without assuming affine automorphisms or imported group orders.

Hamilton paths are counted by an ordinary subset/end-vertex dynamic program, initialized at every vertex. This independently reproduces the 7/11 totals and additionally computes the 19 hostile control. The exact 19 value is not required for the claimed exceptional bridge.

## General proof audit

For an incidence structure with equally many points and distinct blocks, the stated two skew conditions are exactly what is needed to reconstruct a tournament from a point-to-block bijection. Distinctness ensures a point permutation determines a block permutation, and makes the action faithful on points. The matching action beta f alpha^{-1} preserves the skew conditions. Its stabilizer is precisely the reconstructed tournament's automorphism group.

The Hamilton-path construction must be the **tagged disjoint union** of pairs (f,h), not an untagged union of vertex sequences. Under a design automorphism (alpha,beta), transport (f,h) to (beta f alpha^{-1}, alpha h). Incidence preservation proves that the transported sequence is a directed Hamilton path for the transported matching.

Fixing such a pair fixes the entire ordered Hamilton sequence pointwise. Hence alpha is the identity, and distinct blocks imply beta is the identity. The action on the tagged union is therefore free. When the matching set is transitive, every orbit meets the canonical matching fiber, and two points in that fiber are in the same design-group orbit exactly when they are in the same matching-stabilizer orbit. This proves the claimed induced-set and quotient descriptions.

Numerically, 8*189/168=189/21=9 and 12*95095/660=95095/55=1729. The matching enlargement cancels exactly the larger group denominator. Directly dividing the fixed tournament's path count by 168 or 660 is invalid: the remainders are 21 and 55. No stronger divisibility condition or new attainable/missing raw Hamilton-path values arise from this transport alone.

## Source and producer review

Reviewed the delivered `paley/REPORT.md`, `paley/verify.py`, the provenance synthesis, and the current MISTAKE-651 entry via the repository object. MISTAKE-651 corrects the Kantor reference and requires Sylow conjugacy to bridge difference-set equivalence and design isomorphism. It does not license identifying tournament and design groups. The producer explicitly respects this distinction and labels the PSL identifications inherited, rather than inferring an isomorphism type from orders 168 or 660. This independent audit establishes orders/actions only; it does not re-prove Kantor's classification or the cited PSL identifications.

The producer's claims for all 32 order-11 circulant connection sets, rotational controls, and endpoint-refined counts have their own verification artifacts. This audit independently covers the Paley matching fibers at 7/11/19, their lifted full incidence automorphisms, faithfulness/2-transitivity at 7/11, total Hamilton-path counts, and the transport proof. No untested universal statement is promoted from these finite checks.
