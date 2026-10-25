# Held-out Paley certificates: P19, P23 and P31

Status: exact finite certificates satisfying the proved sufficient criterion in CRITERION.md; independent audit requested. All three bounded targets succeeded. No other larger orders were explored.

The criterion proves the complete minimum-blocker classifications:

- P19: beta=18, rank134, exactly361 minima (171 sources-pair,171 sinks-pair,19 stars).
- P23: beta=22, rank208, exactly529 minima (253,253,23).
- P31: beta=30, rank404, exactly961 minima (465,465,31).

These are NOT estimates from sampled blockers and do not require any large blocker census. They follow from the all-size sufficient theorem after checking every hypothesis on each named finite tournament. No all-Paley-order rank or classification claim is made.

## Certificate discovery versus verification

Discovery used a fixed seed540000+n and randomized DFS to find a nonconstant Hamilton cycle, then its orbit under the square-affine group. It retained independent cycle differences by elimination modulo101. A bounded maximum of five seed cycles and200000 DFS visits per search was allowed; the FIRST seed orbit reached full rank for each target, with total observed discovery time below three seconds. A forward0-to1 Hamilton path was also found in each case. Discovery was heuristic; this affects how witnesses were found, not their exact verification.

The verifier trusts none of those discovery ranks. It checks every listed cycle and endpoint path, constructs the selected square minor afresh, and computes its determinant modulo101 by Gaussian elimination. The determinant residues are59,44,38 for P19,P23,P31, respectively, all nonzero. Nonzero modulo101 implies a nonzero integer minor and hence rational independence. The universal row-column upper bound then gives exact ranks134,208,404.

It also verifies regular degree, all common in-neighbor counts, the square multipliers and arc orbit, bipartite-support connectivity with a spanning tree, and both strict expansion inequalities across their entire relevant finite ranges. P11's copied rank certificate is rechecked as a control and gives determinant99 mod101 and count121.

## Meaningful hostile

A constant-difference Hamilton cycle and its square-affine orbit are insufficient. Their difference ranks are only r-1:4,8,10,14 for P11,P19,P23,P31, respectively, well below34,134,208,404. The verifier reconstructs this failure, preventing an unsupported claim that any Hamilton cycle orbit automatically supplies full affine span.

The certificate does not identify or count every Hamilton cycle. Only35,135,209,405 listed cycles are needed for P11,P19,P23,P31. The small optional modular probe enumerates all cycles only at P7 and P11, as explicitly requested, to compare ranks modulo2 and101; it does not extend the larger-order search scope.

## Prior-work search and novelty scope

Targeted current-repository searches for P19/P23/P31 Paley blocking numbers, HP blockers, and cycle-difference rank claims located other Paley invariants but no prior minimum-blocker classification of these orders. In particular, unrelated Betti-number uses of beta are not HP-blocking evidence. This is a repository-level claim, not a claim of literature priority.

## Files and replay

- CRITERION.md: reusable theorem with every sufficient hypothesis and proof, including the separate rational-rank and integer-potential arguments.
- certificate-19.json, certificate-23.json, certificate-31.json: explicit cycle, pivot, seed and endpoint witnesses.
- inputs/P11-rank.json: copied earlier control certificate.
- verify.py: exact independent verification of the supplied witnesses and theorem hypotheses.
- discover.py: bounded seeded witness finder; not run during replay. If rerunning discovery, use a separate copy because it writes the certificate files.
- small_modular_probe.py: optional complete P7/P11 rank comparison, with no torsion inference.

Run replay.py --self-test from any directory. It verifies exact file-set/digests, normal/-O outputs, and file-add/content-mutation hostiles. Standard library only. Previous frozen packets remain unchanged.
