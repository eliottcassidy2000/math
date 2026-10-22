# Publication audit and incoming-source reconciliation

## Subsequent additive closure (Session25)

A separate argument, now fully written and independently audited in ADDITIVE_PRIME_POWER_PROOF.md, proves s(n+phi(p^r))=(2/p)s(n) modulo p^r for every odd p, r>=1, and n>=2r-1, with the same law for D. Thus the incoming additive conjecture is now proved by this new filtered-Hurwitz/Newton argument. The earlier caution that the multiplicative theorem alone does not settle it remains correct. The new proof also gives arbitrary-unit-base closure and the sufficient connected-graph bound n>=2r. No minimal-threshold or literature-priority claim is made.

The main prime-modulus and cycle-stabilizer proofs survive audit. The general multiplicative prime-power contraction theorem is independently proved in PRIME_POWER_PROOF.md.

Incoming commit 3596c5d0b1990b05ae960e980a863956b47cc6b5 and MISTAKE-651 now own the overlapping canon repairs: the 262-class denominator, connected-graph tail law, prime-square equivalence, zero-class boundary, and corrected Mersenne table. Preserve those repairs and their provenance. Our old proposal is withdrawn, renamed HISTORICAL_SUPERSEDED_CANON_CORRECTIONS.txt, and must not be applied.

The decisive correction is a threshold error: connected labelled graphs have the same prime-modulus law from index 2 onward. Only their initial singleton term fails. Thus ordered condensation is not a distinguishing cause of quasi-periodicity. INCOMING_AUDIT_RECONCILIATION.md gives a direct contraction check consistent with incoming Lemma C and preserves its literature credit.

Our all-prime-power result compares mp^r with mp^(r-1). The incoming higher-power conjecture instead compares arbitrary n with n+phi(p^k). The multiplicative result supplements the incoming work but does not settle that additive conjecture. At squares, two independent proofs now coexist.

The checked residues remain unchanged: 9 passes with residue 7, 63 passes with residue 2, and 341 fails with discrepancy 22. The exact closed-presentation congruence remains a useful connection to individual tournament symmetry. No aggregate-to-spectrum or aggregate-to-wiring inference is asserted.

No repository file was edited or canon ID reserved. Residual singleton/multiplier-class wording qualifications are recorded separately rather than restoring a stale patch. Incoming larger computations and bibliography are source-attributed, not claimed as independently reproduced here.

## Reproduction

Run `python3 audit.py --repo /path/to/math` from any directory. Alternatively omit `--repo` when the script is inside the repository or alongside its `math` checkout; the script discovers the repository from its own location and then working-directory ancestors. The repository must contain immutable source checkpoint `f002fc425336`, retained as the denominator-audit input rather than treated as current canon.

All verification conditions use explicit runtime checks and remain active under `python3 -O`. Normal and optimized runs produced byte-identical output. The portable discovery route was separately checked. `audit.out` and `audit.optimized.out` are the frozen outputs. `SHA256SUMS.json` records artifact hashes and the incoming canon hashes.

The checks cover 24,000 quasi-periodicity instances; independent path, cycle, and reflection enumeration on 62 connection sets; exact residues for 9, 21, 63, 341, and 1023; the census denominator; prime-power grids; and 990 connected-graph tail cases with ten initial-term obstructions. Literature novelty remains unclaimed.


## Full CRT extension

Section8 of ADDITIVE_PRIME_POWER_PROOF.md removes the squarefree/cube-free restriction: for each odd p^r exactly dividing n, reduce n to the representative t in [2r-1, phi(p^r)+2r-2] and test p^r|D(t). The2-primary condition is automatic. Direct independent recurrence agrees with this criterion for every one of the831 composites through1000. This is a theorem about aggregate D(n)=s(n)-c(n), not a metagraph invariant.
