# Verification and publication boundaries

**Final root replay: PASS.** Fourteen principal executions (seven verifiers,
each normal and optimized) ran from independently relocated temporary
packets. Eleven manifests covered138 hash entries. Three artifact-integrity
negative controls rejected a changed source, missing source and stale receipt.
Every normal/optimized receipt hash agrees. `SESSION_REPLAY.json` records
actual receipt hashes and runtimes. Reproduce with `python3 verify_session.py`.

The only requirement is Python3 standard library; neither a source checkout
nor external data nor network access is needed for the delivered checks.
Unpublished exploratory C++ coverage files and binaries are deliberately
excluded. They are not dependencies of the new theorems.

## Executed lanes

- `constructor/verify.py`: actual fresh-block matrices, literal graph/cut
  reconstruction, correct arbitrary-strip coefficients, capped all-word
  closures and mixed-cut repairs. Normal/optimized/fresh-directory producer
  receipts match; independent closure uses a different traversal.
- `descent/run_checks.py`: full low-order insertion and induced-pair controls,
  independent strong Hamilton-cycle cover for the gain-four classification,
  actual infinite-family source and separate witness counts. All-order
  kernel, factorial and ray-density conclusions are proved in the manuscript.
- `two_apex/check_word_recovery.py`: exact coordinate identities, recovery and
  rejection of nonactual profiles, direct subset-DP controls and retained-apex
  SCC compression. The existing rank17 realization is credited, not new.
- `intervals/verify_arc_budget.py`:35,000 direct/fast/knapsack comparisons,
  1,100 labelled tournaments through order five, strong-floor and actual
  sparse-family controls, exact finite six-arc counting gate.
- `run_entropy/check.py`:95 exact word-count identities,1,198,372 deletion/
  projection checks and10,385 marked tournament kernel controls. Rational
  partial sums certify the elementary numerical inequalities used in proofs.
- `radix_probe/verify.py`: exhaustive products through four operations,
  exact characteristic polynomials, every integer radix root and eigenspace,
  and the actual radix80 survivor. The finite scope is part of the result.
- `independent/run_checks.py`: eleven separate normal/optimized consumers,
  including direct graph/cut checks, independent full capped closures,
  independent core-budget comparisons, alternative Sturm root isolation,
  literal first-iterate radix80 count, and actual induced-pair/family checks.

## Corrections and scope guards

The preliminary arbitrary-strip cover coefficient3^(r-1) was an unpublished
failed extrapolation. The corrected2*3^(r-1)-2^r+1 is proved, independently
checked, and accompanied by a literal hostile. No earlier published positive
repair from that failed formula is used.

The profile decoder preserves all future responses; singleton-SCC compression
preserves scalar H and presentation budget only. The word3333 distinguishes
these interfaces. The all-length strip and fixed mixed-chain grammars are
distinct; neither is promoted to universal acceptance. Positive eigenrays are
not assumed realizable, signed intercepts are allowed, and zero-digit cases
are not silently excluded by the eigenvalue1 necessary condition.

Classical OCF, Moon and Busch results remain CITED inputs. Independent
computation corroborates the all-order proofs; it does not replace them.
The global odd-H spectrum completeness conjecture remains OPEN.
