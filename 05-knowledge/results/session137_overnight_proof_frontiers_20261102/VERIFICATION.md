# Session137 verification and reproduction

## What the evidence establishes

All-order claims are supported by mathematical proofs and independent review, not extrapolation from the finite checks. Finite classifications retain their exact universes and completeness reductions. A normal/optimized comparison, compiler optimization comparison, sanitizer run, or relocated replay of the same program is a reproducibility control; it is not an independent algorithm.

| Evidence | Scope |
|---|---|
| [Global proof](discriminant/GLOBAL_PROOF.md), [first independent audit](independent_global/GLOBAL_DETERMINANT_AUDIT.md), [further referee](global_referee/REPORT.md), [cycle-side audit](cycle_carriers/DETERMINANT_SIDE_AUDIT.md) | All-order nonempty determinant inequality, equality cases, weighted extension, and proper-cycle defect; each audit states its checked scope |
| [Counting construction](global_referee/COUNTING_GAPS.md), [independent audit](independent/COUNTING_GAPS_AUDIT.md) | Four particular #P functions; implementation-under-test audit computes expected counts independently |
| [Polynomial stability audit](independent/WEIGHTED_POLYNOMIAL_AUDIT.md), [cofactor audit](independent/WEIGHTED_STABILITY_AUDIT.md) | All-order stated coefficients, n=2 counterboundary, n=3 identity, n=4 sharpness, and rounding normalization |
| [Combined lane audit](independent/AUDIT.md), [width universality audit](independent/WIDTH_UNIVERSALITY_AUDIT.md), [second width review](width_phase2/independent/AUDIT.md) | Exact contexts and all-modulus positive-context construction, constructors, carrier reduction, coverage, and harmonic family; no exact-spectrum closure |
| [Bounded-interior audit](independent/BOUNDED_INTERIOR_AUDIT.md), [coverage continuation audit](independent/COVERAGE_CONTINUATION_AUDIT.md) | Finite-base normal form, canonical card reconstruction, selected literal-cut row, and formal-capacity firewall |

## Frozen-source integrity

The remaining lane integration preflight checked **nine hash manifests containing 242 entries**, with no missing files or hash mismatches. Entries can overlap between manifests; this is not a count of distinct proofs or files. It also checked **34 reviewed/inherited source pins**: 13 independent reviewed files, 10 carrier sources, 6 constructor sources, and 5 harmonic inherited copies. The 233-file remaining lane selection preserves the existing directory layout and includes its manifests.

The final weighted-stability manuscript was reviewed separately after that preflight:

```text
discriminant/WEIGHTED_STABILITY.md
SHA256 6bd43447a41e092721d0705069721e46f164eed2428c1655333209eb5fd5ce7c
```

The first checkpoint remains preserved at commit
`252dd4a839b945f0b2be77b4f181dc6241aa6868`; its 46 published blobs were verified against intended bytes. Its mathematical and replay scope is recorded in [CHECKPOINT1.md](CHECKPOINT1.md) and [REPLAY_CHECKPOINT1.json](REPLAY_CHECKPOINT1.json). Publication-byte verification and repository route checks are distinct from mathematical proof review and from CI.

The root `SHA256SUMS` is the aggregate packet manifest. Lane manifests retain their frozen receipts and source identities. `independent/PUBLIC_FILES.txt` describes its payload; `PUBLIC_FILES.txt` and `SHA256SUMS` themselves are also included. `width_phase2/probe.json` is explicitly **exploratory**, not a theorem dependency. Compiled binaries are reproducible run products, not proof artifacts.

## Fresh portable replay status

The initial portable validation completed **26/26 executions**: thirteen final checkers in normal and optimized Python modes, each byte-identical to its frozen expected receipt. It covered modular contexts, local universality and its CRT referee, the three constructor checks, cycle carriers, three harmonic checks, and the two weighted-stability checks plus independent width universality. This is additional to the first checkpoint's own fourteen principal replay executions.

The expanded portable replay has completed the following blocks:

| Block | Completed executions | Result |
|---|---:|---|
| Coverage Python: seven programs, normal and `-O` | 14 | PASS; frozen stdout matches |
| Independent Python: twenty-one programs, normal and `-O` | 42 | PASS; frozen stdout matches |
| Independent C++: four programs, optimized and UBSan | 8 | PASS; frozen stdout matches |
| C++ builds for those routes and the full k=7 producer | 9 | PASS |
| Regenerated coverage JSON sidecar comparisons | 12 | All byte-identical |
| Full k=7 portable producer census | 1 | PASS; byte-identical same-source replay, not an independent census |

The expanded portable replay is **terminal PASS: 65 scientific replays and nine builds**, with no stderr or UBSan diagnostics. All 233 frozen originals remain unchanged; all 242 manifest entries also verify in the portable tree after replay. The initial 26 and expanded 65 contain overlapping checker coverage and must not be reported as distinct mathematical test universes. The full k=7 same-source run took 1,038.263 seconds, used 345,696 KiB peak RSS, and reproduced the complete 301,724-byte receipt exactly. This does not make the finite producer census independently implemented.

[Complete replay report](verification/FULL_REPLAY.md) · [Terminal summary](verification/full-replay-summary.json) · [Execution ledger](verification/full-replay-results.json) · [Postflight](verification/full-replay-postflight.json) · [Initial 26 replays](verification/replay-results.json) · [Source-integrity report](verification/RECONCILIATION.md) · [Final manuscript review](verification/WEIGHTED_MANUSCRIPT_AUDIT.md)

## Essential finite scope boundaries

- **Global determinant:** the root checks cover 33,867 ordinary and 4,325 weighted matrices. A separate ordinary implementation checks 34,380 tournaments, 42,532 root insertions, and 4,748,109 literal paths. The independent cycle-side route covers 7,412 complete classes through order eight and 470 literal parent/cut recurrences. These are distinct universes, not additive headline counts.
- **Weighted polynomial stability:** 856 rational matrices; 3,627 weak and 3,621 strong deletion-root inequalities; all 27 ternary order-three identities; coefficient recurrences through n=200. These controls supplement the symbolic all-order proof.
- **Width-three producer:** 1,950 actual polynomial gate checks, 226 complete positive-word replays, and 24,688 simultaneous modular landings. Independent literal-factor and positive-context controls are separately documented.
- **Natural carriers:** exact augmentation through source order eight covers all carriers with at most fifteen cycle-vertices by the proved source bound. The finite source bound, not the observed class counts alone, makes the 29-type classification complete.
- **Coverage k=6:** the full producer quotient census uses 35 classes, 210 deletion-orbit rows, 860,160 profiles, and 27,525,120 represented cuts. The independent selected-row route checks every one of 131,072 literal cuts for one order-seventeen parent of order-eighteen children; it is not a full independent quotient census.
- **Coverage k=7:** the producer **FINITE-EXACT** census covers 353 strong quotient classes, 2,471 deletion-orbit rows, 40,484,864 profiles, and 2,591,031,296 represented cuts, with zero nonoriginal H-restoring cuts. A full same-source portable replay does not make this independently replicated.
- **Paley7 control:** all 212 literal-child counts reconstruct the 211 constant/linear/quadratic coefficients using forward-extension DP. This is producer cross-route evidence. The selected `-O3`/`-O0` comparison is an optimization control, not an independent algorithm.

See [coverage audit scope](coverage/AUDIT_STATUS.md) for all-order versus finite statements. No full `-O0` k=7 census is claimed. All-order strong-quotient scalar zero-rigidity remains open.

## Reproduce without overwriting frozen receipts

Run from the root of this packet. First check the aggregate manifest and lane manifests, then copy the packet before running programs that generate sidecars:

```sh
sha256sum -c SHA256SUMS
for d in width width_phase2 constructors coverage cycle_carriers independent; do
  (cd "$d" && sha256sum -c SHA256SUMS)
done
(cd harmonic_owners && sha256sum -c MANIFEST.sha256)

W=$(mktemp -d /tmp/s137-replay.XXXXXX)
cp -R . "$W/packet"
cd "$W/packet"
```

The checks used Python 3.12.14 and g++ 14.2.0. Primary lane programs use the Python standard library and/or C++17. The historical independent discriminant checker additionally requires SymPy; its recorded environment used version 1.14.0. Python proof gates use explicit exceptions and remain active under `-O`.

For a typical Python producer, write new receipt names and compare both modes with the frozen expected output:

```sh
(cd width &&
 python3 -B context_residue.py > replay.normal.json &&
 python3 -B -O context_residue.py > replay.optimized.json &&
 cmp replay.normal.json normal.json &&
 cmp replay.optimized.json normal.json)
```

Coverage's seven Python routes should run in the documented order because some regenerate JSON sidecars:

```sh
(
  cd coverage
  for s in audit phase2 fullstar deck_leaf_audit symmetry_edit zero_closure_controls card_rigidity_audit; do
    python3 -B "$s.py" > "$s.replay.out" || exit 1
    python3 -B -O "$s.py" > "$s.replay-O.out" || exit 1
    cmp "$s.replay.out" "$s.out" || exit 1
    cmp "$s.replay-O.out" "$s.out" || exit 1
  done
  sha256sum -c SHA256SUMS
)
```

The full finite k=7 producer replay is:

```sh
(
  cd coverage
  g++ -O3 -std=c++17 -Wall -Wextra -Werror zero_escape.cpp -o zero_escape &&
  ./zero_escape 7 quotients7.txt > zero_escape7.replay.out &&
  cmp zero_escape7.replay.out zero_escape7.out
)
```

The original run took about sixteen minutes with roughly 335 MB peak array storage; these are environment observations, not runtime guarantees. This storage figure is memory, not required output-file size.

Complete lane commands and expected receipts are in [width verification](width/VERIFICATION.md), [phase-two verification](width_phase2/VERIFICATION.md), [constructor verification](constructors/VERIFICATION.md), [carrier verification](cycle_carriers/VERIFICATION.md), [coverage audit status](coverage/AUDIT_STATUS.md), [harmonic theorem and reproduction](harmonic_owners/README.md), and the [independent audit packet](independent/AUDIT.md). For the global routes use [first-checkpoint instructions](CHECKPOINT1.md) and the [counting verifier instructions](global_referee/COUNTING_GAPS.md).

## Claim boundary

The packet asserts neither formal proof-assistant certification nor external priority. It does not give an efficient exact Hamiltonian-path counter, full spectrum completeness, all-order scalar zero-rigidity, or LRC(14). Source hashes certify identity; exact computation certifies the stated finite tests; the manuscripts carry the all-order arguments.
