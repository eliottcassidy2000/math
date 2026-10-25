# Integral observers, cancellation spectra, and root arrangements

2026-10-25. **Status: first audited research checkpoint.** Scoped proofs and
finite controls are separate. This packet does not claim historical novelty,
an F1 model, a universal torification classification, an integer-square solution, or an
LRC/HN consequence. Proposed unifying theories are not promoted by these checks.

## Read the results

- [Cancellation spectra and explicit 50-vertex identification](cancellation/RESULT.md)
  · [independent audit](audits/cancellation/REPORT.md).
  The regular elementary-abelian spectrum and witness are **recovered concurrent
  results**, not new discoveries: [incoming spectrum](../unprotected_ports_cancellation_20261025/valuation/README.md)
  and [incoming witness](../unprotected_ports_cancellation_20261025/counting/POSITIVE_JET.md).
  The negation isomorphism is checked entrywise. The all-odd-action,
  integral-precision and synchronized-action interfaces retain their stated scope.
- [Open, singular and resolved twist mass](twists/RESULT.md)
  · [root-stack presentation](twists/ROOT_STACK.md)
  · [fixed-family zeta and effectivity boundaries](twists/ZETA_AND_EFFECTIVITY.md)
  · [independent audit](audits/twists/AUDIT.md).
  In characteristic different from 2 and 3, the three fresh means are respectively
  `q²−6q+12`, `q²+q+1` and `q²+4q+1`. Characteristic 3 has separate formulas;
  characteristic 2 is excluded. Fixed twists over extensions are not fresh twists.
- [Concept/status board](CONCEPT_BOARD.md) · [reproduction and scope](PUBLICATION_CHECKS.md)
  · [source-copy receipt](FIRST_SOURCE_VERIFICATION.json) · [exact manifest](MANIFEST.json).

## Portable path map

The audited proof and review files are copied **byte for byte**, including their
original scratch-directory references and historical handoff wording. Within
this packet use the following mappings; no scratch path is needed to run tests.

| Original directory | Packaged directory |
|---|---|
| `session63-cancellation` | `cancellation/` |
| `session63-twists` | `twists/` |
| `session63-recovery` (finite twist checker and outputs only here) | `recovery/` |
| `session63-audit-cancellation` | `audits/cancellation/` |
| `session63-audit-twists` | `audits/twists/` |

Repository-relative source citations still refer to the repository root.
Source commit pins remain historical pins; they are not silently advanced.
The broader recovery map and proposed framework are reserved for a later,
separately audited checkpoint, not claimed as part of this one.

## Reproduce

From the repository root, set `P=05-knowledge/results/session63_integral_observers_root_arrangements_20261025`:

```sh
python3 -B "$P/cancellation/check.py"
python3 -B "$P/cancellation/verify_source_recovery.py"
python3 -B "$P/audits/cancellation/audit.py"
python3 -B "$P/audits/cancellation/provenance_audit.py" --repo "$PWD"
python3 -B "$P/twists/verify.py"
python3 -B "$P/recovery/check_twist_histograms.py"
python3 -B "$P/audits/twists/independent_checks.py" --reference "$P/recovery/twist_finite_normal.json"
python3 -B "$P/verify_packet.py"
```

Repeat with `-O` after `-B`. Mathematical checkers print deterministic JSON;
the corresponding saved normal and optimized outputs must match exactly.
All these checks use Python's standard library and explicit runtime guards.
The proof scope is stated in the notes, not inferred from a finite census.
