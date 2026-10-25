# Integral observers, cancellation spectra, and root arrangements

2026-10-25. **Status: audited results plus a separately proposed framework.** Scoped proofs and
finite controls are separate. This packet does not claim historical novelty,
an F1 model, a universal torification classification, an integer-square solution, or an
LRC/HN consequence. Proposed unifying theories are not promoted by these checks.

## Start with the synthesis and the proposed theory

**[Integrated synthesis](SYNTHESIS.md) · [Proposed unifying framework](PROPOSED_UNIFYING_THEORY.md)**

**PROVED / RECOVERED:** the concrete, scoped compatibility results and hostile
examples below. **PROPOSED FRAMEWORK:** integral carrier + operation/context
bank + realizability condition, reviewed for mathematical coherence, not proved
as a universal equivalence or classification. Its named lemmas have their own
proofs; its broad organizing program and next steps remain proposals.

![Keep the construction and distinguish its observers](figures/integral-observer-map.svg)

[Current source reconciliation](CURRENT_SOURCE_UPDATE.md) supersedes the old
search boundary without rewriting the [historical recovery map](recovery/README.md).
[Late concurrent intake](LATE_INCOMING.md) records preserved subsequent work
without importing its broader unaudited claims.

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
- [Integral specialization and signed-orbit base change](specialization/RESULT.md)
  · [Fourier lattice and histogram effectivity](fourier-effectivity/FOURIER_EFFECTIVITY.md)
  · [independent integral/Fourier audit](audits/specialization/REPORT.md).
  These retain the signed-permutation hypothesis, distinguish invariant modules
  from multiplication (including the C3 hostile), and separate an integral
  Fourier lattice from positivity and geometric realizability.
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
| `session63-recovery` | `recovery/` (historical source map/pins plus finite checks) |
| `session63-audit-cancellation` | `audits/cancellation/` |
| `session63-audit-twists` | `audits/twists/` |
| `session63-specialization` | `specialization/` |
| `session63-root` | root synthesis/theory; `fourier-effectivity/`; `figures/` |
| `session63-audit-specialization` | `audits/specialization/` |

Repository-relative source citations still refer to the repository root.
Source commit pins remain historical pins; they are not silently advanced.
The broader recovery map is now archived with its current reconciliation.
The proposed framework is prominently linked above and remains distinct from
proved concrete results; final readbacks certify scope and coherence, not a
universal classification theorem.

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
python3 -B "$P/specialization/check.py"
python3 -B "$P/fourier-effectivity/check_fourier_effectivity.py"
python3 -B "$P/audits/specialization/verify.py"
python3 -B "$P/verify_final_presentation.py"
python3 -B "$P/verify_packet.py"
```

Repeat with `-O` after `-B`. Mathematical checkers print deterministic JSON;
the corresponding saved normal and optimized outputs must match exactly.
All these checks use Python's standard library and explicit runtime guards.
The proof scope is stated in the notes, not inferred from a finite census.

The source Fourier note retains its historical audit-request header. The
[completed independent audit](audits/specialization/REPORT.md), section 6,
records its subsequent proof readback and normal/optimized replay. This is
status reconciliation, not a change to the mathematical source.

The figure generator requires Matplotlib. Give `MPLCONFIGDIR` and
`XDG_CACHE_HOME` explicit writable directories outside this packet before
running `python3 "$P/figures/make_observer_map.py"`. Its default output is
only the published SVG; `--png` additionally makes a visual-QA preview. The
preview PNG was inspected but is intentionally not a published packet file.
