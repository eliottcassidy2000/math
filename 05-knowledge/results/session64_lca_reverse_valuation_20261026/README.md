# Guarded ancestry, exact endpoint repairs, and smaller counting carriers

2026-10-26. **Final independently audited Session 64 packet.**
[Research synthesis](session64-root/SYNTHESIS.md) ·
[exact connection contracts](session64-root/CONNECTIONS.md) ·
[independent synthesis audit](session64-audit-synthesis/REPORT.md).

## Current results and boundaries

1. **Sharp local endpoint theorem: n >= 8.** Every directed closed blue marked
   Hamiltonian path with endpoint outdegrees `(n-5,4)` has three black literal
   K moves to its exact labelled target, without clone or other source filters.
   The six order-12/14 menu failures are completely repaired; independent
   finite audits extend the full source families through orders 11, 10, 9 and 8.
   The order-seven target with nine Hamiltonian paths, all blue, proves the
   threshold is sharp for this local statement—not whole black-core connectivity.
2. **Production reverse valuation for ordinary A000568.** The actual unrestricted
   ordinary DAG retains every structurally reachable state and exact weighted
   merger. A suffix-minimum potential computes normalized binary residues
   without first retaining the extra absolute baseline bits. Arithmetic-budget
   filtering and forced unit-zero-cost skip contraction have separate proofs.
   At N=100 the direct contracted graph has 10,547 input states (excluding the
   common terminal) and 47,147 edges;
   normalized K=2 uses 64 forward updates, excluding graph/potential setup.
   No asymptotic or wall-clock speedup is asserted.
3. **Smaller signed-core observation bank.** Newton coordinates at `2,4,8,...`
   expose precision losses `0,1,3,6,10,...` in the full evaluation image.
   This is an exact absolute-residue mixed-precision presentation for the
   signed clone-excess lane; the ordinary lane was already scalar. It is not
   permission to divide a truncated residue by powers of two.
4. **A precise LCA connection.** A Newton coordinate's precision loss equals
   the sum of low-digit-trie ancestor depths with preceding evaluation nodes.
   Actual unit residues and allowed transfer operations remain necessary.
5. **Exact Antihydra comparison.** The accelerated maps share valuation coding
   on correctly punctured odd 2-adic domains, but positive integrality is not
   preserved. Integer-ray, chronological-debt and barrier conditions remain
   explicit. Collatz convergence and blank-tape Antihydra halting remain open.

## Proofs, certificates, and independent audits

| Lane | Producer result | Independent check |
|---|---|---|
| Six order-12/14 repairs | [Proof and recipes](session64-repairs/README.md) | [Complete order-12/13/14 audit](session64-audit-repairs/AUDIT.md) |
| Order 11 | [Full source family](session64-repairs/next11/README.md) | [224-cell audit](session64-audit-eleven/AUDIT.md) |
| Orders 8–10 | [Order 10](session64-repairs/next10/README.md), [orders 8/9](session64-repairs/boundary/README.md) | [40/80/140-cell audit](session64-audit-lower/AUDIT.md) |
| All-size sufficient bridge | [Degree and port criterion](session64-repairs/bridge/README.md) | [Independent bridge audit](session64-audit-bridge/AUDIT.md) |
| Ordinary reverse potential | [Implementation](session64-counting/REVERSE_IMPLEMENTATION.md) | [Reverse audit](session64-audit-counting/REVERSE_AUDIT.md) |
| Forced-skip contraction | [Exact path bijection](session64-counting/SKIP_CONTRACTION.md) | [Contraction audit](session64-audit-counting/CONTRACTION_AUDIT.md) |
| Signed-core Newton bank | [Mixed precision](session64-counting/NEWTON_QUOTIENT.md) | [Independent Newton audit](session64-audit-newton/REPORT.md) |
| Newton/LCA precision | [Proof and exact controls](session64-root/newton-lca/REPORT.md) | [Newton audit](session64-audit-newton/REPORT.md) |
| Valuation coding and native guards | [Dynamics proof](session64-reverse-valuation/REPORT.md) | [Independent dynamics audit](session64-audit-dynamics/REPORT.md) |

Literal ports, complete tournament rows, anti-reversal witnesses and exact
clone block intersections are retained. Real clone mergers and splits occur
in the lower-order certificates; no quotient path is silently lifted.
The bridge criterion is sufficient, not an existence theorem for its word.
The separate order-ten old-menu diagnostic is explicitly diagnostic only;
the later full lower-order certificate supplies the repair theorem.

## Portable reproduction

Set `P=05-knowledge/results/session64_lca_reverse_valuation_20261026` from the
repository root. Every source is copied byte for byte under its original
session-directory name. The nested `session64-repo/` contains only pinned
reference inputs, not another Git checkout. Thus all mathematical and source
seal checks below run entirely from this packet's bytes.

```sh
python3 -B "$P/session64-repairs/check.py"
python3 -B "$P/session64-audit-repairs/check.py"
python3 -B "$P/session64-repairs/next11/check.py"
python3 -B "$P/session64-audit-eleven/check.py"
python3 -B "$P/session64-repairs/next10/full_check.py"
python3 -B "$P/session64-repairs/boundary/full_check.py"
python3 -B "$P/session64-audit-lower/check.py"
python3 -B "$P/session64-audit-bridge/bridge_check.py"
python3 -B "$P/session64-audit-bridge/n10_check.py"
python3 -B "$P/session64-counting/reverse_ordinary.py" --repo "$P/session64-repo"
python3 -B "$P/session64-counting/contracted_ordinary.py"
python3 -B "$P/session64-audit-counting/audit_reverse.py"
python3 -B "$P/session64-audit-counting/audit_contraction.py"
python3 -B "$P/session64-counting/newton_modular.py" --repo "$P/session64-repo"
python3 -B "$P/session64-audit-newton/check_independent.py" --repo "$P/session64-repo" --producer "$P/session64-counting/newton_modular.py"
python3 -B "$P/session64-audit-newton/check_literal.py"
python3 -B "$P/session64-root/newton-lca/check.py"
python3 -B "$P/session64-reverse-valuation/check_valuation.py"
python3 -B "$P/session64-audit-dynamics/audit.py"
python3 -B "$P/session64-audit-synthesis/verify_sources.py"
python3 -B "$P/verify_packet.py"
```

Repeat with `-O` after `-B`. Normal/optimized outputs and the designated saved
JSON receipts must agree exactly. The provenance checker verifies source
seals; it does not replace the separate native mathematical replays.

## Provenance, links, and publication scope

[Exact final source-copy pins](THIRD_SOURCE_VERIFICATION.json), the original
native manifests/checksums, and [the file manifest](MANIFEST.json) retain byte
identity and scope. [Original-context link index](SOURCE_LINK_INDEX.md) supplies
working repository links for relocated historical/source-relative prose;
those internal archived links are not represented as nested live navigation.
The original counting checksum list is now complete, including Newton inputs.

[Publication checks](PUBLICATION_CHECKS.md) distinguish proofs, finite controls
and transport integrity. Earlier checkpoint receipts record the n >= 12 and
n >= 11 milestones; this final packet supersedes those thresholds only through
the separately audited lower-order extension. Original producer notes remain
unchanged as source history. Full black-core connectivity, positive Collatz
convergence, Antihydra blank-tape nonhalting, and a globally minimal ordinary
state recurrence are **not claimed**.

## Concurrent source credit

The [concurrent LCA/ports/reverse packet at its pinned commit](https://github.com/eliottcassidy2000/math/tree/2a9a9c342e5940f0e50b97fba87a7cfb7a24de7b/05-knowledge/results/lca_ports_reverse_20261025)
independently reaches the shared n >= 12 repair closure, ordinary production
reverse gauge and punctured 2-adic coding/guard conclusions. No exclusive or
historical-priority claim is made for those overlaps. The sharp n >= 8 head,
forced-skip contraction and signed-core Newton presentation are separately
scoped here. Alternative full-cell arc selectors, the Mahler translation and
the fixed-congruence reset obstruction are recorded in the
[bounded incoming comparison](session64-incoming-review/REPORT.md).
That report inspects source/proof scope and pins all 41 incoming files; it is
not a claim that the incoming computational engines were rerun in this session.
The concurrent packet is preserved unchanged and is not silently added as a
new dependency of the already sealed proofs.
