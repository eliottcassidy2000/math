# Proof and computation audit

**Verdict: PASS for the stated scoped results. Full prime/7p realization
remains OPEN.** This packet used separate agents for realization,
attachment, and prime-graph work, followed by cross-agent review. The
coordinator proved B1–B2 and the general R0 extension, reviewed every final
argument, supplied the all-odd-square modulo16 obstruction, and replayed
the integrated checks. No external literature-priority claim is made.

## Proof obligations and independent paths

| Result | Independent review and decisive controls |
|---|---|
| R0 arbitrary-block sorting | Spectrum agent checked coordinator algebra; prime-graph agent independently counted1,098 actual profiles and1,584 named contextual swaps. Both signed commutators, suffix bounds, strictness and empty/full equality boundaries checked. |
| R1–R2 vertex-transitive increment and strong dyadic families | Prime-graph agent's [audit](audit_realization.md) and [literal checker](audit_realization.py) retain the same vertex names, compare every changed arc, and count by DFS and all permutations. r=0 is explicitly outside the formulas. |
| B1–B2 boundary divisor | Spectrum agent's [audit](boundary_independent.md) uses80 source/cut cases,69,342 actual paths and104,316 rotation images; includes a non-ordered cyclic module source and a source rotation broken globally upstairs. |
| B4 strong sharpness | Spectrum agent supplied the new-edge congruence and strong repair; coordinator, prime-graph agent and attachment agent checked the proof. The attachment agent's [separate audit](attachment_boundary_audit.md) uses different vertex labeling and DP. Both independent graph checks reach17 vertices. |
| Triangle two-birth detector | Coordinator checked the coefficient convention, color-transposition argument, eight-case reduction, deletion-deck distinction and linear-rank upper bound; the realization agent supplies the additional [independent audit](audit_attachment.md). |
| Five prime groups and square-carry map | Prime-graph agent enumerated edges both by pairs and by square roots, recovered actual prime witnesses, and counted the H49/H343 anchors by DP and all permutations. Coordinator checked SCC cancellation, the square reflection direction, the special carry at7, and the modulo16 obstruction at2. |

The square-sum graph diagram was opened and visually inspected after
generation. Every vertex, group, isolated vertex, bipartition count and
caption agrees with the exact receipt. Prime-support group-order statements
are CITED, with primary ATLAS checks recorded in `prime_groups.md`; the
numerical graph and witness claims are independently computed here.

## Integrated executable replay

From the repository root, run each listed script with `python3 -I -B` and
again with `python3 -I -B -O`. Compare stdout with its frozen receipt:

| Script in this directory | Frozen stdout |
|---|---|
| `realization_probe.py` | `realization.json` |
| `audit_realization.py` | `audit_realization.json` |
| `attachment_probe.py` | `attachment_receipt.json` |
| `audit_attachment.py` | `audit_attachment.json` |
| `boundary_budget.py` | `boundary_budget.json` |
| `boundary_independent.py` | `boundary_independent.out` |
| `prime_groups_probe.py` | `prime_groups.json` |

The coordinator's integrated replay found byte-identical normal, optimized
and frozen output for every row. `prime_groups_probe.py` also regenerates
the graph's JSON and SVG, and its PNG when Pillow is available. The other
mathematical checks use only the standard library. All checks use explicit
exceptions, not optimization-sensitive assertions.

Attachment controls enumerate all33,866 cuts of all1,099 labelled parents
through order5 and all128 two-birth triangle attachments. Realization
controls compare the full state on723 short block words and2,408 contextual
swaps in addition to the independent profile audit. Boundary controls
enumerate all cuts to one, two and three ordered triangles. These are
declared finite universes and control the proofs; none is treated as an
unbounded realization census.

## Validity boundaries retained

- Actual source states are never replaced by arbitrary formal nonnegative
  coordinates. A target prime needs a supplied source and compatible cut.
- A pair of signed transfer matrices is a difference of constructions;
  its signed linear combination is not itself declared a tournament.
- Equal one-step H responses, equal child H, and isomorphic children are
  distinct assertions. Ordinary isomorphisms must transport named ports.
- Symmetry groups act on untouched path strata even when they are not
  automorphisms of the extended graph. Old induced arcs are fixed in B1.
- Strong sharpness proves exact valuation, not primality: H7009=43*163
  is the explicit hostile at the first two-apex repair.
- The five prime colors are neither equitable classes nor supplied
  automorphism orbits. Their quotient edges cannot be lifted freely.
- The square-carry map preserves the secondary square-sum relation and
  exact path products; it loses strongness and primality of the H labels.
  The24-odd-prime embedding does not extend the prescribed dilation at2.
- Incoming fixed-ray and fixed-wrapper obstructions are scoped to their
  operation languages. The finite adaptive-root success is not induction.

No correction to an established theorem or new theorem-ID reservation was
needed. Existing order-five and weighted-cut results are explicitly
credited. The stronger all-size claims here are mathematical proofs;
finite controls and inherited numerical banks are labeled separately.
