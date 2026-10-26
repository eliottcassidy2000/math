# Source receipt and inheritance boundary

Sources read 2026-10-26. This note's twin restriction and line-graph transport
are elementary deductions from the displayed factorizations. No literature
priority is claimed. Source pins are in `SOURCE_MANIFEST.json`.

## Primary external texts

1. Amitabha Tripathi, *On Pythagorean Triples Containing a Fixed Integer*,
   Fibonacci Quarterly 46/47 (2008/09), 331–340,
   [official full text](https://www.fq.math.ca/Papers1/46_47-4/Tripathi.pdf).
   Read introduction and Theorems 1–3, printed pages 331–333. The introduction
   states the coprime, opposite-parity Euclid parameter conditions. Theorem 2
   counts primitive triples with odd specified leg by distributing its prime
   factors between two coprime factors. The present two-node semiprime fibre
   is that classical theorem's specialization, with its factor splits and
   height recorded explicitly. The paper does not provide the twin-filtered
   graph transport claimed here.

2. Lucia Janičková and Evelin Csókási, *Metric properties in Berggren tree
   of primitive Pythagorean triples*, arXiv:2304.05230v1 (2023),
   [full text](https://arxiv.org/pdf/2304.05230).
   Read Section 2, matrices, Lemma 2.4, and its proof (PDF pages 2–3).
   Their matrices A,B,C are this note's U,A,D. Their repeated-A formula
   agrees with the inherited U spine. This is source confirmation of the
   branch convention, not a new theorem of the present packet.

   **Parity reading guard:** the displayed Theorem 2.3 in this version lists
   `m>n` and coprimality but does not explicitly list opposite parity in the
   extracted statement. Those two conditions alone do not suffice:
   `(m,n)=(3,1)` gives `(8,6,10)`. We use the explicit opposite-parity
   condition from Tripathi and the repository's audited chart, and prove
   the needed primitive chart directly. No inference drops that condition.

## Repository sources and what is inherited

- `THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md`:
  ordered Gaussian lift, parity/content boundary, Farey versus Berggren
  distinction. The present argument does not claim those graphs coincide.
- `THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md`:
  U matrices, consecutive Euclid parameters, exact parabolic spine.
- `THM-3756-odd-square-ordinal-berggren-affine-descent.md`:
  complete two-root chart and unique inverse. The chart is reproduced for
  clarity and explicitly not claimed as a new result.
- `THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md`:
  twin-center arithmetic and its finite/global ancestry boundary. Its
  additive cospan is not transported into the Berggren parent operation.
- `aliquot_power_124_20261026/aliquot/README.md`, Theorem A:
  the independently audited unique perfect target of a prime square.
  This supplies the unique-six corollaries; the present packet does not
  reclassify arbitrary perfect-number incoming fibres.

The routed `CORE-PAPERS` files, relevant current frontier and active
guardrails were read first. The targeted mistakes pass recovered the
September odd-square warnings about lost roots and non-preserved legal
edges, plus the primitive-normalization warning at MISTAKE-418. No
chronological log was used as mathematical authority.

## Reproduction and review lineage

The producer runs with only Python's standard library and explicit gates,
not `assert`. Both `python3 -B check_tree.py` and
`python3 -B -O check_tree.py` print the frozen receipt.

The first draft gained one additional proof corollary and one additional
gate before audit freeze: the unique edge of the target graph consists
exactly of the two twin images with the common square `25` on opposite
boundaries. Initial checks numbered 198,534; the frozen producer numbers
198,535. The mathematical graph and subset statements were unchanged.

Independent audit [PASS](../arithmetic/AUDIT_TREE.md) reviews the frozen
proof and checker. It reproduced all 198,535 producer gates and independently
passed 2,892 gates using direct side enumeration and factors of the odd leg's
square. The resulting status/link-only promotion does not change the proof
or checker. Producer computation is FINITE-EXACT; the unrestricted statements
are the algebraic proofs in `RESULT.md`.
