# Session 65 independent integration and fifteen-kernel audit

2026-10-26. **Verdict: PASS in the stated scope. No unresolved blocker.**

This audit reviews the final `session65-root/SYNTHESIS.md` and `CONNECTIONS.md`
against the sealed producer/audit siblings and independently verifies the new
`fifteen-kernel` corollary. It does not repeat their long native runs, redo
external literature audits, or promote source-inspected formal proofs to local
builds. The repository, Git refs, and every producer file were read-only to this
lane. Only this audit directory was written; the parent owns publication.

## 1. Corrections adopted before the final pins

Two narrow prose issues were reported and corrected by the parent:

1. The kernel report originally said “all 105 unordered-pair orbits.” There
   are 105 unordered **pairs**, partitioned into 65 orbits for duads under a
   transposition and 57 for synthemes. It now correctly says all unordered-pair
   orbits covering the 105 pairs. No formula, code or receipt changed.
2. The two root summaries now explicitly retain `r>=2` in the multiperfect
   arc-family statement. The identity itself also works at `r=1`, but its
   even-family obstruction does not: `D_1=1`. The detailed, previously audited
   multiperfect report already had the correct hypothesis.

All source pins and computations below cover the corrected root bytes.

## 2. Independent finite-carrier and polynomial proof

### Complete carriers and their action

The checker imports no producer implementation. It generates all 15 duads of
six labels, then inspects all `C(15,3)=455` triples of duads; precisely 15 triples
cover each of the six labels once. These are all synthemes, not a chosen list
of fifteen conveniently matched objects.

For an independent incidence control it constructs the even-subset quotient
by the full six-set. Its 16 classes form the binary symplectic four-space;
the 15 nonzero classes are exactly the duads. Enumerating all isotropic
zero-sum triples of nonzero classes gives exactly the 15 syntheme planes,
with 45 point-plane incidences. This verifies the finite geometry used by
the report; the arithmetic identification with Richelot kernels remains the
**cited interpretation** in the inherited Session 52 source, not a new
calculation of curves or isogenies. The kernel carrier is the syntheme/
Lagrangian carrier, not the duad/nonzero-vector carrier.

For the labelled transposition `(01)`, the duad fixed points are `{0,1}`
and the six pairs in the other four labels, giving `(f,b)=(7,4)`. A fixed
syntheme must contain `{0,1}`: otherwise the partners of 0 and 1 are different
other labels and swapping 0/1 changes those pairs. The remaining four labels
have three matchings, giving `(f,b)=(3,6)`. Different fixed-point counts
exclude every identity-`S6`-equivariant bijection between these two carriers.

### Exact signed-orbit factors

For an involution with `f` fixed vertices and `b` swapped vertex pairs, its
unordered-pair orbits have three types:

| Orbit type | Count | Local defect polynomial |
|---|---:|---|
| Fixed endpoints | `C(f,2)` | `2` |
| Swapped endpoints | `b` | `2z` |
| Positive doubleton | `m=b(f+b-1)` | `2(1+z^2)` |

Here a defect counts a failed comparison for one unordered pair. Every swapped
endpoint pair forces one failure. Every doubleton has either zero or two;
there is no coupling between distinct pair-orbits. Multiplication therefore
gives exactly

`2^[C(f,2)+b+m] z^b (1+z^2)^m`.

Substitution proves the reported formulas and their sharp minima:

| Actual action on 15 vertices | Positive singletons | Negative singletons | Doubletons | Orbits | Polynomial | Minimum |
|---|---:|---:|---:|---:|---|---:|
| `(01)` on duads | 21 | 4 | 40 | 65 | `2^65 z^4 (1+z^2)^40` | 4 |
| `(01)` on synthemes | 3 | 6 | 48 | 57 | `2^57 z^6 (1+z^2)^48` | 6 |
| `(01)(23)(45)` on duads | 3 | 6 | 48 | 57 | `2^57 z^6 (1+z^2)^48` | 6 |

In every row, singleton lengths plus twice the doubleton count equal 105.
Evaluating each polynomial at 1 gives `2^105`. Its lowest nonzero coefficient
counts minimizers: `2^65` or `2^57`, respectively.

The new checker independently builds the actual edge permutation and endpoint
reversal bits, decomposes every unordered pair into an orbit, and enumerates
**only each orbit's local assignments**. These have length one or two; precisely
210 local assignments are checked per action, 630 in total. It then multiplies
the resulting coefficient arrays and compares all 106 coefficients, including
zeros, to the closed form. This is complete algebraic enumeration via a proved
factorization, **not an enumeration of `2^105` whole tournaments**.

### Literal certificates and hostiles

For all three actions the supplied 15 adjacency-bitmask rows are separately
checked for bounds, zero diagonal, and exactly one directed arc on every pair.
Every failed invariance comparison is recomputed from those literal rows and
the independently derived vertex permutation. The reported defect lists match
exactly, consist of the swapped vertex pairs, and have sizes 4, 6 and 6.

The auditor also constructs different literal minimizers from local assignments;
all their rows, edge-orbit inventories and complete coefficient arrays are saved
in the receipts. Reversing one edge in a positive doubleton increases the defect
by exactly two, providing an explicit local hostile to an invariance claim.
Since a transposition alone already forces positive defect, no tournament on
either carrier can preserve the full natural `S6` action.

The triple-transposition duad character matches the transposition syntheme
character, as stated. This one matched character is **not** a proof of the full
outer automorphism or outer-action equivalence. The stronger construction is
explicitly inherited: Session 52 builds six pentads and maps a syntheme to its
two containing pentads. Its proof and the distinction between ordinary and
pentad actions were inspected in the saved README; its historical checker was
not rerun here. No identification with the Ogg-prime set follows.

## 3. Integration claim-to-evidence ledger

| Root claim | Checked supporting packet | Boundary retained |
|---|---|---|
| Supplied macro orbit halts; 47,176,870 native steps and 4098 ones | `session65-audit-weakcollatz/AUDIT.md`, native receipts | One exact run, not every-seed halting or a new maximality proof |
| Binary 3-adic Cantor full shift, transitive tournament shadow, density-one integer halting | Same audit, sections 2–4 | Integer residue counting, not a mere Haar-measure inference; exceptional seeds are not excluded |
| Total four-state morphology tables, guarded word macros, tape/head/state/clocks | `session65-audit-cryptids/REPORT.md` | Total-table nonhalting is elementary; no unproved asymptotic morphology theorem |
| Tape-window metric code, actual clone changes, time `2k+1` at distance 1 | Same audit, sections 2–3 | Ordered finite words plus required sidecars; not a configuration-graph isomorphism, K-move simulator, or endpoint-distance runtime bound |
| Frobenius minimum defect, p=37 test, and cited modular-curve genus equality | `session65-audit-fifteens/AUDIT.md` | Optimizes over all orientations; does not supply a canonical arithmetic tournament or reprove the Ogg classification |
| General permutation minimum equals number of even vertex cycles | Same audit, section 2 | Signed unordered-pair action, not permutation sign alone |
| Compatible all-ones ternary branch and ordinary-integer realization failure | Same audit, section 4 | Parity and 3-adic exponent coordinates both retained; neither digit conjecture proved |
| Historical 25+ cross-field atlas | `session65-twentyfive/REPORT.md`, sealed Goldbach audit provenance | Historical mixed-evidence method atlas, not 25 proved transports or a uniquely recovered 25-card enumeration |
| Goldbach wiki/Lean table renaming and halting iff counterexample | `session65-audit-goldbach/goldbach/REPORT.md` | Finite literal isomorphism differs from cited semantic equivalence; 2016 had 27 states, later source 25, newer 21 has separate status |
| Reachable padding changes program size but preserves computation at macro boundaries | Same audit | Description size is not semantic invariant; native clocks change |
| Multiperfect arc identity and odd valuation ledger | `session65-audit-goldbach/multiperfect/REPORT.md` | `r>=2`; necessary odd conditions are not sufficient; restricted even family cannot encode odd witnesses |
| Directed atlas and logical implication examples | `session65-audit-reductions/REPORT.md` | Fixed sentences, uniform languages, halting, beeping, simulation and resource reductions remain different types |
| Deterministic graph and automorphism obstructions | Same audit, graph section | Literal outdegree obstruction and faithful-symmetry obstruction, not impossibility of generic graph-isomorphism reduction |

The root claims do not exceed those scopes. In particular: NP-complete languages
are decidable; a shared halting target gives no converse reduction; the JC/DC
chain is `JC_(2n) => DC_n => JC_n`, not fixed-rank equivalence; and no universal
classification making all touched open problems mutually dependent is claimed.
The machine's fixed-word/tape carrier, unbounded configuration carrier, halting
predicate and resource-bounded language are not silently identified.

Formal-source inspection and observed upstream CI are explicitly distinguished
from local Lean/Coq builds and axiom audits. The inherited source receipts do not
become new execution evidence here. The native-counter clock gap, unknown
integer realization, historical-method status and open converse arrows remain
visible rather than being filled by numerology or untyped analogies.

## 4. Provenance and actually executed verification

`SOURCE_PINS.json` pins **88** reviewed root files and supporting reports,
receipts and seals. The inherited Session 52 README is

`session65-fifteens/inputs/05-knowledge__results__session52_saturation_residues_metagraph_20261024__fifteen-kernels__README.md`,

SHA-256 `747ab7952914b4a9002ecc33ac98d6714596f9f1eb7cfc1cc9ace240930a8c07`.
It matches the fifteens source manifest already pinned by the completed audit.
The inherited repository source commit is
`6a8ad5d3bf01a232da4de41bcbff68390a9337ee`; the moving publication checkout was
not consulted or mutated.

The independent checker also validates **230 existing seal rows** covering the
five producer/audit lanes and the two root corollaries. This includes producer
report identity links to the completed audits and their corrected freezes.
These are integrity checks, not repetitions of the computations described there.

Actually executed successfully:

```sh
python -B session65-audit-synthesis/independent_check.py > session65-audit-synthesis/normal.json
python -O -B session65-audit-synthesis/independent_check.py > session65-audit-synthesis/optimized.json
cmp session65-audit-synthesis/normal.json session65-audit-synthesis/optimized.json
python -B session65-root/fifteen-kernel/check.py > session65-audit-synthesis/root-corollary-normal.json
python -O -B session65-root/fifteen-kernel/check.py > session65-audit-synthesis/root-corollary-optimized.json
cmp session65-audit-synthesis/root-corollary-normal.json session65-audit-synthesis/root-corollary-optimized.json
cmp session65-audit-synthesis/root-corollary-normal.json session65-root/fifteen-kernel/normal.json
```

All comparisons are byte-identical. Every gate in the independent checker uses
explicit exceptions and remains active under `-O`. The root's small authored
checker was read before replay. No downloaded external software, historical
Session 52 implementation, long simulation, formal build, or network program was
executed. No permission denial or access blocker occurred.

**Stopping boundary:** integration consistency, source identity, exact new finite
corollary, literal witnesses and optimization-safe replay are established. The
missing converse reductions, integer realization and headline open conjectures
are mathematical frontiers, not incomplete steps silently certified by this audit.
