# Audit, universes and reproduction

2026-10-28. **PASS, same-agent proof review and alternate implementations.**
No independent researcher or Lean verification is claimed. Statements labelled
PROVED have arguments in their named files; finite tests have the scopes below.

## Declared computations

| Receipt | Universe and controls | Explicit guards |
|---|---|---:|
| ports.json | Literal A1,A2,C1,C2,C3 groups; six signed weight assignments each; all ordered two-port source/target tuples through depth2; sampled triples at depth3; all24 run orders and all4 internal cuts in the inherited eight-vertex hostile; actual16/18-vertex complement contexts; twelve W2 Hamilton paths with all81 source automorphisms; equal cycle-index hostile;81/27 and normal-basis orientation guards | 16043 |
| paley.json | Full automorphism searches and all candidate proper nontrivial modules for Paley3,7,11; every base marking triangle; seeded depth1..3 character/control maps; exact moment recurrences through depth8 and order6 | 3555 |
| counting.json | Unit-tail kernel through n60; ordered ports r0..8, all distinct-port channels; separate direct odd-partition evaluation through26; published A000568 prefix through19 | 839 |
| arithmetic.json | All positive odd residues modulo2^K for K2..14; both signed sheets; every valuation word of length1..4 over1..3 with three integer lifts; four supplied cycles; same-source signed offspring law | 36448 |
| crt.json | Every legal cyclic word at n9,12,15,21,27; actual reader, rotation, reflection and projection fibres; local roots and unit-clock formula | 1431212 |
| audit.json | Separate meet-in-the-middle Hamilton counts with brute-permutation controls; all83349 literal Paley7 depth-two centralizer permutations; independent divisor-Burnside chart counts; all64 O/(8) elements; full rational distributions for two spine laws through depth3; four-vertex operation-type and literal81/27 triple controls | 83374 |

The total 1571471 counts explicit guard executions, not independent theorems
or independent samples. The 81/27 structural definitions and oracle are
inherited; they were not rewritten and silently renamed as new work.

`ports.py` uses the inherited Hamilton dynamic-programming oracle. `audit.py`
constructs the sixteen-vertex complement directly arc by arc and counts paths
by joining two disjoint half-paths. It independently gets 629677591 and
634905785. This tests the consequential isomorphism separation, rather than
only an intermediate endpoint statistic.

The Paley7 enumeration verifies nine multiplier fibres of size9261 inside
the order83349 centralizer. Identifying the identity fibre with all marking
holonomy is the all-size group proof in PORT_LAWS P5. The enumeration does not
claim to list all elements of the much larger source automorphism group or
all of its marking loops. The Paley11 depth-two order is likewise a proved
formula, not a twelve-billion-element enumeration.

## Proof audit and equality boundaries

- **CLOSURE C1:** sufficiency and necessity both use the same fixed aligned
  pair and the declared named-position cut menu. Original lengths and the cut
  word are retained. The singleton decoder recovers the source and original
  forest, not merely an unrelated copy of the partner. General unnamed-menu
  minimality remains OPEN.
- **CLOSURE C2:** the old E2 theorem explicitly excludes refinement, so the
  counterexample does not retract it. All source children in the actual
  sixteen-vertex context are strong; the isolated child partner B need not be.
  Sixteen is not asserted minimal.
- **PORT_LAWS P1:** the hostile compares two actual subgroups, not two full
  centralizers. Their equal full fixed-set distributions and unequal transport
  are independently checkable from six explicit permutations.
- **PORT_LAWS P3:** full wreath action versus full tournament automorphism
  group is distinguished. Criticality requires a transitive root action;
  the general reversal law records whether the root anti-involution has zero
  or one fixed vertex. No free independence is inserted between paired copies.
- **PORT_LAWS P5:** normality of <K> is proved by twisted commutators. The
  root-translation quotient argument proves the reverse inclusion in the
  commutator subgroup, upgrading a lower bound to an equality. Surjectivity
  from the centralizer to H^k is constructive. The p=3 exception is handled
  by its inherited theorem rather than an invalid nonidentity multiplier.
- **CRT R1-R3:** odd readers are bijective in the checked universes; the even
  reader is not. O/(8) is characteristic8, not F64. Cardinality27 does not
  turn Z/3 x Z/9 into a cyclic group. The translation-lift obstruction uses
  unequal full fibres, not just their total number.
- **CLOCKS:** Haar source cylinders are not integer-orbit frequency claims.
  Signed sheet, endpoint, odd-step count and shortcut clock remain distinct.
  Four supplied cycles are verified; completeness is OPEN.
- **COUNTING:** the unit-tail core and ordinary Burnside formula are inherited
  or classical. Named-port factors do not supply actual metagraph transitions,
  and these n60 channels do not supersede the earlier ordinary n300 evaluation.

The all-period golden-reader theorem and its earlier single-toggle repair are
credited in CRT. The March tribonacci priority correction MISTAKE-645 and the
October n11/n12 refutations are preserved in RECOVERY. The incoming selector's
global-converse/native-complement operation-name error was repaired in place
and logged in MISTAKES under this packet's slug. Its correctly typed proof and
counts survive; the four-vertex witness is now part of the alternate audit.

After rebase onto1925bc3385, the hash check correctly flagged three changed
inputs: ENDPOINTS, endpoint_kernel.json, and FULL_HOLONOMY. Their diffs were
read. E2 and the eight-vertex data were unchanged; E6 and the improved ternary
bound were integrated. The full optimized replay was repeated before resealing.
The corrected incoming canonical checker also passed its2619 checks.

## Reproduce and pin

All six normal runs completed. One aggregate `python -O` replay then matched
every saved JSON receipt exactly. The explicit exception guards remain active
when assertions are disabled. No network, external solver or random seed from
the environment is needed for the computation.

Run from the repository root:

```powershell
python -X utf8 -B 05-knowledge/results/complement_port_crt_20261028/seal.py
python -X utf8 -B 05-knowledge/results/complement_port_crt_20261028/replay.py
python -X utf8 -O -B 05-knowledge/results/complement_port_crt_20261028/replay.py
```

The seal pins this packet, loaded inherited Python dependencies, the literal
endpoint witness, and principal proof inputs. It excludes the changing global
results index. Updating a dependency requires a deliberate re-audit and reseal;
do not silently treat a hash failure as a new verified run.

External literature inspected: Davis's formula and sequence attribution via
[OEIS A000568](https://oeis.org/A000568), and the classical independent-branch
spine construction in
[Lyons, Pemantle and Peres, Section2](https://arxiv.org/html/math/0404083v1).
No novelty claim is made for Burnside, size bias, CRT, group abelianization,
or finite-state observational equivalence. The concrete applications and their
limits are the work recorded here.
