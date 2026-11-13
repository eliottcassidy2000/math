# Scope and independent verification

**PROVED scoped deductions / CITED inputs / FINITE-EXACT receipts.** No Lean
formalization or proof of Chan–Pak is claimed. The target-dependent selector,
universal seed-preserving preparation, and the earlier473 repairs remain open.

## Universes and controls

|Receipt|Declared universe|Decisive control|
|---|---|---|
|triangle_probe|Every rooted strong class of orders3–7; all2^(2(n−1)) cross matrices for a triangle sharing its root|2,618 sources,396 ordinary classes,10,014,992 matrices; each minimum respects H'≥3H|
|local_obstruction|16 C3/C3 matrices and256 C3/R4 matrices|Literal permutations, endpoint DP, and the odd-cycle collection formula agree; surplus S and orphan count O separately verify H'=H(A)H(B)+S+O|
|atomic_ears|All2,618 rooted strong classes through7|Physical root survives each induced parent restriction; every atom is strong and nearly transitive; independent Python integer path DP verifies each count descent|
|peeling|Triangle-ear parsing of all396 ordinary strong classes through7 and all their triangle exits|Accepts1,1,6,35,352 classes by order; literal witnesses for all1,436 odd targets9–2881 except21|
|neutral_entry|The explicit25→25 replacement, two births to79, and first full renewal to4769|Actual graph and marked profile through14 vertices;25 arithmetic iterations check the inherited guard, whose induction is proved separately|
|audit|All1,436 selected target graphs and all2,618 parse trees|Independent forward Python path DP recounts every target; literal permutation/two-cover counts check both25 profiles and79|

The source bank is inherited from
[rooted_growth](../multicore_metagraph_20261113/README.md), with each old root
in the final vertex position. Its strong filter is equivalent to a positive
Hamilton-path endpoint count at every vertex: a strong tournament has a
Hamilton cycle, whereas the terminal SCC is the only possible endpoint SCC.
The inherited bank enumerates all ending-rooted classes; hence this filter
retains every rooted strong class. Counts by order are1,4,24,194,2395.
The checksum manifest pins the actual bank and the code used to interpret it.

The complete negative source statements below H45 use Busch's strong-order
minimum f(8)=45 and monotonicity at larger orders. They are not extrapolated
from the seven-vertex bank. They quantify over **strong** source graphs,
arbitrary same-H strong replacements, and arbitrary root choices. They do
not forbid the output integer through a different construction. For example,
count15 cannot exit to47 or49 by one such triangle, but our explicit witnesses
reach both from count9. The stronger inherited H35 unit-readiness obstruction
has its own all-tournament classification and is not inferred from this probe.

## Hostile examples and quotient losses

The seven-vertex count25 chain is the first failure of triangle-only **class**
completeness. It is the unique rejected class in the seven-vertex strong head.
Its count still occurs at order6, enabling the physical neutral replacement.

The inherited exact order9 spectrum in
[THM-4097](../../../01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md)
contains1,482 values. Attaching triangles to all strong order7 sources gives
1,481 values and misses89. Restricting those sources to the triangle-ear
grammar gives1,476 values and misses75,85,89,91,95,103. All six occur in our
grammar at smaller orders. This separates exact-order spectrum, graph-class
completeness and unrestricted value coverage. No enumeration speedup or new
finite coverage record follows from the ten-million-matrix census.

The arbitrary-atom grammar retains a complete graph and a physical root.
It is unbounded in core order and cross-matrix size. Its interface is not a
fixed-dimensional sufficient statistic. The descent preserves the root's
identity, not its start count, deletion count, cover count or readiness.
The neutral25 macro preserves the completed H and incoming port, not every
intermediate edit or the entire cover polynomial. Its uniform context theorem
is restricted to ordinal suffixes.

## Reproduction

From the repository root, with Python, NumPy and Numba available:

```text
python3 -B -X utf8 05-knowledge/results/ternary_core_grammar_20261113/reproduce.py
python3 -B -X utf8 05-knowledge/results/ternary_core_grammar_20261113/seal.py
```

The full command runs each of six producers normally and with Python -O,
then requires their deterministic compressed receipts to be byte-identical.
The long triangle census was run in both modes in this session; the remaining
five pairs were then replayed with `--remaining`. The normal/optimized check
is a guard-control test, not an independent algorithm. Independence comes
from the permutation, cycle-packing and Python-DP controls described above.
All decisions use explicit `need` exceptions rather than removable asserts.

The triangle engine's largest graph has9 vertices; its signed64-bit counts
are safely bounded by9!. The larger14-vertex profile check uses Python
integers. Software versions for these receipts were Python3.13.14,
NumPy2.3.2 and Numba0.64.0.

`manifest.json` contains content hashes for this packet and its directly
used repository inputs. `seal.py` verifies those hashes and local document
links. Published references keep their theorem ID **and filename**, since
the repository has legacy ID collisions.

## Repairs before promotion

The first atom recognizer tested only whether reversing one edge made a
graph transitive. A transitive triple also passes that test, by reversing an
adjacent comparison. The candidate atom must independently be **strong**.
The missing guard was caught by the count-descent audit, repaired, and both
normal/optimized runs repeated. No false atomic certificate was promoted.

The inherited THM-343 bibliography incorrectly paired Grinberg–Stanley with
arXiv2412.10572. Its two affected citations now distinguish their paper
arXiv2307.05569 from Irving–Omar's explicit Corollary20 in2412.10572. This
repairs attribution only; no Lean proof or axiom dependency was changed.
