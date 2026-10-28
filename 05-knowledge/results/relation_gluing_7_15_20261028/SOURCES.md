# Sources, scope and reproduction

Initial clean worktree commit:
`83b601158076669edc8729ce461f150da1c3877e`.
Incoming Session85 integration was read from
`dd908f4f88e7583a4391b6f1ad208f777848784c` before rebasing.
The corrected THM-4839 and MISTAKE-693 were read from its ancestor
`919e32c58b`. Later incoming work through1b86e4eadc is recorded in INCOMING.md.
Final source content hashes and full commit pins are in
[source_pins.json](source_pins.json). Source claims retain their own scope;
reading a packet is not a new independent audit of its dependencies.

## Primary literature consulted

- Alex Fink and Luca Moci, [Matroids over a ring](https://arxiv.org/pdf/1209.6571),
  definition 2.1 and equation 2.1: subset quotient modules for an integer
  configuration. Used to name the classical framework; all numerical
  Smith calculations here are explicit.
- G. David Forney Jr., [Codes on Graphs: Fundamentals](https://arxiv.org/pdf/1306.6264),
  sections 5.5 and 6.3: exact cycle-free code realizations and their state
  spaces. The Fano nongraphic obstruction excludes ordinary cycle-matroid
  realizations, not graphs decorated with arbitrary parity constraints.
- John Cremona and Sage contributors, [number-field Selmer groups](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/selmer_group.html):
  the unit/class-group exact sequence. Consulted as primary implementation
  documentation; Sage was not executed in this session.
- David Dummit and John Voight, [The 2-Selmer group of a number field and
  heuristics for narrow class groups and signature ranks of units](https://arxiv.org/abs/1702.00092),
  section 3 and Theorem 6.1. The elementary sequence is proved in T7;
  the Hilbert-symbol isotropy theorem is a cited direction only.

These four sources were checked on 2026-10-28. No external source establishes
literature priority for the packet's explicit bank or finite counts.

## Recovered repository material

The [atlas](ATLAS.md) gives theorem IDs with exact paths, the map, preserved
predicate, lost information, necessary extra data and cheapest hostile
test for each transferred claim. The source manifest records the files
actually used, rather than treating the session log as current truth.
Important inherited corrections include the 14-versus7 Fano circuit count,
the actual2/4 collision, the norm/prime guard in the Berggren semiprime
map, THM-4839's corrected statement, and the class-grid audit MISTAKE-690.

Incoming Session85 contributes independent row-space LCA accounting and
finite Boolean recovery certificates. The present packet credits those
results and adds the column-kernel ring extension, exact projective tree
width/count, integral-module symmetry census and unrestricted sigma fusion.

The random recovery used Python's Random(715) on the sorted eligible paths
in the initial worktree. Eligibility was a top-level results Markdown file
containing `Fano|seven|fifteen|\b7\b|\b15\b`. The candidate count1772,
six chosen paths, and their treatment are recorded in
[random_recovery.json](random_recovery.json) and the atlas. This discovery
sample supplies no probabilistic claim about mathematical importance.

## Exact computational universes and controls

| Calculation | Complete declared universe | Independent comparison or hostile |
|---|---|---|
|Arithmetic code7/code15|128/32768 prime subsets|Integer square roots versus parity; character-product weight polynomial|
|Minimal dependent supports|All zero-sum subsets of size<=r+1, r=3,4|Proper-subset tests; rank circuit-size bound|
|Gluing|All4096 binary3-by-4 matrices, all 16 cuts|Actual kernel fibres versus image intersections; zero/repeated columns included|
|Fano tree width|All10395 rooted unordered labelled binary trees|Closed combinatorial count693; every cut rank checked|
|Projective width15|All32768 cuts and an explicit decomposition|All balanced cuts satisfy the lower bound; proof covers every r>=2|
|Coefficient change|729 integer2-by-3 matrices with entries{-1,0,1}, characteristics 2,3,5|Smith factors versus modular elimination; signless triangle|
|Arithmetic integral types|All127 nonempty subsets of the seven columns|Separate gcd-of-minors implementation, including the empty subset|
|Symmetry|All168 GL(3,2) permutations, all 128 subset types|Separate enumeration of all 5040 point permutations; explicit integral involution|
|Projective geometry|All35 lines, all 56 spreads, all 6 isotropic spreads|K6 perfect matchings/edge factorizations; explicit4-colour and 7-spread certificates|
|Relations become points|All256 old group-law pairs and 455 new triples|Direct divisor-list sums in the separate audit|
|Fusion|6400 products1<=n,m<=80;864 local cocycles|All shared primes retained; explicit formal/actual520 hostile|
|Prime-power presentation|1296 exponent choices0..5 on{2,3,5,7}|All inherited FUS shapes and off-shape controls|
|Number-field examples|d=5,13,37; two unit quotients|Displayed ideal identity and complete norm2 bounds; general statements proved, not inferred from a search|

Prime-bank construction tries positive a in increasing order for
p=d*a²-1 and stops at the first prime. All displayed primes have a<=9;
primality is tested through the integer square root. This verifies the
least representative in each stated square class, not prime infinitude.

## Reproduction

Python3.13 and SymPy1.14 were used. Only check.py uses SymPy; audit.py and
verify.py use the standard library. From this packet directory run:

```text
python3 -B verify.py --replay
```

The replay runs check.py and audit.py in ordinary and optimized modes,
compares their bytes with the saved receipts, checks the manifest and
local links, and fails on an explicit exception. Checks use ordinary
conditionals, not Python assertions that vanish under -O. The producer
and separate audit are different implementations by the same agent;
this is not an independent human or multi-agent referee report.

The saved outputs are normal.json, optimized.json, audit-normal.json and
audit-optimized.json. Their byte equality and PASS status are recorded
by verification.json. The projective-width and module theorems have
all-size proofs; bounded receipts are supplementary implementation controls.

[make_diagram.py](make_diagram.py) uses Matplotlib to draw the overview
from the saved receipt. The figure is explanatory, not evidence for a
claim omitted from the exact calculations. The PNG hash uses raw bytes;
text hashes use LF-normalized bytes. The packet's .gitattributes fixes
text line endings for reproducibility across Windows checkouts.

## Boundary of this research claim

No global theorem IDs were reserved. The packet is a self-contained,
routed research result with scoped elementary proofs and explicit finite
certificates. It claims no solution of LRC(14), graceful trees, square-sum
Hamiltonicity, Collatz, Rule30, the planar Jacobian conjecture or BSD.
The cited supplied BSD manuscript is not a proved dependency here.
Neither every binary configuration's prime realization nor an efficient
universal arithmetic representation is asserted.
