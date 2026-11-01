# Square-sum graphs: a 14-fraction predicate and literal finite Turing certifiers

**Session135. PROVED:** elementary machine semantics and uniform termination,
conditional only on the stated input conventions. **FINITE-EXACT:** the ranges
in the reproducible receipts. **CITED:** Gerbicz's all-order existence theorem.
No minimality of machine size, new Hamiltonicity threshold, or new general
square-testing algorithm is claimed.

## 1. Inheritance before computation

The objects are the simple graphs `G_n` on the actual labels `1,...,n`, with
an edge `uv` exactly when `u != v` and `u+v` is a square. A cycle has at least
three vertices. A path on one vertex is permitted; empty graphs are excluded
from these machines' accepted language.

The closest prior mechanism is Session131's reflection-root compiler:
`R_k(v)=k²-v`, with an exact alternating affine walk, interval restrictions,
and collision sidecar. It already supplies the even/odd-length closure split.
That is recovered structure, not a new reflection theorem here.

Session124 recovered Gerbicz's parity-compatible two-path invariant and its
49-fold affine lift. A finite table of universally checked transitions,
complete base interval, and a decreasing quotient give a genuine constructor
for every order. The separate protected-prefix constructor preserves a
particular infinite ray; it does not prescribe an arbitrary path. Neither is
replaced by our witness checker.

Session125's rooted-forest bijection concerns an **ordered path**, not arbitrary
square-sum graphs. Its leaf-at-end/leaf-at-middle hostile gives forest counts
21 versus20 from the same old count8. Its Green-matrix sidecar answers forest
queries, not Hamiltonian continuation. Session131's boundary matching retains
actual labels, endpoint incidences, path pairing, and connectedness; cutting a
selection of local cycles is not a complete parametrization of global cycles.
Our certifier retains the full named-vertex sidecar rather than a determinant,
root histogram, or unlabeled port count.

Primary attribution was checked on2026-11-01: Robert Gerbicz's2018-01-21 author
comment in [OEIS A090461](https://oeis.org/A090461) gives paths for every `n>=25`
and cycles for every `n>=32`. Session124's independently replayed construction
is the repository dependency for the actual proof mechanism. This packet does
not redistribute Gerbicz's original C source or require it to run.

The live concept board is:

| Representation | Preserved question | Needed sidecar / hostile |
|---|---|---|
| Reflection channel `k²-v` | Exact adjacency and graph arrival | Bounds and distinct labels |
| Odd-subtraction counters | Is an integer sum a square? | Control state and fraction order |
| Unary vertex witness | Hamiltonian path/cycle verification | Named U/X marks and closing edge |
| Root histogram / computation time | Additive edge-root weight | Forgets cyclic order and collisions |
| Gerbicz paired lift | Uniform constructor | Shared positional parity and complete bases |

## 2. Graph inclusion reduces to a square shell

The new neighbors of vertex `n` are exactly

```
{k²-n : floor(sqrt(n)) < k <= floor(sqrt(2n-1))}.
```

Indeed `1<=k²-n<=n-1` is precisely `n+1<=k²<=2n-1`. Distinct roots
produce distinct neighbors, so the arrival degree is

```
d_new(n)=floor(sqrt(2n-1))-floor(sqrt(n)),  n>=1.
```

This is an exact generator without testing all old vertices. It is an elementary
restatement of the inherited reflection channels, not evidence that graph
arrival transports a particular Hamiltonian cycle. Session125's32-to33 edit
hostile already separates these operations.

For an actual edge query `u,v`, use the square recognizer below on `s=u+v`,
with separate preconditions `1<=u,v<=n` and `u!=v`. The recognizer alone does
not enforce these graph conditions.

## 3. The literal fourteen-fraction program

The exact ordered list in `square.fractions.txt` is

```
43/14, 31/7, 66/43, 65/33, 17/11, 19/26, 37/13,
11/19, 23/34, 31/17, 58/23, 123/145, 99/29, 29/41.
```

For every integer `s>=0`, start at `7*2^s`. Data registers are exponents
`N=v_2`, `O=v_3`, and `A=v_5`. The other primes are one-hot controls.
The halting controls are31 (square) and37 (nonsquare). Then:

* `s=0`: final integer31, in exactly one fraction step.
* `s=k²>0`: final integer `31*5^(2k-1)`, in exactly `5s+2` steps.
* `s=k²+r`, `1<=r<=2k`: final integer
  `37*3^(2k-r)*5^(r+1)`, in exactly `3s+2k²+4k+4` steps.

Thus the outcome is read from a control prime, not from whether the terminal
integer itself is a square. Positive accepting runs also retain the odd number
`2sqrt(s)-1` as a prime exponent. Zero has its separate clean branch.

### First-eligible semantics, without a compiler assumption

Each reduced fraction denominator has exactly one source control, and every
numerator exactly one different target control. No control factor cancels;
there is no self-loop fraction. With a one-hot input, fractions belonging to
other states cannot apply. The data-dependent rule always precedes its fallback.
The state actions are:

| Control | Nonzero operation / next | Zero operation / next |
|---:|---|---|
|7| decrement N, go43 | accept31 |
|43| increment N and O, go11 | — |
|11| decrement O, increment A, go13 | go17 |
|13| decrement N, go19 | reject37 |
|19| go11 | — |
|17| decrement N, go23 | accept31 |
|23| increment N, go29 | — |
|29| decrement A, increment O, go41 | increment O twice, go11 |
|41| go29 | — |

The temporary decrement at17 tests positivity without illegal FRACTRAN
cancellation of a guard factor. Its restoration at23 is load-bearing.
At the start of the `j`-th subtraction, `O=2j-1`, `A=0`, and
`N=s-(j-1)²`. Consuming an odd number moves it into A while removing the
same number from N. If N is zero after the full odd is consumed, accept.
If positive, restore O from A and add two. If N becomes empty before the
odd is exhausted, reject. This proves the exact terminal forms. Every loop
consumes a finite counter; the cumulative subtractions `1+3+...+(2j-1)=j²`
force a terminal branch for every permitted input.

A completed nonterminal subtraction of odd `d` costs `5d+4` steps. A terminal
exact subtraction costs `3d+2`. The positive-input setup costs two steps.
Summing these gives `5k²+2` for squares. For nonsquares, k completed stages
cost `2+5k²+4k`, followed by `3r+2` steps to failure. These are the formulas
above. Execution is linear in the unary exponent s, not a claimed fast binary
algorithm. Fraction order is part of the program, not a cosmetic listing.

## 4. Actual finite six-tape Turing transition artifacts

There are two machines, not a general-purpose virtual machine:

| Machine | Nonhalting states | Disjoint symbolic rows | Fully expanded delta rows |
|---|---:|---:|---:|
| Path |25|59|50,625|
| Cycle |36|84|72,900|

`path.tm.json` and `cycle.tm.json` contain finite disjoint wildcard abbreviations.
The builder-generated `path.tm.tsv.gz` and `cycle.tm.tsv.gz` expand **every** state/read-symbol tuple,
including all rejecting defaults. Their rows are ordinary transitions
`state, six reads, next state, six writes, six L/R/S moves`.
The tape alphabets have sizes5,5,3,3,3,3, so each nonhalting state has2025
explicit entries. Accept and reject are halting states and have no successors.
`build_tm.py` verifies the abbreviations are disjoint before expansion.
The compact publication retains the JSON clauses, deterministic expander, and
frozen expanded SHA256/counts in `TM_BUILD.json`; binary gzip files are derived
locally, not publication dependencies. `check.py` materializes them if absent
and checks every row, unique key, total count, and frozen expanded hash. A row
with wildcards abbreviates a family of **single ordinary TM steps**, not a
loop or an arithmetic macro.

The input tape starts at cell0; all work tapes are blank with heads at0.
The input is the finite nonblank word

```
^ 1^n # 1^v1 # ... # 1^vm $
```

without spaces. The `^` and `$` are literal sentinels, not exponent notation.
Every numeric field is nonempty; the only other symbol is `#`. The cell immediately
after the terminal `$` is blank; as usual, the finite input word contains no
blank symbols. Extra nonblank data after `$` is rejected.

The exact accepted language is:

* **Path:** `n>=1`, `(v1,...,vm)` is a permutation of `1,...,n`, and every
  adjacent sum `v_i+v_(i+1)` is a square.
* **Cycle:** the same, `n>=3`, and `v_m+v_1` is also a square.

The simulator performs only lookup, symbol writing, and one-cell head movement.
It does not call `isqrt`, sort the input, query a graph, or execute Python
arithmetic as an uncounted transition. `isqrt` appears only in the separate
test oracle and in the human-readable theorem.

### Tape invariants and proof

Tape1 is `^ U...U $`, with one cell for each actual label1 through n.
While consuming a vertex's unary field, its head advances by one named cell
per input1. Reaching the terminal marker rejects an out-of-range label.
At the field delimiter, `U` is changed to `X`; an existing `X` rejects a
repeated vertex. A final scan accepts only if every cell is X. Thus no separate
vertex-count counter is needed: injectivity plus complete marked coverage
forces exactly a permutation of the interval.

Tape2 holds the last vertex in unary; tape3 holds the first. Before reading
another vertex, the old Tape2 content is copied destructively to Tape4.
The new vertex is written to Tape2 and appended to Tape4. Therefore Tape4
contains exactly the adjacent sum, without constructing the graph.

Tape5 holds an odd unary subtrahend. One Tape4 symbol is erased per Tape5
symbol, then the odd is increased by two. Exhaustion exactly at an odd-block
boundary succeeds; premature exhaustion rejects. The preserved left sentinel
makes both conditions literal symbol tests. Successful cleanup empties the
sum and odd work tapes. The cycle machine copies the last and first vertices
into Tape4 and invokes the same test once more before its final coverage scan.
The path machine omits this closing test.

All unrecognized symbols or state combinations go directly to REJECT. On
permitted blank-work-tape starts, every rewind has a previously written left
sentinel. Parsing consumes finite input; copying/rewinding crosses finite
stored words; square subtraction consumes the sum. Hence both machines halt
on **every finite input word** over the nonblank input alphabet, not merely
on successful certificates. The accepted-language equivalences follow from
the invariants in both directions.

### Complexity and exact structural clock

For unary sum `s=k²`, the square subroutine costs exactly `2s+6k-2` elementary
transitions. A later vertex v following p costs
`6p+4v+6sqrt(p+v)+6`, including the separator, marks, copying, and edge test.
Consequently for an accepted path define

```
R_path = sum_(i=1)^(n-1) sqrt(v_i+v_(i+1)).
T_path = 5n²+14n+3-2v_1-5v_n+6R_path.
```

For an accepted cycle, including its closing root,

```
R_cycle = sum_(i=1)^n sqrt(v_i+v_(i mod n+1)),
T_cycle = 5n²+14n+11+3v_1+6R_cycle.
```

The root weight is invariant under rotation and reversal. The `3v_1` term is
the machine's initial-coordinate overhead, not a cycle invariant. Normalizing
`v_1=1` removes this gauge and gives `5n²+14n+14+6R_cycle`. The inherited
canonical32-cycle has root weight180 and takes6,662 transitions. The supplied
15-path has root weight55 and takes1,607.

For every Hamiltonian cycle its root multiset satisfies

```
sum k_i² = n(n+1),    sum k_i = R_cycle is even.
```

The first identity counts each vertex twice; the parity statement follows
from `k² == k (mod2)`. Thus normalized cycle runtime is congruent to
`5n²+14n+14 (mod12)`. This is a necessary clock invariant, not a certificate
of existence or connectivity. By Cauchy–Schwarz,
`R_cycle <= n sqrt(n+1)`; in particular the runtime is `5n²+O(n^(3/2))`.

Let L be the actual length of an arbitrary finite input word. On accepted
permutation inputs specifically, `L=n(n+1)/2+2n+2`. All-input runs take O(L+1)
transitions, including rejection of the empty input, with O(n) work cells
once a well-formed positive n-header is read. Each
copy or rewind is charged to the unary field that created it; each square
pass costs O(p+v). Malformed fields stop rather than extending unbounded
work. For accepted inputs this is O(n²) time and O(n) work space. This is
**linear in the expanded unary input**, not linear in a binary witness or a
claim that Hamiltonian search is easy. No adjacency matrix is stored.

## 5. The strongest hostile and the exact boundary

A moment-matched version of Session131's repeated-reflection hostile is:
at n=8 the word `(1,8,1,8,1,8,1,8)` has eight entries, all within the interval,
label sum36, and every cyclic neighboring sum9. Its root histogram satisfies
`sum k²=72=8*9`, with even root weight24. It passes every local square test
and both scalar moment identities, yet is not a Hamiltonian cycle. Our machine
rejects at the repeated named vertex1. Even length plus total label sum plus
all local square edges cannot replace the U/X sidecar.

Conversely the supplied15-path is a genuine permutation with valid internal
edges, but endpoints8 and9 do not close. The path machine accepts and the
cycle machine rejects. The closing-edge obligation cannot be silently dropped.

The connection ledger is therefore precise:

* **Source:** a named unary vertex word and n.
* **Map:** finite Turing verification, with adjacent sums reduced to odd
  subtraction; the same local square predicate has an actual14-fraction form.
* **Preserved predicate:** being exactly a Hamiltonian path/cycle of G_n.
* **Discarded information at final yes/no:** the witness ordering, labels,
  root list, and alternative solutions. Runtime alone retains only root weight
  and the explicit starting/endpoint terms.
* **Sidecar:** named visited marks, first/last vertex, and cycle-mode flag.
* **Cheapest hostile:** `(1,8)^4` at n=8 and path15 with nonsquare closure.

A verifier gives finite-order decision by exhaustive permutation enumeration,
but this packet emits no such enumerator and makes no search-complexity claim.
For fixed n enumeration is finite; for all n, accepted witnesses are not an
induction. Uniform existence remains the cited Gerbicz construction with its
complete bases and preserved paired invariant. No FRACTRAN universality
appeal is needed anywhere in the explicit machine proofs above.

## 6. Reproduction and provenance

With Python3 standard library, from any directory:

```
python3 /path/to/square_machines/build_tm.py
python3 /path/to/square_machines/check.py
python3 -O /path/to/square_machines/check.py
python3 /path/to/square_machines/simulate_tm.py path 15 8 1 15 10 6 3 13 12 4 5 11 14 2 7 9
```

The builder emits full delta tables and their hashes. The producer checks
compare every concrete row with the finite source specification, every stated
terminal exponent shape, bounded malformed words, inherited positive witnesses,
cycle gauge changes, and hostile mutations. The exact ranges are in the JSON
receipts; no test is hidden behind removable `assert`. The independent audit
uses separate interpreters and is distinct from normal/optimized equality.

`SOURCE_PINS.json` records the inherited commit, exact paths, Git blobs, and
SHA256 values. The32/33-cycle words are copied from Session125's independently
regenerated origin census;80..87 are the already published Session131 witnesses.
They are controls, not new cycle discoveries. The15-path is directly checked.
