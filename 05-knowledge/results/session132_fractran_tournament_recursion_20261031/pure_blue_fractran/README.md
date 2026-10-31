# A concrete 22-fraction stream for pure-blue merged tournament nodes

Session 132, 2026-10-31. **PROVED executable recurrence**, relative only to the
inherited all-order pure-blue classification for its tournament interpretation.
The fraction program and its sequence theorem have a self-contained proof below.
**FINITE-EXACT:** literal rational execution through the first 25 outputs;
independent exhaustive labelled-tournament census through order six.

## Result first

The ordered list in `pure_blue.fractran`, started at **817086270**, generates
checkpoints recognized by divisibility by **23**. At its nth checkpoint,

\[
v_{19}(N)=n,\qquad v_2(N)=p_n,
\]

where `p_n` counts pure-blue tournament isomorphism classes, equivalently
pure-blue converse-merged metagraph nodes. The empty object is included at
`n=0`; orders one and two are count conventions, not assertions about the
usual free-complement metagraph, whose geometric statements use `n>=3`.

The first 25 counts, indexed `n=0,...,24`, are

```
1,1,1,2,1,3,2,4,3,6,4,9,6,13,9,19,13,28,19,41,28,60,41,88,60.
```

The run reaches checkpoint 24 after **10,686 actual FRACTRAN steps**. No
macro, direct recurrence update, external input replacement, or interpreter
shortcut performs those transitions. The 22 fractions are small (largest
numerator 7973). This is not a shortest-program claim.

## 1. Inheritance and exact scope

Inherited snapshot: `fae72243028c39ce766cd9c59df763c9afd883e2` in
`eliottcassidy2000/math`, read from `/tmp/session131-source`.

The relevant proved canon is:

- `01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md`, current
  supersession and the SCC center law;
- `01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md`,
  current independent-completion section, **not** its incomplete historical
  odd-order argument;
- `05-knowledge/results/pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md`,
  the all-order prefix/violated-arc classification of pure-blue strong atoms;
- `05-knowledge/results/merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md`,
  section 3, the all-order Hamiltonian-bypass proof of the rigid strong alphabet;
- `05-knowledge/results/pure_blue_atom_classification_20261020/PURE_BLUE_COUNTS_AND_BRANCH_GEOMETRY.md`,
  the inherited counting law and self-converse gauge.

**Correction precedence matters.** The historical THM-4730 proof incorrectly
inferred nonpalindromicity from nonconstant reversal-partner label sums.
Current canon explicitly replaces its all-order consequence with an independent
prefix/violated-arc proof. The old `n<=125` cutoff is superseded; the stronger
auxiliary group lemma ML is not needed or claimed proved. Our actual-arc
checker retains the order-nine nonaffine anti-automorphism as a hostile test.
This packet does not repair or reactivate the rejected affine inference.

The inherited theorem says that the strong pure-blue alphabet is `K1,C3,T5`,
and the rigid strong alphabet is `K1,C3`. Here **T5 is the regular cyclic
five-vertex tournament** `x -> y iff y-x in {1,2} mod 5`; it is not the
nonregular exceptional tournament sometimes also called T(5), nor a Paley
construction over a five-element field.

A tournament is pure-blue when **every directed Hamiltonian vertex word**
has positional reversal an anti-automorphism of all arcs. This is intrinsic,
not a choice of coordinates. Under the fixed-Hamiltonian-path tiling model it
means every tiling in its fibre is grid-symmetric. Strong components are
uniquely linearly ordered, and the theorem gives the normal form

```
A -> Z -> A^op,
A = an ordered word of K1 and C3 components,
Z = empty, K1, C3, or T5.
```

`A^op` reverses component order and every arc. Since these atoms are
self-converse, it is represented by the reverse component word. Every
pure-blue tournament is self-converse, so merging a class with its converse
introduces **no extra factor of two**. These are unlabelled class counts,
not labelled-tournament counts, path counts, tiling counts, or all of A000568.

## 2. Recursive combinatorial specification and recurrence

Let `W1(T)=K1 -> T -> K1` and `W3(T)=C3 -> T -> C3`. The unique SCC normal
form gives the disjoint recursive specification

```
P = {empty, K1, C3, T5} + W1(P) + W3(P).
```

Why disjoint: the first four objects have no paired outer SCCs; every other
object has a unique outer SCC pair, either two singletons or two directed
triangles. Removing that pair gives the unique predecessor. The constructors
add respectively two and six vertices and preserve purity by the center law.
Hence, setting `p_j=0` for `j<0`,

\[
p_n=[n=0]+[n=1]+[n=3]+[n=5]+p_{n-2}+p_{n-6}.
\]

Consequently

\[
P(z)=\frac{1+z+z^3+z^5}{1-z^2-z^6},\qquad
(p_0,\ldots,p_5)=(1,1,1,2,1,3),\qquad
p_{n+6}=p_{n+4}+p_n\quad(n\ge0).
\]

The half-word count is

\[
c_m=\sum_{k=0}^{\lfloor m/3\rfloor}\binom{m-2k}{k},\quad
c_m=0\ (m<0),
\]

because a word with k weight-three letters has m-3k weight-one letters and
m-2k positions. Thus `c_0=c_1=c_2=1` and `c_m=c_(m-1)+c_(m-3)`.
This is [OEIS A000930](https://oeis.org/A000930), checked 2026-10-31. In
particular `p_(2m)=c_m`, `p_(2m+1)=c_m+c_(m-1)+c_(m-2)=c_(m+2)`.
The binomial formula and explicit normal-form generator provide checks
independent of the six-term stream recurrence.

The connection to tournaments is structural, not a tournament artificially
imposed on prime registers: vertices/arcs belong to the inherited tournaments;
the register machine transports **class counts under the two SCC wrappers**.
It loses metagraph adjacency, individual path counts and automorphism groups.
Recovering those needs the component word and center as sidecars; it does not
follow merely from the numerical sequence.

## 3. Complete ordered FRACTRAN program

```
527/58, 29/31, 37/29,
82/111, 37/41, 43/37,
141/215, 43/47, 53/43,
295/371, 53/59, 61/53,
7973/671, 61/67, 71/61,
803/923, 71/73, 79/71,
1079/1343, 79/83, 23/79,
551/23.
```

Use standard FRACTRAN semantics: scan left-to-right, select the **first**
fraction whose product with the current positive integer is integral, then
replace the integer by that product. Restart the scan. Fractions are reduced.
Do not reorder the work fraction behind its zero-case fraction.

**Input:** the single fixed seed

\[
N_0=23\,2^1 3^1 5^1 7^2 11^1 13^3=817086270.
\]

**Output:** observe, without modifying N, each state divisible by 23,
including the initial state. Its output is `v_2(N)` and its index is `v_19(N)`.
To obtain `p_n`, run until the checkpoint with index n. The program intentionally
streams forever; it does not halt after any specified output. Intermediate
states are not outputs, even if an exponent temporarily equals a sequence term.

The exact checkpoint invariant is

\[
N_n=23\,19^n\,2^{p_n}3^{p_{n+1}}5^{p_{n+2}}
7^{p_{n+3}}11^{p_{n+4}}13^{p_{n+5}}.
\]

There is no scratch factor 17 and no other control prime at a checkpoint.
The output convention is a **checkpoint valuation**, not the pure-power-of-two
output convention used by some other FRACTRAN programs.

## 4. Fraction-level invariant and termination before every next output

Registers are prime exponents: `(2,3,5,7,11,13)` hold `(a,b,c,d,e,f)`, 17
is scratch, and 19 is the index. The control states are checkpoint 23 and
seven main/auxiliary pairs

```
(29,31), (37,41), (43,47), (53,59), (61,67), (71,73), (79,83).
```

Exactly one control prime has exponent one at every reachable state, and all
others have exponent zero. Each denominator contains its own control prime;
no other state's fractions can apply. Fraction 22, `551/23=29*19/23`, enters
the first transfer and increments the index. For a transfer with source r,
destination product d, main q, auxiliary q', and following state s, its
ordered triple is

```
(q' d)/(q r), q/q', s/q.
```

If source exponent is positive, the first fraction decrements it, increments
each destination exponent, and moves to q'. The second returns to q. If the
source is zero, only the third can apply, entering s. No transfer has its
source among its destinations. These are actual fractions, not macro-level
execution rules added to the interpreter.

The seven transfers in order are

| Phase | Source | Destinations | Effect on values |
|---|---:|---:|---|
| 0 | 2 | 17 | put a in scratch; clear old output |
| 1 | 3 | 2 | b becomes the new output |
| 2 | 5 | 3 | c shifts left |
| 3 | 7 | 5 | d shifts left |
| 4 | 11 | 7 and 17 | e shifts left and adds to scratch |
| 5 | 13 | 11 | f shifts left |
| 6 | 17 | 13 | a+e becomes the new last register |

Thus the next checkpoint window is `(b,c,d,e,f,a+e)`, scratch is zero,
and the index has increased by one. Starting from the displayed seed and
using `p_(n+6)=p_(n+4)+p_n` proves the invariant inductively.

Every transfer terminates: its nonnegative source exponent decreases by one
per two FRACTRAN steps and is not replenished inside that phase. Seven such
phases and the initial entry therefore reach the next checkpoint in exactly

\[
8+2(2a+b+c+d+2e+f)
\]

steps. This gives termination **between every two outputs**, not termination
of the deliberately infinite stream. Each checkpoint has an applicable entry
fraction, and every intermediate state has either its work, return, or
zero-case fraction. There is no accidental halt. The proof also works for
arbitrary nonnegative initial register windows, not just the sequence seed.

## 5. Reproduction and independent controls

From this directory:

```sh
python3 verify.py
```

Only Python's standard library is used. The checker:

1. Compiles the seven explicit transfer triples, then asserts their equality
   with the **literal ordered list** loaded from `pure_blue.fractran`.
2. Executes that literal list by integer multiplication and divisibility,
   obtaining 25 checkpoint outputs. At each it checks the complete integer
   invariant, not just the output exponent. All 22 fractions are exercised.
3. Exhausts all 729 initial windows in `{0,1,2}^6`, including all-zero input,
   checking the next exact integer and exact step budget. Zero branches are
   therefore exercised directly rather than presumed from positive seeds.
4. Enumerates all labelled tournaments of each order `0,...,6`, independently
   enumerates their Hamiltonian paths, tests positional reversal on **all
   actual arcs**, and canonicalizes surviving classes over all permutations.
   It compares the full canonical-key set with the generated SCC normal forms
   and verifies that every surviving class is self-converse.
5. Tests positive atom controls and the hostile `T5 -> T5`, whose 225 paths
   include only 75 palindromic ones: `tc=9` and blue tiling count 3. Pairing
   pure atoms does not preserve purity when the outer atoms are nonrigid.
   Its Aut order 25 follows from the uniquely ordered SCCs; it is not a
   claimed all-labelled order-ten census.
6. Preserves the inherited order-nine nonaffine hostile: path
   `(0,1,2,3,7,5,6,4,8)` in `Circ(Z9,{1,3,4,7})` has partner sums `{5,8}`
   but its actual positional reversal **is** an anti-automorphism.

Direct census:

| n | Labelled universe | Labelled pure-blue | Pure-blue classes |
|---:|---:|---:|---:|
| 0 | 1 | 1 | 1 |
| 1 | 1 | 1 | 1 |
| 2 | 2 | 2 | 1 |
| 3 | 8 | 8 | 2 |
| 4 | 64 | 24 | 1 |
| 5 | 1024 | 184 | 3 |
| 6 | 32768 | 800 | 2 |

`verification.json` records all 25 six-register checkpoints, exact rational
step counts, every fraction's usage, census canonical codes and controls.
`first_transition_trace.json` records each literal integer and selected
fraction during the first 30-step transition. All files are below 100 KB.

**Novelty boundary:** the pure-blue alphabet, normal form, GF, and Narayana
connection are inherited mathematics, not rediscoveries claimed here. The
new deliverable is the fully explicit 22-fraction realization, its exact
checkpoint/step theorem, and independent regression bundle. This does not
compute all tournament classes A000568 or solve the mixed-node wiring.

## 6. PRIMEGAME-style orbit and valuation analysis

The frozen fraction list has **23 used primes** and **22 fractions**. Put one
valuation-change column per fraction into the `23 x 22` matrix M. The exact
standard-library calculation in `analyze_geometry.py` gives

```
rank_Q(M) = 22,
ker_Q(M) = {0},
dim ker_Q(M^T) = 1,
ker_Q(M^T) = span(total control exponent).
```

In fact deleting the row for checkpoint prime 23 leaves a square matrix of
**determinant +1**. `geometry.json` includes the full matrix, row-prime order,
conservation vector and exact determinant. This is a finite exact calculation
on the displayed fixed program, not an empirical claim about an infinite
family. Reproduce with `python3 -O analyze_geometry.py`.

The total exponent of the 15 control primes is invariant. On well-formed
inputs it is one, yielding the one-hot state manifold proved above. The
index exponent `v_19` is unchanged by every fraction except the entry fraction,
which increases it by one. Therefore no valid stream orbit can be periodic:
every round terminates by the exact budget above and strictly advances its
checkpoint clock.

The full-column-rank result supplies a stronger, separately justified
obstruction: a finite periodic orbit from **any positive integer** would
supply a nonzero nonnegative vector h of fraction-use counts with `Mh=0`.
There is no such vector. Thus no positive integer has a nonempty finite
periodic orbit under this list, regardless of unused inert prime factors.
This says **nothing further** about termination, checkpoint meaning, or
sequence correctness for arbitrary multi-control seeds. Those behaviors are
not classified here. A seed with no control factor simply halts; malformed
control multiplicities are not valid sequence inputs.

Two semantic graph projections clarify why graph cycles are not integer
periodicity:

- The 15-state control graph is one strongly connected component: the
  checkpoint/phase cycle has a reversible main/auxiliary pair at each phase.
  Repeated visits to a control vertex do not restore the integer, because
  the register exponents or clock have changed.
- The seven-register transfer graph is one strongly connected component.
  Its edges are `2->17`, `3->2`, `5->3`, `7->5`, `11->7`, `11->17`,
  `13->11`, `17->13`. The two destinations from 11 encode the recurrence's
  addition while preserving the shift value. The clock is separate.

These are typed directed graphs of machine operations, **not tournaments**.
They have missing pairs and, in the control graph, reciprocal edges; no
arbitrary tie-breaking is used to turn them into tournament objects. The
actual tournament connection remains the SCC-wrapper combinatorics in §2.

## 7. An exact integer clock for every individual FRACTRAN step

There is a stronger, explicit certificate than a nonzero determinant. Give
the used primes the following integer weights; all other primes have weight
zero:

| Data/index prime | 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Weight | 12 | 10 | 8 | 6 | 18 | 16 | 14 | 8 |

| Control prime | 23 | 29 | 31 | 37 | 41 | 43 | 47 | 53 | 59 | 61 | 67 | 71 | 73 | 79 | 83 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Weight | 0 | -7 | -8 | -6 | -7 | -5 | -6 | -4 | -5 | -3 | -4 | -2 | -3 | -1 | -2 |

For every positive integer N define the signed integer grading

\[
G(N)=\sum_p w_p v_p(N).
\]

**Microstep clock theorem.** For every fraction `a/b` in the displayed list,
`G(a)-G(b)=1`. Consequently whenever a FRACTRAN step is applicable,

\[
G(Na/b)=G(N)+1.
\]

Proof: prime valuations add on multiplication and subtract on exact division.
The 22 dot products `w^T M_j` are each one, checked directly in
`analyze_geometry.py` and recorded in `geometry.json`. For example the first
fraction gives `w_31+w_17-w_29-w_2=-8+14+7-12=1`, and entry fraction 22 gives
`w_29+w_19-w_23=-7+8-0=1`. This certificate does not depend on one-hot controls,
valid sequence seeds, fraction priority, or a termination assumption. The
negative control weights are deliberate; G is a signed clock, not a norm.
It is uniquely determined by `w^T M=(1,...,1)` after setting `w_23=0`, since
the only homogeneous freedom is the total-control conservation vector.

Thus on **any** actual trajectory, after t steps,

\[
t=G(N_t)-G(N_0).
\]

In particular no nonempty finite periodic orbit can occur, now with a short
explicit grading certificate rather than only the rank argument. The earlier
scope limitation is unchanged: arbitrary-seed halting and output behavior
remain unclassified. If a trajectory is infinite its integer values cannot
repeat, and therefore must be unbounded; monotonicity of N is not asserted.

For the valid stream `G(N_0)=108`. At checkpoint n the control weight and
scratch valuation are zero, so the **exact cumulative rational-step count** is

\[
\boxed{t_n=
12p_n+10p_{n+1}+8p_{n+2}+6p_{n+3}
+18p_{n+4}+16p_{n+5}+8n-108.}
\]

At `n=24`, with window `(60,129,88,189,129,277)`, this is **10,686**.
Subtracting the formulas at consecutive n and substituting
`p_(n+6)=p_(n+4)+p_n` recovers the seven-loop budget
`8+2(2p_n+p_(n+1)+p_(n+2)+p_(n+3)+2p_(n+4)+p_(n+5))`.
The checker verifies the formula against all 25 actual checkpoint runtimes
and G against each of the 31 recorded integers in the first transition.
This is an exact all-order time/readout certificate, not a fitted runtime law.

The fraction list is unchanged by this extension. `SHA256SUMS` freezes the
portable source, literal program, proof and verification artifacts (excluding
the manifest itself and interpreter bytecode caches).
