# A typed reduction atlas: machines, tournament interfaces, and the limits of isomorphism

2026-10-26, Session 65. Source snapshot: `6a8ad5d3bf01a232da4de41bcbff68390a9337ee`.

**PROVED below:** elementary compilation, observational-equivalence, graph-shape,
and tournament-code statements. **CITED:** explicitly identified external
theorems. **RECOVERED:** the pinned repository continuation results.
**FINITE-EXACT:** the separate codec and hostile checks in `check.py`.
**PROPOSED:** the organization of the repository by these typed arrows.
No famous conjecture is solved, made undecidable, or proved complete here.

## 1. The useful version of the proposed classification

Turing machines are a useful common *execution language*. Tournament modules,
named ports, and continuation signatures are useful *interfaces*. Neither fact
by itself makes the underlying mathematical problems equivalent.

The strongest supported organizing object is a **typed directed atlas**. A node
must specify its input family, encoding, domain, predicate/output, and resource
measure. An arrow must specify a map, its direction, preserved predicate,
initialization, overhead, and the information it discards. Record mutual
arrows only when both have been proved in the same category.

The current genuinely shared mechanism is not an unexplained fractal. It is
**continuation sufficiency**: retain exactly what every permitted future
operation can still observe. Self-similar composition becomes meaningful only
after specifying a substitution grammar, admissible interfaces, base objects,
and the theorem transporting the target through substitution. A general
transition table does not automatically have a useful tournament-module
decomposition.

### Inheritance and concept board

* Closest proved mechanism: Session 33 continuation rows, native involution
  memory, and sharp tournament-module context budgets.
* Canonical hostile: equal current Hamiltonian count but unequal response after
  two external vertices; equal node output but unequal actual complement node.
* Corrected near miss: Session 64's genuine punctured 2-adic conjugacy does not
  preserve positive integers or the native Antihydra barrier.
* Least-used sidecar: the quantified *family of allowed continuations*, together
  with its legal states, outputs, guards, and initialization.

Live concepts: uniform language / quantifier alternation / tape configuration /
typed simulation / tournament interface / continuation congruence / realization.

## 2. Six different kinds of arrow

| Arrow | Required proof | What it does not say |
|---|---|---|
| Logical implication `P => Q` | A theorem in a stated background theory | No algorithm, same-rank claim, or cost bound follows automatically |
| Computable many-one `A <=m B` | One total computable `f` with `x in A iff f(x) in B` | Not polynomial-time; not an isomorphism |
| Polynomial many-one `A <=p B` | The same, with output size and time polynomial in input bit length | A fixed sentence is not an NP-complete language |
| Turing/oracle reduction `A <=T B` | An algorithm making oracle queries to B | May need many adaptive queries; weaker than a specified many-one reduction |
| Simulation/semiconjugacy | An encoding intertwining lawful steps, possibly after a proved finite macro clock | Generally one-way; spurious target states may remain |
| Computable conjugacy | A bijection with computable inverse, exact step/guard/output correspondence | May fail on the positive-integer/native subdomain; without cost bounds it says little about feasible computation |

Graph isomorphism is yet another object. Specify whether the graph consists of
finite control states, all tape configurations, one orbit, macro configurations,
or an observation quotient. Those are not interchangeable graphs.

For NP-completeness the comparison is between *uniform languages*, with an
NP-membership proof and reductions from every NP language, conventionally via
an established NP-complete language. Mere TM expressibility supplies neither.
NP-complete languages are themselves decidable: finite decidability alone
neither establishes nor rules out NP-hardness. It does rule out a computable
many-one reduction from the undecidable HALT language to that decidable target.
Cook's original paper provides the foundational resource-bounded viewpoint;
its original oracle-reducibility formulation should not be silently described
as the modern many-one definition. [Cook, 1971](https://www.cs.cmu.edu/~15455/resources/Cook1971-complx-thm-proof.pdf).

### Proposition 1: unrestricted computable reductions collapse all decidable nontrivial languages

If A is decidable and B has a known yes-instance b1 and no-instance b0, then
`f(x)=b1` when `x in A`, and `f(x)=b0` otherwise, is a computable many-one
reduction A to B. Consequently any two nontrivial decidable languages are
mutually computably reducible. The proof is the displayed algorithm.

This is why a useful atlas must retain resource bounds, structural conditions,
or constructive witness transport. Otherwise finite tournament connectivity,
bounded polynomial algebra, and enormous finite LRC checks collapse for an
uninteresting reason. This proposition gives no practical method to decide A.

A single unresolved sentence is even less suitable as a complexity-class
member: treated as a language with one possible input, it is a constant yes/no
problem in the extensional sense. Not knowing which constant is correct is
mathematical ignorance, not an undecidability theorem for that singleton.
One can instead compare proof strength over a fixed theory, or specify a
uniform family of related predicates.

## 3. Exact compilation and the quantifier boundary

Here predicates called decidable include their terminating algorithms. Merely
giving the code of a program that is *promised* total is a different promise
problem; testing totality is not part of these constructions.

### Proposition 2: counterexample compiler

Let R(e,n) be a uniformly decidable predicate, and put
`P(e) := forall n R(e,n)`. Effectively produce the blank-input program

```
for n = 0,1,2,...:
    if not R(e,n): halt
```

It halts iff `not P(e)`. Thus the false instances reduce to HALT and the true
instances reduce to NONHALT. This is an exact truth-preserving compilation,
not a simulation of the mathematical object's native dynamics. It need not
have a short runtime or preserve witnesses except for the enumerated n.

Sharing this compiler does **not** give `P <=m Q` for two particular families
P,Q: a common target is not a reverse arrow. The *whole class* of universal
decidable predicates has the usual Pi^0_1-complete universal nonhalting
problem; an individual natural subfamily need not be complete.
For a concrete control, the decidable even-number language and HALT both
reduce to HALT, but HALT cannot reduce computably to the even-number language.

### Proposition 3: totality uses an output process, not the same halting test

Let R(e,n,t) be decidable and
`P(e) := forall n exists t R(e,n,t)`. Search t=0,1,... for n=0; on success
emit one beep and advance to n=1, then repeat. Each successful search is
finite, hence:

```
P(e)     iff the program beeps infinitely often,
not P(e) iff the program beeps only finitely often.
```

If n is the least failed search, the program continues searching forever;
there is no finite signal that this n will never succeed. This is the precise
reason the naive full-Collatz counterexample detector does not halt.

For an arbitrary deterministic beeping machine,
`FINBEEP = exists T forall t>T [no beep at time t]` is Sigma^0_2, and
`INF BEEP = forall T exists t>T [beep at time t]` is Pi^0_2.
The preceding compiler applied to complements of arbitrary Sigma^0_2
predicates proves completeness of FINBEEP, using the standard effective
normal form. Aaronson's definition and classification are in section 5.10.
[The Busy Beaver Frontier](https://www.scottaaronson.com/papers/bb.pdf).

There cannot be a computable many-one reduction of all FINBEEP instances to
HALT. An elementary lower-bound witness suffices: simulate M silently; if M
halts, start beeping forever. This has finitely many beeps iff M does not halt.
A FINBEEP semidecision would semidecide NONHALT, contradicting the usual
diagonal halting theorem. This argument does **not** prove that one fixed
Collatz sentence has no alternative Pi^0_1 formulation. A natural Pi^0_2
presentation is an upper bound, not a strictness/completeness theorem for that
single sentence.

Ordinary halting and eventual cessation of a designated beep are therefore
different target predicates. A machine that runs silently forever is already
an elementary hostile to conflating them. Machine-state counts are also
model-dependent: symbols, halt conventions, output convention, and tape
initialization must remain part of every record.

## 4. Principal repository lanes, with exact quantifiers

The hierarchy descriptions below are *effective presentations/upper bounds*,
not claims of strict hierarchy membership or completeness of fixed statements.
`Delta^0_1` means decidable. The table is a principal-lane atlas, not a claim
to have classified every historical continuum problem in the repository.
In the table the ordinary positive Collatz map is `C(n)=n/2` for even n
and `C(n)=3n+1` for odd n. Section 6.5 separately uses its accelerated odd
section, with the valuation clock explicitly retained.

| Lane and exact object | Effective classification | Honest relation/status |
|---|---|---|
| An integer speed tuple `v=(v1,...,vk)`, with `M(v)=max_{0<=t<=1} min_i ||vi t||`, and question `M(v)>=1/(k+1)` | Decidable by finitely many rational interval/linear pieces | No generic NP-hardness or halting-completeness follows |
| LRC through r total runners, r variable | **Decidable by an explicit finite check**, using the cited inductive bound below | LRC(14) remains repository-open as a proof/computation target; finite decidability is not feasibility |
| All-runner LRC | `forall r [the finite check through r succeeds]`, hence Pi^0_1 | A counterexample-search TM exists; no reverse reduction from arbitrary NONHALT |
| The actual order-n tournament black-core graph | Finite, effectively constructible; connectivity of the precisely declared eligible core is decidable | Quantify over n for a Pi^0_1 all-order assertion; do not include blue-only isolated nodes silently |
| A000568(n), number of S_n-orbits of labelled n-tournaments | A total computable counting function, already with exact formulas/algorithms | It is not itself an unresolved yes/no conjecture. Efficient enumeration, exact native recurrences, or black-core connectivity are separate tasks; no #P-completeness asserted |
| Positive Collatz from one input n reaches 1 | `exists t C^t(n)=1`, Sigma^0_1 presentation | Simulation semidecides success; divergence has no generic finite certificate |
| Full positive Collatz | `forall n>=1 exists t C^t(n)=1`, Pi^0_2 presentation | Infinite-beeping compiler; falsity is finite-beeping/quasihalting in that compiler, not ordinary halting |
| Weak Collatz as used by Cryptids | No positive periodic orbit outside `1,4,2`; equivalently `forall n outside {1,2,4} forall t>=1 C^t(n)!=n` | Pi^0_1. Full Collatz implies weak Collatz; weak does not exclude divergent nonperiodic orbits |
| Antihydra blank orbit | Native guarded machine halts, or its Pi^0_1 negation | Session 64's unguarded 2-adic map is a different object; never drop the counter/guard or replace the marked seed by a generic measure claim |
| Strong Goldbach | `forall m>=2 exists prime p<=2m [2m-p prime]` | Inner search is bounded, so Pi^0_1; a counterexample is an even integer with its finite prime-pair search exhausted |
| Erdős / Sloane ternary digit conjectures | Universally quantified n with a finite base-3 computation of 2^n | Pi^0_1. The excluded small n and which digit is required are essential; these are different predicates |
| RH | A Pi^0_1 presentation obtained from a strict finite-precision violation of Lagarias's criterion, detailed below | This is a special effective theorem, not a general rule that arbitrary analytic root questions are decidable |
| `JC_n`: every complex polynomial map C^n->C^n with determinant 1 is invertible | Bounded-degree versions are decidable; `forall d JC(n,d)` gives a Pi^0_1 presentation | Stable JC/DC implications are genuine named theorem arrows, not same-rank isomorphisms |
| `DC_n`: every endomorphism of complex Weyl algebra A_n is an automorphism | Bounded-degree versions are decidable by the inverse-degree bound below; `forall d DC(n,d)` is Pi^0_1 | n is Weyl rank, hence 2n generators. Do not silently replace DC_2 by DC_1 or JC_2 |
| Finite polynomial systems over C or real closed fields, with rational coefficients | Decidable algebraic feasibility; inequalities require the real-field formulation | Changing to integer or rational points changes the problem |
| Integer-square/shared-root existence with polynomial equations and decidable nondegeneracy constraints | Sigma^0_1 existence; its no-solution assertion Pi^0_1 | Finite-field point counts, torifications, and root stacks do not decide integral-point existence |
| General integer-polynomial solvability, polynomial input unrestricted | Sigma^0_1-complete by MRDP | This genuine complete family does not make each fixed square-pattern problem complete |
| Rule 30 specified computable column is not eventually periodic | `forall N,p>=1 exists t>=N [a_t != a_(t+p)]`, Pi^0_2 presentation | No completeness follows; the frequency-1/2 statement instead has a Pi^0_3 presentation |
| ABC in its usual positive-integer formulation | `forall k>=1 exists K>=1 forall coprime a+b=c: c^k <= K^k rad(abc)^(k+1)`, Pi^0_3 presentation | This is a quantifier upper bound, not a resolution of ABC or IUT claims |

Mahler's strict real-tail problem, Euclidean Kakeya, and the global Hopf
problems need their own exact real/set coding and quantifier analysis. A finite
grid model or a polynomial chart is not an automatic equivalence with those
headline problems. They are deliberately left **UNCLASSIFIED HERE**, rather
than assigned ordinary halting merely because some computations touch them.
The repo's NC2/GMC(2) lane is recorded as closed in its current canon, not
silently listed as an open prize. Its algebraic and formalization claims are
not re-audited by this atlas.

### 4.0 Pinning the black-core target, not only its computational type

Let G_n have tournament isomorphism classes as vertices. For every tournament
T and directed Hamiltonian word P, add the line from [T] to [K_P(T)], where
K_P reverses precisely the arcs outside P's consecutive path arcs. The line
is blue iff reversal of P is an anti-automorphism of T; otherwise it is black.
Multiplicity does not affect connectivity. A pure-blue vertex is one for
which every Hamiltonian word is blue. The class-level target is that the
black-induced graph on non-pure-blue vertices is connected, for every n>=4.
All this data is computable by finite enumeration. A separate merged target
identifies each class with its converse before testing connectivity. The
class-level assertion implies the merged one, but a quotient path is not
automatically a lift; this atlas does not equate the two targets. This pins
the two conventions appearing in `HYP-9346-self-lines-connectivity-and-growth-laws-of-the-tiling-line-metagraph.md`
and Session 33's merged-node observer. Pure-blue vertices must not be included
in a supposedly connected black core.

### 4.1 The important LRC finite bound

Theorem A of Malikiosis–Santos–Schymura assumes LRC for n runners and proves
the next runner case for primitive positive integral velocities whenever the
sum of subset gcds exceeds `B(n)=binom(n+1,2)^(n-1)`. Since the singleton
terms include each velocity, a failure not inherited from a lower runner
case has all velocities at most B(n). Check successively n=1,...,r-1,
terminating immediately with false if a lower case fails. This decides the
statement “LRC holds through r runners.” No shifted-LRC claim is imported:
that extension has an additional Lonely Vector hypothesis.
[Theorem A, arXiv:2411.06903v2](https://arxiv.org/html/2411.06903#S1.SS1).

For an individual tuple the finite decision can be seen directly: for each
positive vi, its closed good-time set in [0,1] is the finite union
`[(j+1/(k+1))/vi, (j+1-1/(k+1))/vi]`, over j=0,...,vi-1.
Intersect these rational intervals. A nonempty intersection decides yes and
provides an exact rational witness; otherwise no. This naive algorithm is
pseudopolynomial in numerical speed sizes, not a claimed polynomial-time
algorithm in their binary encoding. Negative relative speeds can be replaced
by their absolute values, and repetitions can be retained.

One must not call this algorithm “enumerating to an unknown Busy Beaver
cutoff”: its finite bound is a mathematical theorem with an explicit formula.

### 4.2 Bounded JC/DC algebra is genuinely decidable

Fix n and degree d>=1. For JC write unknown coefficients for F, require
`det JF=1`, and write unknown coefficients for a candidate inverse G of degree
at most `d^(n-1)`. Equality of `F o G` and `G o F` with the identity is a
finite coefficient system. For DC use PBW-normal forms in A_n, require the
Weyl commutator relations, and seek an inverse of degree at most `d^(2n-1)`.
The inverse-degree bounds for actual automorphisms are the load-bearing
theorems, not an assumption that the conjectures are true.
[Bavula, Theorems 1.1–1.2](https://arxiv.org/pdf/math/0512215);
[Adjamagbo–van den Essen, Theorem 2](https://arxiv.org/html/math/0608009#S3).

Our deduction is then elementary effective algebra: construct the finite
first-order sentence over algebraically closed characteristic-zero fields
asserting that every coefficient tuple satisfying the input relations has
such inverse coefficients. Effective quantifier elimination decides that
sentence. Alternatively enumerate algebraic coefficient counterexamples and
test the finite inverse equations. Algebraic closure and the degree bound
prevent uncountable complex coefficients from being an unaddressed coding
gap. With unbounded d, there is no known global degree cutoff furnished here.

The established arrows are `JC_(2n) => DC_n => JC_n`, hence equivalence of
the *all-dimensional* conjunctions. This is not `JC_n iff DC_n` for each n,
not a reduction to a fixed-rank Poisson shadow, and not a polynomial-time
many-one theorem about the bounded coefficient languages.
[Belov-Kanel–Kontsevich, Theorem 1 and introduction](https://arxiv.org/html/math/0512171).

Status boundary: the repository cites its own exact audits of higher-dimensional
counterexamples and [Gao's 2026-07-31 preprint v1](https://arxiv.org/abs/2608.00222),
which explicitly uses the classical complex constant-Jacobian definition.
This session neither silently retracts that newer record from old background
knowledge nor independently validates its proof or claims consensus. Its live
JC(2)/DC(2) targets and the quantifier analysis above are separate from that
status question.

### 4.3 RH: avoid the equality-of-computable-reals trap

Lagarias proves RH equivalent to
`sigma(n) <= H_n + exp(H_n) log(H_n)` for every n>=1. His proof's reverse
direction shows that false RH yields a *strict* violation of this weak
inequality. Rational interval algorithms approximate the right side to
arbitrary certified accuracy. Dovetail pairs (n,precision), halting when an
upper interval endpoint is below the integer sigma(n). This halts iff RH is
false. It never needs a total decision procedure for equality of arbitrary
computable real numbers. Equivalently there is a decidable relation on finite
interval certificates and RH says no such violating certificate exists.
[Lagarias, Theorem 1.1 and section 3](https://arxiv.org/pdf/math/0008177).

### 4.4 Arithmetic roots: the base ring is part of the problem

The Session 59 shared-root variety is explicitly
`B_r=Proj k[c,a_j,b_j]/(a_j^2+b_j^2-2c^2)`. A split toric coordinate change
over a field containing i, rational point existence, a primitive integral
point with prescribed inequalities, and a finite-field point count are four
different consumers. A solution of the displayed equations with no
nondegeneracy conditions may even be trivial. The exact requested arithmetic
predicate must be supplied before assigning a machine to it.

The general comparison theorem MRDP says recursively enumerable predicates
are Diophantine, and yields a computable reduction from HALT to unrestricted
integer-polynomial solvability. That theorem supplies a *genuine* completeness
benchmark, not a shortcut to completeness of this fixed root surface.
[Matiyasevich, original 1970 paper record](https://m.mathnet.ru/php/archive.phtml?jrnid=dan&option_lang=eng&paperid=35274&wshow=paper).
The precise representation by an existential integer-polynomial equation is
also stated by Matiyasevich–Robinson in section 3, equation (3), of their
[author-posted primary paper](https://arxiv.org/pdf/0802.1052), which was
successfully read here; the older original-record full-page fetch timed out.

### 4.5 Cryptid source status is not proof status

The live Cryptids index distinguishes ordinary halting encodings from the
17-state full-Collatz **beeping** encoding. It currently lists Goldbach25 as
Lean-verified and a newer Goldbach21 as not independently verified; weak
Collatz14 has a linked verification, while Sloane(6,4) and full-Collatz17 have
explicit independent-verification caveats. These are source-reported labels,
not builds performed by this lane. The sibling machine lanes own the exact
tables, simulator/macro audits and any formal build evidence.
[Cryptids, checked 2026-10-26](https://wiki.bbchallenge.org/wiki/Cryptids).

## 5. What a Turing-machine graph actually is

A finite control table has states Q, a finite tape alphabet Gamma, and
instructions `delta(q,a)=(q',b,d)` or halt. It is a sparse directed labelled
multigraph after drawing an arc per instruction. Loops, two symbols with the
same target, opposite-direction arcs, and absent pairs are normal. The
labels include read symbol, write symbol and head movement, not just target.

A full configuration is `(q,h,tape)` with a head position and a finite-support
tape (or another explicitly declared tape model). Its one-step graph is
generally infinite. A deterministic machine has at most one successor per
configuration. The finite table does not contain the tape or its guards.

### Proposition 4: literal tournament obstruction

A finite tournament on N vertices has `N(N-1)/2` arcs. A simple directed
graph with outdegree at most k can be a tournament only if `N<=2k+1`.
Indeed `N(N-1)/2 <= Nk`. For a deterministic configuration graph k=1,
every induced tournament has at most three vertices; therefore no graph with
at least four configurations and outdegree at most one is itself a
tournament. The same argument applies to any four vertices of an infinite
configuration graph. Loops and parallel instruction labels are independent
additional obstructions. If a designated halt vertex has outdegree zero,
the improved bound `N<=2k` follows for N>1 and k outgoing arcs at every other
vertex. These are necessary conditions, never sufficient criteria.

### Proposition 5: one-run isomorphism is almost empty information

The visited-configuration graph of one deterministic computation, rooted at
its start, is a finite path ending in halt, a finite tail entering a cycle,
or a one-way infinite ray. Proof: until the first repeat there is a unique
new successor each time; once a repeat occurs determinism repeats the same
cycle forever. Hence all nonrepeating nonhalting runs have isomorphic
unlabelled rooted orbit graphs. If time instants rather than configurations
are used as vertices, *every* infinite run is a ray.

Conversely, a one-step halter and a two-step halter have the same truth value
but nonisomorphic rooted orbit paths. “Halts iff a counterexample exists”
is strictly weaker than graph isomorphism. Even finite-run speed changes or
idle padding alter the graph without changing the mathematics.

### Proposition 6: no symmetry-faithful ordinary tournament carrier for all control graphs

The odd-order automorphism fact is already used in repository
`THM-643-blue-black-line-parity-structure.md`, T1; the short proof is repeated
to expose exactly what the encoding obstruction needs, not as a novelty claim.

A nonidentity involution of a finite tournament cannot exist: it would swap
some pair u,v and hence send the unique arc u->v to the nonexistent arc v->u.
Thus every tournament automorphism group has odd order (Cauchy's theorem).
The same holds for vertex-colored tournaments because their automorphisms
form a subgroup. A finite labelled directed control graph consisting of two
identical loop states has an order-two automorphism. If a distinguished
start is demanded, add a fixed halting start and leave these two states
unreachable; their swap remains a control-graph symmetry.

Consequently there is no faithful functor on *all* such control graphs into
finite (even vertex-colored) tournaments: faithfulness would inject the
order-two automorphism into the target automorphism group. In particular no
fully faithful symmetry-preserving encoding exists in that category.

This does **not** rule out a graph-isomorphism many-one reduction, which need
not inject automorphism groups. Nor does it obstruct the subcategory of
ordered/symmetry-broken presentations. Rooted deterministic control graphs
in which every state is reachable by its labelled input words already have
trivial label-preserving automorphism groups; the displayed symmetry hostile
does not apply to that narrower category.

### Proposition 7: an honest ordered incidence-tournament codec

Let D be an ordered directed labelled multigraph with n states and m
instructions, allowing loops and repeated endpoints. Create disjoint classes
V, P, Q, of sizes n,m,m. State vertices retain start/halt/output decorations;
port p_e retains the complete instruction label. Define:

```
v -> p_e  iff v is the source of instruction e;
q_e -> v  iff v is its target;
p_e -> q_f iff e=f.
```

For every as-yet undefined pair (within each of V,P,Q), orient by the fixed
input order. Exactly one direction has now been assigned to each distinct
pair: the result is a vertex-typed tournament on `n+2m` vertices. The unique
P-to-Q matching, unique V-in-neighbor of p_e, and unique V-out-neighbor of q_e
recover every instruction, including loops, parallel arcs and its full label.
Thus decoding composed with encoding is the identity on the ordered input.
This proves injectivity and lossless presentation recovery.

Encoding/decoding takes polynomial time in the explicit presentation; the
dense tournament has `O((n+m)^2)` arc bits. For a k-symbol deterministic
table with at most kn instructions the vertex bound is `(2k+1)n`, counting
any halt vertex in n. A tape, head, initial state, guard convention and
observable are still required to execute it. Filler tournament arcs have no
computational meaning. One may decode once and run the original table; that
is a codec plus interpreter, not a computation moving freely along all
tournament arcs. No native black-core move or Hamiltonian-path predicate is
claimed to decide its halting.

The order/type decorations are essential to the *stated* claim. We do not
claim that arbitrarily relabelled source graphs have isomorphic encoded
ordinary tournaments, nor that an untyped tournament uniquely identifies
these port classes. Labels do not magically evade Proposition 6; this codec
deliberately works on ordered presentations rather than preserving every
original symmetry.

## 6. The genuine tournament-composition mechanism already in the repository

### 6.1 Exact modules: path-cover profiles, not a scalar or a universal machine

`THM-3121-path-cover-walk-content-substitution-kernel.md` proves, for a
tournament quotient Q and uniform modules S_i,

```
H(Q[S_i]) = sum_c [x^c] W_Q(x) product_i (c_i! pc_(S_i)(c_i)),
W_Q(x) = 1^T X (I-A_Q X)^(-1) 1.
```

Its proof splits each Hamiltonian path at maximal module runs, records the
quotient walk, and bijects the distinct factor-path components with its
occurrences. This is an exact composition theorem preserving a counting
predicate. It is not a statement about arbitrary tape computations.

Session 33 `tournament/BOUNDARY_OBSERVABILITY.md` sharpens this to an iff:
two factors have equal Hamiltonian counts in every uniform-module context
with at most m exterior vertices iff their ordered path-cover counts
`F(c)=c! pc(c)` agree through `c=floor(m/2)+1`. Each return to a uniform
module costs at least two exterior vertices. Triangular C3 contexts recover
the coordinates, with nonzero diagonal `(j!)^2`.

The first realized hostile has order five: both factors have H=15, but
F(2)=78 versus 90. A two-vertex context returns 123 versus 135. Thus equal
present output does not imply equal computational interface. Conversely,
full profile equality is enough for this *restricted* family of contexts;
named sparse wiring may still distinguish the factors.

### 6.2 Sparse ports require their actual boundary pairing

Session 33 `ports/PORT_CONNECTION_ALGEBRA.md` proves that distinct perfect
pairings of 2k named terminals are distinguished by some near-perfect
continuation. For an edge ab in P but not Q, order and orient Q's fragments
to form one path from a to b. The same continuation leaves P's edge ab
isolated from the remaining vertices. Equal degrees, component counts, or
electrical response are therefore insufficient.

For response matrix M(s,c), the coarsest deterministic observer quotient
consists of distinct rows; the least signed-linear observer dimension is
`rank M`. These are different minima, and rank over F_p preserves modular
counts rather than ordinary existence. Machine instruction interfaces need
their own guard/port output tensor, not an imported Hamiltonian tensor.

### 6.3 Native tournament dynamics: one operation is not every operation

Session 33 `metagraph-memory/NOTE.md` uses actual fixed-path tilings x,
complement of free arcs tau, and merged tournament node o(x). Since
`tau^2=id`, the least continuation state for the sole action tau is exactly
`(o(x),o(tau x))`; include color if color is observed. At order five there
are 38 deterministic/positive states but only 19 signed-linear dimensions.
The all-free-arc-flips extension instead distinguishes every raw tiling.
This is an operation-relative state lower bound, not a computational lower
bound against faster A000568 algorithms.

### 6.4 General continuation theorem, including guards

For states S, actions A with partial deterministic transition, and output o,
adjoin an explicit illegal sink with a distinct legality output. Define
`s ~ t` iff for every finite word w in A*, the entire observed future output
and legality agree when w is applied. Then ~ preserves current output and
is a right congruence: for every a, `s~t` implies `a(s)~a(t)`. Proof: testing
the future word w after a is exactly testing aw before it. Conversely any
output-preserving right congruence is contained in ~, by induction on word
length. Hence this is the unique coarsest exact deterministic observational
quotient. For a finite explicit system, partition refinement computes it.

In a noncommutative composition monoid, right contexts prove only a right
congruence; two-sided composition requires two-sided contexts. This is the
existing MISTAKE-249 boundary, not a new universal category theorem.

For an infinite TM configuration space this equivalence may have infinitely
many classes and need not be decidable. A finite observed prefix gives only
a finite-budget equivalence. No TM minimization algorithm follows from the
abstract characterization.

### 6.5 Session 64 supplies the decisive conjugacy hostile

Its `session64-reverse-valuation/REPORT.md` proves a valuation-word conjugacy
between accelerated Collatz and the twisted map H on correctly punctured
odd 2-adic domains. The inverse codes are

```
X_C(k) = -sum_(j>=0) 2^(A_j)/3^(j+1),
X_H(k) = (1/3) sum_(j>=0) (2/3)^(A_j).
```

Constant word 2 maps the positive Collatz fixed point 1 to the nonintegral
twisted fixed point 3/5. Constant word 1 maps Collatz -1 to twisted 1.
Moreover the native Antihydra odd-return process carries the guard/debt
coordinate `a_n=a_0+2 sum_(i<n) k_i-3n`, with all chronological guards
required. The unguarded 2-adic isometry preserves neither positive integer
realizability nor this native barrier.

This is an actual example where **a proved dynamical conjugacy is still not
an equivalence of the original open problems**. The same caution applies
even more strongly to a common tournament code or a shared valuation tree.

## 7. Finite controls and what they establish

Run from this directory:

```
python3 -B check.py > normal.json
python3 -B -O check.py > optimized.json
cmp normal.json optimized.json
```

The checker uses explicit runtime guards, no disabled Python assertions.
It exhausts all two-symbol partial deterministic target tables on 1–3 states
(4+81+4096=4181 tables), round-trips every ordered incidence codec, checks
loops/parallel labels, and rejects every single-arc corruption of a designated
canonical codec. It independently enumerates all labelled tournaments
through order five and verifies the involution obstruction. Two machines
whose source/read/target skeleton is identical but whose write bit differs
give a repeat versus a three-step halt; their exact simulations are hostile
controls for dropping tape-write labels. Exact rational evaluations check
the Session 64 constant-word domain mismatch.

These controls verify the finite constructions and examples. They are not
evidence for any all-input machine halt, hierarchy separation, LRC bound,
or other theorem beyond the written proofs and cited inputs.

## 8. What would count as real progress toward the proposed classification?

1. Choose one exact uniform source family and one *native* tournament target
   predicate—not “tournament structure” in general.
2. Produce one direction with code size, runtime/macro overhead, initialization,
   legal-state invariant, and a proof excluding spurious accepting behaviors.
3. Test positive and hostile instances independently at the raw configuration
   level. Preserve tape, guards, owners, labels and endpoint obligations.
4. Seek a reverse arrow in the *same* resource/semantic category. If it fails,
   identify the lost coordinate rather than filling the gap by analogy.
5. Only then state a complete-class theorem, an exact conjugacy, or a semantic
   interface equivalence. Keep these three verdicts separate.

The best immediate reusable conclusion is therefore a **reduction-and-interface
atlas**, not one giant isomorphism class of open problems. A finite machine
table can be packaged in tournament data, but the hard mathematical work is
proving which native predicate and which lawful continuations survive.

### Three carriers that must remain separate in the integrated session

* **Tape presentation:** an ordered join of transitive/cyclic three-vertex
  blocks can encode a finite binary tape word. Its elementary one-arc write
  operation is not a native K move; head and program state remain sidecars.
  An ordered finite join is a code, not intrinsically a fractal or an
  isomorphism of whole computation graphs. The sibling lane owns that theorem.
* **Arithmetic address geometry:** an iterated 3-adic/Cantor or related
  self-similar construction can preserve valuations and branch addresses.
  Integer realizability, orientation, native stopping guards and clocks must
  be independently preserved; Section 6.5 gives a concrete failure boundary.
* **Uniform decision predicates:** a compiler can preserve truth of a
  quantified family while changing its carrier completely. Its input size,
  quantifier level, target observer, converse and overhead are what matter;
  the number of vertices in a finite control table is not that classification.

## Source integrity

`SOURCE_PINS.json` records exact Git paths, source commit and byte hashes for
the inherited mathematical inputs copied under `inputs/`. The repository was
read only by this lane; no repository edits, commits or pushes were made.
The parent session owns integration and independent audit. External pages
were checked on 2026-10-26; no external formalization build is claimed here.
