# Injective binary observer assembly and its sharp coupled blocks

23 October 2026. **PROVED:** exact binary-record assembly, a constructive
sharp coupled-move theorem, and native multisource families with arbitrarily
large necessary simultaneous support. **FINITE-EXACT:** the declared
Boolean controls and original tournament seeds. **OPEN:** unrestricted
native wiring and general packet families with more than two full states.

## 1. Inheritance, changed question, and the actual interface

The closest proved mechanism is the
[Session44 shared-port and one-hot compiler](../../session44_guarded_native_moves_20261023/coupled-moves/THEOREM.md).
It retains local kernels, actual owner equalities, literal port capacities,
and every supplied guard; bounded conformal kernel moves connect its
integer fibre. Its fixed-owner specialization uses alternating matchings.
The present result identifies a smaller interface with explicit component
blocks and an exact optimal support budget. General implication-graph
logic is standard; no priority claim is made for Boolean satisfiability.

The canonical hostile is separate feasible source matchings that require
incompatible values of the same actual owner. The corrected near miss is
equating phase equality or an aggregate black capacity with available
distinct ports. The least-used sidecar is the **full source-labelled port
inventory of each local state**, including converse partners when required.

The concept board is: binary full records / shared actual coordinates /
distinct resources / implication blocks / endpoint-fixed padding. The
anchor is an effective assembly theorem, the niche is sharp move cost,
and the wildcard is a native family in which two changing shared variables
force arbitrarily many source slots to move together.

For each packet `i`, supply one or two authenticated full states `s_i(0)`
and, if present, `s_i(1)`. Packets own disjoint labelled assignment slots.
Their states contain:

- every actual source, path, complement partner and required converse copy;
- all local guard and pin certificates;
- the value `o_(i,v)(b)` at every occurrence of each shared owner coordinate;
- a set `R_i(b)` of distinct literal port/resource identities used by that
  packet state.

The set of coordinate occurrences is fixed before choosing a state.
Each state is internally legal and injective. Two states of a genuine
binary packet are distinct full records, not two names for one quotient
record. Every retained kernel coordinate is included in those records;
collapsing a larger fibre to two quotient names does not meet this input
contract. Singleton packets may be encoded by forcing their bit to zero.

Supply the already occupied literal inventory `U` of the unchanged global
remainder. Global compatibility requires equality of the actual values at
all occurrences of each shared coordinate, pairwise disjoint `R_i`, and
disjointness from `U`. Optional fixed state pins are allowed. All other
guards must either hold in every local state with their required per-slot
meaning, or be included as an additional authenticated binary constraint.
An arbitrary higher-arity electrical or owner predicate is outside this
interface unless its reduction is separately proved.

## 2. Exact assembly is a two-clause Boolean problem

For each forbidden joint choice `(b_i=a,b_j=b)`, add the clause

```text
not(b_i=a) OR not(b_j=b).
```

A pair is forbidden if its shared actual owner values disagree or its two
resource sets overlap. Forbid a single state when it meets `U`, violates a
fixed state pin, or is the absent state of a singleton packet. Comparing
each shared occurrence to one nominated occurrence suffices; comparing
all pairs is also exact. Constants and partial intersections of owner
alphabets need no special approximation: test their actual four state
pairs. Equality/complement relations are a useful special case, not an
extra hypothesis.

**Theorem 1.** The satisfying Boolean assignments are in bijection with
the declared compatible, injective full observer assemblies.

**Proof.** A bit vector selects exactly one complete local state. Each
forbidden-pair clause enforces exactly its literal two-record condition.
Pairwise owner equality is global equality; pairwise disjoint finite
resource sets are a global injection. The unary clauses handle the fixed
remainder and pins. Conversely every legal assembly makes each such clause
true. No fractional matching, owner marginal or resource count replaces
the literal records. QED.

This is an explicit finite algorithm polynomial in the supplied record
and clause lists. It does not bound the cost of obtaining those complete
native records. Resource identities can be indexed so that clauses are
generated from actual collisions rather than every pair of unrelated
packets.

### A minimal conditional-matching failure

Let two sources have distinct resource alphabets. Source A has one slot
choosing `p` at shared owner value zero or `q` at value one, plus a fixed
slot using `p`. Source B has one slot choosing `r` at owner value zero or
`s` at value one, plus a fixed slot using `s`.

Each source separately has a complete injective assignment. A forces the
common owner to one; B forces it to zero. Every per-source matching/Hall
test succeeds, but the assembled packet does not exist. The compiler
retains both unary consequences on the same actual owner coordinate and
rejects it. This is an abstract hostile for the interface, not a claimed
realization of these four lists inside a tournament.

## 3. Implication components give the exact simultaneous move budget

Write literal `(i,b)` for the assertion that packet `i` has state `b`;
its complement is `(i,1-b)`. A forbidden pair gives the implications

```text
(i,a) -> (j,1-b),     (j,b) -> (i,1-a).
```

Every implication is accompanied by its contrapositive. Let the strongly
connected components of this directed graph be called literal components.
The formula is feasible iff no component contains a literal and its
complement. For completeness, necessity follows from mutual implication.
For sufficiency, topologically order the component DAG and make a literal
true iff its component lies later than its complementary component. If a
true literal implied a false one, that edge and its contrapositive would
contradict the two strict ordering comparisons.

In a feasible formula, complementation pairs the literal components.
For a pair `C,bar C`, let `I_C` be all packet indices represented in it.
All these packet bits must change together whenever the truth of `C`
changes. The sets `I_C` partition the packet indices; their correlated
values need not all have the same sign.

Give packet `i` a positive additive flip cost `w_i`. Examples are the
number of changed owned path slots or the number of changed actual
slot-clock occurrences between its two full records. Singleton packets
never flip. Define

```text
W_C = sum_(i in I_C) w_i.
```

A component pair is **free** when each of its two truth choices occurs in
some satisfying assignment. Equivalently there is neither a path
`C -> bar C` nor a path `bar C -> C`. A path forces one truth choice.
Conversely adding the unary requirement that `C` be true adds the edge
`bar C -> C`; by the feasibility criterion it creates a contradiction
exactly when the original graph has `C -> bar C`.

**Theorem 2.** In a nonempty accepted assembly family, allow simultaneous
packet endpoint moves whose cost is the sum of weights of the changed
packets. The sharp uniform budget connecting the entire family is

```text
B* = max{W_C : C,bar C is free},        max(empty set)=0.           (1)
```

Moreover any two accepted assemblies admit a path that flips one
complementary component pair at a time, flips each changed packet exactly
once, and has total additive cost equal to their weighted Hamming distance.

**Proof of the constructive upper bound.** Fix accepted current and target
assignments. Their true literal components are upward-closed sets in the
implication DAG, each containing exactly one member of every complementary
pair. Let `D` be the components true currently and false at the target.
If `D` is nonempty, choose a component `C` with no predecessor in `D`.

There is no implication from a retained true component into `C`. A
currently true predecessor outside `D` would be true at the target, which
would force `C` true there as well; a predecessor inside `D` is excluded
by the choice of `C`. Thus turning `C` false breaks no implication from
the other retained true components.

Turn `bar C` true simultaneously. It cannot imply `C`, since the accepted
target makes `bar C` true and `C` false. If it implied any other component
`Y` that remains false, contraposition would give `bar Y -> C`. Apart from
the harmless consequent `bar C` itself, `bar Y` is a different currently
true predecessor of `C`, already excluded. Thus this addition also
preserves every implication.
The new state is accepted. The entire pair now has its target truth value,
so this step strictly removes one pair from the difference and never
revisits it. Iterate. Every changed pair is free because the two accepted
endpoints witness both choices. Each step costs at most (1). Every packet
that changes does so once, giving the asserted total cost. QED.

**Sharpness.** For a free pair choose two accepted assemblies that differ
on it. Any path between them must at some step change its truth value.
The equivalences within the two literal components force every packet in
`I_C` to change at that same step, at cost at least `W_C`. Taking the
largest such block proves the matching lower bound. QED.

This handles mixed positive/negative implications and frozen variables;
it is not a rule to flip positive-variable components indiscriminately.
The free component pairs are not independent choices. A single implication
`b0=1 -> b1=1` has three accepted states, not four; the component DAG must
remain attached to the block list. Likewise a genuine three-bit even-parity
guard has unrestricted pairwise projections but only four of eight global
states. That higher-arity predicate cannot be recovered from those
projections and is not part of this pairwise-only interface.
Costs of globally identified shared variables need not be additive over
packets and must be counted separately. There is also no claim about a
different prescribed elementary-move graph: moves here are simultaneous
authenticated local endpoint replacements with all assembled predicates
retained. Accepted intermediate single-owner permutations are not required.

## 4. An arbitrarily large block on genuinely different native sources

Use the original source1238, ordinary1763, canonical code29558785, and

```text
P0 = (3,0,6,7,1,5,4,2)
P1 = (3,1,6,7,0,5,4,2)
g  = (4,0,2,3,1,5,6,7).
```

The inherited direct map `g` is an isomorphism from `K_P0(T)` to
`K_P1(T)` and fixes both endpoint owners three and two. The core is
strongly connected: its actual arcs contain the directed spanning cycle

```text
0 -> 4 -> 1 -> 5 -> 2 -> 6 -> 7 -> 3 -> 0.
```

Pin all core clocks except one and four. The only remaining owners are
zero and one. There are exactly two possible permutations, namely `P0,P1`,
and both are Hamiltonian. This is a complete eligibility argument from
literal owner pins; it does not assume that padded target-class fibres
remain the same size. The core has trivial automorphism group, as already
proved by the inherited exact census and independently replayed here.

Choose integers `k>=1`, `n>=2k+7`, and `L=n-8`. On the common actual owner
set `{0,...,n-1}`, define `k` source tournaments

```text
T_r = R_r -> T -> S_(L-r),       r=0,...,k-1,
R_r owners = 8,...,8+r-1,
S_(L-r) owners = 8+r,...,n-1,
```

with transitive orders as printed. Both complete path states are the
unique prefix, `P_b`, and unique suffix. Pin every actual clock except
`r+1,r+4`. Again exactly those two records remain. The
[endpoint-fixed padding lemma](../../session44_guarded_native_moves_20261023/stable-padding/THEOREM.md)
extends `g` by the identity on new owners, proving equality of each slot's
full target isomorphism class. Therefore every target-dependent electrical
character at that new order is retained in each slot, at every modulus.
No numeric phase value from order eight is presumed to persist.

The ordered strongly connected component profile of `T_r` is `r`
singleton components, the unique strong eight-vertex core, and `L-r`
singletons. Hence different `r` give nonisomorphic sources. Converse
reverses that profile. Since `r+s<=2k-2<L`, none of these source classes
is the converse of another selected one. Their literal ports, including
converse partners, are consequently distinct without any resource reuse.

The inverse target endpoints also remain injective and disjoint from these
source endpoints. The old target `K_P0(T)` contains the actual strong cycle
`0 -> 2 -> 3 -> 4 -> 6 -> 7 -> 1 -> 5 -> 0`; its isomorphic state-one
target is strong as well. In every padded partner, a prefix vertex reaches
the core along the retained prefix path, and the core reaches the prefix
through a reversed cross arc and that same path. Similarly the core
reaches the suffix through its retained junction and suffix path, while
the suffix has a reversed arc back to a core owner other than the last
owner. The strong core supplies the internal connections. Thus every
padded target is strongly connected, whereas each selected padded source
is decomposable because `L>=1`. No source endpoint can equal a target
endpoint, even modulo converse.

Moreover the normalized inverse endpoint word is the bitwise complement
of the normalized source word; this involution makes inverse injection
equivalent to source injection. The trivial source automorphism groups
ensure the two different path states also have different normalized source
words. Altogether the two endpoint states, their converse copies, and both
ends of each port occupy exactly `8k` distinct normalized endpoint words.
This full reservoir, not just the forward source paths, must be avoided
by an unchanged background inventory.

Identify the actual free owner at clock `r+1` across all slots, for example
by a chain of equalities. The other hole at clock `r+4` contains the other
member of `{0,1}` and is then automatically identified as well. The fixed
actual arcs `(6,7)` and `(5,4)` remain at clocks `r+2` and `r+5`; these are
transported alignment coordinates, not unchanged printed clock numbers.
Exactly two global packets remain: all state zero and all state one.

**Corollary.** At every `n>=2k+7` this supplies a selected native packet on
`k` different source classes, even modulo converse, with sharp simultaneous
source-slot budget `k`. It changes `2k` actual slot-clock occurrences but
only two distinct identified shared owner variables. Keeping converse
copies as separately charged slots doubles the first two counts. All
target-dependent per-slot guards and the fixed shared arcs survive.

The bound is unbounded in `k`, despite the global feasible family having
only two states. State complexity, marking cost, and move support are
different quantities. Each source has trivial automorphism group, so
there is no omitted source-automorphism kernel hiding this support. The
construction uses many pins and decomposable padded sources; it does not
assert a low-cost marking theorem or a regular/primitive-source family.

As a native inconsistency control, use three such sources and identify the
first free owner of source zero with that of source one, then source one
with source two, but the first free owner of source two with the **second**
free owner of source zero. The induced state equations are
`b0=b1`, `b1=b2`, `b2=1-b0`. Each separate relation has legal native pairs;
their joint assembly is empty. The retained actual owner coordinate detects
the obstruction before any move or resource count is invoked.

## 5. Loss ledger, reproducibility, and stopping boundary

Map: authenticated full local binary records to their state bits, then
to complementary implication components. Preserved: every declared shared
owner equality, literal injection, fixed remainder, local guard, and the
exact additive coupled-move threshold. Lost from bits alone: actual source
and target ownership, paths, clocks and kernel alignment. Restoration
sidecar: the full two-record tables and explicit occurrence/resource maps.

The native construction can extend into an unchanged selected remainder
only when that remainder is supplied as `U` and the resulting formula is
satisfiable. The complete all-path metagraph already owns its literal
presentations; this theorem does not alter its complement involution or
manufacture unused ports in that atlas.

The checker exhausts all256 subsets of the eight distinct nontrivial
unary/binary forbidden clauses on two variables, all4,048 formulas with at
most four such clauses on three variables, and480 seeded controls on four
through six variables. Literal accepted-state enumeration and threshold
graph traversal independently test implication satisfiability, frozen
blocks, optimal budget, and the complete proposed monotone paths. Another
400 full-record packet instances test the lossless resource/owner compiler.
Native controls independently enumerate all40,320 core permutations,
verify the strong cycle and trivial automorphism group, and check36 padded
hosts for `k=1,...,8`, including all path/partner/converse arcs and pins.
They also verify strong padded targets, decomposable padded sources, and
the full `8k` converse-closed source/inverse-endpoint reservoirs.
The all-size statement is proved above, not inferred from that range.

```text
python -B 05-knowledge/results/injective_native_blocks_20261023/assembly/check.py
python -O -B 05-knowledge/results/injective_native_blocks_20261023/assembly/check.py
```

The exact source hashes and deterministic output are in `receipt.json`.
No floating comparisons or optimization-sensitive assertions are used.
The methods used are controlled forgetting with a sidecar and exposing the
obstruction before choosing a larger scale. This is a new native assembly
application and sharp support family, not a new general theory of 2-SAT,
nor an all-source wiring classification or Collatz proof.
