---
source: codex-padic-zeta-tournament-20260825
status: >
  RESEARCH SYNTHESIS. THM-4088--4091 and THM-4093 are proved at their stated scopes; the
  external 22-value p-adic-zeta theorem and arXiv:2608.13306v1 remain
  AUTHOR-CLAIMED / UNREFEREED. LRC(14), JC(2), the Sun-hole classification,
  and the named tournament inequalities remain OPEN.
external:
  - https://github.com/octonion/p-adic-zeta-irrationality/commit/b46a1770901551961710e155d775aae7c5ea39e7
  - https://arxiv.org/abs/2608.13306v1
---

# Orientation is not arithmetic type: a p-adic, logic, and tournament session

## Verdict

The two supplied sources do not form one mathematical method. Long's
manuscript is an arithmetic-geometry preprint claiming 22 p-adic-zeta
irrationalities; arXiv:2608.13306v1 is Chen--Roşu's matching-logic preprint.
Their only lawful common theme is **controlled forgetting**: an observable may
retain enough information for one local decision while deleting the global
coordinate needed by the target theorem.

The combined concurrent session produced five scoped advances:

1. [THM-4088](../01-canon/theorems/THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density.md)
   proves that tie-free order and approximation-quality tournaments are blind
   to rational, algebraic-irrational, and transcendental type.
2. [THM-4089](../01-canon/theorems/THM-4089-hybrid-padic-zeta-margin-optimization-and-next-case-obstruction.md)
   globally optimizes the external displayed margin, verifies its 22 formula
   signs, proves four immediate next cases impossible by parameter retuning,
   and separately proves an all-radius `p=13,s=3` one-power no-go.
3. [THM-4090](../01-canon/theorems/THM-4090-two-sort-matching-logic-global-completeness-obstruction.md)
   sharpens the matching-logic preprint's three-sort counterexample to two
   sorts by an independently audited finite argument.
4. [THM-4091](../01-canon/theorems/THM-4091-integral-coordinate-change-lcm-depth-boundary.md)
   proves that cumulative LCM clearing survives every integral formal
   coordinate change, while literal coefficient depth `e>=2` already fails
   at output degree three.
5. [THM-4093](../01-canon/theorems/THM-4093-rational-edge-diagonal-gauge-and-padic-tournament-zeta-tangent.md)
   identifies the exact diagonal gauge preserved by determinant/Bowen--Lanford
   zeta and the directed-triangle divisibility boundary of its first p-adic
   tangent.

No headline conjecture is claimed solved.

## Inheritance pass and concept board

The closest proved mechanism was THM-4057: reduced rationals carry coprime
arcs, while ordinary order supplies only a transitive completion. The
canonical hostile was AP13, whose LRC margin is exactly zero. The corrected
near miss was the recurrent temptation to identify a shared clock or sign
pattern with a target-preserving map. The least-used sidecar was HYP-3114's
arithmetic-type packet for strict witness intervals.

The board stabilized at five objects:

| live object | representation | invariant kept | operation / loss |
|---|---|---|---|
| rational approximants | ordered vertices | orientation sign | forgets determinant, height, residual rate |
| p-adic approximants | valuation-ranked vertices | error-depth order | forgets adelic height and product-formula contradiction |
| holonomy certificate | prime/pivot weighted incidence | additive local cost | total orientation destroys ties and magnitudes |
| matching-logic models | localization / sort-flow graph | reachable consequence | projection loses an unreachable carrier |
| LRC witness set | union of exact intervals | margin and endpoints | arithmetic-type label forgets the boundary geometry |

Comparing every pull against this board led to one stable rule: signs and
reachability can expose a dependency order, but arithmetic conclusions require
a quantitative sidecar that survives local-to-global assembly.

## Why arithmetic types do not make a tournament hierarchy

Given distinct rational approximants, orienting the smaller toward the larger
always produces a transitive tournament. THM-4088 gives explicit increasing
sequences converging to `1/2`, `sqrt(2)`, and a Liouville number with the
identical labelled tournament `i->j iff i<j`; every finite increasing prefix
has continuations of all three types. Ranking by a prescribed strict p-adic
error-valuation schedule has the same collapse for every p-adic target.

So there is no intrinsic cyclic pattern such as

```text
rational -> algebraic irrational -> transcendental -> rational.
```

The categories are not pairwise observables. A forced orientation merely
encodes a chosen complexity score. A faithful arithmetic object must retain

```text
sign + determinant magnitude + denominator height
     + target residual + denominator clearing + decay rate
     + degree/height lower-bound fence.                         (1)
```

Nontransitivity could become meaningful only after choosing an intrinsic
antisymmetric comparison with target content—for example, competing
approximants evaluated under different local places and a preserved adelic
budget. The p-adic manuscript's prime/pivot data are instead a weighted
bipartite incidence structure with ties and hyperedges, not a tournament.

## Anchor: what changed for LRC(14)

For integer speeds `S`, the function

```text
F_S(t)=min_(s in S)||s t||
```

is Lipschitz. Every strict component `F_S(t)>1/14` is an open interval, hence
contains rational, algebraic-irrational, and transcendental times. Arithmetic
type cannot distinguish an LRC counterexample from a safe speed family in the
interior. This closes one speculative lane cleanly.

The survivor is exact and narrower:

```text
F_S(t)=1/14,       margin=0.                              (2)
```

At `(2)`, approximation stability disappears. The next analytic object should
not be a tournament of named constants; it should be the interval-local Stern
sum

```text
S_I(q)=sum_(a mod q, a/q in I) chi(a) e(h a/q),           (3)
```

with endpoint ownership and denominator height retained. THM-4071 controls a
complete all-odd packet; it does not automatically control the truncated
interval `(3)`. A decisive hostile is AP13, where every putative positive
component has collapsed to the boundary.

## Niche: the exact p-adic margin frontier

The external certificate really does verify positive values of its displayed
margin for the claimed 22 parameter pairs. THM-4089 separates that finite
numeric fact from the manuscript's new geometry. It proves an explicit global
minimizer in the cutoff variable and strict concavity with a unique maximizer
in the analytic-radius variable. Exact tangent enclosures then give negative
global optima for

```text
(p,s)=(2,31),(3,13),(5,7),(7,5).                         (4)
```

Thus `(4)` is a mechanism-level stopping result: more precision or a finer
grid cannot extend any of the four rows by one odd weight. The cost theorem
itself must improve.

The orthogonal source-level extrapolation fails too. The first-chamber
stationary point crosses `xi=1` at `p=11`; a complete chamber calculation at
`p=13,s=3` gives idealized `tau>=613/288` and margin `<-37/72` for every
continuation radius, with actual formal margin `<-67/56`. This does not decide
the arithmetic type of `zeta_13(3)`; it rejects one displayed architecture.

The cheapest serious audit is an explicit `(p,s)=(5,5)` worked case with a
named frame, torsor chart, divisor and exceptional primes, source bases, one
small-prime Hasse matrix and kernel, and one large-prime pole-grade pivot.
Until then, the 22 irrationalities remain **AUTHOR-CLAIMED / UNREFEREED;
SPECIALIST AUDIT OPEN**.

## Wildcard: matching logic and the two-sort wall

Chen--Roşu's preprint is not an arithmetic paper. Its localization mechanism
did, however, suggest testing whether the reported three-sort failure was
minimal. THM-4090 proves it is not: two sorts `a,b`, a unary `f:b->a`, and

```text
Gamma={forall x:b forall y:b. f(x and y)}
```

already force the `b`-carrier to be a singleton semantically, while the
`a`-sorted hypothesis cannot feed a `b`-sorted proof. Independently this proves
two-sort failure; together with the preprint's claimed one-sort theorem, it
would make the sort-count boundary sharp for the displayed calculus.

Sort flow is a preorder, not a tournament. The result is useful as a
formalization firewall: a complete calculus for an encoding does not prove
the semantic premise, validate the encoding, or construct the mathematical
witness.

## Typed connections and cheapest decisive tests

| source | target | map | preserved predicate | destroyed information | needed sidecar / cheapest test |
|---|---|---|---|---|---|
| rational approximants | arithmetic type | order sign | pairwise order | gap size, heights, residual rate | determinant-height ledger; three-type continuation hostile |
| p-adic approximants | irrationality | valuation rank | local error depth | global height, nonzero integer form | adelic/product-formula certificate |
| holonomy LCM depth | Apéry framework | prime valuations | denominator budget | recurrence, nonvanishing, root decay | explicit recurrence and cleared linear form |
| binomial atoms mod `p^a` | Sun 2-4-6-8 | bounded lift tree | congruence with height | unrestricted equality and cancellation | height-aware branch-and-bound; universal-solubility hostile |
| strict LRC component | finite grid | density in interval | witness inequality | endpoint owner at zero margin | interval-local Stern estimate; AP13 hostile |
| matching localization | repo formalization | semantics-preserving encoding | consequence if encoding is valid | original object and quantitative witness | translation soundness plus independent semantic proof |
| `L_5=60` | AP/Fibonacci/triangular clocks | common modulus only | residue address | recurrence owner and phase dynamics | clock-specific state transition, not numerical equality |

## New frontier portfolio

### Anchor — LRC(14)

1. Prove an endpoint-sensitive bound for interval-local Stern sums `(3)` or
   construct a boundary hostile showing square-root cancellation is
   insufficient after truncation.
2. Attach endpoint owner, denominator height, and exact zero-margin flags to
   the existing HYP-3114 observer-gluing packet.
3. Keep the `60`-clock split strict: `L_5=60`, the Fibonacci Pisano clock, and
   the AP/triangular phase law share an address modulus but not a transition
   map or owner.

### Niche — p-adic zeta and Sun

1. Build the explicit `(5,5)` geometric audit before trusting the external
   arithmetic theorem.
2. In the Sun problem, use prime-power information only in bounded-height,
   role-labelled lift trees; THM-4027 rules out a fixed local obstruction.
3. Search for a stronger p-adic margin ingredient guided by the exact deficits
   in THM-4089, rather than retuning `xi,Y`.

### Wildcard — logic, JC, and tournament invariants

1. Test whether the two-sort obstruction persists under the smallest added
   `a->b` return symbol; this is the first sort-flow boundary where the proof
   above genuinely breaks.
2. For planar JC, retain the current algebraic frontiers: the supplied sources
   provide no Keller map, Newton edge, degree owner, nonproperness curve, or
   termination argument. A logical encoding is not a JC advance.
3. For tournaments, seek weighted or decorated invariants whose forgetful map
   to orientation is explicit. Any claimed arithmetic correspondence must
   first defeat the three-type continuation hostile of THM-4088.

The session's durable synthesis is therefore not that irrational numbers
secretly form tournaments. It is that a tournament is often the **sign
quotient** of a richer arithmetic comparison, and the open theorem lives in
the information that quotient destroys.
