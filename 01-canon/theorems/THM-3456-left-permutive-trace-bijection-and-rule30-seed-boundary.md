---
id: THM-3456
title: "Left-permutive trace bijections and the Rule 30 seed boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every finite-alphabet radius-one cellular
  automaton permutive in its left input, the center trace together with the
  positive initial half-line is a homeomorphic coordinate system: every fixed
  right boundary compiles each finite target trace to one and only one left
  boundary.  Thus free-input trace cylinders are the full shift and are
  exactly uniform, while the distinguished single-seed trace is selected only
  by the discarded inverse-boundary word.  Rule 30 and Rule 60 have isomorphic
  free trace-incidence structures, but Rule 60's single-seed center is constant.
  No Rule 30 prize problem and no LRC(14) statement is solved.
source: root-rule30-260813291-20260815
depends_on: []
related:
  - THM-3395-small-sheet-typed-cover-star-cochain
  - THM-3453-global-literal-half-twist-cap-seven-support-classification
external:
  - Artem Chernikov, "SOP2 = SOP3", arXiv:2608.13291v1, https://arxiv.org/abs/2608.13291 (CITED VERY RECENT PREPRINT; model-theoretic corollary/context only)
  - Stephen Wolfram, "Announcing the Rule 30 Prizes", https://writings.stephenwolfram.com/2019/10/announcing-the-rule-30-prizes/ (2019; CITED for the problem statements, Boolean local law, sideways two-column discussion, and reported finite-size observations only)
  - Wolfram Rule 30 Prizes, https://rule30prize.org/ (CURRENT OFFICIAL LISTING checked 2026-08-15; problem and prize status only)
script: 04-computation/left_permutive_trace_torsor_rule30_thm3456.py
output: 05-knowledge/results/left_permutive_trace_torsor_rule30_thm3456.out
---

# THM-3456 -- left-permutive trace bijections and the Rule 30 seed boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

Let `A` be a finite alphabet of size `q>=2`, and let

```text
F:A^Z -> A^Z,
(F(x))_j=f(x_(j-1),x_j,x_(j+1))                       (1)
```

be a radius-one cellular automaton.  Assume that `F` is **left permutive**:
for every fixed `(b,c) in A^2`, the map

```text
a |-> f(a,b,c)                                        (2)
```

is a permutation of `A`.  Write

```text
tau_F(x)_t=(F^t(x))_0,                    t>=0         (3)
```

for the center trace.

Then the map

```text
Theta_F:A^Z -> A^N x A^N,
Theta_F(x)=(tau_F(x),(x_1,x_2,...))                    (4)
```

is a homeomorphism.  More precisely, for every `n>=1`, after fixing
`x_0,...,x_n`, the triangular map

```text
(x_-1,...,x_-n)
    |-> (tau_F(x)_1,...,tau_F(x)_n)                    (5)
```

is a bijection `A^n -> A^n`.

Consequently:

1. every length-`n+1` trace word occurs on exactly `q^n` initial blocks
   `x_-n,...,x_n`;
2. under the uniform Bernoulli product measure on `A^Z`, the center trace is
   an independent uniform `A`-valued process and is independent of the
   positive initial half-line;
3. all left-permutive rules on the same alphabet have isomorphic free
   trace-cylinder incidence structures;
4. the two-sorted incidence structure of configurations and finite trace
   words has a formula with `SOP_2` and, in fact, a definable strict cylinder
   order witnessing `SOP` and hence `SOP_3` directly.

For Rule 30,

```text
f_30(l,c,r)=l XOR (c OR r),                            (6)
```

so the theorem applies.  It also applies to Rule 60,

```text
f_60(l,c,r)=l XOR c.                                  (7)
```

The two rules therefore have identical free trace languages and isomorphic
trace-incidence structures.  Nevertheless, from the single-one initial
configuration `delta_0`, Rule 60 has the constant center trace `1^omega`.
Thus no invariant which factors only through the free trace language, its
uniform cylinder counts, or the resulting `SOP`/`SOP_2`/`SOP_3` theory can
prove any of the three Rule 30 prize assertions for the distinguished
single-seed orbit.

The missing coordinate is explicit.  Given a proposed trace `y` and a fixed
positive boundary `r`, let

```text
beta_F(y,r)=(x_-1,x_-2,...)                            (8)
```

be the unique left boundary supplied by `Theta_F^(-1)`.  The Rule 30
single-seed trace is exactly the unique `y` satisfying

```text
y_0=1,            r=0^omega,            beta_F(y,r)=0^omega.    (9)
```

This inverse-boundary word, or an equivalent distinguished-seed predicate,
is the mandatory sidecar.

## 2. Inheritance and connection contract

Wolfram's Rule 30 page isolates three questions about the center column from
one nonzero seed: absence of eventual periodicity, limiting balance, and a linear
computational lower bound.  It also records the sideways identity behind
Rule 30's two-column arguments.  Chernikov's very recent preprint proves
`SOP_2=SOP_3` while retaining witness--parameter pairs through a mixed partial
type dichotomy.  THM-3395 independently imported that typing move into its
elementary small-sheet star-cochain theorem.

The present connection is exact but deliberately negative:

| field | value |
|---|---|
| source | a left-permutive cellular automaton with a full initial configuration |
| target | its center trace, or the incidence structure of trace cylinders |
| map | `x |-> tau_F(x)`; faithfully, `x |-> (tau_F(x),x_(>0))` |
| preserved by the bare trace | every finite center word and prefix consistency |
| destroyed by the bare trace | the right boundary, the uniquely compiled left boundary, and membership of the named single seed |
| required sidecar | `x_(>0)` together with `beta_F`, or the complete initial configuration |
| cheapest hostile | Rule 30 versus Rule 60: identical free trace-cylinder theory, constant Rule 60 single-seed center |

The model-theoretic source and the Rule 30 target therefore meet at a precise
guardrail: a theory of all compatible branches may have the strongest tree
behaviour while saying nothing about one named branch.  This does not
interpret an LRC sheet system in a Rule 30 spacetime diagram, and it supplies
no physical clock, owner, phase, or LRC current.

## 3. The triangular finite trace compiler

For `t>=1`, the value `tau_F(x)_t` depends only on the initial coordinates
`x_-t,...,x_t`.  We claim that, with all other coordinates in this interval
fixed, it is a permutation of the extreme coordinate `x_-t`.

For `t=1` this is exactly (2).  Suppose the claim holds at time `t-1`.  Then

```text
(F^t(x))_0
 = f((F^(t-1)(x))_-1,
     (F^(t-1)(x))_0,
     (F^(t-1)(x))_1).                                 (10)
```

The last two arguments in (10) do not depend on `x_-t`.  By shift
equivariance, the induction hypothesis applied at site `-1` says that the
first argument is a permutation of `x_-t`; applying the permutation (2) once
more preserves that property.  This proves the claim.

Now fix `x_0,...,x_n` and a target trace suffix `y_1,...,y_n`.  At stage `t`,
the already selected coordinates `x_-1,...,x_(-(t-1))` determine every input
to `tau_F(x)_t` except `x_-t`.  The claim supplies a unique value of `x_-t`
for which `tau_F(x)_t=y_t`.  Induction on `t` gives one and only one left
word, proving (5).

If the target includes `y_0`, then necessarily `x_0=y_0`.  The positive word
`x_1,...,x_n` remains arbitrary, giving exactly `q^n` preimage blocks for
every target trace prefix.  This proves the exact fibre count.

## 4. Infinite homeomorphism and Bernoulli law

Given `(y,r) in A^N x A^N`, set `x_0=y_0` and `x_j=r_j` for `j>=1`.  Applying
the finite compiler successively constructs compatible values
`x_-1,x_-2,...`; uniqueness at every finite depth makes the resulting
configuration the unique inverse image of `(y,r)` under (4).

Every output coordinate of `Theta_F` depends on finitely many input
coordinates, and every inverse coordinate `x_-t` depends only on

```text
y_0,...,y_t,              r_1,...,r_t.                (11)
```

Thus `Theta_F` and its inverse are continuous in the product topologies,
proving the homeomorphism.

On a cylinder fixing the coordinates in (11), the finite compiler is a
bijection between two sets of size `q^(2t+1)`.  Hence `Theta_F` pushes the
uniform Bernoulli measure to the product of the two uniform Bernoulli
measures.  In particular the free-input trace is exactly iid uniform at every
finite depth.  No limiting or mixing argument is being used.

For two left-permutive rules `F,G`, the homeomorphism

```text
H_(F,G)=Theta_G^(-1) Theta_F                             (12)
```

preserves the positive boundary and the complete center trace.  It is a
boundary recoding, not a conjugacy of the time evolutions.

## 5. The trace-cylinder tree and the SOP boundary

Let `M_F` be the two-sorted first-order structure with configuration sort
`A^Z`, word sort `A^(<omega)`, and incidence relation

```text
E(x;w)  iff  w is an initial segment of tau_F(x).      (13)
```

Choose two alphabet symbols and identify `2^(<omega)` with the corresponding
word subtree.  For each branch `sigma in 2^omega`, trace surjectivity supplies
an `x` satisfying every formula

```text
E(x;sigma|n),                 n<omega.                 (14)
```

If `eta,nu` are incomparable nodes, no trace can have both as prefixes, so
`{E(x;eta),E(x;nu)}` is inconsistent.  Therefore the single formula `E(x;y)`
witnesses `SOP_2` in `Th(M_F)`.

There is also a direct stronger witness.  Define on the word sort

```text
Q(u;v) :=
  [forall x (E(x;v) -> E(x;u))]
  and [exists x (E(x;u) and not E(x;v))].              (15)
```

This is strict inclusion of the corresponding trace cylinders.  The words
`0^i`, `i<omega`, form an infinite `Q`-chain, and strict inclusion has no
directed cycle.  Hence `Q` witnesses the strict order property and in
particular `SOP_3`.

Chernikov's arXiv:2608.13291v1 proves in full generality that `SOP_2` implies
`SOP_3`, using typed witness--parameter packets even when no definable nested
cylinder order is available.  Equation (15) gives a self-contained special
case, not a new proof of the external theorem.  Conversely, the fact that
every left-permutive rule produces this same highly non-structured incidence
theory shows why the classification result cannot distinguish Rule 30's
named seed from Rule 60's.

## 6. The Rule 60 single-seed hostile

For Rule 60, every site strictly left of the origin remains zero when the
initial state is `delta_0`, because (7) uses only that site and its left
neighbour.  The center therefore satisfies

```text
x^(t+1)_0=x^t_-1 XOR x^t_0=0 XOR 1=1                 (16)
```

for every `t>=0`.  Its center trace is constant, so it is periodic, has
limiting one-frequency `1`, and is computable in constant time.

But Rule 30 and Rule 60 share all conclusions of Sections 3--5.  This gives
three separate non-implications:

```text
full free trace language  -/->  single-seed nonperiodicity,
free Bernoulli balance    -/->  single-seed limiting balance,
trace-cylinder SOP        -/->  single-seed time lower bound.          (17)
```

As checked on 2026-08-15, the current official page still lists all three
prizes and states that submissions are accepted until a satisfactory solution
is achieved.  On that evidence all three questions are treated here as open.
The theorem neither answers nor conditionally reduces any of them.

## 7. Rule 30 finite rings: two columns are faithful, one is not

On a cyclic ring `Z/nZ`, write `x^t_j` for the Rule 30 spacetime diagram.
Solving (6) for its left input gives the exact sideways inverse

```text
x^t_(j-1)=x^(t+1)_j XOR (x^t_j OR x^t_(j+1)).         (18)
```

Thus the two adjacent temporal columns at `j,j+1` determine the column at
`j-1`; iterating (18) around the ring determines the complete spacetime
orbit.  If the two input columns are eventually periodic with periods
`p_0,p_1`, then every reconstructed column and hence the whole orbit is
eventually periodic with period dividing `lcm(p_0,p_1)`.  This is the exact
finite-ring form of the two-column mechanism discussed on Wolfram's page.

One column is sharply insufficient.  In coordinate order `0,...,6`, the
three width-seven periodic states

```text
(0,1,0,0,0,0,0),
(0,1,1,1,0,0,0),
(0,0,0,0,1,1,1)                                      (19)
```

all have center trace `(0110)^omega`, but their adjacent columns have the
distinct period-four traces `1101`, `1100`, and `0010`.  Direct exhaustive
audit of all temporally periodic ring states finds no nonconstant
one-column alias in widths `3,...,6`.

Finally, the width-`n` cyclic single seed lifts to the infinite periodic seed
set `nZ`.  Its center trace agrees with the isolated infinite seed for every
`0<=t<n`, because the radius-`t` cone then contains no other lift.  The bound
is sharp as a uniform statement: at `n=3`, the cyclic and isolated traces
first differ at `t=3`.  Eventual cycles of finite rings therefore certify
finite prefixes, not the infinite single-seed tail.

## 8. Exact verification and failure boundaries

The companion is designed to verify independently of the proof:

1. all `16` binary radius-one left-permutive elementary rules;
2. the finite triangular bijection and exact trace-fibre counts on exhaustive
   small horizons, including Rule 30 and Rule 60;
3. the Rule 60 constant-seed hostile and Rule 30 prefix controls;
4. the sideways inverse on all local Rule 30 cells;
5. every periodic state in widths `3,...,7`, including (19), its three
   adjacent traces, and the first nonconstant one-column alias boundary;
6. the cyclic/infinite single-seed prefix equality and the sharp `n=3`
   boundary.

The frozen universe contains all `16` binary left-permutive elementary rules
and all fixed right blocks through depth `6`, with Rules 30 and 60 extended
through depth `8`: `5,568` fixed-right maps, `502,400` forward traces, and
`5,568` constructive inverse replays.  It exhausts periodic states in widths
`3,...,7`, compares seed prefixes through width `64`, and records exact seed
cycles through width `24`.  Normal, optimized, and stored transcripts are
byte-identical.  LF-normalized SHA-256 pins are

```text
script   e7ac9717c2f5618d18cae2225c4fd68391d2e07d98f0e9b1bcb188e5ef7ad993
output   a04fc7d1c6124888b772f10f0ee5ab22d62c43b669b56f15a440de4e51bce9fd
semantic 1fda893b46126cb9653f758a88a9f0266fe93c6e12bee6d4de527fa09f5b42ca
```

The proof requires genuine left permutivity.  If (2) fails, a target trace
symbol can have zero or multiple choices of `x_-t`, and (5) fails already at
the first corresponding context.  Right permutivity gives the reflected
statement with the negative and positive half-lines exchanged.  Larger
radius rules admit the analogous extreme-input compiler, but its boundary
packet has more than one new symbol per time step and is not asserted here.

The `SOP` statement belongs only to the deliberately enriched free-input
trace-incidence structure (13).  It is not a claim about a canonical bare
first-order theory of Rule 30 spacetime, and it is not a property of the one
distinguished trace.  Likewise, exact Bernoulli balance under random initial
conditions is not Problem 2, which asks about `delta_0`.

No LRC statement follows.  A Rule 30 bit word carries neither the odd
`13`-phase translation, the physical common time, nor the owner-labelled
current needed at the present LRC(14) frontier.  THM-3395 remains the nearest
typed-witness analogy: its affine star cochain restores a common phase, while
`beta_F` restores the missing initial boundary here.  There is no map between
their target predicates.

The independent unnumbered
[Rule 30/LRC depth-observer probe](../../07-reflections/lrc14-rule30-depth-observer-typed-phase-no-go-codex-20260815.md)
makes that boundary finite-exact on the current `165` first-depth-one profile
bank.  Sampling the single-seed Rule 30 trace at the three profile depths has
only four fibres, of sizes `36,36,51,42`, and two fibres mix repeated-owner
with strict profiles.  Reattaching the full profile makes the bit field a
deterministic redundant graph.  Separately,
`Hom(F_13^2,F_2^m)=0` for every finite `m`, so an additive/XOR observer cannot
carry the required odd phase; a nonlinear four-bit encoding can carry a
13-cycle only by explicitly restoring thirteen labelled states.  The probe
excludes zero scalar rows and records `LRC(14)` as open.

## 9. Independent audit

Two independent immutable-package audits of commit `1259c28e5521` rederived
the finite compiler, inverse-limit homeomorphism, Bernoulli product law,
fixed-formula `SOP_2`, direct strict-order/`SOP_3` witness, Rule 30/Rule 60
hostile, sideways reconstruction, alias boundary, and periodic-lift prefix
law.  Separate clean-room implementations reproduced every triangular census
row, the periodic-state profiles through width `7`, all canonical alias
representatives and adjacent traces, and the prefix checks through width `64`.
Normal, optimized, and stored outputs and all three declared hashes matched.

A third primary-source audit reread arXiv:2608.13291v1, the 2019 announcement,
and the current prize page.  MISTAKE-403 records the repaired citation scope
and the distinction between a page's active listing and this repository's
dated openness inference.  The external model-theory theorem remains cited
context only, and none of the three prize questions or LRC(14) is claimed.

**QED.**
