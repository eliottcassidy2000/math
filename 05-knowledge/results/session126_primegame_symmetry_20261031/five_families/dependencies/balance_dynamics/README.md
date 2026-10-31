# Divisor balance under the ordered fourteen-fraction program

**Session 126, 2026-10-31.** **INHERITED:** the definition of `F,U,S` and the
classification of balanced integers. **PROVED HERE:** the spectator reduction,
the complete balance-preserving transition classification, its exact restricted
orbit consequence, and the failure of the unlabelled observation quotient.
**FINITE-EXACT:** all representative transition tables and independently checked
counts. No novelty claim is made for the inherited divisor classification or
for Conway's program. This is not a general orbit classification of FRACTRAN.

## Results first

| Balanced-input universe | Inputs | Balanced after one step | Unbalanced after one step | Unit exit |
|---|---:|---:|---:|---:|
| Supported on the first ten primes | 380 | 24 | 354 | 2 |
| Arbitrary prime support, modulo renaming outside primes | 538 | 39 | 497 | 2 |
| First ten primes, valid control carrier | 98 | 15 | 81 | 2 |
| Arbitrary support, valid control carrier, modulo outside renaming | 181 | 29 | 150 | 2 |

These are complete one-step universes, not bounded integer searches. The two
unit exits are `7 -> 1` and `17 -> 1`. The unit is **outside** the stipulated
`N>=2` divisor-balance domain, not an integer to which `tau(N)-2` may be applied
as a divisor count.

There is also a complete answer to the narrower question “which orbits remain
balanced forever?” They are exactly

    11 <-> 13,
    11*a^2*b <-> 13*a^2*b,

where `a,b` are distinct primes in `{2,5}` or outside the first ten primes.
Every other balanced starting integer exits the balanced set within three
steps; in the valid control carrier it exits within two. This proves nothing
about its unrestricted orbit after that exit, nor does it rule out a later
return to balance.

For the distinguished **seed 2 orbit**, the reconstructed divisor-loop proof
gives a stronger, separate conclusion: there are **exactly 16 balanced visits
in the entire infinite orbit**, the last being `60` at update 70. Section 6a
proves the tail exclusion by a token-count lower bound, not extrapolation
from the first ten prime outputs.

## 1. Inheritance, carrier, and exact program

Read-only source snapshot: `8f76b84d2960f9e14ad62a154b8533c77593702f`.
The exact source paths and byte hashes are recorded in `sources.json`:

- `05-knowledge/results/trijection_reset_20261020/arithmetic.md`, §1:
  definitions, the inherited classification, and the warning that `U` is a
  subset of `S`, not a disjoint summand.
- `05-knowledge/results/session84_guarded_unification_20261028/FUS_DIVISOR_FRONTIER.md`,
  §§1,5: exact counts and the labelled-exponent observer criterion.

For `N>=2`, let `D(N)={d:d|N,1<d<N}`. The functions count all divisors in `D`,
the prime elements, and the squarefree elements, respectively:

    F = tau(N)-2,
    U = omega(N)-[N prime],
    S = 2^omega(N)-1-[N squarefree].

They are not divisor sums. The inherited theorem is

    F=U+S  iff  N=p, p^3, or p^2*q*r,

where prime letters in the mixed shape are distinct. Write `B` for this set.
In particular the possible sums of exponents `Omega(N)` in `B` are `1,3,4`.
The closest inherited mechanism is the labelled prime-exponent box; the hostile
example to count-only dynamics is replacement of a named prime coordinate.

The active alphabet is

    P = {2,3,5,7,11,13,17,19,23,29}.

The ordered program is exactly

| Rule | Fraction | Consumed prime tokens | Produced prime tokens |
|---:|---:|---|---|
|1|17/91|7,13|17|
|2|78/85|5,17|2,3,13|
|3|19/51|3,17|19|
|4|23/38|2,19|23|
|5|29/33|3,11|29|
|6|77/29|29|7,11|
|7|95/23|23|5,19|
|8|77/19|19|7,11|
|9|1/17|17|none|
|10|11/13|13|11|
|11|13/11|11|13|
|12|15/2|2|3,5|
|13|1/7|7|none|
|14|55/1|none|5,11|

All fractions are reduced. The action `T(N)` uses the **first** denominator
dividing `N`; order cannot be replaced by an unordered rewrite relation.
The last rule makes `T` total on positive integers, including the unit.
The identification with Conway's classical PRIMEGAME and its seed-2 output
theorem are separately sourced and audited in `../fractran/`; they are not
needed to prove the transition classification below.

## 2. Why 538 representatives cover arbitrary primes

Let `O` denote all primes outside `P`, and factor uniquely `N=A*M`, where `A`
is supported on `P` and every prime divisor of `M` lies in `O`. Every program
denominator divides a product of primes in `P`, so for every rule denominator
`b`,

    b | N  iff  b | A.

Consequently the first enabled rule is the same at `N` and `A`, and

    T(A*M) = T(A)*M.                                      (1)

Thus projection `pi(N)=A` is a semiconjugacy `pi T = T pi`, on positive
integers with `A=1` allowed. Every individual outside exponent is unchanged.
Outside-prime renamings commute with the fixed program. They also preserve
the exponent profile and hence preserve `F,U,S` wherever `N>=2`.

The balanced shapes use at most three primes. Therefore three distinct
outside representatives, chosen as `31,37,41`, suffice for **all** support
patterns, with no magnitude cutoff. On the active support there are

    10 + 10 + 10*C(9,2) = 380

balanced integers. Modulo permutations of `O`, primes and cubes each have
11 types. Mixed shapes have

    repeated prime active:   10*(C(9,2)+9+1) = 460,
    repeated prime outside:  C(10,2)+10+1    = 56.

The three summands in each line mean zero, one, or two outside singleton
primes, after taking account of whether the repeated prime is itself outside.
Thus the complete orbit-type count is `11+11+460+56=538`. These are orbit
counts under outside-label permutations, not cardinalities of finite sets
of arbitrary integers. The full 538-row table is `transitions.tsv`.

For comparison, enumerate all balanced integers on the thirteen explicit
primes `P union {31,37,41}`. There are `13+13+13*C(12,2)=884`. Independently
quotienting these by the unordered outside exponent profile recovers the
same 538 types; both checkers verify this identity.

### Renaming all primes requires renaming the program

If `rho` is an injective relabelling of prime coordinates and is applied to
both integers and every program numerator and denominator, then

    rho(T_program(N)) = T_rho(program)(rho(N)).             (2)

Divisibility guards and their order are unchanged by an injective exponent
coordinate relabelling, proving (2). This is a conjugacy onto the relabelled
carrier, not an equivariance claim for arbitrary active-prime permutations
of the **fixed** program. The fixed-machine rule at `11` is rule11, whereas
the rule at `13` is rule10. An arbitrary numerical substitution that maps
distinct primes to the same prime is not covered.

## 3. Complete universal balance-preserving rule classification

The following list is an iff, always with the stated first-enabled rule.
All variables denoted by lower-case letters in the families are primes, and
`a!=b`. In the displayed families the allowed prime sets already enforce
every earlier-rule exclusion.

| Rule | All balanced inputs with balanced output | Output |
|---:|---|---|
|1|`7*13*17^2`|`17^3`|
|2|none|—|
|3|`3*17*19^2`|`19^3`|
|4|`2*19*23^2`|`23^3`|
|5|`3*11*29^2`|`29^3`|
|6|`29^3`|`29^2*7*11`|
|7|`23^3`|`23^2*5*19`|
|8|`19^3`|`19^2*7*11`|
|9|none|—|
|10|`13`; `13*a^2*b`, `a,b in {2,3,5} union O`; `11*13^2*t`, `t in {2,5} union O`|replace one 13 by 11|
|11|`11`; `11*a^2*b`, `a,b in {2,5,7} union O`|replace 11 by 13|
|12|`2^3=8`|`2^2*3*5=60`|
|13|none|—|
|14|none|—|

**Proof of completeness.** Rules 1,3,4,5 consume two distinct primes and
produce one prime, reducing `Omega` by one. A balanced prime or cube cannot
enable them. A mixed balanced input can remain balanced only by becoming
a cube; all surviving tokens must therefore have the numerator's prime
label. This forces exactly the four listed inputs. Their earlier guards
are absent, so the listed rules really are first enabled.

Rule2 changes `Omega` from four to five on any balanced input that enables
it, so cannot preserve balance. Rules6,7,8,12 increase `Omega` by one.
For a balanced input, only a cube can then reach an allowed total, namely
four. Divisibility forces the cube to be the denominator prime, and the
result has the mixed balanced shape, giving exactly the four listed cases.
Rules9 and13 remove one prime token. A prime goes to the excluded unit;
a cube becomes a square; a mixed input still has at least two prime labels
and cannot become a cube. Rule14 increases `Omega` by two: only a prime
could possibly reach total three, but multiplication by the distinct
primes `5,11` cannot give a cube.

Rules10 and11 preserve `Omega`, so the exponents must be examined. At
rule10, the balanced-preserving condition is exactly

    (v13=1 and v11=0) or (v13=2 and v11=1).                 (3)

The first case replaces a singleton label (including `13 -> 11`); the
second transfers the repeated-prime position from 13 to 11. All other
possibilities produce `(2,1)`, `(3,1)`, `(2,2)`, or `(1,1,1,1)` rather
than a balanced profile. Earlier-rule exclusion in the first case forbids
`7,17,19,23,29`; in the second it allows the third prime only in
`{2,5} union O`. This gives the rule10 row. At rule11, priority already
forces `v13=0`, and balance is preserved iff `v11=1`. Earlier-rule
exclusion additionally forbids `3,17,19,23,29`, giving its row. This proves
every necessity and sufficiency in the table.

Rule10 has 7 regular active types, or 14 regular types with outside labels;
its repeated-position transfers add 2 active types or 3 outside types.
Rule11 has 7 active or 14 outside types. Together with the eight exceptional
rows this gives respectively `8+9+7=24` and `8+17+14=39` preserved types.

### Per-rule exact counts

Each pair below is `inputs / balanced outputs`. The complete output exponent
profile distributions are in `classification.json`, not suppressed in the
machine-readable result.

|Rule|Active 380|All 538 types|Active-valid 98|All-valid 181 types|
|---:|---:|---:|---:|---:|
|1|24 / 1|27 / 1|6 / 0|8 / 0|
|2|24 / 0|27 / 0|6 / 0|8 / 0|
|3|21 / 1|24 / 1|4 / 0|6 / 0|
|4|24 / 1|27 / 1|6 / 0|8 / 0|
|5|21 / 1|24 / 1|6 / 0|8 / 0|
|6|95 / 1|124 / 1|13 / 0|22 / 0|
|7|71 / 1|97 / 1|13 / 0|22 / 0|
|8|35 / 1|55 / 1|7 / 0|14 / 0|
|9|17 / 0|31 / 0|3 / 0|8 / 0|
|10|17 / 9|31 / 17|7 / 7|14 / 14|
|11|11 / 7|22 / 14|7 / 7|14 / 14|
|12|11 / 1|22 / 1|11 / 1|22 / 1|
|13|5 / 0|13 / 0|5 / 0|13 / 0|
|14|4 / 0|14 / 0|4 / 0|14 / 0|

## 4. Valid control carrier and the sharper preservation theorem

Let the four data primes be `{2,3,5,7}` and the six control primes be
`{11,13,17,19,23,29}`. A valid state has at most one control token in
total: its control-prime exponents sum to at most one. Data and outside
prime exponents are unrestricted. This is a definition of a carrier, not
a claim that every state in it is reachable from seed 2.

The carrier is forward invariant. A state with a control token uses one
of rules1–11: each consumes that token, produces exactly one control token
or none, and creates no extra one. Every occupied control prime has an
unconditional-denominator rule among1–11, so rules12–14 cannot be selected
while a control token is present. With no control token, rules12–13 keep
none and rule14 introduces exactly one. Outside primes do not affect this
argument by (1).

Active balanced inputs now comprise 10 primes, 4 data-prime cubes, and

    4*(C(3,2)+3*6) = 84

mixed shapes, for total 98. With outside primes up to renaming, there are
11 primes, 5 cubes, and

    repeated data prime:    4*((C(3,2)+3+1)+6*(3+1)) = 124,
    repeated outside prime:   (C(4,2)+4+1)+6*(4+1)  = 41

mixed shapes, for total 181. The two groups in each mixed count have zero
or one control singleton. Two control singletons are excluded.

In this valid carrier, rules10 and11 are pure singleton label changes and
preserve balance whenever they are selected from a balanced input. All
control-multiplicity exceptions in §3 are excluded. The **only** other
balanced-preserving transition is `8 -> 60` by rule12. Hence the preserved
counts are `7+7+1=15` active and `14+14+1=29` with outside types.

## 5. Complete graph restricted to balanced states

Retain a directed edge only when both its endpoints lie in `B`. By §3,
the overlap of the regular rule10 and rule11 families is exactly

    11 <-> 13,
    11*a^2*b <-> 13*a^2*b,  a,b in {2,5} union O, a!=b.    (4)

These are genuine two-cycles of the full program. The other regular rule10
cases contain prime3; after their label change, rule5 is enabled and the
next output is not balanced. The other regular rule11 cases contain prime7;
after the label change, rule1 is enabled and the next output is not balanced.
The rule10 repeated-position transfer lands at `11^2*13*t`; its next rule10
output is `11^3*t`, not balanced.

For the eight exceptional edges, direct symbolic substitution completes
the graph. The rule1 edge ends at `17^3`, which next becomes `17^2`.
The rule3,4,5 edges end at `19^3,23^3,29^3`; each cube takes its listed
rule8,7,6 edge once, after which the profile is `(2,2,1)`, not balanced.
The rule12 edge `8 -> 60` is followed by `60 -> 450`, also profile
`(2,2,1)`. There are no further cases by the iff in §3.

Therefore (4) is the entire recurrent set inside `B`, no other balanced
state feeds into it, and all non-cycle starts exit `B` within three steps.
The exact three-step starts are

    3*17*19^2 = 18411,
    2*19*23^2 = 20102,
    3*11*29^2 = 27753.

All violate the valid-control condition. Thus two steps is the optimal
valid-carrier exit bound, achieved for example by `8 -> 60 -> 450`.

|Universe|Exit at step1|Exit at step2|Exit at step3|States in two-cycles|
|---|---:|---:|---:|---:|
|Active|356|15|3|6 (3 cycles)|
|All outside types|499|20|3|16 (8 cycle types)|
|Active valid|83|9|0|6 (3 cycles)|
|All valid outside types|152|13|0|16 (8 cycle types)|

Here “exit” includes reaching the excluded unit. The six active cycle
states are `{11,13,220,260,550,650}`. The eight outside cycle types are
represented by `(11,13)`, `(220,260)`, `(550,650)`, `(1364,1612)`,
`(8525,10075)`, `(21142,24986)`, `(52855,62465)`, and `(391127,462241)`.
There are infinitely many actual outside-prime cycles, not just eight.

## 6. The exact bridge, and what it does not identify

The shared carrier is the **labelled prime-exponent vector**. The divisor
observation counts marked subsets of its exponent box. The program is an
ordered guarded translation on the same coordinates. This gives a precise
arithmetic bridge while retaining two different evaluations of the carrier.

`(F,U,S)` alone is not a closed deterministic program state. For example,
`2` and `11` both have `(F,U,S)=(0,0,0)`, but `2 -> 15` is unbalanced while
`11 -> 13` is balanced. The enabled rule is different as well. Even the
unlabelled full exponent profile loses the coordinate names needed by the
guards. Projection to the active factor is a true semiconjugacy, but it need
not preserve balance: `8` is balanced and `8*31` is not.

At a prime-power output `2^p` with prime exponent `p`,

    (F,U,S) = (p-1,1,1).

Thus such an output is balanced iff `p=3`. The initial seed `2` is balanced
as a prime, but is not one of the later prime-exponent outputs under that
convention. Balance is neither the program's prime-output test nor an
invariant of its dynamics.

Finally, the literal integers `{7,17,55}` have different arithmetic types:
`7,17` are balanced primes, but `55=5*11` has `(F,U,S)=(2,2,2)` and is not
balanced. The program turns the first two into 1 by different rules; it
does not supply a canonical identification of these three literals with
the three inherited shapes or with three blue symmetry points in a merged
tournament metagraph. Such a comparison needs an explicit map, gauge,
preserved predicate, and sidecar, not equality of cardinalities. Nor is a
ten-prime alphabet automatically the ten-node tournament quotient or a
local translation-invariant cellular automaton rule.

## 6a. Exact finite balance window on the seed 2 orbit

This section uses the independently reconstructed canonical macro words in
`../fractran/FRACTRAN_PRIME_EXPONENT_REPORT.md`, §4. The macro proof is an
explicit dependency; the potential minimum and the complete balanced-visit
consequence are the additional argument here. Let `A,...,N` name rules1–14,
write `C(n,d)=5^n*7^d*13`, and assume `n>=2`, `1<=d<n`. Write `n=qd+r`,
`0<=r<d`. The proved macro prefix is

    W_d^q (AB)^r A, where W_d=(AB)^d J (EF)^d K.

Its endpoint is `2^n*3^r*7^(d-r-1)*17`. If `r>0`, the suffix
`C(DG)^n H(EF)^(r-1)K` ends at `C(n,d-1)`; if `r=0`, the suffix
`I L^n M^(d-1) N(EF)^n K` ends at `C(n+1,n)`.

**Exact token-minimum theorem.** Including both endpoints and all intermediate
states, the minimum of the total prime-exponent count `Omega` in this macro is

    min Omega = n+d-1   if d divides n,
                n+d-2   otherwise.                       (5)

In particular, every state in the block has `Omega>=n`.

**Proof.** At a canonical start, `Omega=n+d+1`. Each `A` or `E` in `W_d`
lowers it by one, and the paired `B` or `F` restores it; `J,K` do not
change it. The partial `(AB)^r A` ends at `Omega=n+d`. For `r>0`, `C`
lowers it to `n+d-1`, each `D` lowers it to `n+d-2`, and `G` restores
it. `H` then raises it, and the final `EF` pairs dip only to `n+d-1`.
The minimum is attained after a `D`. Here `r>0` forces `d>=2`, so it
is at least `n`. For `r=0`, `I` gives `n+d-1`; `L` only increases it;
the `M` moves reduce it only to `2n`; `N` raises it to `2n+2`; and the
final `EF` pairs dip only to `2n+1`. Since `d<n`, the minimum is exactly
the one attained after `I`, namely `n+d-1>=n`. This proves (5).

All balanced integers have `Omega` in `{1,3,4}`. The bootstrap word
`LNEFK` takes seed 2 to `C(2,1)`. The macro endpoints keep the candidate
`n` unchanged or increase it by one, so after first reaching any `C(n,d)`
with `n>=5`, no later state can be balanced. The exact integer/hand-rule
replay reaches `C(5,4)=97540625` at update 145. Direct divisor enumeration
at every state up to and including that boundary gives exactly:

|Update|Balanced state|Update|Balanced state|
|---:|---:|---:|---:|
|0|2|43|380|
|2|825|45|950|
|4|1925|53|1450|
|5|2275|56|850|
|13|156|59|580|
|14|132|62|340|
|16|308|69|8|
|17|364|70|60|

This is a complete all-time list: the prefix is checked exhaustively, and
the infinite tail is excluded by (5). In particular, the existence of
infinitely many prime-power outputs does not give infinitely many balance
events. This conclusion concerns the standard seed 2, not all integer seeds.

`seed_balance_window.py` separately checks the exact minimum formula and
the first-enabled macro words for all 1225 pairs `2<=n<=50`, `1<=d<n`,
covering 267171 transitions; it also compares integer and hand-rule replay
and direct divisor counts in the entire finite seed window. These are
controls on the proof, not a replacement for its quantifiers. The independent
auditor checks the seed window afresh by rational multiplication and bounded-
product divisor enumeration.

## 6b. A square-sum graph observer and the exact last small state

**Dependency status:** the finite integer observations and degree obstructions
below are **FINITE-EXACT**. Their infinite-tail extension is **PROVED** by
(5). The Hamiltonian-cycle conclusion is a **COROLLARY OF THE CITED GERBICZ
EXISTENCE THEOREM**, not a new proof of that theorem. The source is
`05-knowledge/results/session124_gerbicz_reconstruction_20261031/GERBICZ_RECONSTRUCTION.md`
at the pinned repository snapshot; its independently reconstructed classical
certificate proves existence for every order at least32. The original author
statement is linked at [OEIS A090461](https://oeis.org/A090461). No full
small-order classification is needed for the following corollary.

Let `N_t` be the seed 2 integer state after `t` FRACTRAN updates, and let `G_n`
be the simple graph on `{1,...,n}`, with distinct vertices adjacent exactly
when their sum is a square. The finite window gives

    min{N_t : 70<=t<=145} = N_70 = 60,
    {(t,N_t) : 0<=t<=145 and N_t<32}
       = {(0,2),(1,15),(19,4),(20,30),(69,8)}.             (6)

After `C(5,4)` at update145, (5) gives `Omega(N_t)>=5`, hence
`N_t>=2^5=32`. Therefore (6) is also the **complete all-time list** of states
below32, and update69 is the last such time. The bound is not `N_t>=60`:
the already checked prime-power output `N_281=32` is a hostile to that
stronger claim and attains the actual minimum over all `t>=70`.

It follows that

    G_(N_t) has a Hamiltonian cycle
       iff t is not in {0,1,19,20,69}.                    (7)

For the positive direction use the cited existence theorem and `N_t>=32`.
The five negative cases have direct degree obstructions:

|Time|Order|Vertex|All its neighbors|
|---:|---:|---:|---|
|0|2|1|none|
|1|15|8|`{1}`|
|19|4|2|none|
|20|30|18|`{7}`|
|69|8|4|`{5}`|

Every vertex on a Hamiltonian cycle must have degree at least two. Thus these
small explicit certificates prove the negative direction without importing
other finite exclusions. In particular, the graph observer is Hamiltonian
at every time `t>=70`, while it is not at `t=69`. This is not its first
Hamiltonian observation: `N_2=825` already has order above the threshold.

At the same update70 the integer is `60`, the final balanced state from
§6a, and `60=|A_5|=5!/2`. These are verified evaluations meeting at the same
integer. No causal relationship, invariant, group action on `G_60`, or
identification of prime generation with Hamiltonicity follows from them.

### READY-only cellular-automaton consequence

The actual CA's previously computed `../fractran_ca/ready_trace.tsv` gives

    FRACTRAN update69: CA time537, READY integer8;
    FRACTRAN update70: CA time547, READY integer60.

Its proved observer theorem makes READY times correspond to successive
integer updates. Therefore (7) transfers to that observer, and every READY
observation at CA time at least547 has a Hamiltonian square-sum graph. No
claim is made for arbitrary non-READY microconfigurations or for a graph
being stored intrinsically by the CA.

### Literal seeds versus states actually visited

Since every member of `{15,16,17,23,25}` is below32, (6) proves that the
standard forward orbit visits **only15** from that set, exactly at update1.
All five remain valid separate initial integers; being a possible seed is
not being a later state of seed 2. The separate seed analysis is consistent
with this distinction:16 would be a forbidden composite-exponent pure
power on the standard orbit;17 and25 lead to traps; and23 returns to2,
which would make the seed 2 orbit periodic if23 were visited. Those checks
use the proved prime-program semantics, whereas (6) already gives a complete
finite-window-plus-tail proof of the requested intersection.

`cross_observers.py` independently replays the same finite prefix with exact
rational arithmetic, verifies all five neighbor lists, checks the existing
READY trace and its times, and counts the60 even permutations of five letters.
It never instantiates a large graph or searches a longer integer orbit.

## 7. Reproduction, independence, and stopping boundary

Run in this directory:

```sh
python3 verify_balance_dynamics.py > normal.json
python3 -O verify_balance_dynamics.py > optimized.json
cmp normal.json optimized.json
python3 -B seed_balance_window.py > seed_normal.json
python3 -B -O seed_balance_window.py > seed_optimized.json
cmp seed_normal.json seed_optimized.json
python3 independent_audit.py > audit_normal.json
python3 -O independent_audit.py > audit_optimized.json
cmp audit_normal.json audit_optimized.json
python3 -B cross_observers.py > cross_normal.json
python3 -B -O cross_observers.py > cross_optimized.json
cmp cross_normal.json cross_optimized.json
```

The producer reports **5784 always-active checks**, including direct trial-
divisor counts for all 538 inputs and all 538 outputs, independent hand-token
versus integer steps, all-rule reductions, spectator invariance, whole-program
prime renaming, the complete preserving criterion, and the restricted-orbit
bounds. The independent auditor imports no producer code: it constructs the
884 labelled sources separately, uses `fractions.Fraction` multiplication
to select integral outputs, enumerates divisors by bounded products, checks
a different spectator assignment `{101,103,107}`, and reconstructs the
restricted graph from the exact table. With its independent seed-window
replay it reports **8776 checks**. The separate token-potential/window
checker reports **271140 checks**. Normal and `-O` receipts agree byte for
byte for all three scripts. The cross-observer checker adds **305 checks**
and also agrees between normal and optimized modes. Its graph-existence
claim retains the explicit classical dependency in §6b.

The universal reach comes from the inherited balance-shape theorem and the
proof of the outside-prime symmetry reduction, not from extrapolation of a
finite numeric sample. The restricted-orbit theorem is universal only while
balance is retained; the exact 16-visit theorem applies separately to seed 2.
No finiteness, eventual-periodicity, reachability, or universality conclusion
is asserted for arbitrary complete program orbits.
