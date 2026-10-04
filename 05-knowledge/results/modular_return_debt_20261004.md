# First-return budgets have symbolic exit tails and an exact debt partition

2026-10-04. **PROVED** the cut-set, final-exponent-tail, density, and finite
observer obstruction theorems below. **FINITE-EXACT** the complete 223/233
cut compilations and declared independent controls. A descending return
supplies a smaller dependency; it is not by itself a certificate to 1.
Universal Collatz certificate coverage remains **OPEN**.

Artifacts: [script](../../04-computation/experiments/modular_return_debt_20261004.py)
and [output](modular_return_debt_20261004.out).

## 1. Inheritance and the concrete improvement

The incoming [triplet/CRT note, section 5](triplet_crt_223_233_20261004.md)
exhausts all 1023 positive valuation words of total cost at most ten. Its
four 233-return words strictly descend at every positive guarded
233-multiple. Their exact relative density among odd 233-multiples is
`1/128`; the least source is 55221. Its compulsory hostile is the expanding
first return `707155 -> 755153` with word `(1,3,1,1,1,3,1)` and cost eleven.

**THM-4489**,
[Fermat sextics and 223 returns](../../01-canon/theorems/THM-4489-fermat-sextic-223-collatz-return-budget.md),
gives the corresponding budget-eight theorem and cost-nine failure at223.
The [ordered-carry theorem](geometry_collatz_drift_carriers_20261004.md),
section 3, already explains why affine word composition retains more than
letter counts. The [unfinished-observation union](adaptive_observation_union_20261004.md)
already retains checked prefixes and closes only from the grounded root.
The [mod19 observer](mod19_recursive_observers_20261004.md) and
[shared-depth note](mod19_resonance_depth_20261004.md) separate finite phase
compatibility from exact source membership.

The new mechanism is to compile the **first exit beyond the cost budget**
as a whole final-exponent tail. Five of the 233 overflow cells immediately
return to the marked section. Four descend throughout; the fifth splits
exactly at final exponent two. This adds relative density `9/2048` to the
old `1/128`, giving **25/2048**. A newly admitted source is

    47765 -> 2239 -> 3359 -> 5039 -> 7559 -> 11339 -> 17009,
    word (6,1,1,1,1,1), total halving cost 11.

This is the least source of the extended 233 rule. It is smaller than the
previous 55221, but still outside the old source benchmark `n<=10000`.
These are exact cylinder densities, not sampled trajectory densities.

Closest mechanism: a final exponent does not enter a word's carry.
Canonical hostile: the next expanding return must remain debt.
Corrected near miss: a budget failure is not a failed trajectory.
Least-used sidecar: the checked crossing edge at the first overflow.
The live board is **source / ordered carry / marked return / cost cut /
terminal tail / finite phase / root boundary**. The anchor is source-aware
return certification, the niche is finite cut algebra, and the wildcard
is the obstruction to replacing exact cost by joined modular clocks.

## 2. One congruence represents infinitely many final exponents

Let v be a possibly empty positive valuation word with length j, cost S,
and affine data

    F_v(n)=(P_v*n+B_v)/Q,   P_v=3^j, Q=2^S.

Appending one exponent k gives

    F_(v,k)(n)=(P*n+C)/(Q*2^k),
    P=3*P_v, C=3*B_v+Q.                              (1)

Crucially C is independent of k. Fix K>=1. The exact family having prefix
v and actual next exponent k>=K has original-source cylinder

    P*n+C = 0 mod Q*2^K.                              (2)

**Proof of both directions.** Actual legality immediately implies (2).
Conversely, reduce (2) modulo 2Q. Since
`P*n+C=3(P_v*n+B_v)+Q`, it forces
`P_v*n+B_v=Q mod2Q`. This is precisely the exact prefix-v cylinder: its
endpoint is odd, and backward division recovers the prescribed preceding
valuations. At that endpoint, equation (2) makes the next valuation at
least K. For empty v it is simply `v2(3n+1)>=K`.

The exact prefix-cylinder fact has an elementary forward induction. For
one letter a it is `3n+1=2^a mod2^(a+1)`. For at least two letters, split
`B_v=3^(j-1)+2^a_1*B_tail`, where B_tail is odd. Reducing
`P_v*n+B_v=Q mod2Q` modulo `2^(a_1+1)` forces
`3n+1=2^a_1 mod2^(a_1+1)`. The first step is therefore an odd integer with
exact valuation a_1. Divide out that step and apply the remaining-word
congruence inductively. This verifies integrality as well as parity.

Equation (2) is one odd residue modulo `Q*2^K`, since P is odd and C is
odd. It represents the union of all exact cylinders `(v,k)`, k>=K, without
enumerating infinitely many k. On a matching positive n,

    k=v2(P*n+C)-S,
    endpoint=(P*n+C)/2^(S+k).                          (3)

The inequality k>=K is not permission to replace k by K in (3). The
reported actual endpoints in the output retain any additional halvings.

### Marked first returns and all-height descent

Let m>=3 be odd; primality is unnecessary for this part. Restrict to positive
odd sources divisible by m. A word prefix u reaches another m-multiple iff
`m|B_u`, because every power of two is invertible modulo m. Thus (v,k) is a
first m-return for every k>=K exactly when

    m does not divide B_(v[:i]) for 1<=i<=j,
    m divides C.                                      (4)

This criterion is independent of the final exponent. Necessity and
sufficiency hold on every legal source. The conditions are not vacuous:
CRT supplies positive odd m-multiples in every cylinder (2).

First-hit root legality is automatic for a matched first-return tail:
if an earlier vertex were 1, all its later accelerated vertices would
remain 1, whereas the final vertex is divisible by m>=3. This does not
authorize root self-loops in unrelated debt prefixes; section 3 retains
the runtime root stop.

A useful sufficient uniform descent condition for the whole tail is

    m*(Q*2^K-P) > C.                                  (5)

Every source is at least m, and increasing the actual k only increases
the denominator. Hence (5) implies `(P*n+C)/(Q*2^k)<n` for all matched
sources. A suitable K always exists. The code's `uniform_descent_tail`
finds the first K satisfying this **sufficient bound**. In general it is
not a necessary condition: the least source admitted by the exact guard
can be much larger than m. In the two audited compilations below, every
discarded exponent has `Q*2^k<P`, so those discarded branches genuinely
expand at every positive source; the displayed cuts are sharp there.

### Density is an arithmetic progression calculation

Adding `n=0 mod m` to (2) gives one residue modulo `m*Q*2^K`. Its relative
natural density among odd m-multiples is

    2/(Q*2^K) = 2^(1-S-K).                            (6)

This uses coprime CRT and a finite periodic set. It does not rely on
countable additivity of natural density for the individual exact-k cells.

## 3. A bounded-cost tree has an exact overflow partition

Fix H>=0. Starting from an m-multiple, stop a branch at its first marked
return if its accumulated cost is at most H. At every remaining prefix v
of cost S<=H, the first overflowing next valuation satisfies

    k >= K=H-S+1.                                    (7)

Store that whole branch as a tail and continue only the next exponents
`1<=k<K` within budget. No modular return is assumed to occur eventually.

**PROVED cut-set theorem.** The bounded first-return cylinders and the
overflow tails are pairwise disjoint and partition the odd original-source
residues modulo `2^(H+1)`. Every overflow tail occupies exactly one such
residue and has relative mass `2^-H`. A return word of cost A<=H occupies
`2^(H-A)` cells and has relative mass `2^-A`.

**Proof.** At a prefix of cost S, actual next valuations split into the
disjoint exact cases `k=1,...,H-S` and their complementary tail (7).
Equations (1)--(2) show these are exact original-source cylinders.
Stopping at a first marked return prevents any descendant from also being
a leaf. Every non-stopped branch increases cost by at least one, so the
construction terminates after at most H within-budget steps. The tail's
modulus is `2^(S+K)=2^(H+1)`; the return word's exact modulus
`2^(A+1)` divides it. These facts prove both the partition and the counts.

For this arithmetic partition only, use the formal accelerated extension
`U(1)=1`. It is useful for residue classification but is not a first-hit
certificate convention. The source-aware `observe_to_budget` stops at
the **first** 1 and never evaluates an outgoing root edge. A formal
overflow class may therefore contain an input whose runtime record has
already closed at 1. A crossing edge may also itself end at 1, which is
recorded explicitly. Neither case is an unresolved convergence assertion.

The retained debt record contains

    (original source; exact prefix; checked path; cost used;
     current frontier; overflow tail;
     actual crossing exponent and checked crossing endpoint).

The crossing edge is evaluated once and kept even though its halving cost
exceeds H. Budget here measures halvings, not Python calls or odd edges.
A later computation can reuse that exact edge. If the tail satisfies
(4)--(5), it immediately provides a descent witness. Otherwise its actual
prefix and exit remain usable observations, without an invented smaller
dependency. A first return and a first descent are different predicates:
none of these rules requires earlier vertices to exceed the source.

## 4. Complete 233 and 223 exit compilations

For m=233,H=10, the inherited four return leaves occupy 8 of the 1024
odd dyadic cells. There are **1016** remaining overflow leaves. Exactly
five satisfy the next-step modular return condition (4):

| Prefix v | Overflow k starts at | Descending tail k starts at | Carry C | Added relative density |
|---|---:|---:|---:|---:|
| (5,2) | 4 | 4 | 233 | 1/1024 |
| (2,1,1,2,3) | 2 | 2 | 1631 | 1/1024 |
| (1,4,5) | 1 | 1 | 1165 | 1/1024 |
| (6,1,1,1,1) | 1 | 1 | 13747 | 1/1024 |
| (1,3,1,1,1,3) | 1 | 2 | 5359 | 1/2048 |

These rows are disjoint because they refine different leaves of one cut.
Combining the first two with their old bounded siblings gives the simpler
all-family statements `(5,2,k)` and `(2,1,1,2,3,k)` for every k>=1. Their
combined mass alone is `2^-7+2^-9=5/512`. Including all five terminal exit
families gives

    bounded mass      1/128   =16/2048,
    added tail mass           9/2048,
    total                    25/2048.                   (8)

The exact complementary mass `2023/2048` is the domain not certified by
**this modular-return rule**. It includes branches handled by unrelated
rules or by reaching 1; it is not a density of nonconvergent or necessarily
pending sources. The complete direct census over the 2048 odd
233-multiples in period 954368 finds exactly 25 admitted sources.

The expanding hostile is retained without reinterpretation:

    original source 707155,
    within-budget prefix (1,3,1,1,1,3), cost10,
    current frontier 503435,
    first overflowing edge: exponent1, endpoint755153.

Its first return increases the original source. Changing the final exponent
to any k>=2 gives a different exact cylinder, and every source there
descends. This is a branch split, not two legal choices at 707155.

For m=223,H=8, the two bounded return leaves have mass `3/256` and leave
253 overflow cells. Exactly three are immediately returning tails:

| Prefix v | Overflow k starts at | Descending tail k starts at | Carry C | Added relative density |
|---|---:|---:|---:|---:|
| (2,3,1) | 3 | 3 | 223 | 1/256 |
| (1,1,2,1,2,1) | 1 | 4 | 2899 | 1/2048 |
| (1,1,1,1,2,1,1) | 1 | 5 | 6913 | 1/4096 |

The total becomes `3/256+19/4096=67/4096`. A complete direct census of
4096 odd 223-multiples in period 1826816 finds exactly 67 admitted sources.
This second instance checks the general mechanism; no larger prime census
or maximality claim is made.

### Added return coverage is not automatically added stopping coverage

Every one of the five 233 completion families already has a uniformly
smaller vertex in its first one or two odd steps. Prefixes starting with
5, 2, or 6 descend immediately for n>1. The other two start `(1,4)` or
`(1,3)`, with maps `(9n+5)/32` and `(9n+5)/16`; both descend for every
positive odd n>1. Even the expanding first 233-return from 707155 had
already descended at its second odd step.

Consequently (8) is strictly increased **marked-return-rule coverage**,
not proved new generic descent coverage or a measured selector improvement.
This is a decisive application boundary of the 233 example. Its useful
new object is the exact overflow cut and reusable crossing-edge record.

The two longer 223 prefixes above behave differently: every nonempty
proper-prefix multiplier is greater than one. Their refined final tails
therefore give a genuine first descent on odd step seven and eight,
respectively. This is proved by the prefix inequalities and (5), and
independently checked in the stated finite controls. No disjointness from
all previously stored generic descent rules is claimed.

In fact the 223 source filter can be removed from these two first-descent
families. They are precisely the corresponding objects of the concurrent
[frontier-family compiler](frontier_family_compiler_20261004.md):

    nominal word (1,1,2,1,2,1,4):
      source=1703+4096*t, endpoint=oddpart(910+2187*t);
    nominal word (1,1,1,1,2,1,1,5):
      source=6815+8192*t, endpoint=oddpart(5459+6561*t);
    t>=0 in both cases.                               (9)

The actual final exponent is the nominal last exponent plus the valuation
of the displayed endpoint argument. The least sources satisfy
`1703>2899/1909` and `6815>6913/1631`, so even the nominal affine endpoints
are strictly below their sources; all extra divisions preserve that
inequality. Every proper prefix grows. This proves first descent on the
whole dyadic classes, independently of divisibility by223. The marked
section merely selects which of these trajectories also return to223.
There is no claim of new coverage beyond the existing first-crossing atlas.
The script replays t=0..63 independently for both families, including the
seed1703 whose actual final exponent is5 although the nominal exponent is4.

## 5. Ordered carries and joined clocks have different jobs

The multiset of exponents determines P and Q but not C. The minimal
comparison needed here is

    (5,2,1): P=27,Q=256,B=233,
    (2,5,1): P=27,Q=256,B=149.

Only the first returns a 233-multiple to the marked section. Both have
legal positive 233-multiple sources by CRT, in separate exact cylinders.
Thus a commutative inventory cannot decide the modular return predicate.
The inherited full affine-word decoder explains the strongest survivor:
retaining the ordered carry together with P,Q retains the whole word.

Finite residue observations lose a different coordinate: the unbounded
final exponent. There is an exact all-precision obstruction using the
expanding 233 prefix above, whose terminal map is

    F_k(n)=(2187*n+5359)/(1024*2^k).

Let M>=3 be any odd observation modulus coprime to233, and put
`d=ord_M(2)`. The formal maps F_1 and F_(1+d) agree on **every** residue
modulo M. Their carries, lengths, and exponent-clock residues agree.
For any prescribed source residue x modulo M, CRT supplies separate
positive sources satisfying

    n=0 mod233,  n=x mod M,
    the exact cylinder for the chosen word.

The k=1 source expands, while the k=1+d source descends, since d>=1 and
all k>=2 tails satisfy the uniform bound. Both are first 233-returns.
Their observed input and endpoint residues agree. No finite M can recover
the missing real denominator magnitude from those observations.

In particular this holds jointly for **every** finite ternary precision
and mod19 precision: choose `M=3^a*19^b`, a,b>=1. Its exponent period is

    lcm(2*3^(a-1),18*19^(b-1)).                       (10)

The script constructs exact sources for all a=1..4,b=1..2, including
last exponents over one thousand. It compares entire affine maps, not
just one accidental endpoint equality.

This does not make modular pruning useless. Given an independently
specified target residue, an incompatible affine endpoint can reject a
candidate; condition (4) itself is precisely such a filter. Fixed-hub
inverse address constraints retain their shared ternary coordinate from
the resonance theorem. But every nonempty dyadic source cylinder meets
every residue modulo an odd modulus coprime to its marked section. Extra
phase labels alone cannot replace its exact source guard or its actual
last exponent. A constructed source with a requested phase is also not a
certificate for a different, independently supplied integer of that phase.

## 6. Reproduction, interfaces, and limits

```text
python -X utf8 -B 04-computation/experiments/modular_return_debt_20261004.py
python -O -X utf8 -B 04-computation/experiments/modular_return_debt_20261004.py
```

The implementation is self-contained. `compile_cut(m,H)` returns bounded
first-return words and overflow tails. `uniform_descent_tail` supplies a
sufficient symbolic cutoff. `apply_return_tail` checks the original source
and replays its entire exact first-return word. `observe_to_budget`
retains prefix and crossing edge, stopping on the first root. No API marks
a smaller endpoint as already certified.

The finite universe is exactly `(m,H)=(223,8),(233,10)`: all 255/1023
positive compositions through those budgets, all 256/1024 cut cells,
all 4096/2048 odd source multiples in the refined periods, and six exact
last exponents at three source heights for each of the eight terminal
families (144 controls), plus128 unconditional frontier-family sources.
Recursive and cut-set word generators agree;
direct carry sums agree with affine recurrence; literal trajectories
independently agree with every counted source cylinder. All assertions
remain active under optimization.

Hostiles include the incoming expanding return, six invalid source/budget
inputs, ordered-carry loss, and eight exact pairs of finite-observer aliases
with opposite drift. Root controls `3->5->1` and the zero-budget crossing
`5->1` verify that neither a root self-loop nor a discarded completed edge
is silently introduced. In the 233 census one crossing edge reaches the
root; the residual modular-rule mass deliberately does not call it pending.

| Source -> target | Preserved predicate | Lost coordinate / required sidecar |
|---|---|---|
| Guarded words -> a final-exponent tail | Exact prefix and minimum final valuation | Actual k; recover it from the supplied source |
| Cost tree -> first-overflow cut | Original-source partition and exact natural density | Work after the cut; retain frontier and checked exit |
| Ordered word -> P,Q,B | Exact affine action and decodable word | Nothing at word level; integer legality still needs the source guard |
| Word/source -> finite phases | Modular input, output and compatible clock | Actual denominator and source magnitude |
| Descending return -> smaller dependency | Checked path to a strictly smaller integer | A home certificate for that dependency remains required |

The positive result is a reusable all-family completion of selected budget
exits. It neither proves all exits return nor supplies a uniform route from
every remaining debt state to a certified smaller dependency.
