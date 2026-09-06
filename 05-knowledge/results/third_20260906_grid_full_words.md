# Full complement words determine exact marginal grid costs

**Status: PROVED finite optimization reduction; FINITE-EXACT tables;
[independent exact and analytic audit accepted](third_20260906_grid_full_words_audit.md).** This retains the words discarded by the
[divisor-bag compiler](third_20260906_grid_bootstrap.md). Its output is an
exact maximum on a declared arithmetic quotient, not a census of physical
entries. The [refined consumer](third_20260906_grid_refined.md) uses it to
reduce the surviving clock set. LRC(14) remains OPEN.

The closest mechanism is the [full-profile graph theorem](third_20260906_decoder_profile_graph.md),
which used visible subwords to force an actual small sheet-gcd edge. The
new operation applies the same hereditary information to the objective
function itself. The six live objects are clock divisors, marginal costs,
complete complement words, a maximizing state multiset, actual pair ratios,
and physical grid safety. The map keeps every profile involving all six
V labels. It still forgets other physical profiles, actual graph realization,
and the simultaneous placement of the danger arcs. Those are the remaining
coordinates; a maximizing state need not come from a dangerous row.

The canonical hostile is the exact-scalar state
(8,8,9,9,10,60,72), which passes every scalar projection but fails a complete
word. The corrected near miss is to maximize seven costs independently,
then attach an unrelated actual pair. The least-used sidecar is the entire
complement word at each prefix, restored before optimization. A cheap probe
at the former maximal clock23760 lowers its excess bound from252 to168.
That positive signal motivates the complete retained-clock calculation below.

## 1. Exact quotient and the objective

Retain the balanced hypothetical-failure hypotheses and the pinned full
profile supplier of the two linked parent notes. Let P_l be its set of
allowed pairs (core gcd, sorted complement word of length l). At a fixed
clock t define F_t to be all seven-entry multisets D=(d_1,...,d_7) such that

```text
d_i divides t,
gcd(d_1,...,d_7)=1,
for every I with1<=|I|=k<=6:
  (gcd(d_i:i in I), sort(gcd(gcd(d_i:i in I),d_j):j not in I))
       belongs to P_(7-k).
```

Repeated d_i are allowed, as distinct physical speeds can have the same
sheet gcd. The exact maximum is

```text
E_*(t)=max_(D in F_t) [sum_i d_i ceil(t/(7d_i))-t]
      =(1/7) max_(D in F_t) sum_i d_i*((-(t/d_i)) mod7).       (1)
```

Every actual surviving profile row maps into F_t by d_i=gcd(t,u_i).
Consequently its true excess E is at most E_*(t). Conversely, this note
asserts no physical realization for an arbitrary element of F_t. The
all-one multiset is always admissible, so the maximum exists.

## 2. Complete branch-and-bound proof

Each d_i must belong to the42 allowed core values of P_6. Intersect them
with the divisors of t, then order them by decreasing cost in (1), breaking
ties by increasing numerical value. Search nondecreasing indices in this
ordered list; repetitions remain possible. This quotients only permutations
of a multiset, not different state multiplicities.

For a prefix S of length n, each subset I of size k<=min(n,6) has a fixed
core gcd and a visible complement word of length n-k. A necessary condition
for completion is that this visible word occur as a submultiset of an
allowed complete word in P_(7-k) at that core. Precompute all such subwords
from the pinned table. Every valid complete multiset passes every one of
its prefixes, since deleting unchosen entries merely deletes letters from
the complement word. At n=7, the condition is exact full membership.
The global gcd is also checked explicitly.

If a prefix has total cost s and a remaining candidate index j has cost w_j,
then every completion using indices at least j has total cost at most

```text
s+(7-n)*w_j.                                               (2)
```

The cost ordering makes this bound valid even with repeated values. Once
(2) is no larger than the best verified complete cost, neither that branch
nor any later candidate can improve the maximum. An attained complete word
supplies only a lower bound; exhaustion of every remaining potentially
better branch supplies optimality. Every final owner is rechecked against
the literal full profile table, independently of the subword cache.

The same procedure can condition on two fixed states by seeding the prefix
with them and allowing the five remaining indices from the entire ordered
list. It preserves all multiplicities and all prefix requirements. Such a
conditioned maximum must not be substituted for a different selected pair.

## 3. A sharp global residue envelope in the full-word quotient

Every one of the42 allowed core values is coprime to7. Thus if7 divides t,
every summand in (1) is zero and E_*(t)=0. Otherwise the cost of a divisor d
depends only on t modulo7:

```text
w_t(d)=d*((-tau*d^(-1)) mod7),       tau=t mod7.
```

Allowing all42 values gives a uniform upper bound for every clock in the
same residue. The exact maxima, with an attaining full-profile multiset,
are:

| tau | maximum E | Attaining multiset |
|---:|---:|---|
|1|210|64,60,46,36,29,44,17|
|2|270|72,58,64,51,66,45,8|
|3|192|66,45,64,24,58,51,20|
|4|239|72,58,46,51,32,48,27|
|5|197|90,48,33,46,34,32,3|
|6|224|90,48,58,51,44,23,6|

All42 values divide L=5388292800, which is coprime to7. A suitable positive
multiple of L realizes each specified residue and every state divisor in
its row. Therefore these bounds are sharp for the full-word arithmetic
quotient, including the clock divisibility constraint. This is not a
sharpness claim for actual decoder entries or strict failure. The uniform
bound E<=270 improves the previous scalar-cap bound447, while retaining
more of the inherited mathematical object.

## 4. All retained clocks and the boundary that survives

The [baseline compiler](third_20260906_grid_bootstrap.md) supplies exactly
8,301 candidate clocks, all at most23760. The source computes E_*(t) for
every one of them, retaining a maximizing seven-state owner, search size,
and semantic digest of the entire ordered table. No clock outside that
audited baseline can re-enter through this refinement.

The former large-clock controls now have

```text
E_*(23760)=168,     E_*(27360)=164,     E_*(28080)=112.
```

The latter two are explicit additional controls outside the final baseline.
The [refined compiler](third_20260906_grid_refined.md) combines this exact
upper bound with the independent actual-pair marginal bound. It retains
an inequality only when both upper bounds permit the component credit.

A further cheap hostile explains why the new largest scale16704 is not
removed by merely conditioning on the first pair found. Its unconditional
maximum is188. The initially encountered margins (6,12) have conditioned
maximum133, but the alternative actual atlas pair

```text
t=16704, e=4, (p,q)=(3,308),
forced sheet margins=(12,16),
D=(12,16,72,58,64,9,9)
```

passes every full word and has E=188. Its strict actual sum is the inert
prime311, and e*gcd(t/e,p), e*gcd(t/e,q) are exactly12,16. Its individual
component credit is172, so this declared relaxation genuinely retains it.
The first excluded pair is not a universal representative of the atlas.
The missing coordinates are still actual connected realization and the
joint placement of all seven danger sets. Neither this state nor the
failed sufficient bound is claimed to be an unsafe physical entry.

## 5. Reproduction and audit boundary

[Standalone source](../../04-computation/third_20260906_grid_full_words.py)
and [frozen output](third_20260906_grid_full_words.out) use exact integer
arithmetic and import no producer implementation. The input profile bytes
and baseline clock set are pinned separately. Controls include complete
unpruned multiset enumeration for clocks1,...,24; zero septimal costs;
all six global residue cases; former endpoint controls; literal membership
of every maximizing owner; and the excluded/retained conditioned pair.

```bash
python3 -B 04-computation/third_20260906_grid_full_words.py
python3 -B -O 04-computation/third_20260906_grid_full_words.py
```

The independent audit uses different numerical ordering and a bounded
multiplicity-knapsack upper bound, checking complete words after weaker
prefix filters. The [audit](third_20260906_grid_full_words_audit.md) verifies all8,301
maxima, six global maxima, and the boundary cases with1,078,253 always-active
gates over2,980,251 independently ordered search nodes. Its normal and
optimized outputs agree byte for byte; sources and hashes are linked there. The global and finite numerical maxima are load-bearing exact
computations; an all-profile analytic closed formula is not claimed.

Normal and optimized producer outputs agree byte for byte, with1,092,449
always-active gates,127,777 search prefixes over the8,301 retained clocks,
and6,878 unpruned small-clock multisets. Raw LF SHA256:

```text
source 428890519bb54074f35de2b5c2088b5020d48a196b4a3512b69821724e0a652b
output 176d1398ae404ef71a46f7233b6459f2ad3f8c5faae88251eddc8dad171eedfc
maxima ca6b6f562db1fc3632f8b7570b89a16020a981ae8aa130be200dc1bdcb4264ca
```
