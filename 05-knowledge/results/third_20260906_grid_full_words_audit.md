# Independent audit of exact full-word grid excess maxima

**ACCEPTED: analytic reduction and pruning, all 8,301 exact maxima, six
global residue optima, and both conditioned boundary controls independently
verified. Normal and optimized outputs agree byte for byte.** This audits
[the full-word optimization](third_20260906_grid_full_words.md), including
its all-clock residue envelope, all 8,301 retained-clock maxima, and the
conditioned stopping boundary at 16,704. Its domain is the declared
seven-state profile quotient, without a physical realization assumption.

## 1. Reduction and primary search audit

For every surviving actual balanced entry, the seven states
`d_i=gcd(t,u_i)` divide t and inherit the full complement words involving
all six smaller-component labels. Thus they belong to the declared finite
set F_t. For each divisor d, the exact identity

```text
7d ceil(t/(7d))-t = d*((-(t/d)) mod7)
```

holds. Summing seven such identities gives the stated objective, including
the factor seven. All-one states are admissible, so the optimization is
nonempty. Repeated states are retained and the global gcd-one condition is
explicit. Maximizing this necessary arithmetic image gives a valid upper
bound for actual excess; it does not prove that a maximizing row is
physically realized or dangerous.

The producer's subword pruning is hereditary: a prefix complement word is
obtained from a completed word by deleting letters. Its core gcd is unchanged
because its selected subset is unchanged. The descending-cost bound is
valid with repetitions: every remaining choice has cost at most the first
available remaining cost. Pruning equal bounds is legitimate because an
already checked owner attains that value and the goal is the maximum value,
not enumeration of all maximizers. Sorting the index sequence quotients only
permutation symmetry. Seeding a fixed pair and allowing the remaining
indices from the full domain correctly retains all conditioned completions.

## 2. Independent optimization route

The [audit source](../../04-computation/third_20260906_grid_full_words_audit.py)
imports no producer implementation. It first checks all 8,301 supplied owners
against literal full profile membership and the direct ceiling objective.
Their scores are then legitimate incumbent **lower bounds**, independent of
any producer assertion of optimality.

The audit orders divisors by increasing numerical value. Its partial filter
is deliberately weaker than the producer's: individual states retain visible
subwords of their length-six words, while larger selected subsets retain
only their exact allowed scalar gcd projections. These necessary conditions
cannot discard any valid full row. At every complete search leaf, the audit
tests all full words literally.

The independent upper bound is a finite bounded-multiplicity knapsack.
If d occurs r times in a valid row, the gcd of those r positions is d, so
r cannot exceed its maximum permitted multiplicity in the exact scalar
projections. The value one has capacity seven. The audit precomputes the
best possible score for each numerical suffix and each remaining length,
using these capacities and ignoring all joint constraints. This supplies
an upper bound. Already used copies and fixed-pair copies are subtracted
from their respective capacities. The search explores every numerically
ordered completion whose relaxed upper bound exceeds the verified incumbent.
Consequently exhausting the search proves optimality independently of the
producer's cost order, upper bound, and stronger prefix filter.

All 42 allowed states are coprime to seven. Therefore both the admissible
divisor domain and all its costs are determined by

```text
(allowed states dividing t, t mod7).
```

There are 1,440 such patterns in the pinned 8,301-clock input. The audit
optimizes each pattern once, verifies that every claimed value in that
pattern agrees, and transfers the result using this exact identity. It
never identifies clocks merely because their previously claimed maxima
happen to agree.

The input clock set has semantic SHA256
`a25d83f0eeb630bb82e84cdfac4e3cf7312f892f6c426d6affd5239a064e4b58`.
The ordered `[t,E_*(t)]` table has semantic SHA256
`ca6b6f562db1fc3632f8b7570b89a16020a981ae8aa130be200dc1bdcb4264ca`.
The full inherited profile bytes are separately pinned to
`935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`.

## 3. Global residue envelope and stopping boundary

The universal residue bound follows from maximizing over all 42 states,
which contains every divisor-restricted domain. Their lcm is 5,388,292,800,
coprime to seven. Suitable positive multiples realize all states and any
nonzero residue, so the displayed maximizing rows prove sharpness in the
arithmetic quotient. For clocks divisible by seven all costs vanish.
The independent optimizer separately checks all six claimed global maxima
`(210,270,192,239,197,224)` and their owners; the uniform bound is therefore
270 in this quotient. No actual-entry sharpness is inferred.

At t=16,704, fixing states (6,12) gives maximum 133, while fixing (12,16)
gives 188. Both are separate complete conditioned searches. The latter
value is attained by

```text
(12,16,72,58,64,9,9).
```

For the actual atlas ratio (3,308) and e=4, the exact sheet margins are
`(e*gcd(t/e,3),e*gcd(t/e,308))=(12,16)`. Its sum 311 is an allowed
inert prime. An independent rational threshold-wall sweep reconstructs
every open component of the primitive two-runner bad overlap. Applying
the strict open-interval grid ceiling to these literal lengths gives
component credit 172. The first margins (6,12), with e=6 and ratio (22,205),
have credit 174. Thus the first conditioned maximum is below its credit,
while the alternative full-word row still exceeds its credit. This is a
valid stopping reason for the sufficient criterion; neither row is asserted
to be an unsafe actual physical entry.

## 4. Finite scope and reproduction

The audit additionally enumerates every admissible seven-entry multiset
without pruning for clocks 1 through 32; checks the extra maxima at 23,760,
27,360, and 28,080; and includes zero-cost clocks 7,49,343 and seven times
the full-state lcm. The full 8,301-clock table, six global residues, two
conditioned optima, and literal pair component counts are separate declared
universes. Their exact values do not assert a census of physical decoder
rows.

[The frozen audit output](third_20260906_grid_full_words_audit.out) records
one independently optimized representative per exact pattern and all global
and conditioned controls. Reproduce with:

```bash
python3 -B 04-computation/third_20260906_grid_full_words_audit.py
python3 -B -O 04-computation/third_20260906_grid_full_words_audit.py
```

No formula correction was required. The independent audit explores 2,980,251
search nodes across the 1,440 clock patterns and also checks 11,311 completely
unpruned small-clock multisets. Every table optimum, every global residue
optimum, and both conditioned optima agree. Literal component sweeps give
33 components and credit 174 for (22,205) at e=6, and 45 components and
credit 172 for (3,308) at e=4.

There are 1,078,253 always-active exact gates. Raw LF SHA256 values:

```text
source 9de11127743456560770e0f98ee854a783f9186b93dfafdcbb3682f5364550ca
output 5bf9f1042176d246719a35854f93134f5a412d059964f65ecd6c6bff7c54e625
semantic 4a00f34f0aab69b57e4086b1418f2ac08394e98e049e9ea0085713cd1528d0f5
pattern_trace adf9b2568cf9045602dedffeea3595f452369c847bfb9bc263c44e8bb75cd01e
```
