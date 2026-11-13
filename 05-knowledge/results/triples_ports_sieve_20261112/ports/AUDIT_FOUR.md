# Independent review of four-step H neutrality

**Status: VERIFIED; no unresolved mathematical defect found.** This review
supports promoting [FOUR_RETURNS.md](FOUR_RETURNS.md) to PROVED IN NOTE,
with its finite local lemmas explicitly scoped. It does not claim any
all-length mutation invariance or a count-changing return of length five.

## All-order proof review

Deleting adjacent repeated mutations is valid because mutation is an
involution. Up to renaming by first appearance, the five listed reduced
four-letter shapes exhaust the possibilities. Orders one and two are
covered by cancellation or the alternating two-name shape.

Restriction commutes with mutation when the entire mutated set is
retained. An entry depends only on that set and its two endpoints.
Consequently a failure of the asserted intermediate tournament property
would restrict to at most five vertices for the three-name words. The
finite implication checked on every labelled tournament of orders three,
four and five therefore proves the all-order entry implication. This
argument neither restricts nor estimates H itself.

For `0101`, the inherited pentagon theorem was read directly from incoming
commit `e08896a818`, file
`session291_cluster_ports_20261113/cluster/PROOF.md`, theorem 3. Its explicit
exterior-pair calculation establishes the required all-order identity.
Equivariance gives `D=tau(mu_1(B))`; tournament outputs therefore come from
one native source/sink step followed by relabelling. This preserves H.

For `0102` and `0121`, the verified local implications split the return
into a native step and a return of length at most three, whose independently
reviewed theorem supplies H preservation. There is no assumption that the
whole four-step path remains native.

For `0123`, each one-exterior restriction forces the same side of the
oriented four-block. Internal exterior edges remain unrestricted. The
global input is an extreme transitive four-block, and the displayed order
of mutation moves four successive sources or sinks. All four steps are
native.

For `0120`, one-exterior restrictions force a cyclic triple and a uniform
orientation to each exterior vertex; two-exterior restrictions force all
exterior vertices onto the same side. Those same restrictions determine
every output edge, including exterior–exterior edges. They prove the exact
global operation `C3 -> R` becoming `R -> converse(C3)`, or its reverse.
Hamiltonian paths in an ordinal sum concatenate, so both counts are
`3*H(R)`. The separately enumerated three-vertex exceptions exhaust the
case with no exterior vertices. This is a global factorization proof,
not an inference of global H from induced H values.

## Independent controls

Both commands were rerun successfully:

```text
python -B 05-knowledge/results/triples_ports_sieve_20261112/ports/four_local.py four_local_receipt.json
python -O -B 05-knowledge/results/triples_ports_sieve_20261112/ports/four_local.py four_local_receipt.optimized.json
```

They exhaust the eleven declared local lemma universes and the six
three-vertex exceptional returns. Explicit exceptions remain active under
optimization.

A separate audit exhausts all five canonical words on every labelled
tournament through order five where the names exist: **5,474 matrix/word
cases and 188 tournament outputs**. At every intermediate step it compares
the matrix formula with directed-path addition and cancellation. For every
output tournament it compares H by literal permutation enumeration,
independent of the subset-DP observer. All comparisons passed. This extra
finite check is corroboration; the preceding argument supplies the
unbounded quantifier.

The additional check is reproducible from the repository root:

```python
import itertools, pathlib, sys
sys.path.insert(0, str(pathlib.Path(
    '05-knowledge/results/triples_ports_sieve_20261112/ports').resolve()))
import verify as v
words = [(0,1,0,1), (0,1,0,2), (0,1,2,0), (0,1,2,1), (0,1,2,3)]
def hp(B):
    return sum(all(B[p[i]][p[i+1]] > 0 for i in range(len(B)-1))
               for p in itertools.permutations(range(len(B))))
cases = returns = 0
for word in words:
    for n in range(max(word)+1, 6):
        for mask in range(1 << (n*(n-1)//2)):
            B = v.matrix(n, mask)
            C = D = B
            for k in word:
                C, D = v.mutate_paths(C, k), v.mutate(D, k)
                v.need(C == D, 'independent mutation implementations')
            cases += 1
            if v.tournament(C):
                v.need(hp(B) == hp(C), 'literal permutation H')
                returns += 1
v.need((cases, returns) == (5474, 188), 'complete finite universe')
```

The depth-seven scout is not needed for this proof and was not rerun in
this review. Its bounded no-witness result remains separate from the
proved lower bound of five on any H-changing return length.
