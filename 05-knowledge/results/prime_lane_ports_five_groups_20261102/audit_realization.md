# Independent audit: actual block sorting and strong dyadic lane exchanges

**Verdict: PASS.** The R0 arbitrary-block monotonicity, R1 uniform-endpoint
commutator, R2 strong dyadic families, and stated finite lane exchanges in
[realization.md](realization.md) are valid in their declared scopes. This
review imports no producer implementation. The claims do not prove every
intermediate count, every prime, or every required 7p is realizable.

## Proof audit

**R0 works for arbitrary actual blocks, including asymmetric blocks.**
The profile entries are literal counts, not formal variables asserted to be
realizable. With `D=E-H-P`, D counts covers whose nonempty first path ends
inside the marked cut. Likewise `E-H-Q` counts covers whose nonempty second
path starts outside it. Both are nonnegative. Removing the empty-first
covers from Q leaves exactly `Q-s` proper covers. For every cut vertex v,
the second path `(v)` and any Hamilton path of `B-v` give at least one such
cover. The assumption `|B|>=2` makes the other path nonempty. Thus
`Q-s>=|S_B|`. The analogous singleton first path outside the cut gives
`P-t>=|B|-|S_B|`. This checks both strictness witnesses and the boundary.

Direct multiplication verifies the two commutators in R0, including their
negative g-coordinate entries. Their signs cannot be inferred merely from
entrywise positivity. The common suffix handles the negative coordinate
correctly: its first row obeys `alpha>0`, `beta>=0`, and
`0<=gamma<=alpha`. The last bound follows because the incoming g coordinate
can contribute only at the suffix's first block; that block then resets g,
while q and h never depend on the old g. Using actual-state `0<=g<=h`, the
lower bounds in the proof follow with the stated signs. In the zero-cut / empty
block-cut and one-cut / full block-cut cases, the **entire** relevant
commutator vanishes. Otherwise the inequalities are strict. No assumption
that the prefix or suffix is vertex-transitive is used.

**R1 is the vertex-transitive specialization.** The automorphism action makes
Hamilton endpoints and endpoints of the relevant ordered cover slots
uniform. Hence `u=H/n` and `d=(E-H)/n` are integers, and the stated P,Q,s,t
expressions follow. Actual proper two-path covers imply `d>u` for a
nontrivial block. The dependence of `f=H(B+x)` on the detailed cut need not
be uniform; it cancels from the commutator. The nilpotent-square identity is
an identity of signed difference operators, not a positive graph operation.
The producer clarified the phrase “transitive case” to
“vertex-transitive case” after review.

**R2 uses the actual prefix clock `r>=1`.** With r transitive prefix
vertices beaten by x, the state is `(1,2^r-1,1,1)`. Applying the two actual
block orders gives

| Cut size a | U | V |
|---|---:|---:|
| 1 | `15*2^r-7` | `18*2^r-7` |
| 2 | `9*2^r-1` | `15*2^r-1` |

Each child has `r+5` vertices. All apex ports and internal block edges stay
fixed; only the three edges between z and the triangle reverse. Since
r>=1, x reaches the first parent SCC. All parent SCCs reach the last, and
the last reaches x, proving strongness. Parent H remains 3, and its
commutative ordered-cover inventory is unchanged. The apex response is
exactly the information that inventory alone omits.

The `r>=1` condition is essential. At r=0 the actual pairs are `(5,9)`
for a=1 and `(5,13)` for a=2; the displayed formulas do not apply. This
boundary was explicitly tested rather than silently extending the prefix
state to the empty parent.

**Lane arithmetic is correct.** The literal constructions give
`35=7*5 <-> 59`, `71 <-> 119=7*17`, and `287=7*41 <-> 479`.
Trial division independently verifies each named prime. The congruence
branches and exact power-of-two inverse describe only these four families;
they have O(log X) values up to X and omit the already attained prime 13.
Neither the relocation monotonicity nor primality of some examples implies
interval filling or uniform prime production.

## Independent literal computation

[audit_realization.py](audit_realization.py) constructs named directed graphs
from scratch, uses adjacency-constrained DFS to enumerate actual Hamiltonian
vertex orders, and never imports the producer's transfer or DP. It checks:

- **1,098** complete block/cut profiles: every labelled tournament of orders
  1 through 4 with every apex cut. Ordered two-slot covers are enumerated
  as actual pairs of vertex words. Both commutators and strict proper-cover
  inequalities are checked from these counts.
- **1,584** contextual graph swaps: all labelled blocks of orders 2 and 3,
  plus the explicit asymmetric strong four-vertex block with arc mask 4;
  every block cut, both singleton cuts, and every empty or one-singleton
  prefix/suffix with its actual cut. This includes **396** empty/full cut
  controls. The graphs share literal vertex names, and their symmetric arc
  difference is exactly the reversed singleton/block cut.
- All four families through **r=5**: counts and strong connectivity checked
  on actual graphs. Full enumeration of **all** vertex permutations through
  r=3 independently agrees with the valid-word DFS.
- The complete ordered nonempty path-cover profile of each parent through
  r=3, counted by permutations and cut positions. A failed consecutive arc
  must be cut; choosing further cuts gives the exact number of covers.
  Both parent orders have identical profiles and H=3.
- The r=0 failure boundary and the displayed prime factorizations.

The finite profile/context ranges are controls; the all-context statements
come from the proof above. The named labels matter: independently renumbering
vertices by their new block position would not by itself verify that the
same physical apex ports were retained. This checker holds those names
fixed throughout.

Reproduce from the repository root:

```bash
python3 -I -B 05-knowledge/results/prime_lane_ports_five_groups_20261102/audit_realization.py
python3 -I -B -O 05-knowledge/results/prime_lane_ports_five_groups_20261102/audit_realization.py
```

Both modes produce the same [audit_realization.json](audit_realization.json).
Checks use explicit exceptions and remain active under optimized Python.
No correctness objection remains for these realization statements.
