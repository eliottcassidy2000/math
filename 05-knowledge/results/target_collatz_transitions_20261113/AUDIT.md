# Scope, independent checks, and reproduction

**PROVED IN NOTES:** the universal quotient-run identities, primitive
guard, finite target inverses, norm-chamber and fixed-partner obstructions,
physical core handoff and clock boundary, and positive inverse-debt
criterion. **FINITE-EXACT:** numerical controls and bank results below.
**OPEN:** complete Chan–Pak coverage, ordinary Collatz convergence, and
universal literal-extra-graph reachability. No Lean claim is made.

## Proof dependencies and priority

The starting snapshot is06847e4c6e6709585bf1b4cfa6f84713081034ea.
Incoming comparisons are read through d8aae7c5ae, including the late
Session307 symmetry constraints. Every new computation
that reads the inherited tournament inventory pins its Git object rather
than reading an unversioned working-tree bank.

- The symmetric mirror, root-star repair and NR maps are inherited;
  [POLARIZED](POLARIZED.md) retains both ordered cross terms instead of
  substituting the same module twice. Its quotient-run table is exhaustive
  by the at-most-two-run argument, not a polynomial fit.
- The existing boundP_M(2)≥2^|M| is explicitly credited to
  [COMPLETE_INVERSE_AUDIT §1](../target_directed_cores_20261112/ports/COMPLETE_INVERSE_AUDIT.md).
  The path-cut boundP_M(2)≥(|M|+1)H(M) is likewise inherited.
- The Π count formula and its interior one-arc transfer are inherited
  from [CLOCKS](../clone_clocks_negative_collatz_20261113/CLOCKS.md).
  New work supplies the source-frame handoff, complete root response,
  target inverse, zero-site response and refill cost.
- The source ternary fork and module qualifications are inherited from
  [CORE_REPLACEMENTS](../prime_ternary_lift_20261113/CORE_REPLACEMENTS.md).
  The new arbitrary-module law has digit4H(K); onlyH(K)=1 gives digit4.
- Old Gold/Pell, Berggren, Stern and square-sum results are used by their
  typed maps in [HISTORY](history/HISTORY.md), not as proofs of tournament
  response realizability. Current incoming work is separated in
  [INCOMING](INCOMING.md).

No scarce theorem ID was reserved. The scoped statements are located by
their descriptive files and section numbers. No old theorem was silently
retracted or re-promoted. Candidate strategies rejected during this
session were not published as established theorems.

## Exact computational universes

| Script | Universe and mechanism | Positive and hostile controls |
|---|---|---|
| response_probe.py | Pinned rooted_growth.json.gz,53,670 rooted records; forget roots and deduplicate ordinary representatives through8 | Ordinary inventory sizes1,1,2,4,12,56,456,6880; exactH,P2 recounted;1,922 distinct response pairs; all ordered pairings tested for e=±1 after the necessary gcd(h,j)=1 filter |
| polarized.py | Complete fixed quotient-run tables;16 literal ordered pairs fromTT1,TT2,C3,Q4 through output order16 | All source/preparation/NR root states; primitive iff guard; ordered asymmetry; direct permutationH for the two order-eight compensation modules; actual22/24-vertex rows checked as tournaments, not recounted by full24-vertex DP |
| selector_head.py | Every pair from the complete response bank withY<39971 | Sharp bounds prove no larger modules can contribute;205 distinctoutputs,34prime,9seven-prime; no claim that this is the full tournament spectrum |
| targets/handoff.py | Complete finite root-response skeleton, five literal retained-module handoffs, three transfers, exact finite target selector | NontransitiveC3 gives digit12; Q4 full primitive state despite gcd(H,q)=5;1081 rejects nonphysical(h,p)=(1,6) but accepts a different chart; retained C3 allows final tail transfer |
| history/check.py | Exact norm, prime uniqueness, Farey and padding controls; independent subset DP for compensation modules | Singleton191; composite604679 has two arithmetic norm candidates;7p obstruction; fixed-pair padding selector accepts only(0,0) for the nontrivial witness |
| collatz/check.py |9,330 short action/halving words,4,917 with paid debt; exact rational prefix maps |6,342 positive paid traces;1,599 endpoint-positive formal fixed points all fail intermediate positivity;150 guarded source lifts; q=0 and even endpoint scope |
| collatz/clock_boundary.py | Complete zero-site run response, literal module graphs and native refill | Γ coefficients(163,180,60); no omitted high runs by indegree bound; radius-two fixed-skeleton probe28single/378double flips distinguishes matchingH from matchingfullresponse |

The ordinary bank's completeness is inherited from its pinned source;
this packet does not rerun a full labelled order-eight tournament census.
The large compensated source follows the proved universal formula, with
independent counts of its two modules. Arithmetic target candidates are
never counted as graphs before their physical module check.

## Independent mathematical review

The constructive-spectrum reviewer checked POLARIZED §§2–5, including
the exact gcd iff, ordered determinant identity, seven-adic phases,
factor inverse, sharp39971 threshold and compensation pair. The history
reviewer independently checked the module observations and the entire
positive-debt proof. The Collatz reviewer independently checked the
norm/Farey proofs and source handoff. The root reviewed each lane's
universal proof, endpoint domain and failure boundary.

The following distinctions were explicitly checked before integration:

- “root name retained” is not “root star retained” in the core handoff;
- three incidence classes can mean3|K| individual old arcs;
- output primitivity, unit readiness and return-chart membership differ;
- gcd(H,q)=1 is sufficient for full primitivity, not necessary;
- a prime cofactor family is not established by a family divisible by7;
- the clock's last site must be treated by a separate graph, not by
  empty-module substitution into a nonempty-module identity;
- a Collatz affine word needs positive intermediate values and the
  appropriate integer guards; ordinary genuine inverses differ from the
  extra odd inverse; an even endpoint may be a partial final odd return.

## Reproduction

From the repository root, run these in the listed dependency order:

```bash
python3 -B 05-knowledge/results/target_collatz_transitions_20261113/response_probe.py
python3 -B 05-knowledge/results/target_collatz_transitions_20261113/polarized.py
python3 -B 05-knowledge/results/target_collatz_transitions_20261113/selector_head.py
python3 -B 05-knowledge/results/target_collatz_transitions_20261113/targets/handoff.py
python3 -B 05-knowledge/results/target_collatz_transitions_20261113/history/check.py
python3 -B 05-knowledge/results/target_collatz_transitions_20261113/collatz/check.py
python3 -B 05-knowledge/results/target_collatz_transitions_20261113/collatz/clock_boundary.py
```

All seven scripts passed both normal and optimized Python execution.
Receipts are byte-identical across those modes; see[validation.json](validation.json)
for their hashes and captured results. Checks use explicit failures and
are not dependent on Python assertions remaining enabled. Root SHA256SUMS
covers the finished packet; the history lane additionally has its own
manifest. New relative links are checked against working files or the
Git object tree, including files omitted by the sparse checkout.

Repository integration: `agents/check_docs.py` passed after adding the
results-index and EIGHT-NINE routes. Maintained header/status lines remain
at their existing positions. `git diff --check` also passed.
