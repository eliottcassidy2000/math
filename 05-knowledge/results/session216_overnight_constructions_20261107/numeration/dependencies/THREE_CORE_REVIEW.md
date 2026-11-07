# Independent audit: three-core positive arithmetic constructor

Status: mathematical proof and independent computations PASS; input freeze receipt recorded separately.

## Exact scope

The construction uses three fixed cut-labeled strong cores with outgoing rows

- H=5: [8,1,3,6], order4;
- H=9: [2,4,9,19,7], order5;
- H=9: [2,8,3,21,7], order5.

A transitive prefix and suffix surround the core by order joins. Every padding vertex has a freely specified cut bit. Adjoin v and x, with L⇒v→x, and x→u precisely for cut vertices u in L. The strong subfamily requires the first prefix vertex to be in the cut if a prefix exists, or a nonempty core cut otherwise.

All odd targets 23..999999 except47,63 have retained recipes with total padding<=16. This is499987 target values, and the largest retained tournament order is22. Independent standalone tournaments of orders7 and8 attain47 and63. Thus the combined finite bank covers all odd targets in[23,1000000], but the two exceptions are NOT outputs of the strong three-core template. No all-target or all-length completeness follows.

## Independent derivation

For a cut-labeled nonempty L, let h count spanning paths; q count ordered two-path covers allowing either path empty; α count such covers whose second path is nonempty and starts in the cut; γ count those whose first path is nonempty and ends in the cut; τ count covers satisfying both restrictions; and s,e count spanning paths starting/ending in the cut.

In particular τ<=α by actual inclusion of cover sets, not an estimate.

Prepending a source with cut bit ε gives

    (h,q,α,γ,τ,s,e)
       -> (h,2q,α+εq,2γ+εh,τ+ε(γ+s),εh,e).

Appending a sink gives

    (h,q,α,γ,τ,s,e)
       -> (h,2q,2α+εh,γ+εq,τ+ε(α+e),s,εh).

All right-hand sides use old coordinates. To prove these, place the new vertex in either ordered path. A new source changes the first-vertex condition on its path, and a new sink changes the last-vertex condition. A singleton first/second path accounts for the h,s,e correction terms. This derivation does not assume the numerical state table or a recurrence fitted to observations.

The singleton-root star closure has H=2α+h−τ. This also follows directly without module machinery: remove x from a Hamilton path. If v precedes x and the following portion is nonempty, the two L-portions contribute α; if x is last there are h paths. If v follows x, the first L-portion must be empty or end outside the cut, so the contribution is α−τ. The sum is2α+h−τ.

The strongness test is exact. Every L vertex reaches v→x. If a prefix exists, its first vertex has no incoming edge from L or v, so x must beat it; that condition also lets x reach all L. With no prefix, x must beat some core vertex, after which core strongness reaches the entire core and suffix. A suffix alone cannot repair an empty core cut.

## All-length exclusion of47 and63 within the strong template

This is an elementary inverse bound, not an extrapolation of a depth16 scan.

For N=H(closed L)<=63, τ<=α gives α<=N−h. Let A be the initial core α, q the initial core q, and P the integer encoded by the successive prefix bits (least significant first in the prepend procedure). After the prefix,

    α0=A+qP.

If the prefix length a is positive, strongness forces its last prepended bit to be1, so P>=2^(a−1). Each suffix step satisfies α_new>=2α_old, so after b steps,

    2^b α0<=N−h.

For either H9 core, q=84>N−h, so a=0. A nonempty core cut has A>=9, hence b<=2. For the H5 core q=32 and N−h<=58, so P<=1 and a<=1. If a=0, the nonempty core cut has A>=5 and b<=3. If a=1, α0>=32 and b=0. These positive lower bounds on A are verified from all actual core-cut states, each independently computed by literal subset DP.

Consequently any strong template attaining47 or63, regardless of padding length, would already occur with total padding<=3. The independent literal enumeration of that full tiny universe contains neither value. This proves the two all-length template gaps. It does not rule out other templates; indeed the explicit standalone exception witnesses attain both values.

## Independent implementation and hostile controls

`audit.cpp` does not import the producer's implementation. It recomputes all80 core-cut states from endpoint-resolved vertex-subset Hamilton-path DP and sums over complementary subsets to count ordered covers. It reconstructs all padded graphs with total padding<=3, with every core cut and every prefix/suffix word, including nonstrong configurations as controls.

Results:

- 3920 literal padded-state checks;
- 3920 literal closed-tournament Hamilton-path counts;
- strongness checked by reachability, agreeing exactly with the stated restriction;
- omitting the τ subtraction produces3690 detected failures;
- all499987 nonzero recipe records decode within the stated bit/range limits and give their exact target values, both by the state compiler and independently implemented closed binary formula using the direct selected-bit pair sum;
- retained recipe use by core:823,290666,208498;
- maximum retained order22;
- an independent depth-first enumeration over prepend and append choices visits exactly88735747 legal recipes, giving precisely the same target support and gaps47,63 as the producer's loop-based scan.

The displayed ordinal composition law, binary formula(4), prefix/suffix bit conventions, and elementary all-length finite membership algorithm in the producer proof were checked algebraically. The literal controls include every configuration needed for the all-length gap proof. The depth16 scan certifies finite coverage only. The recurrence proof, not the scan, gives an all-length construction interface.

`check_exceptions.py` independently enumerates all7! and8! vertex permutations for the two standalone witnesses. The counts are47 and63, and each supplied Hamilton cycle is valid, certifying strongness.

## Reproduction

From a workspace containing this audit directory and the producer packet:

    g++ -O2 -std=c++17 session216-three-core-audit/audit.cpp -o /tmp/three-core-audit
    /tmp/three-core-audit session216-arithmetic/three-core-release/CORE_STATES.txt session216-arithmetic/three-core-release/RECIPES.bin
    python3 session216-three-core-audit/check_exceptions.py session216-arithmetic/three-core-release/EXCEPTIONS.json

Paths are command-line inputs, so the audit is portable. No external packages are needed. The producer's star theorem is consistent with the direct closure derivation above; no external order-minimum theorem is needed.
