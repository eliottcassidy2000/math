# A bounded q5 escape: one strong four-vertex terminal cap

2026-11-03. **FINITE-EXACT complete fixed-carrier census.** For the specified
deleted graph D=TT8->C4, none of its4096 one-apex stars supplies H1409 or1427,
outerroot99 end defect0, strongness and a start defect1 or2 together.
This is not a classification of arbitrary q5 graphs or other C4 placements.

The closest mechanism is the exact [deleted-joint star formula](../q3/PROOF.md).
The q3 terminal/initial-cap obstructions motivate replacing the three-cycle
with the strong four-vertex tournament. The corrected near miss is count-only
acceptance: one new graph attains1409 but fails both port requirements. The
least-used sidecar is the actual end-zero gauge bank, retained alongside
the start-entry bank. The live board is q5 deletion / exact star counts /
fixed root99 / joint ports / count-preserving compensation.

## 1. One declared actual carrier

The retained inner mark is physical0. The twelve other names are arranged as

    TT(99,10,11,2,1,3,4,12) -> C4(200,201,400,401),

where the C4 outgoing rows in that displayed order are [8,1,3,6]. Its five
Hamilton paths are checked by all24 literal vertex permutations; its
two-slot cover count is32. Therefore H(D)=5. This actual carrier preserves
all thirteen names from the original children. Relative to their common
deleted graph it changes exactly three non-star pairs:

    {200,201}, {200,400}, {201,401}.

All other changes in the census are edges incident with0. No native program
or strongness of intermediate edits is asserted.

For x_i=[0->i], the fresh deleted E_i,S_i,J_ij coefficients give

    H=sum_i(1-x_i)E_i + sum_j x_jS_j
                       +sum_(i,j)(1-x_i)x_jJ_ij.

As in the proved parent formula, the three terms place the apex last,
first, or between its two path segments. It classifies every star without
an order cutoff or an isomorphism quotient.

Strongness is equivalent to0->99 and at least one cap vertex->0. Necessity
comes from the source singleton and terminal SCC of D; these two links also
give paths between the apex and every vertex. Independent forward/reverse
reachability confirms the criterion for all4096 cuts. There are1920 strong
states. In a strong state, root99 has sole in-neighbor0, so

    E_99(2)-H=S_0(T-99,1)=5[0->10].

The root guard therefore leaves960 states.

## 2. Exact outcome and its failed implication

| Desired H | Strong exact-count states | With root99 end-zero | With any start1/2 | Jointly compatible |
|---:|---:|---:|---:|---:|
| 1409 | 1 | 0 | 0 | 0 |
| 1427 | 0 | 0 | 0 | 0 |

The sole strong1409 state has physical incoming neighbors of0

    11,2,3,400

and outgoing neighbors

    10,12,1,4,99,200,201,401.

Its exact root99 end defect is5. More strongly, its full endpoint profile
has no end-zero vertex anywhere and no start defect1 or2 anywhere. Thus
merely changing the root's label does not rescue this count-only witness.
Reversing its edge {0,10} restores the root-end guard, but changes H to975.
That last value is an exact second word in the already complete census;
it is not a further search or a proposed count-neutral repair.

`witnesses.json` retains the unique actual graph, all names and changed
pairs relative to its actual H1409 source, incoming/outgoing cuts, every
port defect, the full E1/S1/E2/S2/J profile and fresh golden/F49/F343
observers. `census.json` records all4096 cuts and their exact counts and
strongness. The observers are passive sidecars, not acceptance tests.

## 3. Verification and scope

The deleted graph and C4 have fresh literal subset profiles. The strong
target-count survivor is independently fully profiled. Twenty-six fixed
weak/strong/guard-boundary star controls and both original children are
also literally counted. Together these are30 literal control profiles
plus the target profile, and the cap's independent24-permutation check.
No Python assert is used, so optimization cannot remove the checks.

Reproduce from repository root:

    python 05-knowledge/results/compensation_returns_20261103/star/q5/reproduce.py

Normal and `-O` outputs and artifacts must be byte-identical. Every earlier
packet remains pinned and read only. The result closes exactly this next
carrier; moving C4 elsewhere, changing the deleted strong component, or
varying more than the apex star is outside its universe. No broader q5
impossibility or return closure is inferred.
