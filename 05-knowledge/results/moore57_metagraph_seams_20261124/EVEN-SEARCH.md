# Even Gewirtz transport leaves one sharp-set seed

**PROVED reduction with finite inventories.** This applies to a Moore seam system using the degree-56 Gewirtz action, not to every sharp 56-set in the full group. The [parity gluing theorem](moore-tournament-holonomy.md) confines every normalized seam to K=G∩A56. This is a normal subgroup of order 40,320. The collineation-only subgroup fails the projective incidence denominator, so this is the one surviving ambient Gewirtz group.

Let e(g) be the number of labels moved by g to an adjacent label in the Gewirtz graph. Among even derangements the exact counts are

| type | e=8 | e=14 | e=20 |
|---|---:|---:|---:|
| collineation | 2520 | 5760 | 0 |
| correlation | 10080 | 0 | 2520 |

For each g in a sharp set S, the previous [relative count](../moore57_gewirtz_exact_cover_20261124/SEARCH.md) gives n14(g)+2n20(g)=20 among the other rows of g⁻¹S. At least 28 members of S have one geometric type; any such member has a distinct same-type partner with e(g⁻¹h)=8, as there can be at most 20 with e=14. Translate by g⁻¹ within K. The 2520 even collineation derangements with e=8 form one conjugacy class under the *full* group G. Conjugating by G is legal since K is normal.

Existence of a sharp set in K is therefore equivalent to an exact cover with identity, one fixed seed and 11,212 legal other rows (3400 collineations and 7812 correlations). Exactly 54 must cover all 3024 remaining ordered cells. The [model](even-seed.py) checks the counts without a solver; optional MILP verifies any successful cover before saving. The earlier two-branch test can use a witness if found.

The 75-second HiGHS run of this reduced search reached its time limit at the root LP node, with no incumbent or infeasibility certificate. That is not evidence of nonexistence.
