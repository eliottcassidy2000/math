# Bridge–sink star kernel

**PROVED (direct combinatorial argument; two independent audits passed).** This is a sufficient singularity mechanism, not a converse. No strong-connectivity or unique-Hamiltonian-cycle assumption is required.

## Statement

Partition a tournament into R and Q, with |R|≥3 and |Q|≥2. Assume:

1. In R, vertex a has exactly one in-neighbor s.
2. The only R→Q arc is q→r, where q∈{a,s} and r∈Q.
3. Q−r has a universal sink c, and c→r.

Let M2 count ordered spanning covers by two nonempty directed paths, rows indexed by ordered starts and columns by ordered ends, singletons allowed. Then

    M2 ((e_a−e_s)∧e_r) = 0,

where e_i∧e_j denotes the ordered-pair vector e_(i,j)−e_(j,i). In particular M2 is singular.

## Auxiliary identity inside R

Write H(t,x) for Hamilton-path counts in R, and S((u,t),(x,y)) for the antisymmetrized two-path-cover entry (zero for repeated starts or ends). Write h(t) for Hamilton paths in R−a from t to s, extended by zero at t=a, and L=Σ_t h(t).

Because a's sole in-neighbor is s,

    H(t,a)=h(t),        H(t,s)=1_(t=a) L.

The latter holds because a must be the first vertex of any Hamilton path ending at s, and a dominates every possible initial vertex of a Hamilton path in R−a ending at s. Here |R|≥3 ensures such a path does not begin at s.

A two-path cover ending at s and a must have a as a singleton, since s is already the other endpoint. Consequently

    Σ_u S((u,t),(s,a)) = 1_(t=a)L − h(t)
                       = H(t,s)−H(t,a).                    (1)

In particular Σ_t H(t,a)=Σ_t H(t,s)=L.

## Rowwise proof

Set P=Q−r. We compare antisymmetric endpoint columns (a,r) and (s,r). In any cover counted in these columns, r is an endpoint, so the only R→Q arc q→r can only occur at the end of a path. Therefore every vertex of P occurs in an initial Q-segment of one of the two paths; once a path enters R it cannot return to P.

**Both starts in R.** Impossible, because P is nonempty and cannot be reached.

**Starts (z,t) with z∈P and t∈R.** All P lies on a single initial path segment starting at z. Since c is a universal sink of P, that segment ends at c. Let F(z) be the number of Hamilton paths in P from z to c. For endpoint x∈{a,s}, the antisymmetric entry equals

    F(z) [ Σ_u S((u,t),(x,q)) − H(t,x) ].                  (2)

The negative term is the case where the first path ends c→r and the second spans R. The sum is the case where the first enters R and one of the R paths is extended q→r. If x=q the sum is zero. Equation (1), with antisymmetry in x,q, makes (2) equal for x=a and x=s, for either q=a or q=s.

**Both starts in P.** Fix an ordered cover of P by its two initial Q-segments. When both paths also contain nonempty R-segments, their starting vertices u,v in R may be chosen freely: every vertex of P dominates every vertex of R. Their contribution is proportional to

    Σ_(u,v) S((u,v),(x,q)) = 0,

by antisymmetry in u,v. When just one path contains R, that path spans R and the contribution is a Q-dependent signed coefficient times Σ_u H(u,x). This is equal for x=a,s. These cases exhaust the row.

**One start is r.** Its path is the singleton r. The other path covers P and R, necessarily as a Hamilton path of P ending c, followed by a Hamilton path of R. Its contribution is a coefficient independent of x times Σ_u H(u,x), again equal. Reversing row order only changes the sign, so all ordered rows have been covered.

Thus the two antisymmetric endpoint columns agree, proving the claimed nonzero star kernel.

## Actual grammar: birth and persistence interface

At a `be` birth, the A-state maximum r was inserted immediately after an old maximum a's predecessor s. The new top C3 is r→b→c→r, with bridge s→r, and a retains indegree1 in R with sole in-neighbor s. Thus q=s.

At an `ee` birth, before the two letters the B-state maxima are a→r, and s is a's predecessor on the unique Hamilton cycle. Removing r leaves a with sole in-neighbor s. The new top C3 is r→b→c→r, with bridge a→r. Thus q=a.

### All-order persistence after the last bad pair

Choose the last occurrence of `be` or `ee` and freeze its partition R,Q and vertices a,s,q,r,b,c. Initially Q=(r,b,c) is the directed triangle r→b→c→r, and the state-B maxima are b→c.

Every R vertex has indegree at least2, since both b and c dominate R. Therefore no R vertex is a maximum. The only R→Q arc is q→r, and r has indegree at least2 (from q and c). Therefore no R vertex is the Hamilton-cycle predecessor of a maximum either. Since every legal parent is a maximum or the predecessor of one, every future parent stays in Q. This argument remains true inductively: added vertices dominate all R, and the one-bridge pattern is preserved.

The vertex c cannot be the predecessor of a maximum: all its out-neighbors are in R∪{r}, and each already has indegree at least2. If the suffix is nonempty, its first letter cannot be e, since the birth ended in e and that would create a later bad pair. Its first letter is therefore c or d, neither of which chooses c as parent. The new vertex dominates c, so c itself now has indegree at least2. Before any putative first future choice of c as parent, every added vertex dominates c; hence its out-neighborhood remains R∪{r}. Thus c is neither a maximum nor the predecessor of one, a contradiction. Consequently c is never selected.

Every future vertex therefore dominates c. The frozen c stays a universal sink of the growing Q−r, while c→r, the R tournament, its indegree-one vertex a, and the bridge q→r all remain unchanged. The bridge–sink lemma supplies the same fixed nonzero star kernel throughout the suffix.

**Conclusion (PROVED, independently audited): every legal grammar word containing `be` or `ee` has singular native M2, at every order.** The converse is not established here. No monotonicity of the entire kernel or nullity is asserted; the witness is chosen at the last bad pair.

## Exact checks and limits

`check_birth_star.py` exhausts arbitrary R with prescribed a←s, |R|=3,4,5, both q choices, and Q=C3: 148 cases, all pass. `check_star_persistence.py` exhausts all legal words through length7/order11. Every forbidden word has a two-term star, and every extension not ending be/ee preserves at least one identical old star: 861 forbidden words across lengths2..7. These finite checks are corroboration, not the all-order proof.

`check_arbitrary_root_sink.py` independently exhausts all orientations of Q with the prescribed sink c of Q−r and c→r, |Q|=3,4,5, taking R=C3 and both bridge choices: another148 cases, all pass. `check_last_birth_certificate.py` constructs a concrete star from the last bad pair rather than searching the final kernel. It directly verifies the same fixed star against native M2 for all861 forbidden words through order11 and verifies parent≠c throughout each suffix. Its receipt is `last_birth_receipt.json`.

## Independent reviews

Both [review 1](../audit/ALL_ORDER_SPECTRAL_PROOF_AUDIT.md) and [review 2](../audit2/REPORT.md) passed. The complementary [acyclicity theorem](../triangularity/AVOIDING_WORD_TRIANGULARITY.md) proves the converse and all ordered-k observer unimodularity on avoiding words.
