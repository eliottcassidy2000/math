# Root-module defect slots and the limits of fixed transitive multiports

2026-11-04. All-order evaluator theorem, finite independent controls, and an explicitly inherited obstruction with a ternary-branch corollary. No closed two-branch grammar, native legality certificate, or new global spectrum claim.

## 1. Source map and scope

This packet starts from S180 at commit `c5333ffd05e1fadf9ea5147c1e1dd9155abd5252` and its `ports/PROOF.md`, `older/PROOF.md`, and synthesis. The physical zero-ready counts 170807, 512423, and 1537271 are inherited; two successive +2 returns are not closure under both branches.

Relevant older exact paths:

- `01-canon/theorems/THM-4145-rooted-homogeneous-pair-expansion-two-defect-formula.md`: root pair expansion, good insertion slots, and zero/one/two deletion defects. The general module formula below extends its explicit defect decomposition.
- `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`: arbitrary quotient substitution with complete path-cover profiles. The theorem below is a root-specialized, defect-indexed evaluator of this existing general substitution law, not a claim that general substitution is new.
- `01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`: the complete ordered path-cover profile is a jet, and the scalar Hamilton count loses necessary information.
- `01-canon/theorems/THM-4162-rooted-pair-mixed-two-ear-tensor-and-enumeration-free-johnson-cosets.md`: reconstructs complete one-birth cut capacities and computes exact fixed-cardinality cosets in O(N^3). This is a necessary arithmetic filter, never a claim that every allowed lattice point occurs. Its child tensor is not determined by the parent's scalar count.
- `05-knowledge/results/session175_tournament_spectrum_20261104/clone_obstruction/PROOF.md`: already proves the transitive-inflation growth and polylogarithmic bound used in Section 4. We claim only its application to the S180 binary ternary-branch tree.

No startup helper or `01-canon/MISTAKES.md` was read or scanned.

## 2. Arbitrary root-module formula

Let (Q,r) be a rooted tournament of order q>=2. Replace r by a nonempty tournament B, retaining every other physical name and making all edges between B and the outside uniform with the former adjacency of r. Write Q[r<-B] for this substitution. Let F_B(l) count ordered covers of B by l nonempty directed paths, with F_B(0)=0. Let P_B(t) count ordered covers by t slots, allowing empty slots; P_B(0)=0 because B is nonempty. Then

P_B(t)=sum_l binom(t,l) F_B(l).

For each permutation pi=(u_1,...,u_m) of Q−r, number insertion slots 0,...,m. Let D(pi) be the internal slots i with u_i not pointing to u_(i+1). Let G(pi) include slot 0 if r->u_1, slot m if u_m->r, and internal i if u_i->r->u_(i+1). If D is not a subset of G, this permutation contributes zero. Otherwise put d=|D| and g=|G|. Its exact contribution is

C_B(g,d)=sum_(l=d)^g binom(g−d,l−d) F_B(l)
         =sum_(j=0)^d (−1)^j binom(d,j) P_B(g−j).       (1)

Consequently, with N_r(g,d) the number of eligible permutations with that pair,

H(Q[r<-B])=sum_(g,d) N_r(g,d) C_B(g,d).                (2)

### Proof

Delete B from a Hamilton word of the substituted tournament. The surviving outside permutation is pi. Every backward adjacency must have had a nonempty B-run inserted, and any occupied slot must satisfy the root adjacency on both existing boundaries. Thus all defect slots must be occupied and all occupied slots lie in G. Conversely a choice of occupied slots containing D, together with an ordered path cover of B assigned in the order of those slots, reconstructs one and only one valid Hamilton word. For l occupied slots the choices are binom(g−d,l−d); the path cover has F_B(l) choices. This proves the positive formula. Inclusion-exclusion on the d mandatory occupied slots gives the second expression. Deletion and reconstruction are inverse, proving (2).

The dependency on B is its complete path-cover profile through l<=q. It cannot generally be replaced by H(B). If B itself is later changed, its profile must be updated exactly. The root-conditioned signature N_r is not asserted to be a sufficient multiport endpoint state; physical endpoint/guard information must remain elsewhere.

### Arbitrary transitive-clone depth

For B=TT_k, each vertex independently chooses one of t slots and appears in increasing order within that slot, so P_B(t)=t^k. Therefore

C_TT_k(g,d)=sum_(j=0)^d (−1)^j binom(d,j)(g−j)^k
          =sum_(l=d)^min(g,k) binom(g−d,l−d) l! S(k,l). (3)

Here S(k,l) is the Stirling number of the second kind. The second expression is nonnegative and counts surjections onto the chosen occupied slots. It is zero for d>k, displaying the exact growing defect budget. For k=2 it gives g^2, 2g−1, and 2 for d=0,1,2, recovering THM-4145. One must not retain only the zero-defect term: the directed triangle's pair expansion is 5, whereas zero-defect squaring gives 4.

Once N_r is known, all k can be evaluated by powers of integers between 0 and q, rather than a new exponentially growing vertex DP. The sequence H(Q[r<-TT_k]), k>=1, is an integer linear combination of a^k for 1<=a<=q and therefore satisfies the constant-coefficient recurrence with characteristic polynomial product_(a=1)^q (x−a). Vanishing coefficients may lower its order. This recurrence is for one fixed rooted core only and does not imply ternary branch closure.

Strongness is preserved and reflected for q>=2: nonstrong quotient cuts lift; for a strong quotient, directed quotient routes connect all outside/module endpoints, and a cycle through r supplies a return route between any two module vertices regardless of B's internal direction. This extends THM-4145's same argument.

## 3. Concrete construction/search use

For a fixed rooted core, compute N_r once. Search replacement modules through their true ordered path-cover profiles, and use (2) as an exact count filter. This admits arbitrarily large, internally nontransitive modules and arbitrary defect budgets; it is not confined to two-ear repairs. For TT_k, (3) is a cheap all-k exact evaluator. A candidate still needs the physical joint endpoint transfer, directly recomputed zero-ready root, strongness if required, and literal native guard witnesses. S180's nilpotent named endpoint matrix is a sufficient unbounded sidecar for those endpoint transfers; independent marginal profiles are not.

The existing THM-4162 capacity tensor can then reject impossible birth-count cosets before enumerating actual cuts. Passing a coset test does not produce a cut and is not a proof of realizability.

### Complete finite reduction for an arbitrary-module target

Let Q be a fixed strong quotient of order q>=3, and let all modules B_i be nonempty and otherwise arbitrary. Choose a directed triangle through quotient vertex i, and one physical vertex from each of its other two modules. The induced tournament is C3[B_i,K_1,K_1], whose exact count from THM-3121 is

3H(B_i)+F_(B_i)(2)=H(B_i)+P_(B_i)(2).

Every induced sub-tournament has a Hamilton path, so the subset-split formula for P_(B_i)(2) has 2^{|B_i|} positive terms. Insertion monotonicity therefore gives

H(Q[B_1,...,B_q]) >= 2^{|B_i|}+1 for every i.           (4a)

Consequently every solution of a fixed target N has |B_i|<=floor(log_2(N−1)). This extends the S175 transitive-module size bound to arbitrary module interiors, using the inherited THM-3121 triangle formula and S180's positive path-cover bound. Exhausting all labelled module tournaments through this bound and evaluating the exact substitution kernel is a terminating complete decision procedure for the fixed-skeleton target problem. In the one-root case (2) provides the evaluator. The enumeration remains enormous; no efficient unrestricted-spectrum algorithm follows. Unlike transitive modules, the number of arbitrary interiors through this bound is not polylogarithmic, so this does NOT extend the polylogarithmic obstruction to arbitrary modules.

## 4. Inherited growth bound, new branch-tree application

S175 already proves: if Q is a fixed strong quotient of order q>=3, then

H(Q[TT_(k_1),...,TT_(k_q)]) >= 2^(k_i)+1 for every i.

The proof uses a directed triangle through i, the exact triangle count 2^(k_i)+1, and Hamilton-path insertion monotonicity. Thus through X>=3 there are at most floor(log_2(X−1))^q parameter tuples. For K fixed strong templates of order at most d there are at most

K floor(log_2(X−1))^d                                  (4)

distinct attainable counts. This is precisely the strong-template part of S175, not a new growth lemma.

Combine (4) with S180's branch injectivity: after n maps H->3H+2 or H->3H+4 from H_0, there are 2^n distinct outputs, all at most X_n=3^n(H_0+2)−2. Therefore any representation of all depth-n outputs in this template bank would require

2^n <= K floor(log_2(X_n−1))^d.                        (5)

The right side is polynomial in n, so no such fixed finite transitive-module bank closes both branches. Unlike the S180 fixed-suffix-clock statement, this corollary also allows multiple independently unbounded transitive modules in each fixed strong skeleton. It is an application of S175 plus S180, not an obstruction to arbitrary changing suffixes or nontransitive module families. Arbitrary ordered products are NOT included in (4) or this corollary.

If K_n strong templates of order at most d_n are used at depth n, (5) also gives the necessary information-growth condition

d_n log_2(floor(log_2(X_n−1))) + log_2 K_n >= n.

For bounded K_n this forces d_n=Omega(n/log n). This measures the amount of genuinely growing template structure needed in that restricted transitive-module representation.

## 5. Independent finite controls

`verify.py` imports no repository code and uses only the Python standard library. It independently computes literal child counts by subset endpoint DP, evaluates (1) through slot-cover convolution, and enumerates outside permutations to construct N_r. It checks every rooted labelled quotient of orders 2 through 4 against every labelled module of orders 1 through 4. It separately checks every transitive inflation k=1,...,6 of those rooted quotients, plus strongness preservation/reflection and the inherited exponential lower bound on strong rows. The directed-triangle control has counts 3,5,9,17,33,65.

The all-order results are proved above; finite checks are controls, not extrapolation. Run `python verify.py` and `python -O verify.py`. Exact receipts record the check totals.
