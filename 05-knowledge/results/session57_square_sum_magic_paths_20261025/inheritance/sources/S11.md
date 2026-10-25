# Correlated overlap gluing with original choice identities

25 October 2026. **PROVED** elementary statements; **FINITE-EXACT** replay controls. Classical relational/provenance and join-tree mechanisms are cited, not claimed new. No claim settles LRC(14), arbitrary native shared-pin realizability, or constant-time arbitrary overlap queries.

## 1. Correct object, exact map, and the hostile it repairs

A component identity i refers to one physical choice, not one occurrence in a query expression. Let S_i be its finite set of actual local states (including separate state IDs when different arrays have identical collision masks). A partial assignment a on scope I chooses a(i) in S_i for each i in I. Two rows a,b are compatible precisely when a(i)=b(i) for EVERY i in I intersection J. Define a join b to be their unique union when compatible and failure otherwise.

A block relation R_I is a SET of partial assignments on its fixed scope I. Its natural join with R_J is the set of compatible unions. There is a bijection

    R_I join R_J <-> {u on I union J : u|I in R_I and u|J in R_J}.

Indeed restriction recovers the unique ordered pair, and compatible union is its inverse. Consequently join is associative and commutative, and R_I join R_I = R_I. The empty-scope singleton is an identity and the empty relation is absorbing. Triple associativity follows because either parenthesization succeeds exactly when all three assignments agree wherever they overlap.

For independent CRT components, define R_I as the full product of the local sets S_i, indexed by original identities. Then R_I join R_J = R_(I union J), regardless of overlap. A chosen local state has a fixed-cell collision mask c_i(s). The observable is

    C(u) = intersection over i in domain(u) of c_i(u(i)).

CRT with pairwise-coprime moduli gives a bijection from the full assignment rows to actual full-modulus arrays, and their collision masks equal C(u). Repeated occurrence of an identity is removed BEFORE evaluating the observable or counting. Repeated masks belonging to different actual choices are NOT identified.

In the frozen hostile [2,5,3], the three local state counts are 8,7,2. The [2,5] block has 56 rows and [5,3] has 14. Their compatibility join has exactly 112 rows and reconstructs exactly the 112 valid arrays modulo 30, none with nine distinct cells. Independent pairing has 784 rows, including 128 fabricated zero collision masks. Merely requiring equality of the F5 collision masks is still wrong for counting: its seven actual states produce seventeen same-mask ordered pairs. Shared state identity supplies the exact repair implemented here; the seventeen-versus-seven check refutes the naive same-mask pair-counting rule, not every possible correlated mask-level repair.

### Sufficiency, not universal minimality

Full original choice IDs are a SUFFICIENT exact interface, not a claim of a universally minimal or necessary representation. A certified observational quotient may merge states when every relevant future continuation factors through that quotient. Its weights must be aggregated correctly and its repeated occurrences must use the same correlated diagonal choice, evaluated once, rather than independently pairing all representatives. In particular, when masks alone are sufficient for all future consumers, a properly weighted correlated diagonal on mask classes can preserve the desired counts without retaining each original array ID. The F5 seventeen-versus-seven test rejects only the naive same-mask independent-pair weighting. If an actual array witness or a later constraint distinguishes representatives, the quotient must additionally retain suitable witness information or be refined. This is the inherited continuation-equivalence criterion, not a universal requirement to keep every physical label forever.

## 2. Weights once, no division, and the arbitrary-law distinction

Attach a nonnegative scalar weight w_i(s) to each original local choice. For an assignment u set

    W(u) = product over distinct i in domain(u) of w_i(u(i)).
    H_R(z) = sum over u in R with C(u)=z of W(u).

This is counting when each state has weight 1 and an independent base law when each local weight vector sums to 1. Constraints can restrict R; their total weight need not be 1. Conditional normalization, when requested, is a separate operation requiring positive total mass.

This evaluation uses the union assignment, not W(a)W(b). Every shared weight is evaluated once. It therefore works with zero weights without dividing by a shared marginal or shared product. Scalar block totals, or even scalar block-weighted mask histograms, do not retain enough information to perform this operation.

For ANY finite law f = sum_s w_s delta_s on masks, including non-Dirac laws, the correlated copy is

    Delta(f) = sum_s w_s delta_(s,s),

not f tensor f = sum_(s,t) w_s w_t delta_(s,t). For meet m(s,t)=s intersection t,

    m_*(Delta(f)) = sum_s w_s delta_(s intersection s) = f.

This uses only idempotence on chosen states. It applies equally to finite signed formal measures. It does not contradict the independent-meet probability result that idempotent probability laws on this finite semilattice are Dirac: that result concerns m_*(f tensor f), a different input coupling. Correlated copying is linear in f; independent self-copying is quadratic. If latent choices have the same mask, the displayed mask-law equation remains true, but it does not restore their discarded provenance.

For a genuinely correlated global distribution, one may instead retain the full global state ID and its weight, or a correct explicit factorization with identified factors. Independent local product weights are then not silently assumed. Overlapping marginal distributions alone do not generally determine the global law.

## 3. Formal choice projectors explain the safe order of operations

Introduce symbols x_(i,s), commuting across all identities, with relations

    x_(i,s) x_(i,t) = 0 if s != t;
    x_(i,s)^2 = x_(i,s).

For a partial assignment a put X_a = product_i x_(i,a(i)). Multiplication gives X_a X_b = X_(a union b) for compatible rows and zero otherwise. A concrete model is the free module on partial assignments, with this multiplication and zero for incompatible pairs; associativity follows from the preceding compatibility proof. No relation sum_s x_(i,s)=1 is imposed here.

On FIXED-SCOPE set relations, P_R = sum_(a in R) X_a satisfies P_R P_T = P_(R join T), since restriction makes each compatible union have a unique pair of preimages. In particular P_R squared = P_R. For arbitrary sums of mixed scopes or coefficients this last conclusion need not hold; multiplicities are not magically removed. For example (2 X_a)^2 = 4 X_a.

Evaluate each canonical monomial once by E(X_a)=W(a), extending E linearly. E is generally NOT a semiring homomorphism: if w_i(s)=1/3 then E(x_(i,s)^2)=1/3, whereas E(x_(i,s)) squared=1/9. For different states of one component, E(x_(i,s)x_(i,t))=0 even if both separate weights are positive. The correct order is formal compatibility/union first, scalar evaluation afterward. Evaluating first reproduces the independence error. These orthogonal CHOICE projectors must not be confused with collision-mask meet: two different choices can have intersecting or identical collision masks, but their same-component projectors multiply to zero.

This is a concrete constrained-provenance construction, not an assertion that ordinary polynomial provenance automatically has these relations. Green, Karvounarakis and Tannen provide the classical K-relation/provenance setting; the explicit quotient and once-only evaluation above are specified and proved here.

## 4. General constraint gluing: exactness is not marginal reconstruction

Let each block relation encode its stated local constraints. The natural join of all block relations is EXACTLY the set of assignments satisfying the conjunction of those constraints. This claim holds with cycles too; it does not promise efficient computation or existence. Any omitted global constraint must still be retained as another factor/relation.

A stronger and generally FALSE claim is that arbitrary projections of an unknown global relation recover that relation. For example let R be the four even-parity assignments to three bits a,b,c. Both projections R_ab and R_bc are full two-bit relations. They agree completely on b, but their join contains all eight assignments. Retaining every shared assignment does not repair constraints already discarded during projection. Exact reconstruction from projections needs a lossless-join hypothesis or an equivalent conditional factorization. Likewise a separator division formula for arbitrary probability marginals is not justified without conditional-independence assumptions and careful support treatment.

This boundary is important for the CRT application: its block relations really ARE products of independent local choices, so the required losslessness is proved, not inferred from matching histograms.

## 5. Classical join-tree sufficient criterion, with an explicit proof

Let a finite family of scopes I_v be arranged on a tree. Assume running intersection: for every variable i, the vertices whose scopes contain i form a connected subtree. Write K_uv=I_u intersection I_v for an edge separator. Suppose finite set relations R_v are nonempty and are edgewise fully consistent: every row on either endpoint has a compatible row on the other endpoint (equivalently, their separator projections agree).

**Theorem.** Every row of every R_v extends to a global satisfying assignment. Hence the global join is nonempty.

Proof. Root the tree at a selected row. Choose compatible child rows and recurse. When a child subtree is added, any variable it shares with the already chosen part must occur on the connecting separator, by running intersection. Agreement on that separator therefore suffices; no nonadjacent branch can introduce a new conflict. This constructs the extension. The argument works rooted at every row, giving the stronger assertion.

Full edge consistency can be obtained without changing the global join by deleting unsupported rows (semijoins): first child-to-parent bottom-up, then parent-to-child top-down. Bottom-up removes precisely rows unable to extend into descendant branches; top-down removes those unable to extend toward the root. The same induction proves that the survivors are exactly the projections of the global join. If any relation empties, the global join is empty. This is a classical acyclic-database mechanism, not a new realizability theorem for arbitrary cyclic wiring.

The running-intersection condition matters. On bits take a=b on scope ab, b=c on bc, and a!=c on ac. Every relation is nonempty, every pair has compatible extensions for every row, and every pairwise separator is fully consistent. Yet the global join is empty. The three scopes cannot be organized as a running-intersection tree.

### Exact weighted separator recursion

Root a valid join tree and assign each original variable's local weight to ONE owner node containing that variable; never to all occurrences. Let O_v be those owned variables. For node v with parent p and children t, define

    M_(v->p)(k) = sum over a in R_v with a|K_vp=k of
        [product_(i in O_v) w_i(a(i))]
        [product_(children t) M_(t->v)(a|K_tv)].

At the root omit the separator restriction and sum the same expression over its rows. This equals the total weight of the full join. Proof: conditioned on the separator values, descendant branches have no unexposed shared variable by running intersection; their choices combine bijectively, and each original weight occurs exactly at its owner. Apply induction from leaves. General factors can likewise be assigned once to containing scopes. No division, positive-weight assumption, or prior semijoin reduction is required for this sum-product formula. Keeping a collision-mask channel in each message and intersecting channels gives the full weighted collision histogram, at additional signature cost.

LCA may organize the scope tree or find where information meets, but it does not replace the actual separator assignments or prove that running intersection holds. The inherited Session25 named-port theorem already proves exact two-fragment colouring gluing and necessary boundary information. The present lane specializes that information-retention principle to original CRT component choices and states a classical tree-shaped sufficient extension; it is not a replacement for native geometric/physical realizability checks.

## 6. Interface invoice and loss ledger

A complete separator on k shared independent choices has up to product_i |S_i| assignments, or d^k when each has d choices. This is an explicit-table worst-case size, NOT a universal lower bound for every compressed encoding. Full block relations can be larger still. A nested-loop join of relations with A and B rows makes A*B compatibility tests, each inspecting shared IDs, plus union/output work. Hash indexing can avoid incompatible cross-products but does not remove output size or separator-state size. Separator DP tables and arithmetic depend on their actual support and coefficient bit lengths.

Thus one relational operation per overlapping sparse-table query is not O(1) machine time or space. There is no free arbitrary-overlap improvement over the earlier disjoint-query repair. Efficient special cases may exploit small separators, special factorization, or compressed representations, each requiring its own proof.

Connection ledger:

- Source objects: fixed labelled collision masks, independent CRT choices, Session25 named boundary states, classical relation/provenance operations.
- Map: original component assignment -> compatible union assignment -> once-only weighted collision-mask histogram.
- Preserved: every actual compatible assignment, its component identities, exact multiplicity/weight, CRT witness, and distinctness existence.
- Lost by independent pairing: identity correlation; hostile 112 becomes 784 and invents 128 distinct rows.
- Lost by mask-only shared keys: different actual choices with equal masks; F5 seven states become seventeen same-mask pairs.
- Lost by premature marginalization: cross-block constraints, as even parity demonstrates.
- Sufficient sidecars used here: original state IDs, complete future separator assignments, fixed cell labels, unique weight/factor ownership, and CRT coprimality when interpreting arrays. Certified continuation-equivalent quotients may compress this interface; raw physical IDs are not universally minimal.
- Scope limit: exact abstract relation joining is not itself a physical realizability classifier; tree consistency is sufficient only under the displayed running-intersection hypothesis.

## 7. Replay

Run `python verify.py --verify` and `python -O verify.py --verify`. The verifier uses explicit failures, not removable assertions. It checks 19,683 partial-assignment associativity triples; all 4,096 triples of binary relations on ab,bc,cd, including empty relations; exact semijoin global-extension and weighted DP identities with zero weights; 80 rational diagonal laws; the full mod-30 histogram and actual CRT array bijection; a separately reconstructed line-equation support control; mask-only counting failure; parity projection failure; and cyclic pairwise-consistency failure. Finite checks supplement, not replace, the all-size proofs.
