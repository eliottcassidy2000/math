# A finite representative bound for every rooted four-state

Status: proof supplied; PENDING INDEPENDENT AUDIT. Finite exact controls are distinct from the all-order proof.

## Conventions and theorem

A tournament is finite. Directed paths are ordered, with one empty path when a bin may be empty. Set H(empty)=1. For a nonempty tournament A rooted at a, write

    x(A,a)=(D,E,U,d),
    D=H(A), E=# Hamilton paths ending at a, d=H(A-a),
    U=# ordered spanning two-path covers (P,R), with P possibly empty
      and R nonempty and starting at a.

These are exactly the Session225 graft-recovery conventions. Set V=U-d.

**Theorem.** Every physical four-state has a representative of order at most

- 1, if d=1 and U=d;
- d+1, if d>1 and U=d;
- floor(log_2 V)+2, if d=1 and V>0;
- d+floor(log_2 V)+2, if d>1 and V>0.

Consequently physical four-state membership is decidable by a finite exhaustive search with an explicit bound computable from the input coordinates. No coverage statement about all odd Hamilton-path counts follows.

The proof uses the standard theorems that every tournament has a Hamilton path (Redei existence), and every nontrivial strongly connected tournament has a Hamilton cycle (Camion). The combinatorial recurrence and simultaneous compression are proved below. No Busch bound is needed.

## 1. SCC profile and exact recurrence

List the strongly connected components of A-a in their unique topological order B_1,...,B_k: every edge between B_i and B_j points from B_i to B_j for i<j. Write H_i=H(B_i). Let Q_i be the count of ordered spanning two-bin directed-path covers of B_i, allowing either bin to be empty. Let alpha_i count such covers (L,S) for which S is nonempty and its first vertex lies in N^+(a) intersect B_i. The other bin L may be empty.

Let A_i be induced by a and B_1,...,B_i, and d_i=H(A_i-a), V_i=U(A_i,a)-d_i. Then

    d_0=1, V_0=0,
    d_i=H_i d_{i-1},
    V_i=Q_i V_{i-1}+alpha_i d_{i-1}.                 (1)

Indeed every Hamilton path on the deletion graph passes consecutively through its SCCs, proving the d formula. Exactly d_i of the U covers have R=(a), so V_i counts covers in which R leaves the root. In a V-cover already active before B_i, both paths can append their portions of B_i, independently of the earlier endpoints: all earlier deletion vertices beat B_i. The R path is already nonempty beyond a, so it cannot return to the root, and each path's B_i part is a single final segment. This gives Q_i choices per earlier active cover, including empty segments. If R first leaves a in B_i, P contains a Hamilton path on all preceding blocks, with d_{i-1} choices, and its B_i portion together with the nonempty R suffix is precisely an alpha_i cover. These disjoint cases exhaust the V-covers and have unique restrictions, proving (1).

The following inequalities and equivalences hold:

    Q_i >= 2H_i >= 2,
    Q_i >= 2^{|B_i|},
    alpha_i >= 0,
    alpha_i>0 iff N^+(a) intersect B_i is nonempty.

The first inequality includes the two covers with one empty bin. For the second, each ordered bipartition of B_i admits at least one path in each part. If a beats v in B_i, choose S=(v) and a Hamilton path L on B_i-v, showing alpha_i>0. The converse follows from the definition.

In particular V=0 iff a is a sink. In that case D=E=U=d.

## 2. Exact simultaneous compression

Suppose C is an initial sequence of deletion SCCs to which a loses every edge, and R is the remaining induced rooted tournament including a. Then

    A = C => R,
    x(A,a)=H(C) x(R,a).                            (2)

Hamilton paths and paths ending at a concatenate uniquely, giving D and E. The deletion count factors in the same way. For U, the a-starting path cannot visit C, so the other path contains all of C, followed by its possibly empty part in R. Its C restriction is a Hamilton path on C. Conversely this concatenation constructs every U-cover uniquely, including the case the second portion is empty. Thus U also scales by H(C).

Every singleton SCC in this initial inactive sequence has H=1. Delete all such singleton SCCs. The remaining sequence has the same H(C), and the rooted suffix R is unchanged. Formula (2) therefore preserves **all four coordinates exactly**, not only d or U. This deletion introduces no new edges and changes none of the retained SCCs.

If V=0, all deletion blocks are inactive, so all singleton blocks can be deleted. If V>0, let j be the first active block (alpha_j>0). The retained tournament consists of nontrivial inactive blocks, B_j, its subsequent blocks, and the root. In the elementary proof below it is harmless to count all nontrivial blocks separately.

## 3. Vertex budgets and proof of theorem

For each nontrivial strong block B_i, Camion gives a directed Hamilton cycle; deleting each of its edges gives |B_i| distinct Hamilton paths. Thus

    3 <= |B_i| <= H_i.

Let h_1,...,h_t be the H_i of the nontrivial blocks. Their product is d because singleton blocks contribute 1. If d=1, there are no nontrivial blocks. If d>1, each h_i>=3, and

    sum |B_i| <= sum h_i <= product h_i = d.        (3)

The sum-product inequality is equality for t=1 and follows inductively from ab>=a+b for integers a,b>=2.

When V=0, only the root and these nontrivial blocks remain, giving order <=d+1 (or order 1 for d=1).

When V>0, equation (1) gives V_j>=1. Every later block multiplies V by at least 2. Hence there are at most floor(log_2 V) blocks after B_j. All surviving singleton blocks are B_j itself, if it is singleton, or subsequent blocks. There are therefore at most 1+floor(log_2 V) surviving singleton blocks. Adding the root and using (3) proves the bounds stated above.

There is a stronger useful resource statement: since Q_i>=2^{|B_i|},

    sum_{i>j} |B_i| <= floor(log_2 V).              (4)

Thus the entire suffix after first activation, not only its block count, is controlled. The stated theorem does not need this strengthening.

## 4. Terminating physical-state and numeric-inverse oracles

For an arbitrary integer tuple x=(D,E,U,d), reject immediately if D,U,d<1, E<0, E>D, or U<d. If U=d, reject unless x=(d,d,d,d). Otherwise compute the appropriate bound above. Enumerate every tournament on 1 through that bound with vertex 0 fixed as root, calculate all four statistics exactly, and accept precisely if x appears. Every rooted tournament can be relabeled to have root 0. The representative theorem proves that exhaustion is a valid negative answer, not merely failure of a bounded heuristic.

The script implements `find_state(x)` for tuple input and returns a labeled witness or None. It is a termination theorem, not an efficiency claim: enumeration can be astronomically expensive. Exact integer bit_length evaluates floor(log_2 V) without floating-point error.

If K is any specified invertible rational rooted-graft transfer matrix, the set of physical four-states in its numeric image is decidable: calculate y=K^{-1}x, reject nonintegral y, then apply the physical-state oracle to y. A positive answer produces a rooted representative whose graft has state x. This does **not** decide whether a particular tournament with state x is graph-isomorphic to a graft, nor does it give a size-decreasing recursive inverse, nor prove eventual odd scalar coverage.

## 5. Comparison with prior canon and optional quantitative strengthening

The current repository already contains scalar H-value decidability by deleting all singleton SCCs and applying the Camion sum-product estimate: `05-knowledge/results/camion_busch_gaps_polyhedra_collatz_20261002.md`, section 4. That result gives an H=m representative on <=m vertices. The same section records a Busch-based O(log m) scalar-core bound. Thus finite scalar H decidability and its Camion core argument are **not new claims here**.

The new scoped content here is the rooted-state issue: singleton deletion is allowed only in the inactive prefix; its preservation of D,E,U,d is simultaneous; and U-d controls the remaining singleton/SCC suffix. The proof does not assume arbitrary singleton deletion preserves the state.

Session197 already proves exact numerical-state membership and an O(K^3) response inventory under U<=K for its restricted physical source class: rooted tournaments having exactly one Hamilton path starting at the root. See `session197/source-signature/PROOF.md`, sections 1-3, recovered and checked directly with `git show HEAD:session197/source-signature/PROOF.md` (the path is tracked but sparse-unmaterialized in this checkout). That result uses the complete strong start-one grammar and weak transitive tails, and also gives an injective three-response signature. The present theorem applies to arbitrary rooted tournaments, without the exactly-one-starting-path restriction, by a different SCC compression argument. It does not replace or improve the efficient specialized grammar/inventory of Session197. No priority claim is made.

OPTIONAL, CITED-INPUT refinement only: the above Camion/Busch note records Busch's minimum strong-tournament path bound

    H(B)>=min {3^a 5^b : 2a+3b=|B|-1, a,b>=0}.

This implies |B|<=3 log_3 H(B): for a+b>=1,

    3 log_3(3^a 5^b)-(2a+3b+1)
      = a+(3 log_3 5-3)b-1 >=0,

since 3 log_3 5-3>1 (equivalently 125>81). Summing over nontrivial SCCs gives a total nontrivial order <=3 log_3 d. Conditional on that cited theorem, replace d in the d>1 bounds by floor(3 log_3 d), obtaining an O(log d + log V) representative bound when V>0. This is a consequence of the cited minimum theorem, not a new proof of it. The self-contained compression proof and exact oracle use only the elementary Camion bound and do not depend on this optional refinement.

## 6. Reproducible controls

Run `python verify.py` and `python -O verify.py`. Both must generate the same JSON receipt. No assertions are used for correctness checks, so optimized execution retains them. The code is standalone and imports no repository helpers.

It exhausts all 33,867 labeled tournaments with root fixed at 0 through order 6, covering every rooted isomorphism type. Dynamic-programming counts are checked against independent literal path/permutation enumeration through order 5 (1,099 cases). Each case checks the SCC recurrence, alpha positivity, both Q floors, the suffix vertex bound, all four coordinates after compression, the elementary representative bound, the Camion core budget, and the sink diagonal. These finite checks support implementation correctness; they are not the all-order proof.
