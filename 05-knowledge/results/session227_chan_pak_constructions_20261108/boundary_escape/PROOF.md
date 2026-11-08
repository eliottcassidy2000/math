# Exact graft escape from the three terminal faces

**Status: PROVED, independently audited PASS, with FINITE-EXACT controls.** No scalar-spectrum separation, eventual coverage threshold, or Chan–Pak proof is claimed. No literature-priority claim is made.

## 1. Inheritance, notation, and the loss boundary

For a nonempty rooted tournament `(A,a)`, write

    D = H(A), E = # Hamilton paths ending at a,
    U = # ordered spanning two-path covers (P,R), R nonempty starting at a,
    d = H(A-a), S = # Hamilton paths starting at a.

The unmarked bin P may be empty; an empty path counts once, and H(empty)=1. Put

    f = U-D-d+E.

The Session225 graft recovery proves the four-coordinate operators and the general graft law. Session198, `session198/bilinear-descent/PROOF.md`, Sections 1–3, already proves the positive source grammar, its exact inverse, and the terminal-face condition for **start-one sources**. Its faces are E=0, d-E=0, f=0. THM-4875(A1), `01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`, explicitly attributes to THM-4873(ii)(D) the inequality

    U-D >= (d-E) + 2 sum_{ {u,v} subset N+(a) } H(A-{a,u,v}).       (1)

Thus the defect inequality below is **inherited**, with a self-contained combinatorial proof supplied for the exact application. It is not a new endpoint bound. The new conclusion here is the exact classification of which arbitrary rooted grafts escape all three faces, the resulting separation from every nonempty O/Z/J history on every physical predecessor, and minimality of a two-vertex gadget for this escape.

Closest mechanism: the inherited face condition and graft closure. Hostile: two different rooted states can have the same Hamilton count. Corrected near miss: full rank is sufficient for universal operator separation, but is unnecessary for individual-state escape. Least-used sidecar: which end of a directed-edge gadget is marked. Concept board: missing-cover defect; sink tests; terminal faces; two-vertex delayed rooting; root-start multiplicativity; the integer inverse on the edge-graft image plane.

## 2. The defect counts missing joins

**Lemma (exact interpretation, with inherited lower bound).** Let O=N+(a). Then f counts ordered spanning covers `(P,Q)` of A-a by two nonempty paths such that the last vertex of P and the first vertex of Q both belong to O. In particular

    f >= 2 sum_{ {u,v} subset O } H(A-{a,u,v})
      >= |O|(|O|-1),                                        (2)

and

    f=0 if and only if outdegree_A(a)<=1.                    (3)

**Proof of the interpretation.** The U covers with reserved path R={a} number d. Every remaining cover has R=(a,Q), where Q is nonempty and starts in O. Such a cover reconstructs a Hamilton path by joining P to a,Q exactly when P is empty or its last vertex beats a. The number of these joinable covers is D-E: the omitted E Hamilton paths are exactly those ending at a. The nonjoinable covers therefore number U-d-(D-E)=f, and are precisely the claimed covers.

For a self-contained proof of (2), fix an unordered pair {u,v} in O and a Hamilton word `w_1,...,w_k` of A-{a,u,v}. For 1<=i<=k put

    a_i=[w_i -> u], b_i=[w_i -> v],
    a_0=b_0=1, a_{k+1}=b_{k+1}=0.

At cut i, the paths

    (w_1,...,w_i,u), (v,w_{i+1},...,w_k)

are valid exactly when `a_i(1-b_{i+1})=1`. Swapping u and v gives another valid cover exactly when `b_i(1-a_{i+1})=1`. Empty prefixes and suffixes are allowed; the two resulting bins themselves remain nonempty. For each cut,

    a_i(1-b_{i+1}) + b_i(1-a_{i+1})
       >= a_i+b_i-a_{i+1}-b_{i+1}.

Summing over i=0,...,k gives at least two covers. This also works at k=0. A produced cover recovers its endpoint pair, the cut, and the word obtained by deleting its distinguished last/first vertices and concatenating the remainders. Hence different choices do not collide. Summing proves the first inequality. Every tournament, including the empty one under our convention, has a Hamilton path; this proves the second. If |O|<=1, the exact interpretation permits no cover because its two distinguished vertices are distinct. This proves (3). ∎

We also use the elementary physical inequalities E<=d and D-E<=U-d. Both follow from deletion and the same joining injection. Moreover D-E>0 iff a is nonsink. For the nontrivial direction, choose a Hamilton path; if it ends at a, move a immediately before the first vertex on that path which a beats. The resulting path is valid and does not end at a. Also U-d>0 iff a is nonsink: use the reserved path (a,u), where a->u, and any Hamilton path on the complement.

## 3. The face condition has no source restriction

The inherited operators are

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

On **every physical input**, not only a start-one source, O has output E'=0, Z has output E'=d', and J has output f'=0. Thus the output of every nonempty word lies in

    B = {E=0} union {d-E=0} union {f=0}.                     (4)

The final letter alone proves this. There is no bound on word length, no fixed seed bank, and no requirement that the word be chosen independently of its input. The empty word is excluded: any state can of course be its own empty-word seed.

## 4. Exact escape criterion for arbitrary grafts

Let B be a rooted tournament `(B,r)`, with at least two vertices. Form G=G_B(A,a) by keeping both induced tournaments, putting B->a and A-a->B, and marking r. Let q=H(B-r), e count Hamilton paths of B ending at r, and delta count ordered two-bin spanning covers of B with one specified bin nonempty and ending at r. Its other bin may be empty.

**Theorem.** For every nonempty rooted A and every rooted B with |B|>=2:

1. `d_G-E_G >= d*q > 0`.
2. `E_G>0` if and only if a is nonsink in A.
3. `f_G>0` if and only if r is nonsink in B.

Consequently G lies strictly outside B in (4) if and only if **both marked vertices are nonsinks in their respective input tournaments**. In that case its four-coordinate state is not the output of any nonempty O/Z/J word from any physical rooted tournament.

**Proof.** For (1), choose any Hamilton path of A-a and any Hamilton path of B-r, and concatenate them followed by a. They form a path of G-r, and there are d*q choices, even when A-a is empty. Each ends at a, which loses to r, so none can be extended by appending r. Deleting r injects all E_G paths ending at r into the remaining deletion paths. This proves (1).

For (2), set p=D-E and t=U-d. The inherited graft law gives

    E_G = e*p + (delta-e)*t.

There are delta>=q covers with singleton marked bin r. Also delta-e>=0, since the e Hamilton paths occur among the delta covers by using an empty other bin. If a is nonsink, t>=p>0, whence E_G>=delta*p>0. If a is a sink, p=t=0, giving E_G=0.

For (3), the outgoing neighbors of the new root are its outgoing neighbors inside B together with the old root a. It loses to every old nonroot. Thus its outdegree is `1+outdegree_B(r)`. Apply (3). The three conclusions and (4) prove the equivalence and separation. ∎

The singleton-gadget exception is explicit: it is J, which always has f'=0. The theorem makes no claim that all points outside B are physically realizable, or that an escaping state has a new scalar Hamilton count.

## 5. The smallest escape: mark the source of an edge

Take B to be the directed edge with outgoing rows `[2,0]`, marked at its source 0. Its matrix is

    L = [1 -1 3 -2]
        [0  0 1 -1]
        [2 -1 4 -1]
        [1 -1 1  0].                                      (5)

This is obtained from the inherited graft law with

    (H,Q,s,e,alpha,delta,T,q,b)=(1,4,1,0,2,1,3,1,2).

Directly,

    E'=U-d, d'-E'=D-E+d, f'=U+E.                          (6)

For a nonsink source all three are strictly positive. Hence two vertices suffice. One vertex never suffices, since its graft is J; this proves minimality **in gadget order for this state-escape property**, not minimal order for any scalar count.

The other rooted directed edge, marked at its sink, is exactly the word ZJ (words act left to right). For L, apply two Z births retaining the old root, then mark the first of the two newcomers. Thus the choice of marked endpoint, not determinant rank, causes the first escape.

Four vertices are also the smallest order of any physical state outside (4). At order at most three, a root of outdegree at most one has f=0; any other root is a global source and has E=0. The four-vertex t.L below attains strict escape.

The matrix L has rank 3. Its image is exactly the rational hyperplane

    D=2E+d,                                               (7)

and its kernel is spanned by `(-2,-1,1,1)`. For example its minor on rows 1,2,3 and columns D,E,U is -1, proving rank at least three; (7) proves rank at most three. Individual-state separation therefore does not require the full-rank four-vertex gadget from Session225.

### Explicit collision and first closing witness

Start from t=[2,0], rooted at 0, with state `(1,0,2,1)`. Then

    t.L = (5,1,9,3), while t.OJ = (5,1,5,1).

They have the same scalar H=5; the first lies outside (4), the second lies inside. This is the required hostile to scalar-spectrum separation.

An additional J applied to t.L gives a five-vertex strong tournament with rows

    [18,28,9,17,4], root 4,

and H=13. The intermediate t.L has rows `[2,12,9,1]`, root 2. These rows specify every arc and are checked by literal subset DP and permutation enumeration. H=13 is a known attainable count; this is a concrete witness, not a new spectrum value.

## 6. A multiplicative sidecar and an exact three-coordinate inverse

**Root-start law.** If |A|>=2, then for every nonempty gadget B,

    S(G_B(A,a),r) = alpha_B * S(A,a),                      (8)

where alpha_B counts ordered two-bin covers of B whose first bin is nonempty and starts at r. If A is a singleton, instead `S(G_B(A,a),r)=s_B`, where s_B is the number of Hamilton paths of B starting at r.

**Proof.** A root-start Hamilton path begins in B. Its only possible transition into A enters a, and it cannot return into A again. Thus its old vertices form one Hamilton path of A starting at a. When |A|>=2 this path ends at an old nonroot, which may be followed by an optional B suffix. Splitting B into its nonempty r-starting prefix and optional suffix gives exactly alpha_B choices, independently of the S(A,a) choices. For singleton A, the path ends at a and the B suffix is forbidden, leaving s_B choices. ∎

All O/Z/J moves preserve S on nonsingleton rooted tournaments, by the same path insertion argument; this preservation is already present for S=1 in Session198. In particular L doubles S, whereas the Session225 four-vertex gadget multiplies it by 5. This is a structural sidecar, not a claim that a fifth coordinate is necessary to predict future graft Hamilton counts.

For the L image plane (7), use coordinates `(E,U,d)`. The induced matrix and its inverse are

    M = [0 1 -1]         M^-1 = [-3  2 -5]
        [3 4  1]                [ 2 -1  3]
        [1 1  1],               [ 1 -1  3].               (9)

They are integral and det(M)=-1. Therefore the rank-three graft has an exact unimodular inverse **on its image plane**. These linear identities do not imply that an arbitrary nonnegative candidate has a physical inverse, or that this plane is covered by physical states.

The seed t itself lies in (7). For `x_m=L^m t`, let h_m be its D-coordinate, and let n_m be the Hamilton count after an extra closing J. All these graphs are explicit. The root-start count is `S(x_m)=2^m`. The characteristic identity

    M^3-5M^2+M+I=0

gives, for every m>=0,

    h_{m+3}=5h_{m+2}-h_{m+1}-h_m,
    n_{m+3}=5n_{m+2}-n_{m+1}-n_m,

with

    h_0,h_1,h_2 = 1,5,25,
    n_0,n_1,n_2 = 3,13,61.

Thus h begins `1,5,25,119,565,2681,12721`, and n begins `3,13,61,289,1371,6505,30865`. Every x_m with m>=1, and every closed x_m.J with m>=0, is strong: the initial root reaches all old vertices, and the graft strongness criterion from Session226 applies at every step. Alternatively, use the old nonroot -> B -> old root paths directly. Every x_m with m>=1 escapes (4), while its closing J returns to the f=0 face. These explicit families are not asserted to be new to the literature or to cover an arithmetic progression.

## 7. Verification and open boundary

Run `python3 verify.py` and `python3 -O verify.py`. Their deterministic JSON receipts must agree byte for byte. The checker uses only the Python standard library and explicit exceptions. It does not use the exploratory scripts in this directory or any external data.

The exact finite universe includes all 1,099 labelled tournaments through order 5 with root fixed at label 0; fixing the root loses no rooted isomorphism type. Their profiles and the defect-cover interpretation are independently counted by permutations and subset DP. All literal O/Z/J transitions are checked. Every pair of labelled rooted tournaments through order 4 (75 sources by 75 gadgets, 5,625 pairs) is grafted literally; the general operator, exact escape conditions, strict deletion gap, and root-start law are checked, including singleton and sink boundaries. The edge matrix, its other rooting, image plane, inverse, recurrence, scalar collision and H=13 witness receive separate controls. Larger recurrence terms are algebraic consequences, not claimed literal enumerations of large graphs.

The closest existing C1/root-change censuses already record the familiar 29 missing scalar values in their bounded ranges. Exploratory rediscovery of that list is not a new result and is not used here. Adding edge and triangle grafts gives promising small finite repair witnesses, but no inverse-coverage theorem for prescribed integers follows. Chan–Pak completeness remains OPEN.
