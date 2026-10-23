# Exact incidence collisions survive padding and strong cyclic lifting

Status: PROVED transformations and counting formulas; the named order15 seed is an independently checked finite certificate. This gives an all-size boundary to the successful order8 signature, not a refutation of authenticated port reconstruction.

## 1. The signature and its recoverable scalar

For a tournament T, let c_v count cyclic triples through v and c_uv count cyclic triples through the unordered pair{u,v}. Define

E(T)=multiset_v multiset_(v->w)(c_w,c_vw).

In each cyclic triple through v there is exactly one outgoing edge from v. Consequently

c_v=sum_(v->w)c_vw.

Thus E retains the multiset of all c_v, despite not explicitly storing the root's c_v alongside each row. It also retains every vertex's outdegree as the length of its row.

The current incidence packet proves FINITE-EXACT separation of all6880 ordinary order8 classes by E, and by its residue version E_(4,3). The exact E field is nevertheless not a complete invariant at every order.

## 2. Ordered source/sink padding

Adding a dominating source creates no new cyclic triangle. Every old E row is unchanged. The new source row is the multiset{(c_v,0):v∈T}, reconstructed from the original E rows by their second-coordinate sums. Therefore equal E signatures remain equal after adding a source. H is unchanged.

Adding a sink similarly appends(0,0) to every old row and adds an empty row, again preserving equality and H. Iterating proves that any H-distinguishing E collision propagates to every larger order, without claiming the padded tournaments remain strong.

## 3. A strong lift at every added size at least two

Let r>=1 and let L_r(T)=C3[T,{a},P_r], with block directions

T -> a -> P_r -> T,

where P_r is the transitive tournament of size r. Write n=|T|. Every L_r(T) is strong for nonempty T, even if T itself is not strong.

THEOREM. E(T)=E(U) implies E(L_r(T))=E(L_r(U)) for every r>=1.

Proof by explicit update. Old internal edges retain c_vw; each old vertex gains exactly r cyclic triples, one using a and each P_r vertex. Its outgoing row is obtained by shifting every first coordinate by r and appending(nr,r), corresponding to its edge to a.

Vertex a has c_a=nr and r outgoing edges, each with tag(n,n). Every P_r vertex has c=n. If it has j later vertices in its transitive block, its outgoing row consists of the global multiset{(c_v+r,1):v∈T} and j copies of(n,0). This global multiset is already determined by E(T). No other cyclic triples occur: a or a P_r vertex meets each pair of old vertices uniformly, and P_r itself is transitive. These formulas determine the whole new E signature. QED.

This closure is proved for exact E equality. It does not assert that an arbitrary residue E_(p,q) alone supports the same update without additional coordinates.

## 4. Exact Hamiltonian response

Let F2(T)=2pc_T(2), where pc_T(2) counts unordered spanning covers by two nonempty directed paths, allowing singleton paths. Then

H(L_r(T))=(2^r+1)H(T)+(2^r−1)F2(T).

Proof. The inherited cyclic substitution kernel (THM-3121) allows only one or two T-runs and one or two P_r-runs, because the a block has one run. The one-run path factor of P_r is1, and its ordered two-path cover count is2^r−2: choose the nonempty proper vertex subset of the first increasing path. The kernel's all-one triple has coefficient3; the other three run choices each have coefficient1. Hence

3H + H(2^r−2) + F2 + F2(2^r−2),

which is the formula. Equivalently the first r1 case is3H+F2. This is a consequence of the already proved path-cover interface, not a new general substitution identity.

For two seeds with differences ΔH,ΔF2, the lifted difference is

2^r(ΔH+ΔF2)+(ΔH−ΔF2).

If ΔH≠0 it vanishes at most once as r ranges over positive integers. Thus every H-distinguishing exact-E collision yields strong collisions at all but at most one order n+2,n+3,... . The explicit seed below has no exceptional r.

## 5. Explicit regular strong seed and every-order consequence

The incidence packet and independent audit give two tournaments on Z/15 with arcs x->x+d for

S1={14,13,12,11,5,6,8},
S2={14,13,12,4,10,9,8}.

Each contains step14, hence a directed Hamiltonian cycle and is strong. Both are regular of score7. Every vertex has c_v=28 and outgoing edge-incidence multiset(3,3,3,4,4,5,6). They therefore share exact E and every residue version E_(p,q).

Independently certified values, in the stated S1,S2 order, are

H1=196264905, H2=196277235,
F2_1=4857720210, F2_2=4857619470.

Thus ΔH=H2−H1=12330 and ΔF2=−100740. The strong-lift difference is

−88410·2^r+113070<0  for every r>=1.

Consequently:

- there are exact-E/every-E_(p,q) collisions with different H at every order n>=15 by source padding;
- there are STRONG such collisions at order15 and at every order n>=17 by the cyclic lift;
- this argument leaves the strong order16 case unspecified.

At r1 (order17) the two lifted counts are5446514925 and5446451175, differing by63750. The sign reverses relative to the seed counts; preserving the signature does not preserve which seed has larger Hamiltonian count.

These are equality-of-signature counterexamples, not alternative full metagraph wirings at the same order. They prohibit treating E's order8 class separation as an all-order classification. Capacity and authenticated path data may still supply more information.
