# Exact equality at the physical Z-defect floor

Status: PROVED / independently audited, with the explicit singleton-boundary correction.

## 1. Scope and definitions

For a nonempty tournament U rooted at v write h=H(U), s=start_v(U), e=end_v(U), and d=H(U-v), with H(empty)=1. A directed path is counted by its ordered vertex sequence. Let L and R count ordered spanning two-path covers whose designated first path is nonempty and respectively ends or starts at v. The second path may be empty. Let b count two labeled spanning paths of U-v, both allowed empty, and put

    B=h+L-s,   D=R+b-e-B.

The physical move Z adds z with z->v and every vertex of U-v beating z, retaining root v. For nontrivial U its exact inherited count identity is B(ZU)=B(U)+R+b-e, so D=B(ZU)-2B(U). For the singleton U, the formula-defined D is zero whereas B(ZU)-2B(U)=1; the identification with the physical difference is not asserted at that boundary. The relevant all-input transfer is published in `session203/rooted-cover-closure/PROOF.md`. The argument below reproves all inequalities it uses and does not assume that a numerical state is physically realizable.

**Main theorem.** For every positive integer t, the following are equivalent:

1. A nontrivial positive-start rooted tournament has s=t and D=t.
2. A strong rooted tournament has s=t and D=t.
3. Some nonempty tournament T has H(T)=t.

This is an exact all-order characterization of attainable equality values, not a characterization of every equality graph up to isomorphism. An exact test on individual graphs is given in Section 3. Neither F(8,137) nor an eventual odd-spectrum threshold is settled.

Inheritance board: the physical Z transfer; the nonnegative cover defects; the exact bounded profile census; the moving-gate start-two grammar; and the ordinary Hamilton-path spectrum. The tempting false implication L_a(A)=H(A) => start_a(A)=1 is explicitly avoided: the correct conclusion is start_a(A)=H(A-a). The latter may exceed one.

## 2. Two elementary cover identities

Let A be any nonempty tournament rooted at a, and write h_A=H(A), d_A=H(A-a), and s_A=start_a(A). Cutting a Hamilton path immediately after a injects h_A paths into the L_a(A) covers. Taking the designated path to be the singleton a gives d_A more covers. The intersection consists exactly of the s_A Hamilton paths starting at a. Consequently

    L_a(A) >= h_A+d_A-s_A.                     (1)

Deleting the initial a injects its starting Hamilton paths into paths of A-a, so s_A<=d_A, including a singleton A. Therefore

    L_a(A)=H(A)  =>  start_a(A)=H(A-a).         (2)

For a nontrivial rooted U let A=U-v. Deleting the terminal root from an L cover gives a labeled two-cover (Q,P) of A for which Q is empty or its last vertex beats v. This is a bijection with those covers. Its complement has a nonempty first path ending at an outneighbor of v. Thus exactly

    b-L = sum_{a: v->a} L_a(A).                (3)

Every summand is at least H(A)=d by cutting after a. If s>0 and U is nontrivial, v has at least one outneighbor. Hence b-L>=d.

Set p=d-e and f=R-h+e-d. Root deletion gives p>=0. Cutting each Hamilton path immediately before v and separately taking singleton designated root paths gives two subfamilies of R of sizes h and d; their intersection has size e. Thus f>=0. More precisely, f counts those covers (Q,P) in R for which Q has at least two vertices, P is nonempty, and v beats the last vertex of P.

An exact algebraic rearrangement now gives

    D-s = f+2p+(b-L-d) >= 0.                  (4)

The positive-start and nontrivial qualifications are essential. For the singleton, s=1 but D=0.

## 3. Exact equality test and necessity

**Equality test.** For a nontrivial positive-start rooted U, equality D=s holds if and only if all three conditions below hold:

- v has exactly one outneighbor a;
- L_a(U-v)=H(U-v);
- end_a(U-v)=0.

Indeed, (4) forces every nonnegative summand to vanish. By (3) and d>=1, b-L=d forces exactly one outneighbor a, and its L count equals d. If a is the unique outneighbor, then

    s=start_a(U-v),   p=end_a(U-v).            (5)

The second identity follows because Hamilton paths of U-v fail to extend by v exactly when they end at a. Also f=0 automatically: any nonsingleton Q starting at v contains a, whereas a cover counted by f would require the disjoint path P to end at a. This proves both directions of the equality test.

Suppose D=s. With A=U-v, identity (2) and (5) yield

    s=start_a(A)=H(A-a)=H(U-{v,a}).             (6)

The residual tournament is nonempty. Otherwise U is the rooted directed edge v->a, and end_a(A)=1 contradicts equality. Thus any positive equality value s is the Hamilton-path count of a nonempty tournament. This proves (1) implies (3); (2) implies (1) is immediate because a strong equality source cannot be a singleton.

## 4. A physical strong source for every attainable value

Conversely take any nonempty tournament T, put t=H(T), and let c=C(T) count its two labeled spanning paths, allowing empty bins. Add v,a and orient

    v->a,   a->T,   T->v,

with the orientation inside T unchanged. Root the resulting U at v. It is strong: v reaches every vertex through a, every vertex reaches v directly, and these routes also connect a and all old vertices. Direct path and cover decomposition gives

    h=t+c,   s=e=d=t,
    L=2c-t,  R=t+c,  b=2c,                    (7)
    B=3c-t,  D=t.

For completeness, Hamilton paths starting at v or a contribute t each. Those starting in T consist of a nonempty first T path, then v,a, and an optional second T path; they contribute c-t. For L, the singleton v contributes t. If the designated path is nonsingleton, a is either in that path before a nonempty T segment or at the start of the other path; each case contributes c-t. For R, the singleton v contributes t and a nonsingleton path starting v,a contributes c. In a two-cover of a->T the vertex a starts either labeled bin, giving b=2c. The s,e,d counts follow directly. These exhaust disjoint cases, proving (7).

In particular U is a strong physical equality source. This proves (3) implies (2) and completes the theorem.

## 5. Consequences and limits

The lower bound D>=s is sharp at exactly the ordinary Hamilton-path spectrum. Declaring equality freely attainable at every requested start value is therefore an additional spectrum claim, not a consequence of the cover cone. Conversely, absence of an ordinary count t excludes equality D=t at start t at every source order.

The construction supplies B=3C(T)-t. Its inherited four-core readout is 3B+2s=9C(T)-t. This is a formula on actual sources, not an all-target coverage theorem. The theorem does not assert that every equality graph has the cyclic form in Section 4, or that the source can always be reduced while preserving B. The reduction (6) preserves the start count as an ordinary Hamilton count and deliberately forgets B and the remaining cover coordinates.

For example, T a directed triangle has t=3 and c=12, so the cyclic construction has L=21 and h=15, while its root has s=D=3. A separate illustration of the incorrect start-one inference is obtained by taking A=a->T: L_a(A)=H(A)=3 but start_a(A)=3, not one.

`verify.cpp` enumerates every labeled tournament through order 6, every root, and checks (1)-(6), including the singleton boundary and the exact individual-graph equivalence. It also physically constructs the cyclic source for every such nonempty T and checks every coordinate in (7). Physical Z children are directly recounted for all positive-start rooted inputs through order 5. These are finite controls of the all-order proof, not a substitute for it.
