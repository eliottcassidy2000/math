# Unit binary boundary cancellation forces a physical twin

Session 197. An all-order structural inverse for rank-one clone contexts and exact fixed-mask binary returns. This does not assert scalar completeness or exclude adaptive changes of graph or mask.

## 1. Definitions and statement

For a rooted tournament (T,r), let O=N^+(r), I=N^-(r). Its crossing graph consists of arcs O->I, regarded as bipartite edges. Assume its matching rank is exactly one. Let a_j count collapsed clone-run words with j occurrences of r, no adjacent occurrences, and all exterior vertices appearing once; set a_0=H(T-r). Then only a_1,a_2 can be nonzero at positive indices, a_1=H(T), and a_2>0. Put

    A=a_2, beta=a_1-2a_2,
    R_T(z)=a_0+a_1 z+a_2 z^2.

For positive transitive clone multiplicity m,

    H(T[r:Tr_m])=A*2^m+beta.

Call distinct vertices r,v transitive twins if they have the same orientations to every vertex outside {r,v}. Their mutual arc can have either orientation; the pair is an actual transitive module.

**Theorem.** If beta is +1 or -1, then

    R_T(-1)=0  if and only if r has a transitive twin.   (1)

The implication is physical: a twin can be deleted to halve the clone amplitude while preserving the same unit exponential constant. Zero boundary by itself does not force a twin.

## 2. The rank-one endpoint product

A rank-one bipartite crossing graph has a one-vertex cover c. Taking a converse if necessary, suppose c belongs to O. Set A0=O-{c} and B=I. All arcs between B and A0 point B->A0.

The reviewed rank-one boundary theorem gives

    R_T(-1)=-s_c(T[O]) * sum_(b in B, c->b) e_b(T[B]).  (2)

Here s and e count Hamilton paths with the indicated start and end. One direct explanation is useful: summing all root-gap insertions at z=-1 leaves exterior words starting in I and ending in O. Their internal factor is ordinary adjacency except on an I-to-O step, where it is minus the reverse crossing edge. Since every crossing edge meets c, a surviving word consists of a Hamilton path through B ending at some b with c->b, followed by c and the rest of O. It has one negative crossing and the product factors as (2).

This formula is valid for either available cover center. Its precise dependency is `session197-prime-atoms/rank-one-boundary/PROOF.md`, SHA256 b7935c6f082ff3b6840b71e414a109a15509da3daa8fe91f336796cbac8ebc3b. The structural argument below is additional.

## 3. Vanishing exposes a strong homogeneous block next to r

Every vertex in the source SCC of a tournament occurs as a Hamilton-path start, and no other vertex does. Dually, the possible endpoints are exactly the terminal SCC. For the positive statement within a strong component, a longest path with specified endpoint can always be extended unless the path vertices dominate every outside vertex, which would contradict strongness. Reverse arcs for specified starts; concatenate SCC paths in their transitive order.

Suppose (2) vanishes.

Case I: s_c(T[O])=0. The source SCC C of T[O] does not contain c, so C is a nonempty subset of A0 and

    C -> (O-C),    I -> C,    r -> C.

Therefore {r} union C is a homogeneous module, with internal order r->C. Every vertex of C has the same external incidences as r.

Case II: the second factor in (2) is zero. Let C be the terminal SCC of T[B]. Every vertex of C can be an endpoint of a B-path, so every vertex of C beats c. The other incidences are

    (B-C) -> C,    C -> A0,    C -> r.

Again {r} union C is a homogeneous module, now with internal order C->r, and all its vertices have the same external incidences as r.

These cases may overlap; either choice of C is sufficient. In each case C is a strong tournament. Delete C and call the remaining rooted quotient Q, retaining r. The crossing graph of Q still has rank one: no removed vertex was incident to a crossing edge, and the nonempty crossing edges at c remain.

## 4. The unit constant collapses that block to a vertex

For any rank-one rooted quotient Q and any nonempty tournament D substituted at its root, the path-fragment decomposition gives

    H(Q[D])=a_2(Q) B_D(2)+beta(Q) H(D),             (3)

where B_D(2) counts ordered two-path covers with empty slots allowed. There are at most two D fragments, because two disjoint exterior return fragments would contradict crossing rank one. A single fragment contributes a_1(Q)H(D), while two ordered nonempty fragments contribute a_2(Q)(B_D(2)-2H(D)). This proves (3) for arbitrary D.

Replace r in the original T by Tr_m. The substituted block D is either Tr_m->C or C->Tr_m. Its scalar path count is H(C). Its two-bin cover count is

    B_D(2)=2^m B_C(2).                              (4)

Indeed an ordered two-cover of a transitive join is uniquely the concatenation, in each labelled slot, of the two component covers. Thus these cover counts multiply, and B_(Tr_m)(2)=2^m.

Combining (3)-(4), the original clone sequence has constant

    beta(T)=beta(Q) H(C).                           (5)

If beta(T)=+1 or -1, this forces H(C)=1. A nontrivial strong tournament has a Hamilton path ending at each of at least three vertices, so H(C)>=3. Hence C is a singleton. That singleton is an actual transitive twin of r. The same proof applies after taking a converse when the crossing center was in I.

Conversely, suppose r has a transitive twin v. Contract the pair to a root q of Q. No crossing edge is incident to v, and Q has the same rank-one crossing graph. Replacing r by m clones in T is precisely replacing q by m+1 transitive clones in Q. Hence

    A(T)=2A(Q),  beta(T)=beta(Q),
    H(T-r)=2A(Q)+beta(Q)=A(T)+beta(T).

This is R_T(-1)=0. This converse does not need beta to be a unit.

Deletion of a forced twin therefore gives a genuine smaller graph with the same unit constant and half the amplitude. Repeated deletion terminates because the positive integer amplitude halves at each step. At a terminal graph with no root twin the rank-one boundary is strictly negative by (2). This is a physical reduction of the given context, not a claim that every integer amplitude is realizable.

## 5. Exact fixed-mask binary returns were already literal cloning

Fix any initial tournament S and a subset E of its vertices. Append v_1,v_2,... so that each v_j beats exactly E among the original vertices, loses to all other original vertices, and loses to all earlier newborns. This is a fixed right-birth mask on continuously retained old ports containing E, with every newborn forgotten as a port.

Suppose its scalar sequence, including the original graph at t=0, is not constant and satisfies

    H_(t+1)=2H_t+c for every t>=0, c in {+1,-1}.     (6)

After the first birth, regard v_1 as a clone root. Positive t is exactly transitive cloning of that root. The nonnegative run expansion has largest exponential base nu+1, where nu is the crossing matching rank. The exact matching/run equality, including positivity of the top run coefficient, is the reviewed minimum-cover factorization theorem in `session197-prime-atoms/cover-amplitude-factorization/PROOF.md`, SHA256 2a0451a8473cd17f38a9140c18788fca50325d6f07d8a2716281252f871fbfd6. Nonconstant growth with base two forces nu=1. Solving (6) gives H_t=(H_0+c)2^t-c. The deletion value H_0 therefore equals the clone amplitude plus its exponential constant, so the root-run boundary is zero.

The theorem forces v_1 to have a twin w among the original vertices. Consequently the whole sequence is already literal transitive cloning of this old vertex w in S. If v_1 beats w, the new clones insert immediately before w; otherwise they append after w. Every old vertex keeps its physical identity and all original arcs stay fixed. No forgotten vertex is re-exposed by this statement: the given fixed mask itself performs the cloning.

Thus a fixed-mask exact unit binary return cannot conceal a genuinely different physical source behind an artificial zero-clone boundary. This does not cover a first transient, changing masks or ports, multiple births per iteration, opposite-side mixtures, count-preserving graph rebuilding, or a nonunit offset.

## 6. Necessary scope and verification

The actual source

    rows [2,32,11,19,7,29], root 0,
    (a_0,a_1,a_2)=(15,27,12)

has rank one and zero boundary but no transitive twin at the root. Its exponential constant is +3. It is a counterexample to dropping the unit assumption, and its root deletion followed by fixed-mask births gives 15,27,51,... with affine offset -3.

The included verifier checks the endpoint product, the exposed SCC block and divisibility (5), actual twins for unit states, and literal small root-gap evaluations. The inherited representative controls through order seven supply 1,053 rank-one roots, 321 zero-boundary roots and 169 unit zero-boundary roots. Completeness of this finite bank is not used to prove (1); the structural argument establishes all orders. Python normal and optimized modes must reproduce the same receipt.

The verifier imports only the included standard-library graph core. The rank-one endpoint theorem and the physical clone-run interpretation are the inherited mechanisms; the unit-forced twin, deletion reduction, and fixed-mask consequence are the new structural implications here. No uniform target-spectrum or inverse-coverage claim is made.
