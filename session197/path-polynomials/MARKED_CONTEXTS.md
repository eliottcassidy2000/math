# Positive whole-module transport of reserved-bin endpoint sidecars

Session 197 extension of `UNIVERSAL_CONTEXTS.md`. **All-order marked-profile matrix positivity and explicit physical rooted transfers.** No literature-priority claim is made. This uses the S193 sidecar formulas and S195 fragment kernel already cited in the companion packet. It does not close every multiport interface: the hypotheses below are essential.

## 1. Exactly which marked profiles are covered

Fix k>=1 distinct named vertices a_1,...,a_k of T and k distinct labelled bins 1,...,k. The k=0 case is the unmarked theorem in the companion packet. Give each bin exactly one requirement: its path either starts at its named vertex or ends at its named vertex. No bin has both endpoints prescribed. All k bins are nonempty and are distinct path components. Other bins may be empty.

Let M_(T,a,rho)(d) count such covers in d labelled bins, and use the reserved-bin shift

    f_(T,a,rho)(x)=M_(T,a,rho)(x+k).

The number of ordinary, unreserved bins is x. Then

    f(x)=sum_(P_1,...,P_k) B_(T - union P_i)(x),                  (1)

where the sum is over ordered, vertex-disjoint actual paths satisfying the specified endpoint requirements. It follows that f has nonnegative integer coordinates in the extended basis {1,q_0,r_0,q_1,r_1,...}. Its constant coordinate f(0) is the number of spanning ordered k-path covers satisfying the requirements. For k=1 this is the root's Hamilton start or end count.

The shift matters. A C3 root has unshifted start polynomial d^2-d+1, whose q_0 coefficient is negative. After reserving its bin, f(x)=x^2+x+1 has nonnegative coordinates. Thus the unshifted untyped cone cannot simply be reused.

The S193 L and R profiles are the k=1 cases. Its mixed K_(q,u), where the first path ends at q and the second starts at u with q!=u, is the k=2 case. Profiles with two endpoints on the same path, shared named anchors, or coincident reserved bins are not covered by this theorem.

## 2. A fixed context acts on the entire marked polynomial

Let (Q,r) be any rooted tournament and replace r by T as one intact module. Keep all named anchors inside T and retain their bin labels and endpoint roles. There is a linear operator A_(Q,rho) on Q[x] with

    f_(Q[r:T],a,rho)=A_(Q,rho) f_(T,a,rho).                     (2)

It depends on Q and the role pattern, not on the internal tournament T or the anchor identities.

**Full-profile linearity.** Cut every global path into maximal T-runs and outside runs. Each reserved path has exactly one distinguished local T-fragment containing its anchor at the prescribed global endpoint. These k fragments are distinct, and each must remain first or last in its own reserved output bin. Every other local T-fragment is unmarked. If there are j unmarked fragments, all possible recombinations depend only on j, the k endpoint roles, and Q, because every cross incidence is uniform on T. Therefore the kernel is linear on the falling-factorial coordinates (x)_j of the shifted marked polynomial. This is full-polynomial transport, not only a scalar endpoint count.

To see uniqueness and obtain a monomial test family, put all end-anchors first in a transitive module, then a free m-vertex transitive chain, then all start-anchors. Other anchors cannot be used by a reserved path because they are reserved to different bins. An end-anchor has no possible predecessor from the free chain; a start-anchor has no possible successor in it. Thus every reserved path inside this test module is its singleton anchor, and

    f_test(x)=x^m.                                             (3)

This works for every m>=0 and every fixed role pattern. The test modules are actual tournaments, not a formal assumption about positive-state realizability. Their monomials span Q[x], so equality on them determines the entire marked operator, including its action on nonphysical basis vectors.

## 3. Positive open-path and odd-cycle skeletons

Write S for x->x+1, X for multiplication by x, and retain

    delta=S-S^(-1),       Delta=S+S^(-1)-2I.

There is one additional positive elementary operator:

    nabla=S-I=(delta+Delta)/2.                                 (4)

The factor 1/2 introduces no fractional matrix entries. Every coefficient in the explicit delta/Delta basis formulas in `UNIVERSAL_CONTEXTS.md` is even, including their constant boundary cases. Hence nabla has nonnegative integer entries. Also nabla(1)=0. The shift S=I+nabla is positive as well.

**Theorem (marked universal positivity).** The matrix of every A_(Q,rho) in {1,q_0,r_0,q_1,r_1,...} is entrywise a nonnegative integer matrix.

**Proof.** Apply (1) to the actual monomial test module (3), after substituting it into Q. Choose the k distinguished global paths first. The remaining vertices are counted by the odd-cycle bin formula, exactly as in the unmarked proof.

Delete every maximal run of free-chain vertices from these distinguished paths and from the crossing odd cycles. Retain each named anchor at the appropriate endpoint of its own labelled path. A distinguished path may be its singleton anchor. The transitive test ordering forbids a free-chain run directly adjacent to that anchor; every such run position in the skeleton is separated from it by exterior vertices. All marked paths, cycle skeletons, and unused exterior vertices are disjoint outside the anchors. Exterior-only cycles must be odd. Canonical rotations and fixed bin labels remove all cyclic and path-label ambiguity.

Suppose there are t_0 marked free-chain run positions in the distinguished open paths. Assigning nonempty sets of free-chain vertices to these distinct positions contributes

    (S-I)^t_0 = nabla^t_0.

There is no parity constraint on a distinguished open path. Each crossing cycle with t_i runs and l_i exterior vertices contributes the same kernel as before,

    L_(t_i,epsilon_i)=(S-I)^t_i+(-1)^epsilon_i(S^(-1)-I)^t_i,
    epsilon_i=1-l_i (mod 2).

Its odd-cycle factor two cancels the parity projector's half. If c_0 cycles stay outside, c_1 cycles cross the free chain, and u exterior vertices remain unused by the selected paths and cycles, the skeleton contributes

    2^c_0 X^(u+c_0+c_1) nabla^t_0 product_i L_(t_i,epsilon_i).   (5)

Again X acts after the shifts. Vertices of distinguished paths, including their fixed anchors, contribute no X factor, because their bins were already reserved. Unused free-chain vertices choose among the x ordinary bins.

The run positions are distinct and their assigned named vertices have unique increasing orders. Cutting and reinserting the free-chain runs are inverse operations. There is no factorial division or additional ordering of reserved paths: those paths already have their specified labels. Thus summing (5) gives A_(Q,rho) on every x^m. The full-profile linearity in Section 2 extends this to all Q[x].

Each L is a nonnegative integer polynomial in Delta, or delta times such a polynomial, by the power-sum recurrence in the unmarked proof. Together with positive X, nabla, delta and Delta, (5) proves matrix positivity. QED.

This theorem does not give arbitrary positive matrices a graph realization. It also does not transport a root to a new vertex or split the uniform exterior neighborhood of T without an additional kernel. For a bin with both its start and end prescribed, cutting may create two distinguished T-fragments rather than one; extra split-fragment states are then required.

## 4. An explicit closed physical transfer for all role counts

Let W(T)=C3[T,K1,K1], written a->T->b->a. Suppose s of the reserved bins prescribe starts and e prescribe ends, with k=s+e and distinct anchors. Then the full shifted marked polynomial transforms by

    A_W^(s,e)=D+kX+kS+se I,                                   (6)

where D=X^2+X delta is the unmarked wrapper operator. In particular,

    one retained root, start or end: D+X+S;
    two bins, one end and one start: D+2X+2S+I;
    two starts or two ends:          D+2X+2S.

**Proof on the physical monomial test family.** If neither added vertex is used by a reserved path, the remaining bin profile is D(x^m). If exactly one added vertex is used in one reserved path, it must be b after a start-anchor or a before an end-anchor. There are k choices and the remaining graph gives X x^m. If both added vertices are in a single reserved path, there are k choices; a possibly empty increasing free-chain subpath can be selected at its free endpoint, giving S x^m. The only possible split between two reserved paths uses a before one end-anchor and b after one start-anchor, giving se choices and x^m. These cases are disjoint and exhaustive. Apply full-profile linearity. QED.

All four operators in (6) are positive in the shifted coordinates. This is an actual retained-root transfer, not merely a positive expansion of the initial sidecar.

At x=0, D and X vanish, giving the exact spanning marked-cover update

    f_new(0)=k f(1)+se f(0).                                   (7)

For a retained start/end root, its new Hamilton endpoint count is the previous two-bin rooted cover response. For the mixed K profile the new spanning mixed two-cover count is 2K_T(3)+K_T(2). These formulas do not assert that a fixed finite number of scalar evaluations closes further iterations: S continues to expose higher bin values.

## 5. New apex roots and the limit of an immediate return

For a nonempty input T with unmarked profile f=B_T, the new apex start profiles, shifted by one reserved bin, are

    f_start_a = Xf+(X+I)nabla f,
    f_start_b = Xf+Sf.                                        (8)

The end profiles are interchanged. For example a path starting at a is either just a, or traverses a nonempty T-path, possibly followed by b; the remaining vertices give (8). A path starting at b either is b alone or proceeds through a and a possibly empty T-path. Both root endpoint counts are H(T), since f(0)=0.

These positive formulas clarify a physical limitation: attaching W and choosing one of its new apex roots does not restore a unit endpoint to a nontransitive input. Nor does the shifted cone by itself supply a literal returning ternary-digit operation. Endpoint constraints and the complete marked state still matter. The mixed-endpoint fragment interfaces in the S195 algebra require additional coordinates beyond this one-endpoint-per-bin theorem.

## 6. Verification

Run `python3 verify_marked_contexts.py` and the optimized version, comparing their stdout files. The checker independently computes marked profiles by choosing disjoint anchor-constrained path subsets and summing the ordinary bin profiles of the induced remainders. It checks the monomial-test-module continuation against actual input contexts, positivity on formal q/r vectors, the explicit W transfers, and the new-apex formulas. The emitted JSON specifies the exact finite universe. These computations control conventions and multiplicities; the all-order result follows from (5), not finite extrapolation.
