# Pair-child separation in staircase prime quotients

## 1. Explicit quotient and all-order primeness

For r>=1, let U_r have physical quotient labels A_0,...,A_(r-1), B_0,...,B_r. Orient

    A_j -> A_i for i<j,
    B_i -> B_j for i<j,
    A_i -> B_j exactly when j<=i.

Thus, with numeric labels A_i=i and B_j=r+j, U_r is the cyclic tournament on 2r+1 names with all edges inside its first r names reversed. This equivalence follows directly from forward cyclic distances; the staircase definition above fixes our convention completely.

U_r is strong: A_(r-1),...,A_0,B_0,...,B_r,A_(r-1) is a directed spanning cycle. It is prime, as follows.

Suppose a module M contains two names. If it contains two A-names A_i,A_j with i<j, B_(i+1) distinguishes them and must belong to M. If it contains two B-names B_s,B_t with s<t, A_s distinguishes them and must belong to M; here s<=r-1. Therefore M contains a mixed pair A_i,B_j.

This mixed pair forces B_0 and B_r into M. For j>0, B_0 beats B_j but loses to A_i; when j=0 it is already present. For j<r, B_r loses to B_j but beats A_i; when j=r it is already present. Every A-name sees B_0 and B_r oppositely, so all A-names belong to M. For each intermediate B_j, 1<=j<=r-1, its incidences to A_0 and A_(r-1) are opposite, so it too belongs to M. There are no intermediate B-names at r=1. Thus M is the whole quotient.

## 2. Explicit modules after every deletion

The following list retains physical quotient names and directions.

- Delete A_i, 0<=i<=r-1. The pair B_i -> B_(i+1) is a module. No remaining B-name lies between them, and the only A-name that distinguishes their staircase incidences is the deleted A_i.
- Delete B_j, 1<=j<=r-1. The pair A_j -> A_(j-1) is a module. No remaining A-name lies between them, and the only B-name distinguishing their incidences is the deleted B_j.
- Delete B_0. The remaining A_0 is a universal sink, exposing the global cut (all other names) -> A_0.
- Delete B_r. The remaining A_(r-1) is a universal source, exposing the global cut A_(r-1) -> (all other names).

For r>=2 these are nontrivial proper modules or global cuts in a quotient with at least four remaining vertices. In the global-cut case the side containing all but the isolated source/sink is itself a proper nontrivial module. Hence U_r is critically indecomposable for r>=2: it is prime and every single-vertex deletion is decomposable. The r=1 quotient is C3 and is treated separately below.

## 3. Actual named response separator

Inflate every quotient name v by an arbitrary nonempty tournament D_v. Compare

    T=U_r[D_v] and T'=U_r^op[D_v],

keeping every physical child name and internal edge unchanged. Suppose the child at z is exactly a pair p->q. Let R=T-D_z. For disjoint unions X,Y of other complete child blocks, the marked-pair identity is

    [a_q a_X c_p c_Y](F_T-F_T') = C_YX(R)-C_XY(R).       (1)

Here C_XY counts ordered spanning covers of R by two nonempty directed paths, first ending in X and second starting in Y. The two actual ports have x->y, and a_v,c_v indicate old v->x, old v->y. Group variables a_X,c_Y are uniform only on the indicated child unions. They do not merge physical names or discard higher-degree terms.

If M=X union Y is a module of R and X->Y, the ordered-module injection gives

    C_YX-C_XY >= H(R/M)
        [q_2(X)H(Y)+H(X)q_2(Y)+q_2(X)q_2(Y)].           (2)

Thus |M|>=3 makes (1) strictly positive. The coefficient is exactly the alternating sum of sixteen actual two-port queries, so at least one query separates T and T'. This is an exact integer distinction, with all cubic and quartic contributions retained.

Apply the list in Section 2:

- If z=A_i, take X=D_(B_i), Y=D_(B_(i+1)). It suffices that their total order is at least three.
- If z=B_j with 1<=j<=r-1, take X=D_(A_j), Y=D_(A_(j-1)). It suffices that their total order is at least three. The order of these two sides is reversed from the numeric labels and is essential.
- If z=B_0, take Y=D_(A_0) and X to be all other exterior children.
- If z=B_r, take X=D_(A_(r-1)) and Y to be all other exterior children.

In the last two cases M=R has at least four physical names for r>=2 even with merely nonempty other children, so no extra size condition is required.

**Theorem.** For every r>=2, the preceding explicit pair-child and neighboring-size conditions imply distinct complete named two-port responses for T and T'. In particular, every inflation of U_r by nontrivial children, with at least one child exactly a pair, is distinguished from its quotient converse with all child interiors fixed. All other child tournaments may have arbitrary size, orientation and ordered-cover profile.

## 4. Smallest-order and scope boundaries

At r=1, U_1 is a directed triangle. Every deletion leaves a two-vertex ordered quotient. Formula (2) is positive if the two exterior children have total order at least three. If both are singletons, the marked-run switch is a bijection and this particular quartic is zero. This is not an exact full-response collision.

For a non-endpoint deletion at larger r, two singleton exposed children similarly give a zero quartic, not a collision. Nontriviality of all other children is a sufficient convenience; the sharper local conditions in Section 3 are the actual hypotheses.

No larger child is silently replaced by a pair. The theorem compares the two displayed quotient orientations with interiors fixed, not all arbitrary original-module XOR families. Arbitrary strong outer ambient lifting is not claimed by this exact-integer argument. The unbounded conclusions follow from the elementary proof, not finite checks.

Dependencies: ../marked_pair_quotient/PROOF.md (exact grouped quartic) and ../ordered_module_skew/PROOF.md (support-preserving cover injection and bound). The separate critical-family synthesis, if used, cites its classification source independently; this proof does not depend on that classification.

## 5. Verification

verify.py checks the staircase/cyclic definitions, strong connectivity, primeness and every displayed deletion module through r=7. It retains all sixteen actual query values on both sides for every root at r=1 and r=2 with pair children, larger quotient controls with minimally sufficient sizes, arbitrary nontransitive child controls, and zero-quartic two-singleton boundaries. Both normal and optimized runs must agree. These are implementation controls for the proof, not a bounded replacement for it.
