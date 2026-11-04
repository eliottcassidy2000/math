# Twin reversal retains a quartic contact response

Status: **PROVED exact all-order contact formula; FINITE-EXACT sharp minimum order five for quartic response.** This is not a full-response collision or an obstruction to the existing twin-separation theorem.

## Scope of the inherited cubic theorem

The current-main `two_port_twelve_principles_20261103/PROOF.md`, T4, proves that `F_T-F_(T^op)` has degree at most three when **every core arc** is reversed. Its proof reverses the entire Hamilton word and exchanges the two ports. It does not apply to reversal of one module while the surrounding tournament remains fixed. The following exact formula and strong example make that boundary explicit.

Let a bit `a_i=1` mean `i->x`, a bit `c_i=1` mean `i->y`, with `x->y`. Vertices and queries are named and fixed. Let `{u,v}` be a homogeneous pair with `u->v` in T. Let T' differ only in the internal edge, now `v->u`. Write Q for the contracted tournament with root z. All definitions below retain the individual attachments of the ports to the other vertices.

## Theorem 1: exact contact formula

Consider every permutation of `V(Q) union {x,y}`. Test all consecutive arcs except those between z and a port against the actual tournament Q and the prescribed exterior attachments. Discard any permutation failing a tested arc. Count the remaining permutations in four disjoint classes:

- A_x: x is adjacent to z in the word, but y is not;
- A_y: y is adjacent to z in the word, but x is not;
- P_xy: the word contains the consecutive block x,z,y;
- P_yx: the word contains the consecutive block y,z,x.

Words with neither port adjacent to z are unused. These four nonnegative integer counts depend on Q and on `a_i,c_i` for exterior old vertices, but not on the four twin attachment bits.

Then exactly

    F_T-F_T' = (a_v-a_u) A_x + (c_v-c_u) A_y
             + ((1-a_u)c_v-(1-a_v)c_u) P_xy
             + ((1-c_u)a_v-(1-c_v)a_u) P_yx.             (1)

**Proof.** Hamilton words not using the internal edge have identical contributions in T and T', so cancel. The remaining words contain either the block u,v or the block v,u. Contract this block to z. Every adjacency to an old exterior vertex is unchanged because the pair is homogeneous. With no adjacent port the two contributions cancel. With only x adjacent, incoming weights differ by `(1-a_u)-(1-a_v)=a_v-a_u`, while outgoing weights differ by `a_v-a_u` as well. This gives the first term regardless of which side contains x. The argument for only y gives the second term. When both ports are adjacent, the two possible orders are precisely the last two classes; their boundary products give the displayed coefficients. Expansion and contraction are inverse, and every class includes possible endpoint positions. QED.

The response has degree at most four, and bidegree at most (2,2), by counting the at most two incidences at each port in any Hamilton word. Formula (1) is an exact **full named response**, not a replacement by the scalar Hamilton count, endpoint marginals, or a port-aggregated run profile. In particular A_x can retain one other incidence at x and two incidences at y. Multiplication by `a_v-a_u` can therefore genuinely contribute total degree four.

## A named quartic coefficient is an exact deletion-cover skew

For distinct old exterior vertices i,j, define C_ab(R) to count ordered spanning covers of R by two nonempty directed paths, the first ending at a and the second starting at b. These are defined whether or not the arc a->b exists. Then

    [a_i a_j c_i c_v](F_T-F_T')
        = C_zj(Q-i)-C_jz(Q-i).                            (2)

**Proof.** Only `(c_v-c_u)A_y` in (1) can contribute this monomial. The A_x term contains a twin a-bit; the P terms either contain such a bit or have degree at most two in exterior bits and hence cannot supply `a_i a_j c_i`. In A_y, the c_i factor requires the other neighbour of y to be i. The a_i factor requires x adjacent to i, and the second a-bit a_j requires its other neighbour to be j. There are precisely two possible consecutive blocks:

    z,y,i,x,j     and     j,x,i,y,z.

Their tested attachment products are respectively

    (1-c_i) a_i (1-a_j),     a_j (1-a_i) c_i,

with cubic coefficients +1 and -1. Delete y,i,x and leave the z,j boundary untested. The first block becomes exactly an ordered cover counted by C_zj(Q-i); the second becomes one counted by C_jz(Q-i). This is reversible, including cases when either cover path is a singleton. QED.

Thus a necessary condition for an internal twin-reversal difference to have degree at most three is that every deletion Q-i has zero cover skew between its root and every other retained vertex. This is an exact necessary condition, not a claimed sufficient criterion. It identifies information genuinely omitted by truncating the named response at degree three.

## Theorem 2: degree four occurs first at five old vertices, even strongly

Take

    T  = (24,17,19,6,8),
    T' = ( 8,17,19,6,9).

Rows are outgoing-neighbor bit masks on names 0 through 4. The pair `{0,4}` is homogeneous. In T its edge is `0->4`; in T' its edge is `4->0`. Both have the strong quotient

    Q=(8,1,3,6), root 0.

For example Q contains the spanning cycle `0->3->2->1->0`. Replacing one quotient vertex by a nonempty pair preserves strongness, so both T and T' are strong.

The coefficient in their unique Boolean multilinear response difference is

    [a_1 a_3 c_1 c_4](F_T-F_T') = 2.                    (3)

This is not merely an unreduced monomial: it is a Boolean finite difference on four **distinct** variables. The exact 16-query certificate is retained in `normal.json`. A separate complete permutation expansion gives the same coefficient. All 1,024 attachment assignments for each of T and T' are independently evaluated by integer Hamilton subset DP. Boolean inversion of their difference agrees with every coefficient of the permutation expansion, not only (3). The difference has 4 nonzero linear, 20 quadratic, 56 cubic and 32 quartic coefficients.

The finite minimum-order claim is exhaustive: for quotient order q=1,2,3,4, all `2^(q choose 2)` labelled tournaments, with fixed named root 0, are expanded by a second twin. Renaming any chosen root to 0 shows this covers every rooted quotient up to relabelling. The degree histograms are

    q=1: one response of degree 1;
    q=2: two responses of degree 3;
    q=3: eight responses of degree 3;
    q=4: all 64 responses have degree 4.

Hence no twin reversal on at most four old vertices has a quartic difference, while a strong example exists on five. The q=4 statement is a complete finite census, not an asserted all-order universal degree-four law.

## What this does and does not imply

This refutes transporting T4's **degree-three reconstruction statement** from whole-core converse to arbitrary internal homogeneous-pair reversal. It does **not** refute named twin separation: `two_port_eight_core_20261103/algebra/PROOF.md`, C2, already separates opposite twins by the positive ordered-cover contrast. In particular, this packet does not claim that a low-degree query battery fails to distinguish these twins, only that it cannot recover their entire polynomial by assuming its quartic coefficients vanish.

The same distinction applies to uniform-run laws such as the S185 gate

    H(Q_j[A])=(27+2j)H(A)+21q_2(A)+4q_3(A).

That scalar substitution response, even for every A, is not the named two-port response. Formula (1) explicitly retains the port contacts which a scalar observer forgets. There is no collision claim, no hidden module-converse candidate, and no universal reconstruction theorem here.

## Reproduction

`python verify.py` and `python -O verify.py` produce byte-identical JSON. The checker uses explicit exception gates. It enumerates all 75 labelled fixed-root quotients through order four, computes their complete response polynomials by literal Hamilton words, and independently evaluates all 2,048 hostile attachment queries by endpoint subset DP. It verifies the deletion-skew identity (2) for all 400 ordered distinct exterior pairs in the quotient census. It additionally verifies (1) on all 1,024 hostile assignments by directly enumerating its four contracted-word classes. This is independent arithmetic/algorithmic cross-checking within the author script; an independent author audit is separate.
