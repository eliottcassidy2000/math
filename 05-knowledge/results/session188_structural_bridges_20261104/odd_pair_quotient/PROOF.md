# Every odd quotient with two-vertex children is visible

**PROVED.** Let Q be any tournament on an odd number q>=3 of fixed quotient labels. Replace each label i by a fixed named transitive two-vertex block B_i. Let T=Q[B_i] and U=Q^op[B_i], retaining every child's interior and every physical name. Then their complete named two-port responses differ already on block-uniform queries. No assumption that Q is prime is required.

Equivalently, simultaneously reversing the whole core and all q pair modules cannot cancel. This is a fully covered laminar family. The theorem is not about arbitrary child sizes, even-order quotients, or arbitrary legal inversion sequences.

## 1. The modulo-eight trace

Attach x->y. A coarse bit a_i means that both names of B_i beat x, and c_i means that both beat y. Write

    Delta(a,c)=F_T(a,c)-F_U(a,c)

with coarse masks expanded onto the fixed physical pair blocks. All coefficient extractions below are in the unique Boolean multilinear polynomial of these actual query counts. Define

    D_ij=[a_i c_j] Delta,
    h_ij=[a_i a_j c_j] Delta,                  i!=j,
    L_j=sum_(i!=j) D_ij + sum_(j->i in Q) h_ij.

We prove the stronger identity

    L_j is divisible by 4 for every j,
    sum_j L_j = 4q modulo 8.                         (1)

Since q is odd, at least one actual coarse query has nonzero contrast. We do not assert that every nonzero individual contrast is four modulo eight.

## 2. Classical odd-cycle expansion

For every tournament R,

    H(R)=sum_C 2^(number of nontrivial cycles in C),       (2)

where C ranges over collections of vertex-disjoint directed odd cycles of length at least three; unused vertices are singleton cycles. Cycles are counted up to rotation, with their direction retained. This is the Hamilton specialization of the tournament odd-cycle power-sum formula in Grinberg and Stanley, [The Redei--Berge symmetric function of a directed graph](https://arxiv.org/pdf/2307.05569), arXiv:2307.05569v1, Theorem 1.39. In the symmetric-function convention the one-color specialization counts complementary Hamilton paths; for tournaments complement is whole converse and reversing a path preserves the scalar count.

Modulo eight, only collections of zero, one or two nontrivial odd cycles matter. We retain the part of their attachment polynomial involving at least one a-variable and at least one c-variable, called the mixed part. This includes every D_ij and h_ij in (1), and every other mixed coefficient; no truncation of the actual response is assumed.

In a cycle containing a port, erase the two physical names of each pair block, retaining its quotient label at each occurrence. If just one name of a block occurs, there are two choices for it. If both names occur but are not consecutive in the same directed cycle, their two assignments are both legal because every intervening incidence is exterior and uniform. If both occur consecutively in the same directed cycle, their internal direction forces one assignment. Call the first two types free blocks. A cycle containing a named port has a distinguished cyclic starting point, so these assignments introduce no rotational division.

## 3. A cycle containing both ports

First consider one nontrivial odd cycle containing both x and y. It contains an odd number of old names. Consequently at least one pair block appears exactly once, so there is at least one free block. Each free block contributes an independent factor two to the number of assignments. The coefficient two from (2) means that modulo eight only skeletons with exactly one free block survive. That block occurs once; every other used block occurs twice consecutively in its forced internal order.

Collapse each consecutive pair to one quotient label. The surviving skeleton is a simple directed cycle in the quotient extension Q+x+y containing both ports. Any one of its old labels can be the singleton block, giving two assignments, while all other labels are forced pairs. Thus this skeleton contributes four times its number of old labels modulo eight. That number is odd exactly when the quotient cycle itself has odd length. Therefore the total contribution is

    4 times the number of directed odd quotient cycles containing x,y,

as an attachment polynomial modulo eight.

If there are two nontrivial odd cycles and both ports lie on the same one, the coefficient in (2) is four. The port-containing cycle still has an odd number of old names, hence a block occurs exactly once on it. If its other name is unused, exchange the used name; if its other name lies on the other cycle, exchange those two physical names. Both assignments are legal and distinct, because the block's two occurrences are not consecutive in one cycle and the cycle containing the ports is distinguished. This is a fixed-point-free involution, so these collections contribute zero modulo eight. Collections not using both ports have no mixed part.

## 4. Two cycles with one port each

The remaining mixed contribution has two disjoint odd cycles, one through x and one through y. Its coefficient is four. If a free block occurs, swapping its two physical labels gives a free involution on the port-distinguished cycle pair. Thus modulo eight only collections in which each used block occurs as one consecutive whole pair survive.

Collapse every pair. We obtain two vertex-disjoint simple directed cycles in the quotient extension, one through each port. Their quotient lengths can be even or odd: the expanded lengths are automatically odd. Conversely, every such quotient-cycle pair uniquely expands into a surviving pair collection. Let B_Q(a,c) count these pairs over the integers.

A collapsed one-port cycle cannot have just one old quotient label: uniform attachment forbids x->B_i->x, and likewise for y. Removing the ports therefore gives disjoint directed old paths P and R, each having at least two old vertices. Let P start at s and end at t, and R start at u and end at v. Their contribution is

    a_t(1-a_s)c_v(1-c_u).

The four endpoint labels s,t,u,v are distinct. Reversing both old paths is a bijection to the opposite quotient. Hence, modulo two,

    B_Q+B_(Q^op)
      =sum_(P,R) [a_t c_v+a_s c_u
                  +a_s a_t(c_v+c_u)+c_u c_v(a_t+a_s)].    (3)

The quartic terms cancel exactly in this characteristic-two contrast. This is a proved cancellation in a specified cycle-pair term, not an assumption about arbitrary module-reversal quartics.

Every summand in (3) has two mixed quadratic terms, so the sum of all its mixed quadratic coefficients is zero modulo two. Its cubics have disjoint a-endpoint and c-endpoint labels, so it has no coefficient of the shared-index form a_i a_j c_j. Therefore the trace functional used in (1) annihilates (3) modulo two, for every quotient orientation Q.

## 5. Reduction to the incoming E1 identity

Let

    delta_Q=(F_Q-F_(Q^op))/2 modulo two

be the whole-converse contrast on the quotient's actual singleton names. By (2) modulo four, its mixed part is the difference of the counts of single directed odd cycles containing both ports. Sections 3–4 therefore prove the coefficientwise mixed identity

    Delta = 4[delta_Q+B_Q+B_(Q^op)] modulo eight.          (4)

In particular all mixed coefficients of Delta, including every D_ij and h_ij, are divisible by four.

The incoming [E1 theorem](https://github.com/eliottcassidy2000/math/blob/76b3b1a77d124100f03026420f6b821533cdec15/05-knowledge/results/two_port_parity_separation_20261104/PROOF.md) states, for each quotient label j when q>=3,

    sum_(i!=j)[a_i c_j]delta_Q
      +sum_(j->i)[a_i a_j c_j]delta_Q = 1 modulo two.

Summing over j, applying (4), and using the trace cancellation from (3) gives sum_j L_j=4q modulo eight. This proves (1) and the theorem.

The odd quotient assumption is used only at the final nonzero trace. At even q the identity gives zero and supplies no separator. No conclusion of equality at even q follows. The q=1 case is excluded: its quotient converse is the identity, and incoming E1 does not apply.

## 6. Explicit query battery and reversal-cube corollary

Every coefficient in the trace is an ordinary finite difference of actual coarse queries. A fixed battery sufficient to compute the trace consists of:

- for each c=e_j: a=0, every singleton a=e_i, and a=e_i+e_j whenever j->i;
- for c=0: a=0, every singleton, and every two-element mask.

There are at most 2q^2+q+1 queries. Expanding each coarse bit onto its physical pair makes each a genuine named attachment. Since the trace is nonzero, this finite battery separates T and U. It does not reconstruct arbitrary response functions.

Independently conversing any of the q pair interiors leaves every coarse query unchanged, exactly by the module run-reversal bijection. The battery therefore identifies whether the quotient is Q or Q^op even when all child bits are unknown. Once that orientation is known, one incoming lifted module separator per pair reads its interior bit and is individually insensitive to every other pair reversal. Thus all 2^(q+1) quotient/pair-converse states have distinct complete named responses. The displayed construction uses at most 2q^2+2q+1 adaptive queries or 2q^2+3q+1 fixed queries. These are nonoptimal upper bounds, not lower bounds.

This corollary is a two-level fully covered family with two-vertex children. It does not settle arbitrary child tournaments or arbitrary legal inversion sequences.
