# Every named transitive-module order survives every scalar context

**PROVED:** In any fixed rooted tournament quotient, two distinct named linear orders of the substituted transitive module have different full named two-port responses. This holds for every module size and every quotient, including strong quotients. The separating queries recover the entire module order. The two compared tournaments need not share any homogeneous pair, so this is not merely opposite-twin separation. **No full-response collision or general reconstruction theorem is claimed.** Historical novelty is not asserted.

## 1. Observer and class

Let Q be any tournament with distinguished root z, and let M be a fixed set of r>=1 physical names. For a linear order pi on M, let T_pi=Q[TT_pi], with the same named exterior vertices and identical root incidences for all pi. Write F_T(a,c) for the Hamilton-path count after adjoining named ports x->y with v->x iff a_v=1 and v->y iff c_v=1. Arbitrary comparison retains every Boolean coefficient of bidegree at most (2,2), including total degree four. Only a comparison with the **whole converse** is governed by the inherited full cubic theorem; nothing here replaces that cubic by selected moments.

For each physical v in M, take the following literal query:

- x beats every member of M;
- y loses only to v in M, and beats the other members;
- both x and y have exactly z's incidences to every exterior vertex;
- x->y, as always.

These prescriptions are independent of pi. They make M union {x,y} a common module. Let i be the zero-based rank of v in pi. The query value depends only on i, and is strictly increasing with i. Thus the r fixed named queries recover pi by sorting their values. In particular, any two distinct orders are separated by at least one identical physical query.

## 2. Exact color polynomial

For a tournament A, let q_k(A) count ordered spanning covers by k nonempty directed paths. Define

    P_A(t) = sum_k binomial(t,k) q_k(A).

For integer t>=0, this counts assignments to t named colors, together with a Hamilton path on each nonempty color class; an empty class contributes 1. An ordered join A->B therefore satisfies P_(A->B)=P_A P_B. This is the same color polynomial used in S186's scalar-context hostile.

Let B_(r,i) be the extended module in the query above. Since x dominates M and y, it is a source singleton, so P_B(t)=t P_(M+y)(t). We claim

    P_(B_(r,i))(t)
      = t^(r+2) + t^(r-i+1) ((t+1)^i-(t-1)^i).             (1)

To prove it, first fix the color of y. If v has a different color, every color class is transitive; these assignments contribute (t-1)t^r. If v has y's color, let k of the i predecessors of v also have that color. That color class has exactly 2^k + 1_[k>0] Hamilton paths. Indeed a path containing v->y chooses any subset of those k predecessors to appear before v, and all remaining vertices then follow y in their original order. These give 2^k paths. A path not using v->y must start at y and exists exactly when k>0. Vertices later than v impose no further choices. Summing colors of predecessors gives

    (t+1)^i + t^i - (t-1)^i.

The other r-1-i vertices and the color of y supply t^(r-i). Adding the different-color contribution simplifies to t^(r+1)+t^(r-i)((t+1)^i-(t-1)^i), and multiplying by t proves (1). The counting argument for positive integer t proves a polynomial identity, also valid at t=0. Exponents are nonnegative; usual polynomial convention includes 0^0=1.

For 0<=i<r-1, the consecutive difference is

    P_(B_(r,i+1))(t) - P_(B_(r,i))(t)
      = t^(r-i) ((t+1)^i+(t-1)^i)
      = 2 sum_{0<=h<=i, h even} binomial(i,h) t^(r-h).     (2)

All monomial coefficients are nonnegative. The identity

    t^m = sum_k k! S(m,k) binomial(t,k)

has nonnegative coefficients. Therefore every q_k(B_(r,i+1))-q_k(B_(r,i)) is nonnegative. For k=1 the difference is 2 when i=0 and 2^i when i>=1, so it is strictly positive. This establishes coefficientwise monotonicity of the **complete** ordered-path-cover profile, not only H.

## 3. Every ambient quotient preserves strictness

The occupied-gap substitution identity, inherited from S176 and used in S185, gives

    H(Q[A]) = sum_k c_k(Q,z) q_k(A),                       (3)

where each c_k is a nonnegative integer counting exterior skeletons with k ordered nonempty module runs. In particular c_1=H(Q)>0: contracting the sole run is precisely a Hamilton path of Q. Positivity follows from the usual insertion proof that every tournament has a Hamilton path, including the one-vertex quotient.

Combining (2) and (3), consecutive rank queries differ by a positive integer, at least 2H(Q) when i=0 and at least 2^i H(Q) when i>=1. Thus strictness holds for **every** quotient and unbounded r, with no strongness restriction and no upper bound on the number of exterior runs. In strong Q the resulting examples are strong as well. No inference is made from scalar equality alone: the scalar quantities in (3) are the Hamilton counts of the actual named port extensions, so each is exactly one full-response table entry.

This covers much more than reversing an already-common adjacent twin pair. For example pi=(0,1,2,3) and rho=(1,3,0,2) have no common consecutive unordered pair. Substituted into the directed triangle at one root, the resulting six-vertex strong tournaments have no common homogeneous pair at all, but this theorem still separates them. Different fixed roots, changed exterior arcs, nontransitive modules, or comparisons to arbitrary unrelated tournaments are not covered.

## 4. Exact S185 three-state gate integration

Use the inherited six-vertex seeds

    Q0=(16,5,33,39,46,3),
    Q1=(16,5,49,39,42,3),
    Q2=(18,4,49,39,42,3),

with root 5. For any A, S185 proves H(Qj[A])=(27+2j)q1(A)+21q2(A)+4q3(A). Since q1=P(1), q2=P(2)-2P(1), and q3=P(3)-3P(2)+3P(1), the exact rank query is

    F_(T_j,pi)(query_v)
      = (2j-3) P_(B_(r,i))(1)
        + 9 P_(B_(r,i))(2) + 4 P_(B_(r,i))(3),             (4)

where P is the explicit formula (1). Although the first coefficient in (4) can be negative, monotonicity follows from the original positive q_k expansion, not from incorrectly assuming positivity of this evaluation expression. Formula (4) recovers arbitrary physical clone order for all r and all three states, rather than only comparing the order to its reversal.

## 5. Exact checks and evidence boundary

`verify.py` is standalone standard-library Python, with exception gates that remain active under `python -O`.

- It computes complete q_k profiles by subset Hamilton DP followed by ordered subset-cover convolution for every r=1,...,6 and every rank i (21 extended modules). The independently computed profiles verify (1) and coefficientwise monotonicity.
- It evaluates 750 literal port-extension Hamilton counts across every labelled quotient on 1,...,4 vertices with root fixed at 0 (covering all rooted cases up to relabelling), every r=1,...,4 and each rank. It directly verifies strictness without using (3).
- Its six-vertex strong example has module orders (0,1,2,3) and (1,3,0,2), and no common homogeneous pair. Literal Hamilton-word expansion computes the **complete** named Boolean polynomials for T, U and T^op. Every one of the 484 bidegree-(2,2) interpolation queries for each polynomial is independently checked by subset Hamilton DP (1452 checks). The arbitrary difference has degree exactly four; the whole-converse difference is retained in full and has degree at most three. Both complete differences are serialized, with omitted coordinates explicitly zero within their full respective bases.
- It independently checks all 45 actual S185 gate queries for j=0,1,2, r=1,...,5 against (4).

Reproduce with `python verify.py` and `python -O verify.py`. Their JSON receipts agree byte for byte. Finite controls do not supply the unbounded claim; the color and occupied-gap bijections above do. Universal full named reconstruction, prime hidden-converse cores, arbitrary changes of module type, and existence of a full-response collision remain unresolved.
