# Unequal cyclic modules: a three-query separation criterion

## 1. Scope and statement

Let A, B, C be nonempty tournaments on disjoint fixed vertex sets. Put

T = C3[A,B,C], with A -> B -> C -> A,

U = C3^op[A,B,C], preserving all three interiors.

This is a fully covered laminar reversal: reverse the whole vertex set and
each of A, B, C. Let q_r(D) count ordered covers of D by r nonempty directed
paths. Set q_0(D)=0 and q_r(D)=0 for r>|D|. For two new named vertices x,y,
always x->y, let F_D(a,c) be the Hamilton-path count with v->x iff a_v=1
and v->y iff c_v=1. For I,J in {A,B,C}, define

Delta_IJ = F_T(1_I,1_J) - F_U(1_I,1_J).

**Theorem.** Suppose that for every k>=1 and distinct blocks I,J,

q_(k+2)(I) q_k(J) <= q_(k+1)(I) q_(k+1)(J).                 (1)

Then Delta_AB+Delta_BC+Delta_CA<0. In particular F_T != F_U.

There is no equal-order or equal-profile requirement. Condition (1) is a
finite, exact cross-profile condition, not an asserted property of all
tournaments. Equal log-concave profiles satisfy it, recovering the separate
equal-profile sufficient theorem in ../covered_response. This packet gives
the unequal-profile extension and an explicit unequal seven-vertex control.

Every three blocks of orders 2 or 3 satisfy (1), so the theorem applies to
all these choices, whether a three-vertex block is transitive or cyclic.
Indeed their profiles are (1,2), (1,6,6), or (3,6,6). Only k=1 can be
nontrivial: q_3(I)q_1(J)<=q_2(I)q_2(J) follows immediately from these
three possibilities. For k>=2 the left side is zero.

## 2. The exact unequal-profile identity

Use formal variables a,b,c marking the number of maximal old-vertex runs
in the three modules, not their numbers of vertices. Put d=1-abc and

G = 1/d * [[a, ab, abc], [abc, b, bc], [ac, abc, c]].

G_ij is the generating function of nonempty directed quotient walks that
start in block i and end in block j. Consecutive labels are different.
Each extra traversal of the directed triangle multiplies the weight by abc.
For the opposite quotient use G^t.

Define the linear functional L on monomials by

L(a^r b^s c^t)=q_r(A)q_s(B)q_t(C).

The definition extends to the rational series below because all but finitely
many monomials are killed by the finite supports of the q profiles. Given a
quotient walk with r appearances of A, its appearances are filled, in their
walk order, by an ordered r-path cover of A. This establishes the run-series
enumeration without any assumption on the interiors.

Let e be the all-ones column and let alpha,gamma be zero-one port incidence
columns on the three blocks. For a matrix G of old walks define

P(G;alpha,gamma) =
 (1+e^t G alpha)(1+(e-alpha)^t G gamma)(1+(e-gamma)^t G e)
 +(1+e^t G gamma)((e-gamma)^t G alpha)(1+(e-alpha)^t G e).

The two summands enumerate x before y and y before x, respectively. The
middle segment may be empty only in the first order, since x->y. Empty
end segments are allowed. Thus the literal two-port count is L(P).

The following polynomial identity is exact:

sum_cyclic [P(G;e_A,e_B)-P(G^t;e_A,e_B)] = S/d^2,           (2)

where

S = a^2 b^2 + b^2 c^2 + c^2 a^2
    -abc(a+b+c)
    +(a^2 b+a^2 c+b^2 a+b^2 c+c^2 a+c^2 b)-6abc.

Here the cyclic sum means the three ordered port pairs (A,B),(B,C),(C,A).
One can verify (2) by multiplying by d^3 and expanding integer polynomials;
the supplied verifier does exactly this with dictionaries of monomials and
integer coefficients, without a symbolic-algebra dependency.

Write x_i=q_(k+2)(I), y_i=q_(k+1)(I), z_i=q_k(I), with i=1,2,3 denoting
A,B,C. Since d^-2=sum_(k>=0)(k+1)(abc)^k, (2) gives

Delta_AB+Delta_BC+Delta_CA = sum_(k>=0)(k+1) E_k,           (3)

E_k = sum_cyclic x_i x_j z_l - sum_cyclic x_i y_j y_l
      +sum_(i,j,l distinct) x_i y_j z_l -6 y_1 y_2 y_3.    (4)

The three cyclic indices are (1,2,3),(2,3,1),(3,1,2), and the second
positive sum has all six permutations. Formulas (3)-(4) are unconditional.

## 3. Proof of the sign criterion

For k>=1, (1) says x_i z_j<=y_i y_j for i!=j. Therefore

x_i x_j z_l <= x_i y_j y_l,

so the three cyclic terms in the first positive sum are dominated term by
term by the first negative sum. Also each of the six terms x_i y_j z_l is
at most y_i y_j y_l. Thus E_k<=0.

For k=0 all z_i vanish and all y_i=H(I)>0, by the existence of a Hamilton
path in every tournament. Hence

E_0 = -sum_cyclic q_2(I)H(J)H(L)-6H(A)H(B)H(C)<0.

All remaining weights k+1 are positive, proving (3) is strictly negative.
The proof uses cross-multiplied inequalities and remains valid when some
q coordinates vanish. It does not divide by any potentially zero entry.

## 4. Unequal fully covered example, complete named data

Take A=TT2, B=TT2, C=C3 on names {0,1}, {2,3}, {4,5,6}. Then H(T)=H(U)=159,
but

F_T(3,12)=1959,   F_U(3,12)=2043.

Thus Delta_AB=-84. The complete three-query sum is -324, as predicted by
(3). The q profiles (1,2), (1,2), (3,6,6) are unequal. This example is not
a relabeling of the equal-block six-vertex example.

The receipts retain all coefficients of F_T-F_U, including the complete
cubic and quartic data, with omitted coefficients explicitly zero. The
polynomial is formed from all Hamilton words of the two-port extension,
and independently checked against literal subset dynamic programming at
every input of Hamming weight at most four. These inputs determine every
coefficient of a polynomial of degree at most four. The degree bound holds
because the two new vertices have at most four incident edges in a path.
All 64 block-uniform queries are also retained and checked independently.

## 5. Boundaries and a useful failed sign extension

Separate log-concavity of the three profiles is not enough to infer the
negative sign in the theorem. Take A=TT2, B=TT3, C=TT5. Their profiles are

(1,2), (1,6,6), (1,30,150,240,120),

all log-concave, but the exact cyclic sum in (3) is +4540. The three actual
two-port counts are checked directly, so this is an integer sign reversal,
not a rounding effect. It is not a collision: a positive sum also separates
the response. It shows why replacing (1) by three independent log-concavity
assumptions would invalidate this particular negative-sign proof.

The general unequal-profile cyclic sum might vanish; no universal nonzero
claim is made. Nor is universal q-profile log-concavity proved here. The
primary Grinberg--Stanley paper, arXiv:2307.05569v1, Theorem 1.39 in the
linked PDF, gives an
odd-cycle power-sum expansion of the Redei--Berge symmetric function; that
expansion alone is not a proof of the cross inequalities (1). We do not
attribute such inequalities to that paper.

This is a separation theorem for the specified cyclic three-block family.
It does not prove injectivity of the full named response for arbitrary
legal module-reversal sequences or arbitrary tournaments.

Primary reference consulted: Darij Grinberg and Richard P. Stanley,
The Redei--Berge symmetric function of a directed graph (2023),
https://arxiv.org/abs/2307.05569 .
