# Two small transitive blocks and one arbitrary tournament

## 1. An unrestricted third-block family

Let C be any nonempty named tournament, and let A and B be disjoint fixed
transitive tournaments. Put T=C3[A,B,C] with quotient A->B->C->A, and
U=C3^op[A,B,C], preserving the three interiors. Write

Delta_IJ=F_T(1_I,1_J)-F_U(1_I,1_J),

where x->y is fixed and each new-port incidence uses the same physical
old names. Let q_k=q_k(C) be the number of ordered spanning k-path covers,
with q_0=0, q_1=H(C), and q_k=0 past the order of C.

**Theorem 1.** If A=B=TT2, then

Delta_AB = -8q_1-14q_2+4q_3.                              (1)

For |C|>=2 this gap is 4 modulo 8, and for |C|=1 it is -8. Thus this single
fixed physical query separates the quotient converse for every C, without
a log-concavity assumption on its profile.

**Theorem 2.** If A=B=TT3, then

Delta_BC+Delta_CA = 242q_1-283q_2-552q_3.                  (2)

This sum is strictly negative for |C|>=2 and equals 242 for |C|=1. Thus
these two fixed physical queries separate the quotient converse for every
C, again without a log-concavity assumption.

These are all-order families with an arbitrary third-block interior, not
just unequal transitive blowups or finite controls. The equalities A=B=TTm
mean the two blocks have that isomorphism type and size; their physical
names remain disjoint and are never identified.

## 2. Exact coefficient identities

The run-series identity from ../unequal_covered_response/PROOF.md is valid
for arbitrary unequal profiles. Its old-walk matrix is

G = 1/(1-abc) * [[a,ab,abc],[abc,b,bc],[ac,abc,c]],

and the opposite quotient uses G^t. With e the all-ones vector and alpha,
gamma the port indicators, its two-port enumerator is

P(G;alpha,gamma) =
 (1+e^t G alpha)(1+(e-alpha)^t G gamma)(1+(e-gamma)^t G e)
 +(1+e^t G gamma)((e-gamma)^t G alpha)(1+(e-alpha)^t G e).

Apply the linear rule a^i b^j c^k -> q_i(A)q_j(B)q_k(C).
For TT2 the first two block profiles are (0,1,2); for TT3 they are
(0,1,6,6). Extracting their finite a,b supports gives the following exact
linear forms in C's profile:

For two TT2 blocks:

Delta_AB = -8q_1-14q_2+4q_3,
Delta_BC = 56q_1-46q_2-98q_3-24q_4,
Delta_CA = -58q_1+11q_2+94q_3+24q_4.

For two TT3 blocks:

Delta_AB = -44q_1-270q_2-120q_3+108q_4,
Delta_BC = 1192q_1+1494q_2-2090q_3-2484q_4-432q_5,
Delta_CA = -950q_1-1777q_2+1538q_3+2484q_4+432q_5.

The verifier expands the integer numerator polynomials, convolves the
denominator (1-abc)^-3 using binomial coefficients, and independently checks
these coefficient vectors. There is no omitted tail: if A and B have
order m, deleting the two ports leaves at most three quotient walk segments.
Within each cyclic segment the number of C-runs exceeds the number of
A-runs by at most one, so C has at most m+3 runs. The verifier extracts
through m+3, where the last coefficient is zero, and checks the numerator
support directly as an additional finite-tail certificate.

Adding the last two TT3 formulas proves identity (2). The TT2 identity
(1) is the first formula directly.

## 3. Parity and a universal cover lower bound

For every tournament C, q_1 is odd and positive by Redei's theorem.
For |C|>=2, unordered bipartitions of its vertices give

q_2/2 = sum_{unordered {S,C-S}, S nonempty proper}
                  H(C[S])H(C[C-S]).

Each summand is odd. There are 2^(|C|-1)-1 unordered bipartitions, an odd
number. Therefore q_2=2 modulo 4. Also q_3 is divisible by 3!=6, since
ordering the three paths acts freely on every unordered three-path cover.
Consequently the three terms of (1) are respectively 0,4,0 modulo 8.
This proves Theorem 1 for every nontrivial C. For a singleton q=(0,1),
(1) equals -8.

For |C|>=3, cut any Hamilton path at any of its |C|-1 internal positions.
The resulting ordered two-path cover determines the original concatenated
path and its cut. This is an injection, so

q_2(C) >= (|C|-1)H(C) >= 2H(C).

For |C|=2 the same final inequality holds directly, since H=1 and q_2=2.
Thus (2)<=-324q_1-552q_3<0 for every nontrivial C. The singleton value
242 is immediate. This proves Theorem 2.

The use of Redei parity in Theorem 1 is classical. No universal cover
log-concavity or real-rootedness statement is used or inferred.

## 4. Full two-level reversal cubes

Each block-uniform query is unchanged by any child converse because every
child remains a module in the actual extension. Thus Theorem 1 determines
the quotient orientation independently of the three child-converse bits
using one actual query. Theorem 2 determines it using the sum of two actual
queries; the sum is not counted as one query.

After that, one lifted isolation query per nontrivial child determines its
converse bit, independently of the other children, by the existing disjoint-
module isolation theorem. This is the same argument as
../unequal_reversal_cube/PROOF.md, with a smaller root battery.

If |C|>=2, all three children are nontrivial, and the sixteen named full
responses are therefore distinct:

* C3[TT2,TT2,C]: four adaptive queries or seven fixed queries suffice.
* C3[TT3,TT3,C]: five adaptive queries or eight fixed queries suffice.

The fixed battery uses the root queries plus three child decoders for each
of the two quotient orientations. A singleton C has no independent child
bit: there are only eight distinct states, with corresponding bounds three
adaptive/five fixed and four adaptive/six fixed.

No optimality or arbitrary strongly connected outer-quotient transfer is
claimed. This theorem concerns the stated cyclic substitutions themselves.

## 5. An exact single-probe cancellation

For T=C3[TT3,TT3,C3] and its quotient converse U, both H counts are 2823,
and the cyclic uniform gaps are

(Delta_AB,Delta_BC,Delta_CA)=(-2472,0,-4284).

With blocks on names {0,1,2}, {3,4,5}, {6,7,8}, the physical query
(a,c)=(56,448) gives 65499 on both sides. The other two cyclic queries give
65499 versus 67971 and 62583 versus 66867. This is a cancellation of one
coarse probe, not a collision of three probes or of the full named response.

## 6. Reproduction and scope

The standard-library checker verifies the exact run-kernel coefficients,
all specified literal queries, the parity/bound controls, and the declared
orientation-cube batteries. Normal and optimized runs must agree exactly.

    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json

The proofs above apply to every C. Finite tests are implementation controls,
not their justification. General cyclic three-block probe collisions and
the unrestricted full-response collision problem remain open here.
