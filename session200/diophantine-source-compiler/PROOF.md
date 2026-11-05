# Physical source words as integral, field-invariant Diophantine witnesses

Session 200. **PROVED (independently audited):** full strong C1 and C2
membership admit explicit quadratic systems with O(log n) variables and
O(log n) equations, constructible in polynomial time in the binary target.
Every characteristic-zero solution is already an ordinary integer solution
and decodes actual tournament constructions. Consequently these two scalar
membership problems are in NP with their target given in binary.

**OPEN:** existence of witnesses for every desired odd target; coverage of the
unrestricted spectrum outside {7,21}; and the question of Diophantine
definability of arbitrary number-field integer rings. None is assumed here.

## 1. Inheritance and the useful change of question

The source is the question whether O_K is existentially polynomial-definable
inside its number field K. Its useful challenge is: what stops a fractional
solution from masquerading as an integral combinatorial witness? We answer
that challenge for the actual physical-source problem, rather than assume a
solution to the much broader number-field question.

Closest proved mechanisms are the full C2 source-pair theorem in
`dependencies/BILINEAR_C2.md`, the physical adjoint and midpoint inverse in
`dependencies/PHYSICAL_ADJOINT.md`, and the sharp source-order growth in
`dependencies/SOURCE_NORM_CAPACITY.md`. The canonical hostile is the
nonphysical integer profile (1,0,6,1), which passes every congruence-membership
test and formally gives C2 count 67 with the triangle partner, as proved in
`dependencies/PROFINITE_ESCAPE.md`. Thus being integral is not enough either.
The least-used coordinate is the complete construction word, with its height,
not another inequality on the four numerical coordinates.

Concept board: integer witnesses; physical words; midpoint height; finite-field
reduction; actual graph decoding. The map below preserves source membership,
the exact scalar, and a graph recipe. Passing only to the terminal profile
forgets the word and graph; the selector sequence is their restoration sidecar.
The cheap hostile tests are a non-Boolean fractional selector, a small positive
characteristic where a naive one-hot equation fails, and the phantom 67.

## 2. Physical grammar and a logarithmic word cap

A nontrivial rooted start-one source has profile s=(D,E,U,d), where D is its
Hamilton-path count, E its root-ending count, U its reserved two-path cover
count, and d the path count after root deletion. The seed is t=(1,0,2,1).
Its complete physical grammar is

    O(s)=(D,0,2U,d),
    Z(s)=(D-E+U,d,2U-d,d),
    J(s)=(D-E+U,D-E,D+U,D).

O appends a sink. Z appends a vertex beating only the old root and retains
that root. J performs the same graph birth and marks the newborn. Every word
is physical; every nontrivial source profile has such a word. Each genuine
letter adds one vertex. Add I(s)=s only as padding, not as a graph birth.

The inherited all-order growth bound for an order-m source is U^3>=5^(m-1).
Thus a genuine word of length l from t satisfies

    5^(l+1) <= U^3.                                      (1)

For an integer B>=2 define L(B) as the largest integer l>=0 with
5^(l+1)<=B^3. Every source with U<=B has a word of length at most L(B).
This is an exact integer definition, with no floating-point logarithm.
It follows that L(B)=O(log B). Also each physical step at most doubles U,
so every padded L-step word from t has U<=2^(L+1).

The complete C1 pairing is K(v,w)=v^T Q w, where

    Q=[[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]].

For any integer n>=3, C1 membership is equivalent to K(v,w)=n for two
physical sources with U_v,U_w<=B1(n)=floor(sqrt(4(n-1))). Hence set
L1(n)=L(B1(n)). If v has chronological word alpha and w has word beta, the
literal C1 graph is the source word beta alpha^dagger followed by J, where
dagger reverses the word and exchanges O,Z while fixing J. Its order is
|alpha|+|beta|+3. The inherited bound K(v,w)<=U_v U_w is useful below.

For C2, n<13 is impossible. The exceptional count 15 has its regular-five
witness. At every n>=13, all remaining C2 counts are N(v,w) or M(v,w), with
both physical sources bounded by B2(n)=floor((n-1)/6). Put L2(n)=L(B2(n)).
Writing v=(D,E,U,d), w=(g,s,A,b), the forms are

    N=(D-E+U)g-(D-E+d)s+(D+3U-d)A-(U+E)b,
    M=(D+2U)g-(U+d)s+(2D-E+6U-d)A-(U+E)b.               (2)

Their inherited positive-coordinate matrices, in x=(D-E-d,U-D,E,d), have
maximum entries 6 and 11 and sum(x)=U. Thus

    N(v,w)<=6 U_v U_w,   M(v,w)<=11 U_v U_w.            (3)

The bilinear theorem supplies actual N/M grafts and completeness, not just
these scalar identities. No arbitrary point of the nonnegative cone is used.

## 3. A quadratic compiler over every field

Fix L>=0. For each of the two sources introduce states s_0,...,s_L in F^4.
Set s_0=t. At each step i introduce four selectors e_(i,a), a in {I,O,Z,J}.
Impose the polynomial equations

    e_(i,a)^2-e_(i,a)=0                         for every a,
    e_(i,a)e_(i,b)=0                           for every a<b,
    sum_a e_(i,a)=1,
    s_(i+1)=sum_a e_(i,a) A_a s_i.                    (4)

Here A_a is the displayed integer linear matrix of operation a. There are
four coordinate equations in the last line. All equations have total degree
at most two and coefficients bounded absolutely.

**Theorem 1 (integrality by construction).** Over any field F, every solution
of (4) selects exactly one operation at each step and has all state coordinates
in the prime subfield. Its states are precisely the reductions of the integer
profiles of the selected literal words. Conversely every two padded words give
a solution. Over every characteristic-zero field, all state and selector values
are ordinary integers, embedded via Z->F.

**Proof.** A field has no zero divisors, so each Boolean selector is 0 or 1.
Pairwise orthogonality allows at most one 1, and their sum forces one. Induct
from the integral initial state through the deterministic recurrence. Deleting
I gives an actual physical source word. The converse is direct substitution.
QED. The pairwise equations are redundant in characteristic zero but deliberate:
in F_2 or F_3 the Boolean-plus-sum constraint alone would also admit three or
four selected operations, respectively. We work over fields, not arbitrary
rings with nontrivial idempotents.

For C1 append K(s_L,t_L)=n to (4), using L=L1(n). For C2 use L=L2(n), introduce
h_N=N(s_L,t_L), h_M=M(s_L,t_L), and three more one-hot selectors c_N,c_M,c_15
with the same Boolean, orthogonality and sum constraints. Append

    n=c_N h_N+c_M h_M+15 c_15.                        (5)

All equations remain quadratic, including (5). The source words may be arbitrary
when c_15=1; delete them and return the fixed regular-five graph.

**Theorem 2 (field-invariant exact recognition).** For every n>=3 and every
characteristic-zero field F, the C1 system has an F-point exactly when n is
a strong C1 count. For every n>=13, the C2 system has an F-point exactly when
n is a strong C2 count. Every point decodes an actual graph witness.

**Proof.** Soundness is Theorem 1 and the physical K/N/M construction theorems.
For completeness choose the inherited bounded-source representation, use (1)
to bound each word length, and pad with I. No inequality U<=B needs encoding:
additional length-bounded words are still actual sources and cannot falsely
represent a target. Equation (5) handles 15 separately. QED.

For a number field K this gives precisely the same answer over Z, Q, O_K,
or K. The claim over O_K follows by viewing an O_K-solution in K. Nothing
assumes that arbitrary K-points of an integer predicate descend to Q, nor that
integral numerical source profiles are physical. In this compiler, all the
witnesses themselves descend because their entire deterministic histories are
retained.

Before the terminal readout, the solution set has exactly 4^(2L) selector
histories, each with a unique state lift. It is a finite split reduced algebraic
set. The C1 terminal equation selects a subset; for C2 the additional gate
selector selects one of three scalar readouts. Field extension creates no new
points in characteristic zero. Padding means histories can describe the same
graph; no injectivity from histories to graph classes is asserted.

For the inherited balanced C1 subinventory, start the second source at
b=(3,1,3,1) instead. The same compiler works, with the safe common bound L(B)
(the sharper bound for words starting at b is one letter smaller). For a fixed
root-change budget r, add r Boolean slack variables z_j and the equation

    sum_(both words,all i) e_(i,J) + sum_(j=1)^r z_j = r.

In characteristic zero this enforces at most r J letters exactly. The inherited
physical adjoint preserves the total J cost. Treat the separate balanced
value 3 independently. We do not extend this slack equation's interpretation
to small positive characteristic without a further range bound.

## 4. Polynomial-bit certificates and a single finite-field check

Each source contributes 4(L+1) state variables and 4L selector variables, and
15L+4 equations. The C1 terminal adds one equation. Therefore the two-source
system has 16L+8 variables and 30L+9 equations. C2 adds only a constant number.
Its coefficient bit size is O(log n), entirely from the target. The whole
compiler is polynomial-time in log n and has O(log n) variables and equations.

A certificate is just the two words and, for C2, the readout choice. All states
are O(L)-bit integers, each step uses bounded integer arithmetic, and the
terminal form is directly computable. This proves membership of both C1 and
C2 spectra in NP for a binary target. It does not give a polynomial-time
search algorithm. The exact scalar count is verified by physical transfer
identities rather than a #P-hard count on an arbitrary input tournament.
Actual graph construction from the word has O(log n) vertices and polynomial
output size. This complexity claim concerns these classified C1/C2 families.

In the t,t compiler (not the balanced b,t variant), at fixed word cap L,
every C1 readout lies between 0 and

    M1(L)=4^(L+1),

and every C2 readout, including the exceptional 15, lies between 0 and

    M2(L)=max(15,11*4^(L+1)).

**Corollary (height-synchronized modular exactness).** Fix the appropriate
n,L and choose any prime p>max(n,Mj(L)). Then the compiled system has an F_p
solution exactly when it has a characteristic-zero solution.

Indeed Theorem 1 still decodes literal integer words from an F_p point. The
readout is congruent to n modulo p. Both its integer value and n lie in
[0,p), so they are equal. The converse is reduction modulo p. A prime above
the bound always exists. Its bit length can be O(log n), for example by
Bertrand's postulate; no efficient deterministic prime-search claim is needed.

For the optional balanced b,t variant a safe corresponding scalar upper
bound is 6*4^L, since its two seed cover counts are 2 and 3; use that bound
instead of M1.

This is a finite, target-dependent exactness statement. It does not say a fixed
prime, a fixed word cap, or unbounded modular source membership decides every
target. Proving absence of an F_p solution may still require large search, so
this is not a coNP certificate. In particular it does not contradict the
all-moduli phantom profile: there the required physical witness heights grow
without bound, whereas here word height is fixed before the prime is chosen.
The system distinguishes actual source words from a free integral quadruple.

## 5. What the number-field question contributes, and what it does not

An existential definition of O_K inside K could enforce integral arithmetic in
an arbitrary encoding. It would not, by itself, enforce membership in the
physical source semigroup, supply a witness height bound, or construct a pair
at each requested Hamilton count. The phantom is already an integer profile.

The constructive transfer is stronger in the restricted setting: use Boolean
construction choices to force integrality, use exact physical transfer laws to
force realizability, and use the midpoint/order estimates to keep the entire
certificate at logarithmic length. This avoids importing an unresolved outside
theorem and exposes the remaining universal coverage assertion as the precise
existence problem for this special family of quadratic systems.

We retain a system of equations, not a sum of their squares. The familiar
conversion sum f_i^2=0 is not sound over arbitrary number fields containing i.
If a single polynomial over a particular number field is required, choose a
finite extension with degree at least the number of equations and linearly
independent coefficients alpha_i; the field norm of sum alpha_i f_i vanishes
exactly when every f_i vanishes. That optional replacement changes the degree
and coefficients. It is not part of our uniform quadratic compiler.

The number of variables grows with n. This is a uniform algorithm producing a
family of systems, not one fixed existential formula defining the entire
unbounded spectrum. For each fixed L its finite image is Diophantine; a single
formula for all lengths would require an additional sequence-encoding theorem.
No such theorem is assumed or claimed.

The useful new problem is therefore specific: can the source grammar and the
positive pairing force solutions of these bounded, field-invariant systems
for every odd target in the desired range? Merely changing Q to a larger
number field cannot help, since it adds no witnesses. New physical words,
physical operations, or exact target-level arguments are required.

## 6. Verification and scope

`python3 verify.py` and `python3 -O verify.py` run independent integer state
recurrences, explicit polynomial-equation evaluations, finite-field selector
checks (including hostile naive selectors), direct graph path/cycle counts on
small pair witnesses, finite complete source banks, and target controls 7,21,67.
The checker also emits an explicit C1 target witness and equation-count data.
Finite validation supports the algebra; completeness depends on the frozen
proved source theorems and cutoffs, not a census extrapolation.
