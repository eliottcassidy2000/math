# Marked fixed-degree series and exact recurrence orders

Status: PROVED from the explicitly inherited inputs; independent mathematical
audit passed, with separate computational audit supplied alongside.
Session23, 2026-10-21.
This extends the audited Session22 fixed-degree rationality theorem. It does
not classify the remaining strong atoms or centers and does not determine
metagraph adjacency, bridges, or connectivity.

## 0. Inheritance and scope

The source checkout was read at `bd1c1b9be938b93d8f1b7dfe4e4cf0d5f58e6593`.
The following are inherited, not new results here:

1. `fixed_degree_metagraph_rationality_20261021/session22-fixed-degree/PROOF.md`:
   finite strong alphabets, fixed-degree rationality, the ordered-factorization
   formula, center-divisor formula, and the positive-real leading pole.
2. `effective_black_defect_order_bound_20261020/PROOF.txt`: for a strong
   tournament with h-b<=d, its order is <16(d+2)^2(d+1).
3. `metagraph_classification_20261021a/README.md`, Sections 3-5: the only
   strong h=3 atom is R5; its signature is (order,h,b,|Aut|)=(5,3,3,5).
4. `THM-2454-pure-blue-alphabet-and-center-law.md`, Section 2 and its
   2026-10-20 supersession: unique self-converse decomposition U -> C -> U^op,
   with h=x^2 a and b=x beta. Empty C has a=beta=1.
5. `THM-2450-rigid-self-converse-classes-are-cyclic-ternary-towers.md`,
   Section 2 ONLY: H and |Aut| multiply across ordered SCCs. Its older tower
   frontier is not inherited; the current rigid alphabet is K1,C3, from
   `merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md`.
6. `deficiency_three_20261020/centers.md`: the exact deficit-divisor law.
7. `deficit_two_centers_20261020/session15-deficit-centers/EXACT_LEAF_SERIES.md`:
   the exact minimal recurrence order nine at black degree one is already
   known. It is a control here, not a new result.
8. `continuation_width_symmetry_20261020/README.md`, Section 3: the rigid
   bivariate series 1/(1-z-v z^3), coefficient binom(n-2c,c), and the
   C3^c automorphism law are already proved. Their specialization below is
   inherited; the extension to every fixed h/degree and the raw-order
   obstruction are the additions.
9. `two_shore_deficit_exclusions_20261021/session22-incoming-center-review/RECONCILIATION.md`:
   the trivalent denominator p^2 q^2 is already proved to be reduced, with
   exact minimal recurrence order eighteen. The r=1 case below is inherited.

All result paths above are relative to `05-knowledge/results/` unless they
start with THM, in which case they are relative to `01-canon/theorems/`.
All objects are ordinary tournament isomorphism classes until the explicit
converse quotient. Include the empty object in h=1 with order and SCC count
zero, H=|Aut|=h=1. A non-self-converse class has b=0. Q is the fundamental
paired-black converse-merged graph, not its doubled-edge version. Its degree
is h on a non-self-converse pair and (h-b)/2 on a self-converse class.

Inheritance board: ordered SCC words / involutive reversal / finite strong
alphabets / additive valuations of multiplicative statistics / irreducible
pole carriers / lacunary failure of raw multiplicative marking. Canonical
hostile: R5 -> R5^op has (h,b)=(9,3), not (9,9). Corrected near miss:
MISTAKE-647 forbids rediscovering the old Narayana series as new mathematics.
The least-used sidecar here is the direct-product automorphism law, rather
than only h multiplicativity. No literature priority claim is made.

## 1. A finite multivariate refinement theorem

For a tournament T let c(T) be its number of SCCs. Fix a finite set of primes
Pi, and put

    W(T) = z^|T| u^c(T) product_(p in Pi) v_p^nu_p(|Aut(T)|).

All exponents are nonnegative integers. The variables v_p are formal markers,
not valuations themselves. The direct-product automorphism law makes W
multiplicative across ordered SCCs. It is invariant under converse. If 3 is
not in Pi, interpret v_3 as 1 in formulas below. Write W^2-substitution for
simultaneously squaring z,u,and every v_p, and define

    D(W) = 1 - u z - u v_3 z^3,
    A(W) = 1/D(W).

For odd j>1 let P_j(W) be the sum of W(C) over the actual strong ordinary
atoms with h(C)=j. Let C_a^+(W) be the analogous sum over strong self-converse
centers with h=a, plus the empty center if a=1. These are finite polynomials.
Let F_j(W) count all ordinary classes with h=j, and let J_j(W) count the
self-converse ones. Then

    F_j(W) = sum_(j_1...j_r=j, j_i>1)
                 A(W)^(r+1) product_i P_(j_i)(W),             (1)
    J_j(W) = sum_(x^2|j) F_x(W^2) C_(j/x^2)^+(W).             (2)

The ordered tuple is empty only when j=1. Every tuple has length at most
Omega(j). In (2), x is positive odd. No division by two occurs in (2): the
first half of a self-converse SCC word is canonical. Squaring ALL additive
markers is essential; substituting z^2 alone undercounts both SCCs and the
prime valuations contributed by the opposite half.

For d>0 even and positive odd x, define

    E_(x,d)(W;y,t) = sum_C W(C) y^(x^2 a) t^(x beta),

where C is empty or strong self-converse, a=h(C), beta=b(C), and
x(x a-beta)=d. For x=1 the positive-deficit bound makes this polynomial finite.
For x>1 the inequalities beta<=a and a>=1 give

    x|d, x(x-1)<=d, a<=d/[x(x-1)],

so E is again a finite polynomial. Only actual centers are counted, never
all integer signatures. The fully refined fixed-deficit series is

    S_d(W;y,t) = sum_x F_x(W^2) E_(x,d)(W;y,t).               (3)

The sum is over x=1 and the finitely many x>1 satisfying the displayed
bounds. It marks order, SCC count, the chosen automorphism valuations, h,
and b. In particular h and b take only finitely many values on each fixed
positive deficit slice, a fact not asserting that |Aut| or H is bounded.

The corresponding fixed-positive-degree merged-node series is

    B_k(W;y,t) = S_(2k)(W;y,t)
                  + 1_(k odd) y^k (F_k(W)-J_k(W))/2.         (4)

Here the second term has t exponent zero. Every statistic being marked is
converse-invariant, so each non-self-converse pair has equal weights. Thus
(4) has nonnegative integer coefficients despite the factor 1/2. Equations
(1)-(4) prove multivariate rationality. Their only nonpolynomial denominator
carriers are D(W) and D(W^2). The inherited pure-blue alphabet also gives
the fully marked degree-zero boundary explicitly (interpret missing prime
markers as one):

    B_0(W;y,t) = [y t(1+u z+u v_3 z^3)
                     + y^3 t^3 u v_5 z^5]/D(W^2).          (4a)

For comparison, the converse-merged series at fixed h=j, without fixing
black degree or marking b, is (F_j(W)+J_j(W))/2, not (4).

The same proof permits any finite list of nonnegative, additive-across-SCC,
converse-invariant integer statistics. In particular prime valuations of H
can be included: H=h|Aut|, or directly by SCC multiplicativity. For fixed k,
there are only finitely many possible primes dividing Aut or H: all
nonrigid atoms and centers lie in a finite union of finite alphabets, and
the only unbounded rigid contribution introduces the prime 3. Hence an
appropriate finite Pi can encode the entire multiplicative data without
loss, although that prime set need not be explicitly computed yet.

Converse-variant statistics require a changed formula: the paired weight is
W(C)W(C^op), not W(C)^2. Consequently merely squaring markers is not valid
for an oriented statistic such as a source-sensitive score coordinate.

## 2. Exact denominator of ordinary and self-converse strata

Put p(z)=1-z-z^3 and q(z)=p(z^2)=1-z^2-z^6. Both are irreducible in Q[z], and
they are coprime. For p, the rational-root test on z^3+z-1 suffices.

For q, suppose f(z)=z^6+z^2-1 is reducible. Let alpha be a root and put
beta=alpha^2. The minimal polynomial of beta is z^3+z-1. Therefore
[Q(alpha):Q] is 3 or 6; reducibility forces degree 3 and alpha belongs to
Q(beta). Its monic minimal polynomial is g(z)=z^3+a z^2+b z+c with integer
coefficients, since alpha is an algebraic integer. The distinct polynomials
g(z) and -g(-z) multiply to f(z). Comparing coefficients gives

    2b=a^2, b^2-2ac=1, c^2=1.

The first equality forces a=2e and b=2e^2. The second then has left side
divisible by four, impossible. This proves irreducibility. Finally if
p(z)=0, then q(z)=2z(1-z), which cannot vanish at a root of p.

Let r(j) be the maximum number of nonrigid strong SCCs in any class with
h=j. If no such class exists leave r(j) undefined. Equivalently r(j) is the
maximum length of an ordered tuple in (1) whose atom polynomials are all
nonzero. In particular r(1)=0, and r(j)<=Omega(j).

**Exact denominator theorem.** If F_j is nonzero, its reduced denominator
over Q[z], normalized to constant term one, is exactly

    p(z)^(r(j)+1).                                         (5)

Indeed (1) supplies this denominator bound. At the positive root rho of p,
all top-length numerator products are positive. Since p has a simple root
at rho, this gives pole order r(j)+1 there. Irreducibility forces the whole
factor p to survive to that exponent, and there are no other factors.

For nonempty self-converse deficit d, let

    s(d) = max {r(x): E_(x,d)(z;1,1) != 0 and F_x != 0}.

Then the exact reduced denominator is

    denominator(S_d) = q(z)^(s(d)+1).                      (6)

The same positivity argument at sqrt(rho) applies, followed by irreducibility
of q. If S_d is empty its series is zero, with denominator one. Likewise,
for J_j nonempty, with

    t(j) = max {r(x): x^2|j, F_x != 0, C_(j/x^2)^+ != 0},

its exact reduced denominator is q^(t(j)+1). Thus the non-self-converse
pair series (F_j-J_j)/2 has exact denominator p^(r(j)+1)q^(t(j)+1) whenever
J_j is nonempty, or p^(r(j)+1) if J_j is empty. Coprimality prevents
cross-cancellation. For even k, (6) is already the exact fixed-degree
answer. For arbitrary odd k, the q poles in S_(2k)-J_k/2 may cancel; no
unproved general positivity argument is used to forbid that.

These are effective classifications in terms of actual finite alphabets,
not formulas pretending unknown atom inventories have been computed.

## 3. Exact minimal recurrence order at every degree 3^r

**Theorem.** For every integer r>=0, B_(3^r)(z) has reduced denominator

    p(z)^(r+1) q(z)^e(r),
    e(1)=2, and e(r)=floor(r/2)+1 for r!=1.                 (7)

Its minimal eventual homogeneous constant-coefficient recurrence over Q
therefore has order 3(r+1)+6e(r). The cases r=0 and r=1 recover the inherited exact leaf
and trivalent orders nine and eighteen. The all-r classification is the
new extension, not either special case. The first orders for r=0,1,... are 9,18,21,24,33,36,45,48,... .

Proof. The inherited strong h=3 classifier gives P_3=z^5. The only possible
length-r factorization of 3^r uses r copies of R5, so r(3^r)=r and

    F_(3^r)(z) = z^(5r)/p(z)^(r+1) + lower-p-order terms.   (8)

Put m=floor(r/2). In the center decomposition for J_(3^r), the unique
largest half exponent is x=3^m. The center is one of empty,K1,C3 if r=2m,
or R5 if r=2m+1. Therefore its exact top-q numerator is

    L_r(z)=z^(10m)(1+z+z^3)     if r=2m,
    L_r(z)=z^(10m+5)           if r=2m+1,                  (9)

meaning J_(3^r)=L_r/q^(m+1) plus terms with denominator dividing q^m.
The lower terms have integer polynomial numerators after clearing q^m.

For d=2*3^r, any feasible odd x divides d and thus equals 3^t. If r=2m,
t>=m+1 violates x(x-1)<=d. If r=2m+1 with m>=1, at t=m+1 the excess is

    3^(m+1)(3^(m+1)-1)-2*3^(2m+1)
       =3^(m+1)(3^m-1)>0.

Hence t<=m for every r>=2; for r=0 only x=1 is allowed. By (1), S_d has
denominator dividing q^(m+1). Write its top-q coefficient as an integer
polynomial M_r, allowing it to be zero. Thus the q-leading numerator of
2S_d-J_(3^r) is 2M_r-L_r modulo q.

It cannot vanish modulo q. Otherwise Gauss's lemma for the primitive monic
(up to sign) q would give q|(2M_r-L_r) in Z[z]. Reducing modulo 2 would
force

    (1+z+z^3)^2 | L_r(z) in F_2[z].

But z is coprime to 1+z+z^3, and (9) contains either exactly one copy of
1+z+z^3 or none. This is impossible. Thus the q exponent is m+1 for
r=0 or r>=2. The factor 1/2 does not affect poles over Q.

For r=1, the inherited exact deficit-six law gives

    S_6(z) = R_6(z)/q(z) + z^10(1+z+z^3)/q(z)^2,
    J_3(z) = z^5/q(z).

The q^2 numerator is positive at sqrt(rho), so the q exponent is two.
Finally F contributes the p exponent r+1, and q terms cannot cancel it.
This proves (7). A rational series with reduced constant-one denominator
has precisely that denominator as its minimal eventual characteristic
polynomial (with coefficients reversed in the usual shift convention):
an eventual recurrence of lower order would supply a lower-degree
annihilating denominator, contradicting reduction.

No atom with h>3 and no unknown positive-deficit center coefficient is
needed to determine these exact recurrence orders. Such data still affect
the numerators, initial values, and recurrence onset.

## 4. A sharp marking boundary: raw automorphism order is not D-finite

Valuation marking is not interchangeable with marking the multiplicative
integer as an ordinary exponent. Already on the h=1 stratum, define

    R(z,t)=sum_(h(T)=1) z^|T| t^|Aut(T)|.

Since every such class is an ordered K1/C3 word, a word with c triangles
has automorphism order 3^c. Summing over any number of singleton components
and their placements gives

    R(z,t)=sum_(c>=0) t^(3^c) z^(3c)/(1-z)^(c+1).           (10)

This is NOT D-finite in t over K=Q(z), and in particular is not rational
in z,t. Here D-finite means satisfying a nonzero linear differential
equation in t with coefficients in K[t]. Coefficient extraction from such
an equation gives a nonzero finite-width polynomial-coefficient recurrence
for its t coefficients a_n, of the form

    sum_(j=0)^L P_j(n) a_(n+j)=0,     P_L not identically zero.

For sufficiently large integer n, P_L(n) is nonzero in K. If L=0 this
forces eventual vanishing. If L>0, any block of L consecutive zeros beyond
that threshold propagates forward forever. But in (10), each coefficient
a_(3^c)=z^(3c)/(1-z)^(c+1) is nonzero and all other coefficients vanish.
There are arbitrarily long zero blocks followed by nonzero coefficients,
a contradiction. The same obstruction holds when marking H instead of
|Aut| on h=1, since H=|Aut| there.

The empty tournament contributes t, not 1, because its automorphism group
has order one. In contrast, valuation marking on this same stratum gives
1/(1-z-v_3 z^3). Both claims are exact. No claim about differential finiteness in z after
a numerical specialization of t is made.

## 5. Connection contract, hostiles, and proof boundary

Source: unique ordered SCC words with involutive converse.
Target: weighted rational series and exact irreducible denominators.
Map: multiply additive-statistic monomials, then pair a half with its
opposite and an actual center. Preserved: ordinary and converse-merged
class counts, h,b, SCC count, selected Aut/H valuations, and black degree.
Lost: labelled Hamilton words, endpoint choices, targets, edges and wiring.
Needed sidecar: literal source/word/target ports from the earlier classifier.
Decisive tests: literal involutive-word enumeration versus independent
polynomial convolution; square all markers; actual R5/S4 controls; q mod2
noncancellation; raw automorphism exponent gaps.

The bridge to older multiplicativity work is exact, not numerical: direct
products become sums under prime valuation. The bridge to the old center
law is a multivariate diagonal-square substitution. The distinction between
ordinary exponent and additive logarithmic coordinates explains the
nonrationality boundary rather than merely recording a failed fit.

Hostiles retained: nonrigid pure-blue half R5 is not rigid padding; empty
center differs from K1; ordinary words are not converse orbits; automorphism
order is not SCC-additive; dual-variant markers cannot be blindly squared;
synthetic atom signatures are algebraic controls and not realizability
claims; q-order cancellation in general odd degree is left open rather
than argued from positivity of the final coefficient sequence.

Computational controls and source hashes are supplied alongside this note.
The main checker compares 16,601 literal formal involutive words through
order 13 against independent weighted polynomial convolution, retaining
order, SCC count, nu_3(Aut), nu_5(Aut), h, and blue count. It checks fixed
degrees 1 through 20, rational denominators for r=0 through 7, and the
mod-two obstruction for r=0 through 40. Its invented higher atom signatures
are explicitly synthetic, not a tournament census. Normal and optimized
runs must agree exactly. The independent reviewer uses separately written
word and polynomial controls.
All computation is supporting exact finite evidence, not a substitute for
the all-order arguments. No probabilistic limit-law claims are made here.
