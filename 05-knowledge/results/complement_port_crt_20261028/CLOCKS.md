# Signed Collatz source laws and clocks

2026-10-28. **PROVED elementary / FINITE-EXACT; same-agent audit.** These
statements concern declared residue-source measures and supplied cycles.
They do not prove global convergence or completeness of positive cycles.
The signed-sheet and forward-negation mechanisms are inherited from
[signed_collatz_pair_lifts](../signed_collatz_pair_lifts_20261014d.md).

For an odd source n and b in {+1,-1}, define

    a_b(n)=v_2(3n+b),
    U_b(n)=(3n+b)/2^a_b(n).

Exclude the single 2-adic zero-numerator point on either sheet when taking
valuations; it has Haar measure zero. All finite tests here use positive odd
integers, so their numerators are nonzero.

## S1. Shared sources give a negative covariance

Sample n from normalized Haar measure on the odd 2-adics. The two numerators
differ by2. Exactly one valuation is1, and for every j>=2,

    P(a_+=1,a_-=j)=P(a_+=j,a_-=1)=2^(-j).

For example n=1 mod4 makes a_-=1, while each further condition on a_+ is one
odd source class at doubled precision. The other parity class exchanges the
roles. The exact joint probability generating function is

    E s^a_+ t^a_- = st(s+t-st)/[(2-s)(2-t)].

Both marginals are geometric with mean2 and variance2, but

    E(a_+ a_-)=3,       Cov(a_+,a_-)=-1.

The repeated number2 does not identify this process with ternary fixed-point
offspring. The ternary variable has mean1, values0/3, and positive shared-copy
covariance. Here a common arithmetic source forces negative correlation
between the two signed valuation clocks. `arithmetic.py` verifies the entire
capped law on odd residues modulo2^K, K=2,...,14. With cap K, its covariance
is exactly -(1-2^(1-K))^2.

## S2. Every finite valuation word has an exact source cylinder

For a positive word (a_1,...,a_r), set A=sum a_i and construct B by

    B_0=0, A_0=0;
    B_i=3 B_(i-1)+b 2^A_(i-1),   A_i=A_(i-1)+a_i.

Then

    2^A U_b^r(n)=3^r n+B,
    n = (2^A-B)(3^r)^(-1) mod 2^(A+1).

The second condition is necessary and sufficient for precisely that valuation
word: backward substitution from an odd final value produces odd intermediate
values, with the prescribed powers of2 removed at each step. It is one odd
source residue class and has normalized odd-Haar mass2^(-A). Thus a finite
word keeps its affine intercept, source class, sheet and accumulated clock.
The checker uses both signs, every word of length1..4 with entries1..3, and
three integer lifts of each cylinder.

Summing a cylinder law over residues differs from sampling a single positive
integer orbit. The recent
[Syracuse finite-section result](../syracuse_transfer_finite_section_20261027.md)
also fixes its Haar/Ulam source law; its rank-one-plus-nilpotent finite matrix
does not establish that integer representatives have the same dynamics.

## S3. The cycle equation keeps the two clocks separate

A proposed periodic valuation word must give

    n=B/(2^A-3^r),

with n a positive integer in its source cylinder and no earlier return if a
primitive cycle is intended. For b=+1, B>0 requires 2^A>3^r; for b=-1, B<0
requires 2^A<3^r. These are necessary arithmetic constraints, not a classification
of all cycles. The supplied familiar cycles check as follows.

| Sheet | Odd states | Odd steps r | Shortcut steps A | 2^A-3^r |
|---|---|---:|---:|---:|
| +1 | 1 | 1 | 2 | 1 |
| -1 | 1 | 1 | 1 | -1 |
| -1 | 5,7 | 2 | 3 | -1 |
| -1 | 17,25,37,55,41,61,91 | 7 | 11 | -139 |

The shortcut map divides even n by2 and sends odd n to (3n+b)/2. The
unaccelerated map using a separate multiplication/addition step instead has
A+r steps around the same odd cycle. In particular 7 and11 in the last row
are two different clocks for one supplied cycle, not two lengths under one map.

Negation gives U_(-b)(-n)=-U_b(n), a forward conjugacy, not a reversal of time.
The shared endpoint U_+(3)=U_-(7)=5 does not merge the signed sources:
U_+(5)=1 whereas U_-(5)=7. Neither this seam nor root ordering proves that
the displayed positive cycles are all the cycles.

## S4. The common framework and its limits

Across ports, Lucas charts and signed Collatz, a legal reduction must preserve
the source, joint transition law, operation clock, and any side data that
later guards can inspect. The tournament source law is uniform group elements;
CRT can use uniform words or uniform ring states; S1-S2 use normalized odd-Haar
measure. These are explicitly different probability spaces.

The constructive ternary repair bound counts legal re-marking edges. CRT
periods count one common rotation. A Collatz valuation sum counts divisions
by2 under the stated shortcut convention. Comparing all three is useful only
after declaring a map that transports the target and charges its clock. No
dynamical conjugacy between the three systems is claimed here.

## S5. Incoming exact tree factor, and its signed shared-source coupling

Incoming [Session94 valuation trees](../session94_recursive_ports_crt_clocks_20261029/dynamics/VALUATION_TREE_CLOCKS.md)
constructs a stronger, precisely typed bridge. On odd 2-adics with the poles
and all their finite preimages removed, a single sign's successive valuation
digits are iid with P(a=j)=2^-j. Put

    b=1[a is even], G=floor((a-1)/2), a=1+b+2G.

Then b has Bernoulli parameter1/3, G has probabilities (3/4)(1/4)^g, and they
are independent. Breadth-first offspring xi=3b therefore gives exactly the
critical ternary fixed-point genealogy. Recovering the original arithmetic
source still needs the ordered residual G marks and unused valuation tail.
This is an exact factor of a source law, not a conjugacy of Collatz iteration
with tournament re-marking. The result and its tree-size/clock recursion are
credited to that incoming packet.

Our S1 law specifies how the **two signs at the same source** couple this
factor. Set xi_+=3*1[a_+ even] and xi_-=3*1[a_- even]. Then exactly

    P((xi_+,xi_-)=(3,0))=1/3,
    P((xi_+,xi_-)=(0,3))=1/3,
    P((xi_+,xi_-)=(0,0))=1/3.

Proof: sum 2^-j over even j>=2 to get1/3 on either signed arm of S1. The
remaining event has the larger valuation odd. The two coordinates cannot
both be3 because one original valuation is1. Thus the joint offspring PGF is

    (1+s^3+t^3)/3,
    E xi_+=E xi_-=1, Var(xi_+)=Var(xi_-)=2,
    Cov(xi_+,xi_-)=-1.

This preserves the critical marginal law while making shared-source dependence
explicit. It is a first-step joint law. After applying U_+ and U_-, the two
sources differ; no independent or identically coupled all-generation pair of
native trees follows without another construction.

Incoming [AFFINE A2](../weighted_crt_clocks_20261028/AFFINE.md) also supplies
the exact periodic Lucas bridge with the correct clock. Negate the positive
3n-1 seven-state cycle to get a negative 3n+1 cycle. Its standard clock is
7+11=18, and phi^18-1=76 phi^9 gives R_18=O/(76), of cardinality76^2=5776.
The quotient to R_9 has76-element fibres; it is not the76-state carrier itself.
Likewise the5<->7 minus-cycle has standard period5 after negation and belongs
to the inherited eleven-state rational-cycle dictionary. These are maps of
specified periodic sources, with sign and clock retained, not convergence
results for arbitrary positive integers.
