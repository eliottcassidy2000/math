# All-modulus context clocks, prime-function states, and the actual 5/11 projector pair

**Status: PROVED (all-modulus theorems and the 5/11 construction); FINITE-EXACT (complete contextual-state realization at primes 5, 7, 11).** The exact integer Hamiltonian-path spectrum, bounded exact spectral width, and realization of every formally allowed contextual state at every prime remain OPEN. No external priority claim is made.

## 0. Inheritance, sources, and the objects retained

The source baseline is `291b4d9757842c02aa451beb65acfade2848609b`. Main dependencies:

* Session137 `05-knowledge/results/session137_overnight_proof_frontiers_20261102/width/PROOF.md`, Theorems1–2: the exact universal context coefficient divisor `d_c=c!(c-1)!`, actual sharp detector `E_(c-1)`, and exact modular contextual equivalence.
* The same packet `width_phase2/PROOF.md`, Sections1–3: the path-cover polynomial, its parity, and actual positive `P,C,D,J` operations. Its all-modulus **scalar** universality theorem is inherited, not reproved or upgraded to integer realization here.
* Session138 `05-knowledge/results/session138_stability_and_actual_contexts_20261102/exact_resources/PROOF.md`: actual order-five H=11/13 reflector pair and full-profile, rather than scalar-only, continuation boundary.
* Older `01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md`: adjacency zeta equals1 on every transitive tournament and is insensitive to singleton SCC padding. This old information-loss boundary motivates, rather than proves, our context clock.
* Older `01-canon/theorems/THM-4093-rational-edge-diagonal-gauge-and-padic-tournament-zeta-tangent.md`: diagonal-gauge blindness; it does not assert contextual equivalence. No RH theorem is imported.

A nonempty tournament T has `p_T(c)` unordered covers by c nonempty directed paths, and

    F_T(x)=sum_(c>=1) p_T(c)(x)_c,       (x)_c=x(x-1)...(x-c+1).

For tournaments T,U, write `T ~_M U` if every homogeneous one-port tournament context Q has equal Hamiltonian-path counts modulo M after substitution of T and U. The identity context is included. The inherited sharp theorem is

    T ~_M U  iff  d_c p_T(c)=d_c p_U(c) mod M for every c>=1.     (1)

This is about all finite homogeneous contexts, not merely finitely sampled ones. The profiles have finite support. We also use

    F_T(-x)=(-1)^|T| F_T(x),     F_(T->U)=F_T F_U,
    F_(L_n)(x)=x^n.

Our board is: actual source/branching resource; operational rather than scalar arithmetic; 5/11 factors; p-adic transients; old zeta loss. The canonical hostile is `H(L_n)=1` and `zeta(L_n)=1` at every n, although contexts detect n. The corrected near miss is replacing (1) by unrestricted ordinary polynomial-function equality at composite moduli.

## 1. A function model valid modulo EVERY modulus

Define the integral factorial transform

    Psi_T(x)=sum_(c>=1) (c-1)! p_T(c) (x)_c.                     (2)

For a nonnegative integer x this counts path covers, a cyclic ordering of their c distinguishable path blocks, and an injective assignment of x colors to those blocks. The cyclic-order factor is `(c-1)!`. This supplies an actual combinatorial interpretation; Psi is not a signed subtraction of tournaments.

**Theorem1 (all-modulus operational function).** For every integer M>=1 and nonempty T,U,

    T ~_M U  iff  Psi_T(x)=Psi_U(x) mod M for every x in Z/MZ.    (3)

**Proof.** In the falling-factorial basis, the c-th forward difference at zero is

    Delta^c Psi_T(0)=c!(c-1)! p_T(c)=d_c p_T(c).

If the functions agree on Z/MZ, their integer-coefficient polynomials agree at every nonnegative integer modulo M, so every forward difference agrees. Conversely, the finite Newton expansion recovers Psi at every nonnegative integer from these forward differences using integer binomial coefficients. Combine with(1). No division by a nonunit is used. QED.

Thus the exact all-modulus context quotient is a concrete finite family of special polynomial functions. This does not assert that arbitrary functions or all formal signatures are realizable by tournaments.

**Composite hostile.** `L4 ~_16 L6` by Theorem3 below, but `F_L4(3)=1 mod16` and `F_L6(3)=9 mod16`. Ordinary F-functions need not coincide on contextually equivalent tournaments. The factorial transform is essential at prime powers; one cannot extrapolate the prime simplification in the next section.

## 2. At every odd prime, ordinary profile functions suffice

**Theorem2 (prime-function theorem and upper bound).** Let p be an odd prime and `d=(p-1)/2`. Then

    T ~_p U  iff  F_T(x)=F_U(x) for every x in F_p.               (4)

Every contextual class belongs to the union of the even and odd function spaces

    V_epsilon={f:F_p->F_p : f(0)=0, f(-x)=epsilon*f(x)},
    epsilon in {+1,-1}.

Consequently there are at most `2 p^d - 1` nonempty tournament contextual classes modulo p.

**Proof.** For c<p, d_c is a unit modulo p. For c>=p, d_c vanishes modulo p. Hence(1) is equality of the first p-1 profile coordinates. At a representative x in `{0,...,p-1}`, every `(x)_c` with c>=p vanishes. The triangular evaluation matrix at x=1,...,p-1 has diagonal c!, all units. Thus those p-1 coordinates are equivalent to the F-function. Parity places the function in V_epsilon. Each space is freely described by its d values at x=1,...,d, so has size p^d; their intersection is only zero, because 2 is invertible. QED.

Parity tags are not additional contextual data at the zero function. Counting both zero tags would incorrectly add one class.

## 3. Exact transitive context clocks for every modulus

Define for a prime power:

    ell(p^a)=(p-1)*p^(a-1)                       if p is odd;
    ell(2^a)=1                                  if a<=2;
    ell(2^a)=2^max(1,a-3)                        if a>=3.

The dyadic value is exactly `ord_(2^(a-1))(3)` when a>=2, with the trivial modulus convention at a=1. For M>1 put

    A(M)=max_(p^a || M) a,
    L(M)=lcm_(p^a || M) ell(p^a).

For M=1 put A=L=1.

**Theorem3 (exact clock, transient, and number of states).** For all positive n,m,

    L_n ~_M L_m
    iff n=m or [n,m>=A(M) and n=m mod L(M)].                    (5)

In particular the first periodic index is exactly A(M), the least eventual period is exactly L(M), and the transitive family realizes exactly `A(M)-1+L(M)` contextual classes modulo M.

### Proof at odd prime powers

The transitive profile is `p_(L_n)(c)=S(n,c)`, a Stirling number of the second kind. The elementary inclusion-exclusion formula gives

    q_c(n):=d_c S(n,c)
      =(c-1)! sum_(j=0)^c (-1)^(c-j) binom(c,j) j^n.             (6)

The convention is n>=1, so 0^n=0. If n>=a and t is a multiple of phi(p^a), then every p-unit j has `j^(n+t)=j^n mod p^a`, while every p-multiple j has both powers zero modulo p^a. Equation(6) proves sufficiency.

Now suppose n<m and the signatures agree modulo p^a, and put t=m-n. For c<p, the factors d_c are p-units, so these coordinates recover `j^n` at j=1,...,p-1 by the triangular falling-factorial evaluation transform. Therefore

    j^t=1 mod p^a for all 1<=j<=p-1.                            (7)

Reduction modulo p forces p-1 to divide t, because all nonzero field residues occur. Using j=p-1, binomial induction on `(1-p)^t` gives

    v_p((p-1)^t-1)=1+v_p(t)   for positive even t.

(The t here is even because p-1 divides it.) Hence p^(a-1) divides t, proving phi(p^a) divides t.

If a=1, n>=a holds automatically. If a>=2, use c=p in(6). For j=1,...,p-1 the power differences vanish modulo p^a by(7), and the j=0 term is zero. The remaining difference is

    q_p(m)-q_p(n)=(p-1)! [p^m-p^n] mod p^a.

It has valuation exactly n when n<a, since p^t-1 is a unit. Thus equality forces n>=a. This proves necessity and the exact transient. Notice that the c=p detector remains available when a>=2; it is not erased at the cutoff.

### Proof at dyadic powers

For a=1 all Hamiltonian-path counts are odd and the context signature is constant; the theorem says A=L=1. Assume a>=2. The c=2 coordinate is

    q_2(n)=2^n-2.

For n<m its difference has valuation exactly n. Thus contextual equality requires n>=a. Once both n,m>=a, the c=3 coordinate is

    q_3(n)=2[3^n-3*2^n+3].

At a>=3 this requires `3^t=1 mod2^(a-1)`, giving the claimed period. For a=2, c=3 is already invisible and the period is1.

For sufficiency, in(6) the c=1 coordinate is constant, and c=2 is constant once n>=a. At c>=3 the factor(c-1)! is even. Every even j has j^n=0 mod2^a. For every odd j, `j^ell=1 mod2^(a-1)`: the exponent of the odd-unit group is1 at modulus2,2 at modulus4, and2^(a-3) at modulus2^(a-1) for a>=4. Multiplication by(c-1)! supplies the additional factor2. Thus all q_c agree. The stated unit exponents and order of3 follow directly by repeated squaring: `v_2(3^(2^r)-1)=r+2` for r>=1, together with the base mod4 case. QED for powers of2.

Finally, congruence modulo M is equivalent to congruence at each prime-power divisor. Their transient requirements take a maximum; their period requirements take an lcm. This proves(5).

### A sharp operator identity on EVERY actual tournament

Write P(T)=K1->T. The transitive clock is the universal eventual clock of source padding on the entire contextual quotient, not just on its transitive subfamily:

    P^(A(M)-1+L(M))(T) ~_M P^(A(M)-1)(T)  for every nonempty T. (7a)

The index A(M)-1 and period L(M) are both minimal for a universal eventual identity. More precisely, if delta(T) is the least exponent with nonzero monomial coefficient in F_T, then the contextual padding sequence is L(M)-periodic from padding exponent `max(0,A(M)-delta(T))` onward. Its individual least period can be a proper divisor of L(M).

**Proof.** Write the integer polynomial `F_T(x)=sum_(k>=1) a_k x^k`. Padding by r multiplies it by x^r. The weighted c-th profile coordinate is therefore the integer linear combination `sum_k a_k q_c(r+k)`, using(6). When r+delta(T)>=A(M), every summand is L(M)-periodic by Theorem3. This is an identity of response polynomials; it does not construct signed or virtual tournaments. K1 has delta1 and gives exactly the sharp transitive clock, so a smaller universal index or period is impossible. QED.

### The real 5/11 interaction

For a,b>=1,

    A(5^a*11^b)=max(a,b),
    L(5^a*11^b)=4 * 5^max(a-1,1) * 11^(b-1).                  (8)

The extra factor5 is forced by11-1, not by a numerical resemblance: the mod11 unit field already carries a period10, hence a factor5. For example M=55 and M=275 have the same period20, but their transients are1 and2. In fact `L1 ~_55 L21`, while `L1` and `L21` are separated modulo275 by the c=5 detector. Explicitly, the inherited sharp nine-vertex context E4 has response row `(117,1368,5448,8064,2880)`; its responses on L1 and L21 are respectively117 and172 modulo275, differing by55. At M=25 alone, the same pair has equal lower unit-power responses but differs in q5 by a p-adic transient. This is a precise map from prime-power unit arithmetic to actual homogeneous-context observations.

**Zeta boundary.** Every L_n has H=1 and adjacency zeta1. Nevertheless its entire modular context behavior records the clock(5), and the full collection over all M distinguishes every distinct pair of orders: for n<m choose any prime-power exponent a>n. Thus scalar H and ordinary adjacency zeta lose an explicitly recoverable order coordinate. This is not a connection to zeros of the Riemann zeta function.

## 4. Literal H=5 and H=11 factors become complementary mod5 projectors

Let

    W5=C3[L2,K1,K1],          W11=C3[L2,L2,K1].

They are actual strong tournaments of orders4 and5, with H=5 and11. The inherited binary cyclic formula (or direct enumeration of their odd cycles) gives

    F_W5(x)=x^4+4x^2,
    F_W11(x)=x^5+8x^3+2x.                                    (9)

At residues x=0,1,2,3,4 modulo5, their value lists are

    W5:  (0,0,2,2,0),
    W11: (0,1,0,0,4).                                       (10)

Thus their nonzero supports on F5* are disjoint, not merely their H values. Because ordinal join multiplies the full F-functions,

    Z5=W5->W11        has order9, H=55, and F_Z5=0 on F5.

**Theorem4 (actual universal divisibility and projectors).** Every homogeneous tournament context containing Z5 as its distinguished module has Hamiltonian-path count divisible by5. Moreover the even-parity factors

    E1=W11->W11,                 E2=W5->W5->W5->W5

act, under ordinal join modulo5, as the idempotent coordinate projectors onto `{+/-1}` and `{+/-2}`, respectively. Their function products vanish and their function sums equal1 on F5* (and0 at0).

**Proof.** The first statement is Theorem2 applied to the zero function: all weighted profile coordinates in(1) vanish modulo5, hence every context response does. The projector statements follow from squaring/fourth-powering(10). They are actual positive tournament factors. The formal sum E1+E2 is only an identity of functions; no addition of tournaments is asserted. QED.

This makes the 5/11 connection operational: H=11 on an actual five-vertex module selects one profile channel modulo5, while the H=5 module selects its complementary channel after a positive power. The scalar H alone would miss this because the second evaluations are32 and100. The construction is unrelated, absent an additional map, to the older arithmetic cubic discriminant `4x(x^2-11x-1)` or its exceptional parameters `(11+/-5 sqrt5)/2`.

## 5. Complete prime state spaces at 5, 7, and 11; a branching boundary

The actual unary operations are

    P(T)=K1->T,
    C(T)=C3[T,K1,K1],
    D(T)=C3[T,L2,K1].

Their exact polynomial operators are

    PF=xF,
    CF=x^2 F+x[F(x+1)-F(x-1)],
    DF=(x^3-2x)F+(2x^2+x)F(x+1)+(x-2x^2)F(x-1).               (11)

Source padding and D reverse order parity; C preserves it. In each prime-function space represent the function by `(epsilon;F(1),...,F(d))`; values elsewhere are recovered using parity and F(0)=0. Formula(11) therefore supplies completely explicit matrices and an exact finite positive-constructor graph.

**FINITE-EXACT theorem.** The all-tournament contextual class numbers at p=3,5,7,11 are respectively

    5, 49, 685, 322101.                                     (12)

All classes have actual substitution-width-at-most3 witnesses. At p=5 and11, unary P/C/D words starting at K1 realize every nonzero function in both parity spaces and no zero function. One ordinal join of actual disjoint-support representatives supplies zero. At p=3,7 the unary words already supply both zero parity tags.

The exact exhaustive state counts, including both zero tags if present, are

| p | dimension d | unary tagged states | allowed tagged states | contextual classes |
|---|---:|---:|---:|---:|
|3|1|6|6|5|
|5|2|48|50|49|
|7|3|686|686|685|
|11|5|322100|322102|322101|

Upper bounds in(12) are Theorem2, not a finite tournament census. Lower bounds come from exhaustive reachability of actual operations on an exact finite quotient; each discovered state retains a finite positive word. Termination and completeness are checked by exhausting the queue.

For p=5,11 the determinants of C in the even/odd sectors, and D in either sector, are

    p=5:   det C_even=3, det C_odd=4, det D=4;
    p=11:  det C_even=7, det C_odd=5, det D=6.

P is invertible too. Hence every unary gate is invertible between its appropriate parity sectors, and the nonzero seed can never reach zero at **any depth**. This proves an all-depth constructor-resource obstruction, not merely failure of a bounded search. In contrast det C_odd=0 at p=3 and7, and explicit unary zero words are C and CDDCC.

At p=11, the chronological words

    CDPPPDCPCDD       (order23),
    PPCCDDDCCP        (order21)

give odd-parity evaluation vectors `(1,0,0,0,0)` and `(0,1,0,0,0)`. Their ordinal join gives zero. At p=5, the much smaller explicit Z5 from Section4 supersedes the generic orbit witness.

This is a strict full-context resource separation at primes5 and11: a grammar allowing only unary P/C/D contexts misses the zero contextual state, whereas one genuine two-factor join attains it. Scalar H=0 occurs already within the unary orbit and is **not** the missing state. All claims use actual constructors. Nothing here proves the analogous complete function-state statement at every odd prime or any composite modulus; those are open questions.

## 6. Replay, hostile controls, and loss ledger

Run from any working directory:

    python3 -B verify.py > normal.json
    python3 -B -O verify.py > optimized.json
    cmp normal.json optimized.json

The self-contained standard-library verifier imports no inherited research code. It checks the exact clocks at all M<=512 and seven additional prime powers (519 moduli, 70,588 direct period checks); 35,168 factorial-function comparisons across all moduli1..32 on actual tournaments; all 1099 labelled tournaments through order 5 for the prime-function signature equivalence; 225 actual graph/subset-partition profile checks for the three operators; exact exhaustive finite-state orbits at 3,5,7,11; literal full-integer-polynomial replay of the zero-factor words; and the small 5/11 projectors. It rejects the false dyadic Euler-period prescription and the false absence of transients. The all-M and all-n quantifiers come from the proofs, not the checked range.

Source -> target maps and losses:

* Actual T -> Psi modulo M preserves **exactly** all homogeneous-context H responses; loses the graph, order, integer magnitude, and nonhomogeneous observations.
* Actual T -> F modulo odd p has the same preservation, with the explicit prime restriction. The L4/L6 mod16 hostile blocks a prime-power extension.
* Transitive order n -> its context clock quotient loses only the precise residue/transient ambiguity stated in(5); all moduli together restore n.
* Actual factors W5,W11 -> support projectors preserve their full prime-function action. The scalar product55 alone does not explain universal contextual divisibility.
* Unary expression -> finite orbit retains an actual positive recipe. Zero-state attainability uses join as an extra constructor; no inverse operation or negative vertex count is hidden.
* Zeta -> context arithmetic is not a factor map: the transitive family has one zeta value but many context states. The same-session weighted cospectral family gives an independent, stronger fixed-order loss, but weighted examples are not counted as binary orbit witnesses here.

The next precise boundary is whether the union V_even union V_odd is entirely realized at every odd prime, and how prime-power compatibility restricts the factorial-transform functions. Scalar all-modulus density alone does not answer either question.
