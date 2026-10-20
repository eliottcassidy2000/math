# Observable quotients, Hankel precision, and a recovered binary clock

2026-10-20. **PROVED structural statements; FINITE-EXACT Smith data for the stated integer sequence.** This extends an underused recurrence packet through the plane-coloring blueprint's finite-quotient mechanism. It does not transfer the blueprint's analytic or coloring conclusions to an LFSR, and makes no historical-novelty claim for realization theory or Smith normal forms.

## 1. The paper's missing bridge to LFSRs

In the attached plane-coloring blueprint, §8 constructs a finite-dimensional map Φ from a functional's Cauchy data, proves that ker Φ is stable under multiplication by T, and then uses the characteristic polynomial on the quotient. Two separate facts are essential:

1. the quotient is finite-dimensional;
2. the forgotten subspace is stable under the actual shift/multiplication.

A low-dimensional list of measurements alone gives neither. The natural finite-state object is the observable quotient, not an arbitrary finite sample of values.

The closest repository precedents are [THM-3484](../../../01-canon/theorems/THM-3484-ternary-weighted-determinant-minimal-recurrence-and-cubic-fourier-degree-drop.md), its [Fourier/Jordan classification THM-3485](../../../01-canon/theorems/THM-3485-periodic-polynomial-fourier-jordan-recurrence-classification.md), and [THM-2810's faithful Hankel obstruction](../../../01-canon/theorems/THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go.md). We use these as linear recurrence and quotient mechanisms; no Gaussian-moment conjecture is reopened.

## 2. A finite exact observable module over every residue ring

Let a=(a_n)_(n≥0) be an integer sequence annihilated for all n by a known monic degree-r polynomial χ(E), where E is forward shift. Put

    H=(a_(i+j))_(0≤i,j<r).

For any integer modulus m≥2 let M_m be the module of sequences generated over Z/m by all forward shifts of a mod m. Then

    M_m ≅ image(H mod m) ⊂ (Z/m)^r.                       (1)

Under this isomorphism the first r values are a complete state observation.

Proof. The monic recurrence expresses every shift E^j a, j≥r, in the span of the first r shifts. A combination of those shifts is zero as a whole sequence iff its first r entries are zero: the same monic recurrence propagates those r zeros forward. The map from coefficient vectors to first r values is exactly H. Its kernel is therefore the true observation kernel, and the first isomorphism theorem proves (1). The kernel is E-stable because vanishing of the whole sequence is shift-stable. This is the explicit finite-quotient step that a finite fit alone cannot supply.

If χ(0) is a unit modulo m, E is invertible on M_m, since χ(E)=0 expresses E^(-1) as a polynomial in E. The finite state graph has one incoming and one outgoing edge at every state. If χ(0) is not a unit, this inference is unavailable; tails can feed cycles.

If the integer Smith diagonal of H has nonzero entries d_1,…,d_s and the rest zero, then

    |M_(p^k)| = p^( sum_(i=1)^s max(k−v_p(d_i),0) ),       (2)
    M_(p^k) ≅ direct sum_i Z/p^(max(k−v_p(d_i),0)).        (3)

A zero exponent denotes a trivial summand. This is the image of a diagonal scalar map modulo p^k, transported by the unimodular row/column operations. The module is not generally free over Z/p^k. Smith coordinates classify its additive structure, not the shift operator itself.

This is the smallest faithful linear state module for the span of the sequence's shifts. Its number of states is not the orbit length of the original sequence. A single orbit can be much smaller than the space of all linear superpositions.

## 3. A concrete old degree-22 word has a binary degree-six model

Normalize THM-3484's integer sequence by a_n=D_n/256. All coefficients below are integers:

    P0(n)=−8n^4+16n^5+32n^6−64n^7,
    P1(n)= 8n^4+16n^5−96n^6−64n^7,
    P2(n)=−n+12n²−60n³+160n^4−240n^5+192n^6−64n^7,
    a_n=P_r(n), r=n mod3.

The inherited globally proved minimal characteristic polynomial over Q is

    χ(x)=(x−1)^8(x²+x+1)^7,  degree22,  χ(0)=1.          (4)

Dividing by 256 does not change that rational minimal polynomial. It changes the integral lattice and therefore must be stated before modular reduction.

Modulo two, P0 and P1 vanish and P2(n)=n. Consequently

    a_n mod2 = 1 exactly when n=5 mod6.                  (5)

Its six cyclic shifts are the six coordinate delta words on Z/6, hence are linearly independent over F2. Therefore its minimal binary shift polynomial is

    x^6+1=(x+1)^2(x²+x+1)^2.                            (6)

The observable module has 64 states. It is conjugate as a linear dynamical system to the rotation of six binary coordinates. Its cycle inventory is

    2 cycles of length1,
    1 cycle of length2,
    2 cycles of length3,
    9 cycles of length6.                               (7)

These follow from Fix(E^d)=2^gcd(6,d) and divisor inversion, or direct enumeration. The original delta word has one six-cycle; it does not visit all64 linear states.

The order22→6 drop is not a numerical coincidence involving six. It is a proved change of observable quotient under reduction, with its exact integer lattice retained.

## 4. Smith precision reveals states invisible in a field reduction

The exact 22×22 Hankel Smith computation is reproduced by the companion script, with a separate local valuation-pivot algorithm. Its prime valuation lists are:

    p=2:  0 repeated6; 6 repeated6; 14 repeated2; 15 repeated4;
          20 repeated2; 21 repeated2.
    p=3:  0 repeated3; 2 repeated3; 4 repeated3; 8 repeated3;
          10 repeated3; 12 repeated3; 19,20,22,22.
    p=5:  0 repeated14; 2 repeated4; 3 repeated2; 4 repeated2.
    p=7:  0 repeated21; 8.
    p=61: 0 repeated21; 7.

All other primes have 22 zero valuations. The determinant is

    det H = −2^206 3^191 5^22 7^8 61^7,                 (8)

as also follows from THM-3484's determinant after dividing each of22 rows by256. The complete Smith diagonal, not merely its determinant, appears in the exact receipt.

Consequences of (2):

- The binary field sees dimension6. Through precision k≤6 the observable group is (Z/2^k)^6. At k=7 six additional order-two directions appear. Higher precision reveals the remaining layers at the displayed thresholds.
- At p=7 the field rank is21, despite determinant valuation8. The missing direction first appears at precision k=9 as an order-seven summand.
- At p=61 the field rank is also21, but the missing direction first appears at precision k=8.

Thus a determinant valuation is neither the number of lost state coordinates nor their recovery depth. The full Smith list distinguishes these. Conversely, Smith data alone do not identify the conjugacy class of the shift; the module action is another sidecar.

### A sharp six-register precision boundary

Through precision64 the entire sequence satisfies the order-six recurrence

    a_(n+6)=7a_n+48a_(n+1)+16a_(n+2)+2a_(n+3)
              +16a_(n+4)+48a_(n+5) mod64.              (9a)

The first6×6 Hankel determinant is odd, so these six coordinates are independent over every2-power ring. Checking the first22 residual values in (9a) is sufficient globally: the residual sequence inherits the known degree22 recurrence (4). Consequently the six-register model is exact for every modulus2^k with1≤k≤6, and invertible because its constant coefficient7 is odd.

At modulus128, the unique candidate six-register coefficients fitting the first six residual equations are

    (7,112,16,66,80,48).

It fails later residual equations, first at n=6. Thus the six-register model cannot simply be lifted by one more binary digit. The Smith state count makes the obstruction intrinsic: M_128 has size2^48, whereas six free coordinates over Z/128 have only2^42 states. A linear realization of the generated shift module must retain the additional directions. This does not forbid a nonlinear encoding of a single periodic orbit.

### The complete minimal-register staircase

The same finite Smith certificate determines the exact minimal degree of a monic recurrence modulo2^k, not just the first failure. The degrees are

    k=1,…,6:     6,
    k=7,…,14:   12,
    k=15:       14,
    k=16,…,20:  18,
    k=21:       20,
    k≥22:       22.                                    (9b)

Here is the extra certificate needed for this conclusion. For r=6,12,14,18,20,22, the rectangular22×r matrix formed by the first r Hankel columns has2-adic Smith valuations exactly the first r values of the full list. At a precision in (9b), every remaining full-Hankel valuation is at least k. Thus the prefix image and the full image have equal size, and the former is a subset of the latter. Those first r shifts therefore generate M_(2^k), giving a monic order-r recurrence.

Conversely, Nakayama's elementary finite-module count gives

    dim_F2(M_(2^k)/2M_(2^k)) = #{i:v2(d_i)<k}=r.

Any monic recurrence of smaller degree would generate the shift module with fewer than r elements, impossible. The constant coefficient of the degree-r recurrence is odd: its induced companion action on this r-dimensional mod-two quotient is the invertible shift. These are exact minimal linear recurrences, though at higher precisions the faithful state module need not be all of the free register ambient space.

A crucial distinction occurs at k=7: M_128/2M_128 has dimension12, whereas reducing its observed coordinates modulo2 has image M_2 of dimension6. The latter reduction is a different quotient. Confusing them would erase the newly visible directions.

## 5. The exact binary arithmetic clock and its completion

For every k≥1 the sequence a modulo2^k has exact least period

    P_k=3·2^k.                                         (9)

Proof. Since each P_r is an integer polynomial, shifting n by3·2^k preserves its residue lane and its value modulo2^k, giving the upper bound. Any period m modulo2^k is a period modulo2, so (5) forces6|m. In the lane n=2 mod3,

    P2(X)=−X+4Q(X),  Q in Z[X].

Hence

    P2(n+m)−P2(n)=m(−1+4R(n,m))

for an integer R. The factor in parentheses is odd, so its2-adic valuation is exactly v2(m). A period must therefore satisfy2^k|m as well as3|m. This proves (9).

The shift has the same exact order on M_(2^k): it fixes every shifted generator after P_k steps, while the original generator itself has no smaller period. Its compatible orbit system therefore has limit

    lim_k Z/(3·2^k) ≅ Z/3 × Z_2,                       (10)

and supernatural clock order3·2^infinity. The orbit closure is a torsor for this phase group, identified after choosing a as phase zero; its embedding in the additive state lattice is not a group homomorphism. Indeed the lattice has no3-torsion, and the orbit does not contain the zero sequence. The closure of the acting shift operators is the actual procyclic group. This is a precise completion of a compatible clock, not an isomorphism between an integer, a state graph and a physical path. Its2-adic component is genuine, but it does not have all prime exponents infinite.

The inverse limit of the entire observable state modules is much larger than this one orbit closure: it is H Z_2^22, a finite-index lattice in Z_2^22, abstractly a free Z_2-module of rank22. The finite index is2^206. The orbit completion and state-space completion are distinct objects.

There is also an explicit completed action. Let C be the coefficient companion matrix for multiplication by x in Z[x]/(χ), and let S be the first22-values shift matrix. The recurrence gives H C=S H exactly. Thus H:Z_2^22→HZ_2^22 is an equivariant module isomorphism. The extra finite-precision information is the invariant filtration

    I_k=H^(-1)(2^k Z_2^22),   M_(2^k)≅Z_2^22/I_k.

It is not merely the standard filtration2^kZ_2^22. Ignoring that distinction would incorrectly give22 binary coordinates at precision one.

## 6. Why a finite recurrence fit does not create the quotient

Every finite binary prefix has both a periodic extension and a non-eventually-periodic extension. A periodic extension repeats the prefix (with the empty prefix handled separately). For the second, append zeros and then place ones at factorially spaced positions. The arbitrarily long zero runs and infinitely many ones rule out eventual periodicity. A sequence over a finite field generated by a finite linear state system is eventually periodic; one generated by an invertible system is purely periodic. Therefore the second extension has no finite linear-state realization.

In particular, agreement with an LFSR on a finite sample never proves global shift stability of a proposed kernel. The plane blueprint proves rationality and kernel stability using its functional/cover argument; those are load-bearing conclusions, not consequences of fitting finitely many Cauchy values.

The positive consumer is (1): once a global monic recurrence is genuinely known, a finite r×r Hankel matrix does classify the whole observable module at every modulus. The distinction is proof of global closure versus measurement of a finite prefix.

## 7. Further precise questions

1. Classify the finite observed quotient actions at the bad primes, equivalently the invariant filtrations H^(-1)(p^k Z_p^22) under the known companion action. The completed abstract action itself is already the companion action.
2. Determine which source-free tournament readouts have stable integral observable lattices with a similarly sharp bad-prime filtration.
3. Relate the blueprint's functional quotient to an explicit nondegenerate pairing when its algebra becomes finite-dimensional; do not assume the whole original field has a finite quotient.

These questions extend recovered recurrence work and the supplied finite-quotient mechanism. They do not establish a new coloring theorem, an arithmetic dynamics conjecture, or the proposed sphaleron isomorphism.
