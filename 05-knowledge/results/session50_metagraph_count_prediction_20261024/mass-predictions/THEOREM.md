# Odd-fibre moment reconstruction of the tournament class count

Status: PROVED finite algebraic reconstruction theorem and exact rational bounds; FINITE-EXACT held-out atlas tests through order eight. This is aggregate-statistic recovery, not a closed formula in n, a faster enumeration algorithm, or an all-order asymptotic. No literature-priority claim is made for factorial divisibility, Newton interpolation, or moment-matrix inequalities.

## 1. Inheritance and the genuinely missing coordinate

THM-781, `01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md`, identifies the fixed-Hamiltonian-path word fibre of a tournament class as h=H/|Aut|. Both factors are odd, and the free action of Aut on Hamiltonian words makes h a positive odd integer. Session48 `observer-moments/THEOREM.md` independently authenticates every order-eight merged fibre and proves the observer-weighted class law. These are inherited facts, not discoveries of this packet. Frozen text snapshots and source commit/path/SHA256 are in `sources/SOURCE_PINS.json`.

For n>=2 let m=(n-1)(n-2)/2 and N=2^m. Let the unmerged isomorphism classes have fibre masses h_1,...,h_A. Then

    h_i positive odd,       sum_i h_i=N,
    K_r = sum_words h(word)^r = sum_i h_i^(r+1),    r>=0.       (1)

The unknown A is A000568(n). The known cube size is K_0=N. K_r is an ordered same-fibre collision count: it counts (r+1)-tuples of fixed-path words lying in one unmerged class. This interpretation does not require A as an input, but evaluating the collision relation can still be computationally difficult.

In the converse-merged atlas, let mu be the word-fibre mass. Self-converse fibres have odd mu; non-self-converse pairs have mu=2h, congruent to 2 modulo 4. Therefore, with epsilon=1 for odd mu and epsilon=2 for even mu,

    h=mu/epsilon,       K_r=sum_merged mu*(mu/epsilon)^r.       (2)

Equation (2) can be read as a port-mass statistic only when the port incidence convention really supplies mu; a loop contributes two incidences. A loop-deleting Laplacian diagonal is not interchangeable with mu. We do not prove a new port-degree identification here. No converse label is needed once the authenticated mass parity is retained. Using mu directly in place of h destroys the odd-mass hypothesis.

The live concept board is: odd fibre partition; collision moments; factorial fixed divisors; reciprocal-mass bounds; strong-component peeling; missing moment-generation algorithm. The hostile control is the even mass 2, which makes the first congruence false. The preserved target is the number of unmerged classes; the lost coordinate is the full fibre distribution. The sidecar is a short vector of aggregate moments, not a supplied class count.

## 2. The 2-adic finite-moment theorem

**Theorem 1 (odd-fibre reconstruction).** For any finite nonempty partition into positive odd masses, and any integer k>=1, set

    Q_k(x)=product_(j=0)^(k-1) (x-(2j+1))
          =q_0+q_1 x+...+q_k x^k,
    e_k=k+v_2(k!)=2k-popcount(k),       M_k=2^e_k.             (3)

Then q_0 is odd and

    A = -q_0^(-1) sum_(j=1)^k q_j K_(j-1)  (mod M_k).       (4)

The inverse in (4) is the ordinary modular inverse modulo M_k. If a certified interval [L,U] contains A and contains exactly one integer in the residue class (4), that integer is A. In particular, if M_k>=N, the universal interval 1<=A<=N determines A uniquely. A zero residue then represents M_k, not zero; if M_k>N a zero residue is incompatible with a genuine nonempty partition.

**Proof.** Write h=2t+1 with t>=0. Then

    Q_k(h)=2^k t(t-1)...(t-k+1)=2^k k! binom(t,k).           (5)

For t<k both sides vanish. Otherwise the binomial coefficient is integral. Thus M_k divides Q_k(h) for every allowed h. Summing this divisibility over classes gives

    q_0 A+sum_(j=1)^k q_j sum_i h_i^j =0 (mod M_k),

which is (4) by (1). The constant term is a signed product of odd integers, so it is invertible. Distinct integers in the same residue class differ by at least M_k; an interval of N consecutive integers has at most one when M_k>=N. The true count supplies existence. QED.

Equivalently, the integer polynomial

    P_k(x)=-q_0^(-1) sum_(j=1)^k q_j x^(j-1) (mod M_k)

satisfies h P_k(h)=1 modulo M_k for every odd h. Summing P_k(h(word)) over words recovers A modulo M_k. This is a finite 2-adic reciprocal-mass formula rather than a Burnside sum over permutations.

**Corollary 1 (all orders, bounded-size aggregate sidecar).** Set

    k(n)=min{k>=1: 2k-popcount(k)>=m}.

The residues of K_0,...,K_(k(n)-1) modulo M_(k(n)), together with N, determine A000568(n). Here k(n)=m/2+O(log(m+2)); the supplied modular sidecar has O(m^2)=O(n^4) bits. This is a size bound on a certificate, not a complexity bound for generating its inputs. The cases n=1 and n=2 both have a singleton word and count one, handled directly.

For n=8, m=21 and the least modulus-only choice is k=12, e_k=22. `twelve_residues.json` contains just twelve residues modulo 4,194,304. With N=2,097,152, `predictor.py` recovers 6,880 without any rational moment inequality, support bound, or atlas access.

## 3. Sharpness of the polynomial divisibility exponent

**Theorem 2 (sharp within the declared certificate family).** Let F be an integer-coefficient polynomial of degree at most k with F(0) odd. If 2^b divides F(h) for every positive odd h, then b<=e_k. Thus Q_k achieves the largest possible power of two for this unit-constant, degree-at-most-k polynomial certificate family.

**Proof.** Newton expansion with step two, evaluated at zero, gives

    F(0)=sum_(j=0)^k Delta_2^j F(1) binom(-1/2,j).           (6)

Every forward difference is an integer combination of F(1),F(3),...,F(2j+1), so is divisible by 2^b. Also

    binom(-1/2,j)=(-1)^j binom(2j,j)/4^j

has, in lowest terms, denominator exactly 2^(2j-popcount(j))=2^e_j. To verify the valuation, use v_2(j!)=j-popcount(j), hence v_2(binomial(2j,j))=popcount(j). The sequence e_j is increasing. If b>e_k, each term in (6) is an even integer, contradicting odd F(0). Equation (5) shows Q_k attains e_k, and Q_k(2k+1)=2^k k! has exactly that valuation. QED.

This is not an information-theoretic lower bound for every possible estimator, for the native tournament family, or for bounded-support moment optimization. Native structure or a better interval can determine A with fewer moments than the universal modulus-only threshold.

## 4. Positive-moment reciprocal bounds with no maximum-mass input

The inherited exact identity is

    A=sum_words 1/h(word).                                  (7)

Equation (7) by itself only relocates the unknown. The new use is to replace the reciprocal by certified finite-degree polynomials and combine the resulting interval with (4).

For any real polynomial p and any x>=a>0,

    1/x - (2p(x)-x p(x)^2)=(1-x p(x))^2/x >=0,              (8)
    [1-2(x-a)p(x)+x(x-a)p(x)^2]/a - 1/x
                 =(x-a)(1-x p(x))^2/(a x) >=0.            (9)

Thus both sides are global on the whole half-line, not just on observed mass values. No authenticated upper endpoint is needed.

Let p have degree at most r-1. For the lower bound define

    b_i=K_i,       G_ij=K_(i+j+1),       0<=i,j<r.

When G is positive definite, p's coefficient vector G^(-1)b optimizes the quadratic expression (8), giving

    L_r=b^T G^(-1)b <= A.                                  (10)

It uses K_0,...,K_(2r-1). For the upper bound define

    c_i=K_(i+1)-a K_i,
    J_ij=K_(i+j+2)-a K_(i+j+1),       0<=i,j<r.

When J is positive definite, p=J^(-1)c gives

    A <= U_r=(K_0-c^T J^(-1)c)/a.                          (11)

It uses K_0,...,K_(2r). The matrices are Gram matrices for the positive word measure, with weights x and x(x-a), respectively. Singular systems do not invalidate (8)-(9): any chosen polynomial remains a certificate. The implementation uses the largest available invertible leading block, falling back to 1<=A<=K_0/a; it never silently inverts a singular matrix. Increasing the polynomial space can only improve the corresponding optimized bound. The exact rational coefficient vectors are retained in the certificate.

**Peeling a known stratum.** The inherited rigid-strong theorem proves h=1 classes are exactly ordered words in strong atoms {1,C3}; H and Aut multiply over the unique strong-component decomposition. Therefore their count is c(n), where c(0)=1 and c(n)=c(n-1)+c(n-3), with c(negative)=0. This is the ordinary, not necessarily self-converse, rigid stratum. Its OGF is 1/(1-z-z^3), distinct from the palindromic self-converse subfamily. Subtract c(n) from every K_r, apply (8)-(11) to the residual with a=3, then add c(n) back. `sources/rigid_strong_theorem.txt` §4 states this exact all-rigid counting consequence. We do not guess the h=3 alphabet.

## 5. The held-out order-eight certificate

The original Session48 atlas supplies these eight aggregate word moments, with no class count or maximum-mass field in the prediction input:

    K0 = 2097152
    K1 = 810198400
    K2 = 351545786720
    K3 = 164084691719200
    K4 = 80658618203874752
    K5 = 41239651789795691200
    K6 = 21752734067701334577440
    K7 = 11769806604752582871641440.

For k=8, Q_8 and its odd constant term give

    A = 6880 (mod 32768).

The degree-seven lower polynomial (r=4) and degree-six upper polynomial (r=3,a=1) give exact rational bounds, shown here only in decimal for readability:

    6698.686237470946 <= A <= 26344.84826269345.

The only integer in this interval congruent to 6880 modulo 32768 is 6880. `blind_prediction.json` stores both full rational endpoints, both coefficient vectors, Q_8's integer coefficients, the modular inverse, and the residue. The executable prediction interface receives only `eight_moments.json`, a list of eight integers. It cannot read the atlas, a class total, support values, H, Aut, or a supplied target. The normal checker subsequently compares the result with the held-out count from the authenticated atlas. The same moments alone determine the interval and residue even for a non-tournament odd-fibre partition.

### A shorter hand-checkable upper certificate

The optimized matrix coefficients are not necessary for the eight-moment conclusion. Set

    p(x)=(x^2-1250x+465000)/50000000,
    P_U(x)=2p(x)-x p(x)^2+(1-x p(x))^2.

Then

    P_U(x)-1/x=(x-1)(1-x p(x))^2/x >=0     for every x>=1.

Using only K0,...,K6 to evaluate this degree-six polynomial gives exactly

    A <= 496641366973194539 / 15625000000000
       = 31785.04748628445... <32768.

Meanwhile the eight K residues modulo32768 are

    [0,9600,352,21536,9664,2752,4384,28000].

The coefficient list of Q8, from constant term upward, is

    [2027025,-4098240,2924172,-1038016,208054,-24640,1708,-64,1],

and 2027025^(-1)=17137 modulo32768. Substitution in (4) gives6880. Positivity, the strict upper bound below32768, and this residue alone force A=6880; no lower Gram bound is needed. `hand_certificate.json` retains the fully expanded rational upper polynomial and all these arithmetic steps. This is a particularly short exact proof from the supplied aggregates, not a fit to be extrapolated to other orders.

The n-only h=1 peeling c(8)=13 improves the interval to approximately [6711.2340967,13145.1288248] with the same eight moments. It does not reduce the first successful moment count in this particular interval-and-congruence hierarchy.

Finite tests from the inherited direct-inverse atlas at n=3..7, plus the Session48 order-eight atlas, give first successful moment counts:

    n:                          3   4   5   6   7   8
    without peeling:            1   2   3   4   6   8
    with n-only h=1 peeling:     1   1   3   4   6   8
    recovered A:                2   4  12  56 456 6880.

“First” here means the tested prefix of this specific polynomial/Gram-bound hierarchy, not an optimum among all mathematical methods. No smaller-n fitting or extrapolation is performed.

## 6. Hostile controls and what is still open

The abstract odd partitions {1:17,5:3} and {1:2,3:10} have the same N=32 and K1=92 but respectively 20 and 12 classes. Their difference eight is exactly compatible with the k=2 modulus eight. Thus low moments need not determine the count before the interval/residue isolation succeeds. These are odd-fibre partitions, not asserted native tournament metagraphs.

Replacing odd h by an unnormalized merged mass 2 gives one class but moment K0=2; the k=1 formula would predict an even count, correctly demonstrating the hypothesis failure. The checker also tests polynomial identities at even positive real-domain control points and masses well beyond the order-eight observed support; the analytic bounds do not rely on its empirical maximum.

What has been established is a finite, non-Burnside conversion from a bounded vector of intrinsic fibre/collision moments to the desired unmerged class count, with a sharp dyadic polynomial certificate and a working eight-moment test. What remains open is the important predictive step: derive or cheaply compute enough K_r from n, local native port operations, or smaller-order metagraphs without already performing a census of the same difficulty as A000568. The inherited first moments of score imbalance and cyclic triangles are not K_r and cannot be substituted for them. This packet makes no claim to have solved that remaining step.
