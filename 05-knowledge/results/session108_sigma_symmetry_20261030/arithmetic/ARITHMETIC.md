# Total divisor sums: arithmetic symmetries, prime ports, and update cocycles

Session 108, 2026-10-29. All integers are positive. Here σ(n)=sum_{d|n}d includes n; the aliquot sum is s(n)=σ(n)−n. These must not be interchanged. Statements labelled proved below have proofs here; the complete prime-power classification uses a cited classical theorem, not a finite scan or novelty claim.

## 1. The arithmetic map and its exact domain

Write R_p(a)=1+p+...+p^a, including R_p(0)=1. Unique prime factorization gives

    σ(n)=∏_p R_p(v_p(n)).

Let [x] denote a positive rational squareclass, represented additively by the parity vector of prime valuations. Define Q(n)=[σ(n)]. The target is the direct sum of one F_2 coordinate per prime. If gcd(m,n)=1 then Q(mn)=Q(m)+Q(n). Without that guard this is false: σ(3)^2=16 is square, but σ(9)=13 is not.

For arbitrary m,n the exact multiplicative defect is

    σ(mn)/(σ(m)σ(n)) = ∏_{p|gcd(m,n)} R_p(a_p+b_p)/(R_p(a_p)R_p(b_p)).

Thus define D(m,n)=Q(mn)+Q(m)+Q(n). This is supported computationally on shared **input prime ports**, though its output primes can be new. It obeys the cocycle identity

    D(l,m)+D(lm,n)=D(m,n)+D(l,mn),

by cancellation of each Q term. This is bookkeeping of a genuine arithmetic failure of multiplicativity, not an additional assumption about a binary code.

## 2. Exact exponent updates and a kernel theorem

**Proved.** If n=p^a u with gcd(p,u)=1, replacing exponent a by b gives

    σ(p^b u)/σ(p^a u)=R_p(b)/R_p(a),
    Q(p^b u)+Q(p^a u)=[R_p(b)]+[R_p(a)].

The right side is independent of u. An update universally preserves squareclass exactly when R_p(a)R_p(b) is an integer square. Updates a→b→c telescope to a→c. Updates at distinct input primes commute; updates at the same prime must compose as state replacements, not as independent attachments.

More generally, choose finite allowed exponent sets A_p, and a base exponent a_p^0 at each of finitely many ports p. Then

    Q(∏ p^{a_p}) = Q(∏ p^{a_p^0}) + Σ_p ([R_p(a_p)]+[R_p(a_p^0)]).

Consequently every globally neutral transition is exactly a cancellation among these local parity differences. This is an exact kernel description on a product of exponent-state sets. It becomes an ordinary F_2 linear kernel only for binary choices per port, where each coordinate is toggled between its two fixed states. Three choices at one port are not three independent bits.

For pairwise coprime attachment blocks b_i, the binary cube n(x)=u∏ b_i^{x_i}, x_i∈{0,1}, with gcd(u,∏b_i)=1, satisfies Q(n(x))=Q(u)+Σ x_i Q(b_i). Its neutral toggles form the linear kernel of the columns Q(b_i). Distinct input support is essential. Zero columns suffice but are not necessary: σ(5)=6 and σ(23)=24 share squareclass 6, so σ(115)=144; similarly σ(19)=20 and σ(79)=80 give σ(1501)=1600. These are cancellations of two nonzero columns.

## 3. Complete classification of zero prime-power columns (cited + elementary)

The primary research paper Michael A. Bennett and Aaron Levin, *The Nagell–Ljunggren equation via Runge's method*, Monatsh. Math. 177 (2015), 15–31, DOI 10.1007/s00605-015-0748-1, states in Proposition 1 that a solution of (x^n−1)/(x−1)=y^q outside its four listed solutions has odd q≥3. In particular its positive-x square solutions with n>2 are (x,n,y)=(3,5,11),(7,4,20). [Author-hosted paper, pp. 15–16](https://personal.math.ubc.ca/~bennett/BeLe-Mona-2015.pdf). This is the classical square case, not a solution of the general perfect-power problem.

Apply x=p, n=a+1. For a≥2 the only possibilities are (p,a)=(3,4),(7,3). For a=1, p+1=y² implies p=(y−1)(y+1), and primality forces y=2,p=3. Therefore

    σ(p^a) is square, p prime and a≥1 ⇔ (p,a)∈{(3,1),(3,4),(7,3)}.

Direct certificates: σ(3)=4=2², σ(81)=121=11², σ(343)=400=20².

The zero-squareclass exponent states are {0,1,4} at p=3, {0,3} at p=7, and {0} at every other prime. This classification says nothing by itself about nonzero local squareclass collisions or cancellation between distinct ports.

## 4. A genuine S_3 × C_2 arithmetic symmetry

Fix u coprime to 21 and consider

    Ω_u={3^a 7^b u : a∈{0,1,4}, b∈{0,3}}.

These six distinct integers all have Q(n)=Q(u). Arbitrary permutations of the three a states and the two b states act faithfully on Ω_u, producing an S_3×C_2 action preserving Q. These are state-replacement symmetries, not multiplication by one fixed integer and not preservation of σ's numerical value. The map is concrete on integers: replace v_3 or v_7 while retaining u. σ takes the multipliers {1,4,121}×{1,400}.

In particular 3 is a free *binary attachment* on numbers coprime to 3; 81 supplies a third alternative at that same port, not another independent bit. Their attempted product 3·81=243 has σ(243)=364, not square. The separate 343 port can be independently attached only when coprime to the base. On the restricted states a∈{0,1},b∈{0,3}, this specializes to a C_2² action.

A complementary local symmetry is replacement 19↔79 at its two-port binary cube: toggling both coordinates preserves Q, but requires the base coprime to 19·79 and each exponent in {0,1}. It does not say multiplication by 19·79 preserves Q on arbitrary integers.

### Observer boundary: raw quadratic characters retain zero data

The action preserves the global positive rational squareclass, not the prime support of σ(n). Consequently it need not preserve raw Legendre symbols at primes dividing an allowed square factor: σ(3)=4 and σ(81)=121 have the same squareclass, but χ_11(4)=+1 while χ_11(121)=0. Thus the 3↔81 state swap is not a symmetry of the raw χ_11 observer. To infer character preservation at a fixed prime q, restrict to the common unit domain where q divides none of the allowed σ(n), or retain zero/support data as an additional observer. On the six-state family, requiring q coprime to σ(u)·2·5·11 is sufficient for every allowed state to lie in that common unit domain. The zero value must not be silently replaced with +1 merely because its integer argument is a square.

### Square-sum aliquot edges and the six neutral products

An aliquot edge n→s(n) has endpoint sum n+s(n)=σ(n). The classification therefore says that among n>1 that are prime powers, exactly three edges have square endpoint sum:

    3→1 (sum 4), 81→40 (sum 121), 343→57 (sum 400).

This reconnects the exact aliquot segment 1001→343→57: σ(1001)=1344, σ(343)=400. It does not say that every edge of this segment has square endpoint sum (1344 is not square), or that the segment proves a general convergence claim.

Products assembled solely from individually neutral prime-power factors, with distinct prime ports, are exactly {1,3,81,343,1029,27783}. This is not the classification of all n with square σ(n): nonzero charges can cancel, for example σ(5·19·29)=6·20·30=3600.

A shared-port hostile sharpens the point: 3 and 81 both start in the zero squareclass, but multiplication by 3 sends them to 9 and 243, with squareclasses 13 and 91 respectively. Their quadratic characters modulo 5 differ: (13/5)=−1 whereas (91/5)=+1. Thus equal starting σ-squareclass is not enough information for a multiplication update that touches an occupied prime port. Retain the exponent state.

## 5. The requested small primes and hostile controls

For the squarefree support {3,5,19}, with bits (i,j,k),

    σ(3^i5^j19^k)=4^i6^j20^k,
    Q(3^i5^j19^k)=j[6]+k[5].

The columns [6] and [5] are independent, since the prime-3 coordinate sees j and the prime-5 coordinate sees k. Exactly the bit at 3 is neutral; kernel={000,100}. The eight sums in lexicographic (i,j,k) order are 1,20,6,120,4,80,24,480. None of 5 or 19 individually gives a zero column. A bit is a choice to include an input factor, not the number of divisors or a binary digit of the integer.

At exponent two σ(3²)=13, σ(5²)=31, σ(19²)=381=3·127. σ(5³)=156=4·39; σ(19³)=7240=4·1810. These demonstrate why a prime's first column does not predict later exponent updates.

Aliquot sums do not obey the same multiplicativity: for gcd(m,n)=1,

    s(mn)=s(m)s(n)+m s(n)+n s(m).

Example s(3)=1 and s(5)=1 but s(15)=9, not 1. Moreover s(3)=1 is square whereas s(81)=40 is not: the σ-neutral 3→81 replacement fails for s. No aliquot/Collatz convergence statement follows.

## 6. Verification boundary

`verify.py` uses only exact Python integer arithmetic and explicit checks (works under `-O`). It independently checks σ by direct divisor enumeration up to 2000; all pairs m,n≤100 for the shared-port formula; all triples l,m,n≤20 for the cocycle identity; bounded prime-power square searches p≤200,a≤16 by integer square root; all eight requested bits; the six-state action and hostile examples. The finite scan verifies examples and implementation, not the infinite classification, which depends on the cited theorem. No novelty claim or universal dynamics claim is made.
