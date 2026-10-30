# Exact conditional source energy at every dyadic depth

**PROVED scoped**, same-agent proof review. E1--E4 are local labels in this
file. These are source-character identities for finite words and conditioned
Haar measures. Collatz convergence, fixed-integer H1 and the all-layer5/8
content bound remain OPEN.

## Definitions and inherited inputs

Fix an odd integer q>=3. For a positive valuation word w=(a_1,...,a_m),
let A=sum a_i, P=q^m, Q=2^A and

    B_w=sum_(i=1..m) q^(m-i) 2^(a_1+...+a_(i-1)),
    r_w=-B_w P^-1 mod 2^A,
    chi_w(R)=exp(2 pi i r_w R/2^A),
    V_m(R)=sum_(length w=m) 2^-A chi_w(R).

Each V_m converges absolutely and uniformly because its coefficient weights
sum to1. The source-address partition, full-Haar orthogonality and odd-Haar
nearest-neighbour energy are inherited from
[source_address_two_clocks_20261004.md](../source_address_two_clocks_20261004.md).
In particular all compositions of cost A biject with odd residues modulo2^A.
These are prefix addresses: the final valuation is allowed to continue beyond
the horizon. The extra bit certifying an odd endpoint is not included.

Fix h>=1 and an odd c modulo2^h. Write E_c for normalized Haar expectation
on the SAME original-source cylinder R=c mod2^h. No independent source reset
is used below. Put zeta=exp(2 pi i/2^h).

For any parity string beta of length h, including the all-zero string, let
r_beta be the unique source residue modulo2^h giving that shortcut parity
string under T_q(n)=(qn+1)/2 for odd n and n/2 for even n. Let |beta| be
its number of ones. Define the finite polynomial

    P_(h,u)(z)=sum_beta z^|beta| zeta^(u r_beta)
              =sum_(ell=0)^h A_ell(u) z^ell.              (1)

Here u is an odd residue. All 2^h strings occur once. The coefficients lie
in Q(zeta), and r_beta is computable by inverse parity or forward replay.
Let the short-cost part be

    F_m(c)=sum_(length w=m, A<=h) 2^-A zeta^(c r_w 2^(h-A)),

which vanishes for m>h.

## E1. A uniform exact energy factorization

For every finitely supported complex sequence z_m (m>=1),

 E_c |sum_m z_m V_m|^2
  = |sum_(m=1)^h z_m F_m(c)|^2
    +4^-h sum_(j>=1) 3^-j
       |sum_(ell=0)^h A_ell(c q^-j) z_(j+ell)|^2.          (2)

All inverses in (1)--(2) are modulo2^h. Only finitely many j contribute
for a finitely supported z. This classifies the entire conditional Gram
matrix for every q,h,c without an unbounded word enumeration.

Proof. On the cylinder c mod2^h, two characters correlate only when their
frequency difference has denominator dividing2^h. Characters of different
costs A>h cannot do this: their common-denominator numerator is odd. Nor
can such a character correlate with one of cost at most h. Thus the short
part contributes its deterministic square, and each cost A>h contributes
separately.

Write A=L+h, L>=1. At that cost two source residues correlate exactly when
they agree modulo2^L, equivalently when the first L shortcut parities agree.
Let that common prefix have j ones; its first parity is1. Its affine map
has multiplier q^j/2^L. Varying the suffix residue s changes the source
residue by 2^L q^-j s modulo2^(L+h), up to a common offset. The offset
cancels in the squared modulus. The suffix contribution is therefore the
square in (2). There are binomial(L-1,j-1) such prefixes, each with squared
weight4^-(L+h). Summing over L gives

    sum_(L>=j) binomial(L-1,j-1)4^-L=3^-j.

This proves (2). One may first truncate costs at T>=h; replace3^-j by
the finite sum over j<=L<=T-h. The finite identity is purely algebraic.
Uniform convergence of the finitely many V_m then justifies the limit.
No switch between odd-step and halving clocks at a fixed seed is assumed.

## E2. Sharp memory growth and exact source recovery

Let G_c(m,n)=E_c[V_m conjugate(V_n)]. Then

    G_c(m,n)=0 if |m-n|>h.                               (3)

The bandwidth h is sharp for every q,h,c, and for every j>=1:

    G_c(j,j+h)=4^-h 3^-j zeta^(-c q^-j r_ones),           (4)

where r_ones is the odd source residue for the all-one parity suffix of
length h. Explicitly r_ones=-(q-2)^-1 mod2^h: this is the fixed point of
the all-one shortcut branch, and inverse parity gives the same residue.
In particular (4) is nonzero. Since r_ones and q are odd units,
its normalized phase determines c modulo2^h uniquely. Thus the2^(h-1)
odd-source cylinders have2^(h-1) distinct conditional Gram matrices.

Proof. Each square in (2) couples only j through j+h. The short block
couples only1 through h. A pair j,j+h has only one contributing square,
with coefficients A_0=1 and A_h=zeta^(c q^-j r_ones). This gives (4).
The invertible exponent proves source recovery. For q=3, r_ones=-1 mod2^h,
so the normalized farthest correlation is exactly zeta^(c3^-j).

This is a precise alignment-memory law. It is a theorem about the ordered
Gram observer, not a lower bound for every possible nonlinear program or a
claim that h alone reconstructs a positive integer.

## E3. Periodic exact remainders and the eighth-root formula

Because beta runs over all residues, P_(h,u)(1)=sum_s zeta^(us)=0.
Define B_b(u)=sum_(ell=0)^b A_ell(u), 0<=b<h. Then

    P_(h,u)(z)=(1-z) sum_(b=0)^(h-1) B_b(u) z^b.          (5)

Thus (2) is a finite-difference energy with transported coefficients. Let
S_M=sum_(m=1)^M V_m. For every M>=h,

 E_c |S_M+1/2|^2
   =3^-M 4^-h sum_(b=0)^(h-1) 3^b
                   |B_b(c q^(-(M-b)))|^2.                (6)

The coefficient of3^-M is periodic in M with period dividing
ord_(2^h)(q). Equation (6) is an ordinary finite-M identity at each fixed
source precision, not an empirical asymptotic fit.

Proof. The short-cost sum is -1/2 on an odd source: the cost-one term is
-1/2, and at every cost2 through h the complete sum of odd characters is
zero. Set z_m=1 for1<=m<=M and0 afterward in (2), and center by+1/2.
Every window fully within1..M vanishes by (5). The only remaining windows
have j=M-b, 0<=b<h, proving (6). The phase period is immediate.

At q=3,h=3 let zeta=exp(pi i/4). Direct inverse parity gives

    P_(3,1)(z)
      =1+(-1-zeta+zeta^2)z+(zeta-zeta^2+zeta^3)z^2-zeta^3 z^3
      =(1-z)[1+(-zeta+zeta^2)z+zeta^3 z^2].               (7)

For the other odd u, replace zeta by zeta^u. Define the real character

    chi_8(c)=1 for c=1,7 mod8, and -1 for c=3,5 mod8.

It is multiplicative on the four units, as its multiplication table verifies.
Since zeta-zeta^3=sqrt(2), equation (6) becomes, for EVERY M>=3,

 E_[R=c mod8] |S_M+1/2|^2
   =3^-M [1/4+(-1)^M chi_8(c) 3 sqrt(2)/64].              (8)

To verify the reduction directly, the three partial sums in (7) are
1, -zeta+zeta^2, zeta^3. Their squared norms are1,
2-chi_8(u)sqrt(2),1. Insert u=c3^(-(M-1)) in the middle term of (6),
use chi_8(3)=-1, and obtain (8). Averaging the four source classes cancels
the character term and recovers the inherited odd-Haar remainder3^-M/4.

The square-root term is an actual real trace of eighth roots, not a numerical
analogy. Its sign records both the source class and the retained odd-step
clock. Forgetting either loses the exact conditional remainder.

For general h, conditioning the inherited odd-Haar remainder also gives
the coarse uniform bound E_c|S_M+1/2|^2<=2^(h-3)3^-M. Fixed-depth L2
convergence follows; a specified integer is still a Haar-null point. The
bound is not uniform as h increases, and no fixed-seed H1 statement follows.

## E4. Complete finite-depth product-transport classification

**Concurrent precedence:** Session115 published this same multiplier
classification at65cec73488 while this packet was being developed:
[WORD_UNITARY_CLASSIFICATION.md](../session115_odd_eighth_classifications_20261030/unitary/WORD_UNITARY_CLASSIFICATION.md).
It also proves the stronger equivalences with real preservation, positivity,
and permutation pullback, both on all residues and the odd coset. The
multiplier-only theorem below is therefore a credited independent derivation
and replay, not a distinct advance. Our stated variant additionally allows
arbitrary odd offsets b,c. E1--E3 are separate conditional-energy results.

Consider two odd multipliers q,p>=3 and odd offsets b,c, including -1.
For the map with odd branch (qn+b)/2, the source address is b r_w(q)
modulo2^A. At a fixed cutoff D, match the SAME positive word between the
two systems at every cost A<=D. Together with the constant character, these
words index all characters of Z/2^D. This defines a character permutation U.

Then U preserves all pointwise products of characters if and only if

    p=q mod 2^(D-1).                                     (9)

In the compatible case it is multiplication of every frequency by the unit

    u=c q/(b p) mod2^D.                                  (10)

For p!=q, the FIRST failure is exactly

    D=v_2(p-q)+2.                                        (11)

At that first depth an explicit product relation changes by the character
1/2, so the two outputs are negatives on every odd seed.

Proof. A product-preserving permutation of the cyclic character group must
multiply frequencies by a unit. The one-letter word (D) has frequency
-b/q modulo2^D, a generator, so that unit is forced to be (10). For the
two-letter word w=(1,D-1), D>=2, its target frequency minus its predicted
scaled frequency is

    [-c(p+2)/p^2] - [-(c q/p)(q+2)/q^2]
       =2c(p-q)/(p^2 q) mod2^D.                          (12)

This proves necessity of (9). For sufficiency, expand every word frequency
as -b sum_i 2^(a_1+...+a_(i-1))q^-i. The i=1 term transports exactly.
For i>=2, the difference after scaling contains p-q and at least one power
of2. Under (9) every difference therefore vanishes modulo2^D, including
the embeddings of lower-cost characters. D=1 is immediate.

At the depth (11), let k=1+2q^-1 mod2^D. The source frequencies obey
frequency(w)=k frequency((D)). Equation (12) has valuation D-1, so its
image differs from k times the image of (D) by2^(D-1). This is the explicit
antipodal product defect. For q=3,p=5 it occurs at D=3, the odd eighths;
for q=3,p=7 it first occurs at D=4. Distinct odd multipliers therefore have
no product-compatible all-depth word-matching transport.

The plus/minus change b=1,c=-1 at fixed q has no such defect: it is exactly
frequency negation at every depth. Its source reflection exchanges positive
and negative integers, so it does not identify the positive basins or cycle
counts of3n+1 and3n-1. This is the preserved sign map and its precise boundary.
