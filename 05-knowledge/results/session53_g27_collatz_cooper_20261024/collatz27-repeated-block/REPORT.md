# Repeating 27's first-descent block: binary fuel and a C2 residue quotient

Status: PROVED elementary specialization of inherited affine-word and periodic-point mechanisms; FINITE-EXACT certificate replay. No novelty claim for the general mechanisms. No Collatz convergence, coverage, integer cycle, total-stopping formula, or universal 153N bound is claimed.

## 1. Exact object and inherited credit

Let T(n)=n/2 for even n and (3n+1)/2 for odd n. Let w be the first 59 shortcut parity bits of 27:

11011111010110111011110100111011011111100111100010101000100

The fixed base certificate, already established in the first packet, has 37 odd steps, endpoint 23, and no earlier value below its source. Put

q=2^59=576460752303423488,
a=3^37=450283905890997363,
b=23q-27a=1100931843921811423,
d=q-a=126176846412426125.

On the guarded class n=27 mod q, the native block is F(n)=(an+b)/q. Every positive member of that class has exact first T descent 59 and exact first ordinary descent 96. The base packet proves this by retaining the original source and checking every intermediate multiplier 3^s_j>2^j for j<59.

Prior mechanism search found:
- collatz_mediant_tree_K_20260930.md, Theorems M1/M2: prefix rational fixed points, observer identity, and descent threshold;
- collatz_coherence_20260930_trichotomy_hexagon_motifs.md: periodic parity words and rational 2-adic periodic points;
- collatz_guarded_pumping_memory_20261004.md: exact prime-adic fuel for repeated words, including the dual binary fuel formula;
- collatz_boundary_fibres_20261011.md, section 3: repeated-word multiplier powers, and why deleting the carry does not authenticate a native path.

These mechanisms are inherited and independently specialized below. The new packet does not re-credit the first packet's word or the general pumping principle as a new theorem. Its incremental content is the explicit 27-block repeated cylinder, exact fuel and mod-7 quotient certificate, with a native-versus-quotient hostile. The first packet remains unchanged.

## 2. All finite repeated cylinders

For k>=1 put Q_k=q^k, A_k=a^k, S_k=(q^k-a^k)/(q-a), and B_k=bS_k. S_k is an integer by the geometric-series identity. Define the least residue

r_k = -B_k A_k^(-1) mod Q_k,  0<=r_k<Q_k,
e_k = (A_k r_k+B_k)/Q_k.

Then the following statements hold for every k>=1 and every integer t>=0:

1. n=r_k+Q_k t is positive and has exactly the prescribed first 59k parity bits w^k.
2. T^(59k)(n)=e_k+A_k t=F^k(n).
3. At each of the k block boundaries, the next block gives first descent below that block's own source exactly 59 shortcut steps later. The corresponding ordinary prefix contains exactly 96k steps, with a first descent at ordinary relative time 96 inside each block.

Proof. Composition gives F^k(n)=(A_k n+B_k)/Q_k. Each parity word of length 59k is realized by exactly one class modulo Q_k; its affine coefficient A_k is odd. The integrality congruence therefore identifies that unique class, and proves statement 1 as well as statement 2. Because w^k begins with w, r_k=27 mod q, so r_k>=27. Every nonterminal block source has parity prefix w and is a positive integer in that same class. Apply the base first-descent result at each such source separately. Each block contains 37 odd shortcut operations and hence 59+37=96 ordinary operations. No comparison threshold is silently reset in the claim: each fresh descent is explicitly relative to its own block source.

These are prefix clocks, not total stopping times. In particular, 27 itself only supports one full copy of this word: its next value 23 is not 27 mod q.

First examples:
- k=1: r_1=27, e_1=23.
- k=2: r_2=209923714958077886127437625648218139, e_2=128083994787268184645564261786926256.

The k=2 source and endpoint are equal modulo 7 but not equal as integers; the endpoint is strictly smaller. Further exact examples k=3,4 are in normal.json.

## 3. The exact finite binary-information budget

The rational fixed point is

x*=b/d=1100931843921811423/126176846412426125.

The displayed fraction is reduced, its denominator is odd, and 8<x*<9. In particular it is not an integer. Since d is odd, the residue formula simplifies to

r_k = b d^(-1) mod q^k,
C_k = {positive n: d n-b=0 mod 2^(59k)}.

Indeed modulo q^k, S_k=(q^k-a^k)/d=-a^k/d, and substitution gives the displayed identity. Consequently C_(k+1) is contained in C_k.

For every positive integer n, define

K(n)=floor(v2(dn-b)/59).

Then K(n) is exactly the maximum number of full consecutive initial copies of w in n's T parity sequence. This is a finite integer: dn-b cannot be zero, since x* is not an integer. The proof is the equivalence n in C_k iff v2(dn-b)>=59k, including k=0.

On a source with K(n)>=1, the native block obeys

dF(n)-b=(a/q)(dn-b),
v2(dF(n)-b)=v2(dn-b)-59,
K(F(n))=K(n)-1.

Thus the rank measures the remaining copies of this specific word exactly. At K=0 the repeated-word certificate ends; it supplies no law for subsequent words. It is not a globally decreasing Collatz rank.

The representatives r_k converge 2-adically to x*, and the intersection of all cylinders in Z_2 is the singleton {x*}. Their intersection among positive integers is empty. This does not deny that x* is a periodic point for the rational/2-adic extension of T: it has w repeated forever and T^59(x*)=x*. It denies an integer cycle with this itinerary.

The two geometries must not be confused: F has real multiplier a/q in (0,1), but its 2-adic multiplier has norm 2^59. The forward block expands 2-adic errors; its inverse contracts them. The nested source cylinders, not a forward 2-adic contraction, produce x*.

## 4. A genuine C2 quotient, not a native involution

Modulo 7 the base guarded block induces f(x)=1-x. Therefore

F^k(n) = n mod7 if k is even,
F^k(n) = 1-n mod7 if k is odd,

for sources in C_k. At the level of affine residue maps, k -> f^k is a semigroup homomorphism onto the order-two group C2. Equivalently, the infinite cyclic group generated by the invertible rational affine map F maps onto C2 by reducing its coefficients modulo 7.

But the real affine iterate is

F^k(n)=x*+(a/q)^k(n-x*).

On every admissible positive source, n>=27>x*, so F^k(n)<n for k>=1. Also F^2 is not the identity affine map because its slope (a/q)^2 is strictly between 0 and 1. Restoration of the residue after two copies discards the magnitude and is not physical return. The common-word dyadic guard remains indispensable.

Connection ledger: source is the guarded rational/integer affine block; target is the two-element subgroup {id,1-x} of affine permutations of F_7; map is coefficient reduction modulo 7; preserved predicate is endpoint residue transport and the corresponding Paley orientation sign; lost data are magnitude, carry beyond mod 7, positivity, and legality of future copies; restoration sidecar is the exact binary fuel K(n). The k=2 example is the decisive hostile against upgrading the quotient involution to a native one.

## 5. Density and quantifier boundary

C_k is one residue class modulo q^k, hence has natural density q^(-k)=2^(-59k) among all positive integers (and twice that relative density among odd integers). Any complete interval of q^k consecutive integers contains exactly one member. The family is nonempty at every finite depth, but its density tends to zero.

Because these cylinders are nested, their union for k>=1 is merely C_1, of density 2^-59. Arbitrarily deep construction proves for every k there exists a source, not that one integer realizes infinitely many copies, and not that every integer ever enters the family. The fuel gives no estimate for the remainder after the prescribed prefix, no total stopping formula, and no universal bound such as 153N. No Hamiltonian square-sum path count is inferred from this residue quotient.

### Exact distribution of the specific-word fuel

Write C_0 for all positive integers. For every integer k>=0,

{n: K(n)=k}=C_k minus C_(k+1),
density(K=k)=(1-1/q)q^(-k).

For every fixed residue z mod 7, CRT gives

density(K=k and n=z mod7)=(1/7)(1-1/q)q^(-k).

These are natural densities of deterministic integer sets. They can be called probabilities only under the explicit model of a uniformly sampled integer from 1..X followed by X tending to infinity. They describe initial repetitions of this one fixed word, not a law of unrestricted Collatz stopping times or independence of arbitrary orbit statistics.

Exact finite counts need no asymptotic assumption. For a class r mod M with 0<=r<M, let

A(X;r,M)=floor((X-r)/M)-floor(-r/M).

This counts positive n<=X in that class for every nonnegative integer X, including r=0. Use r_0=0 and modulus q^0=1. Then

#{1<=n<=X: K(n)=k}=A(X;r_k,q^k)-A(X;r_(k+1),q^(k+1)).

For the joint count, replace each residue r_j by its unique CRT lift R_(j,z) modulo 7q^j satisfying R_(j,z)=r_j mod q^j and R_(j,z)=z mod7. Each progression count differs from X divided by its modulus by at most 1; therefore each difference count has error at most 2 relative to X times its displayed density. The constants are uniform in X, k, and z. In particular fixed-character classes do not change the limiting geometric fuel distribution; no prime-restricted distribution is asserted.

## 6. Portable checks and controls

Run:
python3 check.py > normal.json
python3 -O check.py > optimized.json
cmp normal.json optimized.json

Normal and optimized outputs match byte for byte. The standard-library runner uses explicit raising checks, no assertions. Closed bounded universe: depths 1..32, eight lifts per depth (t=0..6 and 2^80+1), and every source 1..2048 for the exact fuel formula. Exact distribution and joint CRT counts are independently compared with enumeration at X=0,26,27,28,100,2048, levels k=0..3, and all seven residues. For each depth it computes r_k both by the closed inverse formula and by independent one-bit-at-a-time parity lifting. Every lifted source is replayed with direct ordinary operations, checking the word, fresh first descent, clock, endpoint, residue quotient, and exact final fuel exhaustion. Hostiles: 27 fails the second copy; the two-copy endpoint returns only modulo 7; x* is nonintegral; each tested fuel boundary rejects the next whole copy.

The proof covers every depth and lift. The finite replay is a certificate sanity check, not the source of the universal quantifier. The closed manifest includes exact inherited snapshots and the frozen first-packet dependencies. No repository files were edited.
