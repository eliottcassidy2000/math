# Fortunate offsets, 30013, and the fixed-line quadratic-prime problem

10 November 2026. **PROVED elementary identities and obstruction; FINITE-EXACT
controls; CITED literature; OPEN infinitude questions.** No new theorem ID,
literature priority, Fortune proof, Landau proof, or LRC consequence is claimed.

## Inheritance and working board

Closest mechanism: [Session51 prime-offset theorem](../session51_structural_errors_relational_packets_20261024/prime-offsets/THEOREM.md),
Sections 2–3: first-hit definitions, finite membership cutoff, shared two-sign
wheel, and prime-offset square threshold. The corrected near miss is
[THM-5110, primorial-crossing-gap-repair-and-brown-prime-power-inverse](../../../01-canon/theorems/THM-5110-primorial-crossing-gap-repair-and-brown-prime-power-inverse.md):
offset-one primality must be restored before claiming a consecutive-prime gap.
The canonical hostile is 30030+19=151*199. The useful sidecar here is the
actual Gaussian coordinate pair, rather than just its prime norm.

Portfolio: anchor = interpret the supplied lists and Landau connection;
niche = distinguish axis-near from diagonal-near Gaussian primes;
wildcard = explicit composites in every finite congruence class.
Board: first-prime offset; sorted support; local root screen; Gaussian line;
prime factor versus prime value. Each calculation below retains these types.
Method cards used: “Search the statement before the method”, “Type every
analogy and every implication”, and “Classify the target over the exceptional
stratum before closing by density”. No new meta-pattern is proposed.

## 1. Identify the lists without changing their indexing

Write p_k for the kth prime and P_k for their product. Define

    F_k = min{d>1: P_k+d is prime}, k>=1;
    L_k = min{d>1: P_k-d is prime}, k>=2.

The first supplied list is L_k for k>=3: 7,11,13,17,29,23,43,41,73,...;
the full list starts L_2=3 ([A055211](https://oeis.org/A055211)).
The second is the increasing distinct prime values of F_k
([A046066](https://oeis.org/A046066)); the indexed sequence instead begins
3,5,7,13,23,17,19,23,... ([A005235](https://oeis.org/A005235)).
In this context the third list is the odd unfortunate primes; including 2
gives [A045493](https://oeis.org/A045493). General primality of every F_k
or L_k is conjectural, not a consequence of their definitions.

Inherited exact control: P_6=30030, F_6=L_6=17, and 30013/30047 are prime.
But 30029=P_6-1 is prime and 30031=P_6+1=59*509. Thus 30013 is the last
prime strictly below P_6-1, not the nearest prime below P_6; the crossing
consecutive-prime gap is 30047-30029=18. These facts were already proved
in THM-5110 and are independently replayed here, not promoted as new.

## 2. What 30013 actually contributes to the quadratic picture

Exact identities:

    30013 = 122^2 + 123^2 = (245^2+1)/2;
    (245+i)/(1+i) = 123-122i.

Thus 30013 is the norm of a Gaussian prime whose coordinate magnitudes differ by 1.
It is not n^2+1: 173^2 < 30012 < 174^2. For all integers k>=0,

    k^2+(k+1)^2 = 2k^2+2k+1 = ((2k+1)^2+1)/2.

The map is multiplication by 1+i on Gaussian integers, or the displayed
substitution on norms. It preserves the sum-of-squares identity but
multiplies the norm by 2. Primality of the odd norm therefore becomes
“twice a prime”, not primality of n^2+1. The factor 2 and the parity of n
are indispensable. This is a genuine fixed-line analogy, not an implication
between the infinitude questions. Bunyakovsky would predict infinitely many
primes for both n^2+1 and 2k^2+2k+1; no unconditional infinitude is claimed.

For Landau itself the prime values begin 2,5,17,37,101,197,257,...
([A002496](https://oeis.org/A002496)). For positive n, apart from n=1 giving 2,
n must be even.

## 3. Shared obstruction, different local screens

Inherited Session51 lemma: if P_k +/- d is prime and greater than p_k,
then gcd(d,P_k)=1. Hence d>1 is prime whenever d<p_(k+1)^2. Every prime
factor of a composite d would otherwise be too large for their product.
The endpoint-size condition is automatic for F_k and for L_k with k>=3;
L_2=3 is the exceptional small endpoint. The strict threshold is necessary
without the first-hit condition: 30030+17^2=30319 is prime, although 289
is composite. This is not a Fortune counterexample because F_6=17.

For n^2+1 an odd prime divisor ell forces n to have multiplicative order 4
modulo ell, so ell=1 mod4. The forbidden-root count is rho(2)=1 and
rho(ell)=2 or 0 according as ell=1 or 3 mod4. At 2 the root is 1.
In contrast, a prime ell dividing P_k forbids just d=0 mod ell on both
primorial signs. The screens have different residue geometry.

**Elementary finite-congruence obstruction.** For any W>=1 and a>=1 put
A=a^2+1. For every integer t>=1 define n=a+tWA. Then

    n = a (mod W),
    n^2+1 = A * (1+2atW+t^2 W^2 A).

Both factors are greater than 1. Therefore every residue class of inputs modulo
any fixed W contains infinitely many n for which n^2+1 is composite. Choose a positive
representative (a=W for residue zero). If gcd(a^2+1,W)=1, all these composite
values still avoid every prime dividing W. This is the elementary polynomial
congruence mechanism; no novelty is asserted and no prime-density statement
follows. It does not rule out sophisticated sieves with growing cutoffs or
additional analytic information.

Positive/hostile controls:

    6^2+1=37;
    30^2+1=901=17*53;
    30030^2+1=901800901=73*12353437.

Both factors in the last line are prime. Every prime through 13 is avoided.
The two-prime example therefore obeys the local restrictions without being
a Landau prime. Connection contract: map integers to their residues modulo
W; preserve divisibility by primes dividing W; lose larger factors and their
multiplicities; restore these with a true primality certificate or a global
prime-counting argument. The displayed factorization is the decisive test.

## 4. Literature audit and logical boundaries

These are **CITED** results, not independently re-proved here.

- [Iwaniec (1978)](https://gdz.sub.uni-goettingen.de/download/pdf/PPN356556735_0047/LOG_0016.pdf):
  infinitely many n^2+1 have at most two prime factors, counted with
  multiplicity. This does not isolate infinitely many prime values.
- [Harman–Lewis, Theorem 1](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/4541183C02D39FD651979931A1B3B302/S0025579300014388a.pdf/gaussian-primes-in-narrow-sectors.pdf):
  the y=O(r^0.119) statement is supported. Their introduction records
  Ankeny–Kubilius's y=O(log r) under RH for the relevant Hecke L-functions.
  This audit does not certify an exhaustive present-day record claim.
- [Green–Sawhney, Theorem 1.1, v2](https://arxiv.org/html/2410.04189v2):
  for each fixed POSITIVE d=0 or 4 mod6, infinitely many primes r=x^2+d y^2
  have x,y both prime. Positivity excludes d=0 itself. The general-form
  extension is an author remark, not an established theorem for arbitrary
  forms without admissibility conditions. Setting y=1 is not permitted.
- [Grimmelt–Merikoski, Theorem 1.1, v2](https://arxiv.org/html/2505.00493v2):
  for a=h=1 and all sufficiently large X, some n in [X,2X] has
  P^+(n^2+1)>X^1.312. In terms of that n the literal deduction is
  P^+(n^2+1)>2^(-1.312)n^1.312; any exponent strictly smaller than 1.312
  absorbs the constant eventually. Do not silently identify X with n.
  [Li v2](https://arxiv.org/html/2406.07575v2) states 1.317; this is a
  **PREPRINT CLAIM, not independently audited here**, so 1.312 should not
  be advertised as an exhaustive current record.
- [Friedlander–Iwaniec](https://arxiv.org/abs/math/9811185): infinitely many
  primes x^2+y^4, with both variables allowed to vary. There is no deduction
  for the fixed slice y=1.
- [Baier–Zhao, Theorem 5′ and Section 9](https://www.impan.pl/shop/en/publication/transaction/download/product/82297):
  for each epsilon>0, infinitely many primes r=a m^2+1 with
  a<=r^(5/9+epsilon); a varies. The epsilon exponent requires their
  square-modulus Elliott–Halberstam type hypothesis. The authors'
  [preprint](https://citeseerx.ist.psu.edu/document?doi=f1df38d28b9929b7fb3d7e4909ba4d5de11d4675&repid=rep1&type=pdf)
  explicitly states the GRH 1/2+epsilon variant. Neither fixes a=1.

Elementary boundary: for n>=2, P^+(n^2+1)>=n^2 implies n^2+1 is prime,
since a prime factor of a composite N is at most N/2<n^2. This is a
statement in the actual variable n; a dyadic X^2 bound with n in [X,2X]
is not the same claim. Merely approaching exponent 2 does not remove the
cofactor. Likewise, y=O(log r) does not imply y=1. Fermat numbers of index
k>=1 have the form (2^(2^(k-1)))^2+1; the Fermat prime 3 is an exception.

The supplied sieve upper bound, if written A(x)=O(sqrt(x)/log x) for prime
values at most x, yields A(x)/floor(sqrt(x-1))=O(1/log x). Thus almost all
quadratic values are composite. Density zero is compatible with infinitely
many prime values. The missing result remains a suitable lower bound.

## 5. Reproduction and close-out

Run from the repository root:

    python3 05-knowledge/results/fortunate_landau_bridge_20261110/check.py

The standard-library verifier uses independent trial division and an
Eratosthenes sieve on 0..30319; it checks all offsets 2..289 at P_6,
the quadratic sequence for 1<=n<=174, all prime root counts through 97,
the displayed factorizations and Gaussian identities, and the congruence
construction for 1<=W<=30, 1<=a<=W, 1<=t<=5. It does not test an infinite
conjecture. `checks.json` records the deterministic results and script hash.
The general obstruction is proved algebraically above; the finite controls
are independent error detectors. Independent reviews accepted the algebra,
inherited scopes, and literature section; normal and optimized replays both
matched the saved JSON. Prepublication wording was corrected to distinguish
coordinate difference from Euclidean distance, and input residue classes from
output residue classes. Landau and the Fortune conjectures remain
OPEN. No proposed implication between them survives as a theorem.
