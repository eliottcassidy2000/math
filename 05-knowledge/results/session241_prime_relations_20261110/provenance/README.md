# Session241 provenance: which arithmetic does the tournament actually retain?

**Status: source-bound synthesis and independent finite regression.** Chan–Pak
completeness, Carmichael's totient conjecture and the finiteness of Fermat primes
remain OPEN. No literature-priority claim. Source snapshot:
`8a6b203c9b33` (full SHA in `source_pins.json`). Historical packets missing from
the sparse working tree were read with `git show HEAD:<path>`, not treated as
missing mathematics. No old numerical census was rerun wholesale here.

## 1. Inheritance and concept board

- **Closest proved mechanism:** ordered strong components multiply Hamiltonian
  path counts, and the actual integer spectrum reduces to the prime and
  seven-times-prime lanes (THM-4094, full path in the pin ledger).
- **Canonical hostile:** 7 is absent, while its multiples 35,49,63,343 are
  attained. Thus ordinary divisibility of integers does not reflect a tournament
  join factorization. A strong tournament can have composite path count.
- **Corrected near miss:** the shared prefix 14,26,34,38 of nontotients and one
  Camion exponent list is not an equality; exponent50 separates them.
- **Least-used useful sidecar:** the full prime-supplier support of a totient
  fiber, or the full positive word state of a tournament recipe. Equal modular
  order alone drops both.

Anchor: arithmetic of the Hamiltonian-path value monoid. Niche: prime-power
unit actions and actual orbit stabilizers. Wildcard: factorization elasticity
as a precise measure of the gap between numerical and structural primality.
The five live objects are counts H(T), ordered strong components, unit actions,
prime-supplier fibers, and positive source words. Every connection below
specifies which of these is preserved; none is identified with all the others.

## 2. A hierarchy of genuine arithmetic maps

### A. Ordered join: exact multiplication, no arithmetic converse

For tournaments A,B, orient every cross edge from A to B. Every Hamiltonian
path traverses all of A and then all of B, giving a literal bijection
Ham(A=>B)=Ham(A)×Ham(B). Consequently H multiplies. The strong decomposition
of a fixed tournament is unique and ordered, but its list of numerical counts
is not prime factorization of H.

The attained-value monoid has identity1; even numbers never occur. Transitive
padding adds vertices without changing H. A prime H therefore forces exactly
one nontrivial strong component, but a strong component need not have prime
H. The exhibited strong values49,63,343 already show this.

THM-4094 proves that attaining every odd prime p!=7 and every 7p with p!=3,
together with the finite carries49,63,343, is equivalent to global coverage
of all odd values except7,21. It does not prove either infinite lane. Arbitrary
multiplication of already attained integers cannot manufacture a missing prime.

**Preserved:** multiplication and existence of ordered-join realizations.
**Lost:** individual graph structure, order, strong-class multiplicity and
endpoint response. **Sidecar:** actual component list and adjacency witnesses.

### B. Repunits: factorization complexity becomes an asymptotic degree

THM-4749 imports the fixed-value class-count polynomial: for a fixed attained
h, the number N(n,h) of n-vertex tournament isomorphism classes is eventually a
polynomial whose degree w(h) is the largest number of nontrivial strong
components possible at H=h. Its proof counts singleton padding around ordered
strong blocks. This is not the number of strong realizers of h.

The cyclotomic identities are useful because they constrain legal factorizations:
3=Phi_2(2) is a legal triangle count; 7=Phi_3(2) is not. In particular63 is
attained but cannot split nontrivially into attained values; its two strong
order8 realizers give N(n,63)=2(n-7). By contrast4095 can absorb7 into35.
Fermat factors F_i=2^(2^i)+1 are always attained by the reversed-extreme-arc
family, regardless of their primality. Thus the product identity for
2^(2^j)-1 supplies at least j strong factors without assuming Fermat primality.

Session219 already strengthens this to the exact smooth law:

    w(3^a)=a;
    w(3^a*7)=a-1 for a>=2;
    w(3^a*7^b)=a+floor(b/2) for b>=2.

The cases(a,b)=(0,1),(1,1) are the two omissions. In particular
w(63^k)=floor(5k/2) for k>=2, although w(63)=1.

**Preserved:** legal multiplicative decompositions and eventual polynomial
degree. **Lost:** the coefficients of N, individual graph topology and prime
support once only w is retained. **Sidecar:** all realizing strong classes and
their orders. Larger repunit tables in THM-4749 explicitly depend on certified
membership of their prime factors; do not upgrade those conditional rows.

### C. Prime unit actions: genuine Fermat and Wieferich detection

THM-4644 states, for an ODD PRIME p, that every circulant tournament on Z_p has
only translations as automorphisms iff p-1 is a power of2. Here “rigid” is
relative to the obligatory translation group, not a trivial automorphism group.
If an odd prime r divides p-1, the order-r unit orbits avoid their negatives,
so selecting one orbit from each opposite pair constructs an extra multiplier.
The converse at the full automorphism-group level uses the cited prime-degree
group theorem, not merely the multiplier calculation.

For d=ord_p(2), its odd part d', and g=2^(d/d'), the number of circulant
connection sets on Z_(p²) invariant under g is2^E, where

    E=(p-1)/d'                       if p is not base-2 Wieferich;
    E=(p+1)(p-1)/(2d')              if p is base-2 Wieferich.

This is a truly target-preserving arithmetic encoding: lift failure of the
unit order changes the number of paired orbits. It does not say that the
Hamiltonian-path count of any one such tournament is prime or a Fermat number.
In fact THM-4861 proves a vertex-transitive tournament has prime H only for C3.

**Preserved:** the specified unit action, stabilizer and order-lifting test.
**Lost:** any general correspondence between integer primes and H-values.
**Sidecar:** modulus, connection set, multiplier and distinction between
translations and the full automorphism group. **Hostile:** composite moduli;
the prime-degree hypothesis cannot be silently discarded.

### D. Aggregate labelled counts: real Korselt relations, not Chan–Pak

THM-4746 concerns s(n), the number of strong labelled tournaments of order n,
and c(n), the number of rotation-invariant tournaments on fixed Z_n. Neither is
H(T). For an odd prime p and m>=1,

    s(m+p-1)=(2/p)s(m) mod p,
    D(m+p-1)=(2/p)D(m) mod p, D=s-c.

For squarefree odd part, s(n)=c(n) mod n is equivalent to p dividing D(r_p)
for each odd p|n, where r_p is n modulo p-1 in1,...,p-1. D(1)=D(2)=D(3)=0
explains the Carmichael and3p families. The prime-power continuation is stronger
but has its own lower threshold; the original modulo-p proof does not imply it.
The first three free residues are structural data, not random numerology.

THM-4884 retains still more information: the labelled spectrum polynomial
P_n(x)=sum_T x^H(T), reduced coefficientwise modulo a prime p|n, equals the sum
over tournaments fixed by a specified order-p rotation. At n=p, all surviving
circulant exponents are multiples of p. This exact orbit statement constrains
multiplicities of values; it does not construct an individual graph at a target.

**Preserved:** orbit congruences and full CRT residues. **Lost:** exact counts
and graph witnesses when reducing modulo p. **Hostile:** 341 fails the Camion
test although it is a base-2 pseudoprime;1023 passes. No primality iff is claimed.

### E. Totient fibers: support, not superficial doubling

Session218 proves a genuine shared subfamily: if M=2^q-1 and
W=(2^q+1)/3 are both prime with q>5, then2q is a nontotient and2^(2q)-1
is a Camion pseudoprime. The first half uses the exact inverse fiber of2q:
for prime q>3 it is {2q+1,2(2q+1)} if2q+1 is prime, and empty otherwise.
The simultaneous M,W primality excludes that latter prime by divisibility.
This is a one-way conditional family, not equality of the two sequences and
not an infinitude result.

Session219 counts complete fibers A(m)=#{x:phi(x)=m}:

    sum_(a>=0) A(2^a)x^a = (2-x)/(1-x) product_(F_i prime)(1+x^(2^i)).

The missing F5 supplier forces A(2^a)=a+2 through31 and A(2^a)=32 for32<=a<=63.
Thus doubling can DECREASE fiber size33 to32. Its factors641 and6700417 have
odd factors5 and3*17449 respectively in their totients and cannot replace the
missing pure-dyadic supplier. This is the exact coordinate lost by looking
only at their common order ord(2)=64.

Session219 also provides literal tournament recipes which separate divisibility
by those two primes: all <=4-run O/Z recipes force6700417|H =>641|H, while
five runs suffice to escape. The result is an all-length finite phase theorem,
not a global obstruction and not a statement that equal orders force identical
behavior in arbitrary recipes.

**Preserved:** permitted prime suppliers and physical recipe state.
**Lost:** fiber multiplicity under membership-only projection, or recipe
behavior under scalar modular-order projection. **Sidecar:** valuations and
full transfer matrices. No map from phi fibers to strong-tournament classes
has been proved.

## 3. Recovered exact smooth value-monoid elasticity

**INHERITED PROVED, rederived here for clarity; not new mathematics.**
The atom/elasticity result is already in
`05-knowledge/results/session160_ordered_programs_20261103/axioms/NOTE.md`;
`05-knowledge/results/ordered_response_ten_20261103/NORMAL_FORMS.md`
extends it to ten rewriting schemas, catenary degree5, delta set{1,2,3},
exact lengths for3969^k, and every rational elasticity in[1,5/2] on
the actual{3,5,7}-supported spectrum. The following narrower{3,7}
calculation is a reusable proof exposition, not a novelty claim.
 Let M be the multiplicative monoid of attained values
3^a7^b, with a,b>=0. The preceding inherited construction proves

    M={3^a7^b:a,b>=0}\{7,21}.

An atom here means an element of this NUMERICAL MONOID which is not a product
of two nonunits of M. It does not mean a strong tournament. The atoms are
exactly3,49,63,343.

*Proof.* For b=0, use3. For b=1, a>=2 is required, and63 absorbs the lone7.
For b>=2, every b is a nonnegative combination of2 and3, so use49 and343
and the remaining3 factors. The four displayed generators cannot split:
3 is prime; splitting49 or343 would require7; splitting63 would require7
or21. This proves both generation and atomicity.

Every unordered atomic factorization is uniquely a tuple(u,v,w,t)>=0 with

    3^u 49^v 63^w 343^t=3^a7^b,
    u+2w=a, 2v+w+3t=b.

Thus its COMPLETE set of lengths is

    L(a,b)={a+v-w+t:2v+w+3t=b, 2w<=a, v,w,t>=0}.

For a nonunit define elasticity rho=max L/min L. Give the exponent pair the
additive weight W(a,b)=a+b/2. Its weights on the four atoms are1,1,5/2,3/2.
Every length ell therefore satisfies W/(5/2)<=ell<=W, proving rho<=5/2.
At63² the factorizations63*63 and3^4*49 have lengths2 and5, proving the
monoid elasticity is EXACTLY5/2. More sharply rho(63^(2k))=5/2 for every k>=1:
its extremal lengths are2k and5k, by the same weight inequalities.

This supplies an exact obstruction to transferring “unique prime
factorization” into the Chan–Pak value monoid. Each ordinary integer still
has a unique factorization into integer primes; the nonuniqueness appears
because7 is unavailable as a standalone tournament value and must be bundled.
For a FIXED tournament, its strong components remain unique: the two
factorizations above describe different join constructions, not a violation
of unique strong decomposition. Atomic value factors cannot use outside primes,
since all exponents are nonnegative; this is why restricting to the smooth
sector does not hide a factorization through another prime.

`verify_smooth_lengths.py` checks the length formula independently against
multiplicative dynamic programming on all a,b in0,...,12. This is a regression
check, not the proof of the all-parameter statement.

## 4. Correction ledger and stopping rules

- **MISTAKE-658:** Mersenne traces do not make their values forbidden.31,63,127
  are attained. Reversal/vertex addition are not automatically a modular group
  action. A formal packing vector need not be realizable.
- **MISTAKE-651 / THM-4746:** connected graphs also obey the corresponding
  congruence beyond their own threshold; they are not a negative control.
  Non-squarefree Camion cases must not be silently omitted.
- **Session218 / THM-4884:** finite missing multiplicities below625 do not
  refute global Ford-type multiplicity occurrence. Word fibers, numerical
  profile fibers and tournament-class fibers are different objects.
- **THM-4904:** the Carmichael chamber transfers a twin/involution/descent
  architecture, not totient values. Its fixed-depth singleton upper bound is
  O_K(sqrt(X)(log X)^floor(K/2)); a matching lower bound is not proved.
- **THM-4847:** the ALL-tournament rotation count on a Fermat number is blind
  to its primality. This must not be confused with the strong-count test, which
  does fail for various composite Fermat numbers. The three-state counting
  reformulation uses Pépin's criterion; it does not settle future Fermat primes.
- **THM-4940:** varying the packing evaluation base preserves the underlying
  packing polynomial but changes carries.1+3b is absent for integer b>=2,
  while1+10b is attained exactly at b=4,8. Binary7/21 is not a universal
  base-independent numerical obstruction.

The most promising next transfer is target-preserving construction with full
prime support or full endpoint state. Symmetry provides exact divisibility
certificates but, beyond C3, obstructs prime path counts in vertex-transitive
families. Any claimed completeness mechanism must therefore explain how it
leaves those symmetric carriers while retaining controllable arithmetic.

## 5. Complementary divisibility chart: inheritance for the session lanes

THM-4748's three-back-arc fan already realizes2^e-1 for every e>=5.
Its e=2k specialization and the classical gcd law for repunits can therefore
be combined in a tournament-valued divisibility chart; the scalar gcd identity
is classical, and the fan is inherited. THM-4801 already proves that
2^e-1 has power-of-two totient exactly for e in{1,2,4,8,16,32}.
Thus phi(4^k-1) is a power of2 exactly for k in{1,2,4,8,16}; this
classification is unconditional and does not assume no further Fermat primes.
The first composite Fermat factor blocks every longer initial product.

All source pins below identify the latest modifying commit visible at the
snapshot, not necessarily the discovery/first-publication commit.
`source_pins.json` gives exact repository paths, commit SHAs and blob SHAs.
