# Prime-power contraction congruences for strong tournaments

**Status: PROVED; independently audited in Session25.**

For every odd prime p, every r>=1, and every m>=1, let epsilon=(2/p). Then

s(mp^r) = epsilon^m s(mp^(r-1)) (mod p^r),
D(mp^r) = epsilon^m D(mp^(r-1)) (mod p^r),

where s(n) counts strong labelled tournaments, s(1)=1, c(n)=2^((n-1)/2) for odd n and0 for even n, and D(n)=s(n)-c(n).

The prime-square theorem is the case m=1,r=2. Literature novelty is not claimed.

## 1. The order-p subgroup captures every short orbit

Put n=mp^r and N=mp^(r-1). Let the cyclic group G=C_(p^r) act on n labels as m disjoint regular cycles of length p^r. Its unique subgroup H of order p has N vertex orbits, each of size p.

On the set X of all strong tournaments on these labels, a G-orbit has size less than p^r exactly when its point stabilizer is nontrivial, equivalently contains H. Therefore the complement of X^H is a union of full p^r-element G-orbits and

s(n) = |X^H| (mod p^r).

No assumption that H-fixed points are G-fixed is made.

## 2. Exact semicomplete quotient polynomial

An H-invariant tournament restricts on each of its N blocks to one of

a=2^((p-1)/2)

circulant tournaments of prime order p. Each block is strong, since any nonzero connection-set step forms a directed Hamiltonian cycle. Between two distinct blocks there are p edge orbits under H, giving2^p independent invariant orientation patterns. Exactly two are uniform, one for each direction. The remaining

w=2^p-2

patterns give at least one arc in each direction between the two blocks.

Collapsing the strong blocks produces a semicomplete digraph on N labelled vertices. Each pair is assigned a directed arc or a digon. Strongness is preserved in both directions: quotient paths lift by using internal strong paths inside the blocks. If

S_N(z)=sum_(Q strong labelled semicomplete on N vertices) z^(number of digons of Q),

then the exact number of H-fixed strong tournaments is

|X^H| = a^N S_N(w).

Write S_N(z)=s(N)+sum_(j>=1) B_j z^j.

## 3. The digon coefficient valuation

Let K=C_(p^(r-1)) act on the N quotient labels as m disjoint regular cycles. For r=1 this is the trivial group, and all subsequent bounds remain valid.

Consider its action on the finite set of strong semicomplete digraphs with exactly j>=1 digons. If an object has stabilizer of order p^b, that stabilizer acts freely on unordered pairs of distinct vertices:
- every nonidentity element of K is fixed-point-free on vertices;
- an element preserving an unordered pair must either fix both endpoints or swap them;
- fixing either endpoint contradicts semiregularity;
- swapping them makes its square fix a vertex, forcing its square to be identity; as K has odd order, the element was identity.

Hence each orbit of its stabilizer on unordered pairs has p^b elements. The object's digon set is invariant under this stabilizer, so p^b divides j. In particular b<=min(r-1,v_p(j)). The full K-orbit has cardinality p^(r-1-b), divisible by

p^max(0,r-1-v_p(j)).

Summing orbit cardinalities proves

v_p(B_j) >= max(0,r-1-v_p(j)).

As p divides w by Fermat,

v_p(B_j w^j) >= j + max(0,r-1-v_p(j)) >= r.

For the last inequality: if v_p(j)<=r-1, use j-v_p(j)>=1; otherwise j>=p^r>=r. Thus every positive-degree contribution vanishes modulo p^r and

S_N(w) = s(N) (mod p^r).

The zero coefficient B_j=0 causes no issue and may be assigned infinite valuation.

## 4. The scalar lift

Euler's criterion gives a=epsilon(1+pu) for an integer u. The elementary principal-unit congruence

(1+pu)^(p^(r-1)) =1 (mod p^r)

follows by induction: raising a number1 modulo p^k to the pth power gives1 modulo p^(k+1). Since p^(r-1) is odd,

a^N = epsilon^(m p^(r-1)) (1+pu)^(m p^(r-1))
     = epsilon^m (mod p^r).

Combining all preceding identities proves

s(mp^r) = epsilon^m s(mp^(r-1)) (mod p^r).

If m is even then both relevant c-values are zero. If m is odd then N is odd and

c(pN)=a^N c(N)=epsilon^m c(N) (mod p^r).

Subtracting proves the statement for D. QED.

## 5. Consequences and limits

1. Taking r=2,m=1 gives D(p^2)=epsilon D(p) modulo p^2. Thus p^2 passes the strong/circulant test if and only if p^2 divides D(p).
2. More generally, for each fixed m>=1 the sequences
   epsilon^(-mr) s(mp^r) and epsilon^(-mr) D(mp^r), for r>=0,
   are Cauchy in Z_p: consecutive differences are divisible by p^r. Hence they have p-adic limits. This proves convergence, not a formula or nonzero value for either limit.
3. For r>2 one must not replace the criterion at p^r by “p^r divides D(p)”. The recurrence relates level r to level r-1 modulo p^r; earlier-level recurrences hold only at their smaller moduli. Iteration loses precision.
4. Odd p is essential to the unordered-pair action and circulant block argument. There is no asserted p=2 version.
5. This counts tournament objects, not Hamiltonian paths within one tournament. It does not imply individual presentation spectra, spectral gaps, S2/S3 identities, or metagraph wiring.
6. The semicomplete quotient is an exact bookkeeping device for mixed interblock arcs. No assertion is made that those arcs are uniform in both directions.

## 6. Independent tests

The separate audit.py uses a fresh modular Wright recurrence, independent of the orbit proof. It checks both s and D lifting for:
- p=3, r=1..5, m=1..4;
- p=5, r=1..4, m=1..2;
- p=7, r=1..3, m=1..3;
- p=11, r=1..2, m=1..5.

These include even m, m divisible by p, r=1, and p=3. See audit.out for the exact completed checks. The proof does not depend on their finite scope.

## Provenance

The parent/root proposed the full cyclic p-power action, semicomplete weighting, stabilizer valuation, and scalar lift. This independent reviewer checked the short-orbit equivalence, block contraction, stabilizer action on unordered pairs, valuation inequality including r=1, c's even-m boundary, and p-adic normalization, then reproduced the stated finite controls. Novelty against published enumeration literature remains unchecked.


## Reconciliation with concurrent audit 3596c5d0

The incoming audit's Lemma B proves the additive law s(n+p(p-1))=(2/p)s(n) modulo p^2 for n>=3, and already deduces the prime-square passing equivalence. Preserve its proof and credit. The present theorem supplies an independent orbit proof at squares and a multiplicative dilation law at all prime powers. These statements have different domains: this theorem compares mp^r and mp^(r-1), whereas the incoming higher-power conjecture compares n+phi(p^k) and n at arbitrary sufficiently large n. This theorem does not settle that additive conjecture. Incoming Lemma C also proves the connected-graph tail law, so no uniqueness claim for ordered condensation should be inferred from the tournament proof.


### Subsequent distinct additive proof

ADDITIVE_PRIME_POWER_PROOF.md now proves the additive HYP-9358(d′) law by a separate filtered-Hurwitz argument. This does not change the scope of the multiplicative theorem above; the earlier distinction remains essential.
