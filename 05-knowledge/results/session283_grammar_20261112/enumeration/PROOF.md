# Complementary elimination for stratified tournament symmetry corrections

**Status: PROVED identities; implementation and timings are FINITE-EXACT controls.**
This counts ordinary tournament isomorphism classes A000568. It does not classify
all classes individually, and supplies no Hamiltonian-path spectrum theorem.
The cycle-type sectors are not disjoint isomorphism-class families: fixed-point
sets of different permutations may overlap; only their Burnside aggregate is
the class count.

## Inheritance and precise advance

The odd-cycle Burnside gas, its external-field recurrence, and maximal
fixed-mass projective quotient are inherited from Session270
`enumeration/PROOF.md` and Session98's response proof. Session279
`cross/PROOF.md` proves the ground-defect moment certificate with cells
`(M,delta,R)`; its producer explicitly listed all odd partitions. Session279's
efficient fixed jets did **not** produce those stratified cells.

The present new implementation changes the variable-elimination order:
**retain the defect-bearing long core; sum the defect-free backbone exactly
behind its finite response interface.** It produces exactly the same richer
packet, not merely the unstratified moments. This is a missing computational
connection between two inherited theories, not a new Burnside formula or
an assertion of polynomial-time full enumeration.

## 1. The conditional family

Fix n>=0 and a finite backbone B of distinct odd lengths >=3. Put S={1} union B.
Every odd cycle multiset has a unique decomposition lambda=mu disjoint-union nu,
where mu uses lengths outside S and nu uses S. Let

    d=|mu| (sum of its lengths),       N=n-d,
    q_j(mu)=sum_(a in mu) gcd(a,j), j in S.

The envelope weight of mu is

    w(mu)=2^[sum_(a in mu)(a-1)/2 + sum_(a<b instances) min(a,b)] / z(mu).

Here z(mu)=product_l l^m_l m_l!, and instance pairs include equal lengths.
For each pair of instances put e=min(a,b)-gcd(a,b); keep the strictly positive
values in D(mu). These are **exactly** the remaining defects in Session279:
every pair involving nu is already integrated exactly.

Define the short-family partition function

    F_S(N;q) = sum_(nu on S, |nu|=N)
        2^[E(nu)+sum_j m_j(nu)q_j] / z(nu),

where E is the exact Burnside orbit exponent. Then the contribution of all
completions of mu to unshifted moment k is precisely

    w(mu) e_k({1-2^(-e): e in D(mu)}) F_S(N;q(mu)).       (1)

Proof: E(lambda)=E(mu)+E(nu)+sum_j m_j(nu)q_j, z factors
because the permitted length sets are disjoint, and the defect list depends
only on mu. Summing nu gives (1). Neither scalar counts nor orbit types are
assumed multiplicative under tournament substitution: this is the exact
factorization of a **permutation cycle-type sum**.

Thus cores having equal `(N,q)` can merge their entire cell-to-moment maps
before completing them. Cores with unequal fields must not merge merely
because their current mass or moments match. Equal class counts are not equal
continuation interfaces; the same loss discipline matters in Chan–Pak, but
there is no direct realization transfer from this Burnside field.

## 2. Exact finite response recursion

Write S={s_0=1<s_1<...<s_t}. For i>=1, choose m cycles of length s_i. Then

    F_i(N;q)=sum_(0<=m<=N/s_i)
       2^[m(s_i-1)/2+s_i binom(m,2)+m q_i]
       /(s_i^m m!)
       * F_(i-1)(N-ms_i; (q_j+m gcd(s_i,s_j))_(j<i)).     (2)

At i=0, only singleton cycles remain, so

    F_0(N;q_0)=2^[binom(N,2)+N q_0]/N!.                 (3)

All operations are exact rational arithmetic. This is the inherited response
recurrence, used in the **opposite elimination order** from the backbone-first
moment-jet producer. The memoization key is the complete `(i,N,q)`.
Every invocation is a sum over precisely the same finite set with precisely
the same weights, so caching is exact even when reached from different cores.

The cell tag is `(M,delta,R)` with M=#D, delta=min D (or0 if D empty),
and R=#{e in D:e>delta}. It stays outside the continuation recursion because
no short cycle changes it. Summing (1) for k=0,...,4 produces byte-for-byte
identical rational cell entries to the inherited all-partition producer.
The consumer is unchanged: no exact A000568 reference is supplied to it.

### Integer compilation of the short family

The default implementation stores G_i(N;q)=N! F_i(N;q). Its transition weight
is the integer

    N!/[s_i^m m! (N-ms_i)!] * 2^[m(s_i-1)/2+s_i binom(m,2)+m q_i].

The coefficient counts selections of m disjoint s_i-cycles on N labelled
points. The terminal value is2^[binom(N,2)+N q_1]. Thus there are no rational
intermediates inside the cached short recursion; divide by N! only when
rejoining the residual core. The normalization is an inherited factorial-
clearing principle applied to this opposite elimination order. The code checks
every coefficient's exact division, and compares this integer path with both
raw-rational and projective-rational paths on the small control suite.

### The inherited fixed-mass gauge can also be retained

The implementation optionally uses Session270's projective response quotient.
Put beta_j=(j q_1-q_j)/2, with beta_1=0. At fixed residual mass N,

    F_S(N;q)=2^(N q_1) F_S(N;(-2 beta_j)_j).

Multiply the core's moment map by2^(N q_1), group on `(N,beta)`, and use
normalized transitions

    exponent = m(s_i-1)/2+s_i binom(m,2)+m(N-ms_i)-2m beta_i,
    beta'_j=beta_j+m(s_j-gcd(s_i,s_j))/2.

The exponent may be negative; exact rational powers of two are used, never
truncating right shifts. The terminal normalized factor is2^binom(N,2)/N!.
Raw and projective producers are compared directly on all small controls.
This quotient reduces states in measured examples, but rational arithmetic
can erase the wall-time benefit; no per-run speed claim follows from state
reduction alone. The gauge is inherited, while its compatibility with the
full defect-cell packet is (1): the scalar factor multiplies every cell moment
by the same amount.

## 3. Cost and exactness boundary

For B={3,5,...,b}, only residual parts >=b+2 are listed. Their number through
mass n is

    C_b(n)=[x^n] (1-x)^(-1) product_(odd j>=b+2)(1-x^j)^(-1).

The extra factor sums all residual masses <=n; it is not a singleton cycle
in mu. This can be far smaller than the full number of odd partitions of n.
For each core the edge list has at most binom(floor(n/(b+2)),2) entries.
Moment cost is O(K times that edge count) rational operations.

At a fixed backbone, each q_j is a nonnegative integer <=n: processed cycle
mass bounds the sum of gcd responses. Consequently the short recursion has
at most O((t+1)(n+1)^(t+2)) raw states and O(n) outgoing choices per state.
This is a coarse bound, not a claim that it is achieved. The residual-core
count remains superpolynomial for fixed b; increasing b grows the backbone
state dimension. No unrestricted polynomial-time counting theorem follows.

All stratified moments and the resulting intervals are **identical** to
Session279's packet at the same n,B,K. This algorithm changes their cost,
not their mathematical exactness boundary. In particular B through11 remains
exact for n<=69 and first has ambiguity at70. Larger backbones or exact
residues are separate improvements, not hidden in the producer.

## 4. Certified completion and the two different outputs

Given the unchanged consumer's integer interval [L,U], let
b=bit_length(U-L). The inherited Session277 modular engine supplies the exact
residue r=A000568(n) mod2^b independently. Then

    A000568(n) = L + ((r-L) mod2^b)

provided the result is <=U; uniqueness follows from U-L<2^b. The program
checks this, and only afterwards compares against the inherited exact cache.
A tight interval alone does not disclose missing low bits, and a low-bit
residue alone does not disclose the huge integer. Both certificates are
necessary unless L=U.

This is a useful analogy for Chan–Pak's missing absolute anchor: complementary
views must have a proved singleton intersection. It does not turn class-number
recovery into existence of a tournament with a prescribed path count.

## 5. Controls and limitations

`residual_compiler.py` compares the full packet (all cell moments, not just an
interval endpoint) against the direct odd-partition producer for n=0..24 and
four backbones, including the noninitial set{3,7}; it also compares n64/70.
`complete.py` performs no-reference interval/residue completion, then checks
against independently inherited exact counts and measures matched baselines.

Runtime observations are machine-specific. Faster packet production is not
necessarily faster than an optimized **full exact count**, and that distinction
is recorded rather than treating a different workload as a speed benchmark.
