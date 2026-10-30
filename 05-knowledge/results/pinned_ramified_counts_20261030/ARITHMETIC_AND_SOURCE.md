# Count queries, normalization, and the 1001 branch

## Source audit: requested Hexagon paper

**CITED PREPRINT CLAIM / NOT REBUILT HERE.** Rudich and Rousseau,
[*Deterministic Factoring in Time N^o(1)*](https://hexagonmath.org/pdf/2610.00600v1),
Hexagon 2610.00600v1, deposited October 28 and displayed October 29, 2026.
This is a Hexagon identifier, not the unrelated arXiv identifier with
the same digits. Read the six-page PDF text, especially Theorem 1 and
Sections 3--4; this session did not download/build the registered Lean
development or reproduce its axiom report. PDF screenshots were requested
but the tool exposed only references; direct PDF download returned HTTP403.
No visually inspected local PDF or formal audit is claimed.

The paper claims one finite-alphabet deterministic Turing machine with
subexponential bit-length running time. Its described interface counts
factor pairs with prescribed leading bits; internal multiplication variables
have unique completions. An adaptive prefix search obtains the least factor.
Its network engine returns counts scaled by a power of two. The described
zero test uses small-modulus residues and a height bound on that scaled
integer. The paper separately labels its stronger quantitative rate a
conjecture. We use these as possible interfaces, not an imported proof of
ETH failure or a measured practical speedup.

## A1. A scale-independent zero test (PROVED; concurrent Session112 credit)

**Incoming precedence:** while this packet was in progress, mainline
`fea3b07542` published the SAME general cutoff and the SAME explicit dyadic
constant in [Session112 MODULAR_OBSERVERS, Section2](../session112_counting_interfaces_20261030/arithmetic/MODULAR_OBSERVERS.md).
We credit that checkpoint, retain the argument below as a local replay,
and do not count A1 as a distinct new result. Its effective-modulus theorem
also handles nonunit scales more generally than the unit-only account here.

Suppose an exact engine supplies V modulo queried moduli, where

    V=2^E C, E>=0, 0<=C<=2^b,

and the target is whether C=0. Query only odd integers m<=2(b+4).
Then

    C=0 iff V=0 modulo every queried m.                 (A1.1)

The query bound is independent of E. In particular one does not need to
materialize 2^E or even know E to use this zero test, provided the promised
scaling relation holds. This is elementary modular arithmetic, not a new
factoring engine or a claim of literature priority.

**Proof.** Put M=b+4 and L=lcm{odd m:1<=m<=2M}. Since 2 is a unit in
every queried ring, all residues vanish iff L divides C. Each prime-power
factor of binom(2M,M) is at most 2M: in the factorial valuation formula
each summand floor(2M/p^i)-2floor(M/p^i) is zero or one and there are at
most floor(log_p(2M)) summands. Thus the odd part of the central binomial
coefficient divides L. Its 2-part is at most 2M. Therefore

    L >= 4^M / [2M(2M+1)] > 2^b.

The first inequality also uses binom(2M,M)>=4^M/(2M+1). The final
inequality holds at b=0 and its ratio increases with b. Hence 0<=C<L
and L|C force C=0. QED.

More generally, for a known scale whose prime support is a fixed set P
of size h, query m<=2M coprime to every prime in P. The same argument gives

    L >= 4^M / [(2M)^h(2M+1)].                         (A1.2)

Choose M so this exceeds the unscaled count bound. The exponents of scale
primes do not enter. For fixed h, M=O(b+1) suffices. For P={2,p}, any odd
prime p, M=b+8 suffices. If scale residues are known, invert them in each
queried ring to recover C modulo L, hence exactly within [0,2^b].

**Hostiles and limits.** Modulo 2^10, scale 2^20 hides the distinction
C=0 versus C=1. Conversely, any fixed finite list of odd moduli has a
positive common multiple masquerading as zero once the count bound is
dropped. A height guard and unit-scale guard are both necessary. Weighted
signed counts need the corresponding absolute-height bound; the displayed
exact interval recovery is for nonnegative counts.

**Shared possible paper improvement, not an implemented result.** A1 can replace
the zero test described in Section 4 whenever its exact output really is
2^E times the desired count and arbitrary odd moduli are supported. This
removes E from the required modulus range. We have not audited the engine's
interface/costs in Lean and do not assert an improved overall factoring
rate. In particular it supplies no engine for our tournament problems.

**Actual tournament consumer.** Burnside numerators in P1--P2 are |G|
times an integer orbit count, and |G| is a known power of p. Modular
evaluation at moduli coprime to p permits exact division before CRT.
If an additional engine scale 2^E is present, avoid both 2 and p and use
A1.2. This changes the arithmetic evaluation of the SAME integer count;
it does not remove group transport, pin geometry, or the target height.

## A2. Inherit the exact arithmetic observer before seeking new analogies

The [Session78 inverse](../session78_fourth_power_aliquot_20261027/ALIQUOT.md),
Section 2, ALREADY proves that (s(n),s(pn)) reconstructs every positive
integer n for any known prime p, including p dividing n. For n=p^a u,
p not dividing u, its coordinates are

    x=s(n), c=sigma(u), y=s(pn)=p x+c,
    d=y-x=p^a[u+p s(u)], a=v_p(d),
    n=[p^(a+1)c-d]/(p-1).                             (A2.1)

An arbitrary pair still needs the exact realization guard sigma(u)=c.
[Session81](../session81_counterexample_families_20261028/ALIQUOT_OBSERVERS.md)
already extends this to prime-power dilations and exhibits a genuine
mixed-prime multiplier40 collision. None of those results is rediscovered
or claimed new here.

The requested chain has the following shared-prime7 observations:

| n | s(n) | s(7n) | c=sigma(n/7^v7(n)) | recovered v7(n) |
|---:|---:|---:|---:|---:|
|1001|343|2569|168|1|
|343|57|400|1|3|
|57|23|241|80|0|

The number400 has two exact roles: it equals sigma(343)=343+57, and also
s(2401). These are related by the prime-power divisor formula, not by
identifying a multiplication arrow with an aliquot edge.

The complete 19-element fibre s^-1(343), including578 and1001, is inherited
from Session78/82 and [Session111's merge](../session111_scale_connections_20261030/ALIQUOT_MERGE.md).
Our checker regenerates it using the sharper bound29241 (342 is not prime),
then verifies both p=2 and p=7 queries by independent divisor formulas.
For this finite fibre, modulus71 is the smallest positive modulus that
separates all19 doubled outputs; all smaller possibilities are checked.
That number is only a finite receipt, not a purported prime-distribution law.
Exact reconstruction works without choosing a modulus.

## A3. Conditional transfer to a previously proved inverse algorithm

**CONDITIONAL on the preprint's claimed factoring machine.** Fix a known
multiplier K>=2. For integer target x>=2, the complete fibre

    {n>=1 : s(n)=x and s(Kn)=y}

can be found in deterministic x^o(1)+poly(log(2+|y|)) bit operations.
This is a complexity corollary of the INHERITED Session81 algorithm, not
a new inversion algorithm or an unconditionally audited complexity result.

Indeed every source lies below N=(x-1)^2. For every profile a_p=v_p(n)
at primes dividing K, with product p^a_p<=N, compute the exact rational

    R_a=product_(p|K) (p^(a_p+b_p+1)-1)/(p^(a_p+1)-1), b_p=v_p(K).

R_a>K, and the only candidate in that profile is

    n=(y-R_a x)/(R_a-K).

There are O_K((log x)^omega(K)) profiles. Retain only positive integers
within N, with the asserted profile and both divisor sums correct. Before
candidate arithmetic, negative y or y>sigma(K)(N+x) can be rejected by
the elementary inequality sigma(Kn)<=sigma(K)sigma(n). Thus only O_K(log x)
bits of y survive. Conditional subpolynomial factoring validates every
candidate and Kn; multiplying by the profile count preserves x^o(1).
The paper's uniform machine quantifiers are needed for that conclusion.
The exponent-vector enumeration is NOT uniform polynomial-size when K
varies freely. The cases x=0,1 have the separate inherited branches.

## A4. What genuinely transfers

| Source -> target | Exact map/predicate | Lost information and required guard |
|---|---|---|
| Pin geometry -> attachment counts | Pointwise stabilizer, then its cycle index P1 | Full configuration orbit and complement relation are not reconstructed |
| Full pinned group -> rigidity test | Take a nonidentity order-p power, then test free-fibre shifts | Valid for triviality; LCA still required for orbit counts |
| Scaled count -> unscaled zero predicate | Restrict to unit moduli, apply A1 | Needs a bound on the true count; cannot erase physical torsion states |
| Aliquot source -> two dilation responses | A2.1 gives an exact inverse | Need known operation and realization guard; one aliquot value is insufficient |
| Conditional factoring -> inverse pair validation | Factor each inherited profile candidate | Does not supply a faster general tournament quotient or prove convergence |

The shared principle is to derive a missing coordinate from a legal
operation when possible, and retain it when not. The paper's prefix
search, the root group's order-p witness, and the aliquot valuation
inverse are three different maps with this property, not an isomorphism
between their objects. In particular prime-dilation clocks and the
ramified p-power clock are not finite-field Frobenius.
