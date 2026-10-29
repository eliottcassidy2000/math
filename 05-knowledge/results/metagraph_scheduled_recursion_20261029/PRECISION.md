# Precision truncation and a sharp reciprocal correction

**PROVED elementary deductions from the inherited positive Burnside
channel.** Exact finite verification is separate. No claim of priority for
classical Burnside theory, 2-adic arithmetic or partition truncation.

## P1. Source and normalization

For integers r>=1 and s>=0, let F_n(r,s) count r independent tournament
layers on common vertices with s binary vertex attributes, modulo common
relabeling. The ordinary count a_n is F_n(1,0). The aligned consecutive
pair count is F_n(2,1), and F_(N-1)(1,1) counts rooted N-tournaments.
These interpretations and the formula below are inherited from
[Session89](../session89_consecutive_tournament_pairs_20261028/counting/README.md).

For an odd partition lambda of n, write k for its number of parts, m_a for
its multiplicities, and

\[
E(\lambda)=\frac{n-k}{2}+\sum_{i<j}\gcd(\lambda_i,\lambda_j),\qquad
z_\lambda=\prod_a a^{m_a}m_a!.
\]

Then

\[
F_n(r,s)=\sum_{\lambda\ {m odd}}\frac{2^{rE(\lambda)+sk}}{z_\lambda}.
\]

Even cycles cannot fix a tournament: a power reverses an antipodal pair.
Keeping r>=1 is essential. The r=0 attribute-only problem does not inherit
this cycle filter. Occurrence pairs in E include repeated cycle lengths.

The inherited [low-bit theorem](../session89_consecutive_tournament_pairs_20261028/counting/LOW_BITS.md)
proves integral normalization by

\[
B_n=\begin{cases}
0&n=0,\\r+2s-1&n=2,\\r(n-1)/2+s&n\text{ odd},\\rn/2+2s&n\ge4\text{ even}.
\end{cases}
\]

For odd n the quotient is odd. For even n>=4 its first bit is phi(n)/2
modulo2. This packet extends the precision, not those already-proved facts.

## P2. Bounded cycle number at fixed precision

The elementary inequalities

\[
\nu_2(z_\lambda)=\sum_a\nu_2(m_a!)\le\nu_2(k!),\qquad
E(\lambda)\ge(n-k)/2+\binom{k}{2}
\]

give the following lower bounds on the valuation of an individual
normalized summand:

\[
L_{r,s}^{\rm odd}(k)=r(k-1)^2/2+s(k-1)-(k-2),\quad k\ge3\text{ odd},
\]

\[
L_{r,s}^{\rm even}(k)=rk(k-2)/2+s(k-2)-(k-1),\quad k\ge4\text{ even}.
\]

Here nu_2(k!)<=k-2 for odd k>=3 and <=k-1 for even k. Confusing these
bounds fails at k=4, where nu_2(4!)=3. Each L strictly increases when k
is increased by2: the difference is2rk+2s-2>0 in its declared range.

**Theorem.** For b>=1, discard every partition with L(k)>=b. The remaining
sum gives F_n(r,s)/2^B_n exactly modulo2^b. Keep the one-cycle term at odd
n, or all two-cycle terms at even n>=4, and the explicit n=0,1,2 cases.

**Proof.** Work in Z_(2), the rational numbers with odd denominators. Each
discarded normalized summand is divisible by2^b, so their sum is zero
in this quotient. For a retained partition, split z=2^v z_odd and compute

\[
2^{rE+sk-B_n-v}\,(z_{\rm odd})^{-1}\pmod{2^b}.
\]

The exponent is nonnegative in this scope. The odd denominator has a unique
inverse modulo2^b. Sum these residues and use the inherited integrality of
the complete count. No approximation or rounding is performed. QED.

The largest retained k is O(sqrt(b/r)+1), independently of n. Enumerating
nondecreasing odd k-tuples with fixed sum uses at most O(n^(k-1)) visits;
each visit uses O(k^2) gcd interactions and integer arithmetic of precision
controlled by b and k log n. This is polynomial in **numeric n** for fixed
b,r,s, not polynomial in log n and not a bounded-size full-count algorithm.
When b grows to the full count's bit length, the advantage need not remain.

The implementation's recurrence is simply

\[
\mathcal P(n,k,a)=\bigcup_{a\le j\le n/k,\ j\ {m odd}}
  (j)\mathbin{\|}\mathcal P(n-j,k-1,j),
\]

with the one-part base case. It enumerates each admitted partition once.
The remaining gcd sum and multiplicity factorial are retained exactly.
This is a new precision-oriented evaluator built on inherited arithmetic,
not a new count formula or a structural decomposition of metagraph vertices.

### The older divisor quadratic form supplies a collision budget

Put M_d=#{i:d divides lambda_i} and

\[
\Delta(\lambda)=\sum_{i<j}(\gcd(\lambda_i,\lambda_j)-1)
=\sum_{\substack{d\ge3\\d\ {\rm odd}}}\varphi(d)\binom{M_d}{2}.
\]

The equality is the gcd/totient identity used in
[THM-516, the core-kernel phi reframe](../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md).
For odd n the **exact** normalized summand valuation is

\[
r(k-1)^2/2+s(k-1)+r\Delta-\sum_a\nu_2(m_a!).
\]

For even n>=4 replace the first two terms by
`r k(k-2)/2+s(k-2)`. Thus every retained k>=3 odd or k>=4 even partition
also satisfies

\[
\Delta\le\left\lfloor\frac{b-1-L(k)}r\right\rfloor.
\]

In particular each pairwise gcd is at most one plus this bound. The old
quadratic form now gives a **precision-dependent collision cutoff** as
well as the cycle-number cutoff. A partial tuple already exceeding the
collision budget may be rejected, because all later collision terms are
nonnegative and the factorial estimate already budgets every multiplicity.
The current implementation evaluates completed tuples; it does not claim
a benchmark for an unimplemented branch-pruning variant.

This cutoff restricts common divisors, not individual cycle lengths or
prime factors. A single arbitrarily long odd cycle has Delta=0 and must
survive. Its odd denominator still depends on its actual length modulo
2^b. A fixed prime wheel deleting long prime cycles would therefore
destroy precisely the leading term. Source size and denominator residues
remain part of the arithmetic carrier.

## P3. The first correction is sharp at every odd n


Let n>=3 be odd, A=r(n-1)/2+s, delta=2r+2s-1 and U=F_n(r,s)/2^A.
Then

\[
\boxed{\nu_2(nU-1)=\delta.}
\]

Equivalently,

\[
U=n^{-1}+2^\delta\pmod{2^{\delta+1}}.
\]

**Proof.** The one-cycle partition(n) contributes exactly1/n to U.
For k=3, let G be the sum of the three pairwise gcds and let v be the
valuation of the multiplicity factorial. Its normalized valuation is
r(G-1)+2s-v.

If all three parts are distinct, v=0 and G>=3, giving at least delta+1.
If exactly two parts equal a, write the third as c. Then v=1 and
G=a+2gcd(a,c). For a>=3 this is at least5, so the valuation exceeds delta.
For a=1 it equals3, giving valuation exactly delta. At every odd n>=5
there is precisely one such partition, (1,1,n-2), with contribution
2^delta/(n-2). If all three parts equal a>=3 the valuation is also higher.
At n=3 the unique triple(1,1,1) contributes2^delta/3, again with an odd
denominator.

For k>=5, P2's first value is8r+4s-3, at least delta+1, and it increases
thereafter. Thus exactly one non-one-cycle term has valuation delta.
Multiplication by odd n preserves its odd leading coefficient and sends
the one-cycle term to1. No other term can cancel that bit. QED.

The boundary n=1 is excluded: nU-1=0 there. The theorem supplies the exact
valuation of a **correction**, not a higher valuation of F_n itself.

Consequences:

- Ordinary tournaments: `v2(n a_n/2^((n-1)/2)-1)=1`; hence
  `a_n/2^((n-1)/2) = n+2 (mod4)` at odd n>=3.
- Aligned pairs P_n=F_n(2,1): `v2(n P_n/2^n-1)=5`. The single n-cycle
  already gives the first five normalized bits; the sixth is changed by
  the repeated-unit triple. This does not follow by squaring a_n.
- The general number delta of initially correct inverse bits depends on
  both tournament layers and attributes. The binary part of2! is the
  reason the correction occurs one bit earlier than for distinct triples.

## P4. Literal consequences for the root-switch recursion

The inherited parent-only root-switch operator Q at order n has one row
for each pointed class, and entry Q_ij counts vertices yielding pointed
class j after re-rooting class i. Each ordinary class U contributes a block
`1 (m_1,...,m_t)` with m_i its vertex-orbit sizes. Therefore

\[
Q^2=nQ,\qquad\operatorname{tr}Q=n a_n.
\]

See [Session79's root-switch theorem](../session79_tournament_enumeration_20261027/ROOT_SWITCH.md).
At odd n>=3, P3 proves

\[
\operatorname{tr}Q/2^{(n-1)/2}=3\pmod4.
\]

For an even child order N>=4, the number R_N of pointed states is
F_(N-1)(1,1). Removing the distinguished vertex leaves a tournament whose
vertices have a binary direction to the root, modulo its actual action.
Thus P3 also proves

\[
\nu_2(R_N)=N/2,\qquad
\boxed{\nu_2((N-1)R_N/2^{N/2}-1)=3.}
\]

These are cheap all-size checks on the generated root banks and their
weighted reset trace. They do not determine Q's blocks, transport maps,
individual child groups or native black/blue edges. In particular dividing
trace(Q) modulo a power of2 by even n requires retaining nu_2(n) additional
digits; nonunit division is not a free operation.

The checker constructs all labelled tournaments and rooted isomorphism
orbits through5, forms Q literally and checks its entries' odd weights,
quadratic identity and trace. It also compares the existing complete bank
through10 with independent channel counts and the new congruences. This
is not a new regeneration or audit of all ten-vertex representatives.

## P5. What the precision carrier loses

The source is the global symmetry sum, the target is its specified residue,
and the map is projection in Z_(2) followed by reduction modulo2^b. It
preserves exactly those requested digits and discards whole high-valuation
partition contributions. Required sidecars are n,r,s,b, odd cycle types,
gcd incidences and multiplicity-factorial valuations.

The12-vertex81/27 examples have the same labelled cyclic triples but
different joint transporters. Their per-object attachment data is never
present in F_n(r,s), even before reducing precision. P2 consequently cannot
serve as a smaller carrier for arbitrary strong-atom composition. A native
port analogue would need its own valuation lower bound and realization
theorem. The earlier12/22 clock proof similarly separates recovering a
generator from recovering a chosen observer's amplitudes; neither is a
general permission to discard an attachment phase.
