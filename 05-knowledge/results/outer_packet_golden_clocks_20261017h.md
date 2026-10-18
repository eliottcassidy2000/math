# Signed outer packets, averaged carry costs, and two complementary golden clocks

2026-10-17, clock_duality_oct15. **PROVED, scoped, self-audited:** the
finite signed-packet inequality, the improved depth/scale remainder, and
the golden two-generator coverage statement below. **PROVED from a named
inherited input:** uniform o(N) cancellation of this one outer packet.
**FINITE-EXACT:** the declared integer, pair-chain and finite-ring censuses.
**OPEN:** a power-saving signed bound for the full ROOT current, ROOT energy,
source-specific universal entry, and Collatz. No independent-referee or
novelty claim.

## Inheritance and concept board

Anchor: retain the isolated outer packet's signs. Niche: golden v-squared
and quadratic-conjugate geometry, in the [companion](golden_square_contexts_20261017h.md).
Wildcard: the user's Kathara sequence, now explicitly defined as
0,1,1,2,4,8,16,... . Board: **marked source / integer clock / signed triple /
affine carry / norm-labelled golden orbit / shared geometric port**.

Closest mechanism: P1--P4 of [positive ROOT transport](root_positive_transport_20261017g.md).
The hostile is P4's actual isolated mass N/2^d+O(1), which no other word
cancels in its terminal columns. The corrected near miss is replacing a
signed discrepancy by all unmerged mass: [signed resonance, section3](collatz_signed_resonance_pair_mass_20261013b.md)
already measured a million-fold loss in that substitution. The least-used
sidecar here is the **initial translation of the pair chain**, which grows
with outer depth; a bound for the adjacent pair is not uniform in that datum.

The inheritance pass recovered [golden prime clocks, section3](golden_prime_clocks_20261003.md),
which already proves pi(3^k)=8*3^(k-1), and
[unit lift clock, U1](unit_lift_clock_20261015a.md), which is stronger than a
new small orbit census. We use those mechanisms and add the explicit
complementary scalar generator and norm-labelled dictionary. The Fourier
renewal estimates in [five mirrors, section2d](collatz_five_mirrors_20260929.md)
concern a different observable; their conditional rate is not imported here.
The later fetch brought S41's energy ledger in a430c818f9. Section7 integrates
its fibre decomposition and repairs the distinction between coefficient
energy and ROOT-projected current energy, and between primitive and ordinary
Syracuse offsets. Earlier S40's Farey sweep and corrected Paley ranking were
already integrated by21dbf2b400; no ranking heuristic is a premise below.

## 1. The outer packet is a three-point signed difference

Retain the shortcut T and ROOT-supported antifunction f from P1--P4:
f(n)=(-1)^tau(n) if n reaches1, zero otherwise; f(Tn)=-f(n), even on1<->2.
Let S be the positive odd integers other than5 modulo8, w(n)=3[n=2 mod3]-1,
and put

    I_(N,d)=((2N+1)/2^d,(4N+1)/2^d],   H=2N/2^d,
    O_(N,d)=(-1)^d sum_(h in I, 3^d h-1 in S)
                         w(2^d h-1) f(3^d h-1).          (1)

This is the actual isolated outer contribution to L^d H_N, not the whole
ROOT current. Its original source is2^d h-1, terminal3^d h-1, and literal
shortcut clock d. The exact source weight is

    w(2^d h-1)=3[3 divides h]-1.                          (2)

The allowed h are even, excluding h=6 mod8 for d even and h=2 mod8 for d
odd. Thus there are three allowed classes modulo8. Within each one,
partition the integers into triples

    h0-8, h0, h0+8,       h0=0 mod3,

whose weights are -1,2,-1. Each complete triple contributes
2f(t)-f(t-8*3^d)-f(t+8*3^d). At each end of I there can be one partial
triple per class. Its absolute weight is at most3; all partial triples
therefore cost at most18. This is an exact grouping, with no random signs.

## 2. A finite, depth-labelled coupling invoice

Let q_K(e) be the Haar probability that the exact affine pair chain from
(k,e)=(0,e) has not absorbed at (0,0) by K shortcut steps. The definition
counts both absorbing self-edges and all branch multiplicities. It requires
no assertion that any integer reaches ROOT. Swapping the two sources shows
q_K(-e)=q_K(e).

**O1 (PROVED finite inequality).** For every N>=1,d>=0,K>=0,

    |O_(N,d)| <= N/(3*2^d) [2q_K(3^(d+2))+q_K(3^(d+3))]
                       +12*2^K+30.                      (3)

Also |O_(N,d)|<=N/2^d+12 by P4's absolute packet count.

Proof. In the three primitive terminal classes1,3,7 modulo8, the first
three shortcut steps have respectively2,2,3 odd steps. Terminals differing
by8*3^d follow the same three parity bits. Their new translation is exactly
3^(d+2),3^(d+2),3^(d+3), respectively. The complete-triple centre parameter
has an odd affine coefficient after these three steps (a power of3), so
it visits each residue modulo2^K once in every2^K consecutive centres.

For m consecutive centres, the number of nonabsorbed comparisons is at
most m q_K(e)+2^K. On absorption the two f values agree: both clocks have
length3+K, and f is an antifunction on all sources. Each unresolved
comparison costs at most2. There are two comparisons per centre and three
classes. Their contribution is at most4 sum_r m_r q_K(e_r)+12*2^K.
Each m_r<=H/24+1. Add the at-most18 boundary weight and bound
4 sum_r q_K(e_r)<=12. This proves(3). Absorption on a cylinder is an affine
identity; an isolated accidental meeting of two evaluated integers is not
counted as an absorbed cylinder. No first-hit clock has been replaced by a
parity-only label in the source certificate.

**O2 (PROVED from the inherited almost-sure input).** Define, for D>=1,

    C_(K,D)=max(2^(-D),
       max_(0<=d<D) [2q_K(3^(d+2))+q_K(3^(d+3))]/(3*2^d)).

Then, uniformly over every outer depth,

    sup_(d>=0) |O_(N,d)| <= N C_(K,D)+12*2^K+30,           (4)
    sup_(d>=0) |O_(N,d)|/N -> 0 as N -> infinity.          (5)

For(4), use(3) for d<D and the absolute count for d>=D. For(5), first fix
D. [THM-4581, haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely](../../01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md),
statement3, supplies q_K(e)->0 for each of the finitely many translations.
Choose K, then send N to infinity, then increase D. Its sketch-level rate
is not used. Equation(4) survives independently of that inherited theorem.
The simultaneous all-depth statement concerns **only this outer packet**;
it supplies no power rate, and says nothing about entry of a fixed source.

### Why this is not yet the power-saving estimate

The exact pair census at K=24 gives:

| e | numerator of q_24(e), denominator2^24 |
|---:|---:|
|3^2|14676583|
|3^3|15506424|
|3^4|16136471|
|3^6|16695569|
|3^8|16774181|
|3^10|16777216|

A finite uniform bound imported from e=1 would miss this deterioration.
The signed packet is much smaller than that absolute invoice. The complete
census uses N=2^k,4<=k<=20,0<=d<=16, certifying every original source through
4194304. At N=2^20:

| d | signed O | absolute outer mass |
|---:|---:|---:|
|0|1574|1048576|
|1|376|524288|
|2|512|262144|
|4|674|65536|
|8|160|4096|
|12|4|256|
|16|-2|16|

These are integer observations, not a decay fit or uniform estimate.
The remaining helper is a **signed survivor discrepancy** on the two
translations in(3), rather than merely a faster upper bound for q_K.
One possible target is |O_(N,d)|<=C sqrt(N/2^d) times a fixed power of
log(2N), uniformly in d. Even that would close only this packet, leaving
all other word families and the full ROOT-current target to be paid.

## 3. Averaging carry costs extends the uniform approximation horizon

Use A^(d)(N) and gamma_d exactly as in positive-transport P3. Set
m=ceil(log2(4N)), K=m+2d. Its inverse words have initial braid index k>=0,
letters A_i>=1, total cost 2k+sum A_i<=K, and source alpha*t+beta.
The proof there actually gives the sharper intermediate inequality

    |A^(d)(N)-N gamma_d|
       <=36 sum_words (|beta/alpha|+2)
         +N*(5/2)/(1-u^2/4)*(u/(2-u))^d*u^(-K),          (6)

for1<u<2. The two terms are respectively affine threshold movement plus
progression quadrature, and the discarded homogeneous tail. The preceding
note bounded every |beta/alpha| by the largest word's value. We sum them.

For1/2<z<1 define

    F=z/(1-z),        G=3z/(2-z).

An inverse letter has slope2^A/3 and normalized offset -c_A/2^A, where
c_A=(3+2^a-2^A)/3, a=1 or2 according to A's parity. Directly
|c_A|/2^A<=1/2. A prefix multiplies a later normalized offset by the
product of3/2^A along that prefix. The initial braid adds at most
(1/3) times this product over the complete word. Consequently

    sum_words (|beta/alpha|+2)
      <= z^(-K)/(1-z^2) *
         [2F^d + F/(2(F-G))*(F^d-G^d) + G^d/3].          (7)

Proof: weight every word by z^(2k+sum A_i). Summing letters gives F for
ordinary suffixes and G for carry-weighted prefixes. Sum the location of
the offset over the d letters, and sum the braid k geometrically. On the
truncated set, 1<=z^(-K)z^cost. Extending this positive majorant to all
words proves(7). We have not changed a source, merged unequal clocks, or
used cancellation in an error term.

**O3 (PROVED uniform remainder).** Put

    r=(12^12/11^11)^(1/10),    log2(r)=0.49658022036436... .

For every N>=1 and0<=d<=floor(m/10),

    |A^(d)(N)-N gamma_d| <=600 r^m
                         =O(N^0.49658022036436...).       (8)

There is no extra log N factor in this bound. To prove it, set z=11/12:
F=11, G/F=3/13. The bracket in(7), divided by F^d, is

    2+(13/20)(1-(3/13)^d)+(1/3)(3/13)^d <=53/20.

Thus the first term of(6) is at most
36*(53/20)*(144/23)*r^m. Choose u=7/4 in the tail. Since N<=2^(m-2), it
is at most (8/3)*[(8/7)^10*(16/7)]^(m/10) <=(8/3)r^m.
The constants sum to599.953623...<600. The exact integer inequalities

    12^24 < 2^10*11^22,
    (8/7)^10*(16/7) < 12^12/11^11

establish the sub-square-root exponent and tail comparison. This extends
the horizon from m/12 to m/10 and improves the error at shorter depths.
It remains an approximation of the absolute envelope. Its outer lower
bound still requires roughly m/2 depth for an O(sqrt(N)) envelope; (8)
does not close that gap or prove a full signed ROOT bound.

## 4. Kathara, Fibonacci, and the word budget

The natural cumulative recurrence realizing the user's supplied sequence is,
with index0, a_0=0,a_1=1 and

    a_n=sum_(j<n) a_j =2^(n-2) for n>=2.                  (9)

For t>=1, a_(t+1) counts all ordered compositions of t into positive parts.
With d parts there are binom(t-1,d-1) words. This is exactly the letter
alphabet A_i>=1 of positive-transport P1. Its clock law is

    Pr(sum A_i=t)=binom(t-1,d-1)*2^(-t).                  (10)

The binary representation is a separator mask on the t-1 interior gaps.
It preserves both the ordered word and its total clock; retaining the
number of separators recovers d. Restricting parts to{1,2} gives F_(t+1)
compositions. This is the Fibonacci restriction of the same word language,
not an assertion that a deterministic Collatz source has random letters.

The distinction explains a useful obstruction. Summing probabilities over
**all depths** at a fixed positive cost t gives1/2. Summing that again over
all t diverges, because prefixes are being counted repeatedly. An accepted
stopping language must be prefix-free (or retain its renewal multiplicity)
before a Kraft/probability bound is used. The entire alphabet has cumulative
doubling; a two-letter controller has Fibonacci growth but omits long
valuation branches. Universal coverage cannot be inferred by dropping them.
The initial braid in(7) is a separate even-cost register, giving1/(1-z^2).
This is an exact connection of the requested recursion to the error proof.

## 5. Fibonacci modulo9 and the complementary doubling generator

Take Q=[[1,1],[1,0]], acting on(a,b), corresponding to multiplication by
phi on b+a phi. Its norm is b^2+ab-a^2. Exactly

    Q^8=I+3[[11,7],[7,4]],    Q^12=-I mod9,   Q^24=I mod9.

The matrix multiplying3 is invertible modulo3. The inherited lifting proof
therefore gives order8*3^(k-1) modulo3^k, on every primitive pair. At k=2
there are72 primitive pairs and three24-cycles. They are precisely the
three classes of a^2-ab-b^2 up to sign: {1,8},{2,7},{4,5} modulo9.

**O4 (PROVED, all k>=1).** On primitive pairs modulo3^k, the two commuting
actions Q and scalar doubling generate a transitive group of size
8*3^(2k-2). Q alone has3^(k-1) orbits, classified by norm up to sign.

Proof. The quadratic is irreducible modulo3, so primitive pairs are units.
Q has order8*3^(k-1); scalar2 has order2*3^(k-1), by v_3(2^2-1)=1 and
binomial lifting. Their cyclic groups intersect in exactly{I,-I}: if
Q^i=cI then c^2=(-1)^i, and -1 is not a square modulo3; the only solutions
to c^2=1 modulo3^k are +/-1. Both groups contain -I. Their product thus
has8*3^(2k-2) elements, the total number of units, proving transitivity.
For the norm classification, Q changes norm by -1. Each prescribed unit
norm has4*3^(k-1) points: four modulo3, then three lifts per level, since
the gradient of this nonsingular norm is nonzero on primitive pairs.
A pair of opposite norms therefore has exactly one Q-orbit's size.

For k=2 the clocks have lengths24 and6, with common half-turn Q^12=2^3 I.
Doubling sends a norm to four times itself, cycling the three Q-orbits.
This is a literal example where adding a commuting generator repairs a
coverage deficit that one recurrence cannot repair.

There is also a precise finite dictionary for(1): h mod24 -> Q^h(1,0)
is a bijection onto one24-cycle. Translation h->h+8 becomes Q^8, an
order-three lift fixing the mod3 vector; it preserves the three positions
of each signed triple. The character(2) and the primitive mask can be
carried as labels. The map discards the integer quotient floor(h/24), so
that quotient and d must remain sidecars. It does **not** intertwine the
Collatz map with multiplication by phi. O4 is universal coverage in a
finite residue ring, not entry of every integer into ROOT.

This distinction is consonant with
[THM-3487, two-twenty-four-state-fibonacci-bundles-cycle-type-obstruction](../../01-canon/theorems/THM-3487-two-twenty-four-state-fibonacci-bundles-cycle-type-obstruction.md):
its two24-state bundles have different cycle types until a seam class is
changed. We have specified an actual action here, not only a cardinality.
The older [THM-486](../../01-canon/theorems/THM-486-involution-modulus-is-pisano-modulus.md)
also records pi(24)=24; that modulus is different from the present modulus9.

## 6. Recovered number contexts, with operation types retained

| number | recovered context | useful mechanism and boundary |
|---|---|---|
|613|[THM-4094, Hamiltonian matching and two prime lanes](../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md), section7|A finite prefix plus multiplicative joins cannot produce a new prime atom. [THM-4097](../../01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md) subsequently supplies explicit strong witnesses for613: it is not a current gap.|
|248|[Heegner/E8 note, section6.2](ramanujan_heegner7_trivial_cycle_20261006.md)|The coefficient of E4/prod(1-q^n)^8 is240+8=248; cubing yields744. The theta contribution and denominator contribution are distinct data, not one unexplained constant.|
|365|[THM-1131, corrected PG213 gap claim](../../01-canon/theorems/THM-1131-pg213-one-twelfth-gap-refuted.md)|The exact runner family{1,...,12,364} has maximum28/365. It is a hostile to quotient-only global-gap claims; endpoint ownership and actual height matter.|

The identity613=248+365 is exact. The retrieved mechanisms do not turn it
into a relation among Hamiltonian counts, modular coefficients and LRC
maxima. In particular an additive decomposition cannot supply a missing
multiplicative prime generator. As a cheap hostile to a shared24-clock
interpretation, the exact Fibonacci periods are respectively1228,60,740.
No common24-period is present for these three integer moduli. Their useful
shared lesson is to name the operation and the information required to
make a construction reach its target; O4 and(9)--(10) are the concrete
new instances, rather than a numerical identification of the three values.

## 7. Incoming energy ledger: the surviving unconditional target

S41's [THM-4667, energy-ledger-of-the-signed-root-transport](../../01-canon/theorems/THM-4667-energy-ledger-of-the-signed-root-transport-partial-currents-merge-cross-terms-and-the-octave-trace.md)
identifies each P-fibre as sibling pairs {s,2s+1}, s=1 mod8, with opposite
signs, plus at most the singleton(2u-1)/3 when u=2 mod3. Write c_d=L^d H_N
and E_d=sum c_d(u)^2. Its coefficient merge identity is exact on every
finite source vector, whether or not that vector reaches ROOT.

**O5 (PROVED, scope repair).** For the global ROOT-supported f,

    R_u^(d)=sum_(P^d(s)=u) H_N(s) f(s)=f(u)c_d(u),
    E_d=sum_u (R_u^(d))^2 + sum_(u not in ROOT) c_d(u)^2.  (11)

The original incoming identification of E_d with the first sum alone
silently required the second term to vanish. The invariant relation is
f(s)=sigma_d(s)f(u); its inverse sigma_d(s)=f(s)f(u) holds only on ROOT.
An abstract disjoint component with f=0 and nonzero coefficient is the
cheapest type-level hostile; no positive Collatz counterexample is asserted.
The sibling identity INTRA=-2 sum c_d(s)c_d(2s+1) is unconditional; its
partial-current version accounts only for rooted siblings and needs the
unrooted complement in(11).

The target survives the repair. Each rooted member of the **finite** initial
support has finite P-depth. Their maximum is finite, although not supplied
by an algorithm here. Unrooted members never reach1. Thus c_d(1) eventually
equals the full ROOT current R(N), regardless of the unrooted components:

    R(N)^2 <= liminf_d E_d(N) <= sup_d E_d(N).             (12)

The same argument for a finite octave gives
D_X<=liminf_d sum_(X<=N<2X) E_d(N). Hence the proposed uniform/logarithmic
energy budget remains unconditionally sufficient for its stated current
budget. Convergence of full energy to R(N)^2 is not supplied without an additional
vanishing statement. A shell-only bound does not automatically bound the
non-ROOT or nonterminal coefficient energy; the converse direction is open.

**O6 (PROVED outer energy).** The isolated all-one word is the repeated
singleton branch of this fibre decomposition. It has

    E_outer(N,d)=sum_(h in I,3^d h-1 in S) (3[3|h]-1)^2,
    |E_outer(N,d)-(3/2)N*2^(-d)| <=18.                   (13)

There are nine allowed h classes modulo24; the squared weights over a
period sum to18. Counting complete periods plus one cropped period proves
the stated error. No other word occupies these terminal columns (P4), so
E_d>=E_outer. This is a hostile to decay to zero at fixed d, but **not** to
a uniform E_d=O(N) bound. It explains why an energy route can bypass the
absolute-envelope obstruction without claiming that its signs are random.

**Clock/phase repair.** The displayed Fourier recurrence in THM-4667 is
correct. Its branch is P_A(r)=(3r+c_A)/2^A with
c_A=(3+2^a-2^A)/3, A=a+2j. It is not the ordinary branch
S_A(r)=(3r+1)/2^A: at A=3 these are(3r-1)/8 and(3r+1)/8;
at r=0 modulo9 their values are1 and8. The positive-kernel relation already
proved in positive-transport section2 is P=C A, S=A C, hence
P^d=C S^(d-1) A. The boundary operators, translation phases and band register
must be retained before an ordinary Syracuse Fourier estimate is used.
If an actual state is B^(j_i)(s_i), its next odd step has exponent
a(s_i)+2j_i, whereas the canonical homogeneous letter is a(s_i)+2j_(i+1).
Their sums differ by2(j_0-j_d). For example11 ->17 ->13=B(3) uses exponents
1,2 (three shortcut steps), while its canonical letters are1,4 (sum five).
Canonical seed depth and the first-hit integer clock are different registers.

**O7 (PROVED bridge from prior P2).** Source rounding already has an energy
version at **every depth**, without a fresh error invoice at each step.
Let Hhat_N be the rounded row on the same literal seeds from P2 and put
Ehat_d=||L^d Hhat_N||_2^2. Column l1 contraction and ||x||_2<=||x||_1 give

    |sqrt(E_d)-sqrt(Ehat_d)| <=4 floor(log_4(4N)),
    |sqrt(T_d(X))-sqrt(That_d(X))|
                         <=4 sqrt(X) floor(log_4(8X)).   (14)

For N=2^k the stronger exact correction is

    E_d=Ehat_d+2e_k (L^d Hhat_N)(1)+e_k^2,                (15)

where e_k is the known six-periodic ROOT coefficient from P2. Proof:
H_N-Hhat_N=e_k delta_1, and L fixes delta_1. This keeps the literal source,
ROOT port and all transported clocks; it does not compare a discrete vector
with its continuous residue-density model. The unresolved matching problem
is that latter comparison. A scalar verifier checks2870 profiles through
depth40, including both sides of dyadic thresholds, independently of the
previous NumPy source compiler.

This fresh connection changes the board: sibling cancellation pays merged
columns; O1 pays signed values on isolated columns; O3 pays average affine
threshold movement; O4 repairs a finite coverage defect by a genuine new
generator; the geometric escape changes a forbidden signature but still
needs non-Cartesian edges. None of these payments can be substituted for a
different one. In particular(8), a scalar l1-envelope approximation, is
not an l2 matching-error estimate.

## Evidence and next targets

The [manifest](../../04-computation/experiments/outer_golden_20261017h/manifest.json)
contains reproduction commands and hashes. The outer census has289 rows;
the pair computation keeps all multiplicities through depth24 for
translations3^ell,2<=ell<=12, with45056 independent12-bit integer-cylinder
replays. Literal outer paths and1200 three-step controls retain integer
clocks. Small inverse-word enumeration checks the carry generating bound.
Normal and optimized receipts are compared. Proofs of the uniform
quantifiers are above, not extrapolations of those finite checks.

Next: signed survivor correlations for(3), and a merged-jump version of(7)
that pays distinct threshold events rather than every word. The positive
transport already controls norm amplification; the unproved rate is actual
overlap with its marked observer. The geometric companion's signature
obstruction similarly locates a forbidden matching while preserving ports.
Meta-patterns used: **Expose the obstruction first, choose the scale second**;
**Keep observer, recurrence class, and finite head**; **Search the statement
before the method**. No new universal method card is promoted on these two
bounded applications alone.
