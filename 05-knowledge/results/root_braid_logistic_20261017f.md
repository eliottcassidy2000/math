# Three marked current folds and the 2-adic logistic representation

2026-10-17, clock_duality_oct15. **PROVED, scoped, self-audited:** the second
and third current folds, their linear asymptotic envelopes, the exact eighteen-band first-transport tariff, the marked odd
transport identity, and the logistic representation statements below.
**FINITE-EXACT:** native snapshots at k=10,14,18; signed transport at three
specified rows; complete conjugacy tables through 17 bits. **CITED:** the
standard parity-conjugacy theorem. **OPEN:** a uniform sublinear ROOT-current
bound, the dyadic energy target, small-uncertainty entry, and Collatz.
No independent-referee or literature-novelty claim.

## Inheritance and incoming work

Anchor: the uniform ROOT-current estimate. Niche: a contact field's surviving
residue colouring, in [the companion note](joint_contact_field_20261017f.md).
Wildcard: logistic maps with a retained arithmetic source measure.
Board: **primitive source / actual clock / signed transport / contact character /
fixed port / local residue**.

The closest mechanism is [R1--R3 of the preceding note](marked_contacts_root_current_20261017e.md):
an exact odd-source kernel, native ROOT profiles, and refinement brackets.
The canonical hostile f=chi satisfies doubling but has linear current. The
corrected near miss is the proposed nonpositive prefix/rim covariance, refuted
there by native data. The least-used sidecar is the braid index in the inverse
fibre. [The reverse-tree note, Proposition5.1](collatz_mod6_20260917_reverse_tree_pieces.md)
already identifies the primitive seeds and the 4n+1 fibres; its complete
inverse-fibre law and order-three residue cycle are inherited, not new here.
[The procedural-generation bridge, B2](procgen_bridges_20260923_lrc_amm_collatz.md)
already records the 2-adic logistic conjugacy. The new work is its marked
ROOT cycle and its precise use/limitation for this current.

Incoming commit2054a8a273 (S40) was read before integration: its
[full ROOT census](collatz_root_current_exact_20261017.out) covers all rows
through octave2^27, with energy/X^2 between0.2126 and4.5367 in the printed
range. These are finite census statistics, not an independent controller
bound; its floating energy reductions are not promoted to integer receipts.
At k18 the completed census has energy about1.6042 X^2, while the native
stage28 snapshot here is1.634206 X^2: the two have different terminal coverage.
The incoming contact tests agree that mixed edges alone need not constrain
ports. The audited [THM-4665](../../01-canon/theorems/THM-4665-ei17-has-no-faithful-realization-over-the-split-and-conjugation-fixed-primes-of-q-zeta55-tested-so-every-q-zeta55-realization-has-bad-reduction-there.md)
and MISTAKE-625 locate the finite-field implication in conjugation and good
reduction, not in a low count's apparent rarity. This is used in the companion
field proof, without importing a finite embedding exclusion into colouring.

The later fetch through309af1b454 adds a full-current replay of the inherited
I/B fold through octave2^26, with zero row failures. Its covariances again
change sign (about+0.752078 X^2 at k22 and-5.091197 X^2 at k23).
[THM-4666](../../01-canon/theorems/THM-4666-class-numbers-are-ranking-excesses-of-kronecker-tournaments-and-local-optimality-of-the-natural-paley-order-excludes-a-real-zero.md)
uses partial sums of the Legendre walk to express ranking changes and a
Fekete polynomial. Section3c uses the same exact Abel-summation operation on
our three-lift coefficient fibres; no ranking-optimality or L-function-zero
conclusion is transported. HYP-9287 remains a conjecture, and its local-search
"unbeaten" cases are not used as proved optimality.

## 1. A second exact fold, and a sharp hostile

Use shortcut T(n)=n/2 for even n and (3n+1)/2 for odd n. Let
w(n)=3[n=2 mod3]-1, and chi(n)=0,1,-1 on residues0,1,2. For a profile f,

    R_f(N)=sum_(N<n<=4N) w(n)f(n).

Assume f(2n)=-f(n) and, for odd n, f(Bn)=f(n), Bn=4n+1. Both hold for
the actual ROOT function: (-1)^tau on paths reaching stopped ROOT1, zero
on paths that never do. This definition does not assume convergence.
The audited native profiles also satisfy both relations within their source
universe. No ROOT value is assigned to an unobserved native path.

Every positive odd n has a unique representation B^j(s), with
S={positive odd s:s !=5 mod8}. The inherited first fold gives

\[
 R_f(N)=\sum_{s\in S}H_N(s)f(s),\qquad
 H_N(s)=\sum_{j:B^j(s)\le4N}
 \begin{cases}-3\chi(B^j(s)),&B^j(s)\le2N,\\
 w(B^j(s)),&2N<B^j(s)\le4N.
 \end{cases}                                                   \tag{1}
\]

**R4 (PROVED).** For all N and s,

| s modulo3 | possible H_N(s) |
|---|---|
|0 or1|0,-1,-3|
|2|0,2,3|

Indeed, B advances residues by1 modulo3, every complete triple of prefix
coefficients sums to zero, and the rim contains at most one B-child. The
three possible prefix remainders and the optional rim term give the table.
In particular each column has the fixed sign sigma(s)=+1 on s=2 mod3 and
-1 otherwise. One fixed profile, f_*(B^j s)=sigma(s) on odd numbers and
f_*(2^a n)=(-1)^a f_*(n), attains every row's l1 norm:

\[
 A_0(N)=\sum_s|H_N(s)|=R_{f_*}(N)
       ={16\over7}N+O((\log(N+2))^2).                         \tag{2}
\]

Proof of the constant: primitive seeds in each residue r modulo3 occupy
three classes modulo24, density1/8. Partition n uniquely as 2^a B^j(s).
The mean of w(n)f_*(n) is

\[
 {1\over8}\sum_{r=0}^2\sigma(r)
 \sum_{a,j\ge0}{(-1)^a\over2^a4^j}
 w\bigl(2^a(r+j)\bigr)={16\over21}.                           \tag{3}
\]

Split a modulo2 and j modulo3; the geometric denominators are 1-1/4 and
1-1/64. Each progression contributes its density times the cutoff plus O(1).
There are O(log^2 N) relevant (a,j); the omitted geometric tail is O(log N).
Subtract the prefix at N from the prefix at4N to obtain (2).

Changing the whole seed1 fibre to f=+1 gives f_*(1)=1 and changes each
current by at most6. Thus normalization, doubling AND braid invariance still
permit a fixed ROOT-normalized profile with cubic dyadic energy. The missing
condition is the actual odd transition. This is stronger than the previous
hostile, but does not refute the arithmetic current conjecture.

## 2. A paid third fold: the odd transition removes 10N/21

For s=1 mod8, put t=2s+1. The accelerated odd map U(n) is the odd part of
3n+1. Then

\[
 \nu_2(3s+1)=2,\quad \nu_2(3t+1)=1,\quad U(t)=B(U(s)).          \tag{4}
\]

Hence any odd-compatible profile obeys f(t)=-f(s). For native finite
controllers this equality is available only after the relevant paths have
been observed or the exact coalescence certificate has been attached.
It must not silently be applied to arbitrary partial native profiles.
The full ROOT function satisfies it without assuming universal entry.

The pairs (s,2s+1), s=1 mod8, are disjoint: the first coordinates are1 mod8
and the second3 mod16. Replace their two columns in (1) by

    (H_N(s)-H_N(2s+1)) f(s),

leaving all unpaired primitive seeds literal. Since sigma(2s+1)=sigma(s),
this removes exactly

\[
 C(N)=2\sum_{s=1\bmod8}\min(|H_N(s)|,|H_N(2s+1)|)             \tag{5}
\]

from the unsigned envelope. This is cancellation forced by an odd clock,
not a random-sign hypothesis.

**R5 (PROVED).** Uniformly as N grows,

\[
 C(N)={10\over21}N+O(\log(N+2)),\qquad
 A_1(N):=A_0(N)-C(N)={38\over21}N+O((\log(N+2))^2),            \tag{6}
\]

and |R_f(N)|<=A_1(N) for every |f|<=1 satisfying the stated relations.

Here is the complete tariff calculation. Ignore the bounded offsets in
B^j(s)=4^j s+(4^j-1)/3, and write x=s/N. On the band
2/2^b<x<=4/2^b, b>=1, the three values of the minimum in (5), for
s=0,1,2 mod3, repeat every six bands as

    (0,1,2), (1,0,2), (1,0,2), (0,1,2), (0,0,0), (0,0,0).

Each relevant residue class has density1/24 within s=1 mod8. Multiplying
by the factor2 in (5) and the band length gives

\[
 \sum_{b\ge1}{\text{row sum at }b\over6\,2^b}
 ={\tfrac14(1+\tfrac12+\tfrac14+\tfrac18)\over1-1/64}
 ={10\over21}.
\]

The actual thresholds differ from their ideal thresholds by at most a
constant: at most1/3 for H(s), at most2/3 for H(2s+1). There are O(log N)
thresholds, each coefficient is bounded by3, and counting each residue
progression costs O(1). The remaining tiny bands contribute O(1). This proves
the error in C; (2) proves the error in A_1. Equality in the finite envelope
can be attained if only the listed disjoint pair relations are imposed;
further odd relations can improve it.

At N=262144 the exact envelopes are A_0=599188 and A_1=474356; exactly
124832 units cancel. The scaled errors 21C-10N and 21A_1-38N are32 and4.
The proof is a uniform linear improvement, not the square-root scale or a
subcubic energy estimate. The six-band pattern reflects the parity/rim
phase together with the order-three braid phase; it is not an identification
with an unrelated mod6 phenomenon.

## 3. Keep the integer clock after quotienting the braid

Write an odd state as n=B^j(s). Set a(s)=nu2(3s+1), which is1 or2 on S,
and write U(s)=B^ell(P(s)) with P(s) in S. The exact marked transition is

\[
 (s,j)\longmapsto(P(s),\ell),\qquad
 \text{shortcut clock increment}=a(s)+2j.                    \tag{7}
\]

This follows from 3B^j(s)+1=4^j(3s+1). Forgetting j preserves clock parity
but loses integer time. The computation retains every source, its actual
odd endpoint and its accumulated integer clock. Equal braid fibres can
then combine signed weights, while their endpoint/clock arrays remain as
authentication sidecars. Seed1 fibres are stopped only after appending their
actual even path of length2+2j (zero if already at1).

For a fixed N push the initial H_N(s) along these observed transitions.
Let b_d be the coefficient at ROOT and D_d the l1 mass at all other primitive
terminals. The exact interval

    b_d-D_d <= R_f(N) <= b_d+D_d

holds for every compatible bounded f with f(1)=1. Signed merging can only
reduce |b_d|+D_d. This is a prospective controller with explicit observed
odd paths, distinct from the inherited first-reset policy.

| N | initial l1 | D_1 | D_12 | D_28 | D_64 | b_64 |
|---:|---:|---:|---:|---:|---:|---:|
|1024|2340|1827|510|196|15|47|
|16384|37448|29232|8869|3023|408|16|
|262144|599188|467688|143552|52456|7344|-82|

These are three rows, not an octave theorem. The exact one-step transport is
slightly stronger than the disjoint-pair tariff because additional fibres
merge. A next sufficient target is a source-weighted contraction of the
transported rows that is stable under repeated transport. An arbitrary l1
operator contraction cannot work: a singleton unmerged terminal retains its
mass. The special row geometry, opposite-sign overlap, and ROOT extraction
must enter together. This restates the obstruction as a concrete row-class
problem rather than discarding the odd relation.

The inherited native tests are retained separately. At k18 and stages1,4,12,28,
max_N |R|/C_cert is respectively0.000998734,0.001409867,0.001671496,
0.001640304. Here C_cert=sum_s |H_N(s)| |f(s)|; (C_cert+-R)/2 are exact
positive/negative masses. All endpoints and clocks are authenticated, and
normal versus optimized runs agree byte-for-byte. Neither sign balance nor
ROOT energy improves monotonically under refinement.

## 3b. Close the exact first-transport tariff: an eighteen-band law

The full first transport can be computed analytically, including the
additional cancellations beyond the disjoint pairs in section2.
For primitive t, all primitive preimages of P are

    s_j=(2^(a_j) B^j(t)-1)/3,
    a_j=2 if B^j(t)=1 mod3, a_j=1 if B^j(t)=2 mod3;

omit the j with B^j(t)=0 mod3. Its signed coefficient is exactly

\[
 C_N(t)=\sum_j(-1)^{a_j} H_N(s_j).                            \tag{10}
\]

Let h_r(x) be the ideal version of H_N(s) with s/N=x and s=r mod3:
sum the prefix weight -3chi(r+k) or rim weight w(r+k) as 4^k x runs
through (0,2] and (2,4]. It satisfies

    h_r(x/64)=h_r(x) for 0<x<=4,    sum_r h_r(x)=0.

For r=t mod9, let rho_j=B^j(r) mod9, with the same a_j rule, and
r_j=(2^(a_j)rho_j-1)/3 mod3. Define

\[
 c_r(x)=\sum_j(-1)^{a_j}
 h_{r_j}\bigl(2^{a_j}4^j x/3\bigr).                          \tag{11}
\]

It is constant on 3/2^b<x<=6/2^b, b>=0. Its nine coefficients have
**exact period18 in b**. To prove this, split the sum for c_r(x/4^9)
into j>=9 and j<9. The first part is c_r(x), since B has order9 modulo9.
Group the remaining nine terms by j modulo3. In each group a_j is fixed,
the three residues r_j run through0,1,2, and their arguments differ by
factors64. All arguments are at most2 when x<=6. Scale invariance of h and
its zero sum make each group vanish. This proves the period for every
band, not just the finite table.

The sums of absolute coefficients in its eighteen bands are

    4,10,8,14,12,12,12,12,16,16,12,12,14,8,10,4,0,0.

The receipt stores all nine signed entries in each row. Primitive t in each
class modulo9 has density1/24. Multiplying the band lengths by those
densities gives the exact uniform asymptotic:

**R7 (PROVED).**

\[
 A_2(N):=\sum_{t\in S}|C_N(t)|
 ={155896\over87381}N+O((\log(N+2))^2)
 \approx1.78409494055N,\qquad |R_f(N)|\le A_2(N).              \tag{12}
\]

Indeed the coefficient is

    (1/8) sum_(b=0)^17 2^(-b) sum_r |c_r(9/2^(b+1))|
    divided by (1-2^(-18)) = 155896/87381.

For the error bound, each actual source B^k(s_j) equals

    (2^(a_j) 4^(j+k)/3) t + 2^(a_j)4^k(4^j-1)/9 -1/3.

After dividing by its slope, the affine offset has absolute value at most1.
Consequently each actual prefix/rim threshold in t is within1 of its ideal
threshold. The O(log^2 N) relevant (j,k) pairs therefore change the sum of
absolute coefficients by O(log^2 N). Counting ideal primitive progressions
costs O(log N); the remaining tiny bands cost O(1). This supplies the stated
uniform error independently of any ROOT census.

The numbers18 and19 now have an explicit relationship in this calculation:
87381=(2^18-1)/3=3^2*7*19*73. The eighteen binary bands come from the affine
braid's order9 modulo9 and its factor4 in size. This is an arithmetic
explanation for the denominator; it does not identify this kernel with the
separate finite-field generator problem where19 appeared earlier.

At N=262144, the complete first-transport l1 mass is467690, within1 of the
leading term in (12). The inverse-fibre formula (10) and an independent
forward push agree on74240 coefficients for N<=256. This is an improvement
on (6), still linear and still insufficient for the dyadic energy goal.
An iteration theorem requires a new family of row norms that survives the
actual transfer. It cannot be obtained by raising this one-row tariff ratio
to a power without proving closure of that family.

## 3c. The finite-depth recurrence closes, with an explicit conductor cost

The preceding tariff extends to every fixed number d of odd transports.
This is an explicit consumer of the inherited ternary-memory principle in
[THM-4658](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md),
not a new infinite-depth spectral theorem.
Let q_d=3^(d+1), P_d=2q_d, and M_d=4(3/2)^d. On the band
M_d/2^(b+1)<t/N<=M_d/2^b, define the ideal coefficient K_d(r,b), r mod q_d.
Set K_d(r,b)=0 when b<0 and K_0(r,b)=h_r(3/2^b) for b>=0. Then

\[
 K_d(r,b)=\sum_{j\ge0:\,B^j(r)\not\equiv0\ (3)}
 (-1)^{a_j}
 K_{d-1}\!\left({2^{a_j}B^j(r)-1\over3}\bmod q_{d-1},
                 b-(a_j-1)-2j\right).                       \tag{13}
\]

Only nonnegative second arguments contribute, so this is a finite exact
integer recursion. The source coordinate and the signed clock both occur.

**R8 (PROVED).** For every d>=0 and b>=0,

1. K_d(r,b+P_d)=K_d(r,b).
2. The three last-digit lifts sum to zero:
   sum_(i=0)^2 K_d(r+i*3^d,b)=0.
3. If A^(d)(N) is the l1 mass after d complete odd transports, then

\[
 A^{(d)}(N)=\gamma_d N+O_d((\log(N+2))^{d+1}),\quad
 \gamma_d={\displaystyle\sum_{b=0}^{P_d-1}2^{-b}
                  \sum_{r=0}^{q_d-1}|K_d(r,b)|
             \over 4\,2^d(1-2^{-P_d})}.                     \tag{14}
\]

Proof. The base case is the three-term prefix cancellation. For the
last-digit property, lifting r by3^d leaves every a_j unchanged and moves
the predecessor through the three lifts modulo3^d; use the preceding
last-digit property term by term in (13).
For the period, write the sum at b+2q_d as j>=q_d plus j<q_d. The first
part is the sum at b, since the affine map B has order q_d modulo q_d.
Group the extra terms as j,j+q_(d-1),j+2q_(d-1). Their second arguments
differ by the old period; their first arguments run through the three
last-digit lifts. All those second arguments are nonnegative. Thus they
sum to zero. The elementary fact used here is
v3(4^(3^d)-1)=d+1, obtained by cubing1+3k and induction; B^(3^d) is a
nonzero translation by a multiple of3^d modulo3^(d+1).

For (14), each primitive residue class modulo q_d has density
3/(8q_d)=1/(8*3^d). Multiply this by band length M_d/2^(b+1), then sum
one period and its geometric copies. Composing d inverse fibres with the
initial B-fibre gives O_d(log^(d+1) N) relevant exponent tuples. Each actual
affine threshold differs from its homogeneous counterpart by O_d(1): an
inverse-fibre slope is at least2/3, and its offset divided by slope is
bounded independently of the fibre index. Summing these threshold errors
proves the stated fixed-depth remainder. This argument also improves the
coarser logarithmic-square remainder in (2) when d=0.

|odd transports d|ternary modulus|binary-band period|gamma_d|
|---:|---:|---:|---|
|0|3|6|16/7 =2.2857142857|
|1|9|18|155896/87381 =1.7840949405|
|2|27|54|26100667631514128/18014398509481983 =1.4488781081|
|3|81|162|1.2280701518 (exact rational in receipt)|

Every signed table entry, last-digit sum and an independently unwrapped
second period was checked. The finite forward masses at N=262144 are
599188,467690,379820,321954. The last is about22.78 above its leading term,
so the period recurrence is not being misread as an exact no-error formula.

**Conductor factor, and the incoming tournament connection.** Put
P_(d,b)(z)=sum_r K_d(r,b) z^r and q=3^d. The exact last-digit balance gives

    P_(d,b)(z)=(1-z^q) Q_(d,b)(z),
    Q_(d,b)(z)=sum_(r<q) z^r [K_d(r,b)
                     +(K_d(r,b)+K_d(r+q,b))z^q].

This is the three-lift Abel transform, analogous to the full Legendre-walk
Abel transform in incoming THM-4666. Every additive Fourier mode that factors
through modulus3^d is annihilated. The surviving modes have exact conductor
3^(d+1). The identity was reconstructed coefficient by coefficient in every
computed band. This explains why a coarse residue average can look perfectly
balanced while a source-specific discrepancy persists: the marked source
weights need not be constant on those three lifts. The map preserves the
coefficient polynomial, but a bound after contracting it with ROOT signs
still needs their within-fibre differences. No intrinsic complete binary
orientation is supplied here, so these coefficient fibres are not declared
a tournament merely to mimic the incoming result.

This produces a positive family of progressively smaller envelopes and an
exact place to test new controls. It also quantifies the obstruction:
each odd transport triples the residue modulus and the band period, and
the remainder constant depends on depth. Neither gamma_d->0 nor a suitable
rate has been proved. Even a decay theorem for gamma_d alone would not
license d growing with log N without a corresponding uniform remainder
estimate. That coupled depth/source-scale estimate is the next critical
mathematical target.

## 4. Logistic maps: exact conjugacy, singular source measure

Use the **2-adic** logistic map L(x)=x(1-x)/2 on Z_2, not the real map
4x(1-x). If x,y have the same parity, 1-x-y is odd, so

    v2(L(x)-L(y))=v2(x-y)-1.

Each parity ball maps bijectively onto Z_2. In branch coordinates the two
maps are V0(y)=y(1-2y) and V1(y)=-y(2y+1); each is a bijective isometry by
the odd difference factor and induction on binary digits. This is an instance
of [Bernstein--Lagarias, The 3x+1 Conjugacy Map, Appendix B, Theorem B.1](https://websites.umich.edu/~lagarias/doc/bernstein.pdf).
For either map F, its parity itinerary Q_F is an isometry onto binary
sequences. Thus H=Q_L^(-1) Q_T is the unique itinerary-preserving isometry
with L H=H T. This is inherited standard structure, not a Collatz proof.

**R6 (PROVED).** The T-cycle1<->2 becomes the L-cycle

\[
 x^2-3x+6=0,\qquad x=(3\pm\sqrt{-15})/2.                      \tag{8}
\]

Indeed (L^2(x)-x)=(L(x)-x)(x^2-3x+6)/4; H(1) is its odd root and H(2)
its even root. Also H(0)=0, H(-1)=-1, H(-1/3)=1. The inverse branches of L
are the even root ell_0(y) of L(x)=y and ell_1(y)=1-ell_0(y).
On odd inputs the braid becomes

\[
 H(Bn)=\ell_1\ell_0\ell_0 L(H(n)),\qquad
 Q_T(Bn)=4Q_T(n)-3.                                           \tag{9}
\]

Thus B inserts two even clocks after the first odd clock. Its conjugate
contracts distances by1/4 and fixes1, the logistic preimage of0. This is the
exact counterpart of the 2-adic limit B^j(n)->-1/3. It does not say that the
ROOT signs converge continuously at that point.

The logistic basin of (8) is countable, dense and Haar-null: its itineraries
are exactly the eventually alternating binary words. Consequently the
logistic ROOT-sign function is zero Haar-almost-everywhere, although the
marked arithmetic sources can carry nonzero signs. H sends complete binary
residue systems to complete binary residue systems; this gives perfect
coarse uniformity and still says nothing about membership in that null basin.
The needed measure is the weighted pushforward of the original integer
sources, including their clocks, not unrestricted Haar measure. This is a
2-adic statement; it is not the 3-adic Wiener law in
[THM-4658](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md).

The same discriminant-15 occurs in the finite unit-distance pentagon's CM
form in [THM-4663](../../01-canon/theorems/THM-4663-polygon-ladder-of-the-finite-unit-distance-graphs-closed-walks-are-kloosterman-moments-pentagons-level-15-hexagons-level-6.md).
The source here is a quadratic dynamical cycle, the target there an elliptic
modular trace. The common imaginary quadratic field is exact; a map preserving
integer-source current or geometric port constraints has not been supplied.
That missing map prevents using the coincidence as a bound.

## Reproduction and next tests

All receipts and commands are in the
[manifest](../../04-computation/experiments/joint_logistic_current_20261017f/manifest.json).
Controls include 256 complete small current rows, 8320 odd-pair checks, exact
rational geometric sums, all parity words through17 bits, and92422 actual
shortcut-clock steps from2048 positive control sources. The finite source
clocks are checked independently of parity-table conjugacy. No finite census
is substituted for the unproved uniform statement.

Next tests: bound the rate of decay of gamma_d jointly with the growing
ternary conductor and affine threshold error; find a cone of marked rows on which a paired
transport bound survives refinement; and estimate the integer-source measure
of logistic cylinders relative to the particular ROOT basin, retaining the
fine digits. A geometric residue quotient and a dynamical parity quotient
both preserve a useful relation while potentially losing the consumer's
marked data; their useful common mechanism is to make that loss explicit.
