# Large Collatz sources as finite radix fractions, with exact block transport

**Status: PROVED elementary representation and block-transport identities;
FINITE-EXACT controls; universal Collatz coverage and unbounded ROOT energy
remain OPEN.** Self-audit using algebra, literal trajectories and independent
modular readers; no independent-agent review or historical-priority claim.

2026-10-18. The owner's proposal is to regard a large integer as the digits
of a real number between zero and one, in any base. This produces useful
proof compression when the radix, scale, carry and actual endpoint survive.
It does not make the leading digits determine the low-digit Collatz branch.

Artifacts: [checker](../../04-computation/experiments/collatz_radix_families_20261018p.py)
and [exact receipts](collatz_radix_families_20261018p.json).

## 1. Inheritance, and what this session adds

Anchor: explicit families of arbitrarily large certified integers. Niche:
compare radix choices by the terminal value they actually preserve. Wildcard:
the marked binary-expansion endpoint in the parity-number dictionary.

The six-item board is source digits / scale / low-digit guard / affine carry /
terminal odd factor / integer clock.

- Closest proved mechanism: the fixed-head terminal compiler in
  [charged terminal families](collatz_charged_terminal_family_20261012a.md).
  In particular its family `(2^(18j+11)-23)/27 -> 1` is **inherited**, not a
  newly discovered family. We identify its exact positional language, give
  a general arbitrary-base transport with its unpaid odd factor, and provide
  a very large compressed certificate and decoder.
- Canonical hostile: [crossing-phase audit, section 4](crossroads_crossing_20260926_phase_audit.md)
  proves that every finite allowed leading-phase word coexists with
  arbitrarily long pure-growth runs. A leading-digit pattern cannot force
  an extra halving. This is stronger than our two-source illustration below.
- Corrected near misses: [projective receipts, section 4](collatz_projective_receipts_20261009.md)
  retain the native lattice after real normalization; their integer midpoint
  `(445,351,668)` has the wrong valuation. The new endpoint repair in section
  7 concerns a different loss, the two binary expansions of a dyadic rational.
- Least-used sidecar: [golden digit/carry](golden_digit_carry_20261003.md)
  keeps radix position and an exact carry. The new radix shift by `floor(C/p)`
  below serves the same kind of purpose, without identifying the two rings.

Incoming commit `987f76f6cf` was read and integrated before deriving these
rules. Its [marked Gray-map calculation](gray_selector_recursion_20261018e.md),
section 4, supplies an especially sharp analogy: dropping a carry preserves
the entire code and every chart but sends source `(1,1,0)` to the wrong word.
Here dropping an odd factor preserves a plausible repeating digit picture
but changes the actual terminal source. Neither unmarked count detects it.

The later incoming commits `87ca2c2843` and `696cd48bf0` classify
[105 optimal atlases](hamming_matching_atlases_20261018e.md) and retain the
central phase lost by the slope quotient. Their ordered-operation identity
prompted section 5A's explicit test of changing radix generators. The geometric
three-colour-list gate is constructed; forcing those lists remains open.
Our arithmetic generators below have a different, noncentral defect.

We use the existing META-PATTERNS moves of retaining the native operation
and the coordinate lost at the first collision. No new navigation card is
needed. Neither chart counts nor a digit-language count supplies ROOT energy.

## 2. The faithful real representation and its two ends

For base `b>=2`, let an integer `n>0` have length `L`, and write

\[
x=n/b^L=0.d_{L-1}d_{L-2}\cdots d_0\quad\text{in base }b.
\]

The pair `(x,L)` recovers `n=b^L x`. For general integers, x alone loses
trailing zeroes: binary `1` and `10` both normalize to `0.1`. On **odd binary
integers**, the exact reduced denominator already recovers L. The reason
to retain L explicitly there is finite precision, limits and exact clocks,
not a false assertion that exact normalized odd binary fractions collide.

For the odd map

\[
U(n)=(3n+1)/2^{v_2(3n+1)},
\]

put `a=v2(3n+1)` and let `L'` be the next length. Then exactly

\[
x'=\frac{3b^{L-L'}}{2^a}x+\frac{1}{2^a b^{L'}}.       \tag{2.1}
\]

The additive term is small at large scale. However a is determined by
trailing binary digits. For every `K>=3`, the odd sources `2^K+1` and
`2^K+3` have normalized binary values only `2^-K` apart, yet their first
valuations are respectively 2 and 1: one decreases, the other increases.
In every nonempty open interval of leading mantissas, sufficiently large
bands contain every fixed trailing residue class, because its normalized
mesh tends to zero. Real closeness does not couple the two ends.

A useful finite state is therefore an **exact source with a leading interval,
a trailing residue and its length**. The interval and residue are checked
views of the same integer, not independently chosen replacements for it.
Reversing source digits makes low digits visible first but does not preserve
ordinary size order. Reading **orbit parities** instead is yet another map:
it requires trajectory information and is not the source's digit reversal.

Compact limits add a second hazard. The positive integers `2^K-1` tend to
the real normalized endpoint 1 and to the negative integer -1 in the
2-adic metric. Each finite string is positive; its infinite trailing-digit
limit is a different arithmetic source. This is the existing compactness
warning in [carry/descent/transport](carry_descent_transport_branches_20261018l.md),
section 5, not a new convergence argument.

## 3. General radix-block transport, including its exact odd-factor cost

Suppose a positive odd source n follows a supplied valuation word of length k>=1,
`w=(a1,...,ak)`, with positive odd endpoint m. Set

\[
A=\sum_{i=1}^k a_i,\quad p=3^k,\quad
C=\sum_{i=1}^k3^{k-i}2^{a_1+\cdots+a_{i-1}}.
\]

Then

\[
pn+C=2^A m.                                            \tag{3.1}
\]

The integer C depends on the first `k-1` valuations, **not the last one**;
`3` does not divide C. Choose an integer `B>1` with `B=1 mod p`, and write
`B=2^e u`, u odd. Define

\[
P_B(n)=Bn+\frac{C(B-1)}p,\qquad
n_j=P_B^j(n)=\frac{B^j(pn+C)-C}{p}.                   \tag{3.2}
\]

**Proposition 1 (PROVED).** Every `n_j` is a positive odd integer with the
actual word and endpoint

\[
\boxed{n_j\xrightarrow{(a_1,\ldots,a_{k-1},a_k+je)}u^j m.} \tag{3.3}
\]

Thus the exact halving clock is `A+je`, and the ordinary clock that counts
each `3n+1` and each halving separately is `A+je+k`. This is a statement about
the first k odd edges. If the endpoint is not 1, its later route remains an
additional obligation.

**Proof.** Integrality follows from `p|(B-1)`. Also

\[
P_B(n)-n=\frac{(B-1)(pn+C)}p
        =\frac{(B-1)2^A m}{p}
\]

is divisible by `2^A`, since p is odd. Every guard for the first `k-1`
odd edges is determined modulo `2^(a1+...+a_(k-1)+1)`, which divides
`2^A`. Those edges are therefore unchanged. After them, the last numerator
is multiplied by B, so its exact valuation increases by e and its odd part
is multiplied by u. Equivalently `pP_B(n)+C=B(pn+C)`. Positivity is immediate.
Repeat, since C and the earlier guards remain unchanged. QED.

To interpret P as appending positional digits, choose `B=b^ell`. Put
`C=pq+r`, `0<r<p`, and retain the shifted coordinate `y=n+q`. Then

\[
P_B(n)+q=B(n+q)+d,\qquad d=\frac{r(B-1)}p,
\quad 0<d<B.                                         \tag{3.4}
\]

This appends the fixed length-ell base-b block for d to **y**. It appends
directly to n precisely when q=0; otherwise forgetting q is a carry error.
For any base coprime to 3 there is an eligible ell. A base divisible by 3
cannot satisfy `b^ell=1 mod3^k`; this excludes this particular unchanged-C
construction, not other ternary representations or changed heads.

**Exact base distinction.** In this construction the same odd endpoint is
preserved for every repeat iff `u=1`, equivalently B is a power of two.
Hence base `2^h` permits it; a base with an odd prime factor leaves the
explicit factor `u^j`. For binary, the least eligible block length is

\[
\operatorname{ord}_{3^k}(2)=2\,3^{k-1}.
\]

For completeness, 4=1+3 and cubing a number `1+3^s v`, `3∤v`, raises the
valuation of its difference from 1 from s to s+1. Factoring a general
exponent into a power of 3 and a 3-unit gives
`v3(4^t-1)=1+v3(t)`. Any power of 2 equal to 1 modulo 3 has even exponent;
the displayed order follows. No random-digit hypothesis is involved.

If a binary pump begins with a supplied strict first-hit certificate to 1,
the pumped sources inherit it, with the longer final halving. Intermediate
head values increase with the source, so no new earlier visit to 1 appears.
The theorem also allows endpoints other than 1 and non-first-hit words;
first-hit status is not silently attached to those inputs.

## 4. A source with eighteen trillion bits, certified symbolically

The inherited seed certificate is

\[
75\xrightarrow{1}113\xrightarrow{2}85\xrightarrow{8}1.
\]

Here `p=27`, `C=23`, `A=11`. Choose `B=2^18=262144`. Then

\[
d=23(2^{18}-1)/27=223307
  =(110110100001001011)_2.
\]

Consequently the exact regular language

\[
\boxed{\quad (1001011\;(110110100001001011)^j)_2,
                    \qquad j\ge0,\quad}               \tag{4.1}
\]

consists of the sources

\[
N_j=\frac{2^{18j+11}-23}{27}
  \xrightarrow{(1,2,18j+8)}1.                          \tag{4.2}
\]

The first sources are 75, 19,884,107, 5,212,499,568,715, and
1,366,425,486,941,448,267. Each has exactly `18j+7` binary digits. Their
finite radix fractions have the exact limit and error

\[
\frac{N_j}{2^{18j+7}}=\frac{16}{27}
             -\frac{23}{27\,2^{18j+7}}.               \tag{4.3}
\]

This is the proposed “large integer as a real fraction” idea with a full
proof attached to each finite truncation. The limit itself is not doing
the proving: (3.1), the guards and the integer clock are.

For `j=10^12`, this certifies an explicitly specified integer with
**18,000,000,000,007 binary digits**, without writing them all down. Its
odd-map certificate has three edges, its shortcut/halving clock is
18,000,000,000,011, and its ordinary clock is 18,000,000,000,014.
The compressed input is the formula or repeat count; this is not a claim
to read an arbitrary eighteen-trillion-bit input in constant time.

The checker separately materializes `j=1,000,000`, an **18,000,007-bit**
integer, and literally replays its three odd edges. The trillion-repeat
source is not materialized. Two independent modular readers—full-denominator
exponentiation modulo `27q`, and powering the affine recurrence modulo q—
agree at all eleven stated moduli, including powers of 2 and 3.

Membership also has an exact source test: a positive integer belongs to
(4.1) iff `27n+23` is a power of two. Indeed its exponent must be
`11 mod18`, since `2^11=23 mod27` and the order is 18. Positivity excludes
earlier exponents. This decoder tests the supplied integer; it does not
replace an unrecognized integer by a convenient member of the family.

## 5. The decimal analogue exposes the remaining terminal problem

For the same head and seed take `B=10^3=1000`, with

\[
d=23(1000-1)/27=851,\qquad B=2^3\,125.
\]

The decimal integer with initial digits 75 followed by j copies of 851 is

\[
D_j=(75\;(851)^j)_{10}
 =\frac{2048\,1000^j-23}{27},\qquad
\boxed{D_j\xrightarrow{(1,2,8+3j)}125^j.}             \tag{5.1}
\]

For example `75851 -> 113777 -> 85333 -> 125`. Inserting the decimal point
gives `D_j/10^(2+3j) -> 512/675`. Both this and (4.3) look like convergence
of repeating fractions, but only the binary construction leaves terminal
value 1. Equation (5.1) does not ground the entire family `125^j`; that
requires another proof. It is a concrete failure of the inference
“a stable repeating digit pattern preserves ROOT.”

The strongest survivor is the exact reduction of one huge family to a
different, smaller family, with all clocks and the odd factor retained.
A subsequent controller may legitimately discharge that factor. Simply
discarding it cannot.

### 5A. Same-anchor generators commute; changing anchors changes the receipt

Write the affine map as `P_(B,c)(x)=Bx+(B-1)c`, where c=C/p is the retained
rational offset. Direct composition gives

\[
P_{B,c}(P_{D,d}(x))-P_{D,d}(P_{B,c}(x))
       =(B-1)(D-1)(d-c).                              \tag{5.2}
\]

For B,D>1, these maps commute iff c=d. With the same marked head, both
compositions retain that head, add their exact halving clocks and multiply
their odd endpoint factors. In particular the binary and decimal pumps
above commute as maps even though the decimal one leaves a terminal cost.

Changing the head requires more than correcting the final translation:
the next head's guard must hold at the intermediate source. For example,

\[
P(x)=4x+1,\qquad Q(x)=262144x+223307.
\]

P is the one-edge sibling pump with c=1/3. Q is the three-edge pump with
c=23/27 and earlier head `(1,2)`, whose source guard is `11 mod16`. At 75,

\[
\begin{aligned}
P(Q(75))&=79,536,429\xrightarrow{(3,2,26)}1,\\
Q(P(75))&=79,128,651\xrightarrow{(1,2,16)}4075.
\end{aligned}
\]

Their difference is 407,778, exactly (5.2). In the second order, the first
failed transfer is already at `P(75)=301`, which does not satisfy Q's
head guard. Q still defines an integer affine value, but its claimed
endpoint-preserving interpretation is unavailable there.

This is the typed connection to the incoming atlas operators: changing an
operation order can be invisible in a quotient and visible in the marked
output. Their defect is a central binary phase. Our translations are **not
central** under affine dilation: conjugation scales a translation by B.
No identification of the two operator groups is made. The useful shared
test is to compute the ordered action and inspect the first failed guard.

## 6. A local deletion rule and the universal-coverage boundary

When B is a power of two, (3.2) is reversible on its marked arithmetic
domain. If `n=P_B(z)` with a positive odd z following the supplied head,
then n and z have the same endpoint after k odd edges, and `z<n`.
One may delete the block in the shifted coordinate, retain q and check
the head on z. This is a source-specific, strictly smaller dependency,
not a statistical assertion about digit frequencies.

At k=1 it recovers the familiar sibling rule `n=4z+1`; at the 75 head
it deletes the length-18 block above. Different heads have different C,
shifts, guards and block lengths. Their union must retain which head was
actually legal. A program with unbounded integer registers is not a
finite-state recognizer merely because its instruction set is finite.

The existing [guarded-pumping theorem](collatz_guarded_pumping_memory_20261004.md),
sections 3–4, proves a stronger relevant limitation: a regular language
of literal first-hit inverse D/R certificates has bounded odd depth and
only polylogarithmic source coverage below a height bound. Our fixed-three-
odd-step language has only O(log X) sources up to X, entirely consistent
with that result. Huge values alone do not fill the uncovered sources.

**Next precise helper target (OPEN).** Construct and prove coverage of a
source-digit reduction system that, for each positive odd input other than
1, authenticates either an actual smaller forward endpoint or a smaller
same-future block-deletion dependency. Each rule must carry its original
source, actual head and integer clock. Strong induction would then prove
ROOT, because each dependency is strictly smaller and genuinely transfers
the conclusion. The missing assertion is universal acceptance by those
rules, not the induction once acceptance is supplied.

This target can use changing heads and exact carry registers. The older
all-growth and cheap-refill examples rule out deriving it from a fixed
leading-digit window or a linear height-plus-run potential. We have not
proved coverage, a new uniform ROOT-current bound, or a bound on unpaid
tail interactions. The previous all-scale clock-through-16 theorem in
[actual ROOT clocks](root_actual_clock_certificate_20261018o.md) remains
the stronger established energy result; these radix certificates do not
extend its clock-uniform conclusion.

## 7. Endpoint repair in the earlier real parity dictionary

[THM-4634, parity-number dictionary](../../01-canon/theorems/THM-4634-collatz-parity-number-dictionary-irrationality-exponents-of-parity-vectors-and-jarnik-strata.md)
reads orbit parity bits as `rho(x)=sum epsilon_i 2^(-i-1)`. Its interval
doubling statement omitted the dyadic endpoint distinction. The exact law is

\[
\rho(Tx)=2\rho(x)-\epsilon_0.                          \tag{7.1}
\]

At x=-1, rho=1 and the next rho is 1, while the fractional part of 2 is
0. At x=-2, rho=1/2 and the next rho is again 1, whereas the fractional
part of 1 is 0. These are precisely the two exceptions to the literal
interval fractional-part formula: failure occurs iff the next parity
tail is all ones, equivalently `Tx=-1`, whose preimages are -1 and -2.

Modulo one, (7.1) always gives the circle-valued semiconjugacy. On the
interval, retain the marked binary expansion, or state the two exceptions.
The theorem's statement and proof are repaired, and the mechanism is
logged in MISTAKES. Its positive-integer endpoint dictionary is unaffected;
the other Diophantine-stratum assertions are not re-audited here and are
not dependencies of Propositions 1 or the explicit certificates.

## 8. Connection contracts and exact verification

| Source -> target | Map and preserved predicate | Lost data / restoring coordinate | Decisive test |
|---|---|---|---|
| Integer -> real mantissa | `n/b^L`; normalized size | L and trailing branch data under truncation | `2^K+1` versus `2^K+3` |
| Guarded word -> repeated digit block | (3.2)–(3.4); exact finite endpoint relation | Shift q, guard, odd factor u, clock | Decimal block851 yields125, not1 |
| Compressed source -> residue | Work modulo full `pq` before division | Denominator divisibility | q divisible by3 |
| Orbit itinerary -> interval | Binary reading; exact marked shift | Two expansions at dyadic endpoints | Sources -2,-1 |
| Marked Gray source -> codeword | Incoming carry-correct codec | Source when carry is erased | Same unmarked charts, wrong source word |

Reproduce:

~~~
python3 -B 04-computation/experiments/collatz_radix_families_20261018p.py --output 05-knowledge/results/collatz_radix_families_20261018p.json
python3 -B -O 04-computation/experiments/collatz_radix_families_20261018p.py --verify 05-knowledge/results/collatz_radix_families_20261018p.json
~~~

Finite universe: every odd source from 1 through 4095; the first k odd
edges for k=1,2,3,4; bases 2,4,5,7,8,10,11,16; repeat counts 1,2,3.
There are 65,536 configurations and **196,608 literal replay controls**.
Sources whose first k edges already pass through 1 are allowed here;
only the separate 75 family is asserted to be first-hit at its endpoint.
Shifted radix strings are checked against literal base conversion.

Additional universes: binary repeat counts 0..100,1000,10000,1000000;
decimal repeat counts 0..50; the eleven giant-source moduli in the receipt;
126 leading-digit hostile pairs; five signed parity-endpoint controls;
64 same-anchor base-pair compositions and the changed-anchor hostile.
The checker executes **1,133,827 active checks**. Optimized Python disables
none of them. These finite controls validate implementations; the proofs
above establish the all-parameter statements. No sampling claim is used.
