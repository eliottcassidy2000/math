# Collatz crossings, exact binary families, and the limit of finite arithmetic features

**Current synthesis, 2026-09-26. PROVED / INDEPENDENTLY AUDITED /
FINITE-EXACT / CITED / OPEN Collatz.** The principal new theorem is
[THM-4507, finite polynomial-valuation potential obstruction](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md).
No proof of Collatz or of the full 1/3--2/3 conjecture is claimed.

## 1. The strongest outcomes

1. A finite graph really does recognize an arithmetic reduction: an
   eight-state binary automaton selects larger starts that provably merge
   with a smaller start. These form one quarter of the odd integers.
2. Actual dyadic crossing words have an exact finite-height Sturmian
   certificate. A forbidden block forces a visit below an explicit height.
   But every allowed finite block also occurs with pure odd growth.
3. No positive logarithmic height plus an arbitrary nonlinear correction
   from finitely many polynomial valuations can decrease globally at
   every Collatz step. Rational periodic shadows prove this, not statistics.
4. Two different marked units of visible length one repair the proposed
   Fibonacci splitting. The supplied colour sequence remains authoritative;
   its 35th entry rules out the simplest Fibonacci continuation.
5. Gaussian squaring preserves complete primitive-triple hypotenuse fibres,
   but finite rational charts cannot encode unbounded Collatz states into
   a degree-two map. Theta to the fourth power supplies full-support
   arithmetic weights; its masked decay remains an open proof obligation.

## 2. Inheritance and the six-concept board

Anchor: the orbit-coupled crossing potential. Niche: exact fate-preserving
binary reductions. Wildcard: marked Fibonacci units, Gaussian fibres and
theta weights. The repository's broader LRC(14) anchor remains open; no
Collatz-to-LRC transfer is inferred from shared constants.

Closest proved mechanism: THM-4476's reciprocal summability and THM-4483's
expanding-cycle charge obstruction. Canonical hostile: arbitrarily long
positive shadows of a negative rational cycle. Corrected near miss: a
balanced comparison under a uniform extension law need not survive a
fixed-source residue projection. Least-used sidecars: polynomial-root
avoidance at several primes and exact synchronized coalescence time.

| Concept | What it contributes | What still has to be proved |
|---|---|---|
| phase and crossings | bounded crossing discrepancy; finite forbidden-block core certificate | relate them to unbounded halving depth |
| finite arithmetic features | exact local predicates and useful small automata | global descent cannot follow from a finite valuation correction |
| poset extensions | comparison laws under a declared measure | retain nonzero residue Fourier modes and affine carry |
| binary families | exact smaller-start reduction and stopping-time identities | connect the alternating pairs across a whole family |
| Fibonacci boundary | two distinct unit charges and exact carry lift | choose a rule fitting all supplied data, then couple actual orbit order |
| Gaussian/theta arithmetic | whole-fibre bijection and full-support positive weights | preserve ancestry or prove masked mass decay separately |

The principal update to the board is constructive and restrictive at once:
a small graph can certify a particular arithmetic predicate while failing
to contain the information needed for global descent.

## 3. Beyond 27: a binary family with an exact reduction

Use U(n)=oddpart(3n+1). Write

    n=2^r u-1, u odd, r>=1;  s=v_2(3^r u-1).

The [run note](crossroads_crossing_20260926_arithmetic_run.md) proves

    U^(r+1)(n)=U^(r+1)(2n+1)  if and only if s=1.

Equivalently, `3^r u=3 mod4`. This gives a valid smaller-start reduction
for m=2n+1, regardless of whether the common future is lower than m.
The eligible larger m have relative natural density1/4 among odd integers.
This is a source density, not an orbit frequency.

For larger m, read its binary digits from the least-significant end:

    m: 1^R 0 b..., R>=2.  Accept exactly when b=R mod2.

Implicit leading zeros handle finite Mersenne strings. An eight-state
automaton checks this predicate exactly. It still retains the integer
to compute `(m-1)/2`; it is not a finite-state simulation of Collatz.

On a fixed-u tower the paired exponents alternate. If u=1 mod4, pairs are
(1,2),(3,4),...; if u=3 mod4, pairs are(2,3),(4,5),... . Thus u=7 gives

    (27,55), (111,223), (447,895), ... .

For an eligible n>1 that reaches1, the two starts have equal odd-step
counts, while both the shortcut and standard lengths differ by exactly1.
For27/55 these are41/41 odd steps,70/71 shortcut steps,111/112 standard
steps. This is a genuine mechanism involving223; it does not identify233
with223 or explain their difference10. The earlier233 connection remains
the exact Fibonacci/poset count described in the previous session.

The pairs do not join into one infinite convergence chain. Two successive
binary reductions are forbidden by the alternating mod4 condition.
Also31 and63 merge at91, above both. Coalescence preserves eventual fate,
not a height decrease. For Mersenne starts, the first maximal-run reset
descends exactly at exponents2,4,8; all later exponents require more work.

## 4. The crossing potential is exact, but one coordinate is missing

For successive odd states x_j, put

    a_j=v_2(3x_j+1), e_j=a_j-1,
    phi_j={log_2 x_j}, alpha=log_2(3/2),
    c_j=log_2(1+1/(3x_j)),
    u_j=the 0/1 upward dyadic crossing at the odd half-step.

Then exactly

    u_j=phi_j+alpha+c_j-phi_(j+1),
    h_(j+1)-h_j=u_j-e_j.

On a hypothetical nonperiodic positive orbit, THM-4476 gives summable
reciprocals, hence finite total carry C_infinity. Thus

    U_J-alpha J=C_J+phi_0-phi_J

has bounded discrepancy on this one fixed orbit. Its odd-state binary
mantissas are Benford; the full shortcut orbit has a companion base-three
rotation clock. Both statements preserve the source. Neither proves
independence of leading phase and trailing valuation.

There is also an unconditional finite certificate. Sort the rational cuts
`2^ceil(log_2 3^j)/3^j`, 1<=j<=L, together with1,2; let rho_L be the least
adjacent ratio. If every odd source in a block is at least m and

    (1+1/(3m))^L < rho_L,

its crossing word belongs to the L+1 length-L rotation words. All tests
use rational arithmetic. Conversely, a forbidden word certifies a source
below that m. This is not an assertion that an entire late tail becomes
one exact Sturmian word: the start of the valid tail can depend on L.

The hostile is sharp: every one of those finite allowed words occurs at
arbitrarily large positive starts with all a_j=1. No extra halving occurs.
An actual completed excursion27->41->62->31 also returns to its original
dyadic band at a larger value. The [crossing note](crossroads_crossing_20260926_potential.md)
records these results and their independent audit.

## 5. The potential class now ruled out

THM-4507 excludes

    V(n)=c log n+f((v_p(P_i(n)))_(p,i))+g(n),
    c>0, finite prime set and polynomial list, arbitrary f, bounded g.

It also excludes any correction bounded on each fixed feature fibre,
including arbitrary interactions with finitely many colours. The theorem
also covers every-step monotonicity for the maximal-rise reset map C;
compressing each rise and reset does not evade this particular obstruction.
It extends the finite linear-counter part of THM-4483; it does not replace
that theorem's stronger treatment of certain infinite counter banks.

The proof chooses an expanding negative rational cycle avoiding every
selected polynomial zero. CRT produces arbitrarily long positive orbit
segments with the same endpoint features and unbounded height gain.
These refute the universal potential inequality. They do not produce one
positive divergent orbit, because the starting integer changes with the
requested prefix length.

The simple visible instance is41->31, which regenerates five units of
v_2(n+1) from one. Its general reset family remains entirely in the units
modulo30 for suitable exponents. Adding the primes2,3,11,5 does not remove
the obstruction; the full theorem allows any fixed finite prime list.

## 6. What the other seeds contribute

The exact supplied35-colour word disagrees with the Fibonacci candidate
at35. Its compressed word has equal-length factors with blue counts1 and3,
so no binary mechanical slope/intercept fits it. For the near-fit candidate,
blue positions are `2 floor(k phi)+k`; its bounded discrepancy is useful
as a spatial correction but cannot pay unbounded orbit growth.

There is a positive realization of the hidden-boundary idea. For Fibonacci
weights f_0=1,f_1=2,

    f_k=1*+f_(k-1)+f_(k-3)+... .

The marked unit has two possible charges(1,0),(-1,1), both of visible
value1 under(A,B)->A+2B. They alternate with k. The
[colour note](crossroads_crossing_20260926_colour.md) proves this and the
exact multiplication-by-three carry; a small colour quotient loses it.

The prime set{2,3,11} is the recovered odd-square-bracket classification.
The [arithmetic note](crossroads_crossing_20260926_arithmetic.md) also proves
the bijection of primitive-triple fibres P_c<->P_(c^2) under Gaussian
squaring, while retaining its sign and ancestry limitations. Finite
rational charts cannot semiconjugate Collatz into a higher-degree rational
map except with finite image; rational-map degree gives the obstruction.

The supplied theta square misses integers such as3. Theta to the fourth
power has positive coefficient r_4(n) at every n. Its weighted survivor
mass `sum_(n>=2, no descent through k) r_4(n)q^n`, at fixed0<q<1, tends
to zero iff Collatz holds. Positivity preserves every fixed integer, but
the orbit mask destroys the unmodified modular identity. The decay
estimate is still missing.

The [dyadic note](crossroads_crossing_20260926_dyadic.md) audits the pasted
functions: the first has factor-four dyadic jumps; the second integrates
to a strictly increasing reparametrization of height. The printed
fractional-part series is actually a floor series away from integers.
This recovers the B_1 boundary issue exactly, without inferring drift.

## 7. Poset-adjacent proof attempt and concrete next obligations

The useful poset is the containment forest of actual matched crossing
excursions. Its vertices are events with integer, valuation and carry
labels. Disjoint excursions can admit different orders of verification;
their arithmetic maps cannot simply be permuted. Uniform linear-extension
balance therefore needs a separate measure argument.

For a fixed-count two-chain parity poset P, with a odd letters among L,
an exact alternative is a carry polynomial in the cyclic group algebra:
`F_P(z)=sum_(w in LE(P))z^(C_w) mod(z^(2^L)-1)`. The coefficient at
`-3^a n mod2^L` selects0 or1 words for the fixed integer n. Here the
extensions correspond bijectively to binary parity words; this claim
does not automatically transfer to an excursion forest whose extensions
may encode the same parity word more than once. Evaluating at z=1
discards all nonzero residue Fourier modes; a balanced comparison can
become deterministic. This locates the missing information precisely.

A proof attempt worth continuing has three distinct obligations:

1. Use crossing certificates and the binary reducer to select arithmetic
   return blocks without changing the starting integer.
2. Keep a variable-depth address or an unbounded cost inside each fixed
   valuation fibre, and bound reset costs across enclosing excursions.
   Defining that cost by the unknown stopping time would be circular.
3. Prove an endpoint-stable inequality for the actual weighted forest or
   the fixed-q theta mass. The attached poset estimates alone do not give it.

Incoming THM-4506 sharpens the worst landing multiplicity and leaves only
an orbit-coupled average as the possible gain in that counting method.
The parallel crossing note independently recovers the phase rotation and
past/future address relation. Its finite balance and potential claims
received the corrections in the [integration audit](crossroads_crossing_20260926_integration.md).
No all-level balance, uniform-law transfer, or Collatz conclusion survives
merely by naming an analogy.

Reproduction, universes, independent controls and source pins are in the
[session audit](crossroads_crossing_20260926_audit.md).
