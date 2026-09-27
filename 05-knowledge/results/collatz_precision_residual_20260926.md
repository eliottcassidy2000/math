# Coefficient-descent residuals and finite-bank precision: corrected dictionary

**PROVED scoped affine inequalities / FINITE-EXACT residual counts and
bounded representative census / OPEN pointwise coverage.**
Original session: opus, gilbreath6-collatz-precision-20260926.
**Correction and integration, 2026-09-27:** the exact-word modulus, the
finite-bank residual, the meaning of entry, and the all-j cutoff were
repaired during an independent audit. The prior published formulation is
preserved in Git history at4403f0015; it is superseded where it conflicts
with this note. The computed density table survives the repairs.

Inherits [the reset bank](reset_20260926_swaplift.md),
[THM-4495, no-descent count order](../../01-canon/theorems/THM-4495-collatz-no-descent-exact-order-spitzer.md),
and the stopping-time references in
[CORE-PAPERS](../reference/CORE-PAPERS.md). The current certificate theorem is
[THM-4512, corrected exact/coarse cylinders](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md).
The [independent incoming audit](entry_20260927_incoming.md) records minimal
witnesses and exact replays. No new literature-priority claim is made.

## 1. Three different unresolved sets

Let U(n)=oddpart(3n+1). For the ACTUAL first j valuations of a fixed source
n, let A_j be their sum and let S_j be the usual positive affine carry:

    U^j(n)=(3^j*n+S_j)/2^(A_j).

Coefficient descent is 3^j<2^(A_j); actual descent below the ORIGINAL n
requires the stronger inequality

    (2^(A_j)-3^j)*n>S_j.

The no-coefficient-descent set at an odd-step horizon k, the no-actual-
descent set at that horizon, and the complement of a particular finite
certificate bank are different definitions. Actual descent forces
coefficient descent. The finite bank is only a selected collection of
sufficient guards; its complement can contain sources already certified
by other short words.

An exact example is4091. It is outside every row of the inherited171-row,
65-cylinder bank. Its first coefficient descent and first actual descent
both occur at step8, ending at1639. It belongs to the new all-k macro
family in [the recursive note](entry_20260927_recursive.md). Thus the
finite bank's insufficient-precision residual is NOT exactly the
classical no-coefficient-descent residual.

## 2. Exact words need the final oddness bit

For a prescribed positive word w=(v_1,...,v_j), write A=sum(v_i) and
S=S_j(w). There are two useful cylinders:

| Object | Congruence | Actual consequence |
|---|---|---|
| Exact valuation word | n=(2^A-S)*3^(-j) mod2^(A+1) | U^j(n)=(3^j*n+S)/2^A is odd |
| Coarse integrality class | n=-S*3^(-j) mod2^A | First j-1 valuations exact, final at least v_j; U^j(n)=oddpart((3^j*n+S)/2^A) |

The exact class has density2^(-A) among odd integers. That is the weight
used correctly by the residual-density programme. The original note
attached that weight and the exact endpoint equality to the coarser
modulus2^A class; those assertions were incompatible.

For word(1), the exact class is3 mod4 and the coarse class is every odd
integer. At n=1 the formal quotient is2, whereas U(1)=1. One missing
binary bit is enough to invalidate an asserted exact arithmetic path.

When 3^j<2^A, put N(w)=S/(2^A-3^j). On the exact cylinder, descent is
equivalent to n>N(w). On the coarse cylinder it is sufficient, because
extra final divisions only reduce the endpoint. This is why the useful
all-height threshold theorem survives the modulus correction.

## 3. Residual density: finite-exact coefficient counts

Define D(k) as the sum of2^(-A_k) over exact valuation words with no
coefficient descent at any prefix of length at most k. This is a natural
density among odd integers. It is computed by a finite dynamic programme:
each surviving prefix has A_j<j*log_2(3), so its next admissible valuations
are bounded. The independent audit uses integer power comparisons.

    k:     1    2    3     4      5      6      7         8
    D(k): 1/2  3/8  1/4  13/64  19/128  1/8  113/1024  367/4096

    k:    12       16       20       24       28       32
    D(k): .05212   .03092   .01925   .01314   .00875   .00593

    k:    36       40       41       48       52       56       60
    D(k): .00427   .00299   .00265   .00156   .00112   .00082   .00062

In particular,

    D(41)=1530343662856563/2^59,
    D(60)=12185976031772265023015553/2^94.

Thus99.73 per cent of odd inputs have coefficient descent by step41,
whereas the selected reset bank covers37.87 per cent of odd inputs with
actual first-descent certificates. The new macro family improves that
selected bank, not the general residue sieve.

For any fixed horizon, the possible small threshold exceptions form a
finite set: there are finitely many preterminal no-descent words and
sufficiently large final valuations give N(w)<1. Hence coefficient and
actual first-descent coverage have the same density at that fixed horizon,
without asserting pointwise equality of their stopping times.

The inherited T-coded counts and Syracuse counts use different clocks.
The earlier asymptotic discussion converts the exponential cost per map
step to cost per odd step, giving (1-h*)/rho*, approximately0.0793, with
rho*=log_3(2). It must not be quoted as1-h* per odd step. That asymptotic
comparison is not a dependency of any new certificate or correction here.

## 4. The one-member theorem has an explicit finite scope

For any positive valuation word,

    S_j<=2^A*((3/2)^j-1).

The integer test

    3^j-2^j < (2^A-3^j)*2^j,
    A=bit_length(3^j),

holds for every j<=5000. The worst ratio is211/416 at(j,A)=(5,8).
Therefore N(w)<2^A in that range and the coarse cylinder has at most
one member not certified by this threshold, namely its least positive
representative. An exact subcylinder inherits that sufficient result.

The independently reproduced606746 first-coefficient-descent labels of
length at most14 leave only the representative1 unproved by the threshold.
The finite enumeration stops at a bound on TOTAL A, not forty additional
valuation bits for every prefix. Its omitted tail is justified by

    S_j<=j*3^(j-1),  N(w)<1 for j<=14 and A>=41.

The producer separately reports sigma=sigma_inf on every odd source
3..10^7, with maximum stopping time155. This remains a finite computation.
The integration audit verifies the relevant integer cutoff comparisons
but does not claim a new independent census of all ten million inputs.

An effective irrationality estimate supplies a sufficiently large
computable cutoff only after its constants are specified. One must then
check the remaining gap between5000 and that cutoff. The original wording
that any effective measure immediately covers every j>5000 was unsupported.
The all-j extension is not a proved dependency.

## 5. Entry must retain its target

Every positive odd Collatz orbit already reaches1 or a locally descending
state1 mod4: the initial valuation-one run has exactly v2(n+1)-1 steps.
Therefore ordinary moving-orbit entry into a local-descent region is
unconditional, and is not equivalent to Collatz.

The equivalent global target is a ROOT-relative one: for every odd n>1,
some actual iterate is less than that original n. Strong induction then
gives convergence. Alternatively one can prove entry into a region whose
entire future is already certified to reach1. These targets must not be
identified with a locally decreasing step after earlier growth.

The density2^(-L) avoidance bound concerns a moving local predicate:
having L successive valuation-one steps. It cannot pay the original
source's accumulated growth. The [entry note](entry_20260927_entry.md)
provides exact reductions, an abstract grow/drop countermodel, and a
total certificate interpreter with rejection allowed.

## 6. Connections retained after the correction

The affine carry is a precise interface between binary valuation data,
ternary factors, and ordinary integer size. The new
[three-type excursion controller](entry_20260927_board.md) retains this
interface: each repeated(1,2) block updates(k,r,b)->(k-1,r,9b), and a
terminal inequality checks the original source. This is an operational
connection, not a density-to-pointwise inference.

The Langlands/loom analogy in the original session was explicitly
speculative; no reduction or automorphic correspondence was established.
It can motivate retaining both arithmetic sides and their shared carry.
The claim that Collatz asks ordinary integers to be Haar-generic was
incorrect: a convergent Syracuse orbit has valuations eventually all2,
a periodic tail rather than a generic independent valuation sequence.
The useful distinction is measure control versus pointwise arithmetic,
not an assertion of generic coding for convergent integers.

The Gilbreath comparison has a concrete surviving mechanism in
[the carry-wall note](entry_20260927_wall.md): boundary dependence disappears
within one scan, while the exact unread tail remains. An irreversible
phase across successive Collatz iterates has not been proved.

## 7. Reproduction

Original producer scripts and outputs remain preserved:

    python 04-computation/experiments/collatz_precision_residual_20260926.py
    python 04-computation/experiments/collatz_coefficient_stopping_20260926.py

Their historical prose and loop descriptions should be read with this
correction. The independent integer audit is:

    python -X utf8 -B 04-computation/experiments/entry_20260927_incoming.py

Global source coverage remains OPEN. Expanding a selected bank is useful
for a compact proof grammar, but the already strong density sieve does not
remove the pointwise obligation.