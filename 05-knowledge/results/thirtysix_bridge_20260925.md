# Thirty-six: relation matrices, half-step phases, and a signed Pell clock

**Date:** 2026-09-25. **Status:** PROVED scoped constructions and identities;
INHERITED PROVED mechanisms identified by source; FINITE-EXACT controls.
Completeness of the minus-cycle list, positive Collatz convergence, and a
prime-distribution theorem from these structures remain OPEN. No novelty
claim is made for Pell equations, finite permutation roots, or Lucas
divisibility. The new work connects explicitly specified representations.

## Inheritance, portfolio, and concept board

The exact earlier36 is in [level11_short_20260922](level11_short_20260922.md),
under the two gluing and diagonal sections: six signed cycle representatives
have36 ordered pairs including diagonal pairs. The objection that there are
only six representatives does not refute that relation-matrix count. It does
rule out treating36 as a count of the same six objects or as their tournament
arc count, which is15.

Closest proved mechanisms: the signed Pell half-step of
[THM-3742, square-triangular Pell central-sign/projective cycles](../../01-canon/theorems/THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle.md),
the parameter selector of
[THM-3335, square-triangular Pell/Pythagorean selector](../../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md),
the [signed sibling fibres](arithmetic_braids_20260917_collatz.md), and the
[four-color ten-letter construction](duck_decoder_20260925.md).
The Q15 endpoint mechanism is already proved in
[seam_threshold_20260925](seam_threshold_20260925.md).

Canonical hostile: a periodic modular clock is not a global Collatz orbit.
Corrected near miss: conflating the period of a cycle with a count of cycle
representatives. Least-used sidecars: the actual time convention, relative
phase of paired cycles, and the signed Pell norm. The earlier level-eleven
note also contains a real eighth-root recurrence, recovered separately in
the [digit companion](thirtysix_digits_20260925.md).

Anchor: recover an exact36 and a lawful half-step. Niche: the two partitions
of the ten-letter object. Wildcard: transport the decimal17 phase to a signed
Pell orbit, then test its prime and integer-lift boundaries.

| Live concept | Actual object | Retained coordinate |
|---|---|---|
| Signed cycles | permutation of a specified known cycle locus | ordinary/shortcut/odd clock |
| Half-step | a square root of that permutation | phase and extra copies where necessary |
| Multipartite quotient | partition plus uniform crossing relation | part sizes, internal arcs, orientation gauge |
| Decimal sequence | a_k=(10^k-7)/3 | digit index, residue, ordinary magnitude |
| Pell interlacing | norm+1 and norm-1 integer points | both coordinates and norm sign |
| Square-sum path | Q15 with endpoints8,9 | exact square edges and endpoint status |

Used research moves: identify the observed object before comparing counts;
retain the missing phase/quotient coordinate; test an apparent recurrence
through the first return and the first hostile. No new META-PATTERNS card is
promoted on this session alone.

## 1. Several exact36s, with distinct maps

The companion proofs make the following meanings concrete:

| Carrier | Why36 appears |
|---|---|
| six signed minima ±{1,5,17} | 6² ordered pairs, with loops |
| ordinary minus18-cycle through17 | two interleaved copies make a half-step36-cycle |
| minimal square-root completion of all three ordinary minus cycles | 2*18=36 relative phase choices |
| ten letters, partitioned by C3 orbit | K_(1,3,3,3) has36 crossing edges |
| same letters, partitioned by XOR charge | K_(4,2,2,2) has36 crossing edges |
| nine nonfixed letters with all pairs restored | K9 has27 crossing+9 internal=36 edges |
| consecutive integers8,9 | T8=8*9/2=36=6² |

Equal cardinality alone does not equate these carriers. The useful connections
are the phase construction, the six-edge trade, and the Pell chart below.

### The half-step is sensitive to the clock

The gluing sign matters. Negation gives U_+(-n)=-U_-(n), with the analogous
identity in ordinary and shortcut time. Joining positive3n+1 to negative
3n-1 therefore gives two reflected copies of the **positive plus** system.
The gluing with the six known1,5,17 cycles uses positive3n-1 and negative
3n+1: on odd signed n its ordinary branch is3n-sign(n). The two gluings
are different piecewise maps, not different drawings of one fixed rule.

The known positive3n-1 cycles have lengths1,2,7 in odd-only time;1,3,11
in shortcut time; and2,5,18 in ordinary time. The representative17 therefore
names a7-,11-,or18-cycle depending on what was counted as a step.

A finite permutation has a square root iff every even cycle length occurs
an even number of times. Odd cycles can be rooted individually; an even
cycle must be paired with another of its length and interleaved.
Consequently, the minimal square-root completion of the ordinary minus
cycles adds copies of the2- and18-cycles, taking25 states to45. There are
exactly2*18=36 roots on that specified labelled completion. Fully duplicating
all three cycles instead gives50 states and216 roots, since the two5-cycles
can be rooted separately or interleaved in five ways.

There is a global formal half-step too. For any positive map P define
G(epsilon*m)=epsilon*P(m). Set R(m)=-m and R(-m)=P(m), for m>0. Then
R²=G on all signed states. This chooses a phase; it is not an ordinary
arithmetic Collatz step or a proof that any new orbit converges.
See [the signed-cycle companion](thirtysix_signed_20260925.md) for proofs,
the reflection obstruction, and independently enumerated phase counts.

### A different, arithmetic half-step on inverse fibres

For sigma in{+1,-1}, write U_sigma(n)=(3n+sigma)/2^k with exact
k=v2(3n+sigma), n positive odd. On an edge together with its sign, define

    F(n,sigma)=(2n+sigma,-sigma).

Then

    U_(-sigma)(2n+sigma)=U_sigma(n),
    F²(n,sigma)=(4n+sigma,sigma),
    F³(n,sigma)=(8n+3sigma,-sigma).                      (1)

Each F step increases the edge's halving exponent by exactly one. For a
fixed positive odd target y not divisible by3, the entire joint inverse
fibre is one ray indexed by k>=1:

    sigma_k is the sign congruent to2^k*y modulo3,
    n_k=(2^k*y-sigma_k)/3.

F translates k by one. Thus the two sign-specific fibres are literally
alternate levels of one ray, and the eightfold rule is three half-steps.
For y=1 this recovers the inherited Jacobsthal trunk1,1,3,5,11,..., with
the first two1s distinguished by sign. This is a statement about edge
preimages sharing a target, not a conjugacy of forward plus/minus orbits.

## 2. The same ten letters have two36-edge multipartite graphs

Use D=Sym²_set(F2²), the unordered pairs with repetition from{0,R,G,B}.
Its C3 rotation has orbit sizes1,3,3,3; its pair-XOR charges have fibre
sizes4,2,2,2. For either partition, join exactly the vertices in different
parts. The edge count is

    ((sum n_i)²-sum n_i²)/2=(100-28)/2=36.

These relations have the same vertex set and share30 edges. To change the
orbit-part relation to the charge-part relation:

* delete00--cc for each nonzero color c, and the three equal-charge edges
  between0c and the pair of the other two nonzero colors;
* add the triangle on{0R,0G,0B} and the triangle on{RG,GB,BR}.

This exchanges six edges for six. The symmetric difference is a disjoint
K_(1,3) and triangular prism. It changes all ten degree parities: the orbit
graph has degrees9,7,...,7; the charge graph has four6s and six8s. The
operator is C3-equivariant and keeps the repeated-diagonal information.

Deleting00 leaves three classes of size3. Joining different classes gives
K_(3,3,3) with27 edges. Adding each internal triangle supplies9 more,
giving K9 with36 edges. With an explicit cyclic order of the three
channels and their inherited internal rotations, this becomes the regular
tournament C3[C3,C3,C3], with four wins and four losses at each vertex.

The user's four-part quotient observation is exact with its orientation
condition: identifying each of four nonempty parts of a complete four-partite
graph to one vertex gives K4. This is a partition quotient, not an edge-
contraction claim about disconnected branch sets. Orienting each crossing block uniformly according to Q gives the
four-vertex tournament Q. Mixed directions instead give a quotient with
two-way arcs; absence of internal edges does not make the original graph
a tournament. The [multipartite companion](thirtysix_multipartite_20260925.md)
gives the quotient, rank, parity, and nonplanarity controls.

## 3. What8 and9 do in the first square-sum path

The unique Q15 Hamiltonian path, up to reversal, is

    8,1,15,10,6,3,13,12,4,5,11,14,2,7,9.

Vertices8 and9 are its endpoints because they are the graph's two leaves.
In Q14 there was also leaf10; adding15 absorbs that obstruction through
the ear1--15--10. The path's edge sums are

    9,16,25,16,9,16,25,16,9,16,25,16,9,16.

There is no36 edge sum in Q15 (the largest possible sum is29), and8+9=17
is not a square. The triangular36 is a different function of the endpoints.

It nonetheless supplies a precise Pell connection:

    T_n=q² iff (2n+1)²-8q²=1.

At n=8,q=6 this is17²-8*6²=1. The familiar primitive triangle
(8,15,17) gives another exact local chart: the odd-root parametrization
of the actual accelerated plus edge3->5 is

    (s*t,(s²-t²)/2,(s²+t²)/2)=(15,-8,17), (s,t)=(3,5).

The same small integers thus occur in specified square-sum, Pell, and
Collatz-edge constructions. This finite correspondence does not produce a
map between their entire recursive graphs.

## 4. A commuting decimal / signed-Pell clock modulo17

Use the inherited integral half-step

    S(x,s)=(x+2s,x+s),         x²-2s² -> -(x²-2s²).     (2)

Even powers of S applied to(1,0) give the square-triangular branch:

    S²(1,0)=(3,2),       n=1,T_n=1;
    S⁴(1,0)=(17,12),     n=8,T_n=36;
    S⁶(1,0)=(99,70),     n=49,T_n=1225.

THM-3742 studies this half-step and its norm sign, including a mod13
projective analysis. Here a mod17 specialization matches the decimal clock.
Modulo17,6²=2. The coordinates u=x+6s and v=x-6s diagonalize S:

    u -> 7u,      v ->12v,
    7=10^9,       12=10^15 mod17.                       (3)

Since10 has order16, S^8=-I and S^16=I. Let a range over F17 except9,
and put z=3a+7 (nonzero). For eta in{+1,-1}, define

    Psi_eta(a)=((z^9+eta*z^15)/2,
                (z^9-eta*z^15)/12) mod17.              (4)

The divisions use inverses in F17. For A(a)=10a+21,

    Psi_eta(A(a))=S(Psi_eta(a)).                        (5)

Proof: z(A(a))=10z(a);(3) multiplies the two channels by exactly the
required powers. The inverse is explicit: recover u=x+6s, then z=u^9
because9²=1 mod16, and a=(z-7)/3. Therefore each chart is bijective onto
one16-cycle. Its norm is

    x²-2s²=eta*z^24=eta*z^8 in{+1,-1}.                 (6)

The two images are disjoint and exhaust all32 points with norm±1.
Indeed their invariant v/u^7 equals eta, and every point of norm±1 has
one of these two values. Each orbit has eight states on each norm fibre.
This is an actual finite dynamical isomorphism with a specified inverse,
not a claim that the entire Pell conic has only16 points.

For the decimal sequence a_k=(10^k-7)/3,

    Psi_+(a_k)=S^k(1,0) mod17.                          (7)

Here k=0 merely extends the affine clock to a_0=-2; the displayed positive
digit family starts at k=1. In particular

    a_(k+8)=1-a_k mod17,
    17 divides a_k iff k=9 mod16.                      (8)

The decimal row k=1 has residue1 and Pell state(1,1); row k=9 has residue0
and Pell state(-1,-1). These are genuinely eight half-steps apart. The
integer square-triangular event n=8 instead occurs at Pell half-step k=4.
The numeral8 is not a license to identify those two index conventions.

### A surviving prime-index obstruction, and its limit

For x_l+q_l*sqrt8=(3+sqrt8)^l, the Pell traces begin

    x_l:1,3,17,99,577,3363,19601,... .

Modulo17,3+sqrt8 can be read as8 and its inverse15. The element8 has
order8 and fourth power-1. Hence x_l=(8^l+8^(-l))/2 vanishes exactly when
l=2 mod4. The later trace19601=17*1153 is a hostile to universal primality.

More generally x_(ab)=T_b(x_a), where T_b is the integer Chebyshev
polynomial. For odd b>1 it is divisible by x_a, with a larger positive
value. Thus if x_l is prime and l>=1, l must be a power of2. This is only
a necessary index condition; it neither proves all surviving terms prime
nor changes the distinct decimal phase condition(8).

## 5. Decimal spacing and the prime return

The user clarified that the intended numbers are the literal decimal
strings31,331,3331,... . Their exact ordinary laws are

    a_(k+1)=10a_k+21,
    a_(k+1)-a_k=3*10^k,
    a_(k+1)/a_k=10+21/a_k.

Thus each row adds one digit. Neither4.2 nor8.4 is its unscaled gap or
successive ratio. The constant21 belongs to the affine recurrence and
doubles to42; dividing a chosen coordinate by5 would turn those constants
into4.2 and8.4, but that normalization was not supplied. The modulo17
eight-step reflection has a definite formula independent of such rescaling.

The [digit companion](thirtysix_digits_20260925.md) extends the finite
primality audit: the first seven displayed terms are prime, k=9 through17
are composite, and k=18 returns to prime. The latter has a complete
factorization/order certificate, not merely a probable-prime test. The
construction continues through all of these events. Its first two odd
Collatz steps descend for k>=4, but later forced rising steps show that
this does not give a stepwise monotone complexity rank.

## Reproduction, correction, and open boundary

Run from the repository root:

    python3 04-computation/thirtysix_bridge_20260925.py
    python3 -O 04-computation/thirtysix_bridge_20260925.py

The [script](../../04-computation/thirtysix_bridge_20260925.py) and
[frozen output](thirtysix_bridge_20260925.out) check Q14/Q15 directly,
Pell half-steps through24,100 trace indices, all32 signed-norm points and
both inverse charts,65 decimal indices, and all eight triple orientations.
An independent lane audited the finite Pell/decimal commuting square and
both charts' complete32-state image, the inverse maps, clock separation,
and the necessary prime-index condition. Companion files state their own
universes, arithmetic certificates, and independent checks.

During the inheritance pass, a known-local error was found in
[THM-060, bipartite skeleton](../../01-canon/theorems/THM-060-bipartite-skeleton.md):
reversing all arcs of a zero-backbone triple contributes0 for a transitive
triple and2 for a cyclic triple, not always2. The local proof is repaired;
the needed evenness survives. This session does not re-audit its other
global tiling claims. The correction lineage is retained in MISTAKES.

What remains unsupported is the implication from these finite or inverse-
fibre structures to complete cycle classification, a root certificate for
every integer, or generation of the primes. The next decoder must preserve
ordinary value, sign, exact clock, and the intended target predicate through
each map; a finite modular conjugacy alone does not lift them.
