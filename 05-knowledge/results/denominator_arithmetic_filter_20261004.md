# An arithmetic realization gate after golden periodic lifting

2026-10-04 (America/Denver).

**Status:** PROVED elementary fixed-point and parity gate; FINITE-EXACT
for exact golden denominators `q=1,...,40,64,76,81,105`. No claim concerns
all denominators, additional integer cycles outside this finite universe,
or global signed or positive Collatz coverage.

## Inheritance and the missing predicate

[The denominator-fibre theorem](denominator_fibre_recursion_20261004.md)
constructs every periodic golden lift of a nonzero rational phase, with
exact period equal to its matrix-phase period. Its zero phase has three
periodic points: zero and the two-cycle `1 <-> phi^-1`. The latter needs
its explicit boundary tag; selecting only the canonical zero lift loses
the signed integer cycle through -1.

The closest mechanism is the elementary affine composition of a parity
word. The canonical hostile is already visible at golden denominator
q=3: its valid golden cycle realizes the rational cycle through 1/11,
rather than any integer cycle. The corrected near miss is identifying a
complete finite golden phase atlas with a complete integer realization
atlas. The least-used sidecar is the reduced denominator of the rational
arithmetic fixed point, distinct from the golden scalar denominator q.

The board here is a phase orbit, its exact periodic golden lift, the
ordered cyclic word, its affine return map, its rational fixed point,
and the integer-realization predicate. The map from the golden carrier
to the rational arithmetic carrier is explicit and fully checked below.

The script imports only `phase_cycles` and `periodic_lift` from the
denominator-fibre implementation. It reconstructs branch digits by exact
quadratic comparisons and independently performs all arithmetic decoding.
It does not call that implementation's inherited `decode_cycle` routine.

## 1. A necessary and sufficient arithmetic test

Use the ordinary, non-shortcut branches

    T_1(n)=3n+1,   T_0(n)=n/2.

Let w be a nonempty cyclic binary word with no adjacent ones, including
across its cyclic boundary. Let s be its number of ones and e its number
of zeros. It has e>=1. Composition gives

    T_w(n)=(3^s*n+K)/2^e.

Compute K without fractions: start `(P,K,e)=(1,0,0)`; reading a one sends
`(P,K,e)->(3P,3K+2^e,e)`, while a zero increments e. The unique rational
fixed point is therefore

    n_w=K/(2^e-3^s).                                  (1)

The denominator in (1) is nonzero: positive powers of 2 and 3 cannot agree,
and s=0 also gives a nonzero denominator. When s>0 it is relatively prime
to 6. The all-zero word has fixed point zero.

To interpret rational parity, write a rational in lowest terms with odd
denominator and use its numerator modulo 2. This agrees with ordinary
integer parity and defines the two branches on this rational domain.

**PROVED parity legality.** If a cyclic word starts with a zero, every
one contributes to K after at least one zero, so K is even. If it starts
with a one, that first one's contribution is odd; every later one has
at least one intervening zero and contributes an even term. Thus K has
the same parity as the first bit. Since `2^e-3^s` is odd, the reduced
fixed point has exactly that rational parity.

Applying the first branch to (1) produces the unique fixed point of the
cyclically rotated word: conjugate the composed affine maps, or compose
the remaining branches and the first branch in the opposite order.
Apply the preceding parity argument to every rotation. The entire cyclic
word is consequently an actual rational parity itinerary, not just a
formal selection of affine branches.

Along this rational cycle the reduced denominator r is constant and
coprime to 6. An even step divides an even numerator by 2 without changing
its odd denominator. An odd step changes `a/r` to `(3a+r)/r`, whose
numerator remains coprime to r because `gcd(r,3)=1`.

It follows that the following are equivalent:

1. every cyclic arithmetic state is an integer;
2. the single fixed point (1) is an integer;
3. `2^e-3^s` divides K.

This is the requested necessary and sufficient realization gate. The
script checks both single-root and all-state integrality independently
after literal rational replay. If the cyclic word is primitive, the
arithmetic cycle has its exact word length: a shorter state period would
force a shorter period in the deterministic parity itinerary.

The proof applies to every finite cyclic golden word. The finite census
below applies it only to the explicitly listed denominators. It does not
show that every integer's infinite golden itinerary is periodic or belongs
to any rational phase fibre.

## 2. Complete finite census

For q>1, enumerate every primitive phase pair modulo q, grouped into its
exact matrix orbits. Lift each orbit's chosen anchor and follow the upper
golden map for one phase period, comparing exact phases at every step.
For q=1 explicitly retain both zero and the boundary two-cycle.

The universe contains **1,065 golden cycles and 40,730 periodic points**.
Of these, 40,727 points have nonzero primitive phase; the remaining three
are the exceptional zero-phase points. Exactly **five integer cycles**
pass the arithmetic gate. The other 1,060 cycles are rational-only, with
209 distinct reduced arithmetic denominators across the whole census
(including denominator 1).

Here the listed root is the cycle member of smallest absolute value, not
an additional equivalence or geometric label.

| Golden denominator q | Root | Ordinary period | Cyclic parity word |
|---:|---:|---:|---|
| 1 | 0 | 1 | `0` |
| 1 | -1 | 2 | `10` |
| 2 | 1 | 3 | `100` |
| 11 | -5 | 5 | `10100` |
| 76 | -17 | 18 | `101010100101010000` |

Their complete ordinary integer cycles are

    0 -> 0;
    -1 -> -2 -> -1;
    1 -> 4 -> 2 -> 1;
    -5 -> -14 -> -7 -> -20 -> -10 -> -5;
    -17 -> -50 -> -25 -> -74 -> -37 -> -110 -> -55 -> -164
        -> -82 -> -41 -> -122 -> -61 -> -182 -> -91 -> -272
        -> -136 -> -68 -> -34 -> -17.

Every state and every prescribed parity bit is checked with rational
arithmetic. All five have arithmetic denominator 1; this should not be
confused with their distinct golden denominators.

Some useful census rows are:

| q | Golden cycles | Periodic points | Integer cycles |
|---:|---:|---:|---:|
| 3 | 1 | 8 | 0 |
| 4 | 2 | 12 | 0 |
| 11 | 13 | 120 | 1 |
| 64 | 32 | 3,072 | 0 |
| 76 | 240 | 4,320 | 1 |
| 81 | 27 | 5,832 | 0 |
| 105 | 192 | 9,216 | 0 |

At q=105 there are 96 cycles of period 16 and 96 of period 80. None passes
the integer gate. Its 20 distinct arithmetic denominators are retained in
the output, ranging from small values such as 55 to much larger values;
this is a finite enumeration, not a prime-factor criterion.

## 3. Small hostile controls

The first failure in increasing golden denominator is q=3. Its unique
golden cycle has word `10100000`. Here s=2, e=6, K=5, so

    n_w=5/(64-9)=1/11.

Its actual rational orbit is

    1/11 -> 14/11 -> 7/11 -> 32/11 -> 16/11
         -> 8/11 -> 4/11 -> 2/11 -> 1/11.

All eight branches have the correct rational parity, but no state is an
integer. Existence and uniqueness of a periodic golden lift therefore
leave a real arithmetic condition unpaid.

At q=4 the two primitive cycles have roots 1/29 and 5/7, with arithmetic
denominators 29 and 7. They are the binary lifts of the q=2 golden cycle,
but neither lifts its integer realization. At q=11 the integer -5 cycle
coexists with the rational cycle through 1/13, word `10000`. Even the
golden scalar denominator alone cannot decide integer realization. The
phase and full denominator-ideal refinements of the inherited note retain
more information, while equation (1) directly tests the required predicate.

The data-preserving connection is thus:

    periodic golden phase + upper-boundary tag
      -> ordered cyclic parity word
      -> exact affine return (3^s,K,2^e)
      -> rational cycle + reduced arithmetic denominator.

The word retains timing; the rational denominator is the integrality gate.
Discarding either coordinate can merge objects that behave differently
under arithmetic realization.

## 4. Failed integer realizations become guarded rational anchors

The incoming [rational-anchor compiler](rational_anchor_returns_20261004.md)
(commit `952b93f75`) already proves the all-prefix-growth rotation lemma
and the guarded return theorem used here. The new connection is to apply
those mechanisms to the complete phase-filter population above. An object
that fails the integer-cycle gate can still be a valid rational chart for
certifying a different, positive integer source. The chart's rational
cycle itself does not thereby become an integer cycle.

Compress a nonzero rational parity cycle into its odd states and their
exact valuations `w=(a_1,...,a_p)`, where every `a_i>=1`. Put

    A_i=a_1+...+a_i, P=3^p, Q=2^A_p,
    B=sum_(i=0,...,p-1) 3^(p-1-i)*2^A_i, A_0=0.

Its return map is `(P*n+B)/Q`, with fixed point `B/(Q-P)`.
For a negative cycle B>0 implies P>Q. The reduced anchor is
`r=-h/d`, where `h=B/g`, `d=(P-Q)/g`, and `g=gcd(B,P-Q)`.
Both h and d are positive, odd, and relatively prime to 3: B is odd and
is a power of 2 modulo 3, while P-Q is also a unit modulo 6.

**PROVED, inherited rotation mechanism.** Let
`R_i=3^i/2^A_i` for `0<=i<p`. These ratios are distinct: equality at two
indices would equate positive powers of 2 and 3. Choose the index j of
their minimum and rotate w to begin just after that prefix. A nonwrapping
proper prefix now has multiplier `R_(j+k)/R_j>1`. If it ends at p, use
`R_p>1>=R_j`; if it wraps to index k, its multiplier is
`R_p*R_k/R_j>1`, since `R_k>=R_j`. Thus every nonempty prefix expands.
The implementation compares exact Fractions and integer powers, with no
numerical logarithms.

**PROVED canonical-anchor corollary.** The cycle state of smallest absolute
value in a negative cycle is its largest state r_0<0. It must be odd:
an even negative state would be followed by the larger state r_0/2.
Every odd prefix has affine form `r_i=M_i*r_0+C_i`, with C_i>0, and
`r_i<=r_0` by maximality. Consequently
`(M_i-1)*r_0<=-C_i<0`, so M_i>1. Thus the canonical negative roots already
chosen in this census are all favorable anchors. A favorable anchor need
not be unique; this argument singles out one intrinsic choice.

The exact sign census is:

| Arithmetic sign | Integer cycles | Rational-only cycles |
|---|---:|---:|
| Negative | 3 | 48 |
| Positive | 1 | 1,012 |
| Zero | 1 | 0 |

All 51 negative cycles pass the favorable-anchor test. Independently,
all **331 cyclic odd phases** of those cycles pass the minimum-prefix
rotation construction; 231 of these phases require an actual rotation.
The 48 rational-only negative cycles have odd-period histogram
`3:1, 4:3, 5:3, 6:3, 7:33, 8:4, 10:1`. They occur at golden denominators
`11:2, 19:1, 20:1, 29:4, 36:1, 38:8, 40:2, 76:25, 105:4`, where each
entry is `q:cycle count`. The positive cycle through 1, with odd word
`(2)`, is a hostile to extending this growth claim to all rational cycles:
its full multiplier is 3/4<1, unchanged by rotation. Zero has no odd word.

For **each of the 48 rational-only negative cycles**, instantiate the
inherited compiler at the least m>=1 with Q^m>h, then the least t>=1 with

    2^t*(Q^m-h)>P^m-h.

Set `beta=h*(P^m)^(-1) mod2^t`, `K=A_p*m+t`, and

    n=(beta*Q^m-h)*d^(-1) mod2^K.                    (2)

The inherited theorem certifies the entire positive cylinder (2) as
having exact first odd descent at pm steps. The saved output retains
`(word,h,d,m,t,residue,modulus,endpoint)` for all 48 instances. This
script independently replays the least positive source and two successive
modulus lifts in each cylinder: **144 integer-source replays**. It verifies
every exact valuation, every earlier state against the original source,
and the terminal strict descent. Forty-six instances use m=2 and two use
m=1; the largest tested first-descent time is 16 odd steps. These finite
replays are controls for the inherited all-height theorem, not its proof.
No disjointness, ranking, or added coverage relative to other families is
inferred from the list.

The retained q=29 cycle gives an explicit connection between the carriers:

    Theta=(-4+20*phi)/29, raw word 1010100,
    -19/11 -> -23/11 -> -29/11 -> -19/11,
    odd word (1,1,2), P=27, Q=16, h=19, d=11.

The exact golden numerator pair is checked at the phase whose first bit
corresponds to -19/11. For m=2, the least budget is t=2, and (2) gives
`999 mod1024`. Its least positive source has route

    999 -> 1499 -> 2249 -> 1687 -> 2531 -> 3797 -> 89,

with valuation word `(1,1,2,1,1,7)`, so its exact first odd descent is at
step six. This example receives 64 consecutive cylinder-lift controls in
addition to its three general compiler controls. The weakened budget t=1
admits the inherited hostile

    487 -> 731 -> 1097 -> 823 -> 1235 -> 1853 -> 695,

which has a legal six-step shadow but no descent during those steps.
The missing coordinate is the exit budget measured against the original
integer source, not legality of the rational anchor.

The source-to-target map now has two explicit gates:

    golden periodic cycle -> rational arithmetic cycle;
    negative rational cycle + favorable phase + source congruence
      + original-source exit budget -> certified positive first descent.

The first map preserves the full cyclic word. Rotation preserves that
cycle and its valuations while changing its marked phase. The second map
uses its rational affine chart to construct new integer trajectories;
it does not preserve the original cycle's arithmetic states or integrality.
The denominator d, marked phase, repetition count, guard, and original
source are required sidecars. Neither map supplies an arbitrary integer
with a successful guarded representation, so global coverage remains open.

## Reproduction and independent paths

Run:

    python -X utf8 -B 04-computation/experiments/denominator_arithmetic_filter_20261004.py
    python -X utf8 -B -O 04-computation/experiments/denominator_arithmetic_filter_20261004.py

The [saved exact output](denominator_arithmetic_filter_20261004.out)
contains every q-row, golden-period histogram, and reduced arithmetic
denominator histogram, together with all integer cycles and the first
rational hostile examples. Normal and optimized outputs agree.
It also retains the sign counts, exact rotation controls, all 48 compiled
rational-anchor cylinders with three replayed endpoints each, and the
q=29 transfer and weak-budget hostile.

The golden phase enumeration/lifting APIs are dependencies, explicitly
separate from this new arithmetic check. Independent paths are the
integer affine composition, literal Fraction replay, all cyclic parity
tests, exact-period checks, and single-root versus all-state integrality.
No inherited arithmetic decoder, floating-point threshold, or assumed
integer-cycle list is used to generate candidates. The known five cycles
serve only as an explicit final positive-control comparison on this
declared finite universe.
