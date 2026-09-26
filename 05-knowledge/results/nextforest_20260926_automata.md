# The ternary unit tower behind a dyadic Collatz carry forest

**PROVED elementary constructions and scoped obstructions / FINITE-EXACT
controls / OPEN Collatz, 2026-09-26.** No claim of literature priority.
The new content is the complete reachable residual-state classification,
its exact horizontal cycles, and the macro-carry reset census. The
infinite-state conclusion itself is inherited, not rediscovered as new.

## 1. Inheritance and the connection that survived

The anchor is the exact source and open boundary of the
[ordered excursion forest](forest_20260926_excursions.md). The niche is
the native section operation from Rule 30; the wildcard is arithmetic
compactification by a finite ternary fraction. The live concepts are
source cylinders, residual functions, macro carries, synchronization,
ternary denominator depth, and the ordinary-integer end marker.

The closest proved mechanism is C10 of
[the creation decoder](creation_decoder_20260925.md#7-a-second-obstruction-the-itinerary-encoder-needs-infinitely-many-residual-states):
the full parity encoder has infinitely many residual functions, whereas
a single multiplication scan has a three-state controller. The hostile is
the all-ones source prefix, whose unread-tail slope is `3^k`. The corrected
near miss is inferring temporal termination from spatial scan termination.
The least-used sidecar is the *entire residual affine map*, not merely its
carry-state count or output parity.

Three older Rule-30 mechanisms were inspected, with their scopes retained:

- [THM-3471, Motzkin strip and innovation carry](../../01-canon/theorems/THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum.md),
  section 7: exact time blocking requires a growing alphabet. Its maximal
  universal block rank is not a fixed-seed complexity lower bound.
- [THM-4204, de Bruijn reset](../../01-canon/theorems/THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation.md):
  a spatial inverse-product reset can become blind to later temporal data.
- [THM-4210, dyadic block-current Cartier tree](../../01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md):
  the transverse channel is needed for physical reconstruction even when
  a selected observer does not see it directly.

The transfer below is a shared *resource argument about exact sections*,
not a claim that Rule-30 Cartier extraction is the Collatz operation.
Rule-30 rows are over F2; the Collatz carry below is an ordinary integer.

## 2. Every source prefix has one exact residual affine map

Write `T(n)=n/2` for even n and `(3n+1)/2` for odd n. Let Q be the
compatible least-significant-bit-first parity encoding on Z2. A length-L
input prefix is an integer `0<=r<2^L`. If its first L shortcut steps
contain a odd steps, put `b=T^L(r)`. For every unread tail q in Z2,

    T^L(r+2^L q)=3^a q+b,
    residual_Q(r,L)(q)=Q(3^a q+b).                    (A1)

The first equation follows by composing the actual parity branches. It
is valid on the entire input cylinder because those L parities depend
only on its low L bits. Finite parity encodings are bijections: each new
input bit changes the last parity through an odd slope. Hence Q is
injective. Two residual maps in (A1) are equal iff their pairs `(a,b)`
are equal: apply Q-injectivity and evaluate their affine arguments at
q=0 and q=1.

The exact deterministic transition from state `(a,b)`, reading input bit
d, is

    m=3^a, v=b+m d, p=v mod2,
    (a,b) --input d / output p-->
        (a+p, ((1+2p)v+p)/2).                        (A2)

The initial state is `(0,0)`. Its invariant is `0<=b<3^a`. Indeed
`0<=v<2m`; the even update is below m and the odd update is at most
`3m-1`. Once an odd step occurs, b is a unit modulo3: the new odd-step
value is2 modulo3 and later halvings preserve being a unit.

## 3. Complete classification: a tower of ternary unit cycles

**A3 (PROVED).** The reachable residual states are exactly

    {(0,0)} union {(a,b): a>=1, 0<b<3^a, 3 does not divide b}.    (A3)

At fixed a>=1, choose the unique input producing output p=0. Then

    d=b mod2,
    b -> (b+3^a d)/2 = 2^(-1)b mod3^a,               (A4)

where the right side means the least nonnegative residue. This is one
cycle through all `2*3^(a-1)` units modulo `3^a`.

For completeness, the order assertion is elementary. The order of2
modulo3 is2. The binomial expansion gives
`v3(4^(3^j)-1)=j+1`, by induction: cubing `1+3^t u`, with `3∤u` and
t>=1, raises the valuation exactly once. Thus4 has order `3^(a-1)`
modulo `3^a`, and2 has order `2*3^(a-1)`. There are exactly that many
units. Multiplying any one unit by successive inverse powers of2
therefore visits all of them.

The all-ones input prefix of length a reaches `(a,3^a-1)` by the exact
identity `T^a(2^a-1)=3^a-1`. Following (A4) then reaches every unit at
that level. This proves the reverse inclusion in (A3), as well as an
explicit reachability bound: each state at level a is reached by some
input prefix of length at most `a+2*3^(a-1)-1`.

The unique input producing output p=1 instead gives the level lift

    b -> 2^(-1)(3b+1) mod3^(a+1),                    (A5)

always a residue2 modulo3. Equations (A4)--(A5), with the input labels
restored using (A2), give the complete infinite minimal residual
automaton. This sharpens C10's selected infinite family to a complete
description of its states. It does not make the automaton finite.

There are exactly

    1+sum_(a=1)^A 2*3^(a-1)=3^A                     (A6)

states through odd-count level A. A literal identification with the
ternary grid is

    (a,b) -> beta=b/3^a in [0,1),   (0,0)->0.

These are precisely the finite ternary fractions in reduced form. Through
level A they are exactly `{j/3^A:0<=j<3^A}`. The value of a is recovered
from the reduced denominator, so this coordinate does not forget it.
This is an actual arithmetic map behind the ternary-tree count, not a
claim about Berggren ancestry or Pythagorean triples.

**Scope warning.** The cycle (A4) follows selected *output-even* edges.
Their input labels are `b mod2`, which usually vary. It is not an even
Collatz cycle on the positive integers and is not the fixed-source
zero-input path considered below.

## 4. Exact macro-carry cost and synchronization

Fix one parity cylinder from (A1). Computing its endpoint from the
unread tail is the affine multiplication `q -> m q+b`, with `m=3^a`
and `0<=b<m`. Its binary synchronous transducer has carry states

    c in {0,...,m-1}, initial c=b,
    output e=(m d+c) mod2,
    next carry c'=floor((m d+c)/2).                  (A7)

**A8 (PROVED).** This transducer is minimal and has exactly m states.
After k bits representing `0<=u<2^k`, its carry is

    c_k=floor((m u+b)/2^k).                          (A8)

For `2^k>=m`, these values run through every integer from0 through m-1:
the first value is0, the last is m-1, and adjacent values differ by at
most1. Thus all m states are reachable from every initial b. Distinct
states are distinguishable by zero continuation, which emits the
binary expansion of c. This is a residual-function lower bound for
the explicitly stated synchronous scan model, not a lower bound on
every algorithm or the one actual orbit.

There is also an exact reset census. A length-h input u sends the entire
carry bank to one state iff

    (m u mod2^h) <= 2^h-m.                           (A9)

This follows by comparing (A8) for c=0 and c=m-1. If `2^h<m`, no reset
exists. Since m is odd, multiplication by m permutes residues modulo
`2^h`; hence the number of length-h reset words is exactly

    max(0, 2^h-m+1).                                (A10)

The sharp reset threshold is `ceil(log2 m)=ceil(a log2 3)`, with the
empty word allowed when m=1. In particular that many zero input bits
flush every carry to0. There is no uniform reset length over all odd
counts a, despite the three states of one multiplication-by3 scan.

The two state counts `3^a` in (A8) and `2*3^(a-1)` in (A3) refer to
different clocks and machines. (A8) scans the unread binary tail of a
*fixed macrostep*. (A3) consumes another input bit and simultaneously
produces another Collatz parity. Conflating them destroys the exact
source interface.

## 5. Fixed integers, finite end markers, and a contraction that is not descent

Fix n>0 and H=bitlength(n). After its H input bits have been read, the
remaining input is identically zero. Thus in the residual automaton

    b_L=T^L(n),  a_L=#odd steps before L, L>=H,
    beta_L=b_L/3^(a_L).                              (A11)

The source remains the same n at every depth. Ordinary-integer support
termination is fully preserved: the input has a known finite end marker.
Nevertheless the *output* orbit still has unbounded possible duration.
For the zero-input transition, (A2) becomes

    beta_(L+1)=beta_L/2                              if b_L even,
    beta_(L+1)=beta_L/2+1/(2*3^(a_L+1))              if b_L odd.   (A12)

Since b_L>=1, the odd ratio is `1/2+1/(6b_L)<=2/3`; the even ratio is1/2.
Consequently

    0<beta_(H+k)<= (2/3)^k beta_H ->0.                (A13)

**A13 is unconditional, even if Collatz were false.** The normalized
state converges to0 for every fixed positive input. The reduced numerator
b_L, which is the actual integer, may grow while its denominator grows
faster. The real-valued coordinate is faithful but is not a well-founded
rank. This is a particularly cheap hostile to using an exact contracting
compactification as a termination proof.

This is pointwise contraction toward0, not a Banach contraction theorem:
the map has denominator-dependent jumps. Near a finite ternary fraction
with odd reduced numerator, distinct fractions of increasing denominator
have updates tending to beta/2, while its own update contains the positive
term in (A12). The end-marker information has not removed this arithmetic
boundary dependence.

Two controls prevent overreading state complexity. For the fixed source
n=1, after2k steps with k>=1, `(a,b)=(k,1)`: the full macro transducer has `3^k`
states while the actual endpoint remains1. For the changing sources
`n=2^k-1`, after k steps `(a,b)=(k,3^k-1)`: almost the entire carry
capacity is occupied. The second family does not exhibit one divergent
positive orbit; the first refutes universal-state-count as a measure of
fixed-orbit difficulty.

What remains is concrete. Collatz is equivalent to every positive
integer's zero-input path eventually having reduced numerator1 or2.
The whole path already satisfies (A13). An effective estimate must control
the numerator, or compensate the denominator depth, rather than prove
more real convergence of beta.

## 6. Connection ledger and the next genuinely different question

| Source -> target | Map and preserved predicate | Lost information / necessary sidecar | Hostile test |
|---|---|---|---|
| Dyadic source-prefix tree -> ternary unit tower | (A1)--(A5), exact residual output function | Different prefixes can share a residual; keep emitted prefix and source address | r=4 and5 at L=3 share (a,b)=(1,2), but have different past parities |
| Repeated three-state scans -> macro transducer | (A7), exact unread-tail endpoint | A fixed alphabet no longer suffices as a grows | all m states reached and separated by zero tails |
| Rule-30 block alphabet -> Collatz macro carries | transfer the exact-block resource argument | No common local evolution or prize predicate | n=1 has huge universal macro complexity and trivial orbit |
| Carry bank -> reset language | (A9), exact synchronization of one scan | Temporal Collatz updates are a different operation | reset depth grows as ceil(a log2 3) |
| Exact residual -> finite ternary fraction | beta=b/3^a, both numerator and depth recoverable | A real limit forgets relative numerator/denominator growth | (A13) holds without assuming Collatz |
| Fixed source -> zero-input ray | known end marker at H, all future endpoints retained | No output stopping bound follows from finite input | numerator descent is still the explicit target |

The tower gives a sharper question than another finite-state search:
can a boundary estimate distinguish the zero-input ray from the freely
chosen horizontal cycles using its *reduced numerator* and the forest's
ordered carry positions? A successful estimate must read the canonical
reduced numerator, respect (A2), and produce an actual smaller positive
numerator at a finite selected return. The scan reset
census alone cannot do this, because it fixes a before flushing while
the orbit can keep increasing a.

## 7. Reproduction and finite scope

Run:

    python 04-computation/experiments/nextforest_20260926_automata.py
    python -O 04-computation/experiments/nextforest_20260926_automata.py

The [stored output](nextforest_20260926_automata.out) records: all8191
source prefixes through length12 with four exact unread-tail controls
each; all2186 non-root residual states through odd depth7; all possible
initial carries for macro machines `a=0..6`; reset-word counts at four
near-threshold lengths (three at a=0); and16320 exact steps along the
zero-input tails of the fixed sources1..255. Every check uses an explicit
exception and remains active under `-O`. Universal claims have the
elementary proofs above; finite controls do not promote them to Collatz.
