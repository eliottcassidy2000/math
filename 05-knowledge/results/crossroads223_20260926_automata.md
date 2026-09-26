# 223, return languages, and the coordinate lost by finite residue lifts

**Status:** PROVED elementary lemmas + FINITE-EXACT computations; research note,
not canon. **Date:** 2026-09-26. **Scope:** finite residue lifts of half-Collatz
parity languages, their cycle means, and integer periodic-source obstructions.
No Collatz convergence claim follows.

## Inheritance and concept board

- Closest proved mechanism: THM-4485,
  [periodic-edit-price-feedback-sets](../../01-canon/theorems/THM-4485-periodic-edit-price-feedback-sets.md):
  bounded-lookahead periodic edits must hit expanding cycles of the binary
  de Bruijn graph. THM-4486,
  [min-max-cycle-density-game](../../01-canon/theorems/THM-4486-min-max-cycle-density-game.md),
  keeps worst-cycle payoff distinct from stationary/average density.
- Canonical hostile: arbitrarily long all-odd prefixes occur at positive
  integers. In fact `n=223*2^H-1` has first `H` bits all one and every one of
  these states has residue `-1 mod 223`. A fixed auxiliary residue does not
  establish a periodic integer orbit or rule out positive realization.
- Corrected near miss: a residue word that fails to close after one turn can
  close after several turns. That defeats its use as a pruning rule for a
  finite mean-payoff graph. It does **not** defeat the one-turn divisibility
  obstruction for an integer or 223-integral rational periodic source; these
  are different assertions. A rational source whose denominator contains
  223 cannot be reduced into F_223. Sections 3 and 5 separate the scopes.
- Little-used sidecar: the integer affine translation, retained before
  reduction modulo an odd prime. The earlier
  [automata note](crossroads_20260926_automata.md) showed that equal-length,
  equal-weight blocks can have identical multiplicative clocks but different
  rational/irrational sources. Here the lost coordinate is visible already
  for a single periodic word and its powers.

Live board: (1) de Bruijn expanding-cycle feedback; (2) mean payoff;
(3) finite-field affine monodromy; (4) Paley orientation;
(5) primitive words versus lifted returns; (6) integer carry/denominator.
The wildcard integer 223 becomes useful only through these exact maps.

Targeted repo recovery found 223 among arithmetic examples, including the
`14a-1` progression at `a=16`, a `12/q` perfect-number-related packet, and an
eventually periodic Collatz reconstruction control `223/45`. None supplies
a new Collatz premise here. The calculations below are derived directly;
no literature novelty is claimed.

## 1. Exact maps and their losses

For a chronological binary word `w=(b_0,...,b_(r-1))`, put `e=sum b_i`.
The half-Collatz branch composition is

    F_w(x) = (A x + C)/B,
    A=3^e, B=2^r,
    C=sum over j=0,...,e-1 of 2^(i_j) 3^(e-1-j),

where `i_j` are the positions of the ones. The formula remains valid for the
constant signed map `(3x+sigma)/2` after multiplying C by `sigma`.

| Source | Target and map | Preserved | Lost; needed sidecar |
|---|---|---|---|
| Ordered parity word | `(r,e)` | multiplicative slope `3^e/2^r` | order and carry C |
| Integer affine map | `x -> a x+c` over F_223, with `(a,c)=(A/B,C/B)` | finite residue transitions | height, positivity, divisibility at higher powers; retain integer `(A,B,C)` |
| Base directed graph | complete residue lift `(vertex,x)` | every edge weight and maximum cycle mean | which closed lifts are integer periodic orbits; retain rational fixed-point quotient |
| Affine permutation | action on Paley arrows | orientation sign `(-1)^e` | carry, locations, integer parity; orientation is not a descent observable |

Every finite parity word has a unique compatible class modulo `2^r`.
CRT combines it with every residue modulo any odd M. Thus every finite
parity/residue start admitted by the complete lift is realized by positive
integers. The realization depends on the horizon; this gives no single
positive integer realizing a prescribed infinite word.

## 2. What is special about 223, exactly?

Trial division proves 223 prime; direct multiplication and order tests give

    2^37-1 = 223 * 616318177,
    ord_223(2)=37, ord_223(3)=222, 2=3^180 mod 223.

The other Mersenne factor need not be prime for any argument here.
The branch maps over F_223 are permutations and

    a_w = 3^(e-180r).

Consequently `a_w=1` iff `e-180r` is divisible by 222. The first nonempty
resonance with `0<=e<=r` is `(r,e)=(21,6)`. The first expanding one is
`(26,18)`. The only expanding resonances through length 49 are `(26,18)`
and `(31,30)`. The order 37 is a residue clock: the word `0^37` acts as
the identity modulo 223, whereas its actual action divides by `2^37`.

Since 223 is 3 modulo 4, the rule `x -> y iff y-x is a nonzero square`
defines the Paley tournament on F_223. The quadratic character satisfies
`chi(2)=+1`, `chi(3)=-1`. A branch therefore preserves every arrow when
even and reverses every arrow when odd, for either constant sign.
Translations cancel in differences. A word's orientation gauge is exactly
`(-1)^e`; it contains no further carry information. This is a precise
Paley connection and a decisive stopping reason for treating it as a new
Lyapunov coordinate: it recovers only the parity of an already known count.

## 3. Complete finite lifts preserve the maximum cycle mean

**PROVED.** Let G be a finite directed graph with real edge weights and at
least one cycle. Replace each vertex v by a nonempty finite fiber X_v.
For each edge `v -> u`, choose a total function `X_v -> X_u`, and lift the
edge from every fiber point, preserving its weight. Then the maximum mean
weight of a directed cycle in the lifted graph equals that in G.

**Proof.** A lifted cycle projects to a closed base walk. Decomposing that
walk into base cycles bounds its average by the maximum base cycle mean.
Conversely, take a maximizing base cycle. Its fiber return map is a total
self-map of a nonempty finite set, so it has a periodic point. Repeating
the base cycle by that point's period gives a lifted closed walk with the
same average. A cycle decomposition of that closed walk contains a cycle
with at least that average; the first inequality gives equality. QED.

This applies to every complete odd-modulus Collatz residue lift, including
moduli divisible by 3: permutation hypotheses are unnecessary. Iterating
this construction, or using any finite product of such moduli, still
cannot lower the maximum cycle mean. A fixed sign strategy also has this
property if its choices are merely lifted from the base graph.

**Boundary:** allowing the strategy to make new sign/edit choices depending
on the added residue changes the strategy class. The lemma does not prove
equality of the resulting optimized game values. Deleting states by an
additional justified arithmetic condition also leaves its hypotheses.
Nor does equality of maximum mean imply equality of the minimum number of
vertices needed to hit expanding cycles: lifted cycles can have very
different lengths and intersections.

**Two exact extensions.** For vector-valued edge observables, the convex
hull of all cycle-average vectors is also unchanged. Every projected
closed walk is a convex combination of base cycle averages, and every
simple base cycle has a periodic fiber orbit whose lifted cycle has the
same average vector. Thus adding several additive statistics at once
cannot evade the obstruction. This does not cover observables depending
on the new fiber coordinate.

Likewise, give each base path any nonnegative weight W, inherited by its
lifts. A path starting at v has exactly `|X_v|` lifts. If `m_f` and `M_f`
are the smallest and largest fiber sizes, its total weighted path counts
satisfy `m_f Z_L <= Z_L(lift) <= M_f Z_L`. Hence their exponential growth
rates agree whenever they exist, including entropy and exponential weights
of additive observables. A complete finite residue lift therefore supplies
no new exponential rarity of base-word-defined events. Requiring a
specified residue return, or changing allowed choices using the residue,
is an additional condition and must be analyzed separately.

## 4. Primitive resonance words: an exact finite test

For an affine permutation `f(x)=a x+c` over F_p:

- If `a!=1`, it has one fixed point and `(p-1)/ord_p(a)` other cycles,
  each of length `ord_p(a)`: translate by its unique fixed point.
- If `a=1,c=0`, it is the identity.
- If `a=1,c!=0`, it is one p-cycle: a nonzero translation.

Two independent exact counts, a chronological carry DP and a
denominator-cleared formula over combinations, agree on all 223 bins:

| Word length, ones | All words | C=0 mod 223 | Range in each nonzero carry bin |
|---|---:|---:|---:|
| 21, 6 | 54,264 | 231 | 212--272 |
| 26, 18 | 1,562,275 | 7,735 | 6,741--7,266 |

The first expanding resonance illustrates exactly why primitive words
and primitive lifted cycles must be distinguished. For length 26, weight
18, an imprimitive word must be a square of a length-13, weight-9 word.
There are `binom(13,9)=715` such words. Every 13-block has multiplier -1
modulo 223, so its square is the identity. Thus the primitive length-26
words number 1,561,560, of which 7,020 have zero carry. Under rotation this
gives 270 zero-carry primitive necklaces and 59,790 nonzero ones.
Zero carry is rotation invariant because rotating the branches conjugates
the affine return map, and conjugation preserves identity/translation.

At de Bruijn order 25 every primitive 26-necklace is a simple base cycle:
two equal 25-windows would make the full rotated words equal, since their
remaining bit is forced by the fixed total weight. Each zero-carry base
cycle lifts to 223 cycles of length 26. Each nonzero-carry base cycle
lifts to one simple cycle of length `26*223=5798`. All have the same
expanding slope per base step and odd density `18/26`.

**FINITE-EXACT consequence:** a one-turn closure filter removes 59,790 of
these primitive base necklaces from an integer-cycle search, but removes
none of their expanding closed walks from the complete finite-field lift.
These counts concern a specified family, not all cycles in that graph.

## 5. The carry quotient explains the apparent contradiction

For a nonempty word, A and B are different. Its characteristic-zero
rational fixed point is `x_w=C/(B-A)`. Repeating the word d times gives

    A_d=A^d, B_d=B^d,
    C_d=C*(A^(d-1)+A^(d-2)B+...+B^(d-1)),
    C_d/(B_d-A_d) = C/(B-A).

**PROVED:** word repetition cannot repair a nonintegral rational source.
Moreover, a rational point periodic under F_w is already fixed under F_w:
`F_w^d(x)-x = ((A/B)^d-1)(x-x_w)`, with nonzero multiplier. This is why
the necessary congruence `C=0 mod p` when `A=B mod p` is valid for integer
or p-integral rational periodic-source searches. It is not necessary for
arbitrary rational sources: their denominators may contain p. The precise
integrality requirement is `v_p(C)>=v_p(B-A)`. For nonzero C, the difference

    v_p(C_d)-v_p(B_d-A_d) = v_p(C)-v_p(B-A)

is invariant. Keeping only whether a residue orbit returns loses this
relative valuation.

Concrete hostile: `w=1^18 0^8` has

    A=387420489, B=67108864, C=387158345,
    C mod 223=17, v_223(B-A)=1,
    x_w=-77431669/64062325.

Its F_223 action is a nonzero translation and closes after 223 repetitions.
The integer numerator/denominator valuations change from `(0,1)` to
`(1,2)`: the source remains exactly the same nonintegral rational number.
In fact no expanding ordinary `+1` word can give a *positive* fixed point,
because C is positive and B-A negative. Expanding de Bruijn cycles are
finite-horizon adversaries, not candidate positive periodic orbits.

## 6. What to pursue, and what to stop

The useful new connection to THM-4485/4486 is a boundary theorem: complete
finite arithmetic memory cannot erase worst expanding cycles, although it
can reorganize them and affect feedback-set geometry. The integer carry
quotient separates this finite graph question from rational-cycle closure.
A residue-enhanced sign strategy remains a legitimate different question;
the cheapest meaningful test compares the optimized game values while
explicitly permitting sign choices to depend on the new coordinate.

A more targeted open question is whether relative p-adic carry valuations
along a *nested* sequence of prefixes yield a height constraint strong
enough to distinguish one positive integer from separate positive witnesses
for each finite horizon. The finite-lift theorem prohibits replacing this
with any fixed collection of residue states. No sign-definite invariant or
uniform bound has been obtained here. The all-odd family is the hostile
control for any proposed finite-horizon version.

The Paley orientation idea stops at the exact gauge `(-1)^e`. The Mersenne
clock stops at order 37. The return-language idea advances to the complete
finite-lift lemma, primitive/resonant counts, and the exact missing quotient.

## Reproduction

    python 04-computation/experiments/crossroads223_20260926_automata.py

The standard-library probe compares independent complete carry histograms,
checks all Paley arrows for both branches and both constant signs, verifies
the affine cycle classification for every word of lengths 1--8 over primes
7, 31, 223, checks 3,048 exact word-power fixed-point identities, and tests
positive finite-horizon and genuine signed-cycle controls. Explicit check
functions remain active under `python -O`. The adjacent `.out` records the
deterministic result. No randomized tests or unbounded-orbit inference.

## 7. Independent audit of the sine supersolution

The companion [Robin note](crossroads223_20260926_robin.md) was audited
independently on 2026-09-26. No mathematical gap was found in the proposed
all-real-state bound `N_m(L)/2^L <= e D^L exp[-A L/(m+4)^2]`.
The audit checked the exact identities
`D*c*lambda^d=D*d*lambda^(-c)=1/2`, phase cancellation in the interior,
the partially absorbing top factor `2d`, and the endpoint comparison
including the special initial state zero. The top estimate has a positive
coarse margin already from `1-10/216-39/44>0`.
The logarithmic sine expansion has leading bracket
`1-c^3-d^3=3cd` and all subsequent brackets positive, so the coefficient
`A=pi^2 cd/2` has the asserted direction without an asymptotic error.

In the inherited private multi-scale inequality the shifted denominator
is `m+5`, which produces exactly the factor `3^6` after minimizing
`x log3 + AL/x^2`. Counting the scale terms is safely absorbed by
`2(1+e*3^6)(L+1)^2`; the explicit argument applies for integer `L>=8`.
The matching lower bound retains THM-4480's CITED Mogul'skii status.
This audit does not upgrade the private conclusion to a globally
consistent pairing construction, a polynomial comparison with the peak
price, or the full uniform Robin comparison.

## 8. A small q=7 adversary repair that reverses its initial signal

THM-4486 leaves `rho*(7,k)` open beyond the computed levels. Its negative
integer adversary selects the top lift at every pair modulo `H=2^(k-1)`.
A tempting repair is to replace the representatives `-1,...,-B` by
`H-1,...,H-B`, leaving the other top lifts unchanged. This prevents the
small free cycle from closing when `B>=1`.

The new probe computes the exact value of each fixed adversary: Min can
choose either sign. In the resulting graph on H pair indices, take the
strongly connected components having no outgoing edges. Within each such
component, the value is its minimum cycle mean, because Min can reach
every cycle from every vertex there. The value from Max's best starting
point is the maximum of those component values: any other start can reach
a closed component. Each minimum is certified by an integer potential on
all of its edges and a zero-weight cycle, not by a floating approximation.
A path-length computation proposes the rational value; failure of either
exact certificate aborts the probe. Memory is linear in H; the largest
tested graph has 4,096 pair states, with a 2,571-state closed component.

| k | B=0 | B=1 | B=3 | B=7 |
|---|---:|---:|---:|---:|
| 5 | 1/3 | 3/7 | 1/5 | 1/3 |
| 7 | 1/3 | 2/9 | 3/10 | 1/4 |
| 9 | 1/3 | 4/15 | 9/28 | 1/4 |
| 11 | 1/3 | 8/27 | 5/17 | 11/37 |

Additionally, `k=13,B=1` gives `12/41`. These are FINITE-EXACT values for
the specified adversaries, not the optimized values `rho*(7,k)`.

The first positive signal, `1/3 -> 3/7` at `k=5,B=1`, therefore does not
survive a larger level. The simplest hostile is the certified `k=7,B=1`
cycle, written using signed representatives modulo 128:

    -32 -> -16 -> -8 -> -4 -> -2 -> 63 -> -36 -> -18 -> -9 -> -32.

Only 63 and -9 are odd. The escape replaces the even transition
`-2 -> -1` with `-2 -> 63`; then Min uses the minus sign at the two odd
states. It buys an extra long even return and lowers the density to 2/9.
All remaining steps are exactly the chosen lift transitions. The graph
has one closed component, and the exact lower potential plus this cycle
prove that 2/9 is the fixed adversary's value from its best starting point.

**Stopping reason:** moving a small recurrent obstruction to height of
order H changes the cost of the whole return word. Breaking the free
cycle locally is not sufficient to improve a worst-cycle lower bound.
A useful next hypothesis must constrain those new return words or charge
the jump with a level-dependent potential. This also explains why merely
combining a height potential with the finite-residue orientation gauge
from section 2 adds no certificate: the gauge only records odd-count
parity and cannot pay for the long even excursion. No improved all-level
q=7 bound or finite provable q=7 strategy was obtained in this round.
