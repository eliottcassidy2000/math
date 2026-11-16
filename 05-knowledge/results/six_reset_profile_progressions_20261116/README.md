# Six steps, phase lifts, and exact two-count targets

Two distinct tasks are solved here: a deterministic selection for
explicit tournament families, and a quantitative obstruction showing why
those families are not universal. The general Chan–Pak conjecture and
the rational quadratic six-cycle problem are not solved.

## 1. What is special about 63, and what is not

Starting at 0 and iterating D(x)=2x+1 gives D^n(0)=2^n−1. At n=6,
63=7·9; 7 appeared at n=3 and 3 at n=2. The case n=1 is the unit
1, so 6 is the first and only exceptional n≥2 for primitive prime
divisors, by [Bang–Zsigmondy](https://bartmichels.github.io/files/zsigmondy_en.pdf).
Yet six is present at the prime-power level:
ord_3(2)=2, ord_7(2)=3, ord_9(2)=6. Thus 63's square factor
carries a genuine clock that disappears on passing to the prime
support. Equivalently 2³+1=9=3² and Φ₆(2)=3. This is a precise
example of old factors carrying genuinely new phase data.

At c=−29/16 the signed quadratic orbit is

    −7/4 → +5/4 → −1/4 → −7/4.

The suggested +1/4 and −5/4 are preimages, not on the cycle:
they go to −7/4 and −1/4, respectively. For an actual six-cycle
x_i in characteristic zero, the three half-turn differences imply

    (x_0+x_3)(x_1+x_4)(x_2+x_5) = −1.

Indeed each difference is multiplied by its pair sum at the
next step, and after three steps its order reverses. Thus if a
six-cycle with good reduction at an odd prime collapses to
precisely the displayed period-three orbit, the three-cycle
multiplier 35/8 must reduce to −1. This forces p=43, since
35/8+1=43/8. The archive's fixed-parameter census and lift
theorem is stronger at that place: the 21 six-cycles modulo
43² do not lift with period six modulo 43³. Varying c is a
different problem.

[Stoll](https://arxiv.org/abs/0803.2836) excludes exact rational
period six conditional on BSD for a specific Jacobian. Over a
quadratic field there is a real six-step construction from the
same source: c=−71/48 and x=−1+√33/12. Its first three conjugate
coordinate pairs (a,b) with x=a+b√33 are (−1,1/12),
(−1/4,−1/6),(−1/2,1/12); the next three negate b.
Projection to the half-trace a gives an apparent three-cycle
but does not give a rational map a↦a²+c' with any constant c',
because a' reads 33b². This is an exact arithmetic counterpart
of a projection losing the next response, not an identification
with tournament Hamilton paths.

Archive trail: [September Zsigmondy, Catalan and critical orbit
synthesis](../collatz_mod6_20260917_zsigmondy_triad.md),
[period-six versus lifted phase](../quadratic_cycle_phase_20261015b.md),
[denominator signatures and the modulus-43 failure](../rational_quadratic_denominator_cycles_20261019/README.md),
[7 and 21 Mersenne clocks](../seven_twentyone_mersenne_openai_math_20261006.md),
[Mersenne switch and order-21 compression](../mersenne_switch_parity_f21_compression_20261006.md),
and [the Monster/Paley affine quotient](../session32_square_monster_20261023/paley/PALEY_MONSTER_BRIDGE.md).
The last supplies the genuine Paley7 count 189=3·63. In the
old row, 756=4·189=12·63=2²·3³·7 and 33=3·11. In particular
the unshifted first digit has the factorization
756·2^s−33=3(63·2^(s+2)−11); this is not a claim about
extracting a tournament by dividing by three.

## 2. New physical returns: five distinct exact counts, same full interface

[QUINTUPLES.md](QUINTUPLES.md) supplies literal nine-vertex raw
seeds for both exact rows

    167·2^t−75+{0,2,4,6,8}, t≥3,
    601·2^t−315+{0,2,4,6,8}, t≥2.

All five outputs in each row have the same full (g,E,O), not
just the same E. Any common fixed-root sink/near-sink word
therefore preserves all offsets. This repairs the loss that
made an ordinary sink unsafe on the earlier 756 row with
three unequal deletion counts. It is not universal placement.
The raw seeds and first two strong heights were independently
checked by subset DP; the recurrence proves all later heights.

## 3. A concrete exact profile that needs a richer block

From [the old 756 selector](../first25_spread_target_returns_20261116/TARGETS.md),
a nonempty right near-sink block B with two-count profile
(u,v)=(H(B),P_B(2)) gives

    N = 756·2^s v + Iu,  s≥2, I ∈ {−33,−31,−29}.

There is an exact bounded-order search for any supplied N:
v≥2u gives a finite bound on v, and v≥2^|B| bounds the
block order. This can certify a hit or absence in this
architecture. It is not a theorem that it always hits.

At the prime target 821861, choose s=2, I=−29. Literal six-vertex
outgoing masks (2,4,25,19,35,15) give (u,v)=(23,272).
There are two disjoint directed triangles, so no single vertex
can leave a transitive deletion. Direct DP and trial division
verify the profile and the primality. The old attachment law
gives the exact target 3024·272−29·23=821861. There is **no**
representation at that target if every nontrivial strong
component of B is one-apex, regardless of how many components
or what height. The proof is exhaustive but tiny: v≤273
since N≥(3024−33/2)v. Ordinal component profiles multiply
and singleton padding is a dyadic factor; after folding
that factor into s there are only nine candidate profile
products, including the identity. The replay checks all three
intercepts and every possible dyadic exponent.

Stronger: put N_s=756·2^s·272−667. If the one-apex block regime
hits this ray at a different height h≥2 with 2h−s≥3, s≥3,
it faces an impossibility. The bound u<138·2^(s−h)
follows from N_s≥(1512·2^h−33)u. For h≤s,
both leading terms are divisible by 756·2^h, but
|Iu+667|<max(667,4554·2^(s−h))<756·2^h
under 2h−s≥3. For h≥s instead use divisor 756·2^s;
the same strict bound is <4554<756·2^s. Thus Iu=−667,
forcing I=−29 and u=23. A one-apex product with this
prime count can have only one nontrivial atom. Every one-apex
atom at 23 has v=296=8·37, whereas the equation requires
v=272·2^(s−h), whose odd part is 17 whenever integral.
To stay in the old block regime on this ray, the selector
would have to drop below about half the prescribed height.

## 4. A power saving for arbitrary concatenations of one-apex blocks

For a strong one-apex tournament with transitive spine length m,
its address begins 1 and ends 0 (1 means the apex beats that
vertex). The odd cycles all contain the apex. Across any 1 at i
and 0 at j>i, put d=j−i. A cycle chooses an even number of the
d−1 intermediate vertices. This contributes 2^max(1,d−1) to H.
The path-colour odd-cycle formula at two colours gives contribution

    2^(m−d) (3^(d−1)+1)

to P(2), in addition to its empty-pack term 2^(m+1).
Thus just the end contacts give v≥L_m=2^(m+1)+2(3^(m−2)+1).
For H=23, even the elementary end contribution to H forces
m≤6, justifying the tiny profile classification used above.

There are 2^(m−2) strong words of length m. At θ=9/10,

    sum_(atom addresses) v^(−θ) ≤ Σ_(m≥2) 2^(m−2)L_m^(−θ) < 1.

The verifier upper-bounds the first eight terms by
(107+115+115+108+97+84+70+56)/1000; the remaining geometric
tail is bounded by (537/1000)(50/67)^8/(1−50/67), since
2^(−9/10)<537/1000 and 3^(9/10)>67/25. All comparison
checks are exact integer powers.

For an arbitrary ordered concatenation of these atoms,
the two-count profile multiplies. The sum of v^(−θ) over
all words therefore converges by a geometric series.
Transitive factors are absorbed into the dyadic height, which
has another convergent θ-sum. Since outputs are bounded below
by a constant times the dyadic height factor times v,
the total number of distinct targets up to X in this entire
architecture is O(X^(9/10)). This applies with no bound on
the number of one-apex components or their orders. The same
proof works for a fixed finite bank of such affine readers
with positive uniform lower bounds. This is much stronger
than saying there are not enough small individual blocks.

The next exact obligation is a target-driven source class beyond
one-apex components, or a seam whose joint count law is different.
The witness at 821861 shows how an actual two-triangle block can
bypass the obstruction; it does not establish that a richer class
hits all progressions.

Run python dynamics.py, python policy.py, and
python check_quintuples.py. No external packages are needed.
