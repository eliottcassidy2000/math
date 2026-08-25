# Sun minimality, odd functions, and the sixty-phase tail — 2026-08-24

## Truth surface

This session began with three deliberately separate questions:

1. **Anchor:** is `896315812331399` the least Sun `2-4-6-8` hole, and can all
   holes be classified?
2. **Niche:** what survives from the repository's odd-function calculus after
   centering even-rank binomial coefficients?
3. **Wildcard:** why does the consecutive-AP deficit have an exact 60-phase
   tail, and what is genuinely shared with triangular and Fibonacci clocks?

The status after the session is:

- **FINITE-EXACT:** THM-4026 proves that
  `N=896315812331399` is a hole and that both neighbours are represented.
- **CITED:** the public 2019 search covered every positive integer through
  `2*10^12`.
- **OPEN:** global leastness. The full interval
  `[2,000,000,000,001,896,315,812,331,398]` remains undischarged.
- **FINITE-EXACT + VERIFIED-EXACT:** THM-4040 proves that `N` is the least and
  unique hole in the one discovery class `459490 mod 1062347` through
  `1,001,999,999,999,999`.
- **OPEN:** all holes outside that class, all holes above that height, and the
  finiteness, infinitude, or density of the global hole set.
- **PROVED:** THM-4036 gives the global energy/support constraints, including
  average-scale-rich support and a three-higher-role support lower bound.
- **PROVED + VERIFIED-EXACT:** THM-4037 gives the centered parity and singular
  fibre laws.
- **PROVED + FINITE-EXACT:** THM-4038 gives the exact holonomic 60-phase law.
  It does not prove AP extremality or LRC(14).

The requested global minimum and global classification were therefore not
obtained. Replacing them by the selected-class statement is not a matter of
wording: the class sees only one integer in `1,062,347` and THM-4027 proves
that no fixed congruence class can be an exceptional tail.

## Inheritance pass

The closest proved mechanism for the Sun anchor was THM-4028's lattice-volume
main term, together with THM-4027's universal modular solubility. The canonical
hostile was THM-4026 itself: complete local support and a growing global mean
coexist with a bounded integral hole. The corrected near miss was the attempt
to turn sparse CRT factors into a necessary obstruction; universal local
solubility blocks that implication. The least-used sidecars were exact
polynomial differences and centered reflection.

For the AP wildcard, THM-4029 already supplied the twelve rational owners and
the affine last-track counters. The hostile predecessor was the sharp onset at
`m=12`; any tail statement beginning there is false. THM-4035 supplied the
least-used comparison: full Fibonacci state modulo ten and pointed triangular
state modulo thirty are both reversible 60-cycle addresses, but scalar shadows
have collisions.

## Live concept board

The board was kept at six objects:

1. the bounded positive representation box for `a(n)`;
2. complete modular period boxes;
3. the target CRT search class;
4. the energy/collision graph of polynomial atoms;
5. centered `C2` reflection and its odd derivative;
6. the AP rational-owner clock with its unbounded affine height.

Every useful connection preserved some of these coordinates and destroyed
others. The central lesson was that three superficially similar periodic
objects—modular representation boxes, the discovery CRT class, and the AP
60-clock—have different semantics. One is a complete local fibre, one is a
search prior, and one is an exact owner state.

## Result 1: exact finite classification of the discovery class

Put

```text
M = 11*13*17*19*23 = 1,062,347,
R = 459,490,
H = 1,001,999,999,999,999.
```

The repaired block engine partitions `[1,H]` into nine contiguous intervals,
enumerates every higher triple

```text
C(x,4)+C(y,6)+C(z,8) < block high,
```

and solves the remaining triangular congruence in exact index lanes of period
`2M`. Exact `uint128` repairs make the floating square-root seeds
non-decisive. Across all `943,194,644` selected targets, the zero vector by
block is

```text
0,0,0,0,0,0,0,0,1,
```

and the final zero is exactly `N`. An independent implementation checked
small Cartesian universes, interval arithmetic, the nine-block partition,
the last three higher boxes, saturation safety, and the target address.

What is preserved: exact target height, role labels, every higher triple, and
the triangular lift. What is destroyed: every off-class integer. The latter
loss is precisely why this result cannot establish global leastness.

## Result 2: energy turns abundance into logarithmically full support

For a degree-`d` binomial atom, a nonzero difference has at most

```text
tau(d!*|h|)
```

oriented preimages. The proof factors the polynomial difference by `u-v` and
then uses strict monotonicity in the gap `u-v`; the zero difference is kept as
an exact diagonal rather than forced into the divisor estimate. Peeling the
roles in order `2,4,6,8` gives

```text
sum_(n<=X) a(n)^2 <<_epsilon X^(13/12+epsilon).
```

Combined with total mass `asymp X^(25/24)`, this gives at least
`X^(1-epsilon)` represented targets. A second mass split shows the same lower
bound after requiring

```text
a(n) >= theta * (mean multiplicity),   0<theta<1.
```

The classwise version uses each fixed residue class's own positive local
mean. This is stronger than positivity and weaker than positive density.

The discarded packets give a useful new ladder. The support of
`C(y,6)+C(z,8)` has exact logarithmic exponent `7/24`, while

```text
#{h<=X:h=C(x,4)+C(y,6)+C(z,8)} >>_epsilon X^(1/2-epsilon).
```

This is an unconditional global support result for the higher-role sumset.
It does not control short intervals or maximum gaps, so it cannot turn the
triangular role into an eventual-coverage theorem without a new spacing or
lower-tail input.

## Result 3: the lawful odd-function bridge

The source is LEM-004's calculus of even and odd functions under an
involution. The target map is

```text
iota_k(n)=k-1-n,
C(iota_k(n),k)=(-1)^k C(n,k),
Delta C(n,k)=C(n,k-1).
```

Thus an even-rank binomial atom is centered-even and its increment is
centered-odd. This is the precise connection to prior odd-function work; the
atoms themselves are not tournament orientations.

On a complete exact-period box, every even modulus gives a free coordinate
reflection action, hence full `2^s` divisibility of every fibre count. At an
odd modulus each coordinate has one fixed point, leaving exactly one odd
target. For the `2,4,6,8` packet it is

```text
-3453/32768.
```

At `N`, the numerator factors as

```text
32768*N+3453 = 5*7*42667*199447*98610539.
```

Differentiating the centered-even atoms gives centered-odd critical equations.
Their degree-24 additive eliminant exposes eight exact singular/nonlifting
primes. The all-centered resultant branch is the same parity fingerprint,
which is the surprising exact overlap.

The first failed implication is equally important. Centered reflection does
not preserve the positive chamber, bounded height, or a common integral lift.
Complete modular parity and singularity therefore cannot detect the exact
hole. THM-4026 is the hostile witness.

For the triangular role, the missing sidecar has an exact name. Writing
`b=2w-1` gives `8*C(w,2)+1=b^2`, while centered reflection sends `b` to
`-b`. Hence all Sun holes are exactly the higher-sumset values for which
every translate `8(n-h)+1` avoids odd squares `b^2`, `b>=3`. The square
quotient preserves the even atom and destroys the sign; the positive chamber
is precisely the choice `b>=3`. This is a lawful odd-function classifier
criterion, but it does not supply the missing global square-intersection
estimate.

## Result 4: sixty is a clock, not a Sturmian cause

After the THM-4029 owner reduction, for `n>=12`,

```text
D(n+1)=(1/7) sum_(c=0)^5 A_c(n mod 60)/(n-c).
```

The six owner-denominator packets have exact periods `1,2,3,4,5,6`; their
least common multiple is 60. The full coefficient vector has 60 distinct
phases, but its leading moment is constant. Phase first appears at order
`1/n^2`.

The Fibonacci denominators `1,2,3,5` produce only period 30. Even after adding
the `q=6` packet the period remains 30; `q=4` contributes the missing 2-adic
factor. A genuinely golden Sturmian slope covers all seven sectors by `m=11`
and is not a persistent tail owner. These are decisive hostile probes against
a Fibonacci-cause story.

THM-4035 nevertheless proves an exact address connection:

```text
phase r mod 60
  <-> (F_r,F_(r+1)) mod 10
  <-> pointed triangular state mod 30.
```

The first two phase moments also inject the 60 phases. These conjugacies
preserve successor and address, but destroy owners and affine height. The
actual cause remains rational-owner localization.

Multiplying by

```text
Q_6(n)=n(n-1)...(n-5)=6!*C(n,6)
```

clears all six active poles. It is the unique monic common desingularizer of
least degree. The cleared sequence is a degree-five quasipolynomial of exact
period 60, yielding an explicit order-360 polynomial recurrence. Its ordinary
generating function is D-finite but nonalgebraic because of its logarithmic
singularity.

The all-prime refinement THM-4042 makes the causal statement sharper. For
`2<=q<=P-2`, the exact owner-word period is

```text
product over ell^e || q of ell^min(e,v_ell(P-1)+1),
```

with top boundary period `rad(P-1)`. Hence the global phase has a closed
prime-adic product formula. It equals 60 at `P=7`, but collapses from the
naive denominator lcm by factors `2,6,6` at `P=5,11,17`. The Fibonacci and
triangular states remain exact addresses for the exceptional `P=7` clock;
the valuation of `P-1` and the image set `(P^(-1)-1)U_q` are the general
cause.

Block summation of each owner track gives an exact triangular-number formula.
That bridge preserves phase and total arrival time, but destroys reciprocal
costs, sides, max-min selectors, and missing-sector labels. It is useful as a
clock identity, not as a second derivation of the AP deficit.

## Frontier after the session

The most valuable next Sun move is not another selected congruence scan. It is
one of:

1. a short-interval or translated energy estimate for the higher-role sumset;
2. a lower-tail/singular-series theorem strong enough to convert average-scale
   mass into pointwise triangular intersection;
3. an exact rank-labelled Pascal carry invariant that retains the positive
   chamber and bounded height; or
4. a genuinely global block census of the interval above `2*10^12`, which is
   computationally far larger than the one-class result.

For the AP/LRC lane, the next test is to Fourier-decompose the six denominator
packets on `C_60` and ask whether any character statistic survives when the
consecutive AP is replaced by a sparse runner set. The required sidecars are
owner, side, missing sector, gap geometry, and unbounded height. A clock match
without those data is only an address coincidence.

THM-4044 is now **PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED**.  Its
confluent `C_60` observer has exact kernel `((P^60-1)^k)`, and on the forced
pure residual ideal the kernel is `P^2(P^60-1)^kK[P]`.  Thus the proposed
cross-frontier hostile is settled: a fixed finite clock depth loses the
mandatory `[P^2]R` coordinate at degree `60k+2` unless the second boundary
Hasse jet or a degree cap is retained.  The map from an AP phase to the same
cyclic evaluation address preserves only `C_60` translation and destroys
owner words, poles, affine height, and the Keller constraints.  THM-4042's
clock collapses at other primes make a bare period match especially weak
evidence.

## Public anchors

- Zhi-Wei Sun's formulation and the reported all-integer verification through
  `2*10^12`: https://mathoverflow.net/questions/323541/
- OEIS representation convention and historical checks:
  https://oeis.org/A306477
- Discovery transcript for the explicit counterexample and selected CRT
  blocks:
  https://gist.github.com/tadamcz/0c578c8b2b3fb92fe8584bc0725187e3
