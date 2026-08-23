# Mahler `3/2`, LRC(14), prefix trees, and ABC radical sidecars

**Session:** root, 2026-08-23
**Portfolio:** Anchor = LRC(14); Niche = Mahler safe-prefix structure;
Wildcard = ABC/IUT and natural-number schedulers.
**Truth boundary:** the exact finite and all-level results are promoted in
[THM-3848](../01-canon/theorems/THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation.md).
Mahler's problem and LRC(14) remain **OPEN**.  Every ABC statement below is
**CONDITIONAL**; IUT-to-ABC remains a cited disputed claim, not a dependency.

## 1. Inheritance pass

The closest proved mechanisms were:

- [THM-2228](../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md):
  exact separation of strict real tails from ordinary `2`-adic residue
  stabilization;
- [THM-2352](../01-canon/theorems/THM-2352-q-adic-prefix-residue-collision-spectrum.md):
  every indexed plateau abscissa occurs in every unrestricted finite carry
  cylinder, while termination is dense and Haar-null;
- [THM-2047](../01-canon/theorems/THM-2047-phase-height-toric-arrangement-for-lrc.md)
  and [THM-2050](../01-canon/theorems/THM-2050-period14-top-germs-do-not-determine-global-loneliness.md):
  a complete local phase carrier still needs the global exit sidecar;
- [THM-3743](../01-canon/theorems/THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction.md):
  every hypothetical LRC(14) counterexample has one exact relation of
  `l1<=356`, but the relation loses owner, phase, sign partition, and arrival;
- [THM-3833](../01-canon/theorems/THM-3833-abc-conditional-cube-radical-and-hyperbolic-power-finiteness.md):
  ABC is usable only after presenting an actual primitive additive packet;
  radicals forget valuation depth, slots, signs, and common scale.

The canonical hostile examples were already pointing in the right direction:

```text
(100)^infinity: strict-safe real tails, nonterminal 2-adic state;
A=1:              ordinary stabilization, unsafe real tails;
greedy equality:  every finite strict prefix, infinite equality boundary;
AP13 vs lift:     identical local germs, different global maxima;
AP13:             many short relations, safely lonely.
```

The least-used relevant sidecar was the exact scaling from a finite rational
power prefix into a torus row.  Making that scaling explicit, rather than
starting from carry words, unlocked the session.

## 2. Live concept board

Five objects stayed live and were compared after each pull:

1. the mixed-power vector `V_N=(q^N,q^(N-1)p,...,p^N)`;
2. the one-sided torus prefix set `H_N`;
3. the strict safe carry-word language and its weak closure;
4. the centered lonely-runner evaluator;
5. ABC's radical evaluator on genuine carry packets.

The key discipline was to avoid calling all five “the same tree.”  They share
addresses only after naming the map, and their branching laws differ.

## 3. The normalized prefix torus is exactly solvable

For `t=xi/q^N`,

```text
{xi(p/q)^n}={q^(N-n)p^n t},             0<=n<=N.       (1)
```

This makes the finite Mahler-style prefix a standard torus intersection.  Its
common wall grid has denominator

```text
Delta_N=p^N q^(N+1).                                    (2)
```

The pullback `t -> qt` preserves every old phase.  Among its `q` inverse
sheets, the appended `p^(N+1)` phase selects exactly one sheet in the
half-open interval `[0,1/q)`.  Therefore the set-theoretic point lift is
unique and the measure is `q^(-(N+1))`.

The selected sheet changes as a point crosses internal walls.  Refining at
the next common grid shows that each safe *atom* has exactly `p` children.
Thus:

```text
point sheet count:        one selected from q;
open wall-atom branching: p;
level-N safe atom count:  p^N.                          (3)
```

This distinction resolves what initially looked like a binary-versus-ternary
contradiction in base `3/2`.

The child remainder in `{0,...,p-1}` is a reversible path digit.  Shortlex
turns it into a natural number, exactly in the spirit of treating the `n`th
odd square as address `n`.  But the ordinal is a scheduler.  The selected
inverse sheet, level, wall map, and target evaluator remain sidecars.

## 4. A full short relation basis can be very safe

The adjacent rows

```text
p e_n-q e_(n+1)                                         (4)
```

generate the complete integer annihilator of `V_N`.  Their maximal minors
are `q^N,pq^(N-1),...,p^N`, whose gcd is one.  Hence the relation lattice is
saturated, not merely a visible sublattice.

Nevertheless the LRC maximum is exactly

```text
M(V_N)=floor((p+q)/2)/(p+q),             N>=1.         (5)
```

The first adjacent pair supplies the upper bound.  Modulo `p+q`, the congruence
`p=-q` supplies a time at which every coordinate alternates between the two
opposite residues nearest the half-circle.

For `p/q=3/2`, every nontrivial level has `M=2/5`.  At level `N`, there are
exactly two maximizing times `+/-b_N/5`, where

```text
2^(N-1)b_N=1 mod 5,
b_N=1,3,4,2,1,3,4,2,... .                              (6)
```

The level-12 row is the exact fourteen-runner interface:

```text
4096,6144,9216,13824,20736,31104,46656,
69984,104976,157464,236196,354294,531441.              (7)
```

It is primitive, has thirteen distinct speeds, twelve independent norm-five
relations, exactly two maximizers, and `M=2/5`.  It misses every multiple of
five, so it is far outside the live covering hard core.

This changes the LRC concept board in a useful negative direction:

> Short relations, a complete relation basis, a recursive atom tree, and a
> tiny arithmetic alphabet can coexist with a very large lonely margin.

THM-3743's short relation is a necessary structural reduction, not a signal
that the row is near `1/14`.

## 5. The exact destroyed coordinate is phase orientation

At either maximizing time in (6), the ordered phases alternate

```text
2/5,3/5,2/5,3/5,... .                                  (8)
```

The centered LRC distance identifies `2/5` and `3/5`.  Mahler's predicate
does not: `3/5` lies outside `[0,1/2)`.  The failed transfer is therefore not
vague “global information”; it is the explicit phase-side bit.

The corresponding formal carry address is `(01)^infinity`, with

```text
Phi=-2/5,              Y_even=4/5,       Y_odd=6/5.    (9)
```

The canonical residues of `-2/5` never stabilize, and `Y_odd>1`.  The
mod-five LRC extremizer therefore fails both Mahler gates.  The equality
`Phi+Y_even/2=0` is only algebra; THM-2228's reconstruction theorem cannot be
invoked without its safety and stabilization hypotheses.

## 6. The incoming safe-language signal extends to a theorem

An agent proposed counting strict safe carry prefixes.  Independent hostile
review found the exact renewal mechanism rather than merely confirming the
first counts.

Let `d` be the greedy equality expansion

```text
1=sum_(k>=1)d_k(2/3)^k,
d=101000001001001010000000001000000100... .            (10)
```

If `P_m` is the length-`m` language satisfying every strict truncated suffix
inequality, then first disagreement with `d` gives

```text
a_0=1,
a_m=1+sum_(k=1)^m d_k a_(m-k),
A(z)=1/((1-z)(1-D(z))).                                (11)
```

The counts start

```text
1,2,3,5,8,12,18,27,40,60,90,134,201,302,452,... .    (12)
```

The inverse limit of the strict finite languages is not the strict infinite
set.  It is the weak closed shift `K` with every suffix tail at most one.
The difference is exactly the countable set of sequences eventually equal to
`d`.  This boundary is finite-prefix invisible.

The greedy remainder is always a strictly interior odd dyadic.  An eventually
periodic tail would also give it an odd denominator dividing `3^L-2^L`, an
impossibility.  Thus `d` is not eventually periodic.  Its shifted follower
sets are all distinct, so `K` is nonsofic.  Since `D(2/3)=1`,

```text
a_m~1.5510451884...(3/2)^m,
h_top(K)=log(3/2),
dim_H(K)=log_2(3/2),
mu_B(union_(w in P_m)[w])~1.5510451884...(3/4)^m.      (13)
```

Here `mu_B` is fair-Bernoulli cylinder measure.  This is a theorem about the
closed formal safe-tail shift.  The actual Mahler
language still intersects strict safety with ordinary stabilization and may
be empty.  It would be incorrect to advertise “the Z-language is nonsofic.”

There are now three exact but nonidentical carrier profiles:

```text
mixed-power torus safe atoms:       3^m;
formal safe carry words:            a_m~C(3/2)^m;
ordinary terminal residues:         dense but Haar-null globally.          (14)
```

The next operation should be their finite-exact fibre product, not an entropy
identification.

## 7. ABC reaches radical complexity, not the terminal clock

At an odd ceiling step `a -> b=(3a+1)/2`, the primitive packet

```text
3a+1=2b                                                   (15)
```

lets ABC act.  Conditionally, for every `epsilon>0`,

```text
rad(ab)>=(2b/K_epsilon)^(1/(1+epsilon))/6.              (16)
```

Every positive ceiling orbit has infinitely many odd states, so its odd-step
radical exponent has conditional liminf at least one.  But `a_0=2^H` has `H`
consecutive even steps with radical exactly six.  ABC sees the entry packets,
not the length of a smooth zero-carry corridor.

The denominator-19 family makes the stopping boundary even sharper.  For
positive multiples `m` of 18,

```text
alpha_m=9*2^m/19
```

has `m+1` safe phases cycling through `9/19,4/19,6/19`.  The identities

```text
19D_m+1=2^m,
2^m+19Q_m=3^m                                           (17)
```

conditionally force both `D_m` and `Q_m` to be logarithmically
radical-saturated.  Yet the family exists for every horizon; its initial
integer changes.  Radical saturation is compatible with arbitrary finite
safety.

IUT contributes no additional proved arrow.  If the disputed IUT-to-ABC chain
were assumed, it would feed exactly the same ABC consumer and stop at the same
hostiles.

## 8. Radical twins refute an LRC shortcut

Let `P` be the product of primes through 121 and `L=lcm(1,...,121)`.  The
thirteen-speed rows

```text
{1,...,11,13,84P},
{1,...,11,13,84L}                                      (18)
```

have identical componentwise radical vectors.  The last speeds are
respectively `55` and `0 mod 121`.  At `t=10/121`, the first row has minimum
distance `9/121>1/14`; the second has a runner at zero.  Moreover every
denominator through 121 divides its last speed.

This proves the precise loss:

```text
radical support preserves: prime occurrence;
radical support destroys:  valuation depth and residue phase;
minimal repair sidecar:     slotwise valuations plus modulus/time.          (19)
```

It does not compare the two rows' global maxima, and it does not exclude an
LRC row.

## 9. Generated task frontier

The strongest next tasks produced by this session are:

1. **Safe-follower x terminal-clock product.**  Construct the exact finite
   automaton product of `P_m` follower state with THM-2352's standard-digit
   plateau clock.  Test which plateau abscissae remain realizable after safe
   restriction.  Do not inherit unrestricted-cylinder surjectivity.
2. **Equality-boundary sidecar minimization.**  Determine the minimal unbounded
   suffix-match state needed to distinguish strict `S` from closed `K` and
   compare it with the ordinary stabilization clock.
3. **Mixed-power deletion laboratory.**  Delete one speed from (7), compute
   the exact change in maximizing-time multiplicity, and ask which relation
   rows remain sufficient.  This is an LRC carrier stress test, not a route to
   the covering residual unless phase owners are attached.
4. **General `p/q` carry-language renewal.**  Identify the quasi-greedy
   boundary and renewal law for digits `{0,...,q-1}`.  Compare its entropy with
   `log(p/q)` and keep it separate from the `p`-ary torus atom tree.
5. **ABC quality along safe families.**  Compute normalized ABC qualities of
   `(D_m,1,2^m/19)` and `(2^m,19Q_m,3^m)` and locate exact valuation spikes.
   The question is distribution of radical complexity, not horizon finiteness.
6. **Valuation-depth LRC bank.**  Generalize the modulus-121 radical twins to
   the live denominator clocks and determine the least valuation truncation
   that preserves each bounded witness bank.
7. **Ordinal with evaluator firewall.**  Use shortlex ranks for tasks, but
   store `(level, path, evaluator, operation word, terminal sidecar)`.  Test
   round-trip decoding before treating the ordinal as mathematics.

The session's main methodological gain is a reusable diagnostic:

> When the same finite address supports two exact evaluators, solve both on a
> hostile family and name the first coordinate on which their predicates
> disagree.  Here it is phase orientation for LRC versus Mahler, and valuation
> depth for radicals versus LRC.
