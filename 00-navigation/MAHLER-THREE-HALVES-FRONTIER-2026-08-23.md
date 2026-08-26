# Mahler `3/2` frontier

**Current status (2026-08-25): OPEN.**  This is the routed detail behind the
compact entry in [CURRENT-FRONTIER](CURRENT-FRONTIER.md).  Primary-source
records live in
[CORE-PAPERS-MAHLER-THREE-HALVES](../05-knowledge/reference/CORE-PAPERS-MAHLER-THREE-HALVES.md).

## Closest inherited mechanisms

- **PROVED:** THM-2228 separates strict real carry tails from ordinary
  positive-integer stabilization.  The periodic carry control `(01)^infinity`
  decodes to `Phi=-2/5` and never stabilizes.
- **PROVED:** THM-2352 realizes every finite carry cylinder in its unrestricted
  abscissa universe.  It does not preserve the Haar-null safe follower tree.
- **CITED:** Mahler gives countability and a height-counting bound;
  Dubickas--Mossinghoff force any candidate above `2^57`.
- **CONJECTURAL:** the 2026 rational-base normality preprint would exclude
  `Z_(p/q)` numbers for `1<q<p<q^2`; the normality premise is not proved.

## Exact finite prefix structure

[THM-3848](../01-canon/theorems/THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation.md)
proves, for every coprime `p>q>=2`, that the normalized length-`N` safe-prefix
set has a common wall grid, exactly `p^N` safe open atoms, Haar measure
`q^(-(N+1))`, and a reversible full `p`-ary atom refinement.  Pointwise there
is instead exactly one safe lift among the `q` inverse sheets.  These are two
different maps; identifying their branching factors destroys the wall address.

The adjacent relations `p e_n-q e_(n+1)` generate the full integer kernel.
The mixed-power speed row has exact maximum

```text
M(q^N,q^(N-1)p,...,p^N)=floor((p+q)/2)/(p+q),  N>=1.
```

For `p/q=3/2`, the maximum is `2/5` and the only maximizing times are the two
mod-five phases.  At `N=12` this gives a primitive thirteen-speed, hence
fourteen-runner, positive control with twelve norm-five relations.  It misses
every multiple of five and is sieve-trivial, so abundant short relations and
a perfect recursive address do not characterize the LRC(14) hard core.

## Where LRC and Mahler separate

The LRC evaluator retains centered distance `||x||` and identifies the phases
`2/5` and `3/5`.  Mahler's predicate is oriented: it accepts only the lower
half-circle.  Along the mod-five extremizer the phases alternate, so the exact
LRC maximizer fails the Mahler test.  This phase side is the first destroyed
coordinate in the map from an oriented orbit to centered loneliness.

The distinct closed formal safe-tail shift `K` has greedy equality boundary
`d`, renewal law

```text
a_m=1+sum_(k=1)^m d_k a_(m-k),
A(z)=1/((1-z)(1-D(z))),
```

entropy `log(3/2)`, binary-ultrametric Hausdorff dimension `log_2(3/2)`, and
is nonsofic.  The strict set removes the countable backward orbit of `d` and
has the same finite language.  This does **not** prove that the unknown
Z-language is nonsofic: a Mahler candidate must also satisfy THM-2228's
ordinary positive stabilization condition.

## ABC/IUT sidecar

Assuming ABC, actual odd carry-one packets have asymptotic radical exponent at
least one.  The denominator-19 safe-prefix family at positive multiples of 18 is
also logarithmically radical-saturated.  Neither fact supplies a horizon bound:
arbitrarily long even zero-carry runs have adjacent radical `6`, and the
denominator-19 construction changes its initial integer with `m`.  IUT adds no
unconditional input here; only the contested claimed IUT-to-ABC implication
could feed the explicitly conditional ABC consumer.

## Exact fibre product and finite-state obstruction

[THM-4072](../01-canon/theorems/THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md)
closes the former finite task.  The safe follower graph has countable state
`q`: nonrejection is the closed language `K`, while strict safety is exactly
infinitely many reset edges.  Carry and native digits are related by a
rooted-binary-tree automorphism with the online integer sidecar

```text
e=c xor (u mod 2),
u'=(3u+c+3^(m+1)e)/2.
```

Thus a candidate code is exactly a synchronized path that never rejects,
resets infinitely often, and has finitely many but at least one native `1`.
Every finite terminal-prefix test removes **zero** safe nodes: its inverse
limit closes strict safety to `K` and eventual-zero native paths to the full
binary path space.  The rooted pair-prefix language is nonregular, and no
finite synchronous transducer changes carry to native coordinates.  The
greedy boundary, `(100)^infinity`, `A=1`, `(01)^infinity`, and the
denominator-19 tower are exact controls, not candidate evidence.

The same theorem now records the exact clocked-cylinder and reset-skeleton
normal form. If a carry prefix of length `m` sends its least representative
`r_m` to `u_m`, then every nonnegative integer in that cylinder obeys

```text
A=r_m+2^m k  =>  T^m(A)=u_m+3^m k.
```

Because `3^ell` is a unit modulo every `2^j`, every finite prescribed reset
skeleton cuts out exactly one residue class modulo `2^L`. This is a finite
compatibility theorem, not an infinite candidate. The minimal hostile for the
terminal-depth-retaining quotient is `A=8` versus `A=13`: after four steps
they share depth, follower state, height, and the native-one-seen flag, but
the next carry `1` rejects the former while `0` keeps the latter alive. Thus
an all-depth argument must supplement this quotient with output information,
such as the residual integer/cylinder address.

## Live deterministic task

After the last native `1`, the product becomes the deterministic integer
system

```text
(q,u) -> (delta(q,u mod 2), ceil(3u/2)),
```

with rejection when `delta` is undefined.  Starting from every reachable
safe state with `u>0`, prove that the orbit rejects or loses the infinite-reset
condition; alternatively exhibit one reachable nonrejecting orbit with
infinitely many resets. Finite reset-skeleton solvability alone cannot decide
this inverse limit: couple the exact residue classes to an Archimedean or
genuinely all-depth `2`-adic obstruction.

[THM-4074](../01-canon/theorems/THM-4074-mahler-denominator19-postterminal-arbitrary-delay.md)
rules out a bounded version of this task.  A reachable denominator-19 terminal
family has arbitrarily long reset runways and then programs every finite carry
word beginning in `1`; safe words stay nonrejecting and greedy prefixes reach
unbounded follower state.

[THM-4077](../01-canon/theorems/THM-4077-mahler-denominator19-2adic-tangent-full-shift-isometry.md)
closes these programmers to a moving-time odd-`Z_2` tangent isometry, not one
orbit. Parameter/output termination give two different open safe intersections.

## Renormalized cross-scale chart

[THM-4082](../01-canon/theorems/THM-4082-mahler-renormalized-linear-chart-and-exact-bit-defect.md)
rescales the depth-`s` parameter by `3^s` and puts every THM-4077 isometry in
one chart. With `Lambda=(243/152)log(3^18)`, an odd 2-adic unit, the exact law

```text
v_2(H_s(x)-Lambda*x)=s+2+2v_2(x),       x!=0,
```

says the nonlinear and linear carry words share exactly that many initial
bits and differ at the next. The induced near-identity isometry transports
the strict-safe and output-termination fibres exactly. Parameter termination
retains the separate sidecar `P_s=3^s N_odd`: these dense countable loci are
strictly decreasing with empty intersection, while output-only and neither
termination remain dense at every fixed scale. Both strict-safe intersections
and the deterministic orbit task remain open; the chart produces no Z-number.
