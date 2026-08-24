# LRC(14): fixed-radius auxiliary centers, t-sheet response, and the cubic tournament sidecar

**Session status -- 2026-08-23.**  THM-3910 is proved relative to cited
`LRCUpTo13`, verified exactly, and independently hostile-audited.  In the
conditional `t>=U` part of THM-3878's `11+2` branch it closes 41 of 57
scale-one certificate survivors, reducing the 58-type ledger to 17.  It also
proves a body-component filter, an exact integer lift-count/Fourier carrier,
a sharp integer failure tariff, and a native no-go for pairwise tournament
compression.  The `t<U` slice and LRC(14) remain open.

## Inheritance pass

- **Closest proved mechanism:** THM-3878's cyclic-slack and one-auxiliary
  deletion certificate.  It already identified open pair-obstruction
  components as the correct target and the final 58 conditional types as the
  exact hostile bank.
- **Canonical hostile:** the scale-two `(1,9)` quotient and the 57 scale-one
  rows on which scalar safe mass, inversion symmetry, component count, and
  isolated AP11 walls close nothing.
- **Corrected near miss:** MISTAKE-464 distinguishes AP11's 14
  positive-length arcs from its 18 closed components.  The four isolated
  walls are real but have zero Haar mass and zero `lambda_+` response.
- **Least-used sidecar:** THM-3377's first deletion response and THM-3729's
  covariantly transported root.  The analogous LRC object is not a tournament
  score but the rooted sheet multiplicity `N_t(w)`.
- **Hostile operation lesson:** THM-3380 warns that a first skew current may
  vanish.  Here inversion symmetry makes a skew orientation unnatural, while
  the actual consumer is a symmetric cubic mixed response.

## Anchor / Niche / Wildcard portfolio

### Anchor

Close THM-3878's 58 `t>=U` certificate survivors without confusing a
sufficient certificate with a counterexample classification.

### Niche

Treat the longest closed body-safe arc as an operation-compatible response,
retaining compact/open endpoint semantics instead of reducing the body to
scalar measure or a component count.

### Wildcard

Recover the tournament ancestry by asking exactly which incidence order the
target consumes.  This produced a three-vertex response hypergraph and a
native AP11 collision showing that complete pairwise Gram/cut data lose the
load-bearing cubic moment.

## Live concept board

1. `A_pq` or `C_pq`: the literal open pair obstruction.
2. Auxiliary-deep center: `||a w_0||>=1/13`.
3. `lambda_+(u)`: longest positive closed body-safe component.
4. `N_t(w)`: integer number of body-safe sheets above the rooted `t`-cover.
5. `F(z,w)`: rooted multiaffine deletion response with cubic term `M_pq`.
6. Endpoint owner: an atomic sidecar invisible to Haar moments.

The meaningful connections were not analogies between names.  Each new move
was checked against all six objects: what it preserves, what it destroys,
and which literal target it evaluates.

## The response reinterpretation that changed the radius

THM-3878 applies cited lower-runner LRC to

```text
u union {a t}.
```

Its earlier deletion certificate kept the auxiliary speed safe over the
whole perturbation.  In `w=t y` coordinates that gives radius `1/(182a)`.
The first center-depth scout retained that choice and, even after adding the
strong condition `||a w_0||>=1/13`, closed none of the 58 survivors.

Tournament deletion history suggested a different question: which predicate
must survive the operation?  The auxiliary runner is needed only to locate a
marked center.  The target perturbation needs to preserve the eleven-body,
not the auxiliary.  Since `t>=U`, the body alone remains safe on the fixed
closed radius

```text
rho=1/13-1/14=1/182.
```

Thus a counterexample requires an `a`-deep center in the **strict**
radius-`rho` erosion of the open pair obstruction.  This stronger response
closes 41 scale-one rows, always with the intrinsic multiplier `a=p`.

This is the central reusable move:

```text
use an added object to select a center;
delete it before transporting the neighborhood;
preserve only the predicate consumed by the target.
```

The move is neither monotonicity nor ordinary deletion.  It is a marked
response calculation, and the choice of retained predicate changes the
metric scale.

## Why the finite exhaustion is genuinely finite

For a pair component of length `beta`, fixed-radius erosion leaves width

```text
gamma=beta-1/91.
```

An `a`-shallow gap `{||a w||<1/13}` has component length `2/(13a)`.  Once

```text
a gamma>2/13,
```

the largest erosion core cannot avoid the `a`-deep set.  Only 398 exact
multiplier cells remain, with terminal cutoffs 2 through 12.  Two standalone
open-component implementations agree on every cell and on the 41/16 split.

Strictness has visible witnesses.  At `(4,13,a=4)` and `(7,13,a=7)`, the deep
set touches only the erosion boundary.  Such a point would make the closed
source arc touch the complement of the open obstruction, so it does not
survive.  This is the same compact/open discipline that repaired THM-1042.

## The body-shape response

If `lambda_+(u)` is the longest positive-length closed component of the body
safe set, a counterexample must pay

```text
t lambda_+(u)<beta_1(p,q)
```

or its scale-two analogue.  Equality closes: a compact arc of the same length
as an open component cannot be contained in it.

AP11 has the owner-labelled arc

```text
[1/14,13/154], length 1/77,
```

so `U lambda_+=1/7`.  This closes all 57 scale-one AP11 controls, including
the four equality-width pairs, and the scale-two control with `beta_2=2/63`.
The response staircase

```text
1/35 -> 15,  1/28 -> 34,  1/21 -> 46,
1/14 -> 52,  1/7  -> 57
```

turns a geometric sidecar into a prioritized target: any structural body
theorem reaching one of these normalized lengths closes a declared portion
of the bank.

## The integer sheet carrier

For `G=G_u`, define

```text
N_t(w)=#{a mod t:(w+a)/t in G}.
```

Then the scale-one safe mass is exactly

```text
(1/t) integral_H N_t,
```

where `H` is the pair-safe set.  Its Fourier transform retains only the
`t`-divisible modes:

```text
Nhat_t(k)=t 1hat_G(k t).
```

This is the correct natural-number address: the root is the integer `t`
itself, not a cosmetic square, rank label, or averaged orbit.  The sheet
multiplicity is integer-valued, so failure pays more than Cauchy--Schwarz
unless its mean over the obstruction is already integral.  With

```text
m=t mu(G), k=floor(m/b), theta=m/b-k,
```

the exact extra payment is `b theta(1-theta)`.

The BV estimate gives `Var(N_t)<=r^2/3`, hence an all-large-`t` theorem for
each fixed body.  It is deliberately nonuniform: no current result controls
`r/mu(G)` uniformly over the admissible eleven-components.

AP11 exposes the analytic frontier cleanly.  The exact CV gate fails only at
the resonant scales `t=13,26`, closes every pair from `t=27`, while the coarse
BV estimate becomes uniform only at `t=88`.  The next useful analytic object
is therefore the actual `t Z` harmonic energy, not another total-variation
bound.

## Tournament connection contract

```text
source:
  the labelled incidence system (G,D_tp,D_tq), or the t-sheet fibre word;

target:
  positive-measure simultaneous safety;

map:
  push Haar measure to the Boolean response deck, or push G along y -> t y;

preserved:
  exact safe mass, the rooted t-harmonic energy, and integer sheet count;

destroyed by pairwise quotient:
  the cubic mixed moment, individual sheet labels, cyclic cell order,
  endpoint owners, and the fixed root if t is averaged;

needed sidecar:
  M_pq for positive Haar mass; endpoint owner on the zero-mass boundary;

cheapest decisive test:
  AP11, (p,q)=(6,41), t=12 versus t=36.
```

Those two rows have identical complete one/two incidence data but different
cubic response and different safe mass by `17/61992`.  Thus an ordinary
tournament or cut semimetric is provably not faithful here.  The faithful
carrier is a rooted response hypergraph, equivalently the integer sheet
carrier.

## Negative results and stopping reasons

- Preserving the auxiliary over the whole interval, even with deep-center
  data, closes `0/58`; its radius is too small.
- Isolated AP11 walls close `0/58` by themselves.  They have zero positive
  component length and cannot feed a Haar-mass gate.
- Endpoint owner switches cannot repair a one-interval deletion certificate
  away from equality.  The sharp survivor `(4,51,4)` has an open obstruction
  interval of length `1/357` strictly containing the guaranteed closed
  `1/364` source.
- Pairwise tournament/Gram data are insufficient even on a native LRC comb;
  16 AP11 pairwise fibres split at the cubic response.
- The sheet-variance theorem is nonuniform in the body and therefore does not
  by itself close the 17 residual conditional types.

Each stopping reason records the missing coordinate instead of merely saying
that a statistic failed.

## Scale-two return: covering number is a physical runner budget

The `(2,1,9)` obstruction has now been pulled in two opposite directions.
First, THM-4002's fixed family extends exactly from eleven-subsets of `[1,15]`
to every eleven-subset of `[1,21]`: 352,716 bodies and 1,356,147 exact finite
`(E,t)` containment instances, followed by the strict variance tail.  An
independent engine measures `G_E\t^{-1}C` directly on the disjoint max-21
stratum, finds positive rational escape in all 781,184 new finite instances,
and rechecks a physical lift.  This is substantial fixed-body scope, but it
does not decrement the 17 arbitrary-body types.

Second, the least-used sidecar—the covering number of the open obstruction
by auxiliary danger combs—has become exact.  No one or two combs cover

```text
C=(2/21,8/63) union (55/63,19/21),
```

while `D_8 union D_9 union D_10` does.  The proof is genuinely all-
multiplier: tooth mass forces the smaller multiplier to be at most 25, the
largest surviving gap leaves 20 necessary pairs, and exact open-endpoint
union rejects each.  `(8,9,10)` uniquely attains minimum possible maximum 10;
triples of larger maximum are not classified.

The connection contract is therefore unusually clean:

```text
source:       the two open quotient arcs C;
target:       a common physical time for an eleven-speed body plus (t,9t);
map:          auxiliary clock a -> danger comb D_a, then the two scale lifts;
preserved:    complete quotient failure and strict endpoint incidence;
destroyed:    physical runner count, body identity, owner and first arrival;
needed sidecar: runner-slot budget and the THM-3818 crossing-row label;
cheapest test: all two-comb covers, with (8,9,10) as hostile positive control.
```

This explains rather than merely reports the stop.  Cited `LRCUpTo13` makes
one auxiliary clock available beyond an eleven-speed body; separate uses need
not share a centre.  Two common clocks are already the general LRC(14)
assertion.  The three-comb control would need fourteen nonzero speeds safe at
`1/14`, stronger in radius than LRC(15).  In THM-3818's rank-eleven scale-two
two-component `W=V_dec` equality branch, `u_i=at` makes actual speeds `2at,t`
a forbidden crossing row of height `2a<=20`; no wider crossing claim is made.
Low relation rank cannot pay a missing physical slot.  The next move must exploit
body-dependent geometry or a simultaneous multi-clock inequality, not reuse
the three-comb cover as if it were a proof.

## Generated next tasks

1. **Remaining-16 erosion atlas.**  Compute the full arrangement of deep
   erosion cores for the 16 rows, not just nonemptiness.  Look for a joint
   incompatibility among several multipliers `a`, retaining the common body
   sheet word as sidecar.
2. **Resonant Fourier passport.**  Explain the AP11 exceptions `t=13,26` from
   the exact `t Z` Fourier lattice and test whether the THM-3818 decoder
   constrains the same low modes for arbitrary bodies.
3. **Normalized body-component lower bound.**  Seek structural hypotheses
   forcing `U lambda_+(u)>=1/35`, then climb the exact response staircase.
   Owner-labelled endpoints and relation support are plausible sidecars.
4. **Integer tariff plus local support.**  Combine the quantized second-moment
   floor with the radius-`1/182` intervals forced for every auxiliary
   multiplier.  The moment problem should remember prescribed local regions,
   not only the measure of the obstruction.
5. **Weighted scale-two carrier after the three-versus-one stop.**  Use the
   exact lift weight `ell(w) in {0,1,2}` rather than discarding it to
   `1_(ell>0)`.  Seek a body-dependent one-clock certificate or a simultaneous
   inequality that does not spend three physical auxiliary runners.
6. **Cubic response kernel.**  Group the exact Fourier expansion of `M_pq` by
   solutions of `i p+j q=r`; seek a sign-controlled low-mode packet before
   the tail can cancel it.
7. **Separate `t<U` program.**  None of the fixed-radius or compact-response
   arguments should be silently transported across `t=U`.  Build a distinct
   physical-entry/owner route for that slice.

The best immediate anchor is task 1, with task 3 as the strongest orthogonal
route and task 5 as the wildcard.  Any claimed connection should continue to
name its source, target, map, preserved predicate, destroyed information, and
cheapest hostile test.
