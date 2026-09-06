# Two overlap profiles and a complete finite clock sieve

**Status: PROVED analytic reductions; FINITE-EXACT complete compiler;
INDEPENDENTLY AUDITED.** The complete declared sieve leaves **8,301 integer scales**,
all at most **23,760**, for a hypothetical balanced actual decoder failure.
Every surviving scale is recorded explicitly. A survivor is not a failing
row, a realizable decoder entry, or a proof that the sufficient count is
sharp. LRC(14) remains open. No theorem ID or external priority is claimed.

The new quantitative conclusion depends on the independently audited
[small sheet-gcd edge theorem](third_20260906_decoder_profile_graph.md),
which supplies **some** actual edge with `e<=6`. The selected edge's ratio
remains unrestricted in the entire 5,855-pair actual atlas.

## 1. Inheritance, board, and the cheapest hostiles

The [translated-grid theorem](third_20260906_grid.md), with its
[independent audit](third_20260906_grid_audit.md), already gives the
necessary bound `t<=97096` under hypothetical balanced failure. It supplies
exact marginal ceiling costs and pair-overlap credit, retaining strict bad
arcs and weak-safe endpoints. The actual ratio atlas is **THM-3818**,
[scaled inert cube-class support-two packet](../../01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md).
The pair measure is the main identity in **THM-739**,
[pairwise-coprime bad-overlap Bernoulli form](../../01-canon/theorems/THM-739-pairwise-coprime-bad-overlap-exact-bernoulli-closed-form.md),
together with the corrected clipped component geometry in the new grid
theorem. The faulty unclipped microscopic formula is not reused.

The [joint-shadow profile supplier](lrc14_joint_shadow_empty_core_next_sep06.md)
gives the scalar subset caps and exact allowed seven-body gcd values. Its
[JSON](lrc14_joint_shadow_empty_core_next_sep06.json) is pinned by SHA-256:

```text
935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f
```

The lower-dimensional lonely-runner existence input remains **CITED** as
documented in those suppliers. It is not re-proved in this compiler.

The six-concept board is: clock divisors; hereditary gcd capacities; the
measure/component tradeoff of an actual edge; individual arc lengths;
common-scale quotient arithmetic; and the literal surviving clock set.
The operation moves from a uniform scalar ceiling to a finite clock with
its divisors retained. It then restores the actual pair's two marginal
gcds, which the independent cost and edge bounds had forgotten.

Three small hostiles determine the necessary sidecars:

* At `t=57995`, the admissible edge gcds `e=1,5,7` give aggregate overlap
  lower credits `828,825,826`. The worst rounded credit need not occur at
  the numerically largest divisor. Every admissible e must be retained.
* At normalized clock `n=2485`, the pair `(1,355)` has51 open components
  each of length `1/2485`. Each can have zero points by the open-interval
  lower count. Erasing the endpoint subtraction is invalid.
* At `t=23760,e=6`, pair `(25,294)` has individual-component credit246,
  below the corresponding credits for both profiles that minimize the
  aggregate measure envelope. A quotient preserving total measure and
  component count need not preserve the consumer using separate lengths.

Targeted inheritance searches found the prior uniform bound and proposed
divisor refinement, but not this complete clock compilation, two-profile
envelope, or forced-pair marginal consumer. No literature priority follows.

Incoming [actual packet work](continuing4_20260906_lrc_packets.md) independently
uses separate interval ceilings on a fully filtered family. That overlapping
mechanism is acknowledged; the complete clock compilation, full-profile edge
supplier and forced-pair marginal bounds here have their own declared scopes.

## 2. The whole aggregate atlas has only two extremal profiles

For an actual coprime pair `1<=p<q`, put

```text
R=ceil((p+q)/14)-1,        J=2R+1,
mu=[p+sum_(r=1)^R min(2p,p+q-14r)]/(7pq).
```

The declared atlas contains exactly5,855 pairs: `p+q<=356`, every prime
of the sum congruent to2 modulo3, and every exponent at most two. For
every real `x>=0`, its exact lower envelope is

```text
min_atlas (x*mu-J)
 = min(x/70-1, 62*x/3045-51).                      (1)
```

The two displayed profiles are `(p,q)=(1,10)` and `(5,348)`. They cross
at `x=304500/37`, with common value `4313/37`.

This has a compact exact finite certificate. Every pair satisfies

```text
J<=51,              mu>=1/70,
304500*mu >= 4313+37J.                             (2)
```

At x=0 and at the crossing, (2) puts every affine function above the
first active segment of the proposed envelope. Linearity gives the whole
interval between them. Beyond the crossing, the lower slope bound in
(2) proves domination of `x/70-1`. Both bounding lines occur in the atlas,
so equality holds. All three inequalities are verified exactly for every
pair; no numerical hull approximation is used.

For any integer `e|t`, the strongest aggregate-component ceiling uniformly
over all actual pairs is consequently

```text
A(t,e)=e*min(ceil(t/(70e))-1,
             ceil(62t/(3045e))-51).                (3)
```

Indeed `ceil(x*mu)-J=ceil(x*mu-J)` because J is integral, and the ceiling
of the finite minimum equals the minimum of the ceilings. Thus this
specific rounding operation preserves the two-profile envelope. Section4
explains the different operation for which it does not suffice.

## 3. At a fixed clock, only actual divisors pay ceiling costs

Under hypothetical failure, write

```text
tV union gU, |V|=6, |U|=7,
gcd(V)=gcd(U)=gcd(t,g)=1.
```

Let `d_i=gcd(t,u_i)`. Every k-subset of these seven values has gcd at
most the kth entry of

```text
(90,30,9,4,2,1,1).                                 (4)
```

For an allowed value d, define its relaxed multiplicity capacity by

```text
c(d)=number of entries in (4) at least d,
w_t(d)=d*((-(t/d)) mod7),       d|t.
```

The parent theorem's excess `E=B-t` satisfies `7E=sum_i w_t(d_i)`.
For the scalar-only relaxation, retain all divisors `d|t` with `d<=90`;
for the full-profile projection, also require d to lie in the exact
42-value seven-body set from the pinned JSON. Construct a finite bag
containing c(d) copies of `w_t(d)`. Let `E_D(t)` be one seventh of the
sum of its seven largest entries. Then

```text
E <= E_D(t).                                      (5)
```

The bag always has at least seven entries because d=1 has capacity7.
Each retained weight is congruent to `-t` modulo7, so the sum of any
seven weights is divisible by7; no unrecorded rounding is needed in (5).
Selecting the largest values is the exact maximum of this declared
per-value multiplicity relaxation. It does not enforce every joint
subset or complement-word condition.

For the scalar case an actual selected edge has `e|t,e<=30`; the
small-edge supplier gives some actual edge with `e|t,e<=6` in the full
connected decoder domain. Equation (3) supplies weak safety whenever

```text
min_(e|t,1<=e<=cap) A(t,e) > E_D(t).                (6)
```

The strict inequality is essential. The largest aggregate survivor in
the full-profile/e6 case is `t=27360`, where `E_D(t)=252=A(t,6)`.
Replacing `>` by `>=` would falsely delete that retained boundary.

## 4. Keep each component's length before applying its ceiling

The individual component lengths have common denominator `D=14pq`:

```text
2p/D                           once,
min(2p,p+q-14r)/D               twice for each r=1,...,R.
```

For `n=t/e`, let `C(t,e,p,q)` be

```text
e*[ceil(n*2p/D)
   +2*sum_(r=1)^R ceil(n*min(2p,p+q-14r)/D) - J].    (7)
```

Each open component of length ell contains at least `ceil(n*ell)-1`
points of any translated n-grid. Therefore (7) is a valid lower overlap
credit, and is at least the aggregate credit for that same pair. It is
nonnegative because every component is nonempty and every individual
ceiling is at least one. The merged aggregate lower bound may be negative.

The sum of the individual ceilings can exceed the ceiling of the total.
For example `t=27360,e=6,(p,q)=(5,348)` has aggregate credit252 and
individual-component credit294. The compiler therefore repeats the actual
pair scan after the cheap aggregate filter. It retains **every** pair:
the two-profile envelope (1) is not passed to this different consumer.

A clock survives the component relaxation precisely if some admissible
e and some actual pair satisfy

```text
C(t,e,p,q) <= E_D(t).                              (8)
```

This condition describes the compiler's relaxation exactly. It does not
say those component minima can be attained simultaneously on one grid,
or that the other five bad masks cover the remaining points.

## 5. The pair ratio forces two of the marginal sheet gcds

Write the selected physical U-pair as `u=hp,v=hq`, with `gcd(p,q)=1`,
and put `e=gcd(t,h)`. Since `gcd(t/e,h/e)=1`,

```text
d_u=e*gcd(t/e,p),       d_v=e*gcd(t/e,q).            (9)
```

In particular `gcd(d_u,d_v)=e`. This is an exact identity, independent
of the unobserved common pair scale h. It connects the selected pair's
primitive coefficients to the marginal cost bag.

For every `(t,e,p,q)`, reject it if either forced value is absent from
the bag or if reserving both exceeds a value's multiplicity capacity.
Otherwise reserve these two entries and take the five largest remaining
costs. Their seven-term total divided by7 is `E_D(t;p,q,e)`. Then

```text
E <= E_D(t;p,q,e) <= E_D(t).                        (10)
```

The final compiler retains exactly the clocks for which some admissible
e and pair have valid forced margins and

```text
C(t,e,p,q) <= E_D(t;p,q,e).                         (11)
```

This recovered coordinate excludes the endpoint witness at `t=74550`,
`e=30,(p,q)=(1,355)`: although its separate component count is zero,
equation (9) forces `d_v=10650>90`. Other pairs must still be checked;
excluding one pair alone cannot certify a whole clock.

In the full-profile/e6 mode, this coupling removes precisely seven
additional scales after (8):

```text
12425,14872,14910,15390,15504,20520,21240.
```

The improvement is modest and exact. The source retains the complete
all-pair scan proving each exclusion, rather than inferring it from one
illustrative forced-margin failure.

## 6. Complete finite result and physical consequence

The declared universe is **every integer `1<=t<=97096`**, supplied by the
previous proved scale bound. For every t the compiler checks every allowed
clock divisor e and every needed actual ratio. A ratio is skipped only
when its already proved aggregate bound is strictly too large to survive.

| Retained marginal domain and edge bound | Aggregate ceilings: count / max | Individual components: count / max | Forced-pair components: count / max |
|---|---:|---:|---:|
| All d<=90; e<=30 | 34,532 / 88,920 | 32,294 / 74,550 | **32,272 / 74,520** |
| Exact42-value projection; some e<=6 | 9,498 / 27,360 | 8,308 / 23,760 | **8,301 / 23,760** |

The last maximum retains the relaxed witness

```text
t=23760, e=6, (p,q)=(25,294),
(d_u,d_v)=(30,36), E_D=E_D(t;p,q,e)=252,
C(t,e,p,q)=246.
```

This certifies membership in the declared relaxation, not an actual
thirteen-speed counterexample. By contrast, every omitted clock is
excluded by a sufficient weak-safety count valid uniformly over the
selected edge supplied by the inherited theorem.

To transfer the finite result, argue by contradiction from a hypothetical
balanced actual failure. All inherited profiles hold, the small-edge
theorem supplies some edge with e<=6, and that edge must satisfy (9).
If its t is absent from the final8,301-element set, every candidate for
that edge has overlap credit strictly exceeding the total ceiling excess.
A lift of the cited six-body safe phase is then weak-safe on all seven
remaining labels, contradicting failure. Thus every such failure has

```text
t in S,    |S|=8301,    max S=23760,    1<=g<=90.    (12)
```

The scalar row of the table needs only the earlier caps and one actual
atlas pair. The stronger row additionally consumes full inherited profiles
and connected actual U geometry through the small-edge theorem. Their
hypotheses are not interchanged.

## 7. Exact output and reproduction

The standalone [source](../../04-computation/third_20260906_grid_bootstrap.py)
and [full retained output](third_20260906_grid_bootstrap.out) import no
producer code. Each `SCALE_SET` line contains a JSON object with `name`,
`count`, `maximum`, `sha256`, and the full sorted `survivors` array. The
hash is SHA-256 of the compact JSON array encoded in UTF-8. Six complete
sets, including the scalar baseline and both intermediate relaxations,
are retained; no selected-height sampling or probabilistic test occurs.

The verifier reconstructs all5,855 atlas pairs, checks the three exact
envelope inequalities, compares every total length with the independent
Bernoulli formula, and checks selected complete length multisets against
literal rational phase-wall sweeps. It retains the largest-divisor
rounding hostile, strict aggregate equality, interval-containment geometry,
and a consumer for which the two-profile quotient loses needed data.
Literal t/h controls independently verify both forced marginal identities.

```text
python3 -B 04-computation/third_20260906_grid_bootstrap.py
python3 -B -O 04-computation/third_20260906_grid_bootstrap.py
```

Normal and optimized output agree byte for byte after **5,094,668 active
requirements**, including105 literal length-multiset sweeps. Raw LF-byte
SHA-256 values are:

```text
source 1c4070e17aaf1825d07899f8b2e056d0e2f0b05f224e90997f732e079814eb3e
output 1a9bc517544bfae2e6061bdf49e9b93dbd6f0722ccd99c9ac8a1b4655ca4838d
```

The final `profile6_coupled` survivor array has compact-JSON SHA-256
`a25d83f0eeb630bb82e84cdfac4e3cf7312f892f6c426d6affd5239a064e4b58`.
The [independent audit](third_20260906_grid_bootstrap_audit.md) accepts the
full proof and reconstructs every retained array through literal interval
intersections and an independent bounded-knapsack C++ engine. It imports
no producer code. The separate [full-word refined consumer](third_20260906_grid_refined.md)
restores an unconditional full-word maximum and reduces the necessary set
further. Joint full complement-word feasibility conditioned on every
selected ratio, together with all remaining actual ratios, is still not
imposed by (11). Its surviving scales are not an exact classification of
genuine failures.
