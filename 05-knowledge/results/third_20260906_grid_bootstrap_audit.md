# Independent audit of the finite translated-grid bootstrap

**ACCEPTED: PROVED analytic reductions with a FINITE-EXACT complete clock
certificate.** The audit accepts the baseline
[grid bootstrap](third_20260906_grid_bootstrap.md): its complete declared
relaxation leaves exactly **8,301 scales**, with maximum **23,760**, for a
hypothetical balanced actual decoder failure. Every full survivor array
was independently reconstructed. A surviving scale need not be realizable
or unsafe. LRC(14) remains OPEN.

The independent path uses literal interval intersections for every atlas
pair, bounded knapsack for both marginal cost problems, and a separate
integer C++ enumeration of every clock. It imports no repository
producer. The later, separate full-profile cost refinement is not part
of this baseline audit.

## 1. Inputs and actual-entry quantifiers

The audited [translated-grid theorem](third_20260906_grid.md) supplies
both the weak-safety implication and the initial clock range
`1<=t<=97096`. Its proper six-component LRC input remains CITED. The
[audited selected-edge theorem](third_20260906_decoder_profile_graph.md)
supplies some actual atlas edge with sheet gcd `e<=6` for a full-profile
surviving balanced entry. The ratio of that edge remains unrestricted
in the complete strict 5,855-pair atlas. The audit does not replace
this existence assertion by a claim about every edge.

The physical row is `tV union gU`, where `|V|=6`, `|U|=7`, and
`gcd(V)=gcd(U)=gcd(t,g)=1`. For `d_i=gcd(t,u_i)`, every chosen subset
of U, together with all tV, has gcd equal to the gcd of the chosen
`d_i`. This gives the inherited subset ceilings and the exact allowed
42-value singleton projection. Each `d_i` also divides t. These are
necessary conditions; the compiler deliberately retains a larger domain
than the full physical profile constraints.

The audit reads the inherited JSON through
[its frozen computation copy](../../04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json),
whose bytes have SHA256

```text
935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f.
```

This matches the producer's separately located profile supplier. The
42 scalar values are reconstructed from the actual level-six profile
pairs and checked against the stored projection.

## 2. Independent full component geometry and all-real envelope

For each coprime atlas pair `p<q`, use the common integer denominator
`D=14pq`. The intervals of `D_p` have endpoints `(14k±1)q`, and those
of `D_q` have endpoints `(14l±1)p`, in these units. The audit splits
wraparound intervals at zero, intersects the two sorted interval lists,
discards zero-length contacts, and joins the two boundary pieces of
the origin component. This directly recovers every component length
as a positive integer numerator over D.

The atlas itself is reconstructed by a prime sieve and multiplicative
generation of sums: primes are congruent to two modulo three, with
exponents at most two, and the sum is at most 356. The sum two gives
no distinct positive pair and is excluded from the actual ratio atlas.
The result is exactly 94 sums and 5,855 coprime ratios, in the strict
scope of **THM-3818**,
[scaled inert cubeclass pair packet](../../01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md).

For all 5,855 pairs, the literal component multiset agrees exactly
with the claimed origin numerator `2p` once and the clipped numerators
`min(2p,p+q-14r)` twice. Thus containment, wrapping, multiplicity, and
strict endpoint conventions are independently verified for the full
universe, not only at the extremizing pairs.

Writing mu for total length and J for the number of components, the
independent literal geometry also verifies

```text
J<=51, mu>=1/70, 304500*mu>=4313+37J
```

for every pair. Both extremal lines occur, at `(1,10)` and `(5,348)`.
The claimed two-line envelope follows on every real `x>=0`: compare
at zero and `x=304500/37`, interpolate the affine inequalities on that
interval, and use the lower slope bound thereafter. No sampled-x
extrapolation enters this argument.

Integer rounding of the finite minimum commutes with minimum because
ceiling is monotone, and subtracting the integer J commutes with
ceiling. This proves the aggregate formula. The independent compiler
does not need that two-line quotient: it retains the exact minimum mu
separately at every J and minimizes those rounded lines over all
eligible e. This supplies a different exact prefilter.

For the full component credit, each positive-length open component
contributes at least `ceil((t/e)*ell)-1` quotient-grid points, then
lifting multiplies its count by e. These contributions are individually
nonnegative. The draft's sentence allowing negative component credit
was corrected before acceptance: a negative lower bound is possible
for the merged aggregate expression, but not for the sum over positive
individual components. This repair changes no formula or computation.

## 3. Marginal capacities and forced pair reservation

Let `c(d)` be the number of entries in `(90,30,9,4,2,1,1)` that are at
least d. If m of the seven sheet states equal d, their m-subset gcd is
d, forcing `m<=c(d)`. Thus these are necessary individual multiplicity
bounds. They do not enforce all joint gcd or complement-word constraints.

For each retained divisor d of t, the exact weight is

```text
w_t(d)=7*d*ceil(t/(7d))-t=d*((-(t/d)) mod7).
```

The seven true weights sum to `7E`. The maximum over the declared
capacity bag therefore gives a valid upper bound E_D. Every seven-term
weight sum is divisible by seven, so division introduces no rounding.

The audit computes this maximum by bounded knapsack on an exact
seven-slot budget. For each value d, the dynamic program considers
all uses from zero through its capacity, rather than selecting and
sorting copies as in the producer.

For an actual selected pair `u=hp`, `v=hq`, and `e=gcd(t,h)`, write
`t=e*n` and `h=e*k`. Then `gcd(n,k)=1`, so

```text
gcd(t,hp)=e*gcd(n,p),   gcd(t,hq)=e*gcd(n,q).
```

These are the exact two forced marginal states, independently of the
unobserved pair scale h. Their gcd is e because `gcd(p,q)=1`.

The audit rejects absent forced values or capacity overflow. Otherwise
it removes the two reserved uses and solves a fresh bounded knapsack
with exactly five slots. Adding the forced weights gives the exact
maximum for the declared conditioned capacity relaxation. Hence
`E<=E_pair<=E_D`. The two copies must both be reserved when the two
forced states coincide.

Weak safety follows when the uniformly valid component credit exceeds
this upper cost. Consequently a clock survives the relaxation exactly
when there exists some eligible divisor e and some atlas pair with
valid margins and `credit<=E_pair`. Equality is retained. The
existential survivor predicate is the correct necessary condition
for an actual selected edge to remain possible.

## 4. Complete independent clock census

A standalone C++ engine performs exact integer arithmetic at **every
`t=1,...,97096`** in both declared domains. It checks every divisor
`e<=6` in the profile mode and every divisor `e<=30` in the scalar
mode. It never replaces this set by its largest element. At every
retained e, all atlas pairs remain available. A pair is skipped only
if its merged aggregate bound already exceeds the unconditioned E_D;
its component credit cannot then survive either stage.

At a surviving clock, a witness permits an early exit from the
existential scan. At a rejected clock, every eligible e and pair is
excluded by an explicit valid bound or a literal margin check. Thus
the enumeration proves both membership and absence in the declared
finite relaxation.

All six independently generated complete arrays agree exactly with
the producer's frozen arrays:

| Domain | Aggregate count / max | Component count / max | Conditioned count / max |
|---|---:|---:|---:|
| Scalar d<=90, e<=30 | 34,532 / 88,920 | 32,294 / 74,550 | 32,272 / 74,520 |
| Exact 42-value projection, e<=6 | 9,498 / 27,360 | 8,308 / 23,760 | **8,301 / 23,760** |

The seven scales removed only after pair reservation agree exactly:
`12425,14872,14910,15390,15504,20520,21240`.
The complete profile-conditioned array has compact-JSON SHA256

```text
a25d83f0eeb630bb82e84cdfac4e3cf7312f892f6c426d6affd5239a064e4b58.
```

Exact controls retain the mechanisms most vulnerable to lost data:

* `t=57995`, e equal to 1, 5, 7 gives aggregate credits 828, 825, 826:
  numerical divisor order does not control rounded credit.
* `t=74550,e=30,(p,q)=(1,355)` has zero full component credit; open
  components one quotient-grid spacing long can contain no grid point.
* `t=27360,e=6,(p,q)=(5,348)` has component credit 294, exceeding its
  merged credit 252.
* `t=23760,e=6,(p,q)=(25,294)` has component credit 246, below the
  corresponding credits from both measure-envelope profiles. The whole
  atlas must be restored when individual lengths are consumed.

These are exact method controls, not unsafe physical examples.

## 5. Reproduction and outcome

The [audit source](../../04-computation/third_20260906_grid_bootstrap_audit.py)
contains the independent C++ engine and requires a C++17 compiler. It
builds that engine in a temporary directory. The
[frozen audit output](third_20260906_grid_bootstrap_audit.out) retains
all six complete survivor arrays and their hashes. Checks remain
active under Python optimization and C++ optimization.

```bash
python3 -B 04-computation/third_20260906_grid_bootstrap_audit.py
python3 -B -O 04-computation/third_20260906_grid_bootstrap_audit.py
```

Normal and optimized outputs match byte for byte. The run passes
25,018,929 C++ checks and 17,595 Python checks. Raw LF SHA256 values:

```text
source ec841b2a3e9898dcf675d4c0212b78e3d16e21b8d13a05aa54b535a231226cc0
output 36a83b08fecc1187e02c916dc255c77ff4448bc9167b247513b6cdc7770fbc16
semantic 580cc630f4f7eac605aafa727bc199262fb9dda106eccc723f3e84403ceb156f
```

The combination of the inherited clock bound, the actual some-edge
supplier, the uniform translated-grid count, and the complete exact
census proves the claimed necessary 8,301-scale restriction. The
unbounded component shapes and actual joint realization questions
remain separate. No surviving-clock statistic is promoted to a
strict failure, and no LRC(14) closure follows.
