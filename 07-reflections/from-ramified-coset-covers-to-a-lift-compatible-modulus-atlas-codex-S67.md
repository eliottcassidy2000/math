# From ramified coset covers to a lift-compatible modulus atlas

## Scope

LRC(14) remains open.  This reflection records a structural correction, not a
closure claim.  The common-scale computations, the resonance identities, the
finite-denominator obstruction, and the Route-B audit all point to the same
missing object:

> a signed, exact-period, lift-compatible atlas over moduli, whose local charts
> are residue/sheet masks and whose global sections are actual integer speeds.

The phrase “atlas” is deliberate.  No single quotient seen so far preserves
all of the LRC predicate.  The positive results work by changing between
incomparable exact charts; the failed results promote one chart to a global
closure and thereby forget either signs, integer lifts, or untested
denominators.

## The perspectives accumulated so far

The project has looked at the same object through several mathematically exact
lenses.

| lens | local object | what it preserves | what it forgets |
|---|---|---|---|
| danger-comb geometry | safe teeth and interval components | the actual lonely predicate and component width | arithmetic address across other denominators |
| rational sieve | danger masks on `(Z/qZ)^*` | every witness at one exact denominator | all untested denominators and integer lift compatibility |
| common-scale sheets | masks in `Z/cZ` and quotient fibres | literal sheet incidence and owner deficits | information transverse to the chosen divisor |
| Fourier/resonance | annihilator lattice and signed coefficients | the exact dual of honest quotient periodicity | positivity after a triangle inequality |
| Bonferroni ledger | centered joint moments | interaction order and exact low slots | joint alignment if slots are bounded independently |
| Farey/continued fractions | one maximizing rational chart | local active runners and rational descent | Cover14 and witnesses at other moduli |
| additive energy | differences and approximate progressions | genuine additive structure when present | a reason that band containment must create that structure |
| tournaments | pairwise order/cost shadow | compact fingerprints and switch-sensitive telemetry | shared carriers and higher compatibility |

The sharp lesson is not that one lens is wrong.  It is that each lens has a
precise forgetful map, and the hard families live in the forgotten fibres.

## Honest quotients: the positive model

THM-1090 closes common scale `30` by switching between the honest quotient
maps

```text
Z/30Z -> Z/6Z,          Z/30Z -> Z/10Z.
```

The first flag reduces the scalar survivors to `120` all-owner rows; the
transverse flag closes every remaining owner obligation.  Neither quotient
dominates the other.  This is the first complete example of the correct
architecture: a nonnested atlas of predicate-preserving local views.

THM-1096 adds the complementary power-of-two model.  At common scale `32`,
heredity is simply “at least two top-order coordinates.”  A single honest map

```text
Z/32Z -> Z/8Z
```

retains every lower-order thick fibre and already forces at least four failed
owners in every scalar row.  Here no transverse chart is needed.  The contrast
with scale 30 is informative: the atlas should permit several incomparable
charts, but it should not require them when one exact carrier is terminal.

THM-1091 proves the dual statement.  Complete divisor fibres are equivalent
to annihilator support in the finite Fourier transform.  The sheet-mask
hypergraph and the weighted zero-sum frequency hypergraph are two exact
presentations of the same ramification data.  THM-1092 then proves the raw
full-support resonance series absolutely convergent, so THM-1093's finite
mod-14 coset regrouping is legal.

THM-1093 also gives the warning.  On its twelve sampled triples, almost all
observed cancellation is **between** residue cosets.  Refining the series and
then summing absolute values destroys exactly the interaction the proof would
need.  The correct dual atlas must retain signed gluing data, not merely local
absolute masses.

## Fixed global atlases cannot close the problem

THM-566 and HYP-2876 show that one divisor-loaded coordinate kills any fixed
finite denominator list.  THM-1098 strengthens this in three directions:

1. the obstructing row is primitive and covers `2,...,14`;
2. it is independently proved lonely by an explicit later rational point;
3. along `M=lcm(1,...,B)`, its least witness denominator satisfies
   `B<q_min<=2M`, forcing an asymptotic logarithmic height cost.

The geometry separates as well.  The same address `17/41` remains lonely for
`{1,...,11,13,84(41k+1)}`, while every lonely component has width at most
`6/[7*84(41k+1)]`.  Arithmetic address and geometric thickness are independent
coordinates.

THM-1105's first unblocked modulus `q0` is therefore a useful address
heuristic, not a law.  The exact prefix row at `B=23` has

```text
q0=25,          q_min=53,          q_min-q0=28.
```

That already exceeds the incoming adversarial-search record `19` without a
search.  No bounded-excess conclusion follows.

THM-1110 sharpens the one-modulus picture.  The classical condition “`q`
divides no speed” is sufficient precisely through `q=14`, because only then
is the forbidden residue window `{0}`.  For composite `q`, however, numerator
counts must be stratified by `gcd(v,q)`.  If

```text
W_q = {r : 14*min(r,q-r)<q},
g = gcd(v,q)<q,
```

then the exact number of reduced numerators killed by `v` is

```text
K_(q,g) = [phi(q)/phi(q/g)] * #{r in W_q : gcd(r,q)=g}.
```

At `q=90`, unit speeds kill two numerators but gcd-five speeds kill eight;
`{5,25,35}` covers all 24 units and adding `1` makes the row primitive.  Thus
even a very small local owner set can close one modulus.  The global content
must concern simultaneous ownership across several moduli by one shared
integer word.

## Route B: a chart with affine carry, not a quotient

THM-1099 isolates the threshold error.  Write a rational margin as `s/q` with

```text
q=13s+e,       e>0.
```

Then

```text
1/14 < s/q < 1/13    iff 0<e<s,
s/q = 1/14           iff e=s,
s/q < 1/14           iff e>s.
```

The signed fourteen-gap identity

```text
sum_(i=0)^13 (g_i-s) = q-14s = e-s
```

shows why the old small-gap pigeonhole is threshold-reversed: it works in the
already-safe band and changes sign on the counterexample side.

Even the specialization `e=1` does not force the pure packet
`s*{1,...,12}`.  There is a one-carry family

```text
{s,2s,...,6s,7s+1,...,12s+1}.
```

Nor are these “cosets of `sZ` in `Z/qZ`”.  Since
`gcd(s,13s+1)=1`, multiplication by `s` is a permutation of `Z/qZ`; the
putative subgroup is the whole group.  What exists is an affine chart:

```text
((x+y) mod q) mod s = (x+y-wrap*e) mod s.
```

Each wrap contributes the cocycle `-e`.  This is an archimedean carry, not a
ramified quotient fibre.

Finally, a residue `u mod q` represents the integer lifts `u+kq`.  Covering a
second modulus `d` is governed exactly by

```text
exists k, d | (u+kq)    iff    gcd(q,d) | u.
```

The lift `7 -> 112=7+105` makes the loss concrete.  It preserves the complete
residue multiset, gaps, active pair, and strict local maximum at `8/105` of the
death-star-S57 row.  Yet it changes Cover14 and creates the exact global
maximum `3/20`.  One maximizer chart cannot see the global LRC predicate.

The compatibility equivalence is now kernel-checked in
`LRCRationalScaleGuardrails`:

```text
(exists k : Z, d | u+qk)  iff  (gcd(q,d) : Z) | u.
```

The reverse direction uses explicit Bézout coefficients, so the transition
map is not merely an existential CRT citation.  This is the first formalized
gluing law of the proposed atlas.

## Definition of the candidate carrier

For each denominator `q`, let

```text
U_q = (Z/qZ)^*,
D_(v,q) = {a in U_q : 14*min(av mod q,q-av mod q)<q}.
```

A local lonely address is an element outside `union_v D_(v,q)`.  An integer
speed `v` determines the compatible residue family `(v mod q)_q`.  To compare
charts `q` and `d`, pass to `lcm(q,d)` or retain the lift coordinate in
`v=u+kq`.  Attach to each chart:

- the exact-period address fibre `U_q`;
- the runner-danger incidence masks;
- the owner/carrier divisibility incidence;
- the affine carry word of chosen representatives;
- the signed Fourier labels on the dual fibre; and
- rational witness obligations at the other active moduli.

A global section is a single set of thirteen positive integers satisfying all
these local congruences simultaneously.  This is the **lift-compatible modulus
atlas**.  The name is provisional; no sheaf theorem is being claimed.  It is a
finite shared-assignment hypergraph on every finite modulus block and an
inverse system as the block grows.

This extends the earlier ramified-coset-cover object.  Honest divisor maps use
literal complete fibres.  Nondivisor Route-B charts require affine transition
data.  The global theory needs both, plus signed dual labels.

## Tournament Analysis: useful shadow, wrong carrier

Runners are not the natural tournament vertices here.  Better candidates are

```text
cover obligations, denominator packets, rational addresses, lift sheets,
safe components, Fourier obligations, and owner--carrier incidences.
```

One pairwise observable is least lift cost: orient `d -> e` when the least
available lift carrying owner `d` is smaller than that for `e`.  With the
`(cost,label)` gauge this tournament is tautologically transitive: score
histogram `0,1,...,12`, zero directed cycles, thirteen singleton SCCs, and one
Hamiltonian path.  That makes it good telemetry and useless as a proof by
itself.  It forgets whether two owners use the **same** lift.

A more faithful pair shadow records shared-carrier compatibility, with a switch
given by the chosen modulus chart and a deterministic tie path by modulus
label.  Even that loses higher intersections.  The proof-bearing object is the
owner--carrier hypergraph (or its signed resonance dual); tournaments should
report its projections, edge flips, SCCs, and path counts without being
mistaken for it.

The challenged assumption is therefore explicit: tournament vertices need
not be runners, and a binary relation need not be the final quotient.  The LRC
predicate is existential across addresses and conjunctive across thirteen
speeds, so shared higher incidence is intrinsic.

THM-1121 supplies a second exact warning about vertex choice.  THM-1111's
pairwise maximum-spanning-tree correction is strong but leaves an adversarial
margin.  Replacing runner or killer pairs by a **weighted bipartite incidence
hypergraph** produces a separating functional: 35 rational-time obligations
have total weight 505, every finite-branch killer has load at most 84, and

```text
6*84 < 505.
```

This closes the complete `92<=k_i<333` r=6 finite horn without enumerating
sextuples or even conditioning on the seven-speed core.  Its killer-load
tournament is transitive and adds no proof content.  The useful object is the
fractional cover dual on obligations.  This is a concrete success of the
alternate-vertex-set default: proof obligations, not runners, expose the
missing one-unit inequality.

## Corrections to the current inverse-theorem narrative

The exact logical chain is

```text
AP-core inverse theorem
    => compact Cover14 floor M>=1/13
    => compact LRC residual M>=1/14.
```

No reverse implication is proved.  HYP-7555's `q=183` two-gap picture is a
suboptimal local maximum for its example, whose global maximum is at `q=24`;
it refutes an arbitrary-local-picture implication, not a premise conditioned
on the global maximizer.  HYP-7565's continued fraction forces the `13`
carrier at `14/183`, but its `14` carrier comes from Cover14/AP-core, not from
maximality above `1/14`.

Likewise, HYP-7575's useful conclusion is only that broad interval containment
does not force additive energy.  The band `[s,12s+1]` already assumes the
unproved `e=1` safe-side specialization.  Generic BSG/PFR produces a large
structured subset or a containing generalized progression; an exact
twelve-term AP requires a sharp full-set stability and filling theorem.

The incoming S57 continuation makes a compatible correction: exact
“near-tight implies AP” is false (`{1,...,11,24}` is the first displayed
near-tight non-AP core), and the far-element target is better phrased as
cover-gap uniqueness.  Its soft-Weyl and component-width lenses leave a
very-near-tight fragmented residual.  The broader assertions that every
non-AP core has cover-gap exactly `1/2`, or that primitive equality implies
the AP at all heights, remain verified/conjectural beyond the stated finite
ranges; they should not be used as closed suppliers.

THM-1115 subsequently gives the full thirteen-speed warning at equality:
`{1,...,11,13,24}` is exactly tight, primitive, and not an AP.  THM-1120's
essential-region criterion explains the local substitution `12 <-> 24`, but
its assertion that only two tight families exist is a structured search
result, not a classification theorem.  Any inverse theorem must therefore
state its covering/centering hypotheses precisely; “tight implies AP” is not
available globally.

## Concrete next theorems

### 1. Coherent finite-lift atlas

Fix `q`, unit `a`, residues `r_i`, and canonical speed residues
`u_i=a^(-1)r_i mod q`.  For a finite modulus set `D`, put

```text
K = lcm_(d in D) [d/gcd(q,d)].
```

Prove that all cover data of `v_i=u_i+k_iq` depend only on
`(k_i mod K)_i`, and that owner `d` is carried by lift `i` exactly when

```text
g=gcd(q,d) divides u_i,
k_i = -(u_i/g)*(q/g)^(-1) mod (d/g).
```

Adding the rational danger masks makes this a finite, proof-faithful
shared-assignment problem.  It is also a natural Lean extension of
`LRCRamifiedCosetCover`.

### 2. Safe-band affine-carry normal form

Prove the complete packing classification:

```text
x_1=s, x_(j+1)-x_j>=s, x_12<=12s+e
iff
x_j=js+c_j with 0=c_1<=...<=c_12<=e.
```

For `e=1`, the carry word is zeros followed by ones.  Then separately decide
which carry words admit a shared integer lift satisfying Cover14.  This repairs
the S90/S101 packet step without pretending to prove the global inverse
theorem.

### 3. Cross-modulus signed miss certificate

For each live `q`, expand the exact safe indicators on `U_q`.  If no unit
address survives, the positive zero mode must be cancelled by named
nontrivial character tuples.  Aggregate this identity over an **adaptive**
block of moduli after deleting divisor-dead packets.

The target dichotomy is:

```text
some modulus has a live numerator
or
the shared speed section carries a quantitatively large, explicitly labelled
cross-modulus resonance/owner certificate.
```

This is the global analogue of THM-1091.  It avoids THM-1093's local triangle
loss and THM-1100's fixed-atlas error.

### 4. Uniform tail after the r=6 weighted finite atlas

THM-1121 eliminates the estimated `3.64e12` sextuple finite horn completely.
The remaining issue is not computation below 333 but the logical scope of the
measure half: THM-1102's `max T=308.4` came from a width-16 near-bottom scan.
Interiority inside that window does not prove a uniform statement for widely
spaced or scaled quintuples.  The next theorem must give a scale-separation or
cluster normal form showing that every `k_6>=333` row lies in the measure horn,
or else exhibit the missing tail branch.  The weighted finite atlas should not
be used to conceal that analytic gap.

## Frontier statement

Common-scale ramification is now understood as a nonnested exact quotient
atlas through scale `32`; scale `33` is next.  The r=6 finite horn is closed by
a universal weighted obligation atlas, while its unbounded measure tail is
still open.
Fixed global denominator atlases are impossible.  One Route-B residue chart is
not lift-faithful, and the elementary single-modulus count must be gcd
stratified.  The remaining global object is the signed shared-lift gluing
across moduli.  Formalizing that gluing and extracting a cross-modulus
certificate is the next structural problem; LRC(14) itself is still open.
