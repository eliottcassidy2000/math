# The allocation word does not lift to an endpoint `K4`

**Status: FINITE-EXACT PARTIAL SCOUT / NON-CANONICAL SYNTHESIS.**  The exact
companion is
[`lrc14_endpoint_origin_k4_horn_partial_scout_20260803.py`](../04-computation/lrc14_endpoint_origin_k4_horn_partial_scout_20260803.py)
with its frozen
[`output`](../05-knowledge/results/lrc14_endpoint_origin_k4_horn_partial_scout_20260803.out).
The result selects outcome 2 of THM-3285's proposed endpoint test: the middle
endpoint fibre is present, but the outer typed companions and the full-depth
vertical lift fail.  It is not a theorem against new endpoint couplings and
has no reserved theorem ID.

## Inheritance pass

The closest proved mechanism is
[THM-2772](../01-canon/theorems/THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction.md):
at zero source and target harmonics, a **common endpoint atom** may carry the
four states bare/source/target/both.  The endpoint pair `(L,R)`, target
difference `q=L-R`, allocation bits, clock, and semantic word are all
load-bearing.  A common cardinality or a determinant label does not attach
that fibre to a physical carrier.

The corrected near miss is
[THM-2807](../01-canon/theorems/THM-2807-positive-graded-address-two-simplex-and-allocation-lift-boundary.md):
the three positive whole cylinders form an exact address triangle, but it
explicitly stops before endpoint origin and allocation.  Its vertical edge
is `Z2^4079` and collapses modulo `169`.

The canonical hostile is
[THM-2820](../01-canon/theorems/THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go.md):
a target-relative endpoint motion can be nonzero while common source support
is absent, whereas a lawful simultaneous rechart can be pure gauge.  The
least-used relevant sidecar is
[THM-2813](../01-canon/theorems/THM-2813-affine-lift-transvection-and-projective-horn-decoder.md):
every point on the residue-seven horn is fixed by the relative lift, and one
off-sheet normal jet is the minimal coordinate that can recover its label.

The live concept board was therefore:

1. THM-3285's whole-cylinder `R-M-R` word;
2. the six separately typed semantic factors;
3. THM-2772's 169-origin endpoint fibre and four allocation states;
4. raw common co-support versus a post-carrier-DFT central `K4`; and
5. the full-depth `Z2^4079` edge versus its collapsed mod-`169` shadow.

## The lawful partial universe

Fix the canonical horn label

```text
(physical clock,sigma,tau)=(1,0,3)
```

and the three addresses

```text
n0=3454614,   n+=3454627,   na=4143978.
```

The scout first reconstructs THM-3285's exact ambient result:

```text
allocation word: R -> M -> R,
uncut target word: B -> B -> B,
weight:           27581135604,
mass:             60781651775958960/371293,
carry-6 value:    790161473087466480.
```

The endpoint-current constructor takes integral intervals, while the narrow
address cylinders have rational endpoints.  No canonical whole-cylinder
endpoint pushforward is present in canon.  Rather than invent one, the scout
uses the least integer unit interval strictly inside each cylinder.  The
three target selectors are

```text
[142082432180573,142082432180574),
[142088047310093,142088047310094),
[142004622528653,142004622528654).
```

They and their source pullbacks follow all three exact horn translations.
Thus this is a reproducible subatom test inside the unchanged whole-cylinder
universe, not a replacement for it.

## First failure: `R` is not an allocation-absence bit

At every vertex, the raw source and target one-sided carrier twist masks are
both the same delta mask supported at harmonic zero.  In particular, the raw
source carrier is present at the two `R` vertices.  What changes is the typed
semantic section.  In factor order

```text
(E3,clock,q1,q2,c2,c3),
```

the source signatures are

```text
n0: (0,1,1,1,0,1),
n+: (1,1,1,1,1,1),
na: (0,1,1,1,1,1),
```

while every target signature is `(1,1,1,1,1,1)`.  Therefore the outer
right-cofiber symbol means “the fully typed source section is absent,” not
“toggle the source carrier off while retaining the same endpoint atom.”
The first invalid implication is

```text
R in the M/R support decomposition
  does not imply
source-absent in THM-2772's bare/source/target/both square.
```

This is the missing coordinate that a syntax-only `R -> M -> R` to `K4`
dictionary would erase.

## Exact endpoint census

The endpoint reconstruction exhausts all `169` origins in both certified
exact-order fields.  The underlying left/right support counts are

```text
             n0       n+       na
left       169/169  169/169    0/169
right      169/169  169/169  169/169
```

At the middle, explicit origin `(0,0)` has four nonzero central allocation
coefficients in both fields.  For the first field, the
`(bare,source,target,both)` product vector is

```text
(85548272376494745,
 33683794039955122,
 33683794039955122,
 56797376486599908),
```

and the second-field vector is also entrywise nonzero.  Thus outcome 1—an
empty middle endpoint fibre—is rejected sharply.

This does not pay THM-2772's mixed-face invoice.  At the sole primal twist
supporting all four states, the vector is `(w,w,w,w)` and its pointwise
Möbius face is zero.  The displayed nonzero central `D3` is again the
post-carrier-DFT complement sum classified by THM-2806, not a nonflat raw
common atom.

At `n0`, the algebraic endpoint values exist but the missing `E3,c2` gates
kill both source amplitudes.  At `na`, the boundary is even sharper: the
underlying left endpoint bank is zero at all `169` origins in both fields,
and `E3` is separately absent.  The target `R` cylinder remains one positive
whole cylinder throughout.  Hence the minimal hostile is a positive target
whole cylinder with all target factors and no source endpoint companion.

## The full-depth edge also fails

The address arithmetic remains

```text
na-n+=169*4079,
n+=na=98 mod169.
```

All three target endpoint banks have full support.  The scout exhausts the
natural affine endpoint covariance class

```text
F'(x)=c chi_u(x) F(x+s),
s,u in F13^2.
```

The first edge `n0 -> n+` has exactly one covariance in each field: zero
shift, trivial character, and one nonzero scalar.  This is the positive
control.  The vertical `n+ -> na` edge has no such covariance in either
field, and neither does the diagonal `n0 -> na`.  Thus even the surviving
target endpoint bank does not transport through the inherited affine-gauge
class along `Z2^4079`.

This is not a no-go for an arbitrary new correspondence.  Canon supplies no
map from the horn's address/ancestry data to a common `(L,R,q,Delta)` atom,
and THM-2813 already predicts that on-sheet data cannot recover the relative
lift without a normal sidecar.

## Connection contract and stopping reason

```text
source:
  the canonical THM-3285 label and one integral unit selector inside each
  of its three translated whole cylinders;

target:
  all 169 left/right endpoint origins and the fixed-sheet central
  bare/source/target/both products in two exact-order fields;

map tested:
  physical source/target restriction, endpoint relation DFT, then every
  affine shift/character/scalar covariance of the target origin plane;

preserved:
  ambient whole cylinder, address translations, raw carrier twist origin,
  clock, sigma, tau, root, carry, and exact open-boundary convention;

destroyed or absent:
  outer source E3 (and n0 c2), the far left endpoint bank, a common
  `(L,R,q,Delta)` atom, raw nonzero mixed face, and the vertical endpoint
  intertwiner;

verdict:
  outcome 2 -- middle fibre present, endpoint/allocation lift fails.
```

The cheapest next positive test is not another origin choice: all `169` were
already exhausted.  It must add a typed source/target connection that
transports `E3` and the endpoint present selector, or move to an off-sheet
atom carrying THM-2813's normal jet.  Either move is genuinely new sidecar
data.  Without it, further scalar or origin scans only repeat this stopping
mechanism.

No root/Čech map, row exclusion, physical current cycle, or `LRC(14)`
conclusion follows.
