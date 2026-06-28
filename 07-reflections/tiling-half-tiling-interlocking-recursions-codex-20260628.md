# Tiling/Half-Tiling Interlocking Recursions

codex-2026-06-28.  The owner asked to think two interlocking recursions: the
tiling model and the half-tiling model.  The useful clarification is that the
two models are not competing encodings.

The fixed-path tiling cube is the witness lift.  It keeps the path, the flip
bits, the canary arc, and the tail/tip deletion coordinates.  The half-tiling
recursion is the quotient/compression fold.  It sees parent automorphism
incident-word orbits, rooted child classes, and the `K_{k,k+1}` coboundary
laws.  The proof object is the span between them, with sidecars.

The n=4 tables make this concrete.  In the fixed-path chart, `a=(0,2)`,
`b=(1,3)`, `c=(0,3)` give fibers

```text
T:{E}, +:{a}, -:{b}, S:{c,ab,ac,bc,abc}.
```

In the half chart, fixed arcs `2->0, 3->0, 1->2, 3->1` give partial score
`[0,1,1,2]`, and free bits `x=(2,3)`, `y=(0,1)` give the four-state section
`E,x,y,xy -> T,+,-,S`.  That is the whole warning: the half chart is clean
because it threw away the four extra `S` witnesses.  Sometimes that is legal;
sometimes it is exactly the missing proof coordinate.

The exact scout `tournament_tiling_half_tiling_interlock_codex_20260628.py`
checks the larger bookkeeping:

```text
fixed-path cover n=6 = 1024
incident-word/rooted orbit count 5->6 = 296
unrooted U(6) = 56
```

and verifies the fixed-path fiber law `H(T)/|Aut(T)|` through `n=6`.  On the
half-tiling side it recovers the `K_{k,k+1}` rank law: `k(k+1)` line
observations have `GF(2)` rank `2k`, with rectangle redundancy `k(k-1)`.

The proof move this suggests:

1. Build witnesses in the tiling cover.
2. Compress through the half-tiling quotient only with named sidecars:
   `H(T)/|Aut(T)|`, parent-aut word orbit, rectangle/hourglass residue, and
   tail/tip deletion signature.
3. Treat every failed descent as a proof obligation, not as noise.

This plugs into the newer LRC-side recursion notes.  HYP-3216 says the proof
route has a cyclotomic ladder and a 2-adic fold.  HYP-3230 says the cap kernel
has a three-gap/Stern-Brocot recursion.  HYP-3231 records the scale-normal
route ledger, HYP-3232 records the modulus-covariance recursion that breaks at
the apex half, HYP-3233 records the cyclotomic-factor layer, HYP-3234 records
the signed-address chart-change sheaf, and HYP-3235 records the totally-real
cyclotomic Fejer-square packet.  HYP-3236 records the Green conductance face,
HYP-3237 records the Vitali-wall/core-construction split, HYP-3238 records the
crossed even-positive/odd-negative compression bridge, HYP-3239 records the
`D_7`/Borsuk-Ulam sign-isotypic refinement and bimodal-discrepancy diagonal,
HYP-3220 records the even-odd/positive-negative parity wall, and HYP-3219
records the Brouwer/sign factorization sidecar.  Incoming HYP-3240 records the
hard-core dilation witness guardrail on the finer `Phi_{14d}` grid, while
incoming HYP-3241 records the equioscillation saddle-index count `(p-1)/2`
and the shared AP/Goddyn-Wong `Phi_14` core.  The topological degree sidecar
now has a concrete antipodal-pair count.  Incoming HYP-3242 identifies the cap
as a measured Euler characteristic of the danger-cover nerve, and HYP-3243
packages the geometry as a typed proof-route atlas.  HYP-3244 adds the
tournament-side lift/compress recursion: tiling for explicit witnesses,
half-tiling for legal quotients.

The next concrete test is finite: attach this span to HYP-3227's
trap-discharge graph, HYP-3236's Green-resistance traps, HYP-3238's
odd/negative payload audit, HYP-3239's `D_7`/bimodality audit, HYP-3220's
parity/sign sidecar, HYP-3219's sign sidecar, HYP-3237's core witnesses, and
HYP-3218's Fejer/equidistribution proof-push.  Every trap edge should declare
its tiling lift, half-tiling
descent certificate, and failed sidecar if the quotient is not legal.  If that
closes, the current finite trap computations become a small theorem schema
rather than another scalar atlas.

After the HYP-3238 scout, the sidecar test has a concrete warning label:
zero negative covariance has `18` primitive non-AP false terminals, and the
`11` non-AP exchange traps split into `8` negative-leakage-plus-odd-`q3` debts
and `3` odd-`q3` debts without negative leakage.  The tiling/half-tiling
descent should therefore record whether endpoint-bimodality pricing has paid
the odd `q3` debt before accepting any positive/even quotient.

Related: HYP-3244, HYP-3243, HYP-3242, HYP-3241, HYP-3240, HYP-3239, HYP-3238, HYP-3237, HYP-3236, HYP-3235, HYP-3234,
HYP-3233, HYP-3232, HYP-3231, HYP-3230, HYP-3229, HYP-3227, HYP-3220,
HYP-3219, HYP-3218, HYP-3216, HYP-3053, HYP-3149, HYP-3199, OPEN-Q-108.
