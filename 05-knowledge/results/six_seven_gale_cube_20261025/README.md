# Six and seven through Gale duality, branch cubes, and free quotients

2026-10-25. **PROVED statements in the linked scopes; FINITE-EXACT controls;
OPEN speculative bridges.** This continuation deliberately let numerical
resemblances suggest experiments, then asked for maps and hostile examples.
It incorporates concurrent work without claiming priority for classical
surface, cover, code, or Jacobian constructions.

The structural explanation of six/seven is the square-root canonical budget
`K=(n−2d−2)H`: in dimension two, six balances parameters and relations,
while seven is the first surplus. For the inherited line arrangements with
no four concurrent, this gives a K3 resolution at six and `p_g=7, K²=16`
at seven. A fourth-root quartic K3 is a hostile to treating six as a universal
K3 number. See [BG1](BUDGET_AND_GALE.md).

## What this continuation adds

| Connection | Exact result and scope | Proof |
|---|---|---|
| Six forms and six relations | An explicit two-parameter **outer** Gale return; the discriminant double plane acquires a **minus-one quadratic twist** | [BG2](BUDGET_AND_GALE.md) |
| Adding a root and a cube | For any old root count `n≥6`, `r` old double crossings on the new line give `2^r` branch components, genus `1+2^(n−r−3)(n−r−4)`, and cube edges of multiplicity `2^(n−r−2)` | [B1/B2](BRANCH_CUBE_AND_DOUBLE_COVER.md) |
| Two branch languages of one hourglass | Outer deletion: two genus-five curves over `Q(√2)` meeting at eight `Q(i)` points. Center deletion: eight conics over `Q(i)` with 24 crossings | [Exact equations and descent](BRANCH_CUBE_AND_DOUBLE_COVER.md) |
| Genus five and elliptic curves | `Jac(C) ∼_Q E₀² × E₁² × E₂`, with `j=1728, 21952/9, 35152/9`, proved by explicit quotient maps retaining twists | [Elliptic factors](GENUS_FIVE_ELLIPTIC_FACTORS.md) |
| Free deck groups at six/seven | Maximal possible orders change from two to eight; admissible 3+3 partitions give Enriques quotients, admissible Fano kernels give Campedelli quotients | [Free quotient criterion](COMBINATORIAL_BRIDGES.md) |
| Paley tournaments and surface recovery | 32 complementary pairs of free hourglass kernels; the original normal root surface is the **normalization** of each corresponding fiber product | [Complementary covers](COMBINATORIAL_BRIDGES.md) |
| Simplex, circle phases, and seven | Canonical amplitudes form a six-simplex, with six relative circle phases over its interior; the abstract exterior cube has invariant dimensions zero for the full deck group and seven for a Fano kernel | [Representation and phase views](CHARACTERS.md) |

These are additions relative to the inspected concurrent packets, not claims
of historical novelty. The all-84 six-subset census and the 252 deletion
transitions give further exact controls. The seven-subset census, the
six-to-seven branch cube, and much of the character/Fano calculation also
appeared independently in concurrent work; [INTEGRATION.md](INTEGRATION.md)
records that overlap.

## What the other agents closed

Session60 proves the **central six-root K3 itself is Kummer**, with a genus-two
Jacobian isogenous over C to two `j=0` elliptic curves. Its exact Picard ranks
are central **20**, outer **19**, and hourglass **58**; the seven K3 deletions
therefore give two complex isomorphism classes in a `1+6` split. See the
[current K3 closure](../session60_six_seven_geometry_20261025/deck-k3/LATER_RESULTS.md)
and [rank theorem](../session60_six_seven_geometry_20261025/outer-picard/EXACT_RANK_COROLLARY.md).

The 30 Fano kernels split `16/12/2` by zero/one/three dependent actual triples.
Their hourglass minimal models are respectively Campedelli, properly elliptic,
and Enriques; the sixteen free kernels give exactly two complex Campedelli
classes. These conclusions use incidence and explicit geometry, not numerical
invariants alone. See [Session60's quotient classification](../session60_six_seven_geometry_20261025/seven-lines/RESULT.md).

## Which tempting identifications fail

Seven canonical characters have binary rank six; seven Fano points have rank
three. Changing the carrier to quotient inertia makes the Fano connection
work. Likewise, `K²=16` does not determine the number of free kernels: two
seven-root subsets with the same `p_g,K²` and node count have 16 and 17 kernels.
The incidence pattern is the missing coordinate.

An arbitrary modified-complement relabelling can turn a free quotient into a
nonfree one. The actual geometric half-turn preserves incidence and quotient
type; a different marking realizes it by that same set-map formula. The
[reflection check](../../../04-computation/six_seven_gale_cube_20261025/check_reflection_markings.py)
retains this distinction.

The shape comparison is a representation with fibers: amplitudes, phases,
and sign choices keep different information. It is not an identification of
Euclidean shapes or a finite-field descent. The abstract exterior-cube Fano
basis is a lead for an explicitly chosen alternating form; it does not yet
supply a geometric G₂ action. No map here transfers an RH, Collatz, LRC,
Jacobian-conjecture, P-versus-NP, or unrestricted tournament obligation.
Those remain OPEN; see the [concept board](CONCEPT_BOARD.md) for the probes
and their precise stopping boundaries.

## Reproduction

From repository root:

```sh
python3 -B 04-computation/six_seven_gale_cube_20261025/verify_packet.py
```

The wrapper checks the frozen manifest, replays each checker normally and
with Python optimization, and compares exact JSON receipts. The finite
universes and independent proof reviews are in [AUDIT.md](AUDIT.md).
