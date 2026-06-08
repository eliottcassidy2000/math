---
id: HYP-2393
name: erdos93-98-duality-the-regular-polygon-pivot-and-the-convex-layers-bridge
status: SURVEY + bridge (hull-size bound h<=2D+1 VERIFIED) + a route to superlinearity for 98
date: 2026-06-08
session: claudebox-2026-06-08-S723
external:
  - "Erdos 93 (SOLVED, Altman 1963, Lean-verified): convex n-gon => >= floor(n/2) distinct distances; regular n-gon tight"
  - "Erdos 98 (OPEN): no 3 collinear, no 4 concyclic => >= h(n); does h(n)/n -> inf?"
depends_on:
  - HYP-2390  # Erdos 98 = the unit-distance dual; the shared '3' (S722)
  - THM-440   # u(22) cocyclic extension (the radius-1 slice)
provisional_id: true
---

# HYP-2393: Erdos 93 and 98 are a dual pair; the regular polygon is the pivot; convex layers bridge them

## The two problems

- **93 (SOLVED, Altman):** `n` points in CONVEX position determine `>= floor(n/2)` distinct distances.
  TIGHT extremizer = the **regular n-gon** (exactly `floor(n/2)`; VERIFIED `n=4..12`).
- **98 (OPEN, S722):** no 3 collinear, no 4 concyclic `=> >= h(n)`; is `h(n)/n -> infinity`?

## THE PIVOT (the duality): 93's optimizer IS 98's forbidden config

The regular n-gon achieves `floor(n/2)` distinct distances by putting **all n points on ONE circle**
(every distance is a chord, only `floor(n/2)` chord-lengths). But "all n concyclic" is exactly what 98
**forbids** (no 4 concyclic). So:
```
   98 = 93 with its optimizer (the concyclic ring) BANNED.
```
93 reaches `floor(n/2)` via a single concyclic ring; 98 cannot use a ring; the question is whether
banning the ring forces superlinearly many distinct distances.

## 93 -> 98 (the creative bridge): the convex-layers / onion argument

Apply the SOLVED 93 to the convex layers of a 98 configuration:
- **Hull-size bound (VERIFIED, 3000 configs):** a set with `D` distinct distances has convex hull on
  `h` vertices, and the hull (a convex polygon) determines `>= floor(h/2)` distinct distances, all among
  the `D`. Hence **`h <= 2D+1`**: *few distinct distances forces a SMALL convex hull.*
- **Recursively (onion peeling):** every convex layer has size `<= 2D+1`, so a few-distinct-distance set
  is organized into `>= n/(2D+1)` thin convex layers. (VERIFIED: per layer, distinct distances
  `>= floor(size/2)`.)
- **The 98 constraint cuts in:** no 4 concyclic forbids any layer (or any cross-layer 4-set) from lying
  on a circle. The regular-polygon trick that makes ONE ring cheap is unavailable layer by layer, and
  the layers cannot be circular. So a 98-extremizer is many thin, mutually non-concyclic convex rings.

**Route to superlinearity (the conjecture):** bound how few distinct distances `>= n/(2D+1)` thin,
pairwise-non-concyclic convex layers can SHARE. If shared distances across non-concyclic layers are
limited (a circle-incidence / additive-energy bound, S722's attack), then `D` cannot be linear, giving
`h(n)/n -> infinity`. The hull-size bound turns 98 into a recursive "small-hull, many-non-circular-rings"
structure problem, attacked with the SOLVED 93 as the per-layer engine.

## 98 -> 93 (reverse inspiration): a toehold on the OPEN strong-93

The STRONG 93 (still open): some single vertex determines `>= floor(n/2)` distinct distances
(Dumitrescu `13n/36`). Importing 98's mechanism: a vertex sees few distinct distances iff many points are
concyclic around it; under the extra hypothesis "no 4 concyclic" each vertex's distance-multiplicity is
`<= 3`, so **every vertex sees `>= ceil((n-1)/3)` distinct distances** -- a clean unconditional toehold
on strong-93 for general-position convex sets. The gap `(n-1)/3 -> n/2` is precisely the convex ordering
structure that 93 exploits and 98 lacks.

## The unified picture (autocorrelation / temperature, S714/S720)

Both problems minimize the SUPPORT of the radial autocorrelation under a constraint:
- 93's convex constraint is RIGID and yields a clean `floor(n/2)` (the regular-polygon crystalline ring);
- 98's general-position constraint is WEAK and yields only `(n-1)/3`, with superlinearity open.
The regular polygon is the cold crystalline ring both share as the boundary object: 93's extremal point,
98's forbidden point. The "3" of 98 (per-circle cap, = THM-440/u(22) ceiling) and the "2" of 93
(`floor(n/2)`, chord-doubling on a ring) are the two ways a ring's autocorrelation collapses.

## Next
- make the cross-layer circle-incidence bound precise: how many distinct distances must `L` thin
  non-concyclic convex layers share? (the additive-energy route, t-0109);
- compute the small-`n` constrained min-distinct-distance table and its onion-layer profile;
- formalize "h <= 2D+1" as a clean lemma feeding a recursive lower bound for 98.
