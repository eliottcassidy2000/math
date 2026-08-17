---
id: THM-2473
title: "The sporadic Keller map's branch tower: depressed-trisection normal form, even escape parity, non-surjectivity, and wreath monodromy"
status: >
  PROVED for the exact algebraic claims (1)-(6) (symbolic identities and
  exact rational fibers, three companion scripts); VERIFIED for the
  full-wreath monodromy of F o F (Chebotarev census, 10 primes, 6.08M
  targets; subsequently PROVED at every depth by THM-3535); FINITE-EXACT
  for the tower census through depth 3;
  CONDITIONAL (off-spine genericity) for the closed-form tower law (7).
source: opus-2026-07-26
related:
  - THM-1300 (the sporadic Keller/JC>=3 counterexample map)
  - HYP-9030-keller-degree-semigroup (tests (i) and (iii) discharged here)
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-1335-trisection-modulus-master-identities-trace-polynomiality
  - HYP-9029-strong-interval-tiling-law
scripts:
  - 04-computation/keller_FoF_ternary_census_opus_20260726b.py
  - 04-computation/keller_plane_anatomy_opus_20260726.py
  - 04-computation/keller_core_cubic_opus_20260726.py
outputs:
  - 05-knowledge/results/keller_FoF_ternary_census_opus_20260726b.out
  - 05-knowledge/results/keller_plane_anatomy_opus_20260726.out
  - 05-knowledge/results/keller_core_cubic_opus_20260726.out
---

# THM-2473 -- the branch tower of the sporadic Keller map, exactly

Throughout, `F` is the Claude-discovered Keller map (THM-1300;
provenance per MISTAKE-205): `u = 1+xy`,

```text
F = ( u^3 z + y^2 u (4+3xy),  y + 3x u^2 z + 3x y^2 (4+3xy),  2x - 3x^2 y - x^3 z ),
det J_F = -2  (F is etale everywhere: NO ramification exists anywhere).
```

Target coordinates are `(a,b,c)`. Since `F` is etale, the honest analogue
of a "branch locus" is the **Jelonek non-properness set** `S_F` (fiber
drop by escape to infinity). This file computes the entire tower exactly.

## (1) Two contracted divisors and the divisor-level `1+2` law [PROVED]

`F1 = u * (u^2 z + y^2(4+3xy))` and `F3 = x * (2 - 3xy - x^2 z)` exactly.
Hence `F` maps the quadric `{u=0}` into the plane `{a=0}` and the plane
`{x=0}` into the plane `{c=0}`, each **birationally** (explicit inverses:
`(0,b,c) -> (2/b, -b/2, b^2(10-bc)/8)` and `(a,b,0) -> (0, b, a-4b^2)`). Over a generic point of each special plane
the residual fiber is an explicit **2:1 fold**:

```text
fiber over (0,b,c) = 1 (u=0 sheet) + 2 (roots of Q = -b(bc-1)x^2 - 2(bc-1)x - c)
fiber over (a,b,0) = 1 (x=0 sheet) + 2 (x = +-2/sqrt(b^2-16a))
```

**The `1+2` motif of the trisection anatomy `(T-+1)(2T+-1)^2` is the
pullback-divisor splitting of the two special coordinate planes: one
birational (anchored) sheet plus one folded 2:1 residual.**

## (2) The depressed-trisection normal form [PROVED]

The generic x-eliminant of the fiber `F^{-1}(a,b,c)` is
`a * x^6 * E(x; a,b,c)` with **core cubic**

```text
E(x) = L x^3 + (4 - 3bc) x - 2c,      L = 27a^2c^2 - 18abc + 16a + b^3c - b^2,
```

`L` irreducible, and the `x^2` coefficient **identically zero**. Hence
the three fiber x-coordinates satisfy `x_1+x_2+x_3 = 0` for EVERY target
(collision fiber `{0, +1, -1}` included), and the fiber structure of `F`
is a global pencil of **depressed cubics** -- the exact shape of the
Chebyshev trisection `4T^3 - 3T = const`. This is the mechanism behind
klein-S324 / HYP-9030's branch-anatomy observation.

## (3) `S_F = {L = 0}` and the even escape parity law [PROVED]

- Off `{L=0}` every fiber has exactly 3 points; since `F` is etale,
  full count forces properness (local diffeomorphism at each fiber point
  + degree bound), so `S_F = Z(L)`, an irreducible hypersurface
  (Jelonek codim-1 purity honored).
- On `{L=0}` generic (`3bc != 4`): because the `x^2` coefficient is
  identically zero, the eliminant degree drops from 3 straight to 1:
  the escaping pair is the fold `x ~ +-sqrt(-(4-3bc)/L)` and the unique
  survivor is `x = 2c/(4-3bc)`. **Every non-proper fiber drops 3 -> 1:
  escape parity is even; a lone point never escapes.** (Census
  cross-check: zero `n1=2` fibers among 6,080,168 mod-p targets.)
- Restrictions reproduce the in-plane drop loci exactly:
  `L|_{a=0} = b^2(bc-1)` (the `b=0` line and the hyperbola `bc=1`) and
  `L|_{c=0} = 16a - b^2` (the parabola).
- Exceptional stratum: on the rational curve
  `t |-> (4/(27t^2), 4/(3t), t)` (inside `{L=0, 3bc=4}`) the eliminant
  is a nonzero constant: the fiber is **EMPTY**. Hence
  **`F` is NOT surjective**: the sporadic Keller counterexample omits an
  explicit rational curve. (Witness `(4/27, 4/3, 1)`, exact.)
- The square factor of the discriminant,
  `disc_x E = -4 (27ac^2 - 9bc + 8)^2 L`, cuts the surface where the
  x-PROJECTION collides while the fiber stays full (witness
  `(1/27,1,1)`: eliminant `(x+6)(x-3)^2 x^6`, the double root `x=3`
  carrying a conjugate y-pair) -- no etale violation.

## (4) Axis dynamics and the rational collision family [PROVED]

`F(t,0,0) = (0,0,2t)` and `F(0,0,w) = (w,0,0)`: the two axes form a
2-cycle and `F^2` restricted to the x-axis is literal **doubling** --
the "torus doubling lambda -> lambda^2" of the collision census made
exact. The famous triple collision is the `m=1` member of an infinite
rational family: for every integer `m >= 1`,

```text
F^{-1}(-1/(4m^2), 0, 0) = { (0,0,-1/(4m^2)) } u { (+-m, -+3/(2m), 13/(2m^2)) }.
```

(The `z`-coordinate `13/(2m^2)` reproducing the LRC prime 13 is FLAGGED
as unadjudicated numerology, in klein-S324's spirit; not used anywhere.)

## (5) Monodromy of F [PROVED] and of F o F [structural + VERIFIED]

- `Mon_geom(F) = S3`: the core cubic is irreducible over `C(a,b,c)`
  (the v* fiber has three distinct x-coordinates, so the x-coordinate
  generates the degree-3 extension), and the restriction of the cover to
  the line `(a,0,0)` contains the explicit fold `4a x^2 + 1 = 0` with
  discriminant `-16a`, a nonsquare in `C(a)`: a transposition.
  Transitive + transposition + degree 3 = `S3`. The arithmetic group is
  also `S3` (specialization at `P- = (-1,3/2,13/2)`: fiber cubic
  `21119x^3 - 404x - 208`, irreducible, disc `-515429008338944`
  nonsquare).
- `Mon(F o F) <= S3 wr S3` structurally (the tower factors). The block
  kernel contains `S3 x S3` exactly: the two depth-2 fiber cubics over
  `P+` and `P-` generate NON-isomorphic cubic fields (discriminant ratio
  nonsquare). Full wreath (order 1296) is **VERIFIED** by the Chebotarev
  census: the joint distribution of rational-fiber sizes `(n1, n2)` of
  `(F, F o F)` over ALL targets in `F_p^3` for ten primes 61..103
  matches the independent-uniform wreath model in all 13 cells (ratios
  0.92-1.02), with no off-model cell ever observed.
- THM-3535 subsequently upgrades this evidence: THM-3530's newest image
  prime isolates one bottom-block inertia transposition, which forces the
  full wreath product at depth two and inductively at every depth.
- **HYP-9030 test (i) verdict: the explicit degree-9 member F o F is
  maximally NON-atomic** (imprimitive with full kernel).  This remains an
  exact statement about the `F`-submonoid, but THM-3438 later refutes the
  global atom-spectrum conjecture: degrees `4` through `8` are realized and
  automatically atomic, while degree one is the unit rather than an atom.

## (6) The depth-2 ternary tree over v* [FINITE-EXACT]

`v* = (-1/4,0,0)`; level 1 = `{P0=(0,0,-1/4)} u {P+- = (+-1,-+3/2,13/2)}`
(the `1+2` sigma-orbit). Depth 2 has **7 points, not 9**:

```text
             v*                          (on the x-axis; NOT in S_F)
        /    |    \
      P-    P0    P+                     P0 is ON S_F (z-axis drop line)
     /|\    |     /|\
    3   (1 only)   3                     7 = 3 + 1 + 3
 cubic  2 escaped  cubic
 K- = Q[x]/(21119x^3-404x-208)    K+ = Q[x]/(20929x^3+532x-208),  K- != K+
```

The naive self-similar "sigma-equivariant 1+2 iterate tree" of HYP-9030
is **REFUTED in its strong form**: the equivariance TWISTS
(`F o sigma = tau o F` with `sigma = diag(-1,-1,1)`,
`tau = diag(1,-1,-1)`, `sigma != tau`, and no sign-diagonal `rho`
satisfies `F o rho = sigma o F`), the anchored branch `P0` absorbs the
Jelonek degeneration, and the two folded branches live in different
cubic fields. The true iterate grammar is: **anchored spine (the doubling
axis 2-cycle, alternating on/off `S_F`) + generic ternary branching off
the spine.**

## (7) The tower census law [depths 1-3 FINITE-EXACT; closed form CONDITIONAL]

`N_1 = 3, N_2 = 7, N_3 = 21` (depth 3: the axis point `(-1/8,0,0)` is
off `S_F` and contributes `1 + sqrt2-pair = 3`; all six cubic-field
points lie off `S_F`, checked by nonzero norms of `L`, contributing 18).
The spine alternates z-axis-type (drop, 1 child) and x-axis-type
(3 children) targets, so conditionally on all off-spine points staying
off `Z(L)` (verified through depth 3),

```text
N_{k+1} = 3 N_k - 2 [spine node of level k is z-axis type], i.e.
N_k = (3^{k+1} + 3^(k mod 2)) / 4 :   3, 7, 21, 61, 183, 547, ...
```

**HYP-9030 test (iii) verdict: the ternary collision census does NOT
match the naive `1+2` iterate tree (9, 27, ...); it matches the spine-
degenerate law above. The "branch tower" of an etale Keller map is its
Jelonek tower, and the degeneration always sits on the anchored branch.**

## Loss ledger / hostile controls

- In this theorem alone, the mod-p census was Chebotarev evidence and the
  exact kernel lower bound was `S3 x S3`; THM-3535 subsequently proves the
  full wreath group at every depth by localized inertia.
- The closed form (7) beyond depth 3 needs each new folded point checked
  off `Z(L)`; a folded point landing on `Z(L)` would break the law at
  that node (no such point found through depth 3).
- Nothing here bears on `JC(2)`/`DC(2)`; all statements are about the
  fixed sporadic map `F` (and its self-compositions), not the family.
- The `13/(2m^2)` observation is flagged, unadjudicated, and unused.
