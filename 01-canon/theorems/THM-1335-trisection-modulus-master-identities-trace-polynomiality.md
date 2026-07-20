# THM-1335 — The trisection modulus, the master identities, and trace-polynomiality of the Jacobian counterexample

**Status:** VERIFIED-EXACT (symbolic identities + 16-sample 9×-overdetermined exact
fit for Tr(z); marked items pending). **Author:** kind-pasteur-2026-07-20-S128c99
(HYP-8130). Pulls: klein-S324 (u-side cubic, T₃(W)=1), klein-S325 (Smith rule),
mac-mini THM-1315 (surjectivity), my THM-1310 (x-side geometry), boxeph-S142,
death-star THM-1305/1320/1325.

Setting: F = the owner's JC₃ counterexample; targets (a,b,g); u = 1+xy;
L = 27a²g² − 18abg + 16a + b³g − b² (THM-1310: leading coeff, resolvent −L,
Jelonek set {L=0}); Q = 27ag² − 9bg + 8 (x-side fold wall).

## (1) The u-side fiber cubic is ALSO depressed

Res_x(C1, C2) = a·u³·(L·u³ + (b²−12a)·u − 4a) — klein-S324's cubic re-derived
from the THM-1310 conic pair. No u² term: **Tr(u) = 0 alongside Tr(x) = 0** —
both natural fiber coordinates have vanishing trace.

## (2) THE MASTER IDENTITIES (the whole map in two lines)

Substituting F's own components (a,b,g) := (F₁,F₂,F₃) — verified IDENTICALLY on ℂ³:

    L(F)·x³ + (4 − 3F₂F₃)·x − 2F₃ ≡ 0
    L(F)·u³ + (F₂² − 12F₁)·u − 4F₁ ≡ 0

Every fiber statement, the monodromy, the walls, and the traces are corollaries
of these two polynomial identities. They are the sharpest Lean targets in the
thread.

## (3) The Chebyshev modulus and the cusp-shaped Jelonek set

Scaling the depressed u-cubic to trisection normal form 4T³ − 3T = m (u = rT):

    m² = 108·a²·L / (12a − b²)³ ,   and the PERFECT-SQUARE surprise:
    m² − 1 = E² / (12a − b²)³ ,     E := 54a²g − 18ab + b³.

Equivalently the polynomial identity  **108·a²·L = (12a − b²)³ + E²** : the
Jelonek/resolvent quartic L is itself a (cube + square)-discriminant object.
Consequences:
- **Jelonek set as a cusp:** {L = 0} = closure of {P³ + E² = 0} ∖ {a = 0} with
  P := 12a − b². The asymptotic variety is a CUSPIDAL hypersurface in the
  (P, E)-net — the map is **A₂-shaped at both ends** (source ruling cone
  b³ + 27a²c = 0 = cone over cuspidal cubic; target Jelonek = cusp;
  monodromy W(A₂) = S₃ in between). The W-ADE reading is now
  structurally evidenced at d = 3, not just numerological.
- u-side discriminant Δᵤ = −4E²L: same resolvent class √(−L) as the x-side
  Δₓ = −4Q²L ✓; {E = 0} is the u-side fold wall and exactly the m = ±1 locus.
- At the collision target m = 1: the fiber is the TRISECTION OF ANGLE 0 —
  u-values (1, −1/2, −1/2) = cos(0, 2π/3, 4π/3), klein's T₃(W) = 1 as the
  m = 1 specialization. Generic fibers = {r·cos((θ + 2πj)/3)}: **F is a
  √(−L)-twisted pullback of the Chebyshev trisection cover** — JC₃ fails
  polynomially the way angle trisection fails: cos(θ/3) is not rational in
  cos θ.
- m is a TARGET-side function ⟹ invariant under source-conjugation: klein-S326's
  G1 (radical 1+xy+x³) and G2 (Nagata) have identically the same m — the
  **essential-class invariant** for the realization program (T1549).

## (4) TRACE-POLYNOMIALITY (the centroid theorem)

Fiber-traces over the target, computed by multiplication-trace in ℚ(a,b,g)[x]/(N):

    Tr(x) = 0                       (exact, = depression)
    Tr(u) = 0                       (exact, = depression)
    Tr(y) = 3b/2                    (exact symbolic — POLYNOMIAL)
    Tr(z) = −81a²g²/2 + 27abg − 51a + 15b³g/8 − 3b²/4
                                    (exact rational fit at 16 random targets,
                                     9× overdetermined, residuals ≡ 0; weight
                                     (−2) and τ-parity constrain the basis;
                                     symbolic confirmation pending)
    Tr(x²) = 2(3bg−4)/L,  Tr(x³) = 6g/L    — genuine L-poles.

So: **every source COORDINATE is trace-polynomial although ℂ[x,y,z] is not
finite over ℂ[F]** (non-properness = {L≠0}-only) and the trace form as a whole
has L-poles from degree 2 on. The fiber-centroid map

    𝒞(a,b,g) = (0, b/2, Tr(z)/3)

is an honest polynomial map of rank 2 with image in the plane {x = 0}: fibers
are balanced across the plane x = 0 at height b/2 in y. SHARPER-JC CANDIDATE
(open, falsifiable): *for every Keller map there exist source coordinates whose
fiber-traces are polynomial in the target* — the "polynomial centroid" law.
(Anchors: the a-axis value Tr(z) = −51a matches the collision family exactly:
−1/(4t²) + 13/(2t²) + 13/(2t²) = 51/(4t²).)

## (5) Rigidity at the widened ansatz; the refined hunt

Extended k=3 control (κ free in A₂ = κu²x, S free in B₂ = S + κxW, wide R):
the full nonzero-det solution variety is a 5-parameter family that decodes
ENTIRELY into the trivial group — r₁₁ = source/target scaling, κ = target
scaling of the second coordinate, w₂₂ = z-translation, s₀₀/r₀₃ = target
translations — around the owner's map (κ=3, r₁₁=1: W = y²(4+3xy), S = y,
R = xu(u+1), det −2). **Rigidity modulo trivial operations at the widest box
yet.** Refined k=4 (κ free, wide R incl. x²u-stratification), k=4 two-term
middle, and k=5: ALL **EMPTY**, certified by Rabinowitsch + Gröbner over
GF(32003) — GB = [1] for (Keller system + t·det(0) − 1) — with the k=3 system
CONSISTENT as control (GB size 28). Specialization caveat: mod-p emptiness is
evidence-grade for ℂ (good-prime argument), not a ℚ̄-proof. Combined with
c98's ℚ-boxes, boxeph's weight-3 emptiness, and Smith's d=2 exclusion:
**the ORDER-{1,3} CONJECTURE** — a z-affine Keller map of ℂ³ has field degree
1 or 3 — is now the sharpest supported structural dichotomy; degree-4 seeds
(A₄/S₄ by Smith) require z-degree ≥ 2 ("the 2-jet architecture").

## Files
- `04-computation/jacobian_trisection_modulus_kps_S128c99.py` + `.out`
- `04-computation/jacobian_trace_z_fit_kps_S128c99.py` + `.out`
- `04-computation/jacobian_k4_refined_hunt_kps_S128c99.py` + `.out`
