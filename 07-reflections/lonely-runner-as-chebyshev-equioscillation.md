# The Lonely Runner as a Chebyshev equioscillation problem on (ℤ/14)*

*kind-pasteur-2026-06-27-S255. The owner noticed the AP's safety function just touches 1/14 at exactly
six points, equioscillating, in three antipodal pairs mirror-symmetric about t=1/2. That observation
is the whole frame: LRC(14) is a Chebyshev min-max extremal problem, and the touch-points are the unit
group.*

## The observation, identified

For a 13-set S the **safety function** is `f_S(t) = min_{s∈S} ‖s·t‖`, and `M(S) = max_t f_S(t)` is the
lonely-runner gap. LRC(14) says `M(S) ≥ 1/14` for all primitive S, with the AP `{1,…,13}` achieving
the *minimum* `M = 1/14` — the AP is the tightest set.

The minimizer of a max-functional **equioscillates** — this is the Chebyshev / Kolmogorov signature of
an extremal. And the owner saw it directly: `f_AP` touches its max 1/14 at exactly six points. Those
points are not arbitrary — they are exactly the **units of ℤ/14**:

> `f_AP` reaches 1/14 at `t = a/14` for `a ∈ {1,3,5,9,11,13} = (ℤ/14)*`, `φ(14) = 6` of them, in three
> antipodal pairs `(1,13), (3,11), (5,9)` mirror-symmetric about `t = 1/2`.

And each unit is pinned by a **binding complement-pair**: at `t = a/14` the runner `s` with `s·a ≡ ±1`
(mod 14) sits at `±1/14`, and these are the three pairs `{1,13}, {3,11}, {5,9}` — each summing to 14,
i.e. *the apex-7 binding pairs* (HYP-2909). GW equioscillates at the same six units. So "equioscillation
at the unit group" is the common signature of the tight locus, and three independently-found conditions
collapse into it:
- the **±units-cover** (mac-mini: every unit `a` has a runner with `s·a ≡ ±1`),
- the **apex-7 binding pairs** (the antipodal pair at `±1/14`),
- the **complement symmetry** `v ↔ N−v` (the mirror about `t = 1/2`).

They are one fact — the safety function equioscillates at `(ℤ/14)*`, with antipodal (complement)
structure.

## Why this is a frame, not a coincidence: the case split falls out

The equioscillation gives an elementary identity that organizes the entire proof. At a unit `a`, the
orbit `{s·a/14}` lies on the (1/14)-grid, so the nearest runner to 0 is at distance `0` or `≥ 1/14`:

> `f_S(a/14) < 1/14 ⟺ some `s ≡ 0 (mod 14)`.

Three regimes follow immediately:

1. **14-free S** (no multiple of 14). Then `f_S(a/14) ≥ 1/14` at *every* unit, so `M(S) ≥ 1/14`.
   **This is THM-523 (q=14) — recovered as the equioscillation of `f` at the six units.** The
   non-covering case *is* the equioscillation case, and it is proved for free.
2. **Tight (M = 1/14).** The six units are the *global* max: `f` equioscillates at 1/14 with no higher
   peak anywhere. This is the census `{AP, GW}` — equivalently, on the spectrum frame, the deepest
   sink stays pinned to the apex Farey node `1/14` and no Farey child out-sinks it.
3. **Covering** (a multiple of 14 present). The unit points are *killed* (`f < 1/14` there — the
   multiple of 14 sits on the observer), so `M` is achieved *off* the units, in the multi-far
   decorrelation floor (the `R ∪ 14Q` lift; HYP-3132).

So the Chebyshev frame is the **local extremal half** — equioscillation at the unit group hands you the
14-free case (THM-523) and pins the tight value 1/14 — and the **global half** (no higher peak / the
covering floor) is exactly the residual the spectrum and decorrelation frames isolate. The same tight
locus, read three ways: equioscillation (local), binding-scale on the Farey tree (global), divisibility
threshold `q(S) = 14` (arithmetic). The frames agree because they are shadows of one object: a
min-max extremal whose active set is the unit group of the apex ring.

## What it buys toward a proof

It re-derives the proved part (14-free, THM-523) as the equioscillation, and it tells you precisely
what "tight" means in extremal-theory language: the six unit maxima are global, i.e. the three binding
complement-pairs leave **no descent direction** (Kolmogorov). The census conjecture becomes "no 14-free
set other than AP/GW makes the unit equioscillation global" — a finite-dimensional extremal-rigidity
statement on `(ℤ/14)*`, mirror-symmetric, with three active pairs. That is a far more structured object
than "min over infinitely many integer sets."

## The pointer beyond

Chebyshev equioscillation at `n` points is the signature of the extremal in an `(n−1)`-dimensional Haar
system. Here `n = φ(14) = 6`, so the AP is the extremal of a *5-dimensional* effective system — and the
mirror symmetry halves it to **3 free pairs**, exactly the three binding complement-pairs. The
conjecture to chase: the right finite-dimensional Chebyshev/Haar system in which `f_S` lives has
dimension `φ(N)/2`, the extremal is unique up to the complement involution, and the AP/GW are its two
representatives. If the LRC tight locus is the equioscillation set of a `φ(N)/2`-dimensional Chebyshev
problem on `(ℤ/N)*`, then the whole conjecture is an instance of the **uniqueness of the Chebyshev
extremal** — and `φ(N)/2` (here 3) is why the census is so small. The runners equioscillate on the unit
group; the tight ones are the Chebyshev points of the apex ring.

— Related: [[lrc14-thread]], HYP-3246/3247 (this frame), HYP-2909 (binding pairs), HYP-2928 (the spectrum
/ binding scale), HYP-3132 (the multi-far floor), THM-523 (q-witness), THM-530 (caps);
`the-tournament-spectrum-is-the-object.md`, `lrc14-three-engines-lift-and-decorrelate.md`.
