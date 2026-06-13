---
source: opus-2026-06-06-S699 (monodromy + Rado-as-tournament)
status: SYNTHESIS — decodes the two seeds. "A loop of the input swaps the two outputs" = the MONODROMY of the witness map (coeff→root, FTA S699l): looping inputs around the DISCRIMINANT (where two roots/witnesses collide = the worry-set, S699l) swaps them — a Galois transposition / braid half-twist (verified z²−t). The LRC runners form a SPACETIME BRAID (S540o); crossings = pair-clock zeros (S699) = the swaps. "Rado graph as a tournament" = the countable RANDOM tournament (Fraïssé limit), the GENERIC fiber: the extension property = "always a witness" = LRC's loneliness genericity (verified P→1). The Galois/monodromy GRADES the fibers: cyclotomic worry-set = ABELIAN monodromy (the floor); Rado generic = full SYMMETRIC monodromy.
tags: [monodromy, galois, braid, discriminant, worry-set, branch-locus, rado-graph, random-tournament, fraisse, extension-property, FTA, cyclotomic, generic-fiber, LRC, spacetime-braid]
---

# Monodromy, the Rado tournament, and the worry-set as branch locus

**Prompt (user):** "a loop of the input causes a swap of the two outputs"; "consider the Rado graph
as a tournament."

Two seeds, one picture. They decode through the FTA-duality (S699l) into a single statement:
**witness-finding is a covering map; the worry-set is its branch locus; the Rado tournament is its
generic fiber.**

## 1. "A loop of the input swaps the two outputs" = monodromy

The FTA map sends **inputs** (coefficients) to **outputs** (roots) — and it is a *branched cover*.
Looping the inputs around the **discriminant** (where two roots collide) **swaps those two roots**:
a Galois transposition, a braid half-twist. **Verified** (`…s699p.py`): `z²−t`, loop `t=e^{iθ}`,
`θ:0→2π`; the root tracked from `+1` returns at `−1` — *the two roots swapped*. The branch point
`t=0` is the discriminant (roots collide).

> **The discriminant is the worry-set.** S699l/HYP-2275 showed the LRC worry-set = the locus where
> the covering-depth polynomial's roots **collide** (multiplicities / the all-orders cancellation).
> So "a loop of the input swaps two outputs" *around the worry-set* is the monodromy that **detects**
> it: the worry-set is exactly the branch locus where two witnesses merge, and a loop around it
> transposes them. The monodromy group is the **Galois group of the witness map.**

## 2. The LRC runners are a spacetime braid; crossings are the swaps

The LRC runners `v_i t mod 1` over `t∈[0,1]` are **strands** of a braid on the cylinder (the
spacetime braid, S540o). Two runners **cross** when `(v_i−v_j)t ∈ ℤ` — and that is *exactly* the
**pair-clock zero** of the signed LRC (S699). **Verified:** `V=[1,2,3]` crosses at
`t=½ (1,3)` and `t=1 (all)`; the braid length / writhe `= Σ|v_i−v_j| = 4`. **Each crossing is a
transposition** (a braid generator); the worry-set is where crossings **collide** (the
discriminant / branch). So the signed-LRC pair-clocks (S699) ARE the braid generators, and the
"loop swaps two outputs" is a runner crossing.

## 3. The Rado graph as a tournament = the generic fiber

The "Rado graph as a tournament" is the **countable random tournament** `T_∞` — the Fraïssé limit
of finite tournaments, universal and homogeneous. Its defining **extension property**: for any
finite disjoint `A,B`, there is a vertex **beating all of `A`, beaten by all of `B`** — *a witness
with any prescribed relation.* **Verified:** for random finite tournaments,
`P(∃ such witness) = 0.34, 0.66, 0.93, 0.99` for `n=10,20,40,80` → `1`.

> **The extension property IS the LRC's "always a witness."** LRC asserts a lonely time always
> exists (a witness avoiding all danger arcs); the Rado tournament asserts a witness with any
> prescribed beating-pattern always exists. The Rado tournament is the **generic fiber** where
> genericity (witness existence) holds maximally; the **worry-set is the special fiber** where it
> *degenerates* to a measure-zero, forced (collided) witness — the branch point.

## 4. The Galois grades the fibers: cyclotomic (abelian) floor vs Rado (symmetric) generic

The monodromy/Galois group of the witness map distinguishes the fibers:
- **Worry-set fiber:** the witnesses are **roots of unity** (THM-403, the cyclotomic clock), so the
  Galois is **abelian (cyclotomic)** — *low monodromy*, the special branch. (This is the same
  cyclotomic floor as S699o: the torsion unit-norm elements.)
- **Generic (Rado) fiber:** the witnesses are unconstrained, so the monodromy is the **full
  symmetric group** — *maximal monodromy*, matching the Rado tournament's huge homogeneous
  automorphism group.

> **The Galois/monodromy group GRADES the problem:** abelian (cyclotomic) at the worry-set floor,
> symmetric at the Rado generic. The "escape" (S699o: non-torsion unit-norm, the field tower) is
> the monodromy **growing** from abelian toward symmetric — and the Rado tournament is the
> `S_∞` limit. So HN's field-tower (S699o), the cyclotomic worry-set (THM-403), and the Rado
> generic are one ladder: *the Galois group of the witness map, from abelian floor to symmetric
> generic.*

## 5. The unified picture

> **Witness-finding (LRC lonely time / Hamiltonian path / polynomial root) is a branched covering
> map.** Its **monodromy** ("loop input ⟹ swap two outputs") is the Galois/braid group, generated
> by transpositions at the **discriminant = the worry-set** (where two witnesses collide). The
> **Rado/random tournament is the generic fiber** (extension property = always a witness, full
> symmetric monodromy); the **worry-set is the special/branch fiber** (collided witnesses, abelian
> cyclotomic monodromy). The LRC runners realize this as a **spacetime braid** whose crossings are
> the swaps. *Genericity (Rado) ⟷ the worry-set is the branch locus where genericity breaks.*

## 6. Honest status

- **Verified:** the `z²−t` monodromy swap; the LRC braid crossings (= pair-clock zeros, writhe
  `Σ|v_i−v_j|`); the Rado/random-tournament extension property (`P→1`).
- **Established (standard, here mapped):** FTA branched cover & monodromy = Galois; the random
  tournament = Fraïssé limit with the extension property; the discriminant = collided-roots locus.
- **New (the synthesis):** the worry-set = the discriminant/branch locus of the witness cover (the
  monodromy detects it); the LRC pair-clocks (S699) = braid generators; the Rado tournament = the
  generic fiber whose extension property = LRC genericity; the **Galois grading** (cyclotomic-abelian
  worry-set floor ↔ symmetric Rado generic), unifying with the field tower (S699o).
- **Conceptual, not a theorem:** a unifying lens (covering/monodromy/Fraïssé), not a new bound. But
  it suggests a probe: **the worry-set is where the witness-cover ramifies; computing the monodromy
  (Galois) group of the covering-depth PGF over config space would locate the worry-set as its
  branch locus and grade it by the Galois group (abelian ⟹ worry-set, symmetric ⟹ generic).**

**Artifacts:** `04-computation/monodromy_braid_rado_s699p.py` (+`.out`). Builds on S699l (FTA
coeff→root / worry-set = roots-collide), S699 (pair-clocks), S540o (spacetime braids), S699o
(unit-norm / field tower), THM-403 (cyclotomic witnesses), Fraïssé/Rado, Galois/braid monodromy.
New: **HYP-2282**.
