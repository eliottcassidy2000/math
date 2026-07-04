# The tight locus is the arithmetic progression — and the 14-grid repels covering

*kind-pasteur-2026-07-03-S37. Asked to prove the tight-locus rigidity, I could not close it —
it is the LRC(14) extremal-uniqueness conjecture, at least as hard as the bound. But I pinned
down exactly what rigidity has to say, confirmed there is no exotic family hiding in it, and
proved the mechanism that connects it to covering. The crux is now a single, clean
classification with a rigorous seed.*

## What rigidity has to classify

The tight crux reduces (S36) to: *primitive covering ⟹ M > 1/14*, where `M = max_t min_i
‖v_i t‖`. The clean route is rigidity: characterize the **tight locus** — the families with
`M = 1/14` exactly — and show covering families avoid it.

**Confirmed (computational).** Over all APs `{a, a+d, …, a+12d}`, all dilates `c·{1,…,13}`, and
thousands of random 13-speed families to magnitude 30, **every family with `M = 1/14` is a
dilated AP `c·{1,…,13}`.** The list is exactly `{c·{1,…,13}}`, `c = 1, 2, 3, …`, and nothing
else. In particular:

- The **unique primitive** tight family is `{1, 2, …, 13}` (all dilates have `gcd = c > 1`).
- `{1,…,13}` attains `M = 1/14` **on the 14-grid**: at `t = k/14` (`k` coprime to 14), and
  only there.
- `{1,…,13}` is **non-covering** (it misses `q = 14`), hence sieve-closed at `t = 1/14`.

There is **no "GW" family** — no second, exotic tight family for `n = 13`. The tight locus is
just the arithmetic progression and its dilates. That is worth knowing: it means rigidity here
is a *single*-family statement, not a classification of several exceptional shapes.

## The mechanism: covering families are repelled from the 14-grid

Here is the rigorous seed connecting covering to the tight locus.

**Proposition (14-grid repulsion).** Every covering family has a runner divisible by 14
(`q = 14 ∈ {2,…,14}` forces `14 ∣ v_a` for some `a`). Write `v_a = 14c`. Then:

1. For every integer `k`, `‖v_a · (k/14)‖ = ‖ck‖ = 0` — on the 14-grid, runner `a` sits exactly
   on the observer. So every point `k/14` is unsafe.
2. The danger set of `a` is `D_a = {t : ‖14c·t‖ < 1/14}`, a union of `14c` arcs centred at
   `j/(14c)` of half-width `1/(196c)`; the 14-grid is the sub-lattice `k/14 = (ck)/(14c)` of
   those centres, so an **open neighbourhood of the entire 14-grid is unsafe** (`‖v_a t‖ < 1/14`
   for `|t − k/14| < 1/(14 v_a)`).

Hence the safe set of a covering family omits a neighbourhood of the 14-grid, and **the
optimizer of `M` lies off the 14-grid**.

## The proof skeleton, and the one gap

Put the two together:

- The unique primitive tight family `{1,…,13}` optimizes **on** the 14-grid.
- A primitive covering family optimizes **off** the 14-grid (repulsion).

So a primitive covering family is **not** `{1,…,13}` — and, granting the rigidity classification
(`M = 1/14 ⟺ dilated AP`), it is not any dilated AP either (the only primitive one is
`{1,…,13}`). Therefore `M ≠ 1/14`, and — with the LRC lower bound `M ≥ 1/14` — `M > 1/14`, i.e.
`μ > 0`, i.e. LRC(14) for that family.

The gap is exactly the **rigidity implication** `M = 1/14 ⟹ dilated AP`. Everything else is
rigorous or elementary: the repulsion is proved, the AP is manifestly non-covering, and the
classification is confirmed with no exceptions. But the implication itself is the extremal
uniqueness for the lonely-runner bound at `n = 13`, and it is not lighter than the bound: to
know that a family off the AP locus has `M` strictly above `1/14` is to know the bound is not
merely met but *strictly* met away from a measure-zero set of configurations. That is the same
wall the field hits at `n = 14`.

## Why the AP is the extremizer, geometrically

The AP `{1,…,13}` at `t = 1/14` puts the 13 runners on the 13 nonzero 14th-roots — the unique
equally-spaced configuration, the "billiard" extremal. Rigidity is the statement that this
equally-spaced configuration is *isolated*: perturb the speeds off the AP and the best you can
do strictly exceeds `1/14`, because you can no longer place all 13 runners exactly on the roots
at one `t`. The 14-grid repulsion is the covering-flavoured half of this: a runner divisible by
`14` cannot sit on a nonzero root (`k/14`) at a 14-grid time — it lands on the observer instead
— so covering families are structurally barred from the equally-spaced configuration and pay
for it in `M`. The deep well `{1,…,12,182}` is where a covering family comes closest, and the
price it pays, `14/183`, is set by the Eisenstein resonance (`14` a primitive 6th root mod
`183 = Φ₆(14)`, opus HYP-4047) — the nearest a covering family can approach the AP locus without
reaching it.

## mac-mini's THM-610 rigorized the mechanism (concurrent)

The owner dispatched this crux to several machines at once, and mac-mini's S30 **THM-610**
proved the mechanism I describe here, in sharper and more general form — I record it so the
credit is right and the pieces compose:

- **Lemma 1 (covering ⟹ `q* ≥ n+1`).** A runner divisible by `q ≤ n` sits on the observer at
  every `t = a/q`, so all shallow (`q ≤ n`) hiding spots are dead — hiding is *deep*
  (`q* ≥ n+1 = 15`). This is the 14-grid repulsion, generalized from `q = 14` to all `q ≤ n`
  and stated as the deep-hiding modulus bound; it is the elementary dual of THM-523.
- **Lemma 2 (tight ⟹ `n ∣ q*`).** The minimal residue-distance `q*/n` is an integer, so every
  tight family — *any* branch — has `14 ∣ q*` and its runners sit on a `(q*/14)`-dilated
  14th-root configuration. So "tight ⟹ 14th-root config" is a **theorem**, not the principal
  branch I could only assert. `q* = 14` is the AP; `q* = 28` is the even block.
- **Corollary: tight covering ⟹ `q* ≥ 2n = 28`**, and the margin map gives primitive
  covering-min `M/(1/n) ∈ [1.06, 1.11]` uniformly for `n = 7..14` — the looseness, bounded away
  from 1.

So the 14-grid repulsion and "tight ⟹ 14th-root config" are now rigorous (THM-610). What this
note adds on top is the **classification by primitivity**: among *primitive* families the tight
locus is the single family `{1,…,13}` (`q* = 14`); the `q* = 28` even block is `2·{1,…,13}`,
imprimitive, reducing to `{1,…,13}`. No third, exotic ("GW") tight family appears to magnitude
30. So rigidity is a one-family statement.

## Honest placement

Not a proof. With THM-610 the mechanism (repulsion, `tight ⟹ 14th-root config`) is rigorous;
the classification (`primitive tight = {1,…,13}` uniquely, no GW) is confirmed to magnitude 30;
and the residual is the one implication `M = 1/14 ⟹ dilated AP` — LRC(14) extremal-uniqueness
at `n = 13`, not lighter than the bound. The crux is now a single, precisely-stated
extremal-uniqueness question with a rigorous mechanism underneath it, not a vague "force the
measure positive."

---
*Linked: [[the-covering-min-and-the-gcd-refinement]] (S36, the value 14/183 + gcd), HYP-4047/opus
(Eisenstein extremizer), HYP-3551 (covering-min), HYP-4058/opus (measure form). Scripts:
`lrc14_tight_locus_enum_kps_S37.py`. HYP-4062.*
