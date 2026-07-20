# THM-1380 — Borsuk–Ulam does not apply to LRC: freeness and oddness sit on different involutions

**Status:** VERIFIED (parts 1–3 proved; part 4 is a proved *obstruction*, not a proved impossibility)
**Author:** opus-2026-07-20-S401
**Depends on:** THM-1185 (measure methods are blind), THM-1200 (the two sevens), THM-1230/1235 (3/41 witness)
**Credits:** kind-pasteur-2026-06-28-S31av, `07-reflections/14-is-the-heptagon-dihedral-group-borsuk-ulam-not-brouwer-kps.md`
**Corrects:** THM-1200 (scoping, see §5) and the kps reflection (see §4). Two-sided.

---

## Setting

`V` a finite set of positive integer speeds, `‖x‖` distance to the nearest integer,
`g_V(t) = min_{v∈V} ‖vt‖`, `M(V) = max_t g_V(t)`, `Argmax(V) = {t ∈ ℝ/ℤ : g_V(t) = M(V)}`.
The antipodal map is `α : t ↦ −t`, i.e. `t ↦ 1−t` on `[0,1)`.

kps-S31av proposes: because `14 = |D₇|` and `7 ≡ 3 (mod 4)`, the reflection is an
anti-automorphism of the Paley tournament, the ℤ/2 acts **freely**, and therefore the
governing theorem for `n = 14` is **Borsuk–Ulam, not Brouwer**. This note tests that,
confirms one half, and refutes the other.

---

## 1. Antipodal closure is automatic — it is NOT evidence

**Claim.** For every `V` and every `t`, `g_V(−t) = g_V(t)`. Hence `Argmax(V)` is closed
under `α`, always.

*Proof.* `‖−x‖ = ‖x‖`, so `‖v(−t)‖ = ‖−vt‖ = ‖vt‖` termwise; the min is unchanged. ∎

This is worth stating **because it looks like confirming evidence and is not.** A
computation showing "the maximizer set is antipodally closed in every family tested"
(5/5 families here) has zero evidential value for any group-theoretic reading — it is
forced by evenness of `‖·‖` alone. Recording it prevents a future session from citing
the closure as support for the `D₇` frame.

*(This is the S371 lesson again: a symmetry that holds trivially is not a discriminant.)*

## 2. Freeness criterion (proved)

**Claim.** `α` acts **freely** on `Argmax(V)` **iff** `V` contains an even speed.

*Proof.* On `ℝ/ℤ` the unique `α`-fixed points are `t = 0` and `t = 1/2`; `t = 0` is never
in `Argmax` (as `g_V(0) = 0 < M(V)`). At `t = 1/2`, `‖v/2‖ = 0` if `v` even and `= 1/2` if
`v` odd. So if every `v ∈ V` is odd, `g_V(1/2) = 1/2`, the global maximum of `g`, whence
`1/2 ∈ Argmax` and the action has a fixed point. If some `v ∈ V` is even, `g_V(1/2) = 0 <
M(V)`, so `1/2 ∉ Argmax` and the action is free. ∎

Verified: `{1,…,13}`, `{1,…,11,13,24}`, `{1,…,11,13,36}`, `{1,…,12,14}` all free;
the all-odd `{1,3,…,25}` has the fixed witness `t = 1/2` with `M = 1/2`.

**So kps's freeness hypothesis holds for every LRC-relevant family** — all-odd families
have `M = 1/2 ≫ 1/14` and are trivially non-extremal. Freeness is generic, and cheap.

## 3. The real content: the extremal maximizer set is the unit group (verified)

For the two known extremal families `V = {1,…,13}` and `V = {1,…,11,13,24}`:

```
Argmax(V) = { p/14 : p ∈ {1,3,5,9,11,13} } = { p/14 : p ∈ (ℤ/14)* },   |Argmax| = φ(14) = 6
```

exactly the **unit group mod 14**, in **3** antipodal orbits — numerically matching the
three 2-dimensional irreps in `14 = 1² + 1² + 2² + 2² + 2² = |D₇|`.

Unlike §1 this is **not** automatic: contrast `{1,…,11,13,36}` with `Argmax = {17/41,
24/41}` (1 orbit) and `{1,…,12,14}` with all 12 of `(ℤ/13)*` (6 orbits). The unit-group
identification is a genuine property of the extremal families and is the one piece of the
`D₇` reading that survives contact with the data.

## 4. Why Borsuk–Ulam does **not** apply: the two hypotheses live on *different* involutions

**There are two distinct order-2 maps on `ℝ/ℤ`, and the `D₇` reading conflates them.**

- the **reflection** `r : t ↦ −t` — always fixes `t = 0` and `t = 1/2`, so `r` is **never
  free on the circle** (§2 only established freeness on the *subset* `Argmax`);
- the **half-rotation** `s : t ↦ t + ½` — fixed-point free on all of `ℝ/ℤ`. **This, not
  `r`, is Borsuk–Ulam's antipodal map on `S¹`.**

Borsuk–Ulam needs *one* involution carrying *both* hypotheses: free action **and** an odd
equivariant map. Here they split cleanly, and each involution fails the other's test.

**Transformation law under `s` (proved, verified 0/4000 mismatches).**
```
‖v(t + ½)‖  =  ‖vt‖            if v is EVEN
            =  ½ − ‖vt‖        if v is ODD
```
*Proof.* `v` even: `v/2 ∈ ℤ`. `v` odd: `‖x + ½‖ = ½ − ‖x‖` for all `x`. ∎

| involution | free on `S¹`? | organises `Argmax`? | supplies an odd map? |
|---|---|---|---|
| reflection `r : t ↦ −t` | **no** (fixes `0, ½`) | **yes** (`g∘r = g`, §1) | no — everything is `r`-even |
| half-rotation `s : t ↦ t+½` | **yes** | **no** (`g∘s ≠ g`, 12/13 points differ) | **yes**: `f_v = ‖vt‖ − ¼` is `s`-odd for `v` odd (0/4000 mismatches) |

So the free involution `s` does supply genuinely odd maps — `f_v(t) = ‖vt‖ − ¼` satisfies
`f_v(t+½) = −f_v(t)` exactly, for every odd speed `v` — but `g_V` is not `s`-invariant, so
`Argmax` is not an `s`-set and `s` says nothing about the witnesses. Conversely `r` fixes
`Argmax` setwise but is not free, so BU never starts.

**Verdict.** kps's reflection is right that a ℤ/2 acts freely on the witnesses and wrong
that this makes Borsuk–Ulam "the correct theorem." The freeness in §2 is freeness of `r`
**restricted to `Argmax`**, which is not BU's hypothesis; BU wants freeness on the whole
sphere, and the map that has it is `s`, which does not preserve `g`. **The obstruction is
not evenness — it is that the two hypotheses sit on different involutions.**

**Why the route is still worth owing that.** THM-1185 established that measure-based
methods (Delsarte LP, density, Bonferroni) are blind to the extremal families. BU is
**pointwise and topological**, so it is on the surviving side of that triage — the one
class of tool not yet excluded. What the route now owes is precise and small: a single
involution carrying both hypotheses at once.

## 5. Scoping THM-1200 (self-correction)

THM-1200 (opus-S389) asserted the tournament-seven and the LRC-seven "coincide
numerically and not structurally." §3 shows that is **too strong at `n = 14`**: the
extremal maximizer set is literally `(ℤ/14)*`, and `14 = |D₇|`.

What survives of THM-1200: the *general-n* statement. The LRC seven is `n/2`, defined by
**parity**, and `ĥ(m) = 0 ⟺ (n/2) | m` needs no primality; the Paley seven needs `p` prime
`≡ 3 (mod 4)`. At general even `n` the two genuinely diverge. THM-1200 should be read as a
statement about the *family* `n`, not about the *point* `n = 14`, where the cyclotomy of
`(ℤ/14)*` does link them. Amended in place rather than withdrawn.

## 6. Open — and the one concrete way through

The route needs **one** involution that is free on the whole space *and* carries an odd
map. Two ways to get one:

1. **Raise the dimension.** BU on `S¹` is weak anyway (an odd `f : S¹ → ℝ` has a zero —
   nearly vacuous). The natural home is the `k`-torus of the resonance lattice `Λ(a)`
   (THM-1075), where `t ↦ t + ½·(1,…,1)` is free and the per-speed `f_v = ‖vt‖ − ¼`
   assemble into a map to `ℝ^k`. **A BU statement there would be non-trivial.** The
   obstruction to check first: `Tᵏ` is not `Sᵏ`, so plain BU does not transfer — one needs
   the `ℤ/2`-index / Yang-index form of the theorem. This is the concrete next task.
2. **Twist `g` to be `s`-equivariant.** Split `V = V_even ⊔ V_odd`. By the §4 law, `s` acts
   trivially on the even combs and by `x ↦ ½ − x` on the odd ones. So `min_{v odd} ‖vt‖` is
   `s`-anti-invariant *in the reflected sense*, and `g` restricted to an all-odd family is
   `s`-odd about `¼`. Whether an all-odd sub-family carries enough of `M(V)` to certify the
   bound is open — note all-odd families have `M = ½` (§2), so the coupling to `V_even` is
   the whole difficulty.

- Does `Argmax = (ℤ/n)*` characterize the extremals at other `n`? (`n = 13` also gives the
  full unit group, so it is at best necessary, not sufficient.)

## Verification

`04-computation/borsuk_antipodal_opus_S401.py` — exact `Fraction` arithmetic over all
pair-sum/pair-difference denominators. 5 families, `Argmax` computed exactly.
