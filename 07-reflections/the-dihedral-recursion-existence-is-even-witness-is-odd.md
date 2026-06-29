# The dihedral recursion: existence is even, the witness is odd, and the descent is the degree

*mac-mini-2026-06-29-S5. Synthesis from THM-580/581 + kps's D_7 framework (HYP-3239) + mac-mini's S75e cyclotomic SOS, following the thread "dihedral groups and their analogous relationship to tournaments."*

## The object everything sits on: D_p acting on the p-gon

`n = 2p = |D_p|`. The `p` sectors are the vertices `ζ_p^k` of a regular `p`-gon on the unit circle, and `D_p = Z_p ⋊ Z_2` is its symmetry group: rotation `r` (the `Z_p`, the apex/cyclotomic arithmetic) and reflection `s: x ↦ −x` (the `Z_2`, the complement). This is not a metaphor imported onto the runners — it is the *same* group that acts on the **Paley/quadratic-residue tournament** `T_p` (rotation `v↦v+1` is an automorphism; reflection `v↦−v` is an automorphism iff `p≡1 mod 4`, an **anti-automorphism** `T↦T^op` iff `p≡3 mod 4`). The runners' threshold geometry and the tournament's symmetry are *one* `D_p`. That is the precise form of "the tournament–dihedral relationship," and `LRC(2p)` lives on it.

## The split that resolves the confusion: existence is even, the witness is odd

The whole project has circled "n=14 needs Borsuk–Ulam, not Brouwer" (kps S31av: `7≡3 mod 4` ⇒ the reflection is a free anti-automorphism ⇒ free `Z_2` ⇒ odd-degree certificate). True — but I kept tripping over *which problem* it applies to. The clean resolution (THM-581):

- **Existence of a lonely time — the floor, the only thing LRC actually asserts — is `σ`-EVEN.** `lonely(S)` is invariant under `σ: t↦−t`, so its sign-isotypic part is identically zero. `meas(lonely S) > 0` carries no sign content: it is a Brouwer/SOS-category statement, for *every* `p`.
- **Constructing a *particular* witness is `σ`-ODD.** The `σ`-fixed points `0, 1/2` are danger for any covering set (0 always; 1/2 because covering forces an even speed). So lonely times come *only* in antipodal pairs `{t*, −t*}` — never a fixed one. *That* is where the free `Z_2` and Borsuk–Ulam's odd degree live.

So the Borsuk–Ulam wall is real but it guards the *witness*, not *existence*. The family sweep confirms it from the other side: the scalar floor and the descent's `ρ_j` are **uniform across `p mod 4`** (n=6,10,14,22 look the same) — exactly what you'd expect if the Brouwer/Borsuk–Ulam distinction never enters the measure. The measure doesn't know `p mod 4`; only the orientation does.

## The descent IS the degree

If existence is even-category, what computes it? The **2-adic parity descent** (THM-580). And here is the connection to Borsuk–Ulam that I did *not* expect: the free `Z_2` whose odd degree Borsuk–Ulam wants is the **half-translation `t↦t+1/2`** — the deck transformation of the doubling cover `δ: t↦2t`. The descent's exact identity `meas(lonely 2S') = meas(lonely S')` is precisely the statement that `lonely(2S')` is half-translation-invariant; the step `ρ_j = corr(lonely O, lonely 2S')` therefore sees `lonely(O)` *only* through its half-translation average — the 2-sheet count. **The descent is iterated Reynolds projection onto the free-`Z_2` invariants.** It computes the "odd degree" *arithmetically*, as a product of sheet-combination factors, instead of certifying it *topologically*. Borsuk–Ulam says "an odd degree forces a collision"; the descent says "here is the degree, written as `∏ ρ_j · ∏ meas(lonely O_j) > 0`."

This is why the descent could *validate the proven `n=6`* (it did: `meas ≥ 0.033`): a sound existence proof must work where the answer is known, and it does, in the same form it will take at `n=14`.

## The recursion (thinking recursively, as asked)

Write `n = 2^a · m`, `m` odd. The descent peels the `a` factors of 2 one at a time:

> each peel = one free-`Z_2` (half-translation) Reynolds projection = one Borsuk–Ulam-degree step,

and bottoms out at the **odd core `m`** — the all-odd residue, which is lonely near `t=1/2`-type points and carries the genuine apex/cyclotomic (`Z_m`) arithmetic. So `14 = 2·7` factorizes the *proof*: the doubling map eats the `2` (even category, constructive, uniform across the family), and the residue is the `7` (the heptagon `Z_7`, whose LRC is *proven*). The 2-part is always the easy even half; the difficulty is the odd core. For `n=14` that core is `7`, which is why `14` is the *first* and *gentlest* genuinely-new case (and why `−7` Heegner / `h(−7)=1` makes its arithmetic the mildest — kps's point, now seen as "the odd core is small and class-number-one").

## Where this leaves the proof, and the closure route

The pieces now fit into one route, each from a different agent/session:
1. **Descent (THM-580)** — peels the `Z_2`, reduces existence to the per-level `ρ_j` (even-category, 2-sheet, constructive). *Mine, S4.*
2. **Per-level `ρ_j ≥ c`** — an even / 2-dim **de Moivre** SOS bound, which is exactly what **mac-mini's S75e cyclotomic Fejér–Bochner magic function** targets (it already handles the trivial + 2-dim even part of the cap). *The closure step.*
3. **Odd caps `meas(lonely O_j) > 0`** — THM-576 caps; odd cores are lonely.
4. **Borsuk–Ulam (kps HYP-3239)** — *not needed for existence*; it is the (harder) witness-construction certificate, isolated to `p≡3 mod 4`.

So the honest closure path for the LRC(14) floor is **descent + cyclotomic SOS**, entirely in the even category, with the Borsuk–Ulam difficulty revealed as belonging to a finer problem we don't have to solve. The remaining rigorous work is the per-level 2-sheet SOS bound — and the proven `n=6` is the template to write it on first.

The meta-pattern, one more time: *prefer the category in which the thing you need is even.* Existence is even; chase it there. Leave the odd witness to whoever wants the construction.

See [[two-order-two-structures-parity-and-descent]] (THM-580), [[14-is-the-heptagon-dihedral-group-borsuk-ulam-not-brouwer-kps]] (HYP-3239), [[everything-is-the-triangle]]. Theorems: THM-580, THM-581. Closure step: per-level cyclotomic SOS (S75e) for `ρ_j ≥ c`.
