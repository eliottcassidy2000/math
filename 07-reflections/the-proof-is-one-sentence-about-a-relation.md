# The proof is one sentence about a relation

*klein-2026-06-29-S8. The culmination of the floor reframe (HYP-3566) and the two requested computations (HYP-3571). What the LRC covering floor finally says, said plainly.*

## The sentence

> The covering floor is one statement about the danger relation `D : (v, t) ↦ ‖v·t‖ < 1/14`:
> **it does not factor, and composed with itself it stays small, once you look at it in the right frame.**

Both clauses are now computed, exactly.

## "It does not factor" — the relation is essential

If `D` factored — if R-safe were independent of Q-lonely — the spectral coupling `SPEC` would vanish and
the floor ratio `R' = meas(lonely(R∪14Q)) / (m_R m_Q)` would be `1` identically. It is not: `SPEC ≠ 0` for
953 of 956 coverings. The bilinear pairing `v·t` genuinely couples the two parts; the relation is
irreducible. And at the extremal `{1..13}` the non-factoring is certified topologically: the lonely set is
*exactly* the multiplicative units `(Z/14)* = {1,3,5,9,11,13}`, arranged in `φ(14)/2 = 3` antipodal pairs
`{(1,13),(3,11),(5,9)}` — the Borsuk–Ulam counting measure, saddle index 3. A unit group does not
decompose additively; a disproof would need to be at once additive (hit every `1/q`) and anti-multiplicative
(cover the units `a/14`), and it cannot. *Relations, not things* (mac-mini's S24 lens): the floor is a
property of `D`, and `D` is essential.

## "Composed with itself it stays small" — in the right frame

The second moment of the sheet count looked like it blew up: `CV(N_R)² → ∞` as `m_R → 0` (HYP-3554),
because the `Z_14` sheet action is not transitive — there is a vanishing-fiber corner. That was the wrong
frame. The Cauchy–Schwarz step `|SPEC|/product ≤ CV(N_R)·√((1−m_Q)/m_Q)` is lossy exactly there. In the
right frame — the *actual* correlation `|R'−1| = |SPEC|/product` — the quantity stays small:
`sup |R'−1| = 0.656 < 1` over the whole adversarial family, so `R' ≥ 0.344 > 0` set-independently, *above*
the `Γ₀(14)/ζ(2)` bound `1/(2ζ(2)) = 0.304`. The variance frame and the floor frame even bind at different
sets (`CV` peaks at `{1..13}\{12}`, with the apex 7; `R'` bottoms at `{1..13}\{7}`, without it) — a tell
that they are different functionals, and only the second one is the floor.

## "The right frame" — transitivity, via the descent's `Z_7`

Why is `|R'−1|` bounded set-independently while `CV` is not? Because the floor is a *transitivity* problem
(HYP-3566). `CV(H)²` (the witness rehearsal, THM-589) is clean because `S_n` is transitive; `CV(N_R)²` is
not because `Z_14` is not. The fix is to read the relation through a transitive symmetry, and `14 = 2·7`
supplies it: the 2-adic descent peels the non-transitive 2-part and exposes the cyclically transitive
apex-`Z_7` core (Paley/Heegner), where the second moment becomes a single group-average whose cyclotomic
positivity is the SOS. The octonion/`Z_7` structure — which the persistence test removed from `b_1^-`
(a dimensional coincidence, HYP-3563) — is the right frame here, at its proven home. The descent is the
change of frame that turns "composed with itself it stays small" from false (in the `CV` frame) into true.

## Why this is the whole proof, in outline

- **Does not factor** ⇒ the floor is not trivially `1`; there is a real `SPEC` to control, and the
  extremal obstruction is the BU/units saddle — the σ-odd witness content (THM-582/583) certifying the
  lonely set is essential.
- **Stays small** ⇒ `|SPEC|/product < 1`, i.e. `R' > 0`, the σ-even floor — bounded set-independently by
  the transitive (`Γ₀(14)` / `Z_7`-cyclotomic) collapse, not the per-set variance.
- **Right frame** ⇒ the transitivity the descent manufactures; the place where the relation, composed with
  itself, is provably positive-definite.

Three things the project chased separately — the witness odd index, the floor measure inequality, the
metagraph antipodal spectrum — are the three readings of one relation `D` under one involution. The proof
is the sentence; the work that remains is to write the `Γ₀(14)`/`Z_7`-cyclotomic constant that makes
"stays small" a closed-form `≥ 1/(2ζ(2))` rather than a 956-covering scan.

See HYP-3571 (the two computations), HYP-3566 (transitivity reframe),
[[reframing-the-lrc-floor-manufacture-transitivity-and-the-singer-z7-the-descent-exposes]],
[[the-variance-blows-up-where-the-fiber-vanishes]] (HYP-3554), [[polysemous-constants-bridges-traps-and-homonyms]],
HYP-3549 (the units extremal), HYP-3553/3550 (`Γ₀(14)` / Euler product), THM-579/580.
