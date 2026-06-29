# The ±1 eigenspace of reversal is the whole split — tested on both sides

*mac-mini-2026-06-29-S8. Pursuing the owner's two tasks and the organizing principle: the even/Brouwer ⊕ odd/Borsuk–Ulam split is literally the ±1 eigenspace split of the reversal `R`.*

## The principle, made operational

The owner's frame turns a pile of dualities — even/odd, Brouwer/Borsuk–Ulam, SOS/odd-degree, symmetric/antisymmetric, real/imaginary — into one linear-algebra statement: **`R` (the reversal/complement involution) has eigenvalues ±1, and every "two-index" phenomenon in this project is that eigenspace split.** The `ε=+1` (R-symmetric) part is the SOS/Brouwer-provable bulk; the `ε=−1` (R-odd) part is the signed obstruction the bulk can't see; and "going to the half" is the projection onto a fundamental domain that must not silently drop the `ε=−1` coordinate. I tested it on both sides and it held, with one clarifying refinement.

## The witness side (Task 1): `f` lives on the half, exactly

The odd index `f` — the palindromic Hamiltonian-path count (THM-582) — *is* a count on the half-system. A `φ`-palindromic path `(v_1,…,v_m, c, φv_m,…,φv_1)` is its own first half plus the fixed center; and the second half imposes **no new constraint**, because the anti-automorphism identity `arc(u,w) ⇔ arc(φw,φu)` makes every mirror arc equivalent to a half arc (THM-583). So `f` is a Hamiltonian-path count on the `(p−1)/2` pairs, computed by a transfer DP over half the vertices (`f = 1, 9, 185` for Paley `p=3,7,11`; DP states `2,21,150` vs `p! = 6,5040,4·10⁷`). The half-compression is *lossless* — provided you keep `φ`. That is the precise content of the warning "don't drop the `ε=−1` part": the `ε=−1` data is not thrown away, it is *stored in `φ`* and regenerated. Lose `φ` and you lose the second half; keep it and the half is the whole.

## The cap side (Task 2): the obstruction is the R-odd eigenspace

On the cap, `R` acts on the six inner sectors as `(1 5)(2 4)` (fixing 3, 6), and kps's co-emptiness matrix `M` commutes with it, so `M = M_even ⊕ M_odd` with `M_odd` two-dimensional on `span(e_1−e_5, e_2−e_4)`. For *every* binding configuration I tested — `consec_8`, `consec_9`, and the actual minimizers `{1,5,7,8,9}` and `{1,11,12,13}` — the **Perron (bulk) mode is R-even** and a **nonzero R-odd eigenspace contributes negatively** to `S₂` (the pairwise co-emptiness, `−tr(M_odd)/2`). The bulk you can bound by a cosine (R-even) Fejér–Bochner SOS is the `ε=+1` part; the residual the cosine cannot reach is the `ε=−1` part. The eigenspace split is not a metaphor here — it is the literal block structure of the matrix everyone has been staring at, and the obstruction is the `ε=−1` block.

## The refinement: the split is *per obligation*

The one correction the test forces is a good one. Last session (THM-581) I argued the LRC floor is purely R-even — "existence needs no Borsuk–Ulam." That is still true, and Task 2 shows *why it doesn't contradict the owner's principle*: the two LRC obligations sit on opposite sides of `R`.

- **The floor / existence** is the lonely *measure* — a function of `t` that is literally R-invariant, so its `ε=−1` component is identically zero. Pure even. (And a lonely tournament has a *source*, which the converse turns into a *sink*: it cannot be self-converse, so it carries no R-odd witness.)
- **The cap / concentration** is the co-emptiness *matrix* — and there `M_odd ≠ 0` is real and is the obstruction.

So the principle is exactly right about the cap, and the floor was already living on the even side. The "signed certificate everyone needs" is the cap's `M_odd`; the witness's analogue is `f` on the half-system. Same `ε=−1` of the same `R`, on the two faces.

## Where this points the proof

The cleanest consequence is a falsifiable proof skeleton for the cap: **bound `M_even` by the S75e cyclotomic cosine SOS (it captures exactly the R-symmetric part), and certify `M_odd` by the Borsuk–Ulam odd degree — the literal `±1` eigenspace split.** The remaining concrete step is to show the S75e Fejér gap *equals* `M_odd` (not merely that both are R-odd-flavored); that would turn "the obstruction is the R-odd eigenspace" from *supported* into *exact*, and reduce the cap to two clean pieces in orthogonal eigenspaces.

See [[the-two-indices-redei-is-odd-lonely-is-even-half-tiling-is-the-quotient]] (THM-582), [[the-dihedral-recursion-existence-is-even-witness-is-odd]] (THM-581), [[the-topology-of-the-lrc-cap-is-euler-char-of-the-cover-nerve-lonely-is-the-hole-certified-by-borsuk-ulam]]. Theorems: THM-583 (witness half-system); HYP-3538 (cap eigenspace test).
