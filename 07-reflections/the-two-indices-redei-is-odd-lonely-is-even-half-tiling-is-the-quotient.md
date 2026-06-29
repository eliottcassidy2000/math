# The two indices: Rédei is the odd half, lonely is the even half, and the half-tiling is the quotient

*mac-mini-2026-06-29-S6. Synthesis from the thread "palindromic Hamiltonian paths and the half-tiling model," following it into the LRC floor. New lemma: THM-582. Builds on THM-281, THM-549/550, THM-088, THM-581, HYP-3242/HYP-3239.*

## One involution, two indices

The project has a single involution wearing many costumes — complement `T↦T^op`, grid reflection `(x,y)↦(n+1-y,n+1-x)`, path reversal, vertex relabeling `k↦n+1-k`, the LRC complement `t↦−t`, the half-translation `t↦t+1/2`. Call it `σ`, the order-two generator of the dihedral group `D_p` that `n=2p` is the order of. Everything the project computes is a function on a space carrying a `σ`-action, and **every invariant splits into a `σ`-even part and a `σ`-odd part.** The long-running confusion — "does the proof need Brouwer or Borsuk–Ulam, SOS or odd degree, an even index or an odd one?" — dissolves once you see that there are *two* indices and they answer *two different questions*.

**The σ-ODD index = the witness = Rédei.** For a self-converse tournament, reversal-conjugation `ρ(P)=φ(reverse(P))` is an involution on Hamiltonian paths (THM-582), so `H(T) ≡ #{palindromic paths} (mod 2)`, and Rédei's "`H(T)` is odd" *becomes* "there is an **odd** number of palindromic Hamiltonian paths." Verified: Paley `T_7` has `189` paths, `9` palindromic; `T_{11}` has `95095`, `185` palindromic — all odd. This is the same fact THM-281 proves one level down (self-complementary tiling fibers are odd because the grid-symmetric tilings are odd), and the same fact kps's `D_7` Borsuk–Ulam guards (the free-`Z_2` odd degree = the imaginary Gauss sum `i√7` = the sign isotype). Three statements, one odd index — the fixed-point parity of the reversal `σ`.

**The σ-EVEN index = existence = the lonely measure.** The LRC cap is the measure-weighted Euler characteristic of the danger-cover nerve, `meas(lonely) = Σ_S(−1)^{|S|}meas(∩D_p) = χ_meas(nerve)` (HYP-3242, verified; and it is exactly my comb-resonance inclusion–exclusion of THM-578/579). Because lonely sets are `σ`-invariant, `χ_meas` has **no sign component** — it is purely `σ`-even (THM-581). Existence of a lonely time is therefore a Brouwer/SOS-category statement; it does **not** see the odd index at all.

## The lemma settles the "odd index" question — and removes it from the floor

I had been hunting (two-order-two-structures reflection; lrc19-s679's "fixed-point-of-involution residual") for an *odd index* that would force the LRC floor non-constructively, a "Rédei-flavored `H(T) mod 2`." THM-582 **finds that index** — it is the palindromic Hamiltonian-path count — **and simultaneously shows it does not belong to the floor.** The reason is sharp and structural: a *lonely* tournament has the observer `0` as a **source**; the converse turns a source into a **sink**; so a lonely tournament is **not self-converse**, `ρ` is not an involution on its paths, and THM-582 simply does not apply. The floor lives entirely on the even side.

This is the most useful thing the palindromic detour produced: it *closes off* the odd-index route to the floor, with a reason, and so **redirects all floor-closure effort to the even-category route** — the 2-adic descent (THM-580, which is itself the constructive even-category degree) plus the cyclotomic Fejér–Bochner SOS (S75e, which handles the even/2-dim de Moivre part). No Borsuk–Ulam, no palindromic parity, for existence.

## The half-tiling is the quotient where the two indices live

Why was the half-tiling model the right place to look? Because it is the **`σ`-fundamental domain itself** (THM-549/550): `h_n = ⌊(n−1)²/4⌋` cells = (one side of the reflection) + (the fixed diagonal `x+y=n+1`), the Burnside orbit count of `σ` on the staircase. The fixed diagonal is the **self-complementary spine** — the palindromic/odd locus — and everything off it pairs up under `σ` (the even, free part). So the half-tiling literally *is* the bookkeeping of "fixed (odd) ⊕ paired (even)" for the project's central involution. The parity recurrence `h_n = 2h_{n-1}−h_{n-2}` (even `n`) vs `2h_{n-1}−2h_{n-3}+h_{n-4}` (odd `n`) is the combinatorial shadow of the even/odd index split — odd `n` carries the extra `(x+1)` factor (the sign mode), even `n` is degenerate (pronic), exactly as `THM-088`'s `SF` is anti-palindromic for `n≡2,3` and palindromic for `n≡0,1`.

## The bold surmise (freely, as asked)

> **Rédei's theorem and the Lonely Runner Conjecture are the odd and even halves of a single `σ`-equivariant Euler characteristic of a conflict-cover complex.**

On the tournament side the conflict complex is `Ω(T)` (odd cycles) and the invariant is `H = I(Ω,2)`, whose `σ`-odd (palindromic) part carries Rédei's oddness. On the runner side the conflict complex is the danger-cover nerve and the invariant is `χ_meas`, whose `σ`-even part is the lonely measure. The two problems are not analogous by coincidence; they are the two isotypic projections of one equivariant index on (different incarnations of) the same `D_p`. Proving LRC(14) is *computing the even projection is positive* (descent + SOS); Rédei is *the odd projection is odd* (palindromic, done). The half-tiling is the quotient that separates them. If this surmise is right, the "unified theory of parity in tournaments" the project keeps circling is precisely the statement that these two indices are the two halves of one equivariant Euler characteristic — and the remaining mathematical content of LRC(14) is entirely on the even side, where the descent has already reduced it to a per-level cyclotomic SOS bound.

## Where this leaves the proof

- **New & proved:** THM-582 (palindromic odd index) — the path-level twin of THM-281; resolves the long-open "odd index" question.
- **Clarified:** the LRC floor is even-category; the odd index does not apply to it (the lonely tournament is not self-converse). Floor closure = descent (THM-580) + cyclotomic SOS (S75e), even category only.
- **Surmised:** Rédei (odd) and LRC (even) are the two isotypic halves of one `σ`-equivariant Euler characteristic; the half-tiling is the quotient.

See [[two-order-two-structures-parity-and-descent]] (THM-580), [[the-dihedral-recursion-existence-is-even-witness-is-odd]] (THM-581), [[the-topology-of-the-lrc-cap-is-euler-char-of-the-cover-nerve-lonely-is-the-hole-certified-by-borsuk-ulam]] (HYP-3242), [[everything-is-the-triangle]]. Theorems: THM-582, THM-581, THM-281, THM-549/550, THM-088.
