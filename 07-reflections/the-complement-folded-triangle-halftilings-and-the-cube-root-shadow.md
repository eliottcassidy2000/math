# The complement-folded triangle: half-tilings, cube roots, and the one/two/three-far shadow

**Source:** mac-mini-2026-06-20-S4. Dispatch: the human owner's half-tiling framework (mirror over
y=x = reverse all arcs; sizes 0,1,2,4,6,9,12,…; even vs odd recursions A+B−C vs A+B−C+D−E−F+G),
to be connected to the full-tiling A+B+C−D−E−F+G recursion, the one/two/three-far structure, and
the cube roots of unity. Canon: THM-549. Built on codex's HYP-2680/2681 (Φ_s Stirling hierarchy,
Eisenstein cube-root modes) and THM-548 (boundary-value decomposition).

## One fold, named three ways

The deepest thread of this project has been *everything is the triangle* — the staircase
`δ_{n-2}` is a right isosceles triangle, and tournament structure is binary tilings of it. The
owner's half-tiling framework adds the move that was missing: **fold the triangle by its own
mirror.** Reflecting a tiling over the line `y=x` is, exactly, reversing every arc of the
tournament — the complement `T → T^op` (PROVED: `σ(x,y)=(n+1−y,n+1−x)` equals `φ∘(reverse-all)`,
verified over all `2^m` tilings for `n≤6`, THM-549). So the half-tiling is the **fundamental
domain of the complement involution**, the same `Z_2` that the merged metagraph `G_n/Z_2` factors
out — but now realized as a literal mirror, with the self-complementary tournaments sitting on
the literal diagonal `{x+y=n+1}`, the spine.

The fold changes the shape from a triangle (`C(n−1,2)` cells) to a quarter-square
(`⌊(n−1)²/4⌋`), and it splits cleanly by parity: **odd `n=2k+1` folds to a perfect square `k²`,
even `n=2k` to a pronic rectangle `k(k−1)`.** That is the owner's "even tournaments have a
slightly different shape" made exact: the square has three recursive corners (`A`,`D`,`B`) and a
seven-piece decomposition; the pronic has none. The order-4 recurrence `(x−1)³(x+1)` governs all
`n`; even `n` degenerates to order-2. Odd is generic, even is special — the same parity asymmetry
the LRC carries through `14 = 2·7`.

## The seven terms are the one/two/three-far

The owner asks us to see the half-tiling recursion next to the full-tiling
`A+B+C−D−E−F+G`. They are the same skeleton. codex (HYP-2681) had already found that the
full-tiling seven-term recursion *is* the coverage pair-tax shadow `H(1)−2(D+E+F)`, with

```
A, B, C  =  one-far residuals  (Δ_u, Δ_v, Δ_w)
D, E, F  =  two-far curvatures (I_uv, I_uw, I_vw)
G        =  three-far residual (Δ_uvw)
```

Three, three, one. The half-tiling's seven pieces — three corners, three edges, one center —
carry the identical `3+3+1` grading, and the difference between the two sign words
(`+++−−−+` full, `++−+−−+` half) is precisely the trace of the complement fold flipping the
level-`n−2` pair. The combinatorial recursion of the folded triangle and the analytic recursion
of the coverage are one object seen upstream and downstream: **the one/two/three-far hierarchy of
the lonely-runner coverage is the measure-theoretic shadow of the half-tiling's recursive
geometry.** That is why the coverage expansion has exactly three one-far, three two-far, one
three-far, with those signs — it inherits the count and the grading from the folded triangle.

## The cube roots are the three runners

A triangle has a three-fold symmetry, and three-fold symmetry is the cube roots of unity. On the
coverage side this is the symmetric group `S_3` permuting the three far runners `{u,v,w}`; the
Eisenstein modes `S_ω = A+ωB+ω²C` and `P_ω = D+ωE+ωF` (codex) are the `C_3 ⊂ S_3` characters,
and their moduli are `S_3`-invariant (VERIFIED: `|S_ω|` constant across all six permutations).
The cube root `ω` is not decoration — it is the rotation of the three corners of the folded
triangle into each other, and the conjugate pair `{ω, ω²}` is the two equal corners `A,B` (size
`n−1`) against the lone corner `D` (size `n−2`), exactly the way `1` is real and `ω,ω̄` conjugate.
The owner's "A, D, B in each corner" is the cube-root basis of the recursion.

## The order is the apex prime

One more thing fell out, and it ties the grading to the number that has run through everything.
The Newton packet of order `s` (an `s`-far curvature) **vanishes unless `|B| + s ≥ 7`** — you
need seven runners to cover seven sectors (VERIFIED: for `|B|=5`, all one-far packets are exactly
zero; coverage first switches on at two-far). So the far-order at which the recursion becomes
live is set by `7 − |B|`: the apex prime is the *length* of the expansion, just as it is the
modulus of the sectors, the period of the kernel's vanishing, and the denominator `7^{s+1}` of
the `s`-far constant. The half-tiling tells us why the seven-term recursion has seven terms
worth keeping near the cap: the triangle has three corners, its fold has a three-fold symmetry,
and the sector count seven is what truncates the far-order tower.

## What this buys, and what to chase

Concretely proved: the fold is the complement, so every complement-invariant tournament invariant
(cyclic-triangle count, Hamiltonian-path count, `H`, the OCF, the order-2 Walsh spectrum —
`c3(T)=c3(T^op)`, `HP(T)=HP(T^op)` verified at 100%) is computable on the half-region, a literal
`2×` saving with the SC spine handled once. Conjecturally ahead: that the coverage's
seven-term recursion can itself be *folded* — a half-coverage with the `++−+−−+` sign word — to
halve the resonance ledger the LRC proof now turns on; that the odd-square / even-pronic split is
the geometric source of the LRC's odd-skeleton discrimination (`14=2·7`, the seven odd speeds);
and that the cube-root three-fold is a genuine third reduction axis (a Mode C, `n→n−3`, Eisenstein
beside the binary Cayley–Dickson `n→n−2`). The triangle was never just a triangle; folded, it is
the place where the seven sectors, the three runners, and the two parities meet.
