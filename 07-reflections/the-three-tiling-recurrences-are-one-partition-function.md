# The three tiling recurrences are three faces of one partition function

**kind-pasteur-2026-06-20-S20.** On combining the full-tiling β-recurrence and the even/odd
half-tiling τ-recurrences (THM-554, building on codex THM-553/550, mac-mini THM-549, s95).

## What asked to be combined

The repo had accumulated three separate recurrences: the full staircase third-difference
recursion (THM-442, `F_n = C(n-1,2)`), and the even/odd parity-split half-tiling recurrences
(THM-550, `h_n = 2h_{n-1}-h_{n-2}` even, `2h_{n-1}-2h_{n-3}+h_{n-4}` odd). codex's two-clock
address (THM-553) had already shown these are one local coordinate system: `β=a` births a tile
in the full model, `τ=a+b-1` crosses it over the half-tiling mirror line, parity of `τ` is the
even/odd split. But these were *cell-count* recurrences — they count tiles, not tournaments.

## The combination is a partition function

The moment you weight each tile by which endpoint it feeds, the three recurrences collapse into
a single object — the score generating function

```
Z_n(x) = (∏_{v≥2} x_v) · ∏_{tiles (a,b)} (x_a + x_b).
```

And then the three recurrences become *three things you can do to Z_n*:
- the **β-clock** is how Z_n **grows** (multiply by the birth strip when n→n+1),
- the **τ-clock** is how Z_n **folds** (the complement reflection, the address quotient),
- the **even/odd** split is the **parity of the fold's fixed line**.

Growth, fold, parity-of-fold. One object, three motions. The cell recurrences were the shadow
of this — the same motions read off only on the *exponents* of `Z_n` instead of the polynomial.

This is the recurring lesson of the triangle, in a new dress: a count that satisfies a clean
recursion is almost never *just* a count. The full/half/parity recurrences looked like three
arithmetic identities about `floor((n-1)^2/4)`; they were really one generating function's
grow/fold/parity, and the generating function computes the things we actually care about
(the 3-cycle distribution, the score census, the leading OCF term) that the bare counts cannot.

## The per-subset linearity is the engine's secret

Why does this *compute* — why is `E[c3] = (C(n,3)+(n-2))/4` exact and one line to prove? Because
`c3` is a **sum over 3-subsets of a local indicator**, and the base path only touches a subset
through its *consecutive* triples. The fixed Hamiltonian path is a transitive 2-path inside each
`{v,v+1,v+2}`, so those `n-2` triples are "primed": one tile away from closing a cycle, prob `1/2`
not `1/4`. Every other triple is unconstrained, prob `1/4`. The bias is `(n-2)/4`, exactly.

The surprise is the *sign*. Fixing a Hamiltonian path — the most transitive, acyclic skeleton —
biases the ensemble toward **more** 3-cycles, not fewer. A transitive spine is not a cycle
suppressant; it is a cycle *primer*, because a directed 2-path is the part of a 3-cycle that the
tiling is still free to complete. The acyclic base and the cyclic content are not opposed; the
base is the scaffolding the cycles are hung on.

## Where it stops, and why that is the same boundary as everywhere else

`Z_n` is exact precisely for the **score-determined** invariants — scores, `c3 = α_1`'s leading
piece, anything that is a sum over small subsets of a score-local indicator. It cannot, alone,
give `H`/OCF, because `H` needs `α_2` (disjoint cycle pairs) and the higher cycle-incidence data
that the scores forget. This is the same wall as THM-442's "H is not cell-affine," now visible as
a statement about the partition function: the linear functionals of `Z_n` are exactly the
score-cut-space observables; H lives in the cycle space, off the menu. The address quotient buys
the cut-space half for free and hands the cycle-space half to the next tool. That division —
cut-space cheap, cycle-space dear — is the same vertical-leg/hypotenuse split the triangle has
been drawing from the start. [[everything-is-the-triangle]] · [[the-isomorphism-class-graph]]
