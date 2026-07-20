# The third dihedral: the Jacobian fibre is a tournament

*kind-pasteur-2026-07-20-S128c105. Owner: "consider previous work regarding the
relationship between tournaments and dihedral groups thoroughly, there will be
multiple instances of related concepts to explore."*

There are. A full sweep of canon finds **two** established dihedral phenomena in
this repo, joined at exactly one point. I want to argue there is a **third**, that
it lives in the Jacobian thread, and that it has been misfiled as a coincidence.

## The two known ones

**The vertex-level `D_n`.** `⟨r: v↦v+1, s: v↦-v⟩` acting on a circulant
tournament, with rotations as automorphisms and reflections as *anti*-automorphisms
onto `T^op`. This is THM-127 (`Aut ∪ anti-Aut = D_{2p}` for Paley `T_p`, `p≡3 mod
4`), and it is the same group doing arithmetic in THM-586 (Burnside on Hamiltonian
paths, `p | H`, `H/p ≡ f mod 2`), representation theory in THM-131 (`189 = 18ρ₀ ⊕
9ρ₁ ⊕ 27ρ₂ ⊕ 27ρ₃ ⊕ 27ρ₄` for `D₁₄`), and isotypic decomposition in THM-581
(LRC existence is `σ`-even, witness construction is `σ`-odd).

**The tile-level `S₃ = D₃`.** Acting on the staircase `Grid(n)` by `σ(r,c)=(c,r)`
and `τ(r,c)=(r,n-r-c)`, with three reflection axes meeting at the centroid and
`|Fix(reflection)| = 2^{⌊(n-1)²/4⌋}`. Order 6 independent of `n`; its "rotation" is
the barycentric 3-cycle, not a vertex shift. A genuinely different group.

**The one joint.** THM-280: the grid reflection `(x,y) ↦ (n+1-y, n+1-x)` induces
`T ↦ T^op`. So one reflection of the tile-level `D₃` *is* the reflection of the
vertex-level `D_n`, and is exactly the `ℤ₂` whose quotient gives the merged
metagraph `G_n/ℤ₂`. The RED transpose edges in the explorer are simultaneously
`D₃`-reflection edges and `D_n`-anti-automorphism edges. That is the whole of the
known picture.

## The third one, and why it was filed as a pun

A thorough sweep ranks THM-1310's `S₃` — the Galois/monodromy group of the
Jacobian counterexample's fibre cubic — as an **order-6 coincidence, not dihedral**,
on the grounds that there is "no rotation/reflection-of-a-polygon structure" and
"neither touches tournaments." Both grounds are wrong, and the same theorem says
why: THM-1310 proves the generic fibre is a **cyclic 3-tournament**.

Take that literally. The fibre is a 3-element set carrying a 3-cycle `1→2→3→1`.
Then

- `Aut(3-cycle) = ℤ₃ = A₃` — the rotations, arc-preserving;
- the three transpositions are **anti**-automorphisms — arc-reversing, carrying the
  3-cycle to its opposite;
- so `Aut ∪ anti-Aut = S₃ = D₃`.

That is not "a group of order 6". It is *the same rotations-preserve /
reflections-reverse split that THM-127 proves for `T_p`*, at `p = 3`. The Jacobian
fibre carries the vertex-level dihedral structure, on three vertices.

## What the identification buys

It is not decorative — it re-reads a classical theorem into the repo's own idiom.
Campbell (1973): *a Keller map that is a Galois cover is an automorphism.* At cover
degree 3, Galois means monodromy inside `A₃`. So:

> **Rotation-only Keller covers are trivial.**
> Equivalently: **every degree-3 Keller counterexample is reflection-driven** — its
> monodromy must contain an arc-reversing element of the fibre tournament.

And the discriminant character that detects this is exactly the rotation/reflection
indicator — the same parity bookkeeping as Rédei's sign on the tournament side. The
repo had "√(-L) is a non-square, 0/8316 sample points" filed as *verified*. Under
the dictionary it is *forced*, and the 8316 points are a corollary (THM-1375).

Degree 3 is also the only place this works: `A₃` is regular, so transitive-in-`A₃`
implies Galois, whereas `A₄` is transitive on 4 points without being regular. The
character detects Galois-ness at `d=3` and nowhere else.

## The negative that matters more

I reached this by the wrong road, and the wrong road is instructive. THM-1300
records the counterexample as `ℂ*`-equivariant with source weights `(1,-1,-2)` and
target weights `(-2,-1,1)` — related by the transposition `(1 3)`, *a reflection*.
That looked like a second, independent appearance of the same reflection, and I
proposed a reduced Jacobian Conjecture around it: untwisted equivariance forces
invertibility.

It is false in one line. Compose with the target coordinate swap: `σ∘F` is Keller,
non-injective, and untwisted. **The weight twist is a coordinate artifact** — it is
not invariant under composition with a linear automorphism, while "Keller" and
"non-injective" both are.

I have now made this exact error twice in unrelated threads. In the LRC work I
spent four sessions looking for what made the direction `(1,2,3)` special, and the
answer was that nothing did: it was the unique *increasing* representative of a
six-fold symmetric orbit, distinguished only by an ordering convention. Here I took
a distinguished-looking weight vector for structure, and it was a labelling.

The generalisation is worth carrying: **when a structure looks distinguished, first
ask which group is allowed to act.** If the apparent distinction is not invariant
under that group, it is a name, not a fact. Both times the real content was one
level down and coordinate-free — the orbit in the LRC case, the monodromy here.

## What this suggests to look at next

The two reflections in the Jacobian story are *not yet known to be the same
reflection*. The `λ = -1` element of `F`'s torus fixes the collision image and
transposes two of its three preimages (exact, verified) — a reflection of the
special fibre. Campbell forces a reflection in the *monodromy*, a loop in the base.
Whether the torus involution generates the monodromy is open, and it is the natural
next computation: if it does, the counterexample is "reflection = torus" and the
`ℂ*`-anatomy stops being an accident of presentation.

More speculatively: if the fibre-as-tournament dictionary is real rather than
notational, then the repo's tournament invariants should have Jacobian meaning. The
first test is cheap — Rédei's theorem says every tournament has an odd number of
Hamiltonian paths, and boxeph-S146's odd-degree conjecture ("every Keller
counterexample has odd cover-degree") is already sitting there as a Rédei-parity
transfer. Two independent threads reaching for the same parity statement is either
the dictionary working, or the second time this session that I have mistaken a
name for a fact. Worth finding out which.
