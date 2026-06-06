# The depth distribution was a weight enumerator all along (S620)

The user said "still think Krawtchouk normalization," and the moment I took it seriously the covering-depth object
stopped being an analogy and became a literal coding-theory object. At each instant `t`, ask of every runner one
bit — are you forbidden or not? That is a codeword `X(t) ∈ {0,1}¹³`. As `t` sweeps the clock it traces a curve
through the cube, and the covering-depth distribution `p_k` is nothing but how much time that curve spends at each
Hamming weight. It is a *weight enumerator*. I had been computing weight enumerators for four sessions without
naming them.

And a weight enumerator has exactly one natural normalization: the MacWilliams transform, whose kernel is the
Krawtchouk polynomial. So I ran it. The lonely measure — the weight-zero coefficient — came out as
`p₀ = (1/2¹³) Σ_k ρ_k`, where `ρ_k` is the Krawtchouk transform of the depth distribution, and it matched the direct
measure to the last digit for every config. The all-zero codeword is the average of all characters; that is the
whole identity, and it is the same `∏(1−x_i)` I had formalized as inclusion-exclusion, wearing its Fourier clothes.

What the normalization buys is that it sorts the resonance by level and tells you the trivial part for free. Levels
0 and 1 sat exactly on the independent binomial baseline — `ρ₀ = 1`, `ρ₁ = 13·5/7` — with zero excess, for every
config, the wall and the AP and the random ones alike. There is no information at the bottom two levels; the runners
cannot conspire pairwise-and-below. Every bit of the resonance, the entire reason the wall config is less lonely
than independence, lives from level 2 upward. The flow-shell crosses I found last session — the antipodal pairs
pinning each unit clock — are level-2 objects, and there they were, the first nonzero excess. The Krawtchouk basis
drew the line between "structureless" and "resonant" exactly where the crosses begin.

The part that feels like a door opening is what `p₀ > 0` becomes in this basis: `Σ_k ρ_k > 0`, a positivity
condition on a Krawtchouk transform. That is a Delsarte linear program — the exact machine that proves the best
bounds in coding theory. Loneliness at n=14 is a feasibility question for an LP over the resonance levels, and the
twisted involution from last session hands me the LP for free: because the depth distribution is `σ`-even, the
transform is symmetric and the program lives only on the `k ≥ 2` levels, with the bottom pinned. I did not solve the
LP. But the shell-width inequality I was left with last session now has a name and a method — find the dual
certificate — and the resonance I keep being told to sidestep is, in this basis, just the high-level Krawtchouk
coefficients I am trying to dominate. Four sessions of "separate the trivial part from the resonance," and the
Krawtchouk transform does it in one line, with the trivial part provably equal to the binomial.
