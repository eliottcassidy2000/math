# A weighted-automaton obstruction to finite affine gadget banks

Status: elementary theorem about candidate interface systems, not a tournament construction or a new spectrum bound. It is an obstruction to a *proposed* finite-state induction from finite seeds. The count-growth argument is a graph-directed version of the Kraft/Dirichlet-series test; no attribution of novelty is claimed. Its graph-interface formulation follows [`finite_certificate.md`](../strategy/finite_certificate.md) and the prime-lane context is [`word_interface.md`](../strategy/word_interface.md).

Let S be a finite set of states, E a finite directed multigraph of rules e:s→t, each with integers a_e≥2 and b_e, and let the initial set consist of finitely many pairs (s,h). The rule e may be applied to (s,h) only when its guard is satisfied, giving (t,a_e h+b_e). Guards may be arbitrarily restrictive or depend on the entire history: dropping them only increases the generated set. Suppose there exists a real C with

    C ≤ b_e/(a_e−1) for every e, and h+C>0 for every seed.

For σ>0 put M_σ(t,s)=Σ_{e:s→t} a_e^(−σ), summing parallel edges. If the spectral radius ρ(M_1)<1, then for some σ<1 the number of distinct generated values h≤X, over all states, is O(X^σ). Consequently this finite bank cannot generate all sufficiently large odd integers or all sufficiently large primes; nor can it generate all sufficiently large elements of the prime lanes p and 7p described in the cited strategy.

Proof. By continuity of spectral radius, pick 0<σ<1 with ρ(M_σ)<1. For every transition,

    h'+C = a_e(h+C) + b_e−(a_e−1)C ≥ a_e(h+C).

Thus, along any rule history from a seed h₀, (h+C)^(−σ) ≤ (h₀+C)^(−σ)∏a_e^(−σ). Sum this bound over *all* finite paths in the finite state graph, retaining multiplicities and disregarding guards. The resulting sum is bounded by a finite seed-dependent constant times Σ_{k≥0} 1ᵀ M_σ^k 1, which converges because ρ(M_σ)<1. Therefore the corresponding sum over distinct generated values converges as well. For X+C>0, each generated h≤X contributes at least (X+C)^(−σ), giving the polynomial counting bound. The exclusion of all large primes uses Chebyshev's standard elementary lower bound π(X) ≫ X/log X; removing a fixed finite set and replacing p by 7p does not alter the order. This argument never assumes any target has a tournament witness, any rule exists graph-theoretically, or any witness family is closed under a gadget. If the allowable seeds do not satisfy the stated common-shift condition, the theorem does not apply without an independently justified transient/threshold reduction.

## Two-state example: locally complete ternary branch, globally too sparse

Let A and B each start with the single seed 23. Supply A→B rules h↦3h, 3h+2, 3h+4, and the only B→A rule h↦9h. All maps preserve oddness and strictly increase positive inputs. The three A→B images cover the three odd residue classes modulo 6, so the familiar local ternary residue check succeeds *for B*. Nevertheless, with C=0,

    M_σ = [[0, 9^(−σ)], [3·3^(−σ), 0]],
    ρ(M_σ) = 3^((1−3σ)/2).

Hence for every σ>1/3 there are O(X^σ) generated values at most X. The state-switch return path is the obstruction: ternary branching happens only once per multiplicative scale 27. It fails even to cover all sufficiently large primes despite perfect local target-residue coverage for B. This is an abstract arithmetic counterexample, not a claim that any edge is realizable by tournament gadgets.

The test ρ(M_1)≥1 is only necessary, not sufficient: equality or a larger spectral radius does not certify residue coverage, seed reachability, state invariance, or graph operations. The theorem has no bearing on adaptive constructions with infinitely many interface types or non-affine rules unless they are separately encoded by a finite system satisfying its hypotheses. In a candidate graph, compute the characteristic polynomial or prove a positive rational vector v with M_1v≤qv for some rational q<1; either yields an exact finite certificate of failure. The latter bound implies ρ(M_1)≤q by weighted sup norm.
