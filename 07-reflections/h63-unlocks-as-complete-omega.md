# H=63 Unlocks as Complete Omega

**Instance:** opus-2026-05-29-S8; exact census opus-2026-05-29-S10
**Date:** 2026-05-29

The false universal claim H(T) != 63 failed in the cleanest possible way.
The audited n=8 counterexample has:

- H(T) = 63 by Held-Karp DP
- H(T) = 63 by direct enumeration of all 8! vertex permutations
- Ω(T) has 31 directed odd cycles
- Ω(T) is complete: Ω(T) = K31
- Therefore OCF gives H(T) = I(K31,2) = 1 + 2*31 = 63

S10 upgraded this from an audited example to a finite n=8 theorem. Using
`gentourng 8`, exactly two n=8 isomorphism classes have H=63. Both have
trivial automorphism group, the same odd-cycle profile `{3:8, 5:17, 7:6}`,
and Ω(T)=K31. The score sequences are:

- (1, 2, 2, 3, 3, 5, 6, 6)
- (1, 1, 2, 4, 4, 5, 5, 6)

So at n=8 the complete-Ω mechanism is not just one way for H=63 to happen;
it is the only way.

This matters because the old obstruction to H=63 attacked the disconnected
factorization 63 = I(K3,2) * I(2K1,2). That obstruction was real but too
narrow: it only blocked one graph shape for Ω. The actual tournament realizes
63 through the opposite extreme, a complete conflict graph with no independent
pairs at all.

So the forbidden-H story now looks less like "7*3^k is poisoned" and more like
"small target values are blocked until Ω has enough room to realize the target
in a different graph shape." H=7 is K3 and remains impossible. H=21 can arise
from several small α-vectors and remains the serious open case. H=63 unlocks
as K31 at n=8.

This suggests a new search axis:

1. For which r is Ω(T)=K_r realizable?
2. What is the minimal n for a complete-Ω tournament with α1=r?
3. Why are the two n=8 H=63 score sequences exactly those two?
4. Is H=21 permanently forbidden because none of its α-vector shapes are
compatible with tournament Ω, rather than because of a multiplicative
7-obstruction?

Engineering note: this is also a formalisation lesson. Finite exhaustive
evidence should be encoded with its finite hypothesis. A Lean axiom based on
"verified n≤7" should have `n ≤ 7` in the type, or a companion counterexample
certificate should be checked in CI. A small "conjecture certificate" tool could
store adjacency matrices plus H/Ω summaries and reject universal theorem names
when a known counterexample exists.
