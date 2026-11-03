# Independent lattice audit

**Verdict: APPROVED**, including the dense-orthogonal-closure strengthening.
Reviewed complete original theorem.md, density-extension.md and hardened
verify.py before the temporary-workspace reset. Restoration review will be
recorded separately below.

The cyclic odd-order tournament followed by deletion realizes a balanced
near-regular score pattern; relabellings realize every pattern. Degree sums
prove the converse. Adjacent shell vectors have 2cos(theta)=2-8/n. Finite
reflection product makes exp(2i theta), hence exp(i theta), a root of unity;
the rational displayed value must be an algebraic integer. For even n>=6
its only possible integer value is 1, yielding n=8. The n=2,4 controls are
correctly separated.

Shell differences generate A7. Adding a half root has index 2, determinant 2,
and exactly 126 norm-two vectors. Integral pairings ensure reflections
preserve this finite set; all 56 added integer roots arise in one step.
The stated Gram graph is E7. A finite reflection closure would be a finite
spanning set with a faithful group action, proving the infinity assertion.

The density lemma is valid. Irrational adjacent angles give full rotation
circles in the closed group. From u wedge w and so(U), commuting with
(a wedge u), a perpendicular to u, gives ||u||² a wedge w up to sign.
Together these generate U wedge w. Johnson-graph connectivity and spanning
therefore yield so(V). The standard closed subgroup theorem gives SO(V),
and one reflection gives O(V). Orbit density follows. This does not identify
the countable group itself with its uncountable closure.

Initial bare asserts were optimization-disabled. Author hardened to explicit
checks, and auditor reran both normal and optimized scripts with identical
outputs matching verification.json before the reset.

Restored independent checker lattice_independent.py imports no author code.
It recreates 240 ambient roots, all 15,876 reflection pairs, one-step closure,
and every labelled tournament on 2,4,6 vertices to verify score-image equality.
Normal/-O both pass 31,759 explicit checks with identical outputs.
These finite checks do not substitute for the all-n or density proofs.

No global novelty is asserted. Dependencies of the density statement are
standard matrix-Lie-algebra and closed-subgroup facts; the finite theorem
uses only elementary algebraic integers and conventional root-system naming.
No implication to tournament automorphisms or LRC is introduced.

## Restored-file review

Read the complete restored theorem.md, density-extension.md and verify.py at
/workspace/shared/session174/lattices. The proofs retain the audited content.
Reran restored author script in normal and optimized modes; both outputs are
identical and match restored verification.json. **Final verdict: APPROVED.**
For prose clarity, the last density paragraph should say explicitly "K also
contains a reflection", since the identity component SO(V) does not.
