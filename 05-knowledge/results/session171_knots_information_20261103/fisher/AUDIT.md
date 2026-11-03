# Independent audit

**Verdict: ACCEPT mathematical substance.**

Read NOTE.md and check.py. Replayed normal and optimized runs independently; both exactly match the corresponding stored outputs. Replayed output SHA-256: e0a74b76e37361ca86e2948b5963aab4f7aacaa94d9d051d9cfceadb48e95be0.

Checked rank–nullity and uniform fiber conditioning, arbitrary-prior upper bounds and equality cases, Fourier annihilator projection, differentiation of the finite exponential family, the Fisher matrix null-space statement, and identity-tilt fractions. These are correct with the stated fixed-bank, uniform-prior, theta=0 restrictions. For the level-four one-point bank, full=15/256, retained=1/256, lost=14/256; the retained *fraction* is 1/15. Suggested clarifying the prose to avoid calling 1/256 a ratio.

Adaptive scope is correct: for a deterministic known branching/stopping rule, a complete realized transcript imposes exactly its queried linear equations. Every compatible source follows that branch, so no unrecorded selection bias invalidates the affine-fiber calculation. The note explicitly does not extend the fixed-bank Fourier subspace formula to adaptive trees. The all-zero identity-certification lower bound needs K independent binary probes and is properly restricted to binary linear probes rather than arbitrary-valued contexts.

The checker is faithful to its declared test scope: all characters at K=4, a declared subset at K=11, exact Fraction arithmetic, complete RM code enumeration, and one adaptive tree. It does not purport to exhaust all adaptive algorithms or general noise models. No assertions disappear under optimization.

An additional independent check, saved alongside this report as independent_audit.py and independent_audit.json, uses three simultaneous statistics on the affine three-variable code, independently computes full/observed/conditional covariance matrices, and verifies every missing Fourier outer-product entry. Normal and optimized outputs agree. The all-family proof, not this small test, establishes the theorem.
