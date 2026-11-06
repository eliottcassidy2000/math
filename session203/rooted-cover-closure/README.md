# Rooted cover closure and guarded unit-start compiler

**PROVED transfer and order bound; independently audited.** A nine-coordinate physical state closes under the actual O/Z/J operations, retaining B=H(U)+L-start_root(U) and the necessary cover sidecars. A four-coordinate nonnegative cover block gives the same new information.

Removing a final O suffix preserves (B,s). A final Z or J step forces a source representing B<=M to have order at most floor(log2 M)+2, with the seed order retained for an empty word. Combined with the inherited complete unit-start graph grammar, this yields an exact bounded compiler, not an occupancy theorem.

- B=137 has no unit-start profile at any order. This excludes neither its Hamilton-count realization nor a staggered profile with another start or shifted base.
- B=139 has a literal unit-start witness.
- Exactly 330 unit-start bases occur through 997.
- A resource limit returns INCOMPLETE and never licenses an absence claim.

PROOF.md gives the proof and scope. The two verifier pairs have matching normal/optimized receipts. audit/ contains independent permutation-and-cut checks, unpruned word enumeration, and frozen graph-level dependencies. The universal staggered lemma F(8,137) and eventual threshold remain OPEN.
