# Independent positive-five clock and prime-offset review

**PASS, no correction.** This report binds both frozen proofs in `../../session197-prime-atoms/positive-five-clock/`:

- `PROOF.md`: `bb7a19deed8eac1ae65feeffbdaa5ce4efd2b62634c44bc167d6c568045e5a9e`.
- `PRIME_OFFSET_BRIDGE.md`: `bc20f87bc6849c2080d999e732feecdb75f50f0721b88a81435a3195da62d345`.
- Main verifier: `cff0f01be9440fda246d2174338f6b9c4b1c50ad87a836993dce65ae2d9ee453`; main receipt: `f4e9992e89d7168f4ae062499cfc90e8c6fc62a5175e0abeb13cafadcc36fc16`; witness: `ba661cc1db521546aac253d466facde6f19b98f5babc15413855b530c0a5b75e`.
- Bridge verifier: `3a9862c38225ba9820453d84969fb0146d4e2f499f7a87f9845beb8b91ce00b7`; bridge receipt: `49a5280d61c1c17b9a486039a5853c370b199bdba9502a1eb9de26874c28a5e8`.

## All-order proof review

The physical quotient substitution and its inherited arbitrary-module coefficient formula give exactly(A,B,C)=(1152,0,5). The independent root-run certificate is sufficient for all positive clone lengths because the stated fixed singleton cover bounds the root by three runs. The same four retained vertices are exactly the old root's inneighbors; each left birth therefore adds a transitive root clone. The displayed canonical-to-birth mapping preserves every incidence and yields the literal3H-10 return, without changing ports.

The prime-offset theorem correctly uses the primality of the **offset**, not of the whole Hamilton count. The boundary statistic is nonzero, its only negative value is-1 at K1, and positive unit value forces Tr2. Factoring C=q leaves exactly the cases analyzed in the proof. The apparent alternative with delta_P=q reduces, after binary-mode cancellation, to B_P(2)=4H(P) and H(P)=q; the cut bound then forces order at most three, excluding q>=5. The q>=5 restriction is essential and is stated.

The H5 specialization correctly removes only singleton SCC padding around its unique four-vertex strong core. The already named strong-order bound is applied with H=5, excluding order at least five. Ordered-join bin multiplication then gives exactly B_W(2)=32*2^s. This proves the specified-quotient amplitude set1152*2^k without asserting global minimality or a new source for an unrealized prime offset.

## Independent computations

The independent checker imports the earlier released graph constructor, not the audited source constructor, and rebuilds all four adjacency matrices. A weighted outside-word DP with factors1_(u->v)+z1_(u->r->v) independently computes the complete root-run polynomial

    1157+3461z+3456z^2+1152z^3.

This is separate from the source's repeated-root-symbol DP. Literal odd-cycle enumeration and disjoint packings, rather than Hamilton endpoint DP, independently give3461,10373,31109,93317 for the four graphs. All are strong. Eight successive physical port births pass1,964 direct incidence comparisons.

Prime-only trial division independently tests all16 primes through58 for3461 and all30 primes through115 for13331. Thus the explicit prime and7p certificates pass without relying on the source's all-integer remainder arrays.

For the structural bridge, an independent complement-subset bin count checks all33,867 labelled tournaments through order six, yielding39 physical profile types. The H5 inventory is24,240,2160 at orders4,5,6, always with the stated two-bin profile. All2,313,441 four-profile coefficient tuples satisfy the claimed descent whenever the offset is prime at least five;204 pure prime-offset tuples occur in this larger control bank.

Both source verifiers were copied into a separate audit subdirectory before execution, because the main verifier writes its witness. Normal and optimized source outputs and witness bytes match the frozen originals. The independent checker also matches normally and under `-O`.

Independent checker SHA256: `ade60c23a6cff9991e60199ed09fcd2e3dae6d4931c4222698f374ac4bd6ffbb`.
Independent receipt SHA256: `3751beac856f53ed031d52a8e2d28f8e333ca7ab6c5d233e707fa8a0a29dd6ed`.

No source file was edited. No unrestricted prime-spectrum, prime-infinitude, or global least-amplitude claim is added.
