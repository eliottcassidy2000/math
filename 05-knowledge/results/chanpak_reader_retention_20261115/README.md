# Chan–Pak: exact repairs and the cost of retaining them

This follow-up starts from the [fixed-response modular windows](../chanpak_modular_kernels_20261115/README.md). Rather than extending a coverage census, we randomly sampled older work, used it to choose new coordinates, and tested whether small integer repairs are still physical under a subsequent read.

## Exact integer repairs

- **A five-vertex source, two consecutive odds at every height.** [Local reader fork](LOCAL_FORK.md) uses one physical source with two marked sets, giving $2^{t+5}+5$ and $2^{t+5}+7$, $t\ge1$. It shares the preclosure count, two-slot count, marked response, and marked end count. But no common nonempty tail makes the closed two-slot counts agree. Even separate fixed tails cannot make them agree throughout subsequent common transitive padding: the constant terms have different parity. Thus this is a final-unit selector, not yet a repeatable switch on the same interface.
- **Compensation rather than strict equality of the intermediate response.** [Compensating reader](COMPENSATING_READER.md) retains exactly the observable the exterior reads, not the internal two-slot count. It gives a second all-height pair, starting at 29,311 and 29,313, by changing the physical host. It also gives separate same-two-slot pairs of fixed gap 6, and records the empty-bin correction that makes an initially plausible gap-2 example fail.
- **Moving both tapes can recenter a changed reader.** [Dilated recentering](DILATED_RECENTER.md) gives four literal height-independent shifts $-11126,-12118,-13786,-14958$. The common-high-word reciprocal dilation model admits exactly these four bounded shifts, even after arbitrary finite refinement. Their gcd is 2, but integer combinations are not yet certified physical compositions. The two-sided tape move is essential: holding one tape fixed gives no exact positive-row recenter, though an isolated count change of +2 occurs at row 84,608.

## What survives a continuation

[Matching and retention](MATCHING_AND_RETENTION.md) gives a structural reason for the missing coordinate. With one inserted module, an $r$-slot read can demand exactly $r+\nu$ visits, where $\nu$ is the maximum number of vertex-disjoint directed triangles through the port *apart from the port itself*. A cyclic reader therefore cannot universally transport the same truncated linear slot state. In a common triangle wrap, a same-order pair whose two-slot count stays equal needs $H$-gap divisible by 30, or 120 at odd order. The gap-6 pairs cannot pass that test. We also give a finite-depth observability theorem and a literal order-eight pair sharing both $H$ and the two-slot count, separated by one triangle wrap.

[Mirror continuation](MIRROR_CONTINUATION.md) identifies the earlier mirror switch’s whole slot difference: it is twice the slot profile of an actual ordinal tournament. Under every unmarked substitution the count gap is exactly twice that auxiliary tournament’s count in the same reader. The gap stays 2 only for a transitive reader and transitive envelope. For a strong reader it has a sharp larger lower bound, attained by the triangle. Selective contact is therefore the relevant escape to try.

## Analytic contrast

[All-prime twists](ALL_PRIME_TWISTS.md) extends the exact alternating-chain analysis to square-root depth-averaged Fourier decay at every odd prime. At the discriminant prime 113, each chain is exactly equidistributed modulo every power of 113. The chain still realizes only logarithmically many absolute counts up to $X$. This sharpens the distinction between uniform output digits across depth and either branching cancellation at a fixed size or actual coverage at a prescribed integer.

## Which old work mattered

Reproducible draws and precise sources are recorded in the notes. The operative moves came from Walecki repair matchings, Collatz guard and shell ledgers, lonely-runner phase pullbacks, and a p-adic transversality foundry. In each case the useful analogy was a question about what a future operation can still see, not a claim that the two problems are equivalent.

**Next proof target:** combine a source-dependent coarse window with a certified terminal repair wider than the remaining error, or compose compensated readers in a way that keeps the relevant higher-slot data ready. The new pairs close one actual unit interval each; no construction here places such an interval at every admissible integer. Chan–Pak remains open.
