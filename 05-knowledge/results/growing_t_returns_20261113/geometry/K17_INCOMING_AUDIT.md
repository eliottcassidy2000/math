# Independent review of the complete k17 marked exclusion

**Accepted: PROVED reduction plus incoming FINITE-EXACT exhaustive
certificate.** The target53084231=405*2^17+71 has no indegree-one unit-ready
root at any tournament order. Therefore an all-k realization of that
exact405-ray interface is **REFUTED**. The ordinary Hamilton count is
attained; other defects and Chan--Pak remain outside this exclusion.

Reviewed sources are
[TARGET_INVERSE](../../target_directed_cores_20261112/TARGET_INVERSE.md),
[COMPLETE_INVERSE_AUDIT](../../target_directed_cores_20261112/ports/COMPLETE_INVERSE_AUDIT.md),
both C++ enumerators and the binary inverse. The current source hashes
and receipt-layer comparisons are in `K17_INCOMING_AUDIT.json`.

For a canonical unique-start core, the target observer is

    Q(x)=h g(x)+(s+p)x+e parity(x)+s+p-e-f-1.

Since g>=1 and all variable coefficients are nonnegative, Q(x)>=Q(0).
The exact identity Q(0)-(p-f)=h-e+s-1>=0 and bin exchange f<=p/2 give
Q(x)>=p/2. The2^c assignments of c vertices to two path bins each have
at least one path pair, so p>=2^c. Thus p<=106168462 and c<=26 are
complete intrinsic bounds, not guessed strong-order cutoffs.

The recursive pruning is also safe. Inserting one new vertex into either
chosen bin gives2p distinct enlarged covers; deleting that vertex recovers
both the old cover and chosen bin. No core ancestor above the budget can
return below it. The stronger legal-tail bound follows from

    p'=2(L3*S3+B2*p-D2*f),
    S3>=p/2, f<=p/2,
    L3>=2^l+3^(l-1), B2>=2^(l-1), D2<=2^l-1.

For S3>=p/2, cut the root-containing path of a two-bin cover just before
its root and retain its original-bin bit. This injects those covers into
root-start three-bin covers times a two-element sidecar. For a nonempty
legal word, its last1 supplies the stated L3 term; forcing that last
symbol to the before position gives the B2 lower bound. Hence
p'>=(2^l+3^(l-1)+1)p. Empty tails have p'>=2p. These lower bounds are
monotone in tail length, validating both the producer's skips and the
auditor's early break.

Both programs enumerate the singleton and every recursively generated
strong unique-start core, retaining full states and chronological
addresses. Neither deduplicates by the scalar target. Singleton plus
empty tail is correctly excluded as a strong constructor. The complete
unique-start grammar is an inherited proved dependency. The weak-source
reduction to strong ready core followed by a transitive tail supplies
the all-order scope; no additional target values come from weak padding.

The producer's even/odd inverse pullbacks and necessary address intervals
are algebraically exact and do not assume g is monotone. The independent
auditor instead visits every outer integer up to its valid elementary
bound. The two saved receipts agree on every core layer, with5640026
retained cores and no hits; the direct auditor tests74845805 outer
addresses. Root additionally recompiled and reran that independent C++
auditor during this session and recovered the same totals and zero hits.
This geometry review does not claim a third large enumeration.

The local checker independently compares those entire layer vectors and
tests the cover inequalities and vertex-insertion floor on every labelled
nonempty tournament through order4 and every attachment. It tests all
legal final1 words through length10 for the coefficient floors. These
small controls support the reviewed reductions; they do not substitute
for the incoming exhaustive certificate.

The result changes the research target: a universal unit-ready405 return
cannot be completed in this interface. The finite physical three-return
chain and the broader nonunit constructions remain valid. No existing
finite positive certificate is retracted.

Reproduce `python -B .../geometry/k17_incoming_audit.py`, then with `-O`.
