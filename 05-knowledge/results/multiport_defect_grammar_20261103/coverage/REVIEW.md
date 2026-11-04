# Independent review of the growing carrier and strong zero-letter channel

2026-11-03. **REVIEWED mathematical statements and typed implementation;
FINITE-EXACT independent two-birth controls. No blocking error found.**
The current source bytes are pinned by `manifest.json`. This review does not
rerun the transport packet's large existing enumeration or reinterpret its
finite checks as all-order proofs.

## Integral forest carrier

Reviewed [transport/PROOF.md](../transport/PROOF.md) and its
`check_transport.py` implementation. The formal signature count g_p is
correct: the least active name has four unpaired states, or pairs in one
of two directions with one of p-1 other names. The resulting bound is a
coordinate upper bound, not an observer-rank statement.

The edge join is well-defined for distinct active physical names. A missing
name is internal to a spanning component, not an absent vertex. Joining a
component to itself is correctly forbidden. Any finite set of join arcs
either extends the old forest to an acyclic graph with indegree/outdegree
at most one, or cannot be completed in any ordering. This proves commuting
joins, opposite-join annihilation and the support-linear-forest lambda bound.
It also explains why the reverse-oriented flip is I-R rather than applying
the same signed I+R twice. Mixed products remain necessary.

The birth implementation has the correct multiplicity for an anonymous
component. Uniform orientation of all forgotten vertices permits an
anonymous component in only one of the two attachment roles, hence at most
one derivative per birth. It cannot attach both ends of the same component.
Deleting the newborn recovers two separate old components when it was
internal; no shortcut edge is invented. Forgetting pushes an endpoint label
to anonymous and introduces z precisely when both ends become anonymous.
The explicit prohibition on reactivating a forgotten name is necessary.
The two forgotten transitive-pair examples correctly separate its failure.

The rooted E2 readout has no hidden factor2: for a two-component cover the
root-ending component uniquely selects the first ordered slot. The positive
defect identity partitions exactly the covers whose attempted concatenation
fails at r→u. The two/three-slot join formula counts the root either at the
end of the joined path or in the remaining third path, with disjointness
enforced. Coincident-root cases are handled by the singleton/path factors.

For the protected-source star grammar, the start terms vanish because the
protected vertex dominates every other vertex. Peeling w from a w-endpoint
path gives exactly the singleton contributions and one adjacent arm. Thus
the kernel is affine in the w-star bits although the whole H response need
not be. In the actual three-flip example the extra root-inport bit multiplies
that affine kernel; this explains the vanishing third mixed defect term.
The 614 lower bound applies to that fixed Q and protected incidence, with
strongness and other constraints still separately tracked.

The C3 Weyl operator expands to the displayed four coefficients. After m
iterations the only route from initial degree m+1 to degree1 lowers the
degree at every step, giving product k(k-1)=(m+1)!m!. Thus a raw finite
component truncation misses actual future contributions. This is correctly
scoped away from a lower bound on every possible nonlinear encoding.
The subsequent x-star correction has exactly the stated N-then-R arc
orientations and hence the inherited 3H+2b reset. No automatic re-entry
into the b=1/2 input locus follows.

## Strong C2 channel and the two newborn choices

Reviewed [construction/RESULT.md](../construction/RESULT.md) and
`algebra.py`. The constant term of P and the a=w convention Q=z_a are
essential and explicit. The injection proving tau=G-F>=0 distinguishes
whether a lies in the first or second old slot; the two P factors in K
remain ordered even though their polynomials commute.

The C2/C4 path polynomials, transported auxiliary x, old-a and new-y debts
have been checked by a separately implemented forward endpoint DP, versus
the producer's backward subset-polynomial engine. The audit covers138
declared cases, including coincident marks and weak/singleton boundaries.
Selected counts also have direct permutation checks. Both simultaneous
(h,F) tripling and strongness on strong inputs pass.

The alternative y transport remembers ancestor data, not a chosen digit:
its scalar debt is h-e_a and its marked debt is F-J. The latter is at
least e_w>0 on a strong input. The old-a debt has the analogous positive
empty-slot term. The proof of positive e_w by longest-path insertion is
valid and does not require a sampled positivity assumption.

The READY section4 obstruction is correctly dualized with its **strong
input order>=4** hypothesis. The explicit C3 boundary is indispensable:
the H9 image has a2-start-ready y and gives31 at the next C2 step. This
does not contradict the all-order scoped obstruction. Nontrivial strong
inputs have unique zero-start auxiliary x; singleton input is separately
excluded from that uniqueness claim.

The B14 inherited from the1537279 return is independently checked in both
polarities. Its minimum start defect616 and end defect2960 show a real
entry failure at every re-marking of that exact source. This rules out
neither a new source-preserving conversion nor a different digit macro.

## Global consequence remains conditional

The proven strong tripler supplies the zero letter on its accepted chart.
It does not make the six arithmetic seed interfaces or the two positive
digit returns available. The complete six-component scalar induction in
`COVERAGE.md` identifies the missing quantifiers precisely. Source entry,
legal operations and output acceptance remain obligations in addition to
the integer recurrence.
