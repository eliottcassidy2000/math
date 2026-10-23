# Independent native branch-width audit

Verdict: PASS, scoped to the stated contiguous cyclic-cut Hamilton-path family. Audit date 2026-10-23. Source manifest SHA-256: `8676a9f9e6823aa756a9c3be845f2e916991d85fd99a516357b0ca71508a9c7a`.

## General proof review

The host is a genuine tournament: each internal cell is a regular odd circulant tournament; each outer oriented cell pair assigns one direction to every cross pair, reversing all other pairs. For a fixed vertex, an outer-outgoing cell contributes |E| outgoing arcs, and an outer-incoming cell contributes p-|E|. Equal numbers of each give total degree (pt-1)/2. The construction requires odd t>=3 and odd prime p as stated.

The cyclic-cut path in cell i starts at a_i and ends at a_i-1 modulo p. Thus its forward boundary to the next cell has cross difference a_(i+1)-a_i+1. Since D=E+1, the boundary is valid exactly when the phase difference lies in E. Distinct phase tuples give distinct literal paths, since their first vertex in each labelled cell identifies the phase. No assertion about all Hamilton paths is needed.

The linear difference map F_p^t -> F_p^(t-1) is surjective with the one-dimensional diagonal kernel. Every affine coset in an exact union representation must itself be contained in the target set. Its image is an affine subspace in E^(t-1). If the image has a nonzero direction, some coordinate projection is the whole F_p, contradicting properness of E. Therefore all contained affine cosets have size at most p. Cardinality gives the cover lower bound |E|^(t-1) even when overlaps are allowed. Full diagonal fibres attain it. This proof includes arbitrary additive subgroup cosets because all subgroups of the elementary abelian group F_p^t are F_p-linear subspaces. Primality is essential to this argument; the source appropriately withholds an unconditional composite-modulus extension.

For p=3,E={0,1}, the transfer I+Q retains the actual common phase and counts shared owners without expanding all branches. The endpoint Fourier formula agrees with binomial counts modulo 3. The sign r-2delta in the cosine is correct. Each A_r(delta) fixes the initial owner as well as the phase difference; removing that owner pin multiplies by 3. Inserting diagonal masks retains specified intermediate owners. Extra guards require sufficient state and are not automatically bounded-memory or additive.

The integer period-six correction applies to all r>=0; regular odd-t hosts use only even r=t-1. This is compatible rather than contradictory. Converse preserves accepted literal paths after reversing both the path and the host arcs, without any assumption that the converse family has the same named phase coordinates.

## Independent evidence

The checker performs 50,919 explicit RuntimeError-guarded checks. It constructs actual arc predicates independently, verifies every pair and every score, enumerates phase tuples and their literal paths in six hosts, and compares endpoint distributions against integer polynomial coefficient convolution. It also checks reversed paths in converse hosts.

For F3^3 and F5^3, all affine subspaces are independently constructed as points, parametrized lines, solution planes of linear equations, and the full space, rather than by the producer's subgroup-extension algorithm. All 6 and 30 nonempty proper subsets E respectively are tested. Maximum contained affine-coset size is p, and diagonal fibres attain the cardinality bound. The counts of distinct affine cosets are 184 and 1056. The period-six expression is checked against exact binomial sums for r=0,...,100.

Normal and optimized runs are byte-identical. A temporary mock publication layout with sibling `audit-branch-width/` and `native-branch-width/` was explicitly tested in both modes, with the same receipts. The checker also supports the original workspace paths. Every source-manifest file hash and size is checked before mathematical tests.

## Limits

This is a lower bound for explicit affine-coset union representations of this specified family, not for arbitrary algorithms, circuits, Hamilton-path counting, or all native metagraph wiring. The p-state transfer is a positive compact alternative, so there is no computational-hardness conclusion. General claims follow from the proofs, not from the finite regressions. No source packet, earlier Gaussian audit, or Git files were changed by this auditor.
