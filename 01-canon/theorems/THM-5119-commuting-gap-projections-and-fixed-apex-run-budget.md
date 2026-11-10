---
id: THM-5119
title: "Commuting gap projections, square-free tape-letter spectra and fixed-apex run budgets"
status: "PROVED, independently audited. Deductions from the inherited tape compiler and deletion kernel; no coverage theorem."
session: Session245, 2026-11-10
source: 05-knowledge/results/session245_physical_spectrum_20261110/forest/PROOF.md
audit: 05-knowledge/results/session245_physical_spectrum_20261110/audit/AUDIT.md
---

# Every constant tape run has a square-free integer spectrum

Fix a labelled tournament on k apexes; other vertices form a transitive
chain with arbitrary apex-neighborhood letters. Use the inherited
root-permutation/gap automaton of Sessions215/216, restated in Session221.
For one letter a and one gap g, its assignment operator E_(a,g) either
places the next vertex in g or rejects the assignment as the zero vector.
It reads only the occupancy/status of gap g and sets that gap to a fixed
letter-dependent occupied status. Therefore

    E_(a,g)^2=E_(a,g),
    E_(a,g)E_(a,h)=E_(a,h)E_(a,g).

The letter operator A_a is their sum, including multiplicities. Simultaneous
kernel/image splitting over Q proves

    product_(j=0)^(k+1)(A_a-jI)=0.

This applies to arbitrary apex orientations, all neighborhood letters and
all formal gap states. Different letters are NOT asserted to commute.
For any fixed physical prefix u and suffix v,

    H(u a^m v)=sum_(j=0)^(k+1)c_j j^m.

The coefficients may be rational or negative; the count is an integer. For
m>=1 only bases 1,...,k+1 remain. At m=0 the zero-eigenvalue correction must
be retained. Multiple runs give an ordered multivariate exponential
polynomial, not a permission to exchange runs.

For k transitive initial roots and leaf tape (-1)^m 0,1,...,k-1, where -1
loses to every root and i beats just root i, the coefficient of (k+1)^m is
exactly k!. The top projection fills every gap; each final leaf must then
repair the gap preceding its matching root, uniquely for each root
permutation. The two-star scalar example is inherited THM-4748, not new.

## A quantitative limitation of few runs

Fix also an alphabet of q>=2 letters. Let C(N,r) be the number of distinct
counts <=N represented with at most r constant runs. Session221's inherited
finished-word SCC deletion preserves H and leaves a subsequence of length

    L(N)=max{2a+3b:3^a5^b<=N} <= floor(3 log_5 N).

A subsequence cannot increase run count or introduce a letter. Counting
positive run lengths yields

    C(N,r) <= 1+sum_(j=1)^min(r,L(N))
                    q(q-1)^(j-1) binom(L(N),j).

Consequently fixed r gives O((log N)^r) values, and r=o(log N) gives
N^{o(1)} values. A finite set of fixed cores has the same conclusion.
Growing apex sets, growing alphabets or unrestricted run budgets are not
excluded. This is a finished-word count, never a continuation-state quotient.

Independent literal controls include all apex orientations through k=3 and
all length-two full-alphabet words, plus 2730 fixed-apex kernel tests. Finite
censuses remain explicitly restricted. Chan–Pak stays OPEN.

[Projection proof](../../05-knowledge/results/session245_physical_spectrum_20261110/forest/PROOF.md) ·
[Run-budget proof](../../05-knowledge/results/session245_physical_spectrum_20261110/forest/RUN_BUDGET.md).
