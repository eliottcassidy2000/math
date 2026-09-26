# Independent audit of growing-prefix multiplicative rank

**Status:** PROVED audit, with the S-unit input CITED; FINITE-EXACT controls.
Date: 2026-09-26. Reviewed source:
[crossroads223_20260926_geometry_growth.md](crossroads223_20260926_geometry_growth.md).
No mathematical gap was found in its main finite bound, factorial
asymptotic, odd-source normalization, or uniform-interval conclusion.
The no-dip refinement below includes the required source-size condition.

## 1. Local unit groups and determinant identity

For an exact valuation word `(k_1,...,k_N)` with total S_N, the extra
terminal oddness condition requires modulus `2^(S_N+1)`. A source can be
constructed independently from

    n = 3^(-N)(2^S_N-C_N) mod 2^(S_N+1).

Adding only `2^S_N` changes the final valuation, a useful off-by-one hostile.
With this full modulus, the affine nodes have the stated slopes and
integer intercepts. Expanding the two affine forms gives exactly

    A_j B_i-A_i B_j
      = -3^i 2^(S_N+1-S_j) C_(i,j).

The carry is nonzero and coprime to 6. The roots are distinct, although
particular sources can still make two node values equal; those collisions
remain in the exceptional set and require no separate exclusion.

A prime outside a node's incident carry bank cannot divide another node.
Its valuation therefore annihilates that node's coefficient in a proposed
first-slot relation. Every supported node is consequently a unit in its
own bank. At least two nodes are supported because all `3m_i>=3`.

For a supported pair, `x=A_j L_i/D_ij` and `y=-A_i L_j/D_ij` are nonzero,
sum to one, and belong to their respective rational unit groups. Their
product group has rank `s_i+s_j`, including when the prime banks overlap:
the generators `(p,1)` and `(1,p)` are independent. Signs contribute
torsion, not an extra free rank. The parameter maps injectively to x.

The author-hosted first page of
[Beukers--Schlickewei, Theorem 1.1](https://webspace.science.uu.nl/~beuke106/s-units.pdf)
was checked directly. It bounds the number of solutions in the Q-closure
of a rank-r subgroup of `(C*)^2` by `2^(8r+8)`. Therefore the note's
bound `2^(16a_N+40)` is correct after `s_i,s_j<=a_N+2`.

## 2. Factorial inverse and interval normalization

The incident valuation-span budget is at most `(N-1)S_N`, and the
incident length budget is at most `N(N-1)/2`. At `S_N<=3N`, these give
the strict integer bound `Q_i<96^(N(N-1)/2)`. If t distinct primes divide
Q_i, their product is at least t!, including the case t=0. Thus the
strict factorial inverse used by the note has the right inequality.

Writing alpha=log_2 3, Stirling's elementary logarithmic estimate gives

    a_N ~ ((5+alpha)/4) N^2/log_2 N,
    log_2 E_N ~ 4(5+alpha) N^2/log_2 N.

At `N=floor(sqrt(Q log_2 Q)/8)`, `Q=log_2 H`, one has
`N^2/log_2 N ~ Q/32`, yielding exponent `(5+alpha)/8<1`.
The bounded-word exception count is global in the source height, so
substitution into an arbitrary interval costs no location-dependent term.

For an odd source, the odd returns after time zero occur at exactly
times `S_1,S_2,...`. Thus `S_N>L` iff the parities at times 1,...,L have
at most N-1 odd terms. The initial odd bit is additional: the appropriate
Terras modulus is `2^(L+1)`, and the number of classes is
`R=sum_(j<N) binom(L,j)`. Among all integers, this has density R/2^(L+1),
not R/2^L. The term `H/2` in the note is correct. Each class contributes
at most `H/2^(L+1)+1` in any interval, giving the stated uniform remainder.

Full first-slot independence implies paired independence, so paired-rank
failure is a subset of the counted first-slot failures. The direction
of the final consequence is correct. For fixed N, the upper Banach
density statement follows by taking H to infinity before increasing the
valuation truncation L; interchanging these limits is unnecessary.

## 3. No-dip refinement and its precise quantifiers

Let `ell=ceil(alpha N)+1`. Suppose `n>=2^ell` is odd and has no strict
descent through time `floor(log_2 n)`. Then

    S_N<=ell.                                         (A1)

**Proof.** If `S_N>ell`, at most N odd steps have occurred by time ell.
Every encountered odd state is at least n. The exact product of relative
step multipliers therefore gives

    T^ell(n)/n <= 3^N 2^(-ell) (1+1/(3n))^N < 1.

Here `2^ell>2*3^N`, and `n>=2^ell>=N` makes the last product at most
`exp(1/3)<2`. Since ell is within the assumed no-dip horizon, this is a
contradiction. QED. This avoids assuming in advance that the Nth odd
return occurs within the no-dip horizon.

The source-size condition cannot be discarded: `n=1,N=5` never descends
but has `S_N=10>ceil(alpha N)+1=9`.

For words satisfying (A1), use the exact integer budget

    P'_N=2^((N-1)ell-N(N-1)/2) 3^(N(N-1)/2),
    a'_N=max{a:a!<P'_N},
    E'_N=binom(ell,N) binom(N,2) 2^(16a'_N+40).

The same local-bank proof applies to every such word and gives at most
E'_N exceptional integer sources over all heights. Now

    log_2 P'_N=((3alpha-1)/2)N^2+O(N),
    log_2 E'_N ~ 4(3alpha-1)N^2/log_2 N.

Consequently, with the same N(H), the number of odd n in *any* interval
I of H consecutive positive integers that both have no descent through
`floor(log_2 n)` and fail the N(H)-node first-slot rank test is at most

    2^ell+E'_N = H^((3alpha-1)/8+o(1)).                (A2)

The omitted small sources are bounded by `2^ell=H^o(1)` uniformly in
the location of I. In particular (A2) applies to n<=X with H=X.
This is an upper count of failures inside a conditioned set. It is not
a lower count, does not say that all no-dip sources have full rank, and
does not alone give relative density one inside that set in every
translated interval. A separate lower bound on the conditioned set's
size would be needed for that relative-density statement.

## 4. Independent exact controls

The standard-library probe checks 594 complete valuation cylinders of
lengths 2,3,4 and total at most 3N; 1,782 positive integer realizations;
9,711 pair identities; and 19 source-collision controls. Its affine-source
construction is independent of the geometry script. Four full odd-residue
universes verify the L+1-bit normalization, and 36 intervals with different
locations verify the residue-count bound. A direct census below 65,536
finds 2,903 no-dip odd sources and checks (A1) in 18,686 eligible cases.
The small-source hostile is checked explicitly.

    python 04-computation/experiments/crossroads223_20260926_geometry_growth_audit.py

All arithmetic gates remain active under `python -O`; the adjacent output
is deterministic. These finite controls do not evaluate the enormous
asymptotic exception bound at a practical threshold and do not replace
the cited S-unit theorem or the quantifier argument.
