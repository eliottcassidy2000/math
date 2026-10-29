# Independent audit: retry language and distinct rational K3 surfaces

## Retry language

**PASS** the source-preserving repairs and the finite failure transducer.
`check_retry_extensions_independent.py` replays 25 pairs of failure charts on
75 original sources, retaining the intermediate integer port, both division
clocks and affine carry. It checks 400 sources on the repaired half-rays,
the exact uniform inequalities, and compressed Mersenne sources for
2<=t<=300. The proof that precisely t=2,4,8 descend on their first compressed
return is valid: LTE supplies c, and the monotonic comparison for t>=10
is sufficient; no infinite-tail extrapolation is used.

The conditional endpoint map on each odd dyadic cylinder has odd unit slope
on its free 2-adic parameter. Hence the Haar product law, failure mass
(1059/2048)^k, and failure-clock means1173/989 and1249/989 are justified.
They are not probabilities attached independently to the deterministic
successive ports of one fixed integer. The stopping mean counts the exit
attempt, whereas those clock means omit the terminal nonfailure word.

The all-period exclusion uses a different deterministic fact. Write each
failure branch as f_i(x)=lambda_i(x-alpha_i)+alpha_i, lambda_i>1.
A composite's fixed point is a positive weighted average of the alpha_i,
with weights (lambda_i-1) times the product of subsequent slopes, divided
by the total slope minus1. Thus every port of a periodic orbit lies between
-16177/139 and -1. The only candidate integer ports are odd[-115,-1].
Exact guarded edges in this finite interval strictly decrease except
-3->-1 and -1->-1. Therefore -1 is the only integer failure cycle.
This does **not** rule out a positive unbounded aperiodic failure orbit.
The checker independently verifies the interval, guarded finite graph,
source-ray inequalities, clock means and chart compositions.

## Three distinct K3 surfaces over Q

`check_k3_independent.py` is a separate standard-library implementation with
manually transcribed named coefficient arrays. It imports no producer module
or SymPy. At p=233 it recomputes the discriminant polynomials, verifies degree24
and squarefreeness, and counts y-solutions on every finite fibre and the
fibre at infinity. Counts are58509,58404,58374 for273,302,4239.

The geometric deduction is valid subject to standard elliptic-surface facts:
coefficient weights2,4,6,8,12 define the O(2) Weierstrass model; a degree24
squarefree discriminant in characteristic>3 gives only nodal I1 finite
fibres, with smooth total space, and nonvanishing discriminant at infinity.
This produces smooth proper good reduction. The corresponding minimal
elliptic surfaces overQ have chi(O)=2 and baseP1, hence are K3.

A Q-isomorphism would identify their l-adic cohomology and Frobenius traces
at this common good prime, forcing equal point counts. Therefore the three
surfaces are pairwise non-isomorphic overQ. Distinct rational forms can
become isomorphic overQbar, so geometric non-isomorphism has **not** been
proved. Nor do these point counts certify Picard or generic Mordell–Weil ranks.

All checks use explicit exceptions and normal/-O outputs agree. No global
integer descent, Sha finiteness, or blanket surface-identity claim is inferred.
