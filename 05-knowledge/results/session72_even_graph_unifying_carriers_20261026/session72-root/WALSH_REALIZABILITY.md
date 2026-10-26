# Exact spectral coordinates are not probability coordinates

Session72, 2026-10-26. Elementary consequence of the inherited signed Walsh
basis; no priority claim. 'Even' means automorphism-sign-even, not Eulerian.

## 1. The correct positive cone must be transported

Let X be a finite state set of sizeN, and let f_1,…,f_N be any real basis of
functions on X. Write M_{x,j}=f_j(x). A signed state law p has response vector
m=M^T p. The realizable positive-law cone is

    K_M=M^T R_+^X,

not generally the standard coordinate orthant. Normalized laws additionally
satisfy the mass coordinate1 if that observable is included. Invertibility
makes this a simplicial cone with N extreme rays; it loses no information.

An operation with column-acting state-law matrix P transports to

    P_F=M^T P (M^T)^{-1}.

A negative entry of P_F does not mean it sends actual probabilities to negative
probabilities. It means these coefficients are not state masses. Positivity
must be checked on K_M, which the transported operation preserves whenever
P preserves R_+^X.

An invertible linear map carries the full standard orthant onto itself iff
it is a positive monomial matrix: it permutes the extreme rays and scales
each positively. If it preserves normalized probability simplexes, the scale
factors are1. Thus an exact change of observable basis can be a linear
isomorphism without a positive object-by-object interpretation.

This is weaker than the new U2 cycle-cone obstruction. The full Walsh cone is
simplicial and can be returned to point-mass coordinates. A projected terminal-
response cone with more extreme rays than its dimension cannot be made into
an orthant by an invertible change of coordinates of that smaller dimension.

## 2. Four-vertex exact control

Use lexicographic edge order01,02,03,12,13,23 and tournament bit1 for the
increasing arc. Canonical tournament masks are0,2,4,5. Canonical sign-even
graph masks are0,3,7,30: empty, P3 plus an isolate, K1,3, and C4.

For each graph representative g, orient its edges upward and form the signed
orbit sum of Walsh monomials. This fixes an explicit sign convention; reversing
a basis sign is harmless. The resulting matrix, with tournament rows and
sign-even graph columns, is

    M = [1,  4,  0,  1]
        [1,  0,  4, -3]
        [1, -4,  0,  1]
        [1,  0, -4, -3].

For moment coordinates(m0,a,b,c), inversion gives probabilities in row order

    p0=(3m0+c+a)/8,    p1=(m0−c+b)/8,
    p2=(3m0+c−a)/8,    p3=(m0−c−b)/8.

Therefore realizability is exactly the four nonnegative right-hand sides;
normalized laws setm0=1. Equivalently, |a|≤3m0+c and |b|≤m0−c, with the
displayed bounds forcingm0≥0. These labelled inequalities are the missing
positive-law sidecar if graph-indexed moments are mistaken for graph masses.

Full tournament converse exchanges rows2 and4 in one-based indexing and
acts in these coordinates by

    (m0,a,b,c) -> (m0,a,−b,c).

It is a valid deterministic state permutation and preserves the true cone,
but its diagonal matrix has a negative entry. No independent sign choice
for the four Walsh basis elements repairs that diagonal: diagonal sign
changes commute with it. An arbitrary change back to point indicators does
repair positivity, so this is not a universal impossibility result for every
possible basis.

In general converse multiplies graph-indexed Walsh observableΦ_G by
(−1)^{|E(G)|}. The sign-even three-edge star appears first at order4. At
orders1–3, all sign-even graph types have even edge size. Thus order4 is the
first standard-coordinate positivity failure of converse in this basis.
This is full converse, not the distinct literal Hamilton-path complement.

## 3. Connection contract

- Source: tournament type laws and their genuine operations.
- Target: signed sign-even-graph observables with transported multiplication
  and operators, not a probability distribution on graph types by default.
- Map: the invertible Walsh-orbit evaluation matrix.
- Preserved: every invariant response expressible in that basis, including
  its exact operator evolution.
- Not preserved without further structure: positivity in the naive orthant,
  ordinary graph toggles, or a natural equivariant map of underlying objects.
- Repair: retain the transported realization cone and exact operators.

The independent standard-library checker `check_walsh_cone.py` enumerates all
labelled tournaments, graphs and permutations through order4, checks orbit
sign well-definedness, matrix ranks, converse transport, and all sixteen sign
gauges at order4. Normal/optimized receipts agree. Its finite test corroborates
the explicit algebra and does not supply the inherited all-order Walsh theorem.
