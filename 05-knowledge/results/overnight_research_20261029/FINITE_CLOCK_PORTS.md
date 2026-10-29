# Finite clocks retain every invariant rational port coordinate

2026-10-29 overnight continuation. **PROVED elementary / FINITE-EXACT
checks**, relative to the stated query language. Same-agent alternative-path
audit. Classical ingredients are finite-group averaging, determinant
multilinearity and polynomial interpolation; no historical novelty is claimed
for those ingredients. This closes the finite-clock, unrestricted integer
coupling subproblem of [MINIMAL_RATIONAL_PORTS](MINIMAL_RATIONAL_PORTS.md).

## F1. The clock and the legal queries

Fix a finite group Gamma of order h and an integral representation
C_gamma on p named ports. A state consists of a nonsingular integral
alternating matrix M, integer port columns P, and source isometries satisfying

    T_gamma^T M T_gamma=M,       T_gamma^T P=P C_gamma.

The transpose is essential: P carries functionals in the update PEP^T.
The representation on the ports and its named basis are known query data.
Allow every integer alternating coupling E satisfying

    C_gamma E C_gamma^T=E       for every gamma,

with update M -> M+PEP^T and an exact determinant observation. Singular
outputs are observed as zero and terminate that branch. Also allow global
sign reversal and unimodular port re-marking P -> PR, transporting every
C_gamma to R^-1 C_gamma R. No new port outside the retained span is added.

For independent port columns the displayed condition is exactly clock
preservation. For dependent columns it is a declared sufficient legality
rule; deciding the larger physically legal set still requires the embedding.
The theorem below applies to this declared language even with dependencies.

Write delta=det M and G=P^T M^-1 P. Then

    C_gamma^T G C_gamma=G.

This follows from T_gamma M^-1 T_gamma^T=M^-1. The update and re-marking
laws of the preceding theorem show that (delta,G) is a sufficient carrier.
The new assertion is necessity even though ordinary pair probes may be illegal.

## F2. Legal orbit sums recover the entire Gram matrix

For i<j let F_ij=e_i e_j^T-e_j e_i^T and form the integer coupling

    E_ij=sum_(gamma in Gamma) C_gamma F_ij C_gamma^T.

Multiplication by a group element permutes this sum, so E_ij is legal.
Invariance of G and cyclicity of trace give

    tr(G E_ij)=h tr(G F_ij)=-2h G_ij.                  (1)

In particular a zero orbit sum forces the corresponding invariant entry
G_ij to be zero. There is no loss hidden by averaging: the unknown being
tested is already invariant, and the entire orbit sum is a legal operation.

Query the polynomial

    f_ij(t)=det(M+t P E_ij P^T)/delta
           =det(I_p+t G E_ij),       t in Z.

It has degree at most p, constant term1 and derivative

    f'_ij(0)=tr(G E_ij)=-2h G_ij.                    (2)

The derivative identity follows directly by choosing the single nonidentity
column in determinant multilinearity. Values at t=0,1,...,p determine the
polynomial, including when an individual queried output is singular.
Explicitly, with H_p=sum_(t=1)^p 1/t,

    f'_ij(0)=sum_(t=1)^p (-1)^(t+1) binom(p,t) f_ij(t)/t - H_p.

Thus legal one-step queries recover every G_ij. Equality of all legal future
responses is equivalent to equality of (delta,G), relative to the fixed named
clock representation. Unimodular re-marking transports the probes and does
not change this equivalence. The same statement holds for a cyclic clock
of any finite order; its order need not be prime or odd.

This is information minimality, not a lower bound on an arbitrary bit
encoding or on the number of adaptive nonlinear queries. The shared empty
query plus p queries per pair gives an explicit finite reconstruction bound.
The degree can be reduced to rank(E_ij) when computing fewer values.

## F3. Exact coordinate count and the C3 example

Let r be the dimension over Q of the invariant alternating-form space.
The trace pairing between invariant G and invariant E is nondegenerate by
(1). Choose any rational basis of the E-space and clear denominators.
The resulting r legal directions and their first determinant derivatives
recover r rational coordinates of G. All remaining entries are determined
by invariance. A finite clock removes redundant coordinates, but it does
not identify two distinct invariant Gram matrices in this query language.

The standard exterior-square trace identity yields

    r=(1/(2h)) sum_gamma [(tr C_gamma)^2-tr(C_gamma^2)].   (3)

For completeness: the eigenvalues on the exterior square are lambda_i lambda_j
for i<j, giving the half-difference in (3); averaging the group action is an
idempotent onto its invariant subspace, whose trace is its dimension.

For a cyclic permutation whose cycle lengths are l_1,...,l_k, equivalently

    r=sum_i floor((l_i-1)/2) + sum_(i<j) gcd(l_i,l_j).    (4)

Within a cycle, unordered pairs have distance classes; an antipodal class
in an even cycle reverses orientation and contributes zero. Between two
cycles there are gcd(l_i,l_j) pair orbits and no orientation reversal.
This proves (4) without character theory.

For cycle type(1,3), r=2. In the ordered basis with fixed port0, every G is

    [ 0  a  a  a ]
    [-a  0  b -b ]
    [-a -b  0  b ]
    [-a  b -b  0 ].

The two legal orbit probes are the star F01+F02+F03 and the oriented triangle
F12+F23+F31. Each has rank2, and its determinant ratio is respectively
(1-3at)^2 or (1-3bt)^2. Querying t=+1,-1 therefore recovers a and b in
four nonempty queries: each coordinate is the negative-minus-positive
determinant difference divided by12delta. This is a complete two-coordinate
carrier for this rational query language, not a native sign-cut classification.

For cycle type(1,3,1,3), r=10. For a single odd n-cycle r=(n-1)/2;
adding a fixed port gives (n+1)/2. At n=25,27,57 these latter dimensions
are13,14,29. These are rational response dimensions, distinct from the
binary invariant-Lagrangian counts, which vanish at all three lengths.

## F4. Infinite-order clocks can lose an entire coordinate

Finite order cannot simply be omitted. Let U=[[1,1],[0,1]], C=diag(U,U),
and for any integer d put

    G_d=[ 0  1  0  0 ]       M_d=G_d^-1=[ 0 -1  d  0 ]
        [-1  0  0  d ]                    [ 1  0  0  0 ]
        [ 0  0  0  1 ]                    [-d  0  0 -1 ]
        [ 0 -d -1  0 ].                   [ 0  0  1  0 ].

Take P=I and T=C^T. These are integral states with det M_d=1,
C^T G_d C=G_d and T^T M_d T=M_d. Directly solving CEC^T=E gives

    E=[ aJ  Y ];     Y=[[u,v],[-v,0]],       a,b,u,v in Z,
      [-Y^T bJ]

where J=[[0,1],[-1,0]]. For either epsilon=+1 or-1,

    det(epsilon M_d+E)=[(a-epsilon)(b-epsilon)-v^2]^2.   (5)

The entire expression is independent of d. Any finite allowed word,
including transported unimodular re-markings and sign reversals, pulls
back to epsilon M_d+E_total with an invariant E_total. Thus all the states
have identical responses, including singularity decisions at every prefix,
although their G_d differ. Nonlinear determinant queries do not rescue
this coordinate. The smallest example used here is four-dimensional;
no universal minimal-dimension claim is made.

The proof is the displayed formula, not a bounded census. The checker also
tests all a,b,u,v in[-2,2] and d in[-2,2] with both signs.

## Inheritance, hostile boundaries and the next experiment

Anchor: finite-clock minimality. Niche: rational reconstruction under an
infinite unipotent clock. Wildcard: arithmetic joint divisibility invisible
to a chosen modulus. The board is invariant form / legal probe / source
embedding / clock / integral saturation / operation cost.

Closest mechanism: [the preceding minimal response theorem](MINIMAL_RATIONAL_PORTS.md).
Corrected near miss: dual ports transform by T^T P, not TP. Canonical
hostiles: the different integral groups H9 and H3+H3 with equal responses,
and the separate tournament81/27 pair with different attachment orbit counts.
Neither is removed by this rational theorem.

The least-used recovered sidecar is
[THM-2953, cyclic charged corank and Reynolds boundary](../../../01-canon/theorems/THM-2953-cyclic-charged-corank-parity-and-reynolds-averaging-boundary.md).
There, averaging a rank-one observation can manufacture joint injectivity
without a lawful single-source response. Here the unknown is already
invariant and the full integer orbit sum is explicitly permitted. Those
two hypotheses are exactly what pay that earlier physicality debt.
The trace-duality step also echoes
[THM-3280, cyclic differential-module endomorphisms](../../../01-canon/theorems/THM-3280-fc-cyclic-module-simplicity-and-ordinary-basis-repair-through-degree-nine.md),
without transferring its differential-module conclusion.

The META-PATTERNS cards used are "Demodulate, then average a complete
subgroup" and "Audit saturation and basis covariance before naming a
lattice bridge". No new card is promoted from one example.

Source -> map -> target: a clocked alternating form maps to determinant
responses through legal orbit couplings; it preserves exactly delta and G.
It destroys the integral cokernel, ambient lattice embedding and unrelated
tournament transporter observations. Their sidecars remain mandatory.
The cheap decisive probes are the explicit pair orbit sums and the
unipotent family (5).

The next test is the genuinely narrower operation set: cross-block-only
or bounded/native sign couplings. Orbit sums may leave it, so F2 is not
automatically applicable. Determine which invariant Gram coordinates those
actual future words separate. General native-complement closure, arithmetic
Selmer realization and Sha finiteness remain OPEN.

## Reproduction

Run `python -B finite_clock_ports.py` and the same command with `-O`.
The standard-library checker uses direct determinant interpolation, an
independent matching Pfaffian for skew determinants, finite-group invariant
linear equations, cycle pair orbits, and the infinite-clock hostile.
The companion JSON fixes the universes and counts. This is not a Lean or
independent-agent proof audit.
