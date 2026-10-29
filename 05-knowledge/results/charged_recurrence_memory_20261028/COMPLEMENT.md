# Actual complements: a closed memory recurrence and its exact scope

**PROVED scoped identities and minima below; FINITE-EXACT ordinary-class
census through order six.** Re-marking is not included in these minima.

## N1. Retaining the path changes the recurrence

Fix the Hamilton path (0,1,...,n-1). A presentation is a tournament containing
these arcs. Its binom(n-1,2) other arcs are arbitrary bits. Native complement
U flips precisely those bits and leaves the path arcs fixed. Hence U^2=I.
It is not global converse. The source below is uniform on presentations,
so an ordinary isomorphism class i has mass f_i=H(T_i)/|Aut(T_i)|.

Let J have class-indicator columns, D=diag(f_i), and K=J^T UJ. Then
K_(i,j) counts presentations in class i whose native partner is in class j;
K is symmetric and its row sums are f. The conditional class transition is

    M=D^(-1)K.

With the original presentation retained, the observed t-step transition is

    F_t=D^(-1)J^T U^t J = I for t even, M for t odd.          (N1)

Generally F_t is not M^t. The latter is the different experiment that
resamples a uniform presentation in the current class after every step.
This distinction was already exposed by
[THM-796, three-sorted incidence](../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md).

## N2. The exact omitted-state recurrence

Let P=JD^(-1)J^T be projection onto class-constant functions and Q=I-P.
In the decomposition im(P) plus im(Q), write U in blocks A,B,C,E, where
A is represented by M. The general first-return elimination gives

    F_0=I,
    F_(t+1)=M F_t + sum_(j=0)^(t-1) L_j F_(t-1-j),
    L_j=B E^j C.                                            (N2a)

Proof: evolve the visible and hidden coordinates separately, starting
with hidden part zero, and substitute the hidden recurrence back into the
visible one. Equivalently,

    sum_(t>=0) F_t z^t = (I-zM-z^2 B(I-zE)^(-1)C)^(-1).

The involution identities A^2+BC=I and AB+BE=0 imply

    L_j=(-1)^j M^j(I-M^2),
    L_(j+1)=-M L_j.                                        (N2b)

Thus every memory term for repeated fixed-path complement is determined
by M, although its powers alone give the wrong retained-path process.
Substitution yields (I+zM)/(1-z^2), agreeing with N1.

This is an application of elementary block elimination, not a claim that
renewal or Schur-complement theory is new. The nearby repo antecedent is
[THM-3305, visible zero mode and rank-two resolvent update](../../../01-canon/theorems/THM-3305-single-edge-visible-zero-mode-and-rank-two-resolvent-update.md):
its selected extra edge changes what a static observer can see. Here the
off-diagonal blocks record what class compression discards.

## N3. Two distinct exact minima

For complement alone with current class observed, the coarsest deterministic
set carrier is the ordered pair (class(x),class(Ux)). All future outputs
alternate these two entries; differing entries are observed at time zero
or one. The transition swaps the pair. Its size is the number of positive
cells of K, with multiplicity retained separately when transporting source
mass.

The smallest U-invariant linear space containing all class indicators is

    span(J,UJ).

Let r be the class count, c the number of support-graph components of K,
and b the number of its bipartite components, counting a loop as
nonbipartite. Then

    dim span(J,UJ)=2r-c-b,
    rank(I-M^2)=r-c-b.                                     (N3)

Proof: the rows of [J,UJ] are (e_i,e_j) for positive cells K_(i,j).
They are the unsigned incidence rows of a bipartite double cover of the
support graph. This cover has c+b components, so its incidence rank is
2r-c-b. Alternatively projecting UJ off span(J) gives the leakage Gram
matrix D-KD^(-1)K and the second rank. The latter rank formula is inherited
from [boundary_support, proposition 0](../boundary_support_20261020_metagraph.md);
that note uses converse-merged classes. The present census keeps ordinary
classes, so its numbers differ.

| vertices | presentations | ordinary classes | leakage rank | minimal linear dimension | minimal set size |
|---|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 0 | 2 | 2 |
| 4 | 8 | 4 | 3 | 7 | 7 |
| 5 | 64 | 12 | 11 | 23 | 50 |
| 6 | 1024 | 56 | 55 | 111 | 732 |

renewal.py constructs every presentation and uses exact isomorphisms,
checks the class count against the independent odd-partition Burnside
formula, and checks f_i against a separate Hamilton-path dynamic program.
It verifies N2 through eight steps modulo a prime and, at n=4, computes
the hidden-state kernels directly over the rationals through six terms.
For N3 a modular minor plus the combinatorial upper bound certifies exact
rational rank. The ordered-pair set minimum needs no matrix rank assumption.

## What changes under re-marking

When a new Hamilton path is selected, a new complement involution acts on
the source tournament. A class-pair carrier for one path need not determine
the new pair. Nor does the scalar matrix M supply legality or transport maps.
For several allowed operations U_a, the linear candidate is the span of
all U_(a_1)...U_(a_s)J. Its minimal observable quotient is a multi-letter
Hankel problem; computing this on actual strong-atom contexts remains OPEN.

The 81/27 hostile pair is rebuilt here: its labelled cyclic triples agree,
while complete automorphism enumeration gives orders 81 and 27. Any proposed
all-context quotient must retain this distinction. No result in this packet
replaces the full ordered-run attachment carrier by cyclic triples, moments,
or the complement-only minimum.
