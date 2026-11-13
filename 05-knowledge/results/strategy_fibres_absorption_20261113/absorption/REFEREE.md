# Cross-review of growth and singular-word reductions

Status: **ACCEPTED within the stated scopes.** Read-only audit of
[GROWTH.md](../GROWTH.md), [orbits/SINGULAR.md](../orbits/SINGULAR.md), the
inherited core-reduction theorem, the literal thirteen-generator file, and
the associated proof/receipt files. No new search or Git operation was used.
No mathematical correction to G1--G3 or the singular pruning theorem was
required. This is a proof audit, not an independent full rerun of their jobs.

## Growth statements

**G1 accepted.** A newborn increases H exactly when it is in the child's
cyclic core, by the inherited induced deletion kernel. Each old cyclic vertex
retains an old directed cycle under births. Therefore positive births inject
distinct persistent vertices into the final cyclic core. Applying the inherited
K(H) bound proves c(T0)+k<=K(H(Tm)). Arbitrarily many neutral births do not
alter this argument. The resulting exclusion concerns one nested branch with
subexponential count in its positive-step index; it neither excludes branching
coverage nor applies after an old arc is changed.

**G2 accepted.** The independent-set partition function at activity two sends
Z to 3Z for an isolated added vertex and to Z+2 for a universal added vertex.
The backward odd-integer procedure strictly decreases, with at most two
subtractions before a division; the displayed logarithmic size bound is
valid. This is a construction of conflict graphs, not of tournament cycle
inventories. The obstruction at 7 and 21 and the stronger nested-birth
obstruction are stated at their proper transfer boundaries.

**G3 accepted for the ordinary directed cycles n>=3.** Retaining the cycle
arcs puts its n distinct vertices in a single SCC, also when extra vertices
are allowed. The inherited strong minima satisfy f(s)>=5^((s-1)/3) in all
three residue classes, and SCC factors are positive. Thus the exponential
completion lower bound follows. It is a restriction on preserving the given
unary cycle geometry, not on all ways to compile an integer into a tournament.

## Singular-word pruning

**Accepted at all word lengths for the listed patterns.** Construction-order
reversal is handled correctly: the word Xa,C,N corresponds to N C Xa.
The displayed matrix identities are consistent with the literal generator
file. Their outer-product factors carry actual observers and actual initial
rays, which is essential:

* N squared has actual return H(Ns)>=5 and actual output ray N*s0 of H=5.
* N C Xa has actual return H(C Xa s)>=5 and the same actual H=5 ray.
* Xb Xa has actual return H(C Xa s)>=5 and actual ray (1,b,b).
  Every R-free generator's first row dominates (1,2,0), so any nonempty
  continuation of that ray has H>=1+2b>=3.

For p or 7p, a product of two nonunit attained odd counts is impossible:
for a prime this is immediate, and for 7p a proper nonunit factorization
requires an attained factor 7, including p=7. Thus NN and XaCN are excluded,
and XX can occur only at the end. Replacing the final Xb by C preserves the
first row exactly. It cannot introduce another prohibited pattern because
the replacement is terminal. This justifies existence of a reduced
representation, not equality of rooted states in subsequent contexts.

The classification of *all minimal rank-one patterns* is only through length
three, as explicitly stated; the proof makes no extrapolation beyond that
finite classification. The gcd-seven example correctly blocks the false
inference that a coordinate divisor must be an attained tournament factor.
The strong H=25 reset also correctly shows that arithmetic factorization need
not be an ordinal-sum graph decomposition.

## Calibration transferred back to absorption

The core budget makes an existential O((log M)^2) batch transformation between
chosen witnesses for already attained adjacent counts <=M automatic: remove
singleton SCCs, pad to K(M), and reverse all discrepant arcs at once. Therefore
a useful new locality statement must impose a uniform small support, preserve
a needed boundary, or construct a repair without assuming a target witness.
This warning has been added to [PROOF.md](PROOF.md). The measured finite radii
R(6)=2 and R(7)=3 remain nontrivial exact small-support facts, but no coverage
or uniform-radius theorem follows from them.

The terminal-absorption certificate A1 retains a different useful obligation:
control the approximation gap and the actual available residual correction
on the same approximate witness. It does not silently reuse a completed
target witness or require an indefinitely nested additive chain.
