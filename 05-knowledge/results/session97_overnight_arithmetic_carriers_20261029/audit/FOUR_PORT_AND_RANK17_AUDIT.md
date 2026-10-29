# Independent audits: joint-port threshold and rank17 height certificate

## Eight marked length32 codes

`check_four_ports_independent.py` reconstructs the eight codes from their
explicit generator masks without using the producer's factorization, rank,
or code constructors. It forms each 32768-word even subcode and its extended
code, enumerates all620 octads, and checks every point/pair/triple response
is155/35/7. It recovers both four-port histograms exactly and verifies that
all eight full weighted four-port tables differ. The named set {0,1,2,3}
is a zero-response witness precisely in the two-code family.

The multiplier orbits are independently identified by direct codeword
membership: sizes2 and6. Multipliers supply coordinate equivalences within
each orbit; the different four-port histograms obstruct arbitrary coordinate
equivalence between them. Thus there are exactly two S32-equivalence classes
*within this eight-code bank*, not eight. Full order4 tables identify the
marked source; their histogram identifies only its two-class quotient.

The reversal theorem is also sound: reciprocal selection is swapped in every
pair, so reflection has no fixed kernel. In the positive-count case, a unit
acting trivially on every Boolean choice must preserve each irreducible
factor. In particular it preserves the doubling orbit of1, and hence belongs
to <2>. This proves the asserted faithful quotient action. No dihedral-stable
kernel or hidden distinguished zero is manufactured.

## Classic rank17 elliptic K3 certificate

`check_rank17_independent.py` is standard-library-only and uses the retained
primary coefficient arrays as inputs. It reconstructs each polynomial square
root by descending coefficient recursion, rather than SymPy factorization;
checks all seventeen section equations; computes all136 off-diagonal height
pairings from rational-polynomial gcd degrees and the infinity correction;
and matches the complete published matrix. Rational Gaussian elimination
finds determinant948 and positive leading principal minors. All fifteen
named secant identities are independently verified.

The height formula requires the stated smooth elliptic K3 model with no
reducible fibres and integral sections disjoint from O. Those are the proper
geometric guards: raw gcd counts without the model and infinity correction
would not be a height proof. The certificate establishes independence of the
seventeen generic sections of this specific determinant948 fibration. It
must not be transferred to the different named high-specialization surfaces.

Replacing P_i by -P_1-P_i for a chosen subset of the fifteen indices is an
integral unimodular basis transformation. Its determinant is a sign and the
Gram transforms by congruence, preserving determinant and positivity, while
named pairings can change. This is a precise rank-preserving loss of named
information, not a new rank record or an identification of unrelated fifteens.

Both independent checks reproduce identically normally and with `python -O`.
