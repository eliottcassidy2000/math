# When residue loops certify a cutoff

Let S be the set of q=(H(T)-1)/2 for actual finite tournaments. Translation
preserves parity as the relation x-y in 2Z, but an odd/even *label* requires
an origin. The ternary relation x+y=z is covariant under translating x,y
by a and z by 2a, not invariant under simultaneous translation. Its
midpoint cousin x+y=2z is invariant under simultaneous translation.

## 1. Residues with size control

Suppose for a target q and modulus m a physical construction gives
s in S with s=q mod m and |s-q|<m. Then s=q. The size requirement
cannot be replaced by a profinite limiting argument. In particular, if
for every q>N there is any sequence of moduli m_j tending to infinity
with physical lifts s_j=q mod m_j and |s_j-q|/m_j tending to 0, S is
cofinite. Even the much weaker requirement that this ratio is eventually
less than one suffices. No such uniform construction is claimed.

A sharper negative control than sparse artificial sets: the set of odd
composite integers is multiplicatively closed and its closure in the
profinite integers is the whole odd coset. For an allowed residue a
modulo m, choose an odd prime p not dividing m and use CRT to impose
x=a mod m, x=0 mod p, x odd, then take arbitrarily large x. It also
realizes every even difference d: take x odd, x=0 mod3 and x+d=0 mod5,
and make both large. Still it omits all odd primes. Adjoining the singleton value 1 preserves
all three properties and all those omissions. This jointly controls
three known kinds of Chan–Pak evidence, not just residue coverage.

The physical control is stronger still: Session199
`profinite-source-escape/PROOF.md` gives a congruence-dense actual C2
construction with a missing physical 67 and a phantom boundary profile.
Its bounded-height recovery is the same distinction between physical
representative and finite-wheel shadows; this section does not claim
discovery of that repair.

## 2. What a finite machine would have to remember

Fix base b>=2, use canonical most-significant-first digit words, and
consider an exact finite-state *physical* language L of covered q. This
can be a complete language of all realizable q, or merely a certified
subfamily if the goal is a sufficient coverage certificate. Determinize
if needed, complete with a dead state, and intersect the complement with
the regular language of canonical positive digit words. Trim to states
reachable from the initial state and co-reachable to a rejecting terminal
of the original covered language. The missing positive q form a finite
set iff this trimmed directed graph is acyclic. A longest path to a
terminal bounds the digit length; choosing the lexicographically largest
missing word among words of maximal length yields the exact largest
missing q (assuming the set is nonempty). The zero case is separate.

A directed cycle gives words uv^j w with nonempty v, and therefore an
explicit infinite family. If [z] is the integer value of a digit word,
B=b^|v|, C=b^|w|, then its integers are

    C([u]B^j + [v](B^j-1)/(B-1))+[w].

This proves a cutoff precisely at the loss of a productive failure loop,
not by estimating a long prefix. It does not assert that a fixed modulus
observer (Session295), or a finite dyadic bank (Session296), is an exact
physical-language recognizer. Session197's automatic-selector sparsity
rules out a particular fixed-final-gadget route to such a certificate.

## 3. A distinct two-clock sufficient theorem

If the SAME physical S is recognizable in two multiplicatively
independent integer bases, Cobham's theorem makes S eventually periodic.
If S has arbitrarily large members in each residue modulo every modulus,
the periodic tail must then include every residue of its period, hence
S is cofinite. Congruence density for the actual spectrum follows from
the inherited physical constructions (e.g. Session199's C2 family);
one still needs BOTH automaticity theorems for the SAME S. Different
base-dependent path spectra do not satisfy this premise.

One clock and profinite density are not enough, even jointly: the
nonnegative integers other than powers of b form a b-automatic,
congruence-dense set with infinitely many holes. To see density, each
arithmetic progression has linear counting growth, while the excluded
powers have logarithmic counting growth. Even an eventual threshold,
if proved, leaves a finite range for the existing 7 and 21 exclusions
and finite witnesses to settle; it is not itself that classification.
