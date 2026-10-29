# Parameter conjugacy with a monotone integer-source sidecar

**Status: PROVED elementary transfer and obstruction theorems; OPEN convergence.**
No claim of literature novelty for valuation coding. Source inheritance is
Session99 `structural/GEOMETRIC_SLOT_BRIDGE.md` (identical valuation laws do
not fix native integer dynamics), `transfer/CARRY_CLOCK_POLES.md` (retain
carry and source), and Session97 `controllers/PERIODIC_EXCEPTION_CERTIFICATE.md`
(periodic words have unique rational sources). The new-to-this-packet use is
an explicit parameter-transfer domain, an integer-source height sidecar, and
an all-preperiod denominator certificate against the named 5x+1 hostile.

## 1. Exact coding, map, inverse, and lost predicate

Fix an odd integer b>=3. Let X_b be the odd 2-adic sources whose accelerated
orbit never hits -1/b, where the next valuation would be infinite. On X_b set
T_b(x)=(bx+1)/2^a, a=v_2(bx+1)>=1. Let W={1,2,...}^N, with left shift sigma.
For w=(a_0,a_1,...), S_0=0 and S_j=sum_{i<j}a_i, define

    Phi_b(w) = -sum_{j>=0} 2^{S_j}/b^{j+1}   in Z_2.             (1)

The series converges since S_j>=j. Its value is odd. Splitting its first
term proves b Phi_b(w)+1=2^{a_0} Phi_b(sigma w). The tail is odd, so the
valuation is EXACTLY a_0, not merely at least a_0. Consequently Phi_b lands
in X_b and T_b Phi_b=Phi_b sigma. Conversely the finite affine identity
below, followed by 2-adic limiting, reconstructs every x in X_b. Thus
Phi_b is a bijection W -> X_b, with inverse the valuation itinerary.

For odd b,c>=3 the map

    H_{b->c}=Phi_c o Phi_b^{-1}

is a bijection X_b -> X_c, inverse H_{c->b}, and conjugates T_b to T_c.
It preserves every exact valuation word, total division clock, and source
cylinder mass. In fact it is an isometry: if two words first differ after
prefix division clock S, with next digits a!=a', their decoded sources have
v_2 of their difference equal to S+min(a,a'), independently of b.
It does NOT preserve ordinary integrality, positivity, real order, or size.
These are the missing predicates; Haar/slot/pole identities alone cannot
restore them.

## 2. Finite residue lift and a noncircular arithmetic predicate

For a finite word u=(a_0,...,a_{J-1}) put D=S_J and

    C_b(u)=sum_{i=0}^{J-1} b^{J-1-i} 2^{S_i}.
    T_b^J(x)=(b^J x+C_b(u))/2^D.                            (2)

The EXACT source cylinder of u is the single odd class

    r_b(u) = (2^D-C_b(u)) b^{-J} mod 2^{D+1},                (3)

taking its representative in {1,3,...,2^{D+1}-1}. Proof: work backwards
from an arbitrary odd terminal 2-adic y by x=(2^D y-C_b)/b^J. Each inverse
step (2^a y-1)/b is odd and gives exactly valuation a. As y ranges over odd
2-adics its source ranges over exactly the class (3).

For prefixes u_J of an infinite word w, r_J=r_b(u_J) is nondecreasing as an
ORDINARY positive integer, since compatible least residues satisfy
r_{J+1}=r_J+k 2^{S_J+1}, k>=0. This gives the equivalences

    Phi_b(w) in positive odd integers
       iff sup_J r_b(u_J)<infinity
       iff r_b(u_J) is eventually constant.                 (4)

If bounded, this integer sequence stabilizes to n and nested residues imply
Phi_b(w)=n. If Phi_b(w)=n>0, every modulus exceeding n has representative n.
This test is stated entirely in finite carries and modular inverses; it does
not assume convergence or remove words on the basis of their orbit outcome.
It is nevertheless an infinite boundedness condition, not an algorithm
deciding all infinite inputs in finite time.

Define W_b^+ by the boundedness condition (4). Then Phi_b restricts to an
explicit bijection W_b^+ <-> positive odd integers. W_b^+ is shift-invariant
because T_b preserves positive odd integers. A parameter-selective target
for ordinary Collatz is

    Y_{3->5}=Phi_5(W_3^+) subset X_5.                         (5)

H_{3->5} gives a genuine conjugacy from ALL positive odd 3x+1 sources onto
this target: no exceptional integer sources have been discarded. Inversion
uses (1) at b=3, or the eventual value of (3). The target can contain negative
rational sources in the 2-adic 5x+1 system: H_{3->5}(1)=-1 because the word is 2,2,2,... .
Hence an intended convergence theorem on Y must target -1, NOT positive 1.
Collatz is exactly the still-OPEN statement that every w in W_3^+ is
eventually constant 2, equivalently every y in Y_{3->5} eventually reaches -1.
The reformulation is exact and parameter-selective, not itself a proof.

## 3. A finite certificate excludes the full 5-cycle itinerary basin

For a nonempty word u of length J and total clock D the unique periodic
source is

    x_b(u^infinity)=C_b(u)/(2^D-b^J).                        (6)

The denominator is nonzero and odd. The inverse-composition identity gives
an odd rational fixed point realizing exactly u; uniqueness also follows
from (1). For u=(1,1,5):

    b=5: C=39, denominator=3, cycle 13 ->33 ->83 ->13;
    b=3: C=19, denominator=101, cycle 19/101 ->79/101
                                            ->169/101 ->19/101.

This transported orbit is POSITIVE REAL and an odd 2-adic unit. A negative
sign certificate is impossible here: positivity of a b=5 periodic source
forces 2^D>5^J>3^J, hence positivity at b=3 as well. The correct separator is
integrality: gcd(19,101)=1 and 101>1.

**All-preperiod denominator persistence.** If y=p/q is reduced, q odd>1,
then an inverse branch is (2^a p-q)/(bq). Its numerator is coprime to q,
because 2 is a unit modulo q. Therefore q divides the reduced denominator
of every inverse branch. By induction, ANY finite prehistory retains q.
Thus no integer source for T_3 can have an itinerary eventually repeating
(1,1,5), or any rotation. This eliminates the entire hostile itinerary
basin, not just three literal source values, with a finite gcd certificate.
It is an arithmetic exclusion proved from the source lift, not a definition
that excludes known counterexamples. The same lemma handles any specified
periodic word whose b-decoding has reduced denominator >1.

**Quantitative residue witness.** For nonintegral p/q, q>1 odd, its least
residue r modulo M=2^L satisfies q r=p+kM. If M>|p|, r>=1 implies k>=1:
k<=-1 makes p+kM negative, and k=0 contradicts coprimality/nonintegrality.
Hence r>=(M-|p|)/q. Such residues diverge at least linearly in M. In particular
r_3((1,1,5)^k)>=(2^{7k+1}-19)/101, and the exact value is computed by (3).
For any prescribed source bound N, M>101N+19 certifies no positive source
<=N has that repeated prefix. There is no single finite prefix excluding
all positive integer sources.

## 4. Why a finite forbidden-pattern repair cannot do this job

Every finite valuation word u has infinitely many positive integer sources:
its cylinder is r_b(u)+2^{D+1}Z_2, containing all r_b(u)+k2^{D+1}, k>=0.
Thus W_b^+ is dense in W for EVERY odd b>=3. Consequently no rule that
rejects an infinite word upon encountering some finite forbidden prefix
can simultaneously accept ALL positive 3x+1 sources and reject even one
other infinite word. Such a rule would reject the positive representatives
of that prefix. This covers finite-state safety filters, not every possible
finite-state acceptance condition on infinite words.

This explains the required design: the guard must retain unbounded source
arithmetic (such as (4)), use a global acceptance condition, or prove an
additional special theorem. A pattern filter cannot silently replace it.
The 5-cycle passes every finite admissibility test but fails the bounded
residue-source sidecar. Precisely this arithmetic distinction is invisible
to the identical Haar valuation laws and geometric-slot/pole statistics.

## 5. Provenance, controls, and next boundary

Classical conjugacy background: Bernstein--Lagarias, *The 3x+1 Conjugacy Map*,
https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/3x-1-conjugacy-map/6975BB4A8C46CF6842217043AAF9EC13 ;
Edgington, *The autoconjugacy of a generalized Collatz map*,
https://arxiv.org/abs/1206.0553 . These are context, not originality claims;
all accelerated coding and arithmetic statements used here are proved above.

`python3 check_valuation_transfer.py` checks independent inverse-source and
forward-orbit paths, parameter pairs, carry composition, periodic rotations,
arbitrary bounded prehistory controls, monotonicity, and negative/rational
hostiles. Finite tests support implementation; proofs supply all-word scope.

Next genuine obligation: use a special b=3 property to prove that every
bounded-residue ray is eventually 2. The named 5-cycle is rigorously absent,
but excluding arbitrary aperiodic bounded-residue rays remains as difficult
as the corresponding Collatz question; no Haar-zero or finite-prefix claim
is allowed to discharge it.
