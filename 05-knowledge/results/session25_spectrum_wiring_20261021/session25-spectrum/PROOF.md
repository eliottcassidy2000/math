# Strong actual-H spectra have every parity-admissible finite residue

Session 25, 2026-10-21. Status: PROVED by the elementary arguments below, with independent exact controls; separate independent reviewer audit passed the mathematical derivation and independently implemented controls. No literature-priority claim. This packet is written outside the repository; it performs no repository or Git mutation.

## 0. Scope, inheritance, and concept board

H(T) is the actual number of directed Hamiltonian paths, not H(T)/|Aut(T)|. A strong tournament in this packet is strongly connected and has at least three vertices. Its feedback vertex number is the minimum number of vertices whose deletion leaves an acyclic tournament.

Closest proved mechanisms: THM-012b (`01-canon/theorems/THM-012b-insertion-decomposition.md`) and the corrected THM-1900 (`01-canon/theorems/THM-1900-the-insertion-response-calculus.md`) retain orphan paths in vertex insertion; THM-4097 (`01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md`) retains the full oriented insertion response. THM-4094 (`01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`) already reduces the global spectrum conjecture to two strong-prime lanes. Later ear certificates already far exceed small-count coverage. These are inherited, not new results.

Canonical hostile: C3 has H=3, whereas inserting its third vertex into the unique path of its transitive two-vertex deletion supplies only two paths. The third is an orphan. Corrected near miss: scalar H and a cut-size statistic do not determine insertion or continuation. Least-used sidecar here: the binary integer carried by the entire insertion word, which becomes a finite exact affine state.

Concept board: (1) actual versus orbit-normalized H; (2) prime-lane attainment versus residue attainment; (3) orphan-preserving endpoint counting; (4) strong one-feedback-vertex families; (5) finite-state block acceleration. The result connects (3)--(5), places an obstruction on congruence approaches to (2), and never crosses the normalization firewall (1).

## 1. Main theorem

For every positive integer M and residue r modulo M, there exists a strong tournament T with feedback vertex number exactly one and

    H(T) = r (mod M)

if and only if either M is odd, or r is odd. Moreover, infinitely many distinct actual path counts occur in each such residue class, witnessed by this same explicit family.

The necessary condition is Rédei parity. The construction and proof below are elementary and do not assume the conjecture that 7 and 21 are the only forbidden odd counts.

Consequences:

- No finite congruence obstruction beyond parity exists for actual H, even after restricting to strong tournaments with feedback vertex number one.
- In particular, 7 and 21 cannot be excluded by any necessary condition depending only on H modulo one fixed modulus, or finitely many fixed moduli (combine them by their least common multiple).
- For any finite list of odd primes and any prescribed finite nonnegative valuations at those primes, some such strong tournament realizes all those valuations simultaneously. Apply the theorem to a CRT residue that is p^e modulo p^(e+1) at each prime.
- This does NOT realize every permitted integer. A residue witness can be enormous. It neither proves actual spectrum completeness nor produces a prime-valued H for each prime.

## 2. Exact insertion-word counting

For a binary word w=b_1...b_n, define T(w) on vertices x,1,...,n. The old vertices are transitively oriented i->j for i<j, and x->i exactly when b_i=1. Let H(w)=H(T(w)). For the empty word, H=1. For a nonempty word,

    H(w) = b_1 + (1-b_n) + #{i : b_i=0, b_(i+1)=1}
           + sum_(j<i) b_j(1-b_i) 2^(i-j-1).                 (1)

Proof. Paths starting or ending with x contribute the first two terms. If x is internal, the old vertices before and after x must each be increasing. Let i be the last before x and j the first after x. The port conditions are b_i=0,b_j=1. When i<j, all vertices can be covered only when j=i+1, giving the adjacent 01 term. When j<i, vertices below j are forced before x, vertices above i are forced after x, and each of the i-j-1 intervening vertices can independently be assigned to either increasing segment. This gives exactly 2^(i-j-1) paths. These cases are disjoint and exhaustive. In particular the last sum retains the orphan paths rather than deleting them.

For a prefix let

    A = sum_(j=1)^n b_j 2^(n-j),     p=b_n,

with A=p=0 for the empty prefix. Then (1) gives the exact state update

    append 1:  (H,A,p) -> (H,       2A+1,1),
    append 0:  (H,A,p) -> (H+A+p,  2A,  0).                 (2)

Appending 1 creates a new 01 contribution precisely when it removes the previous end-x contribution, so H is unchanged. Appending 0 creates new orphan weight A and increases the end-x contribution by p. Thus (2) is proved independently of fitted data. The bit p is redundant if the exact A is retained, since p=A mod2.

Scalar H alone is insufficient: words 00 and 01 both have H=1, but appending 0 gives H=1 and H=3 respectively. The continuation state is not determined even by the old tournament's isomorphism class: both old tournaments in that hostile are transitive.

## 3. Constant-run block acceleration

Fix L>=1, Q=2^L, and let

    w_(L,t) = (1^L 0^L)^t,      t>=1.

At the end of each block p=0. Starting a block with state (H,A,0), the run 1^L changes A to QA+Q-1 and leaves H unchanged. The following run 0^L therefore gives

    A' = Q^2 A + Q(Q-1),
    H' = H + Q(Q-1)A + D,       D=(Q-1)^2+1.                (3)

The additional +1 in D comes from the first zero's previous-bit term; omitting it gives the wrong count already for C3. In homogeneous coordinates (H,A,1), (3) is an exact affine integer matrix, with rows

    (1, Q(Q-1), D),  (0,Q^2,Q(Q-1)),  (0,0,1).

It is a proved constant-run acceleration of (2), not a statistical recurrence. Iterating from H_0=1,A_0=0 gives

    A_t = Q(Q^(2t)-1)/(Q+1),
    H_(L,t) = 1 + Q^2(Q^(2t)-1)/(Q+1)^2 + 2t/(Q+1).        (4)

The expression is an integer because it solves the integral recurrence (3). It is strictly increasing in t, since D>0 and A>=0.

For L=1 the first eight H values are

    3,9,31,117,459,1825,7287,29133.

Every T(w_(L,t)) is strong: its first bit is1 and its last bit is0, so x reaches every old vertex through vertex1, and every old vertex reaches x through vertex n. These paths also provide reachability between any two old vertices. Deleting x leaves a transitive tournament; strongness rules out feedback vertex number0. Hence the number is exactly1.

## 4. Modular tuning and proof of the main theorem

Write M=2^k m with m odd. Choose an integer L>=max(1,k) such that

    2^L = 1 (mod m).                                       (5)

Such L exists because 2 is a unit modulo m; when m=1 this condition is vacuous. For example, take a sufficiently large positive multiple of the multiplicative order of2 modulo m, interpreting that order as1 for m=1.

Then Q=2^L is0 modulo2^k and1 modulo m. Therefore Q(Q-1)=0 modulo M. Equation (3) shows, inductively,

    A_t = 0 (mod M),       H_(L,t) = 1+tD (mod M).          (6)

Also D=1 modulo m. If k>=1, D=2 modulo2^k; for k=1 this simply says D is even. Consequently

    gcd(D,M)=gcd(2,M)=:g.                                  (7)

The progression 1+tD for 1<=t<=M/g visits exactly all residues congruent to1 modulo g. These are all residues for odd M and precisely the odd residues for even M. This proves sufficiency and gives a constructive algorithm: solve tD=r-1 modulo M, choosing a positive t<=M/g, then output the insertion word above. The witness has at most

    2 L M/g + 1 vertices.                                 (8)

All t in the same solution class modulo M/g give the same target residue. Equation (4) is strictly increasing, so these supply infinitely many distinct H values, not merely infinitely many source/sink paddings of one count. All witnesses remain strong.

### Congruence-complete does not mean numerically dense

The entire two-parameter set B={H_(L,t):L,t>=1} already has all the modular coverage just proved, yet it is extremely sparse. From (4), Q>=2 and Q^(2t)>=4 imply

    H_(L,t) >= [Q^2/(Q+1)^2] [1-Q^(-2t)] Q^(2t)
              >= (4/9)(3/4) Q^(2t) = 4^(Lt)/3.

For X>=3, if H_(L,t)<=X, then Lt<=B_X=floor(log_4(3X)). Hence the number of distinct family values at mostX is bounded by

    sum_(L=1)^B_X floor(B_X/L) <= B_X(1+log B_X)
                               = O(log X log log X).

Duplicates only decrease this count. Thus a zero-density subset of the actual strong spectrum is already congruence-complete at every finite modulus. This supplies a constructive sharp warning against promoting local congruence evidence to global integer coverage. The bound is on this family, NOT on the full actual-H spectrum.

The tuning condition matters. For L=1, modulo3 equation (4) gives H_(1,t)=(4^(t+1)+6t+5)/9; the exact block recurrence gives H_(1,t)=1+2t^2 modulo3, so only residues0 and1 occur. The claim is NOT that every fixed block family covers every modulus. For M=3 choose L=2 instead, as (5) requires.

## 5. A precise connection and its losses

Source: binary transitive-ear words with a marked inserted vertex. Target: affine state (H,A,p), then the block state (H,A,1). Map: exact endpoint/orphan counting (1), followed by the proved run composition (3). Preserved target: actual directed Hamiltonian-path count and its modular residue. Destroyed information: generally the full tournament and its isomorphism class, automorphism group, individual Hamilton paths, and all metagraph adjacency. Required sidecar: A, plus p before complete blocks. Cheapest decisive test: 00 versus01, both H=1 but distinct appending-zero responses.

The resulting stateful acceleration is a legitimate finite dynamical model. A paper using arithmetic state refinements or accelerated runs can be compared at the level of proof architecture, but no transfer of a halting, Collatz, metagraph-wiring, or prime-realization theorem follows. In particular the displayed affine matrix is this tournament model's matrix, not an asserted identification with a paper's S2 or S3 maps.

## 6. Reproducible exact controls and negative scope

Run from this packet's directory:

    python3 verify.py
    python3 -O verify.py

The standalone verifier has no third-party dependencies, no repository writes, and no assertion statements that disappear under -O. It checks:

- all1023 binary words of lengths0 through9 against independently implemented vertex-level subset DP;
- all127 words of lengths0 through6 against literal permutation enumeration;
- all1023 words against the direct slot formula;
- transitive-closure strongness for all255 words of lengths2 through9 with firstbit1,lastbit0;
- the closed form against the bit-by-bit recurrence for200 (L,t) pairs;
- all parity-admissible residues for every modulus1 through512, totaling98432 block-residue instances, with direct exact closed-form checks for moduli through32;
- scalar-state, untuned-block, and orphan-retention hostiles.

The full receipt is `receipt.json`. These finite controls audit the construction; the proof, not the finite tests, supplies every modulus and arbitrary t. Incoming f002fc4/fad6a405 arithmetic/cap225 claims are not dependencies. No global integer-attainment result, density-one result, stronger universal obstruction, or arbitrary strong-atom census is established.
