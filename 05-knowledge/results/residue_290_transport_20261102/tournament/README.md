# Scalar A000568 coefficients and the recovered black-core route

**PROVED R1–R2, independently audited PASS; R3 and the black-core results
below are INHERITED PROVED with their stated scopes. Fresh checks are
FINITE-EXACT.** [Independent review](../audit_tournament/REPORT.md). This reconstructs mathematical work from an
interrupted October26 session whose uncommitted temporary files are no longer
available. Its old gate counts and audit claims are not reused. No literature
priority or compact coefficient algorithm is claimed.

## 1. Inheritance, correction search and concept board

Research snapshot `4c6f0ca135`; exact full revision and source hashes are in
`SOURCE_PINS.json`. Read routes: Research Protocol; Active Guardrails' tournament
and scalar-loss entries; Core Papers' tournament source map. The mistakes
search found the March odd-cycle-filter repair and the historical false
self-converse=A000568 shift. Our carrier is ordinary isomorphism types, with
converse types kept separate, and only odd permutation-cycle types contribute.

Anchor: a recurrence evaluated from earlier A000568 values. Niche: positive
unit-tail response rank. Wildcard: exact terminal routes after complements.
The closest proved mechanism is
[THM-514, unit-tail core](../../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md)
and [THM-516, gcd profile](../../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md).
The [Hopf/Lie packet §5](../../primitive_hopf_lie_bijection_20261020/PRIMITIVE_HOPF_LIE_BIJECTION.md)
already provides the Euler/ghost transform. The hostile is a ghost coefficient
computed from the very a_n it is supposed to produce. The least-used sidecar
is connected interaction between permutation cycles.

The live board is **cycle interaction / scalar history / positive response
rank / finite target degree / named path starts / terminal core**. Incoming
[Session98 projective response](../../session98_tournament_growth_20261029/recurrence/PROJECTIVE_RESPONSE.md)
already gives a genuinely coarser fixed-target state quotient. R2 concerns
all unit-tail observations of one fixed core mass, not that different
fixed-degree observer. These results are compatible.

## 2. R1: a noncircular scalar recurrence

Let a_n count ordinary unlabelled tournaments, a_0=1. For an odd positive
length tuple l=(l_1,...,l_r), define

    w_ij = 2^gcd(l_i,l_j)-1,
    C(l) = sum_(connected simple graphs G on [r]) product_(ij in E(G)) w_ij.

The singleton C is1. Equivalently, using partitions of the indexed cycle set,

    C(l) = sum_(set partitions pi of [r])
           (-1)^(|pi|-1)(|pi|-1)!
           2^(sum_(B in pi) sum_(i<j in B) gcd(l_i,l_j)).        (1)

For an odd integer partition lambda of k with r parts, write
z_lambda=product_j j^(m_j) m_j!, and set

    g_k = k sum_(lambda odd partition of k)
                2^((k-r)/2) C(lambda)/z_lambda.                (2)

Then for every n>=1,

    a_n = (1/n) sum_(k=1..n) g_k a_(n-k).                     (3)

The symbol g is used to avoid collision with Session98's normalized beta_n.
Formula(2) uses finite arithmetic only: gcd, powers, factorials and partitions.
It uses neither a_k nor unknown connected tournament or graph-type counts.
In particular g_n a_0 is evaluated without first knowing a_n.

**Proof.** In the inherited odd-cycle Burnside formula, assign a length-l
cycle the atom weight 2^((l-1)/2)x^l/l. An unordered collection of r indexed
atoms carries1/r!. Each cycle pair contributes2^gcd=1+w_ij. Expanding these
factors enumerates all simple graphs on the atoms, with edge weights w.
Unique decomposition into connected components gives the formal identity

    A(x)=sum_n a_n x^n=exp(B(x)),
    B(x)=sum_(r>=1) 1/r! sum_(l_1,...,l_r odd>=1)
          C(l) product_i [2^((l_i-1)/2)/l_i] x^(sum_i l_i).

Each coefficient involves finitely many terms. Grouping equal lengths gives
g_k=k[x^k]B, and xA'=xB'A proves(3). Partition-lattice inversion of connected
components gives(1): the coefficient of a graph with h connected components
is1 for h=1 and0 otherwise. No convergence assertion is involved. QED.

The g_k are positive integers. Positivity follows from the graph sum, while
integrality follows from the integer coefficients of xA'/A. The inherited
Euler product identifies g_k=sum_(d|k) d c_d, with c_d the connected sign-even
graph-type count. This explains the coefficients but is not their evaluation
algorithm. Initial g_1,...,g_10 are

    1,1,4,9,41,262,2787,51313,1661368,95410021.

One implementation fixes cycle1 and decomposes a graph by its component:

    C(l)=2^(sum_(i<j)gcd(l_i,l_j))
          -sum_(1 in S proper subset [r])
             C(l_S) 2^(sum_(i<j in S^c)gcd(l_i,l_j)).          (4)

Repeated lengths require exact subset multiplicities. The checker groups
those multiplicities; literal connected graphs and(1) provide independent
small controls. Dropping gcd interactions first fails at total mass6, from
the pair(3,3): the resulting a_6 is164/3 instead of56.

**Cost boundary.** Equation(3) meets the literal earlier-sequence-values
criterion. Its coefficients still grow through odd-cycle profiles and a
connected interaction calculation. The all-unit profile1^k already contains
the number of connected labelled graphs on k vertices. This is explicit
arithmetic, not an unknown tournament input, but no complexity improvement
on Burnside is proved. Session98's projective response recurrence remains
the stronger existing algorithmic quotient; its changing auxiliary states
are different from this scalar interface. A cheaper g_k kernel remains OPEN.

## 3. R2: exact unit-tail memory, with a first hostile at mass9

For fixed nonunit core mass m, THM-514 gives positive rational B[m,t] for

    T_m={t>=1:3t<=m, t=m mod2}.

Use T_0={0}; T_1,T_2 are empty. After removing the universal r-unit tail
factor2^binom(r,2)/r!, its response is

    Y_m(r)=sum_(t in T_m) B[m,t] (2^t)^r.

Its minimal constant-coefficient rational linear recurrence has order
|T_m| and characteristic polynomial product_t(X-2^t). The first |T_m|
responses determine all weights; the square Hankel matrix of that size is
positive definite. For m>=3 the order is floor(m/6) when m is even and
floor((m+3)/6) when m is odd, so these orders are unbounded.

**Proof.** Write the Hankel matrix as V diag(B[m,t]) V^T, with distinct-node
Vandermonde V_it=(2^t)^i. All weights are positive. The displayed polynomial
annihilates the response. Any other annihilator p must have
B[m,t]p(2^t)=0 for every t, by linear independence of the exponentials, so it
needs every displayed root. The Vandermonde inverse recovers the weights.
Counting the admissible t proves the order formulas. QED.

At m9, the profiles(9) and(3,3,3) have counts1 and3. Two *altered weighted
models* concentrating unit weight on either profile have equal initial
response1 but next responses2 and8. If lower masses agree, their complete
weighted series agree through degree9 and differ by6 at degree10. This is a
universal weighted-model obstruction, not two possible actual a_10 values.
The canonical sequence has fully specified arithmetic weights.

These are transfer eigenvalues2^t, not adjacency eigenvalues of a tournament.
R2 does not exclude nonlinear recurrences or the finite-horizon/projective
quotients in Session98. It tells exactly what is lost by a single response
when all future unit-tail lengths remain observable.

## 4. R3 and the black core are already published

**INHERITED PROVED R3.** A000568 is not P-recursive. The elementary growth
proof is in
[Session98, limits](../../session98_tournament_growth_20261029/limits/LIMITS_AND_NORMALIZATION.md):
a_n>=2^binom(n,2)/n!, whereas any fixed-order polynomial-coefficient linear
recurrence bounds solutions eventually by C^n(n!)^D. Thus logarithmic growth
would be both quadratic and O(n log n), impossible. This does not rule out
coefficients involving2^n, increasing order, nonlinear recurrences or(3).

The originally pending regular black-distance theorem was independently
reconstructed and audited in
[Session75 black-core geodesics](../../session75_overnight_black_root_spectral_20261027/session75-black-core/BLACK_CORE_GEODESICS.md),
with [independent audit](../../session75_overnight_black_root_spectral_20261027/session75-black-audit/ADVERSARIAL_REPORT.md).
For regular order2m+1>=7, exact black distance to the nonstrong locus and to
{1→C,C→1:C strong on2m vertices} is m. A literal endpoint-current construction
uses two different named path starts; this simultaneously excludes later blue
moves and a second singleton sink. C3 and R5 are pure-blue boundaries.

The current stronger
[THM-4819(D), prescribed black core](../../../../01-canon/theorems/THM-4819-regular-tournaments-under-native-moves-the-transitive-slab-exact-distance-and-the-prescribed-black-core.md)
allows any prescribed vertex and either native role. Its endpoint lemma replaces
the pure-blue classification for the first black word. The full all-colour
non-strong distance is arc-connectivity, from
[THM-4811](../../../../01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md).
These targets must not be replaced by the transitive or forest-frame targets.
Global black-core connectivity remains outside the present work.

This packet's fresh controls enumerate every Hamilton path of three pinned
regular order7 representatives and extend every black first word by a second
word with a different actual start. Every move is replayed by literal nonpath
complement and an alternating converse gauge. Original vertex identities and
word order remain intact. A same-start route gives the hostile
1→R5→1 and contains blue moves. Source measures would push by summing original
occurrence weights; no fresh uniform choice or clone merger is introduced.

## 5. Reproduction and contribution boundary

Run `python3 -I -B check.py` and `python3 -I -B -O check.py` here. The outputs
must match. The checker imports no repository implementation. It compares
connected graph enumeration, set-partition inversion, grouped root-component
recursion, odd-cycle Burnside counts, and literal permutation actions. Its
finite universes, receipts and source pins are explicit. These checks do not
prove all-order results without the written arguments.

The contribution is the reconstructed explicit coefficient formula and its
response-rank explanation, plus a current recovery map for the old black
result. It is not rediscovery of the now stronger published theorem. The
four-marker/support-image guards remain relevant for individual tournament
attachments: scalar g_k counts cannot reconstruct their legal ports.
