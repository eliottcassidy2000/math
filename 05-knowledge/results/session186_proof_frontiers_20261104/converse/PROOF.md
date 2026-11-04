# Homogeneous pairs cannot hide a module converse

Status: **PROVED, inherited all-order theorem and explicit gate integration.** No exact full-response collision is claimed. This packet reconciles the older open homogeneous-pair branch with the incoming exact sign theorem, rather than presenting that theorem as new.

## 1. The current source resolves the stale branch

Sources were read from main at `7af964816f3ae8f41137d2a97c02a87cca0bfcc8`:

- `05-knowledge/results/two_port_eight_core_20261103/algebra/PROOF.md`, C2: exact opposite-twin sign separation, even with different ambient quotients and positive projective scaling.
- The same source, C4: any twelve-vertex hidden module-converse core must be prime, using the separately declared finite residue dependencies.
- `05-knowledge/results/two_port_twelve_principles_20261103/PROOF.md`, T4: the **global converse** response difference is a named cubic.
- `05-knowledge/results/session185_proof_frontiers_20261104/alternative/PROOF.md` is the incoming three-state gate packet.

The older twelve-principles structure document leaves the pair branch open. C2 already closes it. This is an inheritance reconciliation, not a new unbounded reconstruction theorem.

## 2. Exact named observer and proof

Let `F_T(a,c)` count directed Hamilton paths after adding physical names x,y with x→y, and arcs u→x iff a_u=1, u→y iff c_u=1. Keep every named Boolean coefficient. The defining path expansion includes the full three-path contribution, without asserting that its individual tensor entries are separately observable. The inherited linear tomography recovers `C_uv`, the number of ordered covers by two nonempty directed paths, the first ending at u and the second starting at v, from this full response. Thus full equality implies equality of each named C entry; the converse implication is **not** assumed.

Suppose the physical pair {a,b} is homogeneous with a→b. Contract it to z, obtaining Q. For each exterior name w,

    C_bw(T) − C_aw(T) = C_zw(Q),
    C_wa(T) − C_wb(T) = C_wz(Q).

For the first identity, swap the twin names in every cover counted by C_bw unless its first path ends in the consecutive block a,b. All other covers biject with C_aw: a cover ending its first path at a cannot contain consecutive a,b. Covers with the twins in different paths participate in this bijection. The unmatched block contracts to precisely a cover counted by C_zw(Q). The second identity is the analogous initial-block argument in the second path.

Consequently the exact contrast

    L_(a,b;w)(T) = C_bw − C_aw + C_wa − C_wb
                 = C_zw(Q) + C_wz(Q) > 0.

Strict positivity is the inherited positive cut-capacity theorem; its sum is twice a positive integer. If the same named pair is homogeneous but has b→a in another tournament, its L is strictly negative. This does not require the two quotients to agree. Thus the complete named responses cannot agree or be positively proportional.

Now let M be any module in any ambient tournament R[M], and suppose M contains a homogeneous pair. Module transitivity makes that pair homogeneous in the ambient tournament. Replacing M by its converse reverses the pair while retaining its homogeneity. The sign certificate therefore excludes equality of the actual ambient responses. At ambient order two, the literal query a=0,c={a} gives response1 for the arc a→b and response3 for b→a (here a in c={a} denotes the physical vertex name). Thus no exterior w or endpoint tomography is needed in that special case. This covers every nontrivial homogeneous-pair case, without an ambient size bound, strongness assumption, or exact-response descent assertion.

**Conclusion.** A hidden module-converse core cannot contain a homogeneous pair. The pair branch is closed; prime even nonselfconverse cores are not closed by this argument.

## 3. Integration of the three-state gate

Use the incoming physical seed rows

    Q0=(16,5,33,39,46,3),
    Q1=(16,5,49,39,42,3),
    Q2=(18,4,49,39,42,3).

Replace name5 by the increasing transitive module X_r on names5,...,r+4. The incoming occupied-gap proof gives, for every r≥1,

    H(Qj(r)) = 4·3^r + 9·2^r − 3 + 2j,   j=0,1,2.

These are strong tournaments, with successive physical single-edge reversals {2,4} and {0,1}. Their 3^r and 2^r channels are both genuine. This construction is positive count arithmetic, not a collision in the full named response.

For every r≥2 the adjacent names5,6 are a homogeneous pair with5→6. Therefore reversing X_r, or reversing any module containing that pair, is detected by the exact sign. Contracting that pair leaves Qj(r−1), so for every exterior name w the certificate is explicitly

    L_(5,6;w)(Qj(r)) = C_zw(Qj(r−1)) + C_wz(Qj(r−1)) > 0.

The reversed-module tournament has the opposite sign. The hidden-converse obstruction and the unit count gate coexist: neither invalidates the other. In particular an isomorphism obtained by reversing the physical clone order is not equality of named responses.

For an arbitrary substituted module A, the incoming identity is

    H(Qj[A])=(27+2j)H(A)+21q2(A)+4q3(A).

It still gives the step2H(A), but the pair obstruction applies whenever A contains twins; it makes no assertion about arbitrary prime A. The three-state carrier alone still does not supply a reusable ternary digit operation.

## 4. Preserve the cubic without overextending its scope

For a **global** converse comparison, T4 proves

    Δ_T(a,c)=F_T(a,c)−F_(T^op)(a,c)

has total degree≤3, bidegree≤(2,2), and no pure-port quadratic terms. Its complete named coefficients are recoverable from1+2n+n³ prescribed queries.

An ambient module reversal R[M] versus R[M^op] is generally not a global converse. T4 does not by itself justify discarding degree-four terms there. This packet interpolates the complete bidegree-(2,2) difference for its ambient examples, and separately retains the complete cubic for each global converse control. No sorting of coefficients, endpoint-marginal quotient, or renaming of vertices is used.

`verify.py` supplies an independent standard-library endpoint subset DP and an ordered-cover convolution. For all three r=2 gates it records actual adjacency, every C entry, all five exact contrasts, the entire named ambient difference polynomial, and the entire named global-converse cubic. Each interpolation uses every assignment of weight≤2 in each port. Noninterpolation assignments are then checked literally. Coefficients not listed are explicitly zero within the complete bidegree-(2,2) basis. The cubic-degree and pure-quadratic gates are tested on the global converse only.

All finite checks are controls. The all-order result is the combinatorial argument above, credited to incoming C2. Universal two-port reconstruction and the remaining prime-core converse problem remain open.

## 5. Explicit ambient quartic hostile

For r=2 the three gate counts are69,71,73. The exact positive L vectors over outside names0,...,4 are respectively

    (10,18,18,24,8), (18,20,18,24,10), (22,16,14,28,14).

Reversing the internal clone arc negates each vector. The ambient module-reversal differences have degree exactly4, with respectively396,380,388 nonzero named Boolean coefficients. The global converse differences have degree exactly3, with respectively342,350,346 nonzero coefficients. Therefore the distinction in Section4 is necessary on these very gates, not merely a missing argument. The elementary query a=0,c={5} already gives ambient differences−56,−60,−64. All names are physical; a=0 is the actual universal-source placement for x.

A particularly small named quartic coefficient is

    [a_0 a_1 c_0 c_5](F_Qj(2) − F_(clone-reversed Qj(2)))
      = −2, −4, −4  for j=0,1,2.

These are masks (3,33), evaluated by all sixteen terms of the exact Boolean finite difference. The same physical old name0 appears once in each port's attachment variables; this is allowed and is not identifying the two variables.
