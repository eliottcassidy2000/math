# Authenticated port currents and the complete dyadic lifting tower

23 October 2026. **PROVED:** the all-size graph/lattice and port-certification
statements below. **FINITE-EXACT:** independent complete tournament controls
through order six, all valid disjoint two-edge switches of the order-five
paired-black graph, and the full order-eight electrical calculation recorded
separately. **NOT CLAIMED:** all-order connectivity, a fixed-size incidence
signature, or an electrical characterization of actual tournament ports.

## 0. Inheritance and the exact increment

The graph is the converse-merged fixed-Hamiltonian-path **complement-line**
multigraph. It is not an arc-flip graph. Session36's
`DYADIC-ELECTRICAL-METAGRAPH.md` already proves:

- every black line has a distinct reflected twin, including loops;
- writing G = Blue + 2 PairedBlack gives A = A_L + 2 A_B;
- critical 2-rank = c_L - c_G + bicycle_dimension(L);
- the first half-energy form q_1 on blue-harmonic potentials has radical
  dimension equal to the number of Smith factors divisible by four;
- black-component contraction, two-color tree coefficients, and their
  conductance-limit interpretation.

The incoming `metagraph_bucket_obstruction_20261023/dyadic_bridge.py` already
exhibits one order-six counterfeit detected by q_1. The incoming port-rigidity
and two-height results already give all-size exact acceptance when **complete
authenticated path-labelled data** is supplied.

The increment here is (i) an exact criterion for which degree-preserving
black switches q_1 detects, (ii) a complete intrinsic dyadic lifting tower,
(iii) a sharp counterfeit invisible even to the entire tower, (iv) its
full-network order-eight evaluation, and (v) the precise connection between
authenticated target buckets, currents and a nonnegative mismatch energy.
The dyadic mechanism is standard integral linear algebra, not a claim of
priority for a new general Smith theorem.

All component counts include isolated vertices. Electrical Laplacians omit
loops, while graph inventories retain them and count a loop twice in degree.
No statement below assumes that the actual metagraph is connected at all
orders. Choose one root in each full connected component and delete its
row and column when a nonsingular matrix is needed.

## 1. What the half-form sees about a black switch

Let W = ker_F2(A_L), and let V be W modulo potentials constant on each
full component of G. Define

    q_1(u,v) = u_hat^T A v_hat / 2 modulo 2.

Consider the degree-preserving switch in the paired-black multigraph

    remove ab and cd; add ac and bd,

where a,b,c,d are distinct and the removed edges exist. Each displayed
edge represents two literal black lines. Put

    x = e_a - e_d,       y = e_b - e_c.

The exact integer change is

    Delta A_B = x y^T + y x^T,
    Delta A   = 2(x y^T + y x^T).                       (1)

This follows by inspecting the four off-diagonal entries; all diagonal
changes vanish. The blue graph and hence W are unchanged. On W the
half-form changes by

    Delta q_1(u,v)
      = (u_a-u_d)(v_b-v_c) + (u_b-u_c)(v_a-v_d) mod 2.  (2)

If the switch is performed within a full component and preserves that
component partition, (2) descends directly to the same V. In full generality
it is always a valid comparison on the common space W, without assuming
matching component quotients. This qualification avoids silently identifying
quotients when a switch disconnects a graph.

**Switch detection theorem.** Write alpha(u)=u_a-u_d and beta(u)=u_b-u_c
on W (or the common V when applicable). Then Delta q_1 is nonzero if and
only if alpha and beta are linearly independent over F2. When nonzero it
has rank two.

Proof: if one functional vanishes or they are equal, the two terms cancel.
Otherwise (alpha,beta):W -> F2^2 is onto. Pulling back the nonsingular form
with matrix [[0,1],[1,0]] gives rank two. This is both directions, not merely
a sufficient test.

An equivalent computable criterion uses endpoint vectors modulo im(A_L):
alpha and beta are independent exactly when x and y have independent
classes in F2^vertices / im(A_L), since W^perp = im(A_L).

**Blue-isolated consequence.** If a,b,c,d are four distinct blue-isolated
vertices, all four unit potentials lie in W, so the two functionals are
independent. Every such unit switch is detected. In particular this applies
to four distinct non-self-converse merged nodes at every tournament order.

More generally, if i,j are distinct blue-isolated vertices, then

    q_1(e_i,e_j) = -m_B(i,j) modulo 2.                  (3)

Thus this principal block of the half-form records every paired-black
multiplicity modulo two, equivalently each literal black multiplicity
modulo four. It is not just a degree statistic. Blue-component indicator
potentials similarly recover the paired-black quotient Laplacian modulo
two after contracting the blue components.

The diagonal cannot detect the switch: the right side of (2) vanishes
when u=v. The obstruction is intrinsically bilinear.

## 2. A counterfeit beyond the entire first half-form

An independent census of all fixed-path tournaments through order six
rebuilds the same canonical node labels as the Session36 small controls.
At order five, replace paired-black edges

    1--4 and 3--7 by 1--7 and 3--4.

The latter edge already exists, so its multiplicity increases by one.
Every blue edge, every color degree, every blue/black orbit count, full
connectivity, the full mod-two Laplacian, and the **entire labelled q_1**
are unchanged. But exact spanning-tree counts differ:

    actual: 41096 = 8 * 5137,
    switched: 36680 = 8 * 4585.

The actual q_1 is nondegenerate on its three-dimensional space. Since the
switched q_1 is identical, both complete 2-primary groups are (Z/2)^3:
there are no further nonzero dyadic layers in either graph. Thus even the
whole dyadic tower cannot force the actual wiring. Here the odd tree factor
already distinguishes them. This is a synthetic colored-network control,
not another tournament realization.

The checker exhausts all 42 connected switches arising from disjoint
paired-black edges at order five: 25 are detected by q_1, and 17 are not.
Both direct matrix recomputation and formula (2) agree in every case. The
receipt preserves full original and altered edge lists, not only local IDs.

A still stronger boundary is unavoidable: changing the matching of actual
ports while keeping the same class-to-class edge multiplicities leaves the
entire Laplacian unchanged. No electrical invariant of the merged graph
alone can authenticate which individual path-labelled ports were paired.

## 3. The complete intrinsic dyadic tower

Let A be any nonsingular symmetric integer matrix, in particular a reduced
Laplacian. For k>=1 define

    M_k = {z in Z^r : A z is divisible by 2^k},
    V_k = reduction of M_k modulo 2, as a subspace of F2^r.

For u,v in V_k choose lifts u_hat,v_hat in M_k and set

    q_k(u,v) = u_hat^T A v_hat / 2^k modulo 2.          (4)

**Tower theorem.** This is a well-defined symmetric bilinear form, and

    rad(q_k) = V_(k+1).                               (5)

Consequently

    dim V_k = number of Smith factors divisible by 2^k,
    rank q_k = number of Smith factors with valuation exactly k.  (6)

In particular V_1 is precisely the blue-harmonic space modulo full-component
constants for G=L+2B. The inherited half-form is the first stage, not the
whole obstruction hierarchy.

Proof of well-definedness: changing v_hat to v_hat+2z changes (4) by
2 u_hat^T A z / 2^k, an even integer because A u_hat is divisible by 2^k.
The other variable is identical. Sums of qualifying lifts qualify, giving
bilinearity and symmetry.

Proof of the radical formula: first, a lift y_hat with A y_hat divisible
by 2^(k+1) pairs trivially with every u_hat in M_k, so V_(k+1) is in the
radical. Conversely take y_hat in M_k and put b=A y_hat/2^k. Its binary
class lies in rad(q_k) precisely when

    b^T z = 0 modulo 2 for every z in M_k.              (7)

The pairing on (Z/2^k Z)^r is perfect, and for symmetric A,

    im(A modulo 2^k) = (ker(A modulo 2^k))^perp.

This equality follows directly by Smith normal form, or from containment
and equal finite cardinalities. Condition (7) is therefore equivalent to
solvability of

    A t = -2^(k-1) b modulo 2^k.

Then y_hat+2t has the same binary reduction and A(y_hat+2t)=0 modulo
2^(k+1), proving the reverse inclusion. Formula (6) follows from Smith
normal form: a binary coordinate admits a lift at level k exactly when its
Smith factor is divisible by 2^k. Domain changes are unimodular and preserve
binary reduction dimensions. No diagonalization of a symmetric form at
prime two is assumed.

**Important computational guard.** A radical vector at stage k need not
lift by merely adding 2^k times a binary vector to an arbitrarily chosen
previous lift. Lower bits may need correction. The proof supplies the
correct congruence A t = -2^(k-1)b modulo 2^k. Any deeper implementation must
persist and verify actual qualifying lifts, not repeatedly divide a
binary form or treat an uncorrected representative as a lift.

For a specified integer current j, the exact complementary test is

    A p = j modulo 2^k is solvable iff
    z^T j = 0 modulo 2^k for every z in ker(A modulo 2^k). (8)

All dyadic tests characterize solvability over the 2-adic integers, not
integer solvability: odd-primary obstructions remain. For a full Laplacian,
also require zero total current on each full connected component.

## 4. Authenticated ports turn these into wiring constraints

Let f be any integer-valued invariant of merged tournament nodes. Given an
authenticated source port p at node s(p), decode its actual complement by
the inherited signed-curl or two-height transport, and call its target
value t_f(p). Define the intrinsic current

    J_f(v) = sum over ports p owned by v of [f(v)-t_f(p)].

Every literal incident line contributes once at an endpoint; loop ports
contribute zero. The actual graph necessarily satisfies

    A_G f = J_f,
    A_B f = [J_f - A_L f]/2.                           (9)

The first identity is the inherited Kirchhoff law specialized to genuinely
authenticated target data. The second is its exact twin-reduced version.
It is useful because a candidate with correct blue wiring and degrees can
still fail these prescribed currents, including their higher dyadic tests.
It is not permission to compute J from the candidate wiring itself, which
would make the test tautological.

For the switch of section 1,

    Delta(A_G f) = 2[x(f(b)-f(c))+y(f(a)-f(d))].         (10)

Over the integers, a disjoint switch preserves this one current exactly
when f(b)=f(c) and f(a)=f(d). The halved current modulo two uses equality
modulo two. Several selected incidence observables can therefore certify
individual ambiguous switches, but simultaneous switches may cancel and
must not be excluded merely from the single-switch criterion.

### A necessary-and-sufficient nonnegative energy, with its exact scope

Consider nonnegative integer transportation entries X_(p,q) from source
port buckets p to candidate target buckets q, with fixed row and column
margins. Suppose every source row has one authenticated predicted target
signature z_p in Z^d, and each target bucket has signature z_q. Then

    E(X) = sum_(p,q) X_(p,q) ||z_q-z_p||^2 >= 0,

and

    E(X)=0 iff every occupied pair has z_q=z_p.         (11)

This is simply a sum of nonnegative integer squares. It is the Dirichlet
energy of the bipartite candidate-routing graph with prescribed signature
potentials. Fixed margins determine the sum of the squared norms, so the
only routing-dependent term is -2 sum X_(p,q) z_p dot z_q.

If the signature only records a short incidence profile, zero energy means
exactly compatibility with that profile, not tournament realizability.
If it is the complete normalized signed-root-curl signature of the source's
predicted complement and the target's actual port, the inherited decoder
is injective on fixed-path tournaments. Thus zero energy, margins, and the
required involution matching give the exact all-size port acceptance test.
This statement does not construct the inventory from counts, reduce the
state-space information, or prove that a short profile is always complete.
The positive energy certificate rejects one mismatched authenticated port
regardless of whether its degree/parity/electrical summaries cancel elsewhere.

## 5. Reproduction and current finite claims

Run `python check_electrical.py` and again with `python -O`. It independently
enumerates all fixed-path tournaments at n=3..6, checks Smith forms at those
sizes, checks every connected disjoint black switch at n=5, and verifies
(5) by exhaustive residue-lift enumeration for small symmetric matrices,
including even hyperbolic blocks and mixed valuations. Checks use explicit
exceptions and survive optimized Python.

The `--n8 PATH` option consumes the full path-replayed order-eight edge
ledger, with literal black counts divided by two only after checking their
evenness. It computes the binary blue kernel and its full labelled half-form
with sparse incidence arithmetic, retaining isolated nodes and loops in the
source inventory. No huge determinant is used. The next checker stores
actual mod-four lifts of every first-form radical basis vector and directly
verifies their Laplacian divisibility before computing q_2.

The all-size theorems are independent of the finite n8 outcome. The latter
is evidence about this one actual graph, not an all-size exponent bound.

## 6. Full order eight: the exponent-four pattern breaks decisively

The full order-eight ledger contains 3,528 converse-merged nodes and is
connected by direct traversal. Its blue graph has one 176-vertex component
and 3,352 isolated vertices. Its blue bicycle dimension is 11. Independent
bit-packed projection gives

    dim V_1 = 3363,      rank q_1 = 3358,
    dim V_2 = 5,         rank q_2 = 4,
    dim V_3 = ... = dim V_8 = 1,
    rank q_3 = ... = rank q_7 = 0,
    rank q_8 = 1,        V_9 = 0.

Consequently the complete 2-primary critical group is exactly

    K(G_8)_(2) = (Z/2)^3358 + (Z/4)^4 + Z/256,
    v_2(kappa(G_8)) = 3374.                            (12)

There are 164 additional Smith factors with valuation zero; (12) does not
compute their odd parts. The degree-eight dyadic exponent is an exact
finite result, not an extrapolation. In particular the full-network
exponent-at-most-four pattern at n<=7 does not continue at n=8.

The blue half-form has rank four and is nonalternating, agreeing with the
inherited independent blue-only n8 census. In one explicitly saved basis,

    q_2 = [[0,0,0,0,0],
           [0,0,1,0,1],
           [0,1,0,1,0],
           [0,0,1,1,0],
           [0,1,0,0,0]].

The unique surviving binary mode has 1,734 nonzero node coordinates,
including 88 of the 176 self-converse nodes, in this grounding. These support
counts are basis/grounding-dependent descriptive checks, not intrinsic
claims or a proposed short combinatorial formula for the mode.

Every basis vector of V_2 is stored as an actual integer potential with
A u divisible by four. The one-dimensional continuation stores a qualifying
potential at every modulus 8,16,32,64,128,256, together with its entire
integer vector A u / modulus. Therefore an auditor can reconstruct each
energy from raw edges and check its residue without trusting an elimination
routine. The last energy is odd after division by 256, so there is no
remaining 512-divisible factor. There is no unresolved precision tail.

### How the lower-bit corrections are certified

Choose integer lifts U_k of a binary basis of V_k. Then

    M_k = 2 M_(k-1) + span_Z(U_k),
    M_0 = Z^r.                                        (13)

For if z is in M_k, subtract the binary combination of U_k matching z
modulo two. The difference is 2w with A w divisible by 2^(k-1). The reverse
containment is immediate. By induction, M_k is generated by

    2^k e_j and 2^(k-i) U_i for 1<=i<=k.

If a chosen u is in M_k, any correction preserving its binary reduction
and its previous divisibility can be taken from 2M_(k-1). Its possible
next binary images are generated by

    A U_i / 2^i modulo two,    0<=i<k,

where U_0 is the standard basis. Solve the ordinary binary equation for
A u / 2^k in these image generators and subtract the corresponding
integer correction. Then **directly verify** divisibility by 2^(k+1).
This is exactly the correction mechanism implemented, rather than the
invalid naive lift warned against in section 3.

At n8, images from U_0,U_1,U_2 already have rank 3526 of the 3527 reduced
coordinates. They suffice to make every next correction in the surviving
one-dimensional radical until its final nonzero q_8. The checker still
checks existence and every integer image at each step, and has an explicit
maximum precision that fails rather than silently certifying an unfinished
answer. Here it terminates at k=8 before that bound.

Additional reproduction:

    python check_electrical.py --n8 ../session38-wiring
    python check_n8_lift.py
    python lift_final_radical.py
    python complete_n8_tower.py

Repeat with `python -O`. The sequence intentionally retains the intermediate
mod-four and mod-eight certificates as well as the terminal certificate.
For publication in a different directory, pass the full ledger path using
`--base` on all three lift scripts, or retain the documented sibling layout.
