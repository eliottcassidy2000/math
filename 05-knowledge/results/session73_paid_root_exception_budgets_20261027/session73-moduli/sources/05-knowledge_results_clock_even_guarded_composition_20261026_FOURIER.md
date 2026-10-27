# From equinumerosity to guarded response operators

## F1. The inherited basis, with the correct symmetry group

Let V be the binary edge space on the named vertices [n]. Choose the upward
orientation only as a coordinate gauge. For a tournament T write x_e(T)=+1
for an upward arc, -1 otherwise, and chi_F(T)=product_(e in F) x_e(T).
These 2^binom(n,2) characters form the full Walsh basis.

If G is a group of allowed vertex relabelings, its action on the characters
permutes their graph supports and inserts the edge-reversal sign. A support
F contributes one invariant orbit sum exactly when every g in Stab_G(F)
has reversal sign +1. Otherwise the putative coefficient equals its negative
and vanishes. Different graph orbits have disjoint character supports, so
the surviving signed orbit sums form a basis for all G-invariant functions.
This is the elementary monomial-representation argument in
[THM-4721 / twisted Burnside and Walsh basis](../../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md).

For G=S_n the surviving graphs are automorphism-sign-even. This is the
[Royle--Praeger--Glasby--Freedman--Devillers counting theorem](https://arxiv.org/html/2204.01947v2)
in linear coordinates. Their Theorem1.1 equates unlabelled counts at the same n;
it does not identify labelled objects. In particular this is not the binary
cycle space of even-degree graphs. The signed-orbit proof above also applies
to G fixing ordered ports: the relevant graph predicate is then G-even, not
evenness under the larger full automorphism group. At n4, two named ports
give32 basis lines and32 rooted tournament types, compared with four unrooted.

The prior [SCC/Lyndon codec](../tournament_even_bijection_20261011c.md) is a
different exact bridge: it uses explicit rankings of irreducible types. Its
transported ordered join is conditional Lyndon fusion, not ordinary disjoint
union. Its induced graph moves can depend on the source state. This distinction
is required by the [all-order equivariance obstruction](../../../01-canon/theorems/THM-4722-no-equivariant-bijection-tournaments-even-graphs-for-every-n-sylow-hall-violator.md).

## F2. A path legality guard becomes Fourier convolution

For an ordered Hamilton word P on the named vertices, let E_P be its n-1
unordered path pairs, and let epsilon_e be +1 if its directed use is upward,
-1 otherwise. Define

    g_P(T) = product_(e in E_P) (1+epsilon_e*x_e(T))/2,
    tau_P(T) = T with every pair outside E_P reversed,
    (K_P f)(T) = g_P(T) f(tau_P(T)).

Thus illegal words have weight zero. One can instead keep an explicit failure
state when failure itself is an observable. The path edges are unchanged,
so g_P(tau_P T)=g_P(T) and K_P^2 is multiplication by g_P.

**Theorem F2 (PROVED, all n).** On every character chi_F,

    K_P chi_F = (-1)^|F\E_P| / 2^(n-1)
                 * sum_(A subset E_P) (product_(e in A) epsilon_e) chi_(F triangle A).

**Proof.** Complement multiplies chi_F by the displayed sign. Expand the
product for g_P, then use chi_F chi_A=chi_(F triangle A). QED.

This is a precise extension beyond the source-independent convolution
operators in [Session72 Theorem A](../session72_even_graph_unifying_carriers_20261026/session72-even/EVEN_GRAPH_OPERATION_ENVELOPE.md).
Those are diagonal in Fourier coordinates. Multiplication by a nonconstant
guard mixes supports. A named P generally fails to preserve unrooted invariant
functions; retain P's marking or sum over its relabeling orbit before using
an unrooted matrix. Forgetting that qualification is exactly an alignment loss.
General guards have the same formula with their Fourier coefficients replacing
the explicit path coefficients. Bounded-order marginals are closed only if
every allowed operation preserves their span, which the Vitali conjunction
family does not do.

## F3. The actual counted operator and a finite matrix

Let L=sum_P K_P over all ordered permutations. It is relabeling-equivariant,
so it descends to tournament types. Its entry L_(T,U) counts actual Hamilton
paths from the chosen labelled T whose partner has type U. It is not a Boolean
adjacency matrix, and includes loops and repeated targets.

Let M_(T,F)=Phi_F(T) be a signed-orbit evaluation matrix. Then

    B = M^(-1) L M,          M B = L M.

This is an exact response representation; compositions obey
M^(-1)(L_1 L_2)M=B_1 B_2. An analogous rectangular formula applies when the
input and output orders differ, with their respective bases.

At n4, graph-column order is empty, P3+isolated, K1,3, C4. With the source
mask convention bit1 means a downward arc, the source rows are0,2,4,5:

    M = [1,  4,  0,  1]        B = [ 3, -6, 0,    0]
        [1,  0, -4, -3]            [-1/2,0, 0,  1/2]
        [1, -4,  0,  1]            [ 0,  0, 0,    0]
        [1,  0,  4, -3]            [ 0,  2, 0,   -1].

The row2/row4 sign convention is converse to Session72's upward-bit convention.
In particular L1=3-Phi_(P3)/2. Negative entries are not negative probabilities:
for a law p, the response vector is M^T p, with positive cone M^T R_+^types.
The inherited [Walsh realization theorem](../session72_even_graph_unifying_carriers_20261026/session72-root/WALSH_REALIZABILITY.md)
already proves that this cone, not the coordinate orthant, is the correct
realization test. This packet supplies the new path-guard matrix as a control.

**Detailed balance (PROVED).** For labelled T,U, the move indexed by P is
reversible by the same P. Consequently

    |orb T| L_(T,U) = |orb U| L_(U,T).

The normalized operator P_(T,U)=L_(T,U)/H(T) has reversible weights
|orb T| H(T). Thus the counted and normalized operations require different
weights. Presentation-response coefficients additionally divide the source
row by |Aut T|. None of these conventions can be exchanged silently.

## F4. Composition at changing size

For a fixed ordered partition I,J with every cross arc I->J, a raw character
restricts as

    chi_F(T_I -> T_J) = chi_(F restricted I)(T_I) chi_(F restricted J)(T_J),

because each increasing cross sign is +1. With another label ordering retain
the fixed cross-sign factor. Therefore pulling back any signed-orbit sum
along ordered join gives an exact linear combination in the tensor product
of the two component invariant bases. This proves a finite structure tensor
at every pair of sizes; it does not identify graph union with tournament join.

More generally, for a finite weighted gluing rule A from typed input tuples
to typed outputs, the tensor in character coordinates is

    B_A = (M_1 tensor ... tensor M_r)^(-1) A M_out.

Associative substitution contracts these tensors associatively whenever the
underlying gluing has a bijection between the two parenthesizations of its
complete witnesses, preserving weights. Work with named incidences first;
quotienting independently before gluing can discard relative frame or
stabilizer multiplicity. A generic prime quotient can revisit a child, so
its state still needs every ordered run component, not two outer endpoints.
Fourier coordinates encode that retained carrier; they do not remove its data.

## F5. What the finite ranks say

Converse carries each path P to its reverse and commutes with L. Hence all
L^k1 are converse-invariant. The ranks 1,2,3,10,34 at n2..6 equal the respective
converse-orbit counts. This is an exhaustive finite calculation over Q, not
an eigenvalue fit. An independent C++ enumeration extends this to n7: rank272
modulo the verified prime 1000000007 supplies the lower bound over Q, and
converse invariance supplies the upper bound272. Conjecture C in README asks
whether it persists.

Since L is self-adjoint for positive orbit weights, C is equivalent to L
having simple spectrum on the converse-invariant subspace and the constant
vector having nonzero projection on every eigenspace there. This equivalence
is elementary spectral decomposition; neither spectral condition is proved
uniformly. The observed ranks do not imply that every other observer, every
normalized kernel, or the full blue attachment grammar has the same quotient.

## F6. Named native-complement test-return is complete at every size

This is a direct application of incoming
[Session72 response Theorem C](../session72_even_graph_unifying_carriers_20261026/session72-response/REPORT.md),
not a new abstract test-return theorem. It makes the restoration operation
literal for native complement on the **full labelled tournament cube**.

Adjoin an absorbing failure state and let a named operation test P, perform
tau_P if legal, and fail otherwise. F2 gives K_P^2=M_(g_P), multiplication
by g_P. Consequently a sequence of doubled named operations has terminal
survival response

    (K_(P1)^2 ... K_(Pr)^2)1 = product_i g_(Pi).

The actual operation returns to the same labelled source after each successful
pair, so the product is a joint test, not a product of separate expectations.

By inherited THM-4807, the complete Hamilton-word legality signature determines
every labelled tournament. For each T, its signature-cell indicator is

    1_{T} = product_(P legal in T) g_P
              * product_(P illegal in T) (1-g_P).

Expand the second product. Every summand is a terminal response of a doubled
operation word (the empty word gives 1). Thus these responses span the entire
2^binom(n,2)-dimensional function space on labelled tournaments, over every field. This proves
full linear law recovery on the valid states, for every n>=1. The indicator
expansion uses signed linear combinations of experiments; it is not a claim
that one positive survival experiment isolates every tournament.

On a fixed-path subcube, the same observation separates initial states, but
the experiment may pass outside that subcube between the two moves. Closure
therefore requires permitting those intermediate marked states. The incoming
re-marking/inverse protocol stays within its own fixed-path carrier and has
rank 2^((n-1)(n-2)/2); these are different declared universes.

The new finite control constructs the full guard algebra through n4 over F2.
This named all-size statement does not prove Conjecture C: replacing each
named K_P by their sum L averages away the ability to choose a specific test.
It also makes no black-only accessibility claim.
