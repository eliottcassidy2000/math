# The conic square hides the arrow: Pythagorean geometry, Petersen, and a Paley target chart

**Status:** proved elementary dictionary and exact obstruction, with finite graph reconstructions. The geometric ingredients are inherited, and no novelty for their classical identifications is claimed. This does not solve the character target equation or Chan–Pak.

## Inheritance

`math/05-knowledge/results/pythagorean_form_conic_graphs_rp_margin_20261018.md` connects the integral Pythagorean form with its finite-field conic, Petersen, five exterior triangles, K5 contractions and K3,3 subdivisions. It credits Coxeter's classical conic graphs and records its audit corrections. `THM-4602-tournaments-as-chirotopes-pfaffian-criterion-and-legendre-frieze-strips.md` already identifies the Legendre determinant chart with Paley+sink and explains switching under representative rescaling. `paley_target_bridge/PALEY_TARGET_REFORMULATION.md` supplies the bounded universal Paley host and exact integer target equation. The present note isolates the precise square-root information lost between these carriers and its impact on using A5 symmetry for target selection.

## 1. One quadratic form, two levels of information

Work over any field F of characteristic other than2. Put

    Q(X,Y,Z)=X²+Y²−Z²,
    B(v,w)=X_vX_w+Y_vY_w−Z_vZ_w.

For spinors a=(u,v), define

    V(a)=(u²−v²,2uv,u²+v²).

Then Q(V(a))=0, and direct expansion gives the exact polarization identity

    B(V(a),V(b))=−2 det(a,b)².                    (1)

In the affine chart a(t)=(1,t),

    V(t)=(1−t²,2t,1+t²),
    B(V(t),V(s))=−2(s−t)².                       (2)

In particular V(a)=V(−a), but replacing a by−a reverses the quadratic-character sign of every determinant involving a when q≡3 mod4. The same conic vector therefore admits two lifts with opposite tournament contacts.

These formulas are valid for rational Pythagorean parametrization as well as finite fields. They do not identify a tournament Hamilton count with a Pythagorean norm.

For F_q with q odd and q≡3 mod4, let χ be its quadratic character. For distinct projective points in this affine chart, equation(2) implies

    χ(B(V(t),V(s)))=χ(−2),                       (3)

which is **constant**. The Pythagorean conic's polar pairing has squared away precisely the signed determinant that distinguishes the two arc directions.

This is a strict information-loss statement, not merely an unsuccessful example. B is symmetric under exchanging t,s, whereas χ(s−t) changes sign because χ(−1)=−1. No rule using only B on an unordered pair can recover the Paley direction. Independent rescaling of projective conic representatives can introduce factorsχ(λ_t)χ(λ_s), but the resulting pair datum remains symmetric and still cannot orient the ordered pair by itself.

## 2. The constructive repair is a signed lift, not a new numerical identity

Retain a representative spinor a_i=(1,x_i), rather than only its quadratic conic image. Define

    i→j iff χ(det(a_i,a_j))=χ(x_j−x_i)=+1.

This is exactly the physical induced Paley tournament used in the target reformulation. For m distinct coordinates its Hamilton count is

    2^(1−m) sum_{π∈S_m} product_{j=1}^{m−1}
                 (1+χ(det(a_{π(j)},a_{π(j+1)}))).       (4)

The sum is in the ordinary integers. Thus the character target equation is equally a **signed spinor configuration problem on a conic**, provided the representatives are retained. Every contribution corresponds to a literal Hamilton path.

An SL2 transformation preserves all determinants before representatives are renormalized. Returning to the affine chart introduces per-vertex scaling factors. Their square/nonsquare classes implement vertex switching, exactly as in THM-4602. This is why the companion Möbius-switching route is appropriate: conic geometry provides switching operations, not a free identification of all configurations or all counts.

Selecting x_i to make(4)=N is still the unsolved physical target problem. Replacing det by its square would make all off-diagonal character tests constant and would not be an equivalent equation. Real signs also cannot simply replace finite-field quadratic characters: the real rank-two determinant patterns are locally transitive, whereas the universal finite-field host contains every bounded tournament.

## 3. What the A5/Petersen/K5/K3,3 connection really preserves

At p=5, consider projective points of Q over F5. The conic has6 points. Among the nonconic points, the10 with nonsquare Q have orthogonality graph B=0 isomorphic to Petersen. The15 with square Q have graph5 disjoint triangles. These facts are recovered from the prior packet, not newly discovered.

An independent reconstruction in the replay supplies an explicit Petersen isomorphism to the graph whose vertices are2-subsets of five labels, adjacent when disjoint. This identifies the natural A5 action. The replay also verifies:

* its six perfect matchings each contract to K5;
* deleting any Petersen vertex and suppressing its three degree-two vertices gives K3,3;
* the even permutation(01)(23) is a nontrivial involutory automorphism on those2-subsets.

These are actual maps between undirected graphs. They preserve incidence in the specified way, not Hamiltonian path counts under tournament completion. Contracting an edge or suppressing a degree-two vertex does not preserve the number of spanning paths in general.

There is an exact obstruction to keeping this A5 action as tournament symmetry. A nontrivial involution on a finite set exchanges some pair of vertices; if it were a tournament automorphism, it would send the edge between that pair to its reverse. Hence no faithful even-order group acts by tournament automorphisms. The natural A5/Petersen action therefore cannot survive unchanged on any tournament orientation/completion of those ten labels. The same objection applies to any claimed full A5 tournament symmetry, not only this graph.

At p=5 there is an additional character obstruction:χ(−1)=+1, so the determinant Legendre test is symmetric, not an orientation rule. Passing from that conic to a Paley tournament requires a different field with q≡3 mod4 and a chart/sign choice. A repeated occurrence of five cannot remove that distinction.

## 4. A concrete design lesson for the selector

The full conic geometry can organize possible points and transformations, but its quadratic pairing is too coarse to select the count. The useful interface is the **oriented square root**, namely the signed spinor determinants and their normalization factors. A construction that retains only norms, projective incidence, or the A5 orbit loses the data that equation(4) reads.

For an algorithm, this gives a justified separation:

1. Use conic/Möbius operations to generate point sets and correlated switches.
2. Carry the actual signed representatives or their induced arc pattern throughout.
3. Evaluate the exact target residual using(4) or subset DP.
4. Do not quotient by the larger projective or A5 action as though it preserved H.

The universal Paley host ensures existing witnesses can be represented in this language, not that a residual-reducing move always exists. The new question is whether a controlled sequence of signed lifts can change the count to N while preserving the continuation information. Neither the conic identity nor the graph contractions establish that descent.

## Reproduction

Run `python conic_sign_replay.py`. It independently verifies(1)/(2) on affine controls over5,7,11,19,23, builds all projective points at5, finds an explicit Petersen isomorphism, checks all six K5 contractions and ten K3,3 subdivisions, and verifies the displayed A5 involution. Normal and optimized output are identical. The field identities and symmetry obstruction have all-order proofs above; the finite graph identifications are explicitly bounded to p=5.
