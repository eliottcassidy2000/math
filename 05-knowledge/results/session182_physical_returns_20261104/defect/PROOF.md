# Evolving joint ports and the transitive-clone ternary obstruction

2026-11-04. **PROVED all-order transport and a single-branch obstruction; FINITE-EXACT verification.** No repeatable two-digit grammar or global spectrum theorem is claimed.

## 1. Scope and genuinely new boundary

The starting point is Session181's named endpoint state and its multiport response `(h,F,J)`. Earlier zero-ready structure results already supply a general ordered path-cover interface and an affine closure criterion. The explicit simultaneous update below is a specialization/compiler refinement, not a claim to have invented that interface. The new structural boundary is §4: even a single persistent `+4` branch cannot be obtained by repeatedly lengthening one transitive vertex module, regardless of the size of the fixed surrounding tournament. This is stronger in that direction than a finite-template counting obstruction against closing both digits.

All physical old names and incidences are retained in births. A change in a port set is an explicit change of membership, not a scalar reset. If an additional finite-field phase is present, retain/recompute it from its actual defining endpoint data; no phase equivalence is used in this packet. Native complement programs need their own literal legality witnesses, which these birth and clone identities do not supply.

## 2. Simultaneous evolving-port transport

Let C be a tournament on the named set V. In the squarefree ring `Z[z_v:v in V]/(z_v^2)`, let `G_ij` enumerate nonempty simple directed paths from i to j, including singleton paths. Set

    e_j=sum_i G_ij, s_i=sum_j G_ij, P=sum_ij G_ij,
    E=sum_{w in W} e_w, S=sum_{w in W} s_w,
    h=[V]P, F=[V]E(1+P), J=[V]ES, alpha=h-J.

The port set W may contain any actual old vertices; the empty set is allowed algebraically. For the canonical strong rooted realization require C strong and W nonempty. Add a new named vertex x. Let `a_i=[i -> x]`, `q_i=[x -> i]`, and choose `delta in {0,1}` to declare whether x is a new port. Define, entirely on the old state,

    X=sum_i a_i e_i, Y=sum_i q_i s_i, A=1+X, B=1+Y,
    C_W=sum_{i,w in W} q_i G_iw,
    D_W=sum_{w in W,j} G_wj a_j.

Here the first sum means i ranges over V and w over W. The full polynomial updates are

    P'=P+z_x AB,
    E'=E+z_x A(C_W+delta),
    S'=S+z_x(D_W+delta)B.

Consequently the new scalars are exactly

    h'=[V] AB,
    F'=[V]{ A(C_W+delta)(1+P)+EAB },
    J'=[V]{ A(C_W+delta)S+E(D_W+delta)B },
    alpha'=h'-J'.                                      (1)

Proof: a path visiting x has a unique incoming and outgoing old piece, either of which can be empty. For a path ending at an old port the outgoing piece is nonempty and ends there; for a path ending at x it is empty. This gives E'; reversing the role of starts and ends gives S'. Multiplication and `z_x^2=0` yield (1). This proof applies after any finite sequence of updates, with the current graph and current ports.

For completeness the full named matrix updates are

    G'_ij=G_ij+z_x (sum_k G_ik a_k)(sum_l q_l G_lj),
    G'_ix=z_x sum_k G_ik a_k,
    G'_xj=z_x sum_l q_l G_lj, G'_xx=z_x.

Thus (1) is an exact unbounded-state grammar when this full matrix is retained. It is not a closed recurrence on only h,F,J.

### Each individual joint port and inner defect

Define `J_wv=[V]e_w s_v`. Write `C_w=sum_i q_i G_iw` and `D_v=sum_j G_vj a_j`. For old w,v,

    J'_wv=[V](A C_w s_v+e_w D_v B),
    J'_xv=[V]A s_v, J'_wx=[V]e_w B, J'_xx=0.           (2)

The defect at every current physical vertex is recovered with its current actual in-neighbors:

    b'_w=F'_w-h'=sum_{v -> w in C+x} J'_wv.             (3)

Equation (3) is the bad-join bijection: a spanning ordered two-slot cover with first slot ending at w either joins to a Hamilton path or has a second slot starting at an in-neighbor v. In particular the term involving x appears for an old w exactly when x -> w. The newborn's own defect sums (2) over the old vertices beating x. This carries all inner defects together rather than only an aggregate port defect.

### Strongness and exact branch predicates

If C is strong, then C+x is strong if and only if x has at least one in-neighbor and one out-neighbor in C. Necessity follows from a source/sink obstruction. For sufficiency, old vertices reach x through any old in-neighbor; x reaches every old vertex through any old out-neighbor and strongness of C. This also correctly excludes both births from the singleton source.

For the rooted family with a physical root r, `TT_s -> C`, r beating the prefix and C minus W, and W beating r, Session181 gives

    H_s=alpha+2^s F, E_r(2)=h+2^s F, b_r=J.

It follows that a birth operation realizing a clock-uniform digit d must satisfy, and at the response level is equivalent to,

    F'=3F, alpha'=3alpha+d.                            (4)

Test (4) with (1), retain J' and (2), and apply the strongness guard. This provides a precise evaluator and target, not existence of a successful cut. In particular the actual Session181 multiport candidate has `(h,F,J)=(18583,124383,3690)`, so a next uniform branch still requires `F'=373149` and `alpha'=44681` or `44683`. Its J' cannot be discarded.

## 3. Exact unbounded transitive-clone spectrum

Fix any tournament C and one physical vertex x. Replace x by a transitive module of length t>=1, with every clone having x's old incidences to `U=C-x`. For t=1 identify its single member with x. A literal nesting can keep x and add later clones after all existing clones, so no old edge or name need change. Call the result C_t.

In the squarefree ring on U, write I for old in-neighbors of x and O for old out-neighbors. Set

    L=1+sum_{j in I,i in U} G_ij,
    R=1+sum_{i in O,j in U} G_ij,
    M=sum_{i in O,j in I} G_ij,
    A_r=[U] L M^(r-1) R.

Then for all t>=1,

    H(C_t)=sum_{r>=1} A_r r! S(t,r),                   (5)

where S(t,r) is a Stirling number of the second kind. Only `r<=floor(|U|/2)+1` can contribute, and

    A_r >= 0, (r-1)! divides A_r.                     (6)

Proof: a Hamilton path uses r nonempty runs of clones. The old piece before the first run is empty or ends in I; after the last run it is empty or starts in O. Each of the r-1 intervening old pieces starts in O and ends in I, so it has at least two vertices. Their ordered cover count is precisely A_r. Assigning the t individually named clones onto r ordered nonempty runs gives the number of surjections `r! S(t,r)`; within each run the transitive order is forced. Conversely these data reconstruct the path uniquely. Permuting the r-1 distinct disjoint middle paths acts freely on these covers, proving divisibility (6).

If C is nontrivial strong, every C_t is strong: paths in the quotient lift between modules, and two different clones can reach one another through an old directed cycle through x. Alternatively use the strong birth guard at each clone addition. For a singleton C the clone family is transitive and nonstrong at t>1.

Formula (5) is a finite exponential spectrum because

    r! S(t,r)=sum_{k=1}^r (-1)^(r-k) binom(r,k) k^t.

It explicitly exhibits genuine base-3 terms, not just binary growth. For example

    C rows=[20,1,2,7,14], clone physical vertex x=0,
    (A_1,A_2,A_3)=(9,10,2),
    H(C_t)=2*3^t+4*2^t-5.                             (7)

The additional base-2 mode does not vanish merely because a base-3 mode exists.

There are also exact pure-ternary positive controls. With strong rows

    C=[26,24,67,116,36,71,19], clone physical vertex x=6,
    (A_1,A_2,A_3,A_4)=(145,150,50,0),
    H(C_t)=50*3^t-5,
    H(C_{t+1})=3H(C_t)+10.                            (9)

Thus the obstruction below is not a blanket obstruction to an affine ternary clone mechanism. The fixed digit in (9) is 10, not either requested digit. Replacing H by H+4 algebraically would yield a digit2 recurrence, but that is not a physical construction or an admissible count reset.

## 4. No persistent +4 by lengthening a fixed transitive module

**Theorem.** No fixed tournament C and physical vertex x can satisfy

    H(C_{t+1})=3H(C_t)+4

for every t>=1, or even for every t beyond a finite threshold. More generally if the fixed digit d is 2 or 4 and that recurrence holds eventually, then it holds for all t>=1 and necessarily

    A_r=0 for r>=4,
    A_2=3A_3,
    A_1=3A_3-d/2.                                    (8)

The d=4 case is impossible. More generally any fixed integer digit d in such an eventual recurrence must satisfy d congruent to 2 modulo4. Conditions (8) are necessary and sufficient for the algebraic response recurrence with a given d; this statement does not assert that a tournament realizes them when d=2.

Proof: a recurrence solution has form `c*3^t-d/2`. Distinct positive-integer exponentials are linearly independent even on a tail (a finite Vandermonde matrix proves this). If R is the largest index with A_R nonzero, the coefficient of R^t in (5) is A_R. Thus R<=3. With those remaining terms,

    H(C_t)=A_3*3^t+(A_2-3A_3)*2^t+(A_1-2A_2+3A_3).

Comparison gives (8), including when some A_r vanish. By (6), A_3 is even. If d=4 then `H(C_t)=A_3*3^t-2` is even, contradicting the classical Rédei parity theorem that every tournament has an odd number of Hamilton paths. Equivalently (8) makes A_1=H(C) even. For any integer d, the same comparison gives d/2=3A_3-A_1, which is odd because A_3 is even and A_1=H(C) is odd. Hence d is 2 modulo4. The parity theorem is used only at this last step; see Richard Stanley's [Rédei theorem slides](https://math.mit.edu/~rstan/transparencies/iceca.pdf) or the primary research discussion [The Tournament Theorem of Rédei revisited](https://arxiv.org/abs/2510.10659).

This excludes the strategy of repeatedly adding one transitive twin of the same fixed vertex for the persistent +4 branch, even on an arbitrarily large fixed context and even allowing a finite initial transient. It does not exclude changing the module chosen at each step, changing external incidences, growing nontransitive modules, nonlinear reparametrization of module length, or other named operations. Applying the theorem to a whole rooted tournament is legitimate if all other vertices and incidences are fixed; it does not prove (4) impossible for general births.

## 5. Exact controls

`verify.py` has no dependencies or removable assertions. It constructs squarefree endpoint matrices by ordered-tuple enumeration and checks new scalar responses against independent complete-permutation enumeration. It exhausts all labelled tournaments of order1 through4, every birth cut, every old port subset W, and both newborn memberships. It also checks every individual joint coefficient, each actual inner defect, and the full matrix update including new endpoints. Strongness is checked by literal reachability.

An independent subset Hamilton-path DP checks (5) for all these tournaments, every cloned physical vertex, and t=1,...,4. It verifies the divisibility in (6), strongness when applicable, and the nontrivial three-run example (7) through t=7. The pure-ternary +10 example (9) is independently DP-checked through t=6. Normal and optimized runs have identical outcomes. These finite checks audit the all-order proofs; they are not used as evidence of a successful two-digit grammar.
