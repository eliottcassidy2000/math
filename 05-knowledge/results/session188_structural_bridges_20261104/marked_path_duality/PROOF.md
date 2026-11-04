# An exact marked path-duality identity over the integers

## 1. Setting and provenance

Let A be any finite tournament on fixed physical names. Work in the commutative squarefree ring Z[z_v:v in A]/(z_v^2). Let G_pq be the sum of all directed p-to-q paths, including a singleton when p=q, weighted by their vertex supports. Put

    C_p=sum_u G_up,  R_p=sum_v G_pv,  H=sum_(u,v) G_uv,
    E_pq=G_pq-G_qp,  D_pq=C_p R_q-C_q R_p.

Let E_pq^odd retain just the monomials with odd vertex-support size. Then, coefficientwise over the integers,

    D_pq = -2(1+H) E_pq^odd.                             (1)

Every product retains disjoint physical supports. This is not a product of independent endpoint counts on overlapping sets.

The squarefree path algebra and its converse/parity identities are inherited tools in the pinned E1/P2 work and the Session188 response packets. The purpose here is to give a short self-contained derivation of this exact marked integer specialization and one response-kernel application. No historical priority claim is made for the underlying reciprocity or rank-one identity. Equation(1) alone is not a non-annihilation or reconstruction theorem.

## 2. Rank-one converse calculation

Write B for the zero-diagonal tournament adjacency matrix, Z=diag(z_v), and J=11^t. Then B+B^t=J-I and

    G=(I-ZB)^(-1) Z.

The inverse is a finite squarefree series; every repeated vertex vanishes. The matrix (I-ZB^t)^(-1)Z is G^t. Define N=(I+ZB)^(-1)Z=-G(-z), and h=1^t N1. Because Z^2=0,

    I-ZB^t=(I+Z)(I+ZB)-ZJ,
    ((I+Z)(I+ZB))^(-1)Z=N.

The rank-one inverse formula, valid since every denominator has constant term one, gives

    G^t=N+(N1)(1^t N)/(1-h).                            (2)

Summing all entries yields H=h/(1-h). Since h=-H(-z), this proves the exact reciprocity

    (1+H(z))(1+H(-z))=1.                               (3)

Row and column sums of (2) give

    C_p(z)=-(1+H(z))R_p(-z),
    R_p(z)=-(1+H(z))C_p(-z).                            (4)

Subtract the (q,p) entry of (2) from the (p,q) entry. Its left side, after moving N, is E_pq(-z)-E_pq(z). By (3)-(4), its right side is D_pq(z)/(1+H(z)). Therefore

    D_pq/(1+H)=E_pq(-z)-E_pq(z)=-2E_pq^odd,

which proves (1). No division by two in the coefficient ring is used: the odd-support projection is defined directly.

## 3. A fixed-name arbitrary-child quartic specialization

Adjoin two old singleton names u,v to A, with A->u->v->A, and let T be this cyclic substitution. Let W converse only A. Add the actual ports x->y, with a_w=[w->x] and c_w=[w->y]. For distinct physical p,q in A define the disjoint quartic

    Delta_pq=[a_p a_u c_q c_v](F_T-F_W).

Then the complete exact coefficient is

    Delta_pq=[A]{(1-H)D_pq+(2H+H^2)E_pq}               (5)
            =[A]{2(H^2-1)E_pq^odd+(2H+H^2)E_pq}.     (6)

Here [A] extracts the full A-support after squarefree multiplication. H is the path series, not merely the scalar Hamilton count.

For an explicit check of (5), use the complete eight-term disjoint quartic normal form in ../cyclic_child_disjoint/PROOF.md. Enumerating the placements of u,v in its three old path slots gives the source coefficient before full A extraction as

    C_p C_q+C_p R_q+G_qp+(1+H)^2G_pq
                 +H C_q R_p+R_p R_q.

Converse of A interchanges C and R and transposes G while fixing H. Subtracting cancels the symmetric two-column/two-row terms and leaves exactly (5). Equation(1) then yields (6). All six named old/port roles remain physical; no independent relabeling or profile averaging is used.

This corollary is specifically for the two-singleton cyclic exterior. Equation(1) holds for every A, but no arbitrary-exterior positivity, arbitrary ambient sign transfer, or full named-response equality follows from it.

## 4. Useful boundary

For a transitive pair A, all these disjoint quartics vanish. For the seven-name Paley tournament with arcs u->u+d modulo7 for d in{1,2,4}, the retained calculation gives every off-diagonal Delta_pq equal to +8 or -8. Thus a universal nonzero-modulo-eight assertion for this marked kernel is false, despite its exact nonzero integer values. This is only a boundary of this chosen coefficient test, not a full-response collision.

The verifier checks coefficientwise identities on every labeled tournament through five names, additional named larger controls, literal complete fourth differences for small arbitrary children, and a complete named Boolean coefficient bank for the Paley example. All cubic/quartic coefficients and zeros are retained. These finite checks support, rather than replace, the all-order algebraic proof.
