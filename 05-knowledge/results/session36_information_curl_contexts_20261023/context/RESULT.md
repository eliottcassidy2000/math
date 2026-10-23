# Two exact context quotients for Hamiltonian paths

Status: PROVED elementary reconstruction and rank theorems; FINITE-EXACT checks below. Repository-relative advance, not a literature-priority claim. Written during Session36 by recovering the March/August path-polynomial lineage, not extrapolating a finite census.

## Later completion

`CURL-AND-CLAN-COMPLETION.md` proves a stronger exact quotient classification after recovering the primary C3-structure theorem: on labelled tournaments, `W^+` equality, full `W` equality, equality of the cyclic-triangle hypergraph, and equality of every vertexwise-substitution Hamiltonian response are equivalent. A positive full-support coefficient detects each triangle. This supersedes the deliberately cautious unresolved wording in Section6 below, without changing TheoremsA-C or the exact module rank.

## 0. Main result and the distinction that matters

There are two different places a tournament can occupy in substitution.

1. A **module** is one indivisible interface block: every exterior vertex meets all of its vertices uniformly. Its exact Hamiltonian context state is the full spanning path-cover profile. At fixed order `n`, the connection matrix for all one-hole module contexts has rational rank exactly `ceil(n/2)`.
2. A **quotient** supplies the individually addressable slots to be filled by modules. Its exact Hamiltonian context state is the full-support part of its multivariate directed-walk series. Replacing each vertex by a transitive tournament already supplies a complete family of probes for this state, via multivariate Stirling inversion.

Neither state reconstructs labelled adjacency. Neither applies to arbitrary vertex attachments without extra boundary information. The older six-vertex twins in THM-3235 prove that the two states are genuinely different, not competing descriptions of one context category.

This makes the Session35 slogan precise in this lane: information is sufficient only relative to a specified interface and a specified observable.

## 1. Inheritance, provenance, and live concept board

Closest proved mechanisms:

- `math/01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`: full-profile substitution via factorial contraction against `exp(y W_Q)`. This repairs the merely verified/open kernel language in THM-1975.
- `math/01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`: path-cover profile equals the endpoint jet of a parity polynomial `E_T`.
- `math/01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md` and `THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md`: injective path-colour polynomial `Q_T`, multiplicativity, pure order parity, half-grid reconstruction. THM-3181 explicitly makes no sharp minimality claim on realized tournament polynomials.
- `math/01-canon/theorems/THM-059-universal-master-polynomial.md`: much older central-factorial/master-polynomial view. Its `W` is the same centered permutation polynomial as `E`, with notation changed.
- `math/01-canon/theorems/THM-250-flip-formula.md`: the single reversed arc is the canonical cheap near-transitive test. Here its whole path-colour polynomial becomes an explicit full-rank family, instead of only its Hamiltonian value.
- `session35-recursions/CONTEXT-AND-PORTS.md`: square-zero simple-path series and logarithmic odd-cycle identity. This is a transparent route to the coloured odd-cycle formula used below, not a new OCF claim.

Canonical hostiles:

- THM-3134's order-five masks40/76: equal `H=15`, different first derivative; the one-vertex-per-side cyclic context has `H=123/135`.
- `math/01-canon/theorems/THM-3235-transitive-tournament-blowup-saddle-and-hypergeometric-decimation.md`: order-six masks408/1332 have equal complete profile, but different vertexwise transitive blowups. A module invariant cannot be silently promoted to a quotient invariant.
- `math/01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md`: port score/curl census does not fix wiring at order six. Its exact/finite scope is retained; no all-order metagraph reconstruction is inferred here.

Concept board:

- Scalar `H`: closed output, loses run count. Probe by one cyclic context.
- `pc/E/Q`: same module state in three bases. Probe completeness and rank.
- `W_Q`: quotient walk content, resolves vertex slots. Probe by Stirling inversion.
- Square-zero ports: preserve exclusion labels and endpoints; exact but exponential.
- Metagraph wiring: needs transport of a chosen path/ports, not merely a node scalar.

The recovered `complete_invariant_search.out` is historical search provenance only. Its sorted overlap/walk collisions support the methodological warning but are not dependencies of the proofs here. Likewise historical formula-versus-#P rhetoric in THM-1945/1975 is not used as a complexity theorem.

## 2. Definitions and exact context category

For a finite nonempty tournament `T`, let `p_T(d)` count spanning covers by `d` unordered vertex-disjoint directed paths, singleton paths allowed. Let `F_T(d)=d!p_T(d)`, and set it to zero outside `1<=d<=|T|`.

A one-hole module context is built from fixed nonempty tournaments using tournament substitution and exactly one occurrence of a distinguished hole. The identity context is allowed. Thus the varying tournament is always inserted as one module, possibly inside successive larger modules. It never supplies the quotient's individually addressable vertices, and is never duplicated.

For fixed `n`, define the connection matrix

    M_n(T,C)=H(C[T]),

with rows all `n`-vertex tournaments and columns all such one-hole contexts. Rank means rank over `Q`, equivalently the maximum rank of its finite submatrices.

The **exactly one occurrence** condition is important for the linear-rank theorem. If a context uses several copies of the varying tournament, its response can be nonlinear in the profile. Equal profiles still imply equal responses, but the stated rank bound is not asserted for that larger matrix.

## 3. Module contextual equivalence is exactly profile equality

Let `P_m` denote the transitive tournament of order `m`. Define

    C_0[T]=T,
    C_m[T]=C3[T,P_m,P_m] for m>=1,

with the triangle cyclically oriented.

Write `S(m,k)` for a Stirling number of the second kind and

    G_m(k)=k! S(m,k),

extended by zero outside `1<=k<=m`. Every subset of a transitive tournament has exactly one directed increasing path, so its path-cover counts are precisely `S(m,k)`.

THM-3121 gives the quotient walk series

    W_C3(x,y,z)=(x+y+z+xy+yz+zx+3xyz)/(1-xyz).

The coefficient of `x^a y^b z^c`, for positive indices, is three when all indices agree, one when their maximum and minimum differ by one, and zero otherwise. Consequently

    H(C_m[T])=sum_(k=1)^(m+1) a_(m,k) F_T(k),             (1)

where

    a_(m,k)=3G_m(k)^2+2G_m(k)G_m(k+1)+G_m(k+1)^2
             +G_m(k-1)^2+2G_m(k-1)G_m(k).              (2)

In particular

    a_(m,m+1)=(m!)^2 != 0.                              (3)

For `m=0`, use `H(C_0[T])=F_T(1)` separately. Thus the responses to `C_0,C_1,...` recover `F_T(1),F_T(2),...` successively. For example

    H(C_1[T])=3F_T(1)+F_T(2).

**Theorem A.** For finite nonempty tournaments `T,U`, the following are equivalent:

- `p_T(d)=p_U(d)` for every `d`;
- `H(C[T])=H(C[U])` for every one-hole module context `C`;
- `H(C_m[T])=H(C_m[U])` for every `m>=0`.

Proof. The full-profile substitution law is multilinear in the factor profiles. Induction up the unique branch from the hole proves the first implication, including all nested contexts. The second implication includes the listed probes. Equations (1)-(3) prove the last implication by induction on the coefficient index. There is no fixed-order assumption in this equivalence. QED.

Thus the profile is not merely a sufficient sidecar: its equality classes are the largest quotient preserving all the stated observations. The claim is about equivalence classes, not a claim that the usual list of `n` numbers is the shortest encoding.

## 4. Half-size reconstruction from actual Hamiltonian contexts

For a permutation `pi` of `V(T)`, let `b(pi)` count backward consecutive arcs and `g(pi)=n-1-b(pi)`. Put

    E_T(u)=sum_pi (u-1/2)^b(pi) (u+1/2)^g(pi).

Concatenating ordered path-cover components and recording the cuts proves

    sum_d F_T(d) z^d = z E_T(z+1/2).                    (4)

Reversal of the permutation gives

    E_T(-u)=(-1)^(n-1)E_T(u),
    [u^(n-1)]E_T(u)=n!.                                (5)

These are the THM-3134 identities, repeated here to keep the dimension proof self-contained.

Let `r=ceil(n/2)`.

**Theorem B (finite tomography).** If the order `n` is known, the `r-1` integer counts

    H(C_0[T]),...,H(C_(r-2)[T])                         (6)

determine the complete path-cover profile. For `n=1,2`, no measurements are needed because all profiles of that order agree.

Proof. The triangular equations recover `F_T(1),...,F_T(r-1)`, hence the first `r-1` Taylor coefficients of `E_T` at `u=1/2`. If two such polynomials agree on this jet, their difference `D` has a root of multiplicity at least `r-1` at `1/2`. Parity supplies a root of the same multiplicity at `-1/2`. Their universal leading terms cancel, so `deg D<=n-3`. But a nonzero polynomial with both roots has degree at least `2r-2>n-3`. Hence `D=0`. QED.

This differs from the older half-grid theorem in two ways: the measurements are genuine Hamiltonian-path counts of finite tournament contexts, not evaluations at an abstract colour number, and known monicity saves one measurement.

The largest added context in (6) has `2(r-2)` vertices. This is a simple sufficient family, not a vertex-optimal gadget theorem.

## 5. Exact linear rank and an explicit basis of real tournaments

The parity vector space

    V_n={E in Q[u]: deg E<=n-1, E(-u)=(-1)^(n-1)E(u)}

has dimension `r=ceil(n/2)`. Equation (4) identifies it linearly with the span allowed by the path-cover parity identities. Every one-hole response is linear in the varying profile. Therefore `rank M_n<=r`.

To prove equality on actual tournaments, use the injective path-colour polynomial

    Q_T(t)=sum_d p_T(d)(t)_d.

The transform is invertible because the falling factorials form a polynomial basis.

### 5.1 Recovered coloured odd-cycle identity

Let `R=Q[x_v:v in V]/(x_v^2:v in V)` and let `Phi_T=1+sum_(nonempty simple directed paths P) x_V(P)`. For every nonnegative integer `t`,

    [x_V] Phi_T^t = Q_T(t):

choose at most one path for each of `t` labelled colours; square-zero multiplication enforces vertex disjointness, and full support enforces spanning. Both sides are polynomial in `t`, so this holds identically.

The Session35 square-zero determinant identity, equivalently the established OCF mechanism, gives

    log Phi_T = sum_v x_v + 2 sum_(directed odd simple cycles C) x_V(C).

Exponentiating `t log Phi_T` and taking full support yields

    Q_T(t)=sum_(vertex-disjoint odd-cycle collections A)
             2^|A| t^(n-sum_(C in A)(|C|-1)).            (7)

The empty collection contributes `t^n`. Formula (7) is a coloured form of the old odd-cycle collection formula, not a new claim for OCF itself.

### 5.2 The full-rank near-transitive family

Let `B_0=P_n`. For `1<=j<=r-1`, let `B_j` be `P_n` with just the arc `0->2j` reversed. Every directed cycle must use `2j->0`, and then an increasing sequence of vertices strictly between those endpoints. Consequently all directed cycles intersect, and the number of odd cycles of length `2ell+1` is

    binom(2j-1,2ell-1),  1<=ell<=j.

Equation (7) gives the complete profile in a closed form:

    Q_(B_j)(t)=t^n + 2 sum_(ell=1)^j
                        binom(2j-1,2ell-1)t^(n-2ell).   (8)

In the coordinate order `t^n,t^(n-2),...,t^(n-2r+2)`, the `r` polynomials from `B_0,...,B_(r-1)` form a triangular matrix with diagonal `1,2,...,2`. Hence they are linearly independent, with coefficient determinant `2^(r-1)`.

### 5.3 A nonsingular finite connection minor

The `r` contexts `C_0,...,C_(r-1)` recover `F(1),...,F(r)` by triangular inversion. If a polynomial in `V_n` has this jet zero at `1/2`, parity supplies multiplicity `r` at `-1/2`; that would require degree at least `2r>n-1`. Thus those `r` probes are injective on all of `V_n`.

Together with the explicit basis family, this proves

    rank (H(C_m[B_j]))_(0<=j,m<r)=r,
    rank M_n=r=ceil(n/2).                               (9)

The affine hull of the realized profiles has dimension exactly `r-1`, because all are monic in `Q` and the same basis differences are independent.

**Minimality scope.** Equation (9) is an exact linear-realization lower bound. Equation (6) is optimal for linear reconstruction on the entire affine hull. It is NOT a lower bound on the number of integer probes needed merely to distinguish the finite set of realized tournament profiles, on arbitrary encodings, on bit complexity, or on time needed to compute any of these #P observables.

### 5.4 The same carrier is the signed linear-forest polynomial

Let `S=A-A^T`, and let `N_L(T)` be the sum of the products of the signed arcs over all directed linear forests with `L` arcs, including isolated vertices. Expanding the consecutive-arc product in each vertex ordering counts every such forest exactly `(n-L)!` times. Therefore

    E_T(u)=sum_(L=0)^(n-1) 2^(-L)(n-L)! N_L(T) u^(n-1-L). (9a)

Reversing an odd-length component is a sign-reversing involution, so `N_L=0` for odd `L`. This is exactly the forest/master-polynomial dictionary recorded in THM-4757, not a new forest expansion. It proves a useful new consequence of tomography: the measurements (6) recover every coefficient of the signed linear-forest polynomial, and hence its full formal logarithm. The path-cover, endpoint-jet, path-colour and full forest descriptions have the same information content. A truncated spectral one-loop approximation need not.

## 6. Quotient contextual equivalence: the second exact state

Let `Q` now be a labelled quotient tournament on `[q]` and define the ordinary commutative walk series

    W_Q(x)=1^T X(I-A_Q X)^(-1)1=sum_c w_Q(c)x^c.

Unlike the square-zero path series, its vertex variables may repeat: each occurrence represents a different run in a substituted module. Let

    W_Q^+(x)=sum_(c_i>=1 for every i) w_Q(c)x^c,           (10)

the **full-support part**, not the whole walk series.

THM-3121 says that for nonempty factor tournaments `S_i`,

    H(Q[S_1,...,S_q])
       =sum_(c_i>=1) w_Q(c) prod_i (c_i! p_(S_i)(c_i)). (11)

Only (10) is consumed. In particular, putting each `S_i=P_(m_i)` gives

    h_Q(m):=H(Q[P_(m_1),...,P_(m_q)])
       =sum_(1<=c_i<=m_i) w_Q(c) prod_i c_i!S(m_i,c_i). (12)

Let `s(a,b)` be signed Stirling numbers of the first kind, inverse to `S(a,b)`. The transform is lower triangular in every variable with nonzero diagonal, so

    w_Q(c)=1/(prod_i c_i!)
           sum_(1<=m_i<=c_i) h_Q(m) prod_i s(c_i,m_i).   (13)

**Theorem C.** For two labelled `q`-vertex quotients `Q,R`, the following are equivalent:

- all Hamiltonian counts `H(Q[S_1,...,S_q])` and `H(R[S_1,...,S_q])` agree for every tuple of nonempty tournaments;
- they agree for every tuple of transitive nonempty tournaments;
- `W_Q^+=W_R^+` coefficientwise.

Proof. Equation (11), restriction to transitive factors, and inversion (13), respectively. QED.

Thus `W_Q^+` is the exact quotient-side contextual state. The full `W_Q` is sufficient, while formula (11) directly consumes only its full-support part. No strict separation between equality of `W_Q` and equality of `W_Q^+` is asserted here. If deletion/empty slots are admitted as additional contexts, apply the same theorem to every induced quotient: the family of full-support series is then exactly the full `W_Q`. Do not omit this support distinction.

No finite-rank formula for this quotient-context matrix is proved here. Converse reverses every walk and preserves its content, so neither `W_Q` nor `W_Q^+` recovers orientation in general.

## 7. Hostile boundaries and what the objects fundamentally retain

- Equal `H` is insufficient even for module contexts: the reconstructed masks40/76 have `H=15` but responses123/135 to `C_1`.
- Equal full profile is sufficient for all one-hole module contexts, and therefore none of those contexts can distinguish the THM-3235 twins. Yet the twins as quotients have `H(Q[P_2,...,P_2])=421425/411513`, proving the interface category changed.
- The profile counts **how many disjoint runs** can be exposed to a uniform exterior. It forgets which endpoints those runs have, which vertices carry them, and how a vertex-resolved exterior can join them.
- The quotient series counts **how run visits can be threaded through labelled slots**. Its coefficients are walk multiplicities by vertex-content, not labelled path order or adjacency reconstruction.
- The square-zero port matrix restores vertex exclusion and endpoints. It is the right completion for arbitrary labelled gluing, at an explicit exponential cost.
- The fixed-path metagraph has still another observable: transport and compatibility of a chosen Hamiltonian path under its moves. The module rank theorem supplies no automatic metagraph wiring theorem.

The concrete direction is not a universal scalar invariant. It is a hierarchy of interface-relative response states, with actual reconstruction theorems separating sufficiency from minimality and explaining exactly when label information must return.

## 8. Verification

Run from the parent workspace:

    python3 session36-context/verify_context.py

No repository code or data are imported. The verifier separately implements endpoint-subset Hamiltonian DP, least-vertex set-partition path-cover DP, permutation-and-cut enumeration, and odd-cycle collection enumeration.

- All 1,099 labelled tournaments of orders1 through5: path covers by two independent algorithms; coloured cycle formula (7); parity/affine tomography; direct vertex-level context Hamiltonian checks (2,120 total probe instances).
- Deterministic random controls:128 at order6,64 at7,32 at8,16 at9. Each has full-profile reconstruction and a direct Hamiltonian DP for its largest probe, reaching15 vertices.
- The full-rank near-transitive family for every order1 through12: full profile directly computed, formula (8) checked, actual finite connection determinants nonzero.
- The same-H hostile masks40/76 rebuilt directly, with responses123/135.

Raw outputs: `verification.out`, `verification.json`. The proof is all-order; these tests audit implementations and conventions, not establish the quantifiers.

A second independent companion is

    python3 session36-context/verify_quotient.py

It constructs substituted tournaments directly and verifies the Stirling inverse against an unrelated content/last-letter walk DP on1,258 cells: every labelled quotient of order2 with block sizes1..3, order3 with sizes1..3, and order4 with sizes1..2. It also reconstructs the THM-3235 profile twins and their twelve-vertex Hamiltonian counts421425/411513. Outputs: `quotient_verification.out`, `quotient_verification.json`.
