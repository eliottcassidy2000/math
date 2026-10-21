# Distinguished partial paths: exact orbit channels and rational deweighting

Status: PROVED identities and FINITE-EXACT controls. This is a session-local
research packet, not a new scarce theorem-ID reservation or a literature
priority claim. The classical Burnside mechanism and the full-Hamiltonian
path/tiling specialization are inherited; the result here is their explicit
partial-path interpolation, fixed-deficit rational series, and typed symmetry
boundary.

## 1. Inheritance, objects, and lost information

Repository sources read:

- `01-canon/theorems/THM-260-rooted-tournament-layer-decomposition.md`: rooted
  tournament counts and the automorphism-deficit viewpoint. Its displayed
  general mixed-cycle denominator has the narrow error recorded below.
- `01-canon/theorems/THM-283-sc-burnside-formula.md`: the corrected
  anti-automorphism cycle filter, including endpoint-swap monodromy.
- `05-knowledge/results/converse_tournament_count_channels_20261021/`:
  `F_m(a,b)` and the half-order self-converse reductions.
- `01-canon/theorems/THM-644-gridsym-fiber-law-and-anti-redei.md`:
  the full-path tiling fiber and positional-reversal anti-action. We use its
  explicit definitions, not any historical conjectural status as a proof.
- `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`:
  full path covers versus partial marked paths, and the requirement to audit
  when a component action is free. This packet counts one distinguished path
  on a subset, not a spanning cover or a quotient-walk kernel.
- `01-canon/theorems/THM-024-involution-anti-aut.md`: involutory anti-actions
  exist, but anti-automorphisms need not all be involutions.

For `n >= k >= 1`, let `P(n,k)` be the number of isomorphism classes of pairs

    (T, (v_1,...,v_k)),     v_1 -> v_2 -> ... -> v_k,

where T is an n-vertex tournament and the v_i are distinct. An isomorphism
must preserve every path position. A one-vertex path is a distinguished root.
Let `B(n,k)` count the fixed classes of the involution

    (T,(v_1,...,v_k)) -> (T^op,(v_k,...,v_1)).

The converse-merged path count is `(P(n,k)+B(n,k))/2`. Reversing the path is
essential: merely taking the converse while fixing its order does not preserve
the directed-path condition for k >= 2.

For an odd partition lambda of m, write

    K(lambda) = number of parts,
    z_lambda = product_l l^(m_l) m_l!,
    E(lambda) = sum_i (lambda_i-1)/2 + sum_(i<j) gcd(lambda_i,lambda_j),
    F_m(a,b) = sum_(lambda odd, |lambda|=m) 2^(a E(lambda)+b K(lambda))/z_lambda.

The empty partition has E=K=0 and z=1. This is exactly Session22's channel,
not a new definition of a directed-cycle statistic. All generating functions
below are initially formal.

Connection contract: deleting the distinguished path leaves a tournament on
m vertices with a k-bit adjacency color at every vertex, plus the unfixed
internal arcs on the ordered path vertices. This preserves the complete
marked isomorphism class. Forgetting the path loses its embedding and changes
stabilizers; retaining only a count does not recover individual tournaments,
Hamilton-path multiplicities, or metagraph adjacency.

## 2. Ordinary path count and the labelled-scaffold extension

**Theorem 1.** Put m=n-k. Then

    P(k+m,k) = 2^binom(k-1,2) F_m(1,k).                    (1)

More generally, prescribe a fixed oriented simple graph D on the k ordered
marked positions, with e arcs and no oppositely oriented duplicate pair. Count
injective, not necessarily induced, copies of D whose individual positions
remain distinguished. The result is

    2^[binom(k,2)-e] F_m(1,k).                            (2)

Thus, for this rigidly labelled scaffold notion, its shape is irrelevant to
the total: only k and e matter. This is not the count of unlabelled copies
modulo automorphisms of D.

**Proof.** Fix the marked vertices as 1,...,k and quotient only the remaining
m vertices by S_m. There are `binom(k,2)-e` free internal arc bits. Each
unmarked vertex carries its k adjacency bits to the marked vertices, so the
remaining object is a tournament with 2^k vertex colors. Session22's Burnside
formula gives F_m(1,k). A path has e=k-1. This also gives a direct inverse to
the deletion map above. QED.

For k=n this recovers the usual tiling count `2^binom(n-1,2)`. For k=1 it
recovers the rooted count `F_(n-1)(1,1)` underlying THM-260. Ordered k distinct
roots without arc constraints instead give `2^binom(k,2) F_(n-k)(1,k)`.

## 3. The exact fixed-deficit rational series

Let

    L_m = {0},                                      m=0,
    L_m = {K: 1<=K<=m and K == m (mod 2)},            m>=1,
    w_(m,K) = sum_(lambda odd m, K(lambda)=K) 2^E(lambda)/z_lambda.

Every K in L_m occurs: for m>0 take K-1 singleton parts and one part
m-K+1. Therefore every w_(m,K) is strictly positive. By (1),

    A_m(k) := P(k+m,k)/2^binom(k-1,2)
            = sum_(K in L_m) w_(m,K) 2^(kK),         k>=1.

**Theorem 2.** The normalized ordinary generating function is

    G_m(z) := sum_(k>=1) A_m(k) z^k
            = sum_(K in L_m) w_(m,K) (2^K z)/(1-2^K z).      (3)

Its exact reduced denominator over Q is

    D_m(z) = product_(K in L_m) (1-2^K z).                   (4)

Consequently A_m satisfies the minimal constant-coefficient recurrence with
characteristic polynomial `product_(K in L_m)(X-2^K)`. For m>=1 its order is
ceil(m/2); for m=0 it is one. No pole cancels because each corresponding
partial-fraction weight is positive and all poles are distinct.

Examples, for k>=1:

    A_0(k) = 1,
    A_1(k) = 2^k,
    A_2(k) = 4^k,
    A_3(k) = (4*8^k + 2*2^k)/3,
    A_4(k) = (8*16^k + 4*4^k)/3,
    A_5(k) = (128*32^k + 40*8^k + 12*2^k)/15.

For instance `G_3(z)=4z(3-8z)/[(1-8z)(1-2z)]`. This is a rational series
obtained only after removing the Gaussian-sized free-arc factor in (1).

## 4. Converse-fixed partial paths

**Theorem 3.** Write

    k=2r+delta,    m=n-k=2s+epsilon,    delta,epsilon in {0,1}.

If delta=epsilon=1 then B(n,k)=0. In every other case,

    B(n,k) = 2^[r(r-1+delta+epsilon)] F_s(2,k+epsilon).      (5)

The even-path range is r>=1, while the odd-path range includes r=0.

**Proof.** With the marked positions fixed, the converse-twisted Burnside
sum runs over maps g that reverse positions 1,...,k and arbitrarily permute
the m unmarked vertices. Their restriction to the path is r transpositions
and delta fixed points. By the corrected THM-283 monodromy criterion, the
unmarked permutation must have type `2lambda` plus epsilon singleton parts,
where lambda is an odd partition of s, and delta+epsilon<=1. This explains
the vanishing case.

Let K=K(lambda), E=E(lambda). The pair-orbit counts are:

- internal to the marked vertices: r(r+delta);
- internal to the unmarked vertices: 2E+(1+epsilon)K;
- crossing the two sets: (2r+delta)K+epsilon*r.

The k-1 adjacent path pairs form exactly r pair-orbits under reversal.
Prescribing their forward orientations consistently pins those r bits.
The unmarked permutation centralizer has order `2^K z_lambda`. Thus its
contribution after cancellation of that 2^K is

    2^[r(r-1+delta+epsilon)] * 2^[2E+(k+epsilon)K]/z_lambda.

Summing gives (5). This proof permits noninvolutory anti-actions on the
unmarked vertices; restricting them to involutions would be wrong. QED.

At m=0 the formula recovers the full-path grid-reflection fixed count:
`B(2r,2r)=2^[r(r-1)]`, `B(2r+1,2r+1)=2^(r^2)`. At k=1 it gives zero for even
n and exactly C_n for odd n, with the self-converse count C_n of THM-283.

A generalization of the same proof: if an anti-invariant scaffold on k
ordered positions is preserved by their reversal and its prescribed pairs
occupy q pair-orbits, replace the prefactor exponent in (5) by
`r(r+delta+epsilon)-q`. The assignments must be consistent with the
anti-action. A directed path has q=r. This extension does not permit two
fixed vertices or contradictory prescribed arc orientations.

## 5. Rational converse parity series

Fix m=2s+epsilon and delta in {0,1} with delta*epsilon=0. Set r_0=1-delta and

    V_(m,delta)(r)
      = B(2r+delta+m,2r+delta)/2^[r(r-1+delta+epsilon)]
      = F_s(2,2r+delta+epsilon),                         r>=r_0.

Put `v_(s,K)=sum_(lambda odd s,K(lambda)=K) 2^(2E(lambda))/z_lambda`.
Then

    sum_(r>=r_0) V_(m,delta)(r) z^r
      = sum_(K in L_s) v_(s,K) 2^[(delta+epsilon)K]
          (4^K z)^r_0/(1-4^K z).                       (6)

Its exact reduced denominator is

    product_(K in L_s)(1-4^K z).                        (7)

Again every simple-pole coefficient is nonzero. The forbidden parity channel is identically
zero, whose reduced denominator is 1; it must not inherit (7). The recurrence
order is ceil(s/2) for s>=1 and one for s=0. Thus converse symmetry halves
the residual order before the odd-cycle spectrum is used.

## 6. Rooted converse: one central orbit, not one central vertex

An ordered tuple of at least two distinct roots, with converse acting
pointwise on the root slots, cannot be self-converse: an anti-isomorphism
fixing both endpoints would have to reverse their mutual arc. For one root,
its self-converse count is C_n when n is odd and zero when n is even.

There is a stronger objectwise statement.

**Theorem 4.** In every self-converse tournament of odd order, exactly one
Aut(T)-orbit of vertices can be fixed setwise by converse on rooted
isomorphism classes. Equivalently, all vertices fixed by some
anti-automorphism form one Aut(T)-orbit. For even order there is no such
orbit.

**Proof.** Let G=Aut(T) and choose an anti-automorphism a. The anti-coset aG
acts on vertices, and its average fixed-vertex count equals the number of
G-orbits preserved by a. Indeed, a nonpreserved orbit contributes zero; in
a preserved orbit O, for each v in O the number of g with ag(v)=v is
|G|/|O|, so its total contribution after averaging is one. Every anti-action
has exactly n mod 2 fixed vertices, by its allowed cycle type. Therefore the
number of preserved G-orbits is n mod 2. An orbit is preserved iff some
anti-action fixes one, hence each suitably transported, vertex in it. QED.

This does not select a canonical vertex. The directed 3-cycle is a hostile:
its unique central root orbit consists of all three vertices. Rooted ordinary
counts also do not come from differentiating the unlabelled OGF: at n=3 the
rooted count is 4, while n*A_n=6.

## 7. The free-action threshold and its first symmetry correction

Write h_k(T) for the number of directed ordered simple k-vertex paths in T.
The automorphism-weighted mass is

    W(n,k) := sum_[T] h_k(T)/|Aut(T)|
            = 2^[binom(n,2)-k+1]/(n-k)!.               (8)

This is rational in general, not necessarily an object count. The labelled
identity behind (8) chooses an ordered k-tuple, fixes k-1 arc orientations,
then divides the total labelled incidence by n!.

The unlabelled marked count is instead

    P(n,k) = sum_[T] (1/|Aut(T)|)
               sum_(g in Aut(T)) h_k(T[Fix(g)]).        (9)

A permutation fixes an ordered path only if it fixes all its vertices. Every
nonidentity tournament automorphism moves at least three vertices. Thus
Aut(T) acts freely on all directed k-paths whenever m=n-k<=2, and (8) equals
P in this range. For every m>=3 the nonidentity odd-partition contribution
is strictly positive, so P>W. The first failure is n=4,k=1:

    P(4,1)=12,       W(4,1)=32/3.

This is why the familiar full-Hamiltonian fiber identity cannot be blindly
extended to arbitrary partial paths.

For fixed m>=3, the first nonidentity partition is `(3,1^(m-3))`, with
K=m-2. Its weight is `2^[1+binom(m-2,2)]/[3(m-3)!]`. Therefore, as k tends
to infinity through positive integers,

    P(k+m,k)/W(k+m,k)
      = 1 + c_m*4^(-k) + O_m(16^(-k)),
    c_m = [m(m-1)(m-2)/3] * 2^(4-2m).                 (10)

For m=3,4 the remainder is identically zero. This is a rigorous statement
about a weighted symmetry correction, not an asserted exact probability of
a nontrivial automorphism. The full finite expansion is already given by
(3), so no asymptotic limit is required to compute it.

## 8. Formal versus analytic scope

The normalized G_m has convergence radius 2^(-m), including radius 1 at
m=0. Every nonzero normalized converse parity series in (6) has radius
4^(-s), including radius 1 at s=0.

In contrast, for each fixed m>=0 the raw generating series
`sum_(k>=1) P(k+m,k) z^k` has radius zero, even after division of coefficients
by k!. The Gaussian free-arc factor is at least `2^binom(k-1,2)`. Similarly,
the unnormalized nonzero converse parity series has radius zero, even after
division by r!. These raw series remain valid formal generating functions.

They are not D-finite over Q(z). A polynomial-coefficient recurrence of fixed
order forces the coefficients of any formal solution, after finitely many
initial terms, to be bounded by `C^n (n!)^d` for some C,d: solve for the last
term outside the finitely many roots of its leading polynomial and iterate
a polynomial bound on the ratios of coefficient polynomials. The positive
coefficients here grow at least as `2^(c n^2-O(n))`, contradicting that bound.
This argument is about the raw fixed-deficit sequences, not a claim that
zero radius alone rules out D-finiteness. The normalized sequences in
Sections 3 and 5 are rational and therefore D-finite.

## 9. Narrow repair of THM-260's mixed-cycle denominator

THM-260's general displayed layer divides by `product_i k_i * r! * m!`, where
k_1,...,k_r are the nontrivial odd cycle lengths and m is the fixed-point
count. The correct denominator is

    m! * product_(l>=3 odd) l^(m_l) m_l!.

The r! factor is valid only when all r nontrivial cycle lengths are equal.
For cycle type (5,3), n=8, E=4 and no fixed points, the correct rooted-deficit
layer is `2^4 * 8/(5*3)=128/15`, agreeing with THM-260's own table. Its
stated general denominator instead gives 64/15. The repaired full layer is

    D_lambda(n) = (n-m_1) * 2^E(lambda)/z_lambda.

This follows immediately by subtracting the pointed Burnside term
`m_1*2^E/z_lambda` from `n*2^E/z_lambda`. It changes no correct table entry
and requires no alteration to the inherited rooted counts. THM-260's separate
"graded ring" sentence is not used here: no layer multiplication has been
defined, and disjoint permutation cycles still interact through cross-gcd
terms and centralizer multiplicities.

## 10. Exact controls and next questions

Run:

    python3 check.py --receipt receipt.json
    python3 -O check.py --receipt optimized-receipt.json

The script imports none of the repository's producer implementations. Its
literal labelled-mask orbit census quotients the unmarked vertices directly
for all 21 cells `1<=k<=n<=6`. It checks ordinary, reversed-converse and merged
counts; directly verifies the objectwise central-root-orbit theorem on all
76 ordinary classes at n=1,...,6; checks the weighted path mass by a separate directed-path DFS in every class;
and checks the normalized recurrences and leading symmetry correction for
m=0,...,12. Normal and optimized runs both execute 8,954 explicit checks.
The receipts include every count and all denominator/prefix data. Optimization
does not disable tests: checks use explicit exceptions, not assert.

Cheap hostile controls include:
- n=3: rooted count 4 versus `n*A_n=6`;
- k=1: even-order self-converse rooted count zero;
- n=3,k=2: reversal-fixed marked-path count 2, versus zero for two roots
  held pointwise under converse;
- n=4,k=1: marked count 12 versus weighted mass 32/3;
- mixed (3,5) cycle multiplicity: 128/15 versus the erroneous 64/15.

Useful next directions, not claimed here:
1. Mark an unlabelled path or cycle with a nontrivial scaffold symmetry group.
   Then quotient by that group as well; the extra Burnside terms, absent from
   (2), are exactly the missing coordinate.
2. Refine the omitted-vertex tournament by a score/degree or strong-component
   sidecar. A common finite carrier might preserve rationality in k, but one
   must track how the k-bit colors constrain the target property.
3. Extract a stabilizer-size distribution rather than the weighted correction
   (10). Burnside totals alone do not identify how many marked objects have
   each nontrivial stabilizer.
