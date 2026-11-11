# Response-order budgets, ordinal congruences, and a head-balanced source family

2026-11-11. **PROVED:** a sharp run-order filtration, its exact accumulation under triangle wrapping, a closed bounded-observer congruence for ordinal reuse, and a literal all-k family preserving the core lane's zero-tail head and positive tail. **FINITE-EXACT:** new controls in `response_filtration_receipt.json`. **OPEN:** universal target coverage and arbitrary rooted-inflation rearming. This is a focused extension using the other lanes' results as inputs, not an audit of them.

## 1. A sharper rule for how much response a context consumes

Let F[M] replace one vertex of a tournament F by a nonempty unmarked tournament M, uniformly oriented to each of the k exterior vertices. Let C_p count ordered spanning covers by p nonempty directed paths.

**Sharp run-order theorem.** There are nonnegative integers a_(p,r)(F), determined by a finite compressed-skeleton enumeration, such that

    C_p(F[M]) = sum_(r=1)^(p+floor(k/2)) a_(p,r)(F) C_r(M).

The same upper bound applies with exterior endpoint restrictions, or with empty global path slots permitted. It is a bound on the module run order, not on the module's vertex count.

**Proof.** Compress each maximal M-run on each global path. If r runs occur on s paths meeting M, then there are r−s gaps between successive runs within those paths. Every such gap uses at least TWO exterior vertices. A single exterior vertex v between consecutive runs would demand both M→v and v→M, contradicting uniform module incidence. These gaps are disjoint, so 2(r−s)≤k. Since s≤p, r≤p+floor(k/2). Ordered local path covers expand a skeleton in C_r(M) ways. Summing the finite skeletons gives the formula.

This improves the generic k+p run bound without assuming transitivity of M. It uses the fact that M is ONE genuine module; applying it to interleaved nonmodules would be invalid. For several substituted modules, another large module counts as many exterior vertices unless a separate multiport skeleton bound is proved.

The bound is universally sharp. For k=2d, take exterior pairs a_i,b_i with M→a_i→b_i→M. A path can contain d+1 M-runs separated by the pairs; p−1 further paths each contain one M-run, for r=p+d. Choose the remaining exterior arcs to support the specified order. For odd k, add the last exterior vertex at the end of one path. Taking M large enough gives a real expanding cover and a positive leading coefficient. The theorem asserts sharp response order, not algebraic independence of all cover statistics on every restricted physical family.

## 2. Two added vertices can consume one response order forever

Let F[M] be the directed triangle blowup

    M → a → b → M.

The coefficient of C_(p+1)(M) in C_p(F[M]) is exactly p. At maximal run order, a,b must form the unique internal separator on one of the p ordered paths; all other paths contain one M-run. There are precisely p choices for that path.

Consequently repeated triangle wrapping gives

    C_p(F^ell[M]) = p(p+1)...(p+ell−1) C_(p+ell)(M)
                    + a nonnegative combination of lower C_j(M).

This holds for every p≥1 and ell≥0, with the empty product equal to1. It follows by composing the linear response kernels and comparing the unique highest-order contribution. In particular the ordinary Hamilton observer after ell wraps has leading coefficient ell! on C_(ell+1)(M).

Thus an operation may be very small physically while spending an unbounded response-order budget under reuse. A fixed-order skeleton truncation does not remain closed under these repeated attachments. This does **not** rule out a different finite sufficient carrier on a special class, or prove that distinct actual tournaments with every prescribed truncated profile exist.

The new finite certificate gives the first three exact kernels:

    C_1(F[M]) = 3C_1 + C_2;
    C_2(F[M]) = 6C_1 + 8C_2 + 2C_3;
    C_3(F[M]) = 6C_1 + 15C_2 + 15C_3 + 3C_4.

These are certified by the proved finite cutoff and the triangular surjection basis, with extra transitive and nontransitive module controls. They are new tests of this filtration, not a rerun of the inherited Q-source certificate.

For nested single-module contexts with exterior counts k_1,...,k_ell, an initial p-path observer needs at most

    p + sum_i floor(k_i/2)

cover orders at the innermost module. Triangle nesting attains this budget exactly.

## 3. A closed exception: ordinal contexts spend no order

Define Q_r(M) to count ordered spanning r-path covers with empty slots allowed. For nonempty M,

    Q_r(M)=sum_(j=1)^r binom(r,j) C_j(M),
    Q_1(M)=H(M).

The lower-triangular change of coordinates has diagonal1. Hence retaining C_1,...,C_d is equivalent to retaining Q_1,...,Q_d.

For the literal ordinal operation A⇒B,

    Q_r(A⇒B)=Q_r(A)Q_r(B)                    for every r≥1.

Indeed each slot independently concatenates its A path and B path, ignoring empty portions. All A→B arcs exist, and every global slot has a unique restriction to each block. This is a bijection even when one or both restrictions are empty.

**Bounded-observer congruence.** For d≥1, equality of the jet (Q_1,...,Q_d) is preserved by ordinal composition on either side, converse, and repeated ordinal copying. If two pairs of modules have equal jets componentwise, ordinal composition of the pairs again has equal jets. Therefore equal d-jets can be exchanged under arbitrarily long expressions in these operations while preserving every external observer whose module skeleton uses at most d runs.

There is also a useful projective version. Put

    rho_r(M)=Q_r(M)/H(M),       1≤r≤d.

Nonempty tournaments have H>0. Then rho_r(A⇒B)=rho_r(A)rho_r(B), and converse preserves every rho_r. Equality of the normalized d-jet is consequently a congruence for these same operations. Equivalently it is equality of C_j/H through j=d. A homogeneous context response can scale by the appropriate H ratio, so projective equivalence alone must NOT be described as preserving its absolute output.

A finite proof certificate consists of the triangular coordinate conversion and the d multiplicative identities, plus actual initial modules. This gives a closed finite observer interface without the potentially large ideal computation in `PROOF.md`. The polynomial-ideal method remains available when the designated operations are not diagonal in these coordinates.

Internally marked modules need the corresponding endpoint tensor as well; no unmarked Q-jet is asserted to preserve arbitrary interior marks or an arbitrarily chosen cloning root. Here the core lane's mark and the constructor lane's central mark are exterior to every traded atom.

## 4. The core exchange guard becomes an inhabited congruence

The core lane has proved that swapping the modules in its quotient Q preserves BOTH the bare-apex head and every positive-tail output exactly iff

    C_2(A)/H(A)=C_2(B)/H(B),
    C_3(A)/H(A)=C_3(B)/H(B).

This is precisely equality of their normalized third jets. Thus the head-balanced exchange guard is closed under synchronized ordinal extensions:

    (A,B) → (A⇒C, B⇒D)

whenever C,D also have equal normalized third jets. The same is true of synchronized left extensions, converse, and equal positive ordinal powers. These are physical operations, not labels attached to a scalar count. One can change each module presentation repeatedly without spending additional response order.

By contrast, matching only the positive-tail projective coordinate [201H−62C_3 : 201C_2+120C_3] is not being declared a congruence for ordinal operations. The complete third jet is a proved sufficient closed interface; no unproved closure of the smaller observer is used.

Take the actual inherited atoms

    A=(2,4,9,19,39,15),
    B=(2,4,57,3,11,27).

Their full profiles agree, with initial entries

    H=17, C_2=190, C_3=828;
    Q_1=17, Q_2=224, Q_3=1449.

These are the same inhabited atoms used by the constructor lane, not new abstract source points. Let A_k and B_k be independently chosen ordinal words of k atoms, each position chosen from {A,B}, with k≥1. Then

    h=H(A_k)=H(B_k)=17^k,
    q=Q_2(A_k)=Q_2(B_k)=224^k,
    r=Q_3(A_k)=Q_3(B_k)=1449^k,
    C_2=q−2h,       C_3=r−3q+3h.

Every one of the 2k labelled atom positions can be traded independently. The third-jet congruence preserves the Q-source head and positive-tail response. In fact the inherited FULL profile equality and the all-r product law preserve every finite external module observer at each such atom, provided no operation chooses an internal root without retaining its extra response.

The k=0 empty ordinal word is NOT allowed as a substituted Q-module. The actual initial case here is k=1. This differs from the constructor sandwich, whose central marked singleton allows its separate legitimate k=0 initial case.

## 5. Exact source and head formulas for the new all-k presentation

Substitute h,q,r above into the core lane's proved kernels. Straight algebra gives

    K_k = 9h²−4hq+37q²,
    U_k = −8h²+19hq−4hr+10q²+62qr,
    E_k = 7hq.

The negative coefficients arise solely from the change from nonempty C-slots to Q-slots; the original compressed-run kernels are nonnegative and provide literal graph semantics. For all k≥1,

    bare-apex count = K_k+U_k−E_k,
    positive-tail count = K_k+2^t U_k,          t≥1.

Every source has 12k+6 vertices; bare-apex outputs have 12k+7; positive-tail outputs have 12k+7+t. Its 2k atom positions remain actual independent rebase locations, and swapping the two complete module words preserves all displayed responses. These are presentations of a new uniform family within the inherited universal Q constructor, not a claim of new scalar spectrum coverage beyond that constructor.

The genuine k=1 initial source is

    K_1=1,843,881,   U_1=20,596,980,   E_1=26,656;
    bare count=22,414,205;
    first positive-tail count=43,037,841.

The constructor lane uses the SAME atom words in its single-mark sandwich and obtains

    D=289^k,       U=V=3808^k.

Thus one atom library supplies two different source services with exact initial cases: a balanced selected sandwich and a head-balanced asymmetric Q exchange. Their scalar outputs differ, but their reusable local rebase guard is the same complete cover-profile congruence. This is a concrete connection rather than identifying the two constructors' offsets.

The inherited domination-defect inequality still charges the actual initial source: at bare birth its increment is U_k−E_k, and d(G_z)−1≤U_k−E_k. The congruence preserves that exact increment under all certified atom trades. It supplies no independent choice of domination count or offset.

## 6. Where a rooted operation spends the saved information

Use the inherited distinguishing clone ports: vertex1 of A and vertex0 of B, each cloned transitively to length2. New direct cover computations give

    cloned A: (H,C_2,C_3)=(35,606,3726),
    cloned B: (H,C_2,C_3)=(33,546,3426).

Their normalized-second and normalized-third cross differences are

    606·33−546·35=888,
    3726·33−3426·35=3048.

So the full head-preserving exchange guard already fails after this first rooted continuation. The failure is not merely different Hamilton counts: the PROJECTIVE response jets diverge. Replacing one atom by its profile twin may choose a different useful next operation, but it does not certify return after that operation.

This is exactly where the polynomial future-guard compiler must retain rooted response coordinates or add another control state. The ordinal congruence tells us which operations are free with respect to observer order; the sharp run filtration tells us which general attachments consume order; this explicit clone pair shows a different loss, namely root incidence hidden inside the module.

## 7. Reproduction and precise boundary

Run from the repository root:

    python3 -B 05-knowledge/results/session262_interface_supply_20261111/grammar/response_filtration.py

The receipt records sixteen new ordinal Q-jet composition controls, three exact triangle-wrapper kernels and their extra module controls, four all-k formula evaluations, and the two rooted-clone guard failures. The all-size proofs are the run-separator argument, slotwise ordinal bijection, and inherited literal Q kernels; the finite samples are not used to infer unbounded coverage.

**Outcome:** a sharp response budget plus an explicitly inhabited operation-closed equivalence. It explains why indefinite ordinal reuse is cheap, why repeated nonordinal triangle wrapping consumes arbitrarily high cover order, and why rooted cloning needs more than unrooted profile equality. It supplies the actual initial source and all-k rebasing family for the core exchange. Universal exact target coverage remains OPEN.
