# Full named modular cubics have a strict support cutoff at even order

**PROVED:** an all-even-order support reduction for every named mixed cubic of the modular converse response. **FINITE-EXACT:** the radius-four neighborhoods of two specified order-twelve controls and a literal prime twelve-vertex full integer cubic. **OPEN:** the global prime twelve-vertex hidden-converse residual and universal named response reconstruction.

## 1. The new support reduction retains the full named cubic

Use the physical response convention of incoming T1–T4 in `two_port_twelve_principles_20261103/PROOF.md`: x→y; a_u=1 means u→x and c_u=1 means u→y. Put δ=(F_T−F_(T^op))/2 modulo2. The division precedes reduction. In the squarefree algebra over F2, let G_uv enumerate simple directed paths from u to v including singleton paths when u=v, let s_u=sum_v G_uv, U=product_v(1+z_v), p=sum_u a_u s_u, q=sum_u c_u s_u and X=sum_vw a_v c_w G_vw.

Incoming T2 proves that the complete cubic part is

    [δ]_3 = L_odd((U+1)(p+q)X).

This is an inherited identity, not a new derivation of the observer. The following simplification is the increment here.

**Support-cutoff theorem.** If n is even, then

    [δ]_3 = sum over W, |W| odd and |W|≤n−3, of [z_W](p+q)X.    (1)

In particular every named modular cubic at order twelve is computed from two disjoint paths whose union has size3,5,7 or9. No sum is over unnamed isomorphism types; the full original labels remain.

**Proof.** Fix a nonzero squarefree monomial z_W of (p+q)X and put r=n−|W|. Its weight after multiplication by U+1 and application of L_odd is the parity of the number of subsets R of the unused names for which |W|+|R| is odd, plus 1 if |W| is odd. For r≥2 the first number is 2^(r−1), hence even, leaving exactly 1 for odd |W|. For r=1, even n makes |W| odd, and the weight is1+1=0. For r=0, the two identical contributions cancel (both are zero at even n). Thus only odd |W|≤n−3 survive. The two path supports in pX or qX are nonempty and disjoint, so |W|≥2 and the first possible odd size is3. QED.

For u≠v the specific named coefficient κ_(uv;w)=[a_u a_v c_w]δ is consequently

    κ_(uv;w) = sum over disjoint nonempty A,B,
                  |A|+|B| odd and ≤n−3,
                  (S_u(T[A]) G_vw(T[B])
                    + S_v(T[A]) G_uw(T[B])) modulo2.          (2)

Here G_vw(T[B]) counts Hamilton paths on B with exactly those endpoints, and S_u(T[A]) counts Hamilton paths starting at u; absent endpoints contribute zero. The same coefficient occurs at a_w c_u c_v. This equality follows either directly from incoming T2 or from X(a,c)+X(c,a)=pq and (p+q)pq=0 in the squarefree algebra. A repeated name w=u or w=v is allowed. In that case the only compatible middle path is the relevant singleton, so (2) reduces to T2's inherited shared-label even-support formula. Excluding singleton paths would corrupt the full cubic.

At order four the sum is empty, so every exact integer converse cubic coefficient is divisible by four. This is not a claim that those integer coefficients vanish. At order twelve all792 AAC coordinates and their792 mirrored ACC coordinates remain; there is no coefficient sorting, label renaming, or endpoint-only quotient.

## 2. An exact finite algorithm, with the integer observer kept separately

Compute every G_uv(T[B]) modulo2 by the recurrence

    G_uv(B) = sum over t in B−u with u→t of G_tv(B−u),
    G_uu({u})=1.

The implementation stores each fixed-start row as a bitset of possible endpoint parities. Its row parity supplies S_u(B). Enumerate ordered disjoint supports A,B and apply (2), updating the original named coordinates. This takes O(n³·3^n) elementary work in the simple implementation and O(n·2^n+n³) bitset/integer storage. It never substitutes z=1 into a squarefree product and never treats the two paths as independent endpoint marginals without support-disjointness.

The separate integer calculation evaluates every attachment assignment with port weights at most two and total weight at most three. It retains all1885 corresponding integer coefficients, including132 pure-port quadratics checked to vanish, using the inherited exact converse degree-three theorem. Twenty-four additional physical assignments check the resulting polynomial by literal Hamilton-path subset DP. This validates the whole named cubic, not only (2) or a shared-label projection.

Equation (1) is specific to the modular global converse cubic. It does not discard degree-four terms for an ambient internal module reversal. Incoming T1 makes δ=0 a necessary condition for a hidden ambient module converse. A nonzero coordinate of (1) therefore excludes that core in every uniform ambient quotient, but vanishing is not sufficient for equality of integer responses.

## 3. Precisely bounded twelve-vertex neighborhoods

The two source rows are

    A=(1284,1297,1322,1411,1325,1355,1311,1399,1536,1279,2048,1023),
    B=(3706,680,3603,1508,3914,1300,294,1397,519,1256,2370,1002).

A is the earlier iterated cyclic-extension endpoint-balanced control; B is the incoming prime endpoint-balanced control. Enumerate every set of zero through four distinct edges among their66 named old edges and reverse exactly that set. There are768692 presentations per source, counted separately even if the balls overlap. Endpoint parity is computed with exact F2 subset DP. For every endpoint-balanced target, incoming T3 is applied: each entry of (C_uv−C_vu)/2 modulo2 equals the parity of all odd-support paths between u,v in either direction. Every proper subset is checked for modularity when reporting whether an endpoint-balanced target is prime.

The complete receipt is `search.out`. For A the successive cumulative (tested, endpoint-balanced, prime endpoint-balanced) counts at radii0,1,2,3,4 are (1,1,0), (67,7,6), (2212,29,11), (47972,131,90), (768692,1021,703). For B they are (1,1,1), (67,1,1), (2212,3,3), (47972,32,32), (768692,386,386). Every row has zero combined-sieve survivors. Thus the two radius-four balls contain1537384 source/edit presentations and1407 endpoint-balanced presentations, of which1089 are prime; none passes C-skew modulo four. This is a finite neighborhood exclusion by necessary intrinsic obstructions. It does not enumerate every twelve-vertex tournament or settle the remaining prime universe. No filter conflates pair-free with prime.

## 4. A prime endpoint-balanced hostile and its full integer cubic

Reverse just the edge{1,10} of A. The resulting physical rows are

    P=(1284,273,1322,1411,1325,1355,1311,1399,1536,1279,2050,1023).

Every proper vertex subset is checked: P is prime. Its score multiset (2,2,3,3,5,5,6,6,7,8,9,10) differs from its converse score multiset (1,2,3,4,5,5,6,6,8,8,9,9), certifying nonselfconverseness. H(P)=84599. Exact endpoint vectors are

    E=(14994,12527,8034,6540,7043,6183,5902,4542,7191,4105,3769,3769),
    S=(4250,5257,5208,6424,6541,7725,10832,12928,3769,7447,4425,9793).

They agree entrywise modulo two. Forty unordered cover-skew entries fail modulo four, including C_02−C_20=26018. The entire named integer converse polynomial has1738 nonzero coordinates among1885 retained coefficients. All1584 mixed cubic coordinates are retained, of which760 have nonzero half-coefficient modulo two. The direct integer interpolation and the separate truncated two-path algorithm agree on every one of these coordinates.

Two explicit exact separators are

    [a_0 c_0 c_1](F_P−F_(P^op)) = −56686,
    [a_0 a_1 c_2](F_P−F_(P^op)) = −20574.

Both are2 modulo four, so incoming T1 excludes P as a hidden module-converse core in **every** uniform ambient quotient. This consequence is unbounded in the quotient, but applies to this particular core; it is not a global prime-twelve exclusion. The first separator deliberately repeats the physical name0 in different port variables, while the second uses three distinct old names. Independent literal finite differences also give the expected zero quartic at masks(3,3). Twenty-four further assignments in C++ and two in Python check the complete integer cubic away from its interpolation battery.

## 5. Reproduction

    g++ -O3 -std=c++17 -Wall -Wextra -Werror search.cpp -o search
    ./search 4 > search.out
    g++ -O3 -std=c++17 -Wall -Wextra -Werror cubic.cpp -o cubic
    ./cubic > cubic.json
    python verify.py > verification.json
    python -O verify.py > verification.optimized.json

The full integer cubic and all named zeros are in `cubic.json`. The small-control driver `run_controls.py` checks all2 labelled order-two tournaments, all64 labelled order-four tournaments and16 specified order-six presentations; every mixed cubic is compared against literal integer interpolation. All order-four modular cubics vanish as predicted. `verify_neighbors.py` independently checks all134 source/edit presentations through radius one with integer endpoint DP and actual two-cover convolution, reproducing counts(67,7,6) and(67,1,1) and no survivors. Python's separate integer endpoint and two-cover convolution checker verifies the principal literal witness and selected coefficients; the C++ path-bitset calculation checks every modular mixed cubic against integer interpolation.
