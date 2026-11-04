# Fixed-quotient module substitution: exact defect growth and a contextual no-return theorem

2026-11-04. **PROVED with a finite symbolic skeleton certificate:** an explicit growing strong family retaining the address24 old graph, its complete seven-coordinate signature, a sharp transitive all-markings joint-defect bound, and exclusion of both next requested digits for every replacement tournament module, even after arbitrary marked ordered contexts. **FINITE-EXACT:** independent permutation/subset-DP checks and literal full attachments. Reusable two-branch closure, unrestricted spectrum completeness, and native-complement legality remain **OPEN**. Historical novelty is not established.

## 1. Inheritance and precise scope

The source snapshot is `8880585ffd98072e93aa39d97476c7a1d4c9129d` in `eliottcassidy2000/math`. The targeted sources are:

- [S182 induction, seven-coordinate law and terminal-port criterion](https://github.com/eliottcassidy2000/math/blob/8880585ffd98072e93aa39d97476c7a1d4c9129d/05-knowledge/results/session182_physical_returns_20261104/induction/PROOF.md).
- [S185 closure, actual SCC7/8 escape witnesses](https://github.com/eliottcassidy2000/math/blob/8880585ffd98072e93aa39d97476c7a1d4c9129d/05-knowledge/results/session185_proof_frontiers_20261104/closure/PROOF.md).
- [S185 mining, ordered-context barrier and suggested mixed-module inflation](https://github.com/eliottcassidy2000/math/blob/8880585ffd98072e93aa39d97476c7a1d4c9129d/05-knowledge/results/session185_proof_frontiers_20261104/mining/PROOF.md).
- [S185 carriers, constant-bounded fixed-front tail inverse](https://github.com/eliottcassidy2000/math/blob/8880585ffd98072e93aa39d97476c7a1d4c9129d/05-knowledge/results/session185_proof_frontiers_20261104/carriers/PROOF.md).
- [S187 incoming, scalar observer loss and prepared physical returns](https://github.com/eliottcassidy2000/math/blob/8880585ffd98072e93aa39d97476c7a1d4c9129d/05-knowledge/results/session187_named_response_20261104/incoming/NOTE.md).

The closest proved mechanism is the ordered signature law. The hostile is the address24 carrier: arbitrary ordered context cannot continue it, but S185 supplies genuinely different replacement carriers realizing both children. The proposed underexplored sidecar was its mixed-port homogeneous pair. This packet first settles that transitive-inflation route and then extends the exclusion to every replacement tournament module in the same fixed quotient, allowing arbitrary nonempty internal port masks and arbitrary ordered context.

For a named marked suffix `(B,W)`, write

    sigma=(h,p,F,G,J,e,s),     c=h-J.

Here h is the Hamilton-path count; p is the ordered two-slot spanning path-cover count with empty slots allowed; F counts covers whose first slot is nonempty and ends in W; G counts covers whose second slot is nonempty and starts in W; J counts both conditions jointly; e and s are the marked Hamilton end/start counts. The actual graph and actual port subset remain part of the state.

The inherited strong attachment has TT_t -> B, W -> r, and r -> TT_t and B\W, for t>=1. Its count is

    H_t(B,W)=c+2^t F,

and its physical root defect is J. It is strong iff W meets the terminal SCC of B. A clock-uniform d return, with any integer clock shift q, necessarily has target constant `c'=3c+d`. Thus the address24 constant c=-35 requires -103 or -101 for d=2 or 4. Constant failure alone excludes every q; no bounded clock search is involved.

## 2. Actual mixed-module family and strongness

The S185 address24 suffix has outgoing rows

    (2,20,9,3,45,79,31),     W={1,3,5},

and sigma=(35,676,252,267,70,13,14). Its homogeneous pair is {3,0}, directed 3->0, with 3 marked and 0 unmarked.

Collapse that pair to a proxy v. The remaining vertices, in the code order, are `(1,2,4,5,6)`. Replace v by the transitive block

    x_0 -> x_1 -> ... -> x_(m-1),       m>=1,

with every external edge inherited from v=3. Keep the external ports {1,5}; let S be any subset of the m clones and mark precisely `{1,5} union S`. Denote the result B_m(S). The code `carrier.py` constructs the literal outgoing rows and ports.

The quotient contains the directed spanning cycle

    v -> 1 -> 4 -> 5 -> 6 -> 2 -> v.

It is therefore strong. Inflating one vertex of a strong quotient by a nonempty tournament preserves strongness, since routes through the other quotient vertices connect any two vertices in the same module, and quotient routes lift between modules. Thus every B_m(S) is strong, including S empty. The fixed external ports are nonempty, so every full attachment at t>=1 is strong.

For m>=2, identify x_0 with old vertex3 and x_(m-1) with old vertex0. In the first-clone family S={x_0}, increasing m inserts one fresh clone immediately before old vertex0. All old names, all old arcs and their port statuses survive. Every completed direct birth is strong. This is a literal adjacency construction, not a certified program of Hamilton-complement operations.

## 3. Exact seven-coordinate family

For m>=2 and S={x_0}, the full signature is

    h = 9*2^m - 1,
    p = 66*3^m + 28*2^m - 30,
    F = 66*2^m - 12,
    G = 24*3^m + 15*2^m - 9,
    J = 35*2^(m-1),
    e = 13,
    s = 7*2^(m-1).

Hence

    c = -17*2^(m-1) - 1,
    H_(t,m) = 2^(m-1)*(132*2^t-17) - 12*2^t - 1.

This is an explicit unbounded strong carrier with its root defect retained. At m=2 it is the original address24 representative up to the displayed name encoding. One direct clone birth gives

    H_(t,m+1) = 2 H_(t,m) + 12*2^t + 1.

The 3^m channel in p and G cancels out of the actual attachment observer; it is not a ternary return. The family does not satisfy either desired 3H+d map merely by increasing the clone clock.

### Finite symbolic proof

Fix an order of the five exterior vertices. A gap can hold a clone run exactly when its preceding exterior vertex, if any, beats v and v beats its following exterior vertex, if any. Every backward exterior adjacency must have a clone run inserted. If k gaps are available and r specified gaps must be nonempty, assigning m distinct linearly ordered clones contributes

    sum_(a=0)^r (-1)^a binom(r,a) (k-a)^m.

Each clone is assigned to a gap; its order inside that gap is forced. This is a bijection, not an interpolation from sample m. For two-slot covers, also choose the cut in the exterior permutation; all backward adjacencies within the two slots must be occupied. Endpoint restrictions remove the corresponding boundary gap. For a distinguished first clone ending a path, its forced gap cannot contain any other clone; when it starts a path, that gap may contain every greater clone.

`family_coeffs.py` enumerates exactly 5! exterior orders, and the six cuts where applicable. After inclusion-exclusion its power-basis coefficients are:

| Observable contribution | Exact expression |
|---|---|
| h | 9*2^m - 1 - 3*0^m |
| p | 66*3^m + 28*2^m - 30 |
| F, endpoint external | 33*2^m - 16 - 3*0^m |
| F, endpoint first clone | 66*2^(m-1) + 4 - 15*0^(m-1) |
| G, start external | 18*3^m + 14*2^m - 6 |
| G, start first clone | 18*3^(m-1) + 2*2^(m-1) - 3 |
| e, endpoint external | 4 - 4*0^m |
| e, endpoint first clone | 9 - 4*0^(m-1) |
| s, start external | 3*2^m - 0^m |
| s, start first clone | 2^(m-1) |

The convention is 0^0=1. For m>=2, combining the disjoint endpoint contributions gives the claimed coordinates. The next section supplies J. The zero-base terms explain why m=1 should not silently be substituted into the simplified family.

## 4. Exact named singleton defect and the sharp all-markings bound

With no clone marked, the joint defect is

    J_empty(m)=6*2^m-4.

Let clone x_i be the unique marked clone, where i+j=m-1 and i,j>=0. Then

    J_i(m)=6*2^m-4
           +(18*3^i+5)*2^j + 8*2^i
           -1_(j=0)-4*1_(i=0).                       (1)

This retains the position of the marked clone. Permuting a mark to a different clone changes the actual response.

For completeness, `skeleton.py` computes the complete symbolic coefficient table. The both-external contribution is `6*2^m-4`. For x_i at the first path's end and an external marked second-path start, the coefficient terms in the pair of bases for i smaller and j larger clones are

    +18*(3^i*2^j) + 8*(2^i*1^j) - 1*(1^i*0^j).

For an external marked first-path end and x_i at the second path's start, they are

    +5*(1^i*2^j) - 4*(0^i*1^j).

These three disjoint cases exhaust J with one marked clone: the same clone cannot belong to both disjoint paths. Their sum is (1). The same gap bijection establishes each coefficient for every i,j, rather than a bounded fit.

**Sharp bound.** For every m>=2 and nonempty clone-port set S,

    J(B_m(S)) >= 35*2^(m-1),
    c(B_m(S)) <= -17*2^(m-1)-1.                       (2)

Equality holds exactly for S={x_0}.

Proof. At i=0 and m>=2, (1) gives 35*2^(m-1). For i>=1, subtracting this value yields

    [18*3^i+5-23*2^i]*2^j + 8*2^i - 4 - 1_(j=0),

which is strictly positive: `(3/2)^i>=3/2` makes the bracket at least `4*2^i+5`, and the final terms are positive. Adding ports cannot decrease J. More precisely J is a quadratic sum of nonnegative ordered endpoint-pair counts. Each new clone contributes its strictly positive singleton increment over J_empty, so a second mark makes the inequality strict. Finally h=9*2^m-1 gives the constant bound. QED.

The m>=2 boundary matters. At m=1 the sole-clone signature has h=17,J=34,c=-17. That contraction deletes an old module name and is not included in the no-return theorem below.

## 5. No continuation even with arbitrary ordered contexts

**Theorem.** Let m>=2, let S be any nonempty clone-port set, and let X,Y be arbitrary marked tournaments on disjoint names, with either allowed to be absent. Retain the internal arcs and ports of B_m(S), and form

    D = X -> B_m(S) -> Y.

Then `c_D` is neither -103 nor -101. Consequently no such target attachment, with any legal fixed integer outer-clock shift, realizes either `3H+2` or `3H+4` from the actual address24 source.

This covers all transitive inflations of the original mixed pair retaining its marked old vertex, even allowing extra clone marks and arbitrary marked ordered contexts. It does not cover nontransitive module replacement, changes to external arcs or ports, another quotient, or an arbitrary non-ordered enlargement.

Proof. The inherited ordered-sum law, with its nonnegative joint corrections, gives

    c_(X->B->Y) <= h_X h_Y c_B.                      (3)

For absent contexts omit the corresponding factor. Since every tournament has h>=1 and c_B<0, adding such context cannot increase c_B.

By (2), all m>=4 have c_B<=-137, already below both targets. For m=2,3, exact named endpoint counting gives, in increasing nonzero binary module-port-mask order,

    m=2: -35, -59, -143;
    m=3: -69, -107, -273, -171, -337, -375, -611.

Thus only three candidates can reach either target by decreasing their constant:

| m | Module-port mask | c | F | G |
|---|---:|---:|---:|---:|
| 2 | 1 | -35 | 252 | 267 |
| 2 | 2 | -59 | 307 | 250 |
| 3 | 1 | -69 | 516 | 759 |

The receipt retains all seven coordinates, literal rows and literal port masks for every m=2,3 candidate. Each of the three rows satisfies, for both C=-103 and C=-101,

    3c < C < c,       c-min(F,G) < C.                (4)

If either X or Y is nontransitive then its Hamilton count is at least3, so (3) gives c_D<=3c<C. The fact h=1 iff transitive follows from the standard deterministic insertion of a directed triangle's three Hamilton paths; a nontransitive tournament contains a triangle.

It remains that both contexts are transitive. If X contains a port, its G coordinate is at least1 and the ordered joint formula subtracts at least F_B, giving c_D<=c_B-F_B<C. If X has no port but Y has one, the intermediate G remains G_B and the right-hand joint formula subtracts at least G_B, again below C. If neither context has a port, its joint/start/end coordinates vanish and c_D=c_B>C. These are all cases. QED.

The finite step is genuinely exhaustive after the all-order bound; it is not an extrapolation from failed small searches. No F matching or clock congruence needs to be assumed because the required constant is already impossible.

## 6. Verification and stopping boundary

Run from this directory:

    python verify.py
    python -O verify.py

The stored `normal.stdout` and `optimized.stdout` byte-match. `receipt.json` includes hashes of every source verifier and literal small cases. Verification uses explicit exceptions, not assertions removable by optimization.

The exact checks comprise:

- the complete finite symbolic coefficient tables for all-order defect and family formulas;
- 36 named singleton formulas at clone sizes1 through8;
- all 120 nonempty clone-port masks at sizes1 through6;
- complete independent permutation-versus-subset-DP endpoint response arrays at sizes1,2,3, including all directed two-slot joint coordinates;
- 22 literal full attachments at prefix clocks1 and2, with direct Hamilton counts and strongness;
- every small mask surviving the analytic order reduction, checked against the contextual inequalities.

The all-order conclusions follow from the gap bijection and inequalities, not from the finite sample range. No separate independent-agent audit is claimed inside this packet.

This closes the transitive-inflation proposal in S185 mining, including its combination with arbitrary ordered contexts. Section7 below also excludes nontransitive replacement of this pair. Neither result refutes S185's positive whole-suffix replacement witnesses. A productive next constructive move must change a quotient incidence, an external marking or another part of the carrier while controlling the full named response. Merely using p's nonzero 3^m channel does not supply an actual digit return. Global two-branch closure and spectrum coverage remain open.

## 7. Stronger extension: arbitrary nonempty tournament modules

The transitive assumption can be removed from the no-return theorem. This extension supersedes the proposed next step of merely choosing a nontransitive replacement for the same homogeneous pair.

**Theorem.** Replace the collapsed pair by any tournament A with m>=2 vertices, keep all its external arcs uniform as above, keep external ports {1,5}, and mark any nonempty subset U of A. Add arbitrary marked ordered contexts X and Y. The resulting suffix `X -> B[A,U] -> Y` has neither constant -103 nor -101. Consequently this entire fixed-quotient module-substitution architecture cannot supply either requested clock-uniform next digit from address24, for any fixed integer outer-clock shift.

No assumption of transitivity, strongness or a fixed port pattern is imposed on A. B[A,U] is strong by the same quotient-cycle argument. The statement retains the external quotient and external markings. It does not exclude changing a quotient arc, changing an external port, a non-uniform enlargement, or the already supplied S185 whole-suffix replacements. The restriction m>=2 retains both old names; m=1 contraction is outside this theorem.

### Decorated path-cover identity

Fix a marked vertex u in A. Let h=H(A), q2 be the number of ordered spanning covers of A by two nonempty directed paths, and E_k (respectively S_k) count ordered spanning covers by k nonempty directed paths whose **first slot** ends (respectively starts) at u. These are counts of actual directed paths, not scalar surrogates for a future named operation.

The same exterior skeleton contraction, now assigning arbitrary directed module paths to the occupied gaps, gives

    H(B[A]) = 17h + 9q2,
    J(B[A,empty]) = 8h + 6q2,
    J(B[A,{u}]) - J(B[A,empty])
       = 25E_1 + 44E_2 + 18E_3 + S_1 + 5S_2.         (5)

The endpoints in the last line are disjoint cases: u is the first external slot's end, or the second external slot's start. A singleton u cannot serve both endpoints because the two spanning path slots are disjoint.

`module_coeffs.py` directly enumerates the 5! exterior permutations and six cuts, retaining nonnegative occupied-gap counts before any power-basis expansion. Its coefficient table is:

| Exterior condition | k=0 | k=1 | k=2 | k=3 |
|---|---:|---:|---:|---:|
| Hamilton path | 5 | 17 | 9 | 0 |
| Both joint endpoints external | 2 | 8 | 6 | 0 |
| u is first-slot end | 0 | 25 | 44 | 18 |
| u is second-slot start | 0 | 1 | 5 | 0 |

The k=0 terms vanish because A is nonempty. In the last two rows the gap containing the distinguished path is forced; reorder the module cover so that this distinguished path is its first slot, and the remaining paths follow exterior-gap order. This is a bijection and validates (5) for every tournament A, not just transitive ones. It also independently explains the positive coefficients hidden in the signed power-basis formulas of Section4.

Subtracting the joint defect gives the exact one-port constant

    c_single = 9h + 3q2
       -25E_1 -44E_2 -18E_3 -S_1 -5S_2.             (6)

For larger U, nonnegative joint endpoint counts imply `c(B[A,U])<=c_single` for every selected u in U.

### Two injections paying the positive terms

Every Hamilton path of A either ends at u, contributing to E_1, or splits immediately after u into a first path ending at u and a nonempty remainder, contributing injectively to E_2. Concatenation recovers the original path whenever the joining edge is present. Therefore

    h <= E_1 + E_2.                                  (7)

For an ordered two-path cover of A, identify the unique slot containing u and retain its original slot number as one bit. If u is terminal in that slot, place that path first and the other path second: this is an E_2 cover plus the bit. Otherwise split the u-containing path immediately after u and encode `(prefix ending at u, nonempty suffix, other original path)`: this is an E_3 cover plus the bit. The two cases have different numbers of slots, and the bit and concatenation recover the source. Hence

    q2 <= 2(E_2+E_3).                                (8)

Substitute (7)-(8) into (6):

    c_single <= -16E_1 -29E_2 -12E_3 -S_1 -5S_2.     (9)

For m>=2, E_2>=1 and S_2>=1 by taking singleton u as the distinguished path and any Hamilton path on A-u as the other path. Also E_3>=q2(A-u), by taking singleton u followed by an ordered two-path cover of A-u. Every nonempty proper subset of A-u and its complement admit Hamilton paths, so

    E_3 >= 2^(m-1)-2.

Consequently, for every m>=2,

    c(B[A,U]) <= -34 -12(2^(m-1)-2)
                = -10 -12*2^(m-1).                  (10)

This bound is deliberately coarse. At m>=4 it gives c<=-106, already strictly below both -103 and -101. Ordered contexts can only lower a negative constant by (3), so all these module orders are eliminated at once.

At m=2 there is only the transitive module. At m=3 the two isomorphism types are TT3 and C3. The transitive possibilities were exhausted in Section5. For C3 and any singleton u,

    h=3, q2=6, E_1=S_1=1, E_2=S_2=2, E_3=2,

and (6) gives c=-115. Every larger nonempty port set only decreases c. Thus the only remaining candidates are exactly the three transitive rows already eliminated by the arbitrary-context argument in Section5. This proves the stronger theorem.

### Additional checks and resulting frontier

`module_extension.py` independently forms the literal quotient substitution for every labelled tournament A of orders1 through4 and every nonempty U. It computes h,q2,E_k,S_k by subset partitions inside A, checks the decorated identity against direct endpoint response on the full substituted graph, verifies both injections numerically, and checks the full small-order contextual reduction. It also verifies the finite nonnegative skeleton table and C3 control. Normal and optimized executions agree.

An unpublished exploratory search through module order7 found no target constant. That bounded computation motivated the proof; it is not needed for the theorem, and the proof is not an extrapolation from its empty search. Its generation/canonicalization is not independently audited, so no census-minimality claim rests on it.

The remaining constructive move must change something beyond the tournament substituted into this one homogeneous pair: an external marking or quotient incidence, another part of the carrier, or a genuinely non-uniform architecture. This is still a fixed-quotient no-return theorem, not a global spectrum obstruction or a proof that other two-branch mechanisms cannot exist.
