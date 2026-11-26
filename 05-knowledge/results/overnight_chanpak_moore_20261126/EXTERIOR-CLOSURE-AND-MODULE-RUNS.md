# Exterior closure at 128 links, and an exact run formula for spine reversals

2026-11-26. Continuation of the [overnight tandem](README.md) and [four-contact catalogue](LOW-LINK-TRADES-AND-ADJACENT-SPINES.md). All degree-57 claims are conditional separation statements about distinct full completions; they do not construct one. The spine statements are all-order identities within F1.

## 1. Exterior closure is not restricted to low branches

Let k=108+t be the number of changed links between two full charts, with t≤20. The old bound of 133 for four or more untouched branches means there are exactly the three frozen stars, and all 54 exterior branches are touched. Call a branch low if it has exactly four changed incident links and high otherwise, with h≤2t high branches. Let S be the set of **all** exterior branches with at most three high changed neighbors, regardless of their own degree.

Every active sheet changes at least four exits. In a branch of S, one of those exits must go to a low branch; its removed neighbor changes exactly four exits. The [difference identity](../moore57_quadratic_tensors_20261125/PRIMITIVE-HOLONOMY-AND-EDGE-SEPARATION.md) then bounds the active sheet's old inactive degree by four. Its head and exits to the three frozen branches exhaust that allowance. Thus every old exterior neighbor of an active sheet in S is active.

The total incidence into high vertices of the changed-link graph is 4h+2t. Each exterior branch outside S uses at least four such incidences. Therefore s=|S| satisfies

    s ≥ 54−h−floor(t/2) ≥ 54−2t−floor(t/2).

If t≤18, then s≥9. The old matchings among the s branches preserve activity, so they share an active size p>0. If p<56, a centered vector constant with opposite signs on the active and inactive parts within S and zero elsewhere has Rayleigh quotient s−1>7, impossible. Thus p=56, and one fully active S branch forces every exterior sheet active by the matching property and closure. This gives **127 changed links** in general if even one exterior sheet stays fixed. On the actual fixture a low branch has at most 28 active sheets, so this already excludes k≤126 without requiring an extra fixed sheet.

## 2. The actual fixture excludes 127

Take t=19. There are at least 16 low branches and s≥7. The s≥9 case is impossible as above. If s=8, let f be 1 on active sites in S, −p/(56−p) on inactive sites there, and zero elsewhere. This is a centered 7-eigenvector of the exterior adjacency: its Rayleigh quotient is 7, the maximum on centered vectors. But any inactive site in a low branch outside S has negative image, since all its neighbors in S are inactive. Such a site exists by the catalogue cap, and f is zero there. Thus s=7, h=38 and every high degree is five.

The independently replayed histogram of trades per quartet is 4 quartets of multiplicity 1, 4 of multiplicity 6, 141 of multiplicity 7, and 3 of multiplicity 14. A low branch has activity at most 14 except when its changed-neighbor quartet is one of the three exceptional ones. Each exceptional quartet contains two of the three hub branches 1,2,54. The hub branches have together at most fifteen changed incidences when t=19,h=38, so at most seven low branches can be exceptional. Outside S there are at least nine low branches. Consequently two of them have cap 14, and p≤14 by closure and bijectivity. Also p≥s−1=6: from an active sheet in S, its direct route and the s−2 two-step routes through the other S branches to any other S branch yield distinct active endpoints by girth five.

For the above f, ||f||²=56sp/(56−p), and Lf=(s−1)f on S, where L is the old exterior adjacency. The centered spectrum is in [−8,7], so λ²≤56−λ gives

    ||(L−(s−1)I)f||² ≤ (8−s)(s+7)||f||².

For each outside branch write g_y for the number of old neighbors of y in the active portion of S, m=sp. Then Σg_y=m and Lf(y)=(56g_y−m)/(56−p). In an outside exterior branch only active sheets can have g_y>0. If at most a sheets can be active, then Σg_y²≥Q(m,a), where

    Q(m,a)=a floor(m/a)² + (m mod a)(2 floor(m/a)+1).

Indeed, the minimum over at most a nonzero integral coordinates is obtained by spreading visits over a slots as evenly as possible. In the 47 outside exterior branches use two caps 14, seven caps 28, and 38 caps 56. The necessary inequality is

    56(2 Q(7p,14)+7 Q(7p,28)+38 Q(7p,56))−47(7p)² ≤ 14·7p·(56−p).

For p=6,…,14 the differences left minus right are 18620, 13083, 3136, 25627, 43708, 58947, 69776, 83251, 92316. All are positive. Thus **any two distinct full completions sharing the actual frozen stars differ in at least 128 whole branch matchings**, even before fixing the fourteen neighborhoods. The bound on individual moved edges in [the earlier checkpoint](OPPOSITE-ARMS-AND-BOUNDARY-VARIANCE.md) is a separate count.

## 3. An antipodal constraint and ten necessary 128-link lanes

Every active sheet of a branch in S lies on a changed link to a low branch. From the low-branch endpoint, that link exchanges its old target images by the half-turn τ. Because the link remains a permutation, the set of old images of changed sheets is τ-invariant. The activity set of a branch in S is a union of these sets from links to low branches. Thus its active set is a union of antipodal pairs, even when the branch itself is high. In particular p is even.

At t=20 there are l=54−h low branches. The sum of the three hub degrees is at most 15+(2t−h): start from at most five per hub and add the total excess beyond degree five at all high branches. Therefore at most E=floor((15+2t−h)/2) low branches are exceptional. Among the outside exterior branches there are at least L₀=l−s low ones, and at least B=max(0,L₀−E) of those have cap 14, so give L₀−B more branches cap 28.

We can also use the full old adjacency instead of the exterior compression: the same f, extended by zero, has unchanged first moment, and on the centered space the eigenvalues satisfy λ²=56−λ. The frozen branches have the same neighbor-count sum m per branch as the other outside branches and can conservatively be given cap 56. The root and heads contribute zero since the sum on each branch vanishes. This replaces 54−s by 57−s outside branches in the integer inequality with three extra Q(m,56) terms. The s=8 equality argument still excludes that case, since there is a low branch outside S. The [replay](exterior_closure.js) checks all remaining h,s,even p with the bounds above. At k=128 the necessary lanes are:

| high h | closed s | possible p |
|---:|---:|:---|
| 37 | 7 | 8 |
| 38 | 6 | 8,10 |
| 38 | 7 | 8 |
| 39 | 5 | 6,8,10,12,14 |
| 39 | 6 | 8,10 |
| 39 | 7 | 8 |
| 40 | 4 | 4,6,8,10,12,14 |
| 40 | 5 | 6,8,10,12,14 |
| 40 | 6 | 6,8,10,12 |
| 40 | 7 | 8 |

In the tight incidence lanes s=44−h exactly four high incidences go to each branch outside S and none to S. Any low branch inside S must then have all four of its changed neighbors among low branches in S: outside low branches use only high neighbors. Thus if any closed low branch exists there must be at least five, and exactly five would force a K5 whose five quartets violate the catalogue. Consequently all closed branches are high in the (39,5) and (40,4) lanes; in (38,6) and (37,7) there are either no closed low branches or at least six.

These are numerical survivors, not compatible complete charts. They reduce search to small paired activity sets; in particular, s=7 requires exactly four antipodal pairs per closed branch.

## 4. The exact change for a constant-bit spine interval

In F1 take an interval of t≥3 consecutive spine sites with identical root bits. It is a homogeneous module: every outside vertex sees it uniformly. Reverse the first-to-last arc. Let H_j denote the count with that whole module replaced, in the same position, by a transitive j-site block having the same exterior contacts. The exact increase is

    2^(t−2) H₁ + (3^(t−2)−1−2^(t−2))(H₂−H₁).

Here is a path-run proof. Group each full path by the outside order and the positions of maximal runs in the interval, not by which vertices fill the runs. Uniform exterior arcs make validity of the patterns independent of the filling. There are at most two such runs: each separator between two runs must use the root of that spine. In fact, deleting that root removes the special seam bridge; the two sides are now an ordinal join and that spine is ordered outside its interval. A leave-and-return through exterior vertices is impossible in that deletion. Let W₁ and W₂ be the numbers of one-run and two-run patterns. Then H₁=W₁ and H₂=W₁+2W₂.

For one run the full interval gains 2^(t−2) paths: the original increasing path survives, and in a path traversing last→first each interior site can independently go before last or after first, in increasing order within its part. For two runs, only a run containing **both** endpoints and at least one interior site changes its number of fillings (the reversed two-site tournament still has exactly one path). Choose which of the two runs is special and put m interior sites there. For 1≤m≤t−3, there are binom(t−2,m) choices and 2^m extra spanning paths in that run; the other run has one each. This gives 2(3^(t−2)−1−2^(t−2))W₂, proving the formula. The reverse operation has the negative of this increment. The literal DP replay checks 360 cases, including 336 with a nontrivial two-run contribution.

In particular the only unit pair in this class is **1↔3**: the nonnegative formula forces t=3 and H₁=1. The contracted tournament is transitive (any directed triangle gives three order-preserving path extensions by vertex insertion), so the original with its expanded transitive module also has count one. This explains why neither the adjacent edit lanes nor the homogeneous nonadjacent lane can provide unbounded forks. Mixed-bit intervals still need a separate treatment.

Replay with Node: 'node exterior_closure.js'. It needs no presumed full graph. The quartet input itself can be regenerated with 'node minimum_trades.js ../moore57_quadratic_tensors_20261125/period_fourteen_site.txt'.
