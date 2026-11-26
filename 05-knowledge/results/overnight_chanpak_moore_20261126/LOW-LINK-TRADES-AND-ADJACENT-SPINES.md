# Frozen link trades and the two adjacent-spine unit edges

This is the next [overnight Chan–Pak / Moore checkpoint](README.md), after [opposite-arm forks and boundary variance](OPPOSITE-ARMS-AND-BOUNDARY-VARIANCE.md). The Moore claims concern separation between distinct hypothetical full completions with the actual frozen stars, not existence. The spine result is an all-order theorem within F1.

## 1. At least 119 links if one exterior sheet stays fixed

Retain the [frozen-star difference identity and 108-link bound](../moore57_quadratic_tensors_20261125/PRIMITIVE-HOLONOMY-AND-EDGE-SEPARATION.md). A link is one whole branch matching. A branch is touched if any of its links changes; a sheet is active if any incident edge changes. A chart sharing four or more full stars already needs at least 133 changed links, by the [old frozen-link estimate](../moore57_quadratic_tensors_20261125/WILSON-ACTION-AND-RIGID-MOORE-CHARTS.md), so consider exactly the three frozen branches. All other 54 branches are touched. Every active leaf changes at least four exits. Call a branch **low** when it has exactly four changed links and **high** otherwise.

Suppose there are k=108+t changed links, with t≤10. The excess changed-link degree is 2t, so if there are h high branches then h≤2t, and the high side has 4h+2t incidences. A low branch with no low changed neighbor uses up four of those incidences. The number s of low branches with a low changed neighbor is therefore

    s ≥ 54−h−floor((4h+2t)/4) ≥ 54−4t−floor(t/2) ≥ 9.

On any of these s branches, every active sheet changes all four exits; its old removed neighbor in the adjacent low branch also has removed degree four. The inactive-neighbor consequence of the difference identity bounds the active sheet's old inactive degree by four. Its head and the three frozen exits already account for four, so all its 53 old exterior neighbors are active. In particular the old matchings among the s branches preserve the active and inactive sets. The active sheet count p is the same in all s branches, by bijectivity and closure in both directions, and p>0.

If p<56, the active and inactive parts of those s branches are disjoint nonempty (s−1)-regular graphs. A nonzero centered linear combination of their indicators, extended by zero, has Rayleigh quotient s−1≥8; the old chart has nontrivial upper eigenvalue 7. Thus p=56. A fully active branch with this closure is matched into every exterior branch with no inactive exterior neighbor, so *all* exterior sheets are active. The contrapositive proves **k≥119 if even one exterior sheet stays inactive**, in particular for two distinct completions agreeing on the fourteen planted exterior neighborhoods. This is a different restriction from the previous bound on the number of individual edges.

## 2. Exact four-contact trade catalogue on the physical fixture

The three frozen matchings of an exterior sheet give a row of four coordinates: the branch and the three slots (56,0,28). Any two distinct rows agree in at most one coordinate: within a branch the three maps are injective, and two coincident center values in different branches would give a square. A four-contact trade is a pair of disjoint tetrads with identical coordinate sets, with four distinct entries in each coordinate of each tetrad. Four is the smallest possible nonempty side under these hypotheses: each old row must meet a new row in each of the four coordinates, and those four new rows must be distinct by linearity.

The [fixture replay](minimum_trades.js) independently reconstructs the three frozen matchings from the actual seed and both non-periodic repairs. It finds the following unordered trades:

| Board | Contact rows | Trades | Four-branch sets | Rows covered | Max trades per quartet |
|---|---:|---:|---:|---:|---:|
| Full three-star fixture | 3,024 | 1,057 | 152 | 2,478 | 14 |
| Marked contact board | 1,741 | 83 | 33 | 484 | 5 |

**Every** trade flips all three slot coordinates by +28 mod 56 on each of its four branches. Slot matching across the sides, relative to branch matching across the sides, gives a double transposition in each of the three coordinates. Exactly three quartets attain fourteen trades: (1,2,36,37), (1,31,44,54), (2,23,38,54). No five-branch set contains more than two catalogued quartets.

The enumeration uses a direct table for two-coordinate keys. For two disjoint old rows A,B, each new row must meet both of them, because each old row meets four new rows in four distinct coordinates. A new tetrad therefore determines a bijection of which coordinate it hits on A with which it hits on B; it has no fixed point, since A and B used different values in each coordinate. Try the nine derangements, reject repeated values in the new tetrad, then split the two remaining values per coordinate into the other two old rows (eight ways up to swapping). Require A,B to be the two least old row IDs and A to be the least on either side. This counts every unordered trade exactly once. The replay asserts the linearity, half-turn and permutation properties and checks the integer totals.

**In a full chart:** for any active sheet at a low branch, its four removed and four added contacts form precisely such a trade, with quartet equal to the changed-neighbor list of that branch. The difference identity against the frozen centers gives equality of the slot sets, and repeated slots would create a square. An oriented trade cannot be reused by two sheets of the branch, by injectivity on any of the four links. Thus a low branch has at most 28 active sheets. Each changed matching at it acts, on its old image sheets, by a subset of antipodal transpositions, so its active size is even. At k≤118 some low branches necessarily exist, and hence some exterior sheets are unchanged. Therefore k≥119 holds for this *specific* three-star fixture without assuming the optional fourteen fixed neighborhoods.

## 3. If k=119, at most eight active sheets in the coupled low branches

This adds a restriction even after reaching the 119-link threshold. Suppose k=119, so total excess degree is 22. As before write h for the number of high branches, l=54−h≥32 for the number of low branches, and s for the number of low branches with at least one low changed neighbor. The previous counting gives s≥5. If s≥9 the spectral argument forces full activity, contradicting the trade cap. If s=5, necessarily h=22. Indeed if e_LL counts the low-low changed links then

    4l−2e_LL ≤ 4h+22, so e_LL ≥ 97−4h.

With five coupled low branches this forces at least nine internal edges and h=22. Their graph is K5 minus at most one edge; at least three vertices have their entire four-neighbor quartets within this five-set, and those quartets are distinct. The catalogue's five-set maximum of two rules that out.

There is also a second-moment constraint at s=6,7,8. On the s coupled low branches define f to be 1 on active sites, −p/(56−p) on inactive sites, and zero elsewhere, with common activity p as above. Then f is centered (even within each branch), (L f)=(s−1)f on the s branches, and ||f||²=56sp/(56−p), for L the old exterior adjacency. The spectrum in [−8,7] implies

    ||(L−(s−1)I)f||² ≤ (8−s)(s+7) ||f||².

This is the usual second-moment bound from λ²≤56−λ, using mean λ=s−1. For an outside low branch, every inactive sheet sees only inactive sheets in the s branches, by closure. At those sites Lf=−sp/(56−p), and there are at least 28 such sites. Also Lf sums to zero on that entire outside branch, since branch links are matchings. Cauchy–Schwarz bounds its norm contribution below by 56(sp/(56−p))². After summing over l−s branches and dividing by ||f||², we get

    (l−s)sp/(56−p) ≤ (8−s)(s+7).

This excludes s=8. For s=6 and 7 it gives respectively p≤8 and p≤4 (p is even). Combining e_LL≥97−4h with e_LL≤2s gives h≥21, with h=21 possible only for s=7. In particular the 119-link case needs only the following numerical lanes:

| high branches | coupled low branches | common active size in those coupled branches |
|---:|---:|---|
| 21 | 7 | 2 or 4 |
| 22 | 6 | 2, 4, 6 or 8 |
| 22 | 7 | 2 or 4 |

These are necessary conditions, not witnesses. They give a small activity parameter for a search of the equality case, conditional on existence of the two full graphs.

## 4. Adjacent spine contacts give only 3↔5 and 9↔11

In F1, reversing an edge between consecutive spine vertices preserves a transitive spine after swapping their order, but swaps the two address bits. Equal bits leave the count unchanged. An unequal swap in the increasing numeric direction is v_0=4ℓq+ℓ+r → v_1=4ℓq+2ℓ+r, 0≤r<ℓ, with ℓ a power of two. The recurrence gives F(4q+2)−F(4q+1)=q+2−[q odd]. Padding the common suffix adds Z−ε, where Z is its padded reverse complement and ε=Z mod 2 if ℓ>1, zero otherwise. Hence

    ΔH=cℓ+F(c−1)(q+2−[q odd]+Z−ε)−[ℓ=1][c odd].

The last term accounts for the changed root endpoint parity at the lowest pair. The bracket is at least two, so a unit change is possible only at c=1, ℓ=1, q=0 or 1. This gives b=2↔3 and 6↔7 and the two total-count edges stated above. At these endpoints the other arm (address one) can have a unit spoke only at b=2, and there its +2 sign agrees with the adjacent reversal. Thus these edits provide **no opposite-sign adjacent-swap/spoke fork**. Nonadjacent spine reversals remain unclassified.

Replay: from this directory run 'node minimum_trades.js ../moore57_quadratic_tensors_20261125/period_fourteen_site.txt'. It uses only the existing seed, not a presumed full chart.
