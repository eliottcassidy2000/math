# Exact sevenfold squares at equal Jacobian: the third-slot direction

**PROVED** for the slot identities and infinite fixed-Jacobian fibers. **FINITE-EXACT** for the frame census and literal certificates. The fixed-apex-edge/two-contact template is the one in the [seven-site census](TEMPLATE-CENSUS.md); its bounds do not apply to unrestricted edits.

Write \(D,Q,R\) for the ordered one-, two-, and three-slot path-cover readings of a nonempty filler, with empty slots allowed. The first lower read involving \(R\) occurs at **eight fixed spine sites**. Tags are \(\mathtt{000210230}\), filler placeholder index 6, marked contact bits 1 and 5 (the two \(b\) contacts to sites 0 and 2, in the original placeholder indexing). The four literal reads are
\[
 11D+42Q+20R,\quad71D+258Q+116R,\quad23D+90Q+44R,\quad77D+294Q+140R.
\]
In order \(00,10,01,11\), the ordered nonempty run skeleton counts are \((155,102,20),(935,606,116),(335,222,44),(1085,714,140)\). Evaluate these against \(D,Q-2D,R-3Q+3D\). The two orders are strictly increasing for all nonempty fillers. The half-mixed gain is \(d=-3D-6Q\), so the normalized gain Jacobian is constant with diagonals zero and off-diagonals \(d\). It does not depend on \(R\).

In fact the entire square depends on just the height and the mixed gain, not three independent linear readings. If \(h=H_{00}\), then
\[
 (H_{10},H_{01},H_{11})=((29h-12d)/5,\ (11h+2d)/5,\ 7h).
\]
The three-slot variation keeps the Jacobian fixed but changes the affine gain map's intercepts as well as the target heights. This is stronger than matching the determinant alone, but not equality of the two gain maps. In the [older three-slot inverse](../chanpak_induced_gap_three_slot_descent_20261118/README.md), a fourth slot was hidden from the first three readings. Here the next slot is hidden specifically from the differential of a literal two-flip map.

## Literal primes and independent fixed-Jacobian directions

| Filler outgoing masks | \(D,Q,R\) | Heights \(00,10,01,11\) |
|---|---|---|
| \(86,284,392,401,484,79,398,291,33\) | \(2393,38272,317403\) | \(7981807,46862827,17465251,55872649\) |
| \(6,244,344,51,417,69,153,301,107\) | \(2393,38272,317643\) | \(7986607,46890667,17475811,55906249\) |

Both starting counts are prime (trial division through their square roots); their difference is 4,800. Both marked Jacobians have off-diagonal entries \(-236811\). The full hosts have 19 vertices, with marks \((0,18),(2,18)\); the lower hosts are not strong. Ordinary subset/end-vertex DP replays all eight heights in [the literal checker](third-slot-check.cpp).

The least maximum filler order at which equal \((D,Q)\) can hide different \(R\) is eight. The complete [canonical extension replay](third-slot-profile-check.js) through order eight finds exactly one split reading: masks \((94,60,184,112,32,193,150,27)\), \((30,92,216,48,32,199,153,27)\) give \((345,5640,45225)\), \((345,5640,44505)\). Their starts in the eight-site frame are 1,145,175 and 1,130,775, composite, both with \(d=-34875\). The canonical enumeration extends each representative by a node, refines its classes by their outgoing counts to every current class, and minimizes adjacency encoding over permutations left within those classes. It gives 456 and 6880 types at orders seven and eight.

Slot counts multiply coordinatewise under ordinal join because each slot follows the join order. At any fixed product length \(t\) of the prime-example fillers, varying the number of factors of the first sort produces \(t+1\) distinct \(R\)'s, and hence distinct starts, at fixed \(D,Q,d\). Using both twin pairs strictly enlarges this: take \(u\) factors from the order-eight pair and \(v\) from the order-nine pair. The ratios \(45225/44505\) and \(317403/317643\) are multiplicatively independent, as the first has nonzero 67-adic valuation and zero 113-adic valuation while the second has zero at 67 and nonzero at 113. Consequently there are exactly \((u+1)(v+1)\) distinct starts in this product grid on one common order \(8u+9v+10\), all with identical \(D,Q,d\). The [14-pair exact bank](third-slot-bank.py) extends this to fourteen independent directions. Each pair has a prime that divides exactly one of its two \(R\)'s and none in any other pair (verified together with all 28 literal profiles). Take \(t\) factors from each pair: those prime valuations recover the fourteen independent choices of multiplicity, so \((t+1)^{14}\) distinct starts share the same \(D,Q,d\) and host order \(125t+10\). This does not claim that the product starts are prime.

## A structural reason the third slot can be invisible

In a skeleton contributing to the coefficient of \(R\), all three runs have nonempty module parts. A fixed site before the module can then meet an apex only immediately **after** it (at the start of a run); a fixed site after the module can meet it only immediately **before** it (at the end of a run). Thus two marks on the same apex and same side compete for one boundary, and no three-occupied-run skeleton uses both. Their mixed \(R\) coefficient is identically zero before imposing the sevenfold equation. Both primitive frames here use this geometry.

More generally the \(R\)-coefficient is \(2L_{\mathrm{side}}R_{\mathrm{side}}\). On the left of the module the three run assignments only check the first contact after each apex (use the module itself if no left site is assigned to that run); on the right they only check the last contact before each apex. The decisions are independent. Reversing the order in which the apices are read relabels the two constrained runs on each side without changing their assignment counts, giving the factor 2. Consequently when the marks lie on opposite sides, the mixed \(R\)-coefficient is \(2\Delta L_{\mathrm{side}}\Delta R_{\mathrm{side}}\); if the lower \(R\)-coefficient is positive, exact sevenfold imposes the separate ratio product \((L'_{\mathrm{side}}/L_{\mathrm{side}})(R'_{\mathrm{side}}/R_{\mathrm{side}})=7\). This is a direct filter for a future search for third-slot curvature, not a theorem that it always vanishes.

## Template closure and fixed-frame exclusions

The [prefix-sharing scan](template-scan-wide.cpp) uses 12 states for each apex order. As the spine is read, assign each site to the run before, between, or after the apices. The first run remembers whether it is empty or has a valid last contact (these cases share one state), or has an invalid last contact; the middle remembers empty or the validity of its last contact, and its first contact is tested when the run starts. The final run remembers only whether it has started, checking its incoming contact once. An empty middle is legal only in the order of the direct apex edge. Width-one values are shared across all placeholder choices. As before, the first three transitive filler widths determine each triple of coefficients, and the nonempty-slot basis certifies all four strict edge inequalities.

| Fixed sites | Width-one survivors | All-filler ratio and both increasing orders | Nonzero \(R\) in \(H_{00}\) | Nonzero \(R\) in \(d\) |
|---:|---:|---:|---:|---:|
| 7 | 6752 | 1360 | 0 | 0 |
| 8 | 16560 | 3746 | 2 | 0 |
| 9 | 37334 | 8814 | 10 | 0 |
| 10 | 78696 | 18902 | 28 | 0 |
| 11 | 156770 | 37990 | 60 | 0 |

Through eleven sites every \(R\)-dependent lower row also has \(3\mid d\) for all fillers (its \(d\)-coefficients are multiples of 3). At nine sites, the two old-form possibilities are the eight-site form and its singleton pad (which scales \(Q\) by 2 and \(R\) by 3). The new primitive is tags \(\mathtt{0003311301}\), placeholder 7, bits 0 and 2, with four rows of coefficients
\[
 (-7,158,122),\quad(-23,842,610),\quad(-15,500,366),\quad(-49,1106,854).
\]
It has \(d=-9D-39Q\), and again the four reads lie in a plane with \(H_{10}=5h-4d/3,\ H_{01}=3h-2d/3,\ H_{11}=7h\). The finite table identifies no claim of an all-order obstruction to \(R\) entering \(d\).

Width two and three in the eight-site frame give prime starts 359 and 887, respectively, and width two in the nine-site primitive gives 1723. The literal checker with the full scan output also excludes 887 and 1723 from every \(R\)-free frame through nine sites, and 359 through eight sites. For every variable form not already excluded by divisibility, the compiled lower bound is \(h\ge6Q\) using \(D\le Q/2\). For a nontransitive \(n\)-vertex filler \(Q\ge(3/2)2^n\), by placing a directed triangle wholly in either slot (each such subset contributes at least three, rather than one). Thus at cap 1723 the nontransitive modules have order at most seven; enumerate their arcs anchored at a Hamilton path and handle the transitive powers of 2 separately. A constant-width form would have to have coefficient dividing the target prime. At nine sites the constant-width form \(359D\) does give 359.

There is also a uniform scale separation: every \(j\)-slot cover on \(n\) vertices has \(n+j\) ways to split one slot in two, and each resulting cover has at most \(j\) preimages (one per proposed adjacent merge). Hence \(Q_{j+1}\ge (n+j)Q_j/j\). In both primitive third-slot frames \(|d|/h=O(1/n)\) since \(R\ge(n+2)Q/2\). This explains why the height can gain a new scale while the Jacobian does not.

Replay using `g++ -O3 -std=c++17 template-scan-wide.cpp -o scan && ./scan 7 11 > scan.out`, then `g++ -O3 -std=c++17 third-slot-check.cpp -o check && ./check scan.out` and `node third-slot-profile-check.js` and `python3 third-slot-bank.py`.
