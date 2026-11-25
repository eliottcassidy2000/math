# Seating a fourth neighbor set on the 28-trade chart

**FINITE-EXACT.** On the three-star 28-trade seed, the exterior site \((1,0)\) admits **51 of its 53** required new neighbors simultaneously. Including its three assigned star neighbors and its branch head, this fixes 55 of its eventual 57 incident edges. A specified 53-choice reference has just two conflicting pairs. Every exact local completion would need to change the reference choices on at least **eight** outside branches, even allowing every possible slot on those branches. The existence of a complete 53-choice neighborhood on this seed remains open. Changing the seed permits several complete neighborhoods at once. [Completed neighborhoods](PLANTED.md) [Three-star construction](three_stars.py) · [Physical slot certificate](fourth.py)

## Three permutations of the same 53 choices

Use centers \(C=(56,0,28)\), and write \(p_h,q_h\) for the maps from exterior branch \(h\) to branches 0 and 28. The map to branch 56 is the identity. For a proposed new edge from \((b,x)\) to \((h,y)\), the three new two-step endpoints at the centers are therefore the triple \((y,p_h(y),q_h(y))\). At each center forbid the three destinations already reached directly or through another center. All selected triples must be pairwise disjoint in **each of their three coordinates**, with at most one slot selected per exterior branch \(h\ne b\). These conditions are sufficient for inserting the selected edges locally: a triangle through a new edge would require a common old center neighbor at its ends. A square involving one new edge would close through two old centers. These possibilities are exactly the pre-existing forbidden routes. Two new edges have the common endpoint \((b,x)\), and could produce a square only when the two new neighbors share a center. Heads cannot give such a collision because the new neighbors are in different branches. All old triangle and square tests already pass. No matchings between other exterior branches have yet been filled.

For \((b,x)=(1,0)\), `PICKS` in the replay gives one eligible slot in each of the 53 possible branches in increasing order. The only conflicts are the pairs of branches \(\{6,17\}\) and \(\{40,53\}\). Discard choices 17 and 53 to obtain the claimed 51 new edges. To test all repairs on at most seven rows, enlarge the allowed repair support to exactly seven, permitting each selected row to retain its reference value. Every admissible support must meet both conflict pairs. There are exactly 8,493,464 such supports. For each, freeze the other 46 rows, remove the slots they occupy in each coordinate, and recursively try all eligible triples of each of the seven free rows. No support succeeds. This is a Hamming-distance certificate for the specified reference, independent of how it was found. [Complete radius-seven search](fourth_repair.cpp)

## A joint cost in the two voltage coordinates

**PROVED selection rule.** For *every* source branch \(b\) and slot \(x\) on this seed, the two coordinate displacements of a completed neighborhood satisfy
\[
 \sum_{h\ne b}\big(y_h-p_h(y_h),\ y_h-q_h(y_h)\big)
          =(28,28)\pmod {56},
\]
where the sum runs over the other exterior branches. Write \(P=p_b(x),\ Q=q_b(x)\). The forbidden triples of slots, separately in the three center boards, are \(\{x,P,Q\}\), \(\{P,x,Q+28\}\), and \(\{Q,x,P+28\}\), with arithmetic modulo 56. The selected coordinates must exhaust the complements of these three sets. Subtracting their sums gives the identity. It retains the joint pair of changes performed on a single physical row.

Compare a chosen cell against the translation pair of its own row in the untraded starting seed. Each row offers a finite set of displacement pairs. Assign cost zero to its original pair and cost one to any other option. For each fixed source branch add the row options in \((\mathbb Z/56)^2\), maintaining separate sets of reachable pairs for costs zero through seven. The first cost at which \((28,28)\) occurs gives the following exact lower bounds. Allowing all cells in each row relaxes the slot exclusions and distinctness constraints, so the numbers apply simultaneously to every physical source slot.

| Minimum traded cells forced per source slot | Number of eligible source branches |
|---:|---:|
| 3 | 1 |
| 4 | 4 |
| 5 | 45 |
| 6 | 4 |

Branch 38 has the first bound. Branches 2, 12, 35 and 42 have the last. For a complete new star the 56 source slots use each target cell in each included row exactly once, so at least 336 such changed placements are required at each of the last four choices. This is a necessary cost relative to the original seed, with no sufficiency claim for the relaxed row selections.

From this directory run:

~~~sh
python3 -B fourth.py
c++ -std=c++17 -O3 fourth_repair.cpp -o /tmp/moore57-fourth-repair
python3 -B fourth.py --board | /tmp/moore57-fourth-repair
~~~
