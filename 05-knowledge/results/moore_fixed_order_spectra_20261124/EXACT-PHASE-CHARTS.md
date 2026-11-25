# Exact Hoffman–Singleton codes and longer fixed-order prefixes

This extends [the fixed-order census](README.md). These are additional **affine phase encodings**, not exact counts of the native 50-site graph spines. Full spectra beyond order ten are not enumerated.

## An exact 50-site chart, and compression to order 14

Use vertices \(U_{i,x},W_{j,y}\), for \(i,j,x,y\in\mathbb F_5\), with steps \(\pm1,\pm2\) in the respective cohorts and cross edges \(y=x+a_{ij}\). By shifting the labels inside cohorts, normalize \(a_{0j}=a_{i0}=0\). The remaining 16 base-five entries need 38 bits. For an arbitrary payload use a transitive chain with one apex, whose outcontact word is \(1010101010\), the 38 bits, then \(0\). Along the displayed spine (chain followed by apex), all consecutive arrows are forward. The standard table \(a_{ij}=ij\) gives code \(47503389336\).

The apex is uniquely identifiable in the unmarked class: deleting it leaves a transitive tournament, and the prefix has two triangles sharing only the apex, so no other deletion does so. The remaining order is unique. The decoder reads the word in that order, reconstructs the phase table, and checks the Moore condition. These charts are strongly connected since their first bit is 1 and last bit is 0. For the standard word the source is rigid and nonselfconverse: a converse would force equality with the reverse complement, which has a different prefix.

A smaller code is possible when the table factors as \(a_{ij}=u_i v_j\), with \(u,v\) permutations of the nonzero coordinates after the leading zero. The multiplicative gauge \((u,v)\mapsto(cu,c^{-1}v)\) lets us insist \(u_1=1\), after which the pair is uniquely recovered from the table (read the first nonzero row for \(v\), then any nonzero column for \(u\)). This leaves \(6\cdot24=144\) labelled normalized tables, all presentations of the same unlabelled Hoffman–Singleton graph. Encode the two Lehmer ranks together as one byte; use prefix \(1010\) and final \(0\). This takes exactly **14** vertices. Within this fixed-sentinel, fixed-width binary layout, eight payload bits are necessary to distinguish the 144 tables, not a minimum among general tournament encodings.

| Order | Standard phase chart | Exact source \(H\) | Exact marked tile mate \(H\) |
|---:|---|---:|---:|
| 14 | scalar-normalized product phases | 5,119 | 65,647 |
| 16 | product phases before removing scalar gauge | 20,479 | 458,141 |
| 50 | arbitrary normalized affine table | 356,498,466,438,421 | 129,545,487,780,991,322,017 |

The replay evaluates the [word and mate transfers](../moore_merged_count_addresses_20261124/word_code.py) and independently verifies both ends by ordinary subset DP for orders 14 and 16. Across the 144 order-14 sources there are 134 heights (125 singletons, eight doubles, one triple), whereas **all 144 mate heights are distinct**. So the source height alone loses phase information, but in this finite family the **mate height identifies the exact phase table**, not merely the unlabelled graph. The source range is 2663–5365 and the mate range 5533–74985: no scalar overlap. The 144 source classes are distinct by apex-and-order decoding, and the 144 mate classes are distinct by their counts. The mate identification is within this chart family; it is not a decoder for arbitrary tournaments.

## Certified initial intervals and exact ranks

For each row, all positive odd values up to the displayed bound except 7 and 21 occur at that order and at every later order by adding sinks. These are certified **lower bounds on the end of the initial interval**, not claims about the first missing value of the full spectrum.

| order | certified through |
|---:|---:|
| 11 | 15,745 (previous replay) |
| 12 | 23,929 |
| 13 | 51,857 |
| 14 | 137,583 |
| 15 | 471,341 |
| 16 | 2,035,707 |
| 17 | 5,394,339 |
| 20 | 6,164,957 |
| 40, hence 50 | 4,294,967,295 (previous independent certificate, order bound below) |

The new replay generates the complete ordinary class bank through order eight with nauty; for orders nine and ten it adds fixed samples of 8192 and 4096 masks, respectively. A rooted bank record is \((h,e,R)\): total, ending-at-root, and two-path split with the second path beginning at the root. The rooted banks have 254 distinct triples at order six, 2732 at seven, 51268 at eight, 71921 in the order-nine sample, and 40940 in the order-ten sample. The [one-arc join identity](../../../01-canon/theorems/THM-4866-single-cycle-floor-bottom-band-uniquely-hamiltonian-one-arc-lorentzian-joins-and-crevasses.md) gives \(h_xh_y+R_xR_y-e_xe_y\) using the converse of the second bank tournament. Include two one-apex words, one word plus any sampled rooted triple, and joins of two triples with cut sizes restricted as in the replay. Appending a bit to the numeric apex word gives \(g(2v+1)=g(v)\), \(g(2v)=g(v)+v+(v\bmod2)\), with \(g(0)=1\); the apex-rooted triple is \((g(v),[v\text{ even}],v+1)\). The exact sieve records the smallest order among its witnesses. For example, one-apex-based joins already give the whole prefix through 305381 at order 16; including the two finite banks extends it to 2035707. The stopped sample is reproducible and not exhaustive beyond order eight. Every positive mark comes from a constructed tournament.

For the old [independent 32-bit certificate](../../../session203/coverage-recertification/README.md), its F1 joins satisfy \((u+1)(v+1)\le H\). Therefore the two word widths, with a zero word taking zero bits, total at most 33 for \(H<2^{32}\); two apexes give order at most 35. In F2, \(g(v)y_h-[v\text{ even}]y_e\ge0\) and \(y_R\ge1\), so \(v+1\le H\). One apex, at most 32 bits, and at most seven rooted sites give order at most 40. This is an order bound on the previously verified witnesses, not a rerun of the 32-bit sieve.

The entire order-14 phase packet now has exact scalar ranks: by Rédei parity and the universal omissions 7 and 21, each height \(h\) in its ranges has rank \((h+1)/2-2\). For the standard source and mate the ranks are **2558** and **32822**. The corresponding order-16 values have ranks **10238** and **229069**. These ranks stay the same at 50 if each endpoint is padded individually (not a claim that the padded endpoints remain tile mates). In a padded source chart, the strong source code is recoverable as its unique non-singleton strongly connected component.

## Ordinary classes, and both conditional degree-57 endpoints

A Burnside bound transfers the earlier labelled percentiles to uniform **ordinary classes**. For any invariant event, write its labelled probability as \(p\) and the nonidentity Burnside correction for all tournaments as \(\delta\), divided by the total number of labelled tournaments. Its ordinary class share is \((p+\delta_E)/(1+\delta)\) with \(0\le\delta_E\le\delta\), so the change is at most \(\delta\). A nontrivial tournament automorphism has odd order. With support size \(s\) every unordered edge meeting that support has orbit length at least three (an odd permutation cannot swap endpoints while stabilizing an unordered pair). Thus \(\delta\le\sum_{s=3}^n n^s 2^{-s(2n-s-1)/3}\). For \(n=50\) the \(s=3\) term is less than \(1.58\cdot10^{-24}\); for \(s\ge4\) the logarithm of the summand is convex, so bounding by its values at \(4,50\) gives less than \(2\cdot10^{-24}\) in total. The three previous percentile statements therefore also hold for ordinary classes (the strict certified margins exceed that error). The two exact 50-site phase-chart values are both below \(2^{-98}\) times the order-50 labelled mean, putting more than 96.15% of ordinary order-50 classes above each by the same Cantelli estimate.

The degree-57 bound also applies to the **dense tile mate**. Take its converse: the backward graph pairs now are the original sparse pairs together with the two incident spine arcs where present. So at most 59 flips meet any vertex relative to the forward order, and converse does not change \(H\). Inserting \(59\) in the previous elementary Brégman estimate gives a gap exceeding 974 bits (use \(\log 2<0.6932,\log 3250<8.087,\log59>4.077,\log(3250/59)<4.01\)). Thus both conditional endpoint counts are below \(2^{-970}\) of the labelled order-3250 mean. For that order the earlier variance formula has relative variance \(<0.000617\): its \(j=1\) term is \(2(n-2)/(n(n-1))\), and for \(j\ge2\) bound the numerator by \((2(n-2j))^j\) by arranging \(j\) dominoes in the \(n-2j\) gaps between unselected edges; the direction-choice factor is at most \(2^j\). The denominator is at least \((n-2j+1)^j(n/2)^j\). This bounds the tail by a geometric series with ratio \(4/n\). More than 99.938% of ordinary order-3250 classes lie above either hypothetical native height, again by Cantelli and Burnside. The underlying graph's existence and the exact endpoint values remain open.

## Replay

Run python3 hscodes.py. For the prefix, compile order-prefix.cpp with C++20 and optimization, run gentourng -q n > bank-n for n=2,..,8 (default text records), and then run order-prefix 20 1 in the same directory. The optional second argument multiplies the fixed sample sizes; the receipt uses 1. Counts and ranks are exact within the indicated interval, not extrapolations to untested targets.
