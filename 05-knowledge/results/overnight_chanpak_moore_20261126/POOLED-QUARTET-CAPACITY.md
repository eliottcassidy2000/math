# Pooled quartet capacity forces at least 129 changed branch matchings

2026-11-26. Continues [the Hall and exact-orbit checkpoint](HALL-AND-PAIR-ORBITS.md). **Conditional theorem on the actual three frozen stars:** two distinct full degree-57 Moore completions agreeing on those stars differ on at least **129 whole branch matchings**. This is a separation between hypothetical completions, not an existence or nonexistence proof. It is also a different metric from individual removed or added edges.

## 1. Capacity belongs to a quartet, not to a source branch

For a low branch (changed degree four), let $Q$ be its four changed neighbors. Each active source sheet uses an oriented trade from the regenerated catalogue. Its old profile is four actual target sheets. Distinct source sheets, even in distinct branches, cannot share that profile, since they would share at least two old neighbors and hence form a square. A quartet of multiplicity $M$ therefore supports at most $2M$ active source sheets **pooled over all low branches that use it**.

For $m$ visits into at most $a$ sheets, write $Q_{\rm int}(m,a)$ for the lower bound obtained by making integral loads as equal as possible. As before, each outside branch receives $m=sp$ visits from the $p$ active sheets per closed branch; for an outside low branch only active sheets can receive them. The ordinary quartets have $M\le7$. For any of the three exceptional quartets with $M=14$, if $j\ge2$ outside lows use it their total squared visits are at least

    Qint(jm,28) ≥ Qint(jm,14j) = j Qint(m,14).

The equality on the right follows by multiplying both the quotient-remainder counts and the number of slots by $j$. Reusing a quartet therefore never improves on the ordinary lower bound for those branches. Only a singly used exceptional quartet can improve it. Grant at most three of the outside lows a cap of 28 and cap the others at 14 in the lower-bound calculation. First remove any $J$ branch caps already charged at their exact activity size $p$; it is safe to grant all three concessions in the remaining problem.

## 2. Period four defeats the only arithmetic equality at 128

Let $h$ count high branches, $s$ count closed branches and $l$ count closed lows. The earlier closed-set bound and [Hall deficit bound](HALL-AND-PAIR-ORBITS.md), with the pooled variance test, leave only $(h,s,l,p)=(40,6,6,8)$ at 128. Among the eight tuples left by the earlier exact-mask screen, the gaps in increasing lexicographic order are

    +28224, +18816, +9408, 0, +3360, +40096, +18816, +9408.

At equality $m=48$, and all three exceptional quartets must be used singly among the outside lows. Each such outside low must use all 28 sheets: $Q_{\rm int}(48,27)=90>88=Q_{\rm int}(48,28)$. This saturates the associated quartet over *all* sources, so no closed low can use an exceptional quartet.

Every closed low has a changed link to another closed low: its four changed neighbors include at most three highs by closure, and no outside low can have a low changed neighbor. This link transfers its entire active set. In the regenerated catalogue, every ordinary quartet sends its active contacts into one residue class modulo four at each of its four targets, in each of the three frozen coordinates. The orientation change adds 28 and preserves the residue. Thus every closed low sends all eight visits into a single residue class at frozen branch 56, where the matching labels are the identity. Among six such sets, two use the same class. Sixteen visits into its 14 sheets force a repeat. But equality would require $Q_{\rm int}(48,56)=48$ in the frozen branch, so no repeated visits are allowed. This contradiction excludes 128.

The replay regenerates 2,114 distinct oriented profiles in 152 quartets and checks every target-coordinate residue set. The multiplicity histogram is $1:4,6:4,7:141,14:3$, and each coordinate has one residue for ordinary quartets and two for exceptional ones.

## 3. The cutoff layer at 129

Put $k=108+t$ and $t=21$. Apply the existing high-incidence bound $\max(2,54-h-\lfloor t/2\rfloor)\le s\le7$, with $h\le2t$. In the old exact chart, closure implies at least $p$ active sheets in any outside branch. There are at least five outside lows; if $p>14$ they would require distinct exceptional quartets, more than three. Thus run the variance test at $s-1\le p\le14$ with even $p$. The replay uses the deficit and exact-cap formulas in the earlier notes.

The arithmetic leaves 114 tuples, of which four fail by the same saturation and mod-four argument: $(40,6,6,8)$, $(41,6,5,8)$, $(41,7,7,8)$ and $(42,7,6,8)$. In each, there are at least five closed lows with eight visits each; the total visits $m=48$ or 56 still cannot repeat at equality. Both values force saturation because $Q_{\rm int}(m,27)>Q_{\rm int}(m,28)$.

There are **110 necessary tuples across 11 rows**, all with $(h,s)$ of the form $(40,5)$, $(41,3\ldots6)$ or $(42,2\ldots7)$. None has been shown to lift to actual matchings. The replay lists the exact remaining tuples and their variance gaps.

Replay: `node pooled_capacity.js ../moore57_quadratic_tensors_20261125/period_fourteen_site.txt`.
