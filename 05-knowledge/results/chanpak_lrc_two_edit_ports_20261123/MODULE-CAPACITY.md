# An absolute sparsity barrier for glued one-apex fillers

**PROVED** for profile bounds and the infinite profile exclusions below; **FINITE-EXACT** for the indicated counts and explicit primes. This is a limitation on the filler class, not on unrestricted two flips. The [literal scalable frames](TWO-APEX-SCALABLE.md) have starting counts \(6Q+5D\) and \(15Q+D\) where \(D=H(M),Q=Q_2(M)\).

## The two counts see different exponential scales

Strip ordinal singleton pads from a nontransitive one-apex module. Its reduced spine has \(n\ge2\) positions, starting O and ending I. Let \(d=j-i-1\) for an O at \(i\) and I at \(j>i\). In a two-slot cover the apex occupies one of the two slots. Sum the ordinary one-apex height over subsets of the spine assigned to that slot; the other slot has a unique transitive traversal. The constant 1 in each subset height contributes \(2^{n+1}\). For the chosen pair, fixing the endpoints in that slot, the contribution with \(u\) of the \(d\) intervening positions there is \(2^u+[u=0]\). This gives
\[
 Q=2^{n+1}+\sum_{\mathrm O_i,\mathrm I_j,\ i<j}2^{\,n-1-d}(3^d+1).
\]
In particular the extreme pair gives \(Q\ge 2\cdot3^{n-2}\). Also \(Q=O(3^n)\) uniformly: for \(d=n-2-k\) at most \(k+1\) pairs contribute, each at most \(4\cdot3^{n-2}(2/3)^k\). The \(2^{n-2}\) reduced words therefore grow with smaller base than their two-slot scale.

The profile \((D,Q)\) multiplies coordinatewise under ordinal join: in each slot concatenate the chosen paths in ordinal order, and conversely they cannot visit any block twice. Each singleton pad multiplies \(Q\) by 2 and leaves \(D\) fixed. For either lower height \(\le X\), \(Q\le X/6\); at length \(n\) the allowed pad multiplicities number \(O(1+\log(X/(12\cdot3^{n-2})))\). Summing the word counts backwards from the maximal length is a geometric sum with a linear weight, giving \(O(X^\alpha)\) distinct lower counts, where \(\alpha=\log_3 2\). Purely transitive fillers add only \(O(\log X)\).

## Even arbitrary ordinal gluing does not fix the deficit

Take one generator per distinct reduced one-apex profile, plus the singleton and any fixed finite bank of nonempty tournaments. At \(s=\alpha+\delta,\ \delta>0\), form the product \(\prod_g(1-Q_g^{-s})^{-1}\), which counts generator multisets with weight \(Q^{-s}\), overcounting the distinct products but never missing one. The logarithm is \(\sum_g Q_g^{-s}+O(1)\) as \(\delta\downarrow0\): the geometric tail is uniformly summable since \(Q_g\ge2\) and \(\sum_g Q_g^{-2\alpha}<\infty\).

We can halve the main bound by pairing each reduced word with its reverse-complement (the converse, with spine order reversed), which preserves \((D,Q)\). No odd-length word is fixed. At even \(n\) there are \(2^{n/2-1}\) fixed words, whose contribution is summable at \(s=\alpha\). Consequently
\[
 \log\prod_g(1-Q_g^{-s})^{-1}
 \ \le\ \frac{C}{\delta}+O(1),\qquad C=\frac{2^{-\alpha}}{2\log 3},
\]
using the geometric sum \(\tfrac12\,2^{-s}\sum_{n\ge2}(2\cdot3^{-s})^{n-2}\). Counting \(Q\le X\) by comparison with \(X^s Q^{-s}\), and taking \(\delta=\sqrt{C/\log X}\), gives at most
\[
 X^\alpha\exp(2\sqrt{C\log X}+O(1))
\]
profiles (and up to a factor two as many starts across both frames). The same proof permits any additional bank with \(\sum_g Q_g^{-\alpha}<\infty\), including every fixed finite bank. This is \(o(X/\log X)\); gluing cannot cover the primes even though the same class can pass any fixed collection of modular screens.

## Explicit escapes and an infinite joint-profile obstruction

The first scalar escape, across both frames, is already \(689\). Indeed \(Q\le114\) if a start is at most 689, and \(Q\ge2^{|M|}\) by choosing a Hamilton path on each ordered subset pair; thus the module has at most six vertices. Relabel along a Hamilton path and enumerate every remaining arc. The checker compares this universe against the glued closure and gives no escape below 689. Frame A realizes 689 with \(D=13,Q=104\), using module masks \((10,4,9,18,7)\).

For the bank of **just** one-apex strong blocks and singletons, the [exact checker](module-capacity-check.js) generates the entire multiplicative profile closure up to \(Q=825\) from 28 generator descriptions (16 distinct pairs). It finds 50 distinct products including the empty identity. Any profile representing 1747 or 4951 in either frame is within that cap, yet none does. The following actual nontransitive-spine modules nevertheless give both primes as exact two-flip starts:

| Module row masks | \(D,Q\) | Frame | Literal \(00,10,01,11\) |
|---|---:|---|---|
| \((30,60,40,48,4,17)\) | \((23,272)\) | A | \(1747,5831,9915,12229\) |
| \((30,44,24,48,2,21)\) | \((31,328)\) | B | \(4951,6795,25313,34657\) |

In fact \(1747\) is the **first prime** scalar escape across the two frames. For a nontransitive module on \(n\) vertices select a directed triangle. In the \(Q\) sum, the two disjoint classes where it lies wholly in one slot contribute at least 3 rather than 1 (the positive odd-cycle expansion gives \(H\ge3\) whenever a triangle is present). Hence \(Q\ge(3/2)2^n\). At starts \(\le1747\), \(Q\le291\), so nontransitive modules have \(n\le7\), and transitive modules are already glued. The exact checker enumerates all those orders anchored at a Hamilton path and compares with the glued closure; its only prime escape is 1747.

The checker confirms the literal path counts by subset DP and the cyclicity of each one-vertex deletion. It also replays the single-module-with-pads census through one million: A gives 1,159 starts (193 prime), B gives 655 (84 prime), union 1,813.

The paired readings escape the glued class along two **infinite** rays, not just at the initial examples. Because \(D=23\) or 31 is prime, a glued representative would require just one nontrivial component; all other factors are singleton pads. A reduced one-apex word with a given \(D\) has \(D>2^{n-2}\); the finite search gives only \(Q=296\) at \(D=23\) and only \(Q=644\) at \(D=31\). Their odd parts 37 and 161 disagree with 17 and 41. Appending any number \(j\ge0\) of singleton pads to the displayed modules therefore gives \((23,272\cdot2^j)\), \((31,328\cdot2^j)\) outside the glued profile set. This does not rule out the same scalar count being reached by a different profile on those rays.

**Next capacity test.** A two-apex spine with two fixed early positions both apices beat and two fixed late positions both apices lose to has \(4^n\) free contact words on \(n\) interior positions. Also \(Q=\Theta(4^n)\) uniformly: put one apex in each slot, end each before-apex run with its late guard and begin each after-apex run with its early guard. Independently send every free position to any of the four runs, giving the lower bound. Any two-slot cover is determined by the placement and order of the apices and allocation of all spine positions among the resulting four runs, giving the upper bound. Two apices thus escape the exponential **word-capacity** deficit. Distinct profile counts and an absolute inverse in that larger chart are still open.
