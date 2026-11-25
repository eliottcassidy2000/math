# 290-site Moore voltage code by port colours and interval ranking

This refines the [Fourier overlap decoder](FOURIER-OVERLAPS-AND-THREE-DELETIONS.md) and [marked-spine construction](MATCHING-FLIP-ZIGZAGS.md). The degree-57 graph is still hypothetical. The argument applies only to the order-11 branch, with its fixed pentagon, 25 moving port orbits and 270 remote orbits. The code concerns labelled tournament contacts; it does not identify an unmarked tournament from its Hamilton count.

## Sharper overlap bound: average over all nontrivial characters

Retain the three cases \(c=25,100,175\), write \(t=c/5\), \(r=(2065+t)/15\), and let \(m_{ij}\) denote the number of shifts on a pair of distinct moving orbits. The total number of shifts is \(M=8395-c\), and \(\sum_{i<j}\binom{m_{ij}}2=675-c\). The latter follows by [the exact gap census](ORDER11-ORBIT-PACKETS.md): an independent remote orbit uses all five unordered gap classes once, while a cyclic orbit uses three in cross edges. A port spends none on repeated shifts. Every cross gap is used at its two endpoints.

For an off-diagonal Fourier entry \(B_{ij}(\zeta)\), average its squared absolute value over the ten primitive characters. Distinct shifts have mean cross term \(-1/10\) per ordered pair. Hence
$$
\mathop{\rm mean}_{\zeta}\sum_{i<j} m_{ij}|B_{ij}(\zeta)|^2
 =\sum_{i<j}\left(m_{ij}^2-\frac{m_{ij}^2(m_{ij}-1)}{10}\right)
 \le M+2(675-c)-2(135-t).
$$
For the inequality, \(m^2(m-1)\ge4\binom m2\) for integer m. Thus the displayed subtraction can even be replaced by \((2/5)(675-c)\), which is the same \(2(135-t)\); any triple improves it further. The previous determinantal proof's trace, lifted quadratic bound and degree loss have the same values for every nontrivial character. Average the whole basis law over these characters. Some basis has the resulting promised overlap; nonsingularity of its principal minor at that character implies nonsingularity at every primitive character by Galois conjugacy. There is no character index to store.

Replace the previous upper bound on the Fourier error by \(U/225\), where \(U=M+2(675-c)-2(135-t)\) above. In the notation of that proof set \(Q=((295-c)49+268t)/225\), and take
$$
R=57r-\left\lceil\frac{11r^2}{100}-4Q-\frac{268t+U}{225}
                      +\frac{175+72t}{15}\right\rceil.
$$
This gives 6,042, 6,019 and 5,997 stored shifts respectively. Diagonal entries are still recovered from the stored off-diagonal rows.

## Port-colour coding and an arithmetic decoder

For a fixed principal basis, let \(s\) be the number of its port orbits, \(s_i\) those over colour i, and \(k=r-s\) the number of its remote orbits. Write
$$
L=\frac{11}{2}[k(539-k)+s(49-s)].
$$
This is the number of possible shift positions from touching remote–remote and port–port pairs. In counting, dominate the choices on these pairs by \((1+x)^L\). A selected remote has one port contact per colour, with 55 possibilities (port name and phase). An unselected remote has at most one stored contact per colour, offering \(11s_i\) positive choices. The weighted inventory of all raw records is therefore at most
$$
Z(x)=(1+x)^L(55x)^{5k}\prod_i(1+11s_i x)^{270-k}.
$$
Evaluate at \(x=1/55\). For \(s=5a+b\), shifting one unit from an \(s_i\) at least two larger than another increases the product. The worst distribution uses \(5-b\) values a and b values \(a+1\). This gives the rational upper bound
$$
55^R Z(1/55)\le
\frac{56^L(5+a)^{(5-b)(270-k)}(6+a)^{b(270-k)}}
     {55^{L-R}5^{5(270-k)}}.
$$
The normalized bound divides by \(11^{294}\), as justified by the connected touching support in the Fourier proof.

**An actual finite arithmetic algorithm.** First encode only the support and multiplicity shapes, by consecutive subintervals in a mixed-radix probability split. A free pair with m shifts is assigned weight \(\binom{11}m x^m\) for \(0\le m\le3\). A selected remote's port-name options are equiprobable, and an unselected remote's absent or present choices have respective weights 1 and \(11x\) (one \(11x\) option per available port name). Each factor is normalized locally. These weights give exactly the intended raw model except for omitting free pairs with more than three shifts, so its normalizing denominator is no larger than \(Z\).

Then, within the selected shape interval, choose the phases uniformly in mixed radix. For non-tree pairs there are \(\binom{11}m\) possibilities. On the edges of a deterministic spanning tree of a connected support, choose only the canonical translate in each translation class, giving \(\binom{11}m/11\) choices. Because the tree has no cycle, one can normalize every tree edge simultaneously by changing the phase origins. This is exactly the factor \(11^{294}\) used in bounding probabilities. For disconnected formal shapes, simply use unnormalized radices (or reserve and discard their intervals); this never affects a valid record's interval. A normalized record with h stored shifts thus has interval length at least \(11^{294}x^h/Z(x)\). If the bound above is at most \(2^b\), the half-open interval has a dyadic point on the mesh of width \(2^{-b}\). Encode the smallest such point; the decoder retraces the same partitions. If original phases must be preserved, skip normalization and do not divide the bound.

The replay checks all 78 rational comparisons using integer powers, and implements both directions of the interval algorithm for a full-size **synthetic** normalized touching packet with all 6,042 allowed shifts. It tests actual shape storage, connected-tree normalization, phase rank, and the recovery of the same packet, not existence of a degree-57 graph. The [prior exact Fourier replay](fourier-overlap-check.py) gives an actual Hoffman–Singleton reconstruction control for the subsequent Gram completion step.

| Cyclic orbits | Basis size | Stored shifts | Normalized payload bits (before header) | Original-phase payload bits (before header) |
|---:|---:|---:|---:|---:|
| 25 | 138 | 6,042 | 41,846 | 42,863 |
| 100 | 139 | 6,019 | 41,751 | 42,768 |
| 175 | 140 | 5,997 | 41,661 | 42,678 |

With orbit names kept, add a 2-bit stratum and a 295-bit basis mask: **42,143 bits**, fitting 291 sites when phases are normalized; **43,160 bits**, fitting 295 sites if phases are preserved. If names may be permuted within each port colour and within the remote set, label selected sites first in each group. Five 3-bit selected-port counts replace the 295-bit mask (the selected remote count follows from r). That uses **41,863 bits**, fitting **290 sites**. These are upper bounds, not claims of optimality.

## Path guarantees and a reserved spine

On 290 sites the relabelled code leaves 42 free contacts. Use 39 to orient 13 vertex-disjoint directed triangles and one other contact as a marker. Both readings decode and have at least \(3^{13}=1,594,323\) paths by the [odd-cycle collection identity](https://arxiv.org/abs/2307.05569). With names preserved at 291 sites there are 52 spare contacts and 17 such triangles, giving \(3^{17}=129,140,163\) paths at both readings. The marker has no guaranteed displacement.

On **291** sites it is also possible to reserve a Hamilton spine for the relabelled code, since its nonspine contacts supply exactly the old 290-site capacity. Fix two of the eight-site marked-spine gadgets (42 nonspine contacts altogether). The chart and its tiling mate have at least \(24^2=576\) and \(23^2=529\) paths, by the previous concatenation argument. This is a smaller marked-spine variant of the 296-site code; it trades away some of that chart's stronger uniform lower bound. For all 290-site readings \(H<2^{1687}\); on 291 sites \(H<2^{1694}\). These upper bounds replay the previous finite Brégman calculation with the respective orders. None evaluates H of an actual degree-57 graph.

Run `python3 port-colour-code.py` for the capacity and codec replay.
