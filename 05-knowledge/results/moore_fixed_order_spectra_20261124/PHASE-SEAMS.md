# Native +2 seams of the Hoffman–Singleton phase chart: an exact inverse

The [14-site affine phase code](EXACT-PHASE-CHARTS.md) has two readouts: the original tournament and the marked tile complement. A single apex-contact flip commutes with tile complementation when it is off the spine. The following are changes of **one literal arc** with both code words valid; the decoded unlabelled graph is Hoffman–Singleton at both ends.

| Payloads | Source counts | Marked-mate counts |
|---|---|---|
| 11,43 | 3053,3051 | 25443,23899 |
| 13,45 | 3377,3379 | 29695,28437 |
| 19,51 | 3143,3141 | 27937,26135 |
| 21,53 | 3469,3471 | 32685,31169 |
| 67,99 | 3263,3261 | 32155,30353 |
| 69,101 | 3595,3597 | 37061,35545 |
| 89,121 | 3705,3707 | 36969,34943 |
| 93,125 | 3379,3377 | 26055,23185 |
| 6,38 | 4169,4181 | 45047,45049 |

All nine moves flip the bit of weight 32; the first eight are \(\pm2\) at the source, the last is \(+2\) at the mate. The two 3377/3379 edges are genuinely different marked lifts, **but not different pairs of merged source classes**: payloads 13 and 125 become converses, and so do 45 and 93. They have different marked exits because the marked chain-then-apex spine is not the transported spine under that converse operation. All 18 pairs of endpoint counts in the table were checked independently by ordinary subset DP.

## Why only the middle bit works

Let the word be \(1010\), a \(W\)-bit payload, and terminal \(0\). The [weighted inversion formula](../one_apex_join_representation_mass_pi2_over_12_20261106.md) expresses \(H-1\) as twice the sum of \(c_d\) over \(1,0\) inversion pairs at distance \(d\), where \(c_1=c_2=1\), \(c_d=2^{d-2}\) for \(d\ge2\). For a flipped bit with \(\ell\) sites to the left and \(r\) to the right, the \(0\to1\) change is \(2(B-A)\), where \(A\) sums the left 1 weights and \(B\) the right 0 weights. The first, second, third and fourth leftmost bits are fixed \(1,0,1,0\), and the final right bit is 0. Therefore
\[
c_\ell+c_{\ell-2}\le A\le 2^{\ell-1}-c_{\ell-1}-c_{\ell-3},
\quad
c_r\le B\le 2^{r-1}.
\]
If \(r\ge\ell+1\), \(B-A\ge c_{\ell-1}\ge2\). If \(\ell\ge r+1\), \(A-B\ge c_{\ell-2}\ge2\), unless \(\ell=4\), which would make \(W<4\). Thus there are no \(\pm2\) edges for odd \(W\ge5\), and for even \(W\ge6\) the flip must be central.

For equal \(\ell=r\), write \(A=c_\ell+s\). The bounds give \(2\le s\le c_\ell-3\). After the last site, which supplies \(c_\ell\), the right weights \(1,1,2,4,\ldots\) give exactly two bit patterns for every strictly interior sum, by binary expansion with two interchangeable units. Each of \(s-1,s+1\) is interior. This gives four suitable right words per left word, hence exactly \(4\cdot2^{\ell-4}=2^{W/2}\) labelled seams in the full payload cube (half increasing and half decreasing). The replay also covers the small exceptions \(W=3,4\), with one and three seams.

## A bounded-time exact target test, not just a phase test

For a central flip let \(u\) be the integer given by the left half, and \(v\) the integer given by the reversed complement of the right half. Both have \(\ell\) bits. Write \(G(x)\) for the one-apex height function and \(K(x)=\lceil x/2\rceil\). Then \(A=K(u)\), \(B=K(v)\). The involution exchanging \(u,v\) and complementing the central bit is the converse on the tournament (reverse the chain), explaining the two marked representations of the repeated seam.

The lower endpoint of a seam is exactly
\[
F(u,v)=2uv+G(u)+G(v)-1+2\min\{K(u),K(v)\}, \qquad |K(u)-K(v)|=1.
\]
Indeed within each half the weighted inversions give its own height minus one, whole cross inversions separate into the product \(2uv\), and inversions touching the center contribute \(2A\) or \(2B\). This furnishes an **absolute inverse** in this family. For prescribed lower odd \(H\), set \(s=\lfloor\sqrt{H/2}\rfloor\). Let \(a=\min(u,v)\) and \(b=\max(u,v)\). The unit condition implies \(1\le b-a\le3\). The elementary recurrence for \(G\) gives \(1\le G(x)\le x+1\) (for an even \(x\), use the recurrence bound and oddness). Consequently
\[
 2a^2 \le F(u,v)\le 2a^2+9a+5<2(a+3)^2.
\]
Test only \(a=s,s-1,s-2\), and \(b=a+1,a+2,a+3\), with both assignments to the halves. Retain those giving exact equality and satisfying the bit length, prefix \(1010\), and valid-payload restrictions. For the upper endpoint use \(H-2\). This is at most nine scalar tests, or eighteen ordered decoders, independent of the width; each evaluation reads \(O(\log H)\) bits. It also shows each absolute seam height has bounded multiplicity in the whole word family. Since \(u/2^\ell\in[10/16,11/16)\) and \(v=u+O(1)\), the lower counts concentrate in bands
\([25/32,121/128]\,2^{2\ell}+O(2^\ell)\). The next width multiplies the scale by four; these seams leave large gaps and cannot alone solve the full target problem.

## Independent order-eleven context

[Royle's own census report](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament) puts all odds through 80,555 besides 7 and 21 in \(S_{11}\) (pad with sinks if needed). His [maximal-count report](https://mathoverflow.net/questions/415359/maximum-number-of-hamilton-paths-in-a-tournament-on-n-vertices) gives 95,095 as the maximum, attained by Paley \(P_{11}\); this height was independently recomputed. The \(P_{11}\) affine automorphisms are transitive on unordered pairs, since \(-1\) is not a square, so to enumerate a nonempty flip shell it suffices to include edge \(\{0,1\}\). Through five flips the shell sizes in distinct heights are 1,1,14,239,2555,11875. The union supplies 1308 heights above 80,555 (508 through four flips). Thus \(41584\le |S_{11}|\le47546\). This neither reproduces Royle's full generator run nor classifies the entire remaining fringe.

## Replay

Run python3 phase-seams.py to check all source seam locations, the inverse on widths 5–12, the small exceptions and the nine literal chart squares. Compile paley-neighborhood.cpp with C++17 and optimization for the Paley shell counts. All charts in the first section are ordinary Hoffman–Singleton only; no degree-57 witness is assumed.
