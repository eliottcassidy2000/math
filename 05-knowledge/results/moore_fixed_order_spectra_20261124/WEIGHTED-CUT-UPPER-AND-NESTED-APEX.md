# Weighted cuts cut the degree-57 upper exponent; a nested-apex exact guard

This extends the [overnight tandem](OVERNIGHT-TANDEM.md). The Moore claims are conditional on existence, uniform over every allowed complement spine, and refer only to the native tournament and its marked tile complement. Neither Chan–Pak nor Moore existence is resolved.

## 1. The new upper height window

For the native tournament \(T\) and its marked tile complement \(T^\ast\),
\[
 2^{18927}<H(T),H(T^\ast)<2^{27021}.
\]
The floor is the [weighted augmentation result](WEIGHTED-AUGMENTATION.md). The old upper exponents were 28860 and 28863. More precisely the new uniform computation gives \(\tfrac45\,2^{27021}\) as a strict upper bound, excluding \(E_p\pm1\) for \(p\ge13511\) at either native endpoint, in addition to the earlier small-exponent exclusions. Here \(E_p=2^{p-1}(2^p-1)\); no primality assumption is needed.

### Scaling paths without changing their weights

Adjoin a bidirectional vertex \(a\). Each Hamilton path closes uniquely through \(a\); its weight is unchanged if every arc \(i\to j\) is assigned weight \(w_i/w_j\). Every cycle cover also telescopes to weight one, so its permanent is invariant. The direction of scaling matters: smaller displayed positions get smaller \(w_i\), decreasing cheap forward weights and increasing expensive backward weights. The positions of the backward edges, not just tournament scores, now enter the estimate.

For any nonnegative matrix \(W\) of order \(N=3251\), whose \(i\)-th row sums to \(r_i\) and is bounded entrywise by \(m_i>0\), we use
\[
 \operatorname{per}W\le e^{-N}\prod_i(r_i+6m_i).
\]
In row \(i\), sample an independent Bernoulli indicator for each entry, with probability \(W_{ij}/m_i\); the expected permanent is the permanent of the matrix of probabilities. Apply [Brégman's \(0/1\) inequality](https://www.mathnet.ru/eng/dan37755) and take expectations rowwise using independence. An empty row receives factor zero. For a row of size \(k>0\), its factor is \((k!)^{1/k}\le(k+6)/e\): concavity of \(\log\) and the trapezoid bound give \(\log(k!)\le k\log k-k+s\), \(s=1+\tfrac12\log k\). For \(4\le k\le3251\), \(s/k\le1\) and \(k(e^{s/k}-1)\le s+s^2/k<6\), since \(s^2/k\) is decreasing and \(s<5.05\). The sizes 1,2,3 satisfy the inequality directly. Taking expectations proves the claim.

Use 25 consecutive bins of 130 sites each and \(w_a=4\); the exact \(w_i\ge1\) are constant on bins and listed in the replay. The per-row entry bound \(m_i=w_i\) is valid since every \(w\ge1\), including the apex. Allow a second positive multiplier \(z_i\) per bin, \(z_a=1\), and apply AM–GM:
\[
 H\le \operatorname{per}W\le
 \frac{\big(\sum_i z_i(r_i+6w_i)/(Ne)\big)^N}
      {\prod_i z_i}.
\]

### Bound many overlapping cuts simultaneously

For bins \(i<j\), begin with every arc forward. A graph edge flips that arc, changing the weighted sum in the numerator by
\[
 g(i,j)=z_jw_j/w_i-z_iw_i/w_j .
\]
Set \(g(i,j)=0\) for \(i\ge j\). Define the mixed difference
\[
 d(i,j)=g(i,j)-g(i+1,j)-g(i,j-1)+g(i+1,j-1).
\]
Summing over all graph edges gives exactly \(\sum_{i<j}d(i,j)e_G(P_i,S_j)\), where \(P_i\) is the prefix through bin \(i\), and \(S_j\) the suffix starting at bin \(j\). Indeed, a backward arc joining bins \(b<c\) is charged by all rectangles \(b\le i<j\le c\), whose mixed differences telescope to \(g(b,c)\). Graph edges inside a bin have no effect; all within-bin weights are the same.

For disjoint \(P_i,S_j\) of sizes \(x,y\), center both indicators and use \(A+\frac12 I\). It acts by \(57.5\) on the all-ones vector and has norm \(7.5\) on its complement, so
\[
 \left|e_G(P_i,S_j)-57.5xy/3250\right|
 \le 7.5\sqrt{x(1-x/3250)y(1-y/3250)}.
\]
Use the upper endpoint for positive mixed coefficients and the lower endpoint (or zero) for negative ones. Bounding separate bin pairs by the same spectral lemma would discard the larger prefix/suffix constraints; the telescoping is what lets them enter. For the converse of the marked mate, add the exact change from the spine arcs, \(\sum_i g(i,i+1)\), one across each bin boundary. The within-bin spine arcs again contribute zero.

All numerical claims use rational inequalities. The square roots are rounded upward at thousandth precision, \(e>27182818/10^7\) follows from a finite partial sum, and the final exponential comparison is done with integer powers. Taking every \(z_i=1\) already gives \(H<2^{27051}\); the nonuniform AM–GM step lowers this to \(2^{27021}\) at both endpoints.

## 2. A nested apex gives actual +2 seams outside both earlier variable-tail charts

For a strong binary word \(w\) (starts 1, ends 0), let \(U\) be its transitive chain with a single inner apex \(b\). As usual 1 means that \(b\) beats the corresponding site. Before the old six-site marked fork \(F\), insert \(U\) as an ordinal block under the **outer** common apex; the outer apex beats all chain sites of \(U\) but loses to \(b\). The sole qualifying full path of \(U\) that can lead directly into the outer apex ends at \(b\), so \(e_U=1\). After \(F\) append any transitive tail starting with an apex outcontact. The [exact endpoint determinant](ENDPOINT-AND-SIGNED-DECK.md) proves a physical one-arc gap of 2 regardless of the two flank lengths, because the tail has \(h_V=s_V=1\).

The parameters for the unchanged shapes of the two old six-site forks are \(h_F=29,\ Q_F=320\). The lower readings are \(C_F=111\); \((J_F,R_F,s_F)\) is \((162,158,11)\) or \((129,191,16)\). To compute the left-block coordinates, partition the chain among three ordered subsequences \(D,E,N\): \(D\) precedes \(b\) in its lane, \(E\) follows \(b\), and \(N\) is the other lane. In each subsequence, the chain stays increasing. A legal assignment needs \(D\) empty or ending in 0, and \(E\) empty or starting in 1. Let \(Z\) be the number of such assignments, with subscripts indicating restrictions. Then
\[
 h_U=Z_{N=\varnothing},\quad
 R_U=Z+Z_{D\ne\varnothing},\quad
 C_U=Z_{E=\varnothing}+Z_{D\ne\varnothing,N=\varnothing}.
\]
If \(b\) is in the left lane, the right lane is always eligible, and the left meets the outer apex only when \(E\) is empty. If \(b\) is in the right lane, that lane must start before \(b\), and the left lane cannot end at \(b\), so must be empty. These two cases prove both restricted sums. Streaming through the word only needs 18 states: last contact of \(D\), first contact of \(E\) (each empty/good/bad), and whether \(N\) is empty.

Write \(h=h_U,R=R_U,C=C_U\). The lower endpoint for transitive tail parameter \(2u=\operatorname{bin}(\overline{\operatorname{rev}(t)})\) is
\[
 H(w,u)=D_0+L u-58h q(u),\quad
 D_0=h(111-s_F)+29(C-1)+(R-h)J_F+s_F,\quad
 L=2[320(R-h)+hR_F].
\]
This is the [ordinal lane transfer](VARIABLE-TAIL-SELECTOR.md) after aggregating the left block. Two explicit twin-prime pairs, each checked by **independent full subset DP** on the actual tournaments:
 
| Strong word | Fork | \(u\) | Vertices | Literal pair |
|---|---|---:|---:|---:|
| 11110 | A | 0 | 14 | \(31151\to31153\) |
| 10100 | B | 1 | 15 | \(176459\to176461\) |

Both ends of both pairs are absent as either endpoint of either earlier transitive-variable-tail family. This last check exhausts the old chart coordinates using \(q(v)\le v/2\), so it proves strict enlargement of the realized height set, not just a difference of descriptions.

## 3. Complete exact searches in the new family at four exponents

For word length \(\ell\), fix the last (zero) chain site in \(D\), put \(E=\varnothing\), and assign all earlier sites freely between \(D\) and \(N\), except that \(N\) must be nonempty. This gives \(R-h\ge Z-h\ge2^{\ell-1}-1\). Therefore the test \(D_0\ge129(2^{\ell-1}-1)\) bounds the left length for a prescribed target, even with an unrestricted right tail (\(H\ge D_0\)). For each strong word, the right tail is an exact one-coordinate inverse. With \(T=H-D_0\), \(K=29h\), test \(u\in[\lceil T/L\rceil,\lfloor T/(L-K)\rfloor]\); this is complete since \(0\le q\le u/2\). In suffix buckets of dyadic size \(B\), the [reverse-bit formula](TWO-COORDINATE-RIGIDITY.md#3-a-sublinear-exact-integer-inverse-for-that-short-guard-family) makes the height strictly decreasing in suffix reversal whenever \(K(\lfloor\mathrm{lo}/B\rfloor-1)>(L+K)(B-1)\). Binary-search each bucket, then verify exact equality and interval membership. A singleton is always available.

For \(p=7,13,17,19\), the maximum word lengths were 7,18,26,30 respectively. Both forks and all three lower queries \(E_p-3,E_p-1,E_p+1\) were searched. None hits, including the \(536{,}870{,}911\) strong words at exponent 19. This is an exclusion only for the new family; it does not rule out those targets as tournament counts. The exponent-61 full problem is unchanged.

The transfer lesson is to count with boundary data before discarding it. In Moore, a closed path forgets local scaling exactly, but the weighted sum of all possible arcs records large overlapping cuts. In Chan–Pak, collapsing a block to just its height would lose which lane meets the outer apex. Retaining that small state space makes the new exact counts and a complete finite-height search possible.

Reproduce with:

    python3 weighted-cut-upper.py
    python3 nested-apex-controls.py
    c++ -O3 -std=c++17 nested-apex-targets.cpp -o nested-apex-targets
    ./nested-apex-targets 7 13 17 19

The exponent-19 run examines over half a billion words; use the smaller exponents for a quick replay.
