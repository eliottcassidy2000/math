# Two local apices turn all 4096 Golay words into different path heights

**PROVED with a bounded exact certificate.** Replace each of the two marked chain sites by 12 controls, a tail of length \(t\ge1\), and a local apex. Flip the 24 control–apex contacts according to the extended binary Golay code. At **every** fixed \(t\ge1\) the 4096 heights are distinct. The zero word is exactly the shifted-square rail at \(k=t+13\); a nonzero word gives the prime \(139472229779\). At fixed tail the minimum gap in this entire chart is at least **1,023,132**, with equality at \(t=1\). Thus this construction distinguishes codewords but does not give a unit-gap inverse for Chan–Pak.

## One-module reads

Use the [two-slot quotient](UNEQUAL-CHAINS.md) in the order/transposition that gives
\[
 H=4Q_AQ_B+Q_A D_B+3D_AQ_B-D_A D_B.
\]
Explicitly, use the converse copy with outgoing rows \((10,4,9,18,39,15)\), marks \(A=5,B=1\). The old-to-new label map from the unequal-chain note, with all arcs reversed, is \((3,2,1,0,5,4)\); the two roles are exchanged. The replay computes the corresponding four run counts independently.

For each module use an ascending spine, first twelve controls, then \(t\) tail vertices. Every spine vertex is initially directed to its local apex. An O contact reverses the corresponding control–apex arc; all other vertices are I contacts. There are twelve independent physical arcs in each module. With the all-zero word each module is the original ordered chain of length \(t+13\). All expanded tournaments are strong, since the six-site quotient was strong and outside walks connect every module in both directions.

Index the controls \(i=0,\ldots,11\). For one module put \(B=\sum_{i:O}2^{11-i}\), and let \(E\) sum \(2^{j-i-1}+[j=i+1]\) over O–O pairs \(i<j\). Then \(D=1-E+B2^t\). This follows from the [one-apex word count](../chanpak_tiling_viewpoint_20261123/precedence-binary-chart.md): the transitive spine is split into the two ascending pieces around the apex, and an O–I pair with \(d=j-i-1\) intermediate sites contributes \(2^d+[d=0]\). For fixed O at control \(i\), summing over all later sites tentatively treated as I gives \(2^{12+t-1-i}\); subtract E to remove the actual O–O pairs.

Let \(n=12+t\). The second read is
\[
 Q=2^{n+1}+\sum_{i<j:\, O_i,I_j}2^{n-1-d}(3^d+1),\qquad d=j-i-1.
\]
Choose the slot containing the apex and the subset of the spine assigned to that slot. The complement is transitive. For an O–I pair, summing over the \(d\) optional intermediate sites gives \(\sum_u\binom du(2^u+[u=0])=3^d+1\). All outside sites are free and the apex has two slots, giving the multiplier. For the tail sum this becomes \(Q=\alpha 3^t+\beta 2^t-\gamma\), where
\[
 \alpha=\sum_{i:O}2^{i+1}3^{11-i},\quad
 \gamma=\sum_{i:O}2^{i+1},\quad
 \beta=8192+\gamma-\alpha+
 \sum_{i<j\le11:\,O_i,I_j}2^{11-d}(3^d+1).
\]
This gives the exact integer height for unequal tails, too; the all-width distinction below is for equal tails.

## Separation, a nonzero prime, and an exact distance scale

Use \(1+z^2+z^4+z^5+z^6+z^{10}+z^{11}\) and its twelve shifts, with parity coordinate appended. Assign the first twelve coordinates to A and last twelve to B. The first twelve are an information set (the projection of the basis is triangular). So are the last twelve: a systematic matrix \([I\mid L]\) for the self-dual code has \(LL^{\mathsf T}=I\). Every nonzero word therefore activates both blocks; its count has leading term \(4\alpha_A\alpha_B\,9^t>0\), while the zero word grows on the \(4^t\) scale.

At equal widths the count is a linear combination of the six modes \(1,2^t,3^t,4^t,6^t,9^t\). The replay computes distinct six-mode vectors for all 4096 codewords. The \(9^t\) coefficient alone has only three collisions, all pairs, and within each the \(6^t\) coefficients differ. Thus any pair difference has highest mode with base at least 6 and nonzero integer coefficient. The sum of all coordinate spans is \(W=17696218772230\). Every lower-mode ratio is at most \(3/4\); \(W(3/4)^{107}<1\), so the leading coefficient dominates at every \(t\ge107\). For \(t\ge108\), the difference has magnitude at least \(6^t/4\). The finite scan at widths 1 through 107 establishes the minimum gap 1,023,132, at width 1 between masks 0x18cf4b and 0x1e49d9 (heights 3147435342621 and 3147436365753). These statements concern equal-width comparisons, not different widths.

At \(t=1\) mask 0x1007b4 gives the prime \(139472229779\) by exact trial division; zero also has prime height 1073807359 at that width. This contrasts with the [triple/dodecad read](DODECAD-CHAIN.md), in which every nonzero word forces a factor of three. Neither result supplies infinitely many primes in its chart or makes an arbitrary prescribed prime attainable.

## First statistical difference occurs at variance

Each of \(D,Q\) is multilinear of degree at most two in its own twelve switches (also clear from “O at i, I at j”). Hence \(H\) has degree at most four, and its expectation over the Golay code is exactly the full-cube mean by strength seven.

For each pair of contacts in one half, the quadratic coefficients of \(D,Q\) are \(-a,-b\) with \(b\ge a>0\), by the displayed formulas; at separation \(d+1\), compare \(2^{n-1-d}(3^d+1)\) with \(2^d+[d=0]\). Each quartic coefficient of \(H\) supported on two switches per half is therefore positive (\(-a_Aa_B+b_Aa_B+3a_Ab_B+4b_Ab_B>0\)); no other shape supports a quartic term. Write \(h_I\) for this coefficient at a four-subset \(I\). On reducing \(H^2\) to squarefree monomials, only eight-subsets with four coordinates per half can contribute to the expectation difference. Such a projection has the full law unless the subset is an octad; for an octad the all-ones probability is larger by \(2^{-8}\). Thus
\[
 \operatorname{Var}_{\rm Golay}(H)-\operatorname{Var}_{\rm cube}(H)
 =2^{-8}\sum_{\substack{O\text{ octad}\\ |O\cap A|=4}}
   \ \sum_{\substack{I\subset O\\ |I\cap A|=|I\cap B|=2}}
       h_I h_{O\setminus I}>0.
\]
There are 255 such octads in this placement (checked). This is a genuinely two-module effect: [in the earlier linear cycle register](../../../session212/cycle-spectrum/PORTS-AND-GOLAY.md) the first pure moment gap is at degree eight in the *counts*, whereas here it is already visible in the variance. Since all path heights now distinguish the codewords at each width, no nonidentity affine Mathieu action on the marked cube preserves this count readout, although all preserve the labelled edit metric.

Replay: run apex-contacts-check.py with Python 3.
