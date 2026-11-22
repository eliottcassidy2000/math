# Every within-port and between-port edit misses the strong lower rail

**PROVED.** Let \(T_p\) be the three-port strong realization of \(E_p-1\) for odd \(p\ge5\), with ports \((6,4,1)\) of widths \((p,m,m)\), \(m=(p-3)/2\). No single reversal with both endpoints in the union of these ports has count \(E_p+1\). Each reversal strictly increases the Hamilton count except for reversals between consecutive vertices in the same port, which preserve the isomorphism class. The three exact minimum increases across distinct ports are below. [Source, two-return bottleneck and 45 one-run skeletons](../chanpak_support_birth_grammar_20261122/ODD-PERFECT-PORT.md)

## Within a port

Write its width as \(w\) and the distance between reversed vertices in transitive order as \(d\). For \(d=1\), exchange the labels to exhibit the isomorphism. For \(d\ge2\) every new odd cycle in the block must use the reversed arc. Its remaining vertices are an odd, hence nonempty, subset of the \(d-1\) intermediate vertices in increasing order. The cycles pairwise intersect, so the one-run Hamilton count increases by \(2^{d-1}\). The two-run count is a sum over ordered nonempty bipartitions of products of induced Hamilton counts, each at least one. Its increase is nonnegative. In fact its exact increase is
\[
 4\sum_{\substack{1\le r<d\\r\ \mathrm{odd}}}{d-1\choose r}(2^{w-r-2}-1).
\]
Each new cycle on \(r+2\) vertices must lie on one of the two sides, every remaining vertex can choose its side except for the assignment leaving the other side empty, and its packing weight is two. All Hamilton paths in the inflated host still use at most two runs per port, because the bound on returns uses only edges exterior to the port. In particular, changing its orientation preserves the allowed external skeletons and can only increase their filling counts. The 45 skeletons using all three ports once each already give a gain of at least \(45\cdot2^{d-1}\ge90\). [THM-002, odd-cycle packing formula](../../../01-canon/theorems/THM-002-ocf.md)

## The minimum crossing a transitive pair is a corner

The following monotonicity holds for **any** tournament constructed by uniform substitution of transitive blocks. Fix two distinct blocks \(A\to B\). Number their vertices from zero in transitive order. Reversing \(A_i\to B_j\) gains all paths using the exceptional \(B_j\to A_i\) and loses all using the original \(A_i\to B_j\); every other path survives. In each summand of the compressed path grammar, distinguish precisely that adjacency, and let \(r\) be the number of runs in a specified transitive block.

The number of fillings where a specified run ends at index \(i\) is
\[
 D(w,i,r)=\sum_{k=0}^{r-1}(-1)^{r-1-k}{r-1\choose k}(k+1)^i k^{w-1-i},
\]
where \(0^0=1\). Indeed the earlier labels may use the designated run, the later ones may only use other runs, and every other run has to be nonempty. For a specified starting vertex the weight is \(D(w,w-1-i,r)\) by label reflection. An unmarked port of multiplicity \(r\) contributes \(P(w,r)=\sum_{k=0}^r(-1)^{r-k}{r\choose k}k^w\), the number of surjections of its vertices onto its runs.

The endpoint weight \(D(w,i,r)\) is nondecreasing in \(i\). In a filling with designated endpoint \(i\), swap the run assignments of labels \(i\) and \(i+1\), then sort each run. The latter was in another run, and it becomes the new endpoint. This is injective. All skeleton coefficients are nonnegative. Therefore the count lost from the old orientation increases weakly with \(i\) and decreases weakly with \(j\), and the count gained in the new orientation has the opposite monotonicities. Subtracting proves that **the smallest signed change occurs at \((i,j)=(w_A-1,0)\)**. This also applies when one or both widths equal one.

## Exact three-port corner values

Set \(x=2^m, y=3^m\). For \(m\ge2\), the minimum increases in \(T_p\) are

| Source arc | Least change |
|---|---:|
| \(6\to4\) | \(12x^3-2x^2+9x/2-2\) |
| \(6\to1\) | \(40x^4-4x^3+14x^2-x/2-6+y(36x^3-40x^2+9x-10)\) |
| \(4\to1\) | \(32x^4+4x^2\) |

Every entry exceeds two for \(x\ge4\). At \(m=1\) the three values are respectively 94, 860, 426. At exponent 61 the respective minima are 1,856,910,058,351,609,662,461,181,950, 382,324,636,509,376,268,260,372,254,741,437,563,607,100 and 2,658,455,991,569,831,746,960,535,625,167,536,128. These cover all \(61\cdot29+61\cdot29+29^2=4379\) labelled cross-port flips.

Here is the finite derivation of the table. For ordered ports \(a,b\), form the tensor \(C_{ab}(r_6,r_4,r_1)\) of spanning path skeletons with exactly one distinguished \(a\to b\) transition. Every other transition has to follow the source. A mark may occur at most three times. Indeed in a fixed block, the subwalks exiting and returning between consecutive visits are arc-disjoint. At most one uses the distinguished arc, including either exit or entrance. Among the remaining subwalks the original singleton bottleneck allows at most one. These bounds apply also when the exceptional arc has a direction forbidden by the source. The ordinary vertices may occur only once. To compute \(C\), append each possible next vertex to the current end, keeping the ordinary-vertex mask, the three counts, and whether the distinguished step has been used. A step of the old orientation may be read as ordinary or as distinguished; if there are several such occurrences the specified endpoints can belong to at most one of them in a filling. This counts exactly the desired paths.

At the minimizing corner the old endpoint weights are \(P(w,r)/r\), because specifying the run of the largest or smallest vertex makes it the required endpoint automatically. Both new endpoint weights instead become \(P(w-1,r-1)\). Multiply the corresponding tensor entries by these two weights and by \(P(w_c,r_c)\) for the remaining port. Subtract the old total from the new. For \(r=1,2,3\) the three unconstrained factors are \(1,2^w-2,3^w-3\cdot2^w+3\). With \(w\ge2\), the new corner factors are \(0,1,2^{w-1}-2\). Substituting \(2^p=8x^2\) and \(3^p=27y^2\) gives the displayed polynomials. The replay certifies the entire symbolic identity by comparing values at an independent five-by-five grid in \((x,y)\): the degree in each variable is at most four because the three widths have slopes 2,1,1 in \(m\). It also checks the actual \(p=5\) changes against path counts of literally flipped tournaments.

Run `python3 -B 05-knowledge/results/chanpak_fermat_solids_seams_20261122/port_flips.py` from the repository root. This result concerns the specified source and reversals with both endpoints in its three large ports.
