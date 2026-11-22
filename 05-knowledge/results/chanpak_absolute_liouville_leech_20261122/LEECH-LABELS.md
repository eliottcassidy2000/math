# Completing the eight pure-A height labels, and measuring the perfect-flank defect

**PROVED** physical count constructions and elementary arithmetic bounds. Chan–Pak remains open. No map from a Niemeier lattice or the Leech lattice to tournament isomorphism is asserted.

## 1. The two omitted Coxeter labels are now physical

The [earlier Leech recovery](../chanpak_mrb_leech_clock_20261114/MRB_LEECH_SOURCES.md#4-what-the-six-integers-mean-on-the-leech-side) identifies the eight pure-\(A\) Niemeier Coxeter numbers \(2,3,4,5,7,9,13,25\). Indeed the root system has the form \(A_{a-1}^{24/(a-1)}\), so \(a-1\mid24\). The historical six-label traversal did not include \(9,25\). The [six signed-height pairs](../chanpak_six_signed_height_pairs_20261118/README.md) already build both \(p=a2^k+1\) and \(7p\) for the other six labels, except that no lower witness exists for \(p=7\) or \(21\).

**New result.** For \(a=9,25\) and every \(k\ge1\), both \(a2^k+1\) and \(7(a2^k+1)\) are ordinary counts of **strong** tournaments. In particular both rails are now attained for the full set of eight pure-\(A\) Coxeter labels, apart from the two inherited lower holes. The labels specify the finite catalogue; the actual construction is in marked tournament profiles, not in lattice geometry.

Use the inherited [one-return scaffold](../chanpak_prime_frontier_selector_20261118/PRIME_RAIL.md): a fixed core \(D\), contact set \(I\), new vertex \(b\), and arbitrary nonempty incoming module \(X\). Orient \(b\to X\to D\), \(I\to b\), and \(b\to D\setminus I\). Let \(d=H(D)\), let \(G\) count ordered spanning path pairs \(P,S\) in \(D\) with \(P\) nonempty and ending in \(I\), and let \(D_{II}\) count the same pairs in which \(S\) is nonempty and starts in \(I\). The inherited run formula gives the output as
\[
G f_X(2)+\beta H(X),\qquad \beta=d-D_{II},
\]
where \(f_X(2)\) counts ordered two-slot path covers, allowing either slot to be empty. For \(X=\mathrm{TT}_w\) this is \(G2^w+\beta\).

Outgoing masks below are in vertex order, with vertex 0 at the low bit. The contact column is a bitmask in the same convention.

| Label | Rail | Masks of \(D\) | Contact | \((d,G,D_{II})\) |
|---|---|---|---:|---|
| 9 | lower | 136,17,3,22,101,15,47,126 | 81 | (277,1152,276) |
| 9 | upper | 4,1,10,67,15,31,55 | 96 | (29,63,22) |
| 25 | lower | 2,48,3,7,13,29 | 12 | (23,100,22) |
| 25 | upper | 32,1,3,7,15,158,63,95 | 144 | (51,175,44) |

All four cores have the property that every vertex reaches \(I\). This is an exact finite reachability test and ensures the full host is strong for every nonempty \(X\): \(b\) reaches every \(X\) vertex, \(X\) reaches every \(D\) vertex, and all \(D\) vertices reach \(b\) through \(I\). The upper rows give every \(k\ge1\) directly. The lower 9 row gives \(w=k-7\) for \(k\ge8\); the lower 25 row gives \(w=k-2\) for \(k\ge3\). Literal strong witnesses for the remaining finite lower starts are:

| \(H\) | Outgoing masks |
|---:|---|
| 19 | 2,40,3,5,15,29 |
| 37 | 32,17,3,7,13,30 |
| 51 | 64,9,3,5,15,31,62 |
| 73 | 4,65,34,7,15,27,61 |
| 101 | 2,112,67,7,13,29,57 |
| 145 | 66,48,67,7,13,29,58 |
| 289 | 192,5,49,7,11,27,62,126 |
| 577 | 130,196,97,135,15,27,57,116 |
| 1153 | 144,321,131,7,14,31,61,122,253 |

The same pairs work on a shared **arbitrary** nonempty module \(X\), not just transitive widths. In the right-hand scaffold replace \(X\) by the ordinal join of \(X\) with 7 singleton pads for label 9, or 2 singleton pads for label 25. Ordinal pads preserve \(H(X)\) and multiply \(f_X(2)\) by \(2\) each: assign the new singleton to either slot, where its order is forced. The two hosts therefore attain \(p\) and \(7p\) exactly for
\[
p=1152 f_X(2)+H(X)\quad\text{or}\quad p=100 f_X(2)+H(X),
\]
respectively. The replay verifies the four new profile rows and whole hosts at small widths and on a cyclic input, and verifies all finite starts by subset DP. It does not rely on any completeness assertion for the scaffold search.

## 2. Rough physical pairs give almost-neutral, not prime, aliquot flanks

The [inherited perfect-core translation](../perfect_core_twin_edges_20261019/PERFECT_CORE_TWIN_EDGES.md) concerns prime flanks: if \(C\) is perfect and \(p\nmid C\) is prime, \(s(Cp)=C(p+2)\), where \(s=\sigma-\mathrm{id}\) is the proper-divisor sum. For any coprime integer flank \(N>1\), multiplicativity gives the exact defect
\[
s(CN)-C(N+2)=2C\bigl(s(N)-1\bigr).
\]
It vanishes exactly when \(N\) is prime. Otherwise \(N/p\ge\sqrt N\) for the least prime factor \(p\), so this defect is at least \(2C\sqrt N\). Thus relative neutrality does not imply the literal affine transition.

For \(z\)-rough \(N\), with \(z\ge2\), we have
\[
1\le\frac{\sigma(N)}N\le\prod_{p\mid N}(1-1/p)^{-1}
  \le (1-1/z)^{-\omega(N)},\qquad
\omega(N)\le\frac{\log N}{\log z}.
\]
With the simultaneous roughness cutoff \(z=\exp(c\sqrt{\log N})\) this yields
\[
0\le s(CN)-C(N+2)
 \le O_C\!\left(N\sqrt{\log N}\,\exp(-c\sqrt{\log N})\right).
\]
Indeed \(-\log(1-1/z)=O(1/z)\) and the total exponent tends to zero. The same bound applies at \(N+2\) after adjusting constants, and for all sufficiently large heights both flanks are automatically coprime to a fixed \(C\). The proved rough pair family therefore supplies relatively near-neutral **matched** aliquot flanks at unbounded heights. At either endpoint a literal perfect-core \(+2\) transition still demands exact primality, separated by the absolute square-root defect above. No twin-prime or odd-perfect-number inference follows.

There is a literal tournament interpretation after removing the unavoidable dyadic factor. For an even perfect core \(C=2^{r-1}M\), \(M=2^r-1\) prime, \(C\) itself cannot be an ordinary Hamilton count (Rédei parity). The odd \(M\) is realized for \(r=2\) by the directed triangle and for \(r\ge5\) by the inherited signed dyadic rail. Ordinally adjoining one such \(M\)-witness to each member of a strong \(N,N+2\) pair realizes the normalized endpoints \(MN,M(N+2)\) with gap \(2M\). These larger hosts need not be strong. For \(r=3\), the odd factor \(M=7\) is a genuine obstruction to this **factorwise** construction, not a proof that its normalized endpoints never occur by different cores. Reintroducing the arithmetic factor \(2^{r-1}\) does not make the perfect-core step exact except at prime \(N\).
