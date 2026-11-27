# The closing tile as a boundary charge on polygon dissections

This links the [closed/open tiling split](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md), the [cycle seam and path leakage](../../../session212/cycle-spectrum/PORTS-AND-GOLAY.md), and the [dihedral polygon path action](../fermat_dihedral_tournaments_and_lonely_chord_grids_20261127/README.md). The polygon gives a local gluing law for the marked path's closure, and an explicit reason that gluing alone does not recover *how many other paths or cycles* exist.

## 1. Add the right triangle statistic and all shared diagonals cancel

Put a marked directed Hamiltonian path \(0\to1\to\cdots\to n-1\) around a convex polygon. Its closing contact \(\{0,n-1\}\) is one free tile, with bit \(e=1\) when \(n-1\to0\) and \(e=0\) in the other direction. For an arc read in a given boundary direction put \(s_{ij}=1\) if \(i\to j\), and \(-1\) otherwise. On a triangle \(i<j<k\) use its cyclic polygon order, and define the boundary charge
\[
q(i,j,k)=\tfrac12(s_{ij}+s_{jk}+s_{ki}-1).
\]
The four possibilities are \(1\) for a cycle following the polygon order, \(-2\) for a cycle opposed to it, \(0\) for a transitive triangle with circulation \(+1\), and \(-1\) for a transitive triangle with circulation \(-1\).

For **every** triangulation \(\mathcal D\), the exact identity is
\[
\sum_{\tau\in\mathcal D}q(\tau)=e .
\]
Indeed each diagonal occurs along opposite directions in its two adjacent triangles, and thus contributes zero. The boundary circulation is \(n-2+2e\); there are \(n-2\) triangles, and the shift by one on each triangle subtracts exactly this baseline. The same proof works for polygonal dissections if the charge of a face of size \(r\) subtracts \(r-2\): these baselines add under gluing. For a local quadrilateral retriangulation it gives the four-term equality \(q(abc)+q(acd)=q(abd)+q(bcd)\), purely from the common boundary. A diagonal flip between *dissections* is different from an orientation flip of a tile in the metagraph. For the latter, an internal diagonal shifts the charge of its two incident faces in opposite directions, while the closing tile changes just the root face and hence the whole sum by one.

Writing \(t_+,t_-\) for the two types of cyclic face and \(u_-\) for the negatively read transitive faces gives \(t_+=e+2t_-+u_-\). A directed boundary therefore forces a positively oriented cyclic triangle in **every** dissection into triangles, even if reversing one local shortcut cannot peel the cycle. For example, a directed pentagon with every step-two shortcut reversed has no removable vertex preserving the directed boundary, but its signed face charges still add to one. A single literal diagonal of a directed boundary closes exactly one of the two child boundary paths, not both: the children need opposite directions of that diagonal to both follow the original polygon orientation.

## 2. Exactly when all the little polygons can be directed triangles

An entirely cyclic triangulation of a directed boundary has \(u_-=0\) and no transitive faces. Adjacent faces must have opposite signs, and every boundary-touching face must be positive. Each negative face is therefore internal with degree 3 in the dual tree. Now \(t_+=1+2t_-\) and \(t_++t_-=n-2\), hence \(n=3(t_-+1)\). Conversely, in a dissection whose two-colored dual has every negative face internal, orient positive faces along the polygon and negative faces oppositely; shared orientations agree, the boundary is directed, and every face is cyclic. Remaining tournament chords are free. The condition concerns existence of a suitable orientation for a given dissection, not existence in every tournament.

The [rooted planar-tree model](https://ics.uci.edu/~eppstein/pubs/Epp-SoCG-07-talk.pdf) counts such dissections on a fixed labelled polygon. Root at the given closing side. A positive face has two ordered child slots, either empty or occupied by a negative face, and every negative child must have two positive children. With \(x\) marking negative faces the count series \(P\) satisfies \(P=(1+xP^2)^2\). Thus \(P=R^2,\ R=1+xR^4\), and Lagrange inversion gives, for \(n=3(t+1)\),
\[
[x^t]P=\frac{2}{4t+2}\binom{4t+2}{t}.
\]
The first counts at \(n=3,6,9,12\) are \(1,2,9,52\), independently confirmed by enumerating the triangulations. With uniformly random independent orientations of the nonboundary contacts after fixing a directed boundary, the expected number of entirely cyclic triangulations is that count divided by \(2^{n-3}\): each of its \(n-3\) internal diagonals has exactly one prescribed orientation. This expectation gives no claim that a specific polygon or tournament has such a triangulation.

## 3. The metagraph reads both the marked path and the global response

Write \(H(T)\) for all directed Hamiltonian paths, \(C(T)\) for directed Hamiltonian cycles up to rotation, and \(a(T)=|\mathrm{Aut}(T)|\). Each cycle has exactly \(n\) path openings, and a closable path records its closing arc. The fixed-spine tiling class therefore has
\[
\frac{nC(T)}{a(T)}\text{ closed presentations},\qquad
\frac{H(T)-nC(T)}{a(T)}\text{ open presentations}.
\]
Here the group action on ordered paths is free. Toggling the closing tile in the marked cube always exchanges the two presentation states and changes their face charge by \(1\) in the opening-to-closing direction. After passing to isomorphism classes it need not increase \(C\): the move is not just adding a cycle while retaining every old one.

For consecutive endpoint order \(s=0,t=n-1\), let \(T_+\) close the spine using \(t\to s\), and \(T_-\) have \(s\to t\). The exact cycle displacement is
\[
C(T_+)-C(T_-)=N_{s,t}-N_{t,s},
\]
where \(N_{u,v}\) counts all directed Hamiltonian paths with those ordered endpoints in the common network with the closing contact removed. Such paths never use the contact internally. Removing the flipped arc from every affected cycle gives the claimed two terminal counts. In contrast the path displacement counts paths using the flipped arc in any position (two different ordered contractions), so its sign need not match the cycle displacement.

Three exact controls from the fixed-spine cube, reporting \((H,C)\) before and after **closing that very same marked path**:

| order | opening | closing | result |
|---:|---:|---:|---|
| 5 | \((13,2)\) | \((15,2)\) | cycle count unchanged |
| 6 | \((37,4)\) | \((37,5)\) | path count unchanged |
| 7 | \((109,10)\) | \((115,9)\) | more paths, fewer cycles |

The marked witnesses are masks 10, 391, and 1095, where bits reverse the nonconsecutive nonclosing pairs in lexicographic order. The replay scans all fixed-spine pairs up to order seven. It finds no path-neutral cycle-changing case at smaller orders than six, and no opposite-sign change at smaller orders than seven, and validates the terminal formula throughout.

## 4. Why overlap data must be kept

All triangle orientations in a fixed labeled triangulation determine its boundary and internal diagonals but leave \(\binom n2-(2n-3)=(n-2)(n-3)/2\) contacts free. The allowed completions are a literal subcube. For the fan at 0 on five vertices, the closed tilings from masks 0 and 4 have **identical oriented triangles on every face**, and in both cases the same local closing change is seen. Yet their global changes are \((\Delta H,\Delta C)=(8,1)\) and \((12,3)\). So gluing the same overlapping triangle subtournaments does recover the marked closure, but it cannot recover the global path or cycle response without data about cross contacts. This is the sharp distinction between the polygon's shared diagonal, which cancels, and the unused diagonal of the complete tournament, which can change the count.

Run `node replay.js` for the fixed-spine scan, 25,426 direct charge comparisons, and the compatible dissection counts.
