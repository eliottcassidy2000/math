# An exact binary strip of local seven-counts inside a ternary bulk

For any digit $d\in\{0,1,2,3\}$ and $k\ge1$ consider the strong two-hub word $13\,d^k\,00$. Match the first two spine sites to the last two, independently of where the intervening $k$ sites go. Both pairings make two disjoint increasing components with endpoint types $(1,0)$ and $(3,0)$. Each intervening site can be assigned to either component without moving its endpoints. Hence there are exactly $2^{k+1}$ **selected bare forests** (two pairings times $2^k$ assignments); every one has seven gluings. They contribute precisely $14\,2^k$ directed Hamiltonian paths to the ordinary count. This is an actual doubling strip, not 2-adic divisibility.

Yet the ordinary count on this same tournament is

| central digit $d$ | exact ordinary count |
|---|---|
| 0 | $-5+24\,2^k+18\,3^k$ |
| 1 | $-9+25\,2^k+20\,3^k$ |
| 2 | $-9+5\,2^k+40\,3^k$ |
| 3 | $-11+24\,2^k+24\,3^k$ |

Proof of the all-size formula: fix the two-hub order and occupied run set of the homogeneous central block; its interior assignments then depend only on how many of the three slots it occupies. They are surjections counted by $1$, $2^k-2$, or $3^k-3\cdot2^k+3$. Thus all ordinary counts lie in the span of $1,2^k,3^k$. The first three positive lengths give an invertible basis and are independently checked by literal permutation enumeration in the replay below. Strongness holds because the hubs and the first site form a directed triangle; the second site is reached from the hubs and can reach them via a last zero. Every middle site can be reached through an earlier spine site and reaches the hubs via a final zero.

Consequently the fraction of paths in this selected family decreases on the scale $(2/3)^k$ despite doubling their number at every step. This makes the 2/3 competition visible in one literal tournament family. No one of the four central letters removes the ternary bulk. A constructive inverse must change more than the multiplicity of the seven-gluing witnesses. In particular multiplying the local seven by a smooth weight is not a way to prescribe an ordinary count of $7\cdot2^k$.

## Replay

Run **python replay_seven_sheets.py**. With k=1,2,3 the largest check has nine vertices. The replay also reconstructs each selected bare forest directly from all the literal paths, for a check independent of the permanent formula.
