# An exact crossing at exponent 18 and two different pairs at exponent 21

**FINITE-EXACT.** There are two strong tournaments on 32 vertices differing in a single arc with Hamilton path counts **34,359,607,295** and **34,359,607,297**. Their midpoint is \(E_{18}=2^{17}(2^{18}-1)\). At \(E_{21}\) there are at least two different strong lower isomorphism classes on 40 vertices with a single-arc exit to the upper count. The two upper classes also differ. Both exponents are composite. At exponent 61 the corresponding existence question remains open.

Use a core \(D\) with outgoing row masks in vertex order starting at zero. The apex beats precisely the vertices in the tag mask. The left transitive flank takes the bits of \(x\) from high to low and the right flank the complemented bits of \(y\) from low to high. A one means the apex beats the flank vertex. The left flank, core and right flank form an ordinal sum after deleting the apex. Each pair uses identical flanks and contacts and reverses only the indicated core arc.

| Exponent | Lower core \(D\) | Arc | Tags | \(x\) | \(y\) |
|---:|---|---|---:|---:|---:|
| 18 | 80, 21, 73, 131, 44, 15, 58, 119 | (4,6) | 89 | 7 | 749445 |
| 21 | 4, 17, 34, 7, 13, 27 | (2,3) | 22 | 7 | 923014801 |

The exact [two-flank transfer](../chanpak_gap_star_ports_20261122/README.md) uses the core path count \(h\), ordered two-slot count \(q\), and even/odd masses for tagged starts \((s_0,s_1)\) and untagged ends \((e_0,e_1)\). With \(k\) the apex-only count, its linear and constant terms are
\[
 (e_0+e_1)x+(s_0+s_1)y-(e_0-e_1)[x\ {\rm odd}]-(s_0-s_1)[y\ {\rm odd}]+k .
\]
The other terms \(h(f(x)+f(y))+qxy\) agree across each of these flips.

At exponent 18 the two profiles share \(h=413,q=6200,(e_0,e_1)=(1266,1527)\). The start pair changes from \((1042,1253)\) to \((1043,1252)\) and \(k\) rises from 1575 to 1579. Thus the gap is \(4-2[y\ {\rm odd}]\) for **every** choice of both flanks. The indicated addresses give the stated absolute heights. At exponent 21 the gap is \(-2+2[x\ {\rm odd}]+2[y\ {\rm odd}]\). The indicated pair has counts 2,199,022,206,975 and 2,199,022,206,977.

Each endpoint of the new exponent-21 pair has exactly one vertex of outdegree one. Each endpoint of the earlier 40-vertex pair has two. Thus the count pair and even its common order forget which marked tournament gave the crossing. This gives a lower bound of two on the number of marked edges at order 40. The first possible order and the total number of edges at order 40 are separate questions. [Earlier exponent-21 certificate](../chanpak_middle_out_square_sum_20261122/EXACT-SEAMS.md)

Run the linked Python script from this directory. It uses subset path DP on the two small cores to derive both profiles, evaluates the exact transfer, and constructs all full orientations to check strong connectivity, the unique edit and the degree distinction. [Certificate replay](other_seams.py)
