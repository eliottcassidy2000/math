# All-size growth and exact unicyclic range

Let G_n be the simple graph on {1,...,n} with distinct vertices joined if their sum is a square. For every n≥5, the new vertex n has at least one neighbor in {1,...,n−1}. Equivalently, there is a square k² in the interval n<k²≤2n−1.

Proof. For n=5,6,7,8 use k²=9. For n≥9 let r=⌊√n⌋≥3 and k=r+1. Then k²>n, and

(r+1)² ≤ 2r²−1 ≤ 2n−1,

since r²−2r−2≥0 for r≥3. Thus v=k²−n is an older positive vertex joined to n. This is an all-n proof; the finite checks below are supporting tests, not a substitute for it.

Consequently, if G_n is connected for some n≥5, every later prefix graph is connected. Indeed, G14 is connected: G12 has three path components; 13 joins the first and second and 14 joins the second and third. The same edge enumeration gives 13 edges at 14 vertices, so G14 is a tree.

For a connected graph its cycle rank is μ=|E|−|V|+1. Adding a new vertex with d older neighbors changes μ by d−1≥0. Exact older-neighbor sets are {1,10}, {9}, {8}, {7}, {6,17} at n=15,16,17,18,19, respectively. Hence μ=1 throughout G15,...,G18, while μ(G19)=2 and μ(G_n)≥2 for all n≥19. Therefore the connected unicyclic graphs in the entire prefix family are exactly G15,G16,G17,G18: G1 is an isolated vertex, G2,...,G12 have 1 and 2 in different G12 components, G13 still has the component containing 2 separate, and G14 is a tree.

This limits application of a connected-unicyclic Hamilton-path theorem to precisely those four prefix orders, not the later graphs.

`check_growth.py` independently enumerates older square-sum neighbors, checks components of G2,...,G14, verifies the old-neighbor and cycle-rank statements through n=10000, and checks the square-window witness through n=10000. It uses explicit checks that remain active under `python -O`; its standard and optimized outputs agree. The unbounded conclusion follows from the displayed inequalities and induction, rather than from the finite range.
