# A prime-bearing exact ×7 rail with a deletion tradeoff

A background lane found this complementary construction after the [Paley rail](../chanpak_fano_square_compression_20261118/SEVEN_SHEETS.md), which had true deletion compatibility but composite lower counts only. This new rail has genuine free rotations and known prime lower counts, but no general prime or deletion-compatibility theorem.

Let b beat all of A, A beat all of D, and let I⊂D beat b while b beats O=D−I. Write d=H(D), G for ordered splits of D as P,S with P a nonempty directed path ending in I and S an optional directed path, and Q for the same splits with S empty or starting in O. Split a spanning path according to the number of A runs. If there are two, it starts with an A run and the second starts immediately after b, giving G times the ordered nonempty two-cover count of A. If there is one run starting the path, its continuation has the form P,b,S and contributes Q. If that single run is not first, it follows b, giving d+G. Thus for any nonempty A the count is

    (d+G+Q) H(A) + G Q_2(A).

For the lower family take D outgoing masks (42,116,121,114,97,64,1) and I={0,5}, giving (43,284,266). For the upper take circular R on Z/7 with steps 1,2,3, with an ordinal sink c appended inside D, and I={c}. Then the profile is (175,1988,1988): the first split ends at c, and the two optional pieces of R give 1638+2·175=1988. The pair of coefficients is (593,284) and seven times that pair. Both hosts are strong. The upper path set admits actual free rotation of R on the vertex labels; it has seven sheets preserving the allocation of A to nonempty slots. For transitive width a≥1 this gives exactly

    H_lower = 284·2^a+25,  H_upper = 7 H_lower.

Cases include 593→4151 (prime), 1161→8127, and 2297→16079 (prime). This makes no assertion that there are infinitely many primes of the exponential form.

The two-run deletion indicator separates the families. Removing the first A run always gives a valid insertion; removing the second A run is valid iff S is empty or starts in O. At width two each template has two allocations. The lower one-face/both-face counts are (2(G−Q),2Q)=(36,532), while the upper gives (0,3976). Genuine sevenfold replication of counts and compatibility with actual deletion are different conditions even after admitting prime inputs.

`python prime_rail.py | diff -u prime_rail.out -` verifies all the coefficient data and literal host instances.

## First appearance versus the paired presentation

The complete small-order atlas locates 593 first at **order 8**, in 22 isomorphism classes, all with positive arc loads. For example, use masks (196,193,50,7,11,27,60,124). The paired lower rail uses nine vertices instead, so minimizing the source and obtaining the coupling are different optimization problems. The upper paired witness for 4151 uses **order 10**, and that is minimal: the complete order-nine atlas has maximum 3357. It has no unused arcs, and the sevenfold label symmetry compresses its 45 arcs to nine load orbits. With R vertices 0–6, c=7, A={8}, b=9, orbit representatives and loads are:

| Representative arc | Orbit size | Load per arc |
|---|---:|---:|
| 0→1 | 7 | 1050 |
| 0→2 | 7 | 1062 |
| 0→3 | 7 | 978 |
| 0→7 | 7 | 543 |
| 7→9 | 1 | 3976 |
| 8→0 | 7 | 543 |
| 8→7 | 1 | 175 |
| 9→0 | 7 | 259 |
| 9→8 | 1 | 2163 |

The sum of multiplicity times load is 9·4151, as required. The order-nine atlas also supplies a smaller first-born witness for 2297 (masks 320,133,417,71,15,27,54,121,250), whereas this rail's paired witness uses eleven vertices.
