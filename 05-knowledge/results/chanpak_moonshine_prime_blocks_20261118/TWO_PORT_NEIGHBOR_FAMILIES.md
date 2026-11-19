# Exact two-parameter families with Hamilton counts at distance two

**All-order construction.** Two independent transitive module lengths
produce several new +2 pairs and, in one family, both neighbors of
every central odd target. These are sparse families, not universal
neighbor coverage or a Chan–Pak proof.

## A two-port extension of the marked reader

For a tournament with distinct marks 0 and 1, substitute nonempty
tournaments A and B at the marks. Every cross-module contact is
uniform, in the direction of the corresponding source arc.
Let c_rs count legal compressed words visiting all the unmarked
vertices once and marks 0 and 1 respectively r and s times.
A visit counts one maximal run through the corresponding module.
Between two visits to the same mark there is at least one unmarked
vertex: the other mark alone cannot supply both directions of a
return. Thus r,s<=n−1 for a source on n vertices. This two-port
bound would not hold automatically with three cyclic marks.

Compression of a Hamilton path is a bijection to a legal word
and ordered nonempty path covers for the two modules. Hence

    H = sum_(r,s) c_rs Q_r(A) Q_s(B),

where Q_r counts ordered r-slot nonempty covers (the inherited
substitution principle, THM-3121). For A=TT_m, B=TT_k this is
sum c_rs · r!S(m,r) · s!S(k,s). In particular Q_1=1,
Q_2(TT_m)=2^m−2, Q_3(TT_m)=3^m−3·2^m+3.
If source signatures agree except in c_11, switching sources
translates the result by exactly the c_11 difference for *all*
independent lengths m,k. For arbitrary nontransitive modules
it instead scales the offset by H(A)H(B).

## Literal witnesses and exact formulas

Use outgoing masks and mark labels 0 and 1 in every row.
All omitted c_rs are zero. A bit at i in a row represents an
arc from that row's vertex to i.

**Family I (seven-vertex sources):**

    [110,12,88,0,11,30,58]
    [12,61,80,4,9,29,59]

Common coefficients c_12,c_21,c_22 = 14,18,8;
c_11 runs through 31,33. For every m,k>=1 these realize

    H = 8·2^(m+k) + 2^(m+1) − 2^(k+1) ± 1.

Each member is an odd target with the other member at distance 2.
Unlike the earlier doubled reader in which the count depended
only on m+k, fixing m+k here still yields m+k−1 distinct centers:
2^(m+1)−2^(k+1) is strictly increasing in m on that diagonal.
For sum s the centers lie between 7·2^s and 9·2^s,
disjoint across s. This certifies order (log X)^2 distinct
centers up to X, both upper and lower bounds.

**Family II (seven-vertex sources):**

    [18,24,91,1,8,31,59]
    [94,124,72,0,12,29,56]

Common coefficients c_12,c_21,c_22 = 6,9,6;
c_11=15,17. The counts are

    H = 3(2^m−1)(2^(k+1)−1)+6  or  +8.

The lower count here is always a multiple of 3, the upper is not.

**Family III (eight-vertex sources):**

    [80,121,227,5,12,25,56,123]
    [96,121,203,33,13,20,56,123]
    [240,249,203,161,140,20,184,32]

Common coefficients c_12=88, c_21=179, c_22=88,
c_31=48, c_32=24; c_11=177,179,181. For every m,k>=1 the three
counts are N−2, N, N+2 with

    N = 16·2^(m+k) + 3·2^m + 24·2^k·3^m − 16·2^k − 3.

Thus **both neighbors** of every designated central odd N
are realized, and so is N itself. For instance m=1,k=2 gives
353,355,357; at m=k=1 the counts are 177,179,181.
For 1<=m<=k, v_2(N+3)=m: the term 3·2^m is the only
term of that valuation. With m fixed the formula is strictly
increasing in k. This proves genuinely separate count parameters
and a lower bound of order (log X)^2 distinct centers up to X
by keeping both lengths below a sufficiently small multiple
of log X.

**Family IV (an actual single-arc reversal at every height):**

Start with [54,100,64,23,70,92,9]. Mark 1 and 2 (not 0 and 1
for this row). Flip the arc between unmarked vertices 3 and 4
to obtain [54,100,64,7,78,92,9]. The higher coefficients are
c_12=44, c_21=89, c_22=44, c_31=24, c_32=12, and c_11=89,91.
This yields the counts

    H = 2^m + 2^k(8·2^m + 12·3^m − 8) ± 1

at all m,k>=1. The same arc reversal still relates the expanded
witnesses because neither changed endpoint was inflated. The
coefficient of 2^k has 2-adic valuation exactly 2, so the even
center recovers m by its valuation whenever m<k+2. Then it
recovers k by the exact formula. Thus this gives another
order-(log X)^2 family, now with a literal single-arc +2 move
throughout.

## What the verification establishes

Run python3 two_port_neighbors_replay.py. It computes the full
finite marked-word signature (recursive state = used outside
vertices, last symbol, two visit tallies) and independently
counts literal inflated tournaments on a 3-by-3 grid by subset DP.
The coefficient identity proves all heights; the grid is only
a test of the implementation. Eight-vertex witnesses were found
by trying all new-vertex contacts above selected seven-vertex
two-port sources, not by extrapolating numerical differences.
No maximum width or catalog completeness is asserted.
