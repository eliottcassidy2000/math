# Target-placed mixed-scale unit pairs from two independently inflated vertices

This proves an exact strong two-choice repair and a constructive way to
select its absolute baseline in a relative interval around a target.
It is not an exact selection for every target.

Start with the six-vertex tournament of outgoing masks
`[52,41,42,17,6,24]` (orientation code 7002 in the low-bit
lexicographic convention). Replace vertices 4 and 5 by arbitrary
nonempty tournaments U and V, uniformly to the other vertices.
Call the resulting tournament A. Append two ordinal singletons
(buffer and closer), then reverse the closer's incoming arc from
either root 2 (low branch) or root 1 (high branch), but not both.
Every output is strong: the base tournament is strong, uniform
substitution by nonempty tournaments preserves strength, and the
buffer and closer feed back to A through the selected root.

Write p_r=P_r(U), q_s=P_s(V) for the numbers of ordered covers by
exactly r and s nonempty directed paths. Then the two exact counts
are
```text
F = Σ c_rs p_r q_s,     F+2 H(U)H(V)
c = [191 218  92 12
     186 290 172 30
      48 102  84 18].
```
In particular the relative step is always exactly 2 for transitive U
and V, without synchronizing their lengths or order parity.

## Proof of the coefficient table

As in SINGLE_SEAM.md, the closer reads h(A)+2 S_i(A). In h(A)
enumerate the paths after compressing each U/V run to one letter.
In S_i(A) enumerate the ordered pair of paths covering A, first
starting at root i and second possibly empty. Record one formal
variable per visit to U and V. Each gap between old vertices
permits the empty subword, U, V, or their two-letter subword in
the unique direction of the template edge between the modules.
Check the arcs at the ends of the gap and prohibit an initial
module gap in the rooted path. The same test handles an empty
second path and paths' end gaps. Every compression is unique,
and its expansions are counted by p_r q_s. Thus the finite
gap calculation suffices. The component coefficient tables
(rows r=1,2,3; columns s=1,2,3,4) are
```text
H(A)                     S₂(A)                   S₁(A)-S₂(A)
43 36  8 0               74  91 42  6            1 0 0 0
54 78 40 6               66 106 66 12            0 0 0 0
16 34 28 6               16  34 28  6            0 0 0 0
```
All other coefficients vanish. The replay counts these compressed
patterns separately from literal subset DPs on inflated tournaments.

## Moving the baseline from the target

Take U,V transitive of orders m,n≥1. An ordered path cover of
a transitive order is exactly a surjection to its labelled slots;
the compressed coefficients yield the exact exponential formula
```text
F(m,n)=51−20·2^n−24·3^n+6·4^n
       −70·2^m+80·2^(m+n)+16·2^m·3^n−24·2^m·4^n
       +24·3^m−42·3^m·2^n+12·3^(m+n)+18·3^m·4^n.
```
So F=18·3^m·4^n[1+O((2/3)^m+(3/4)^n)]. This is a
macroscopic selector with a simultaneous exact +2 repair.

For any sufficiently large real target X, one can explicitly
choose m,n with
```text
|F(m,n)-X| = O(X/log log X).
```
In particular every fixed relative interval about X contains
both branches for sufficiently large X. To prove the quantitative
version, let α=log₄3 and L=log₄(X/18). Both exponents can be kept
a fixed positive fraction of L. Start m at m₀=ceil(L/(4α)),
and allow an additional o(L) search range. Take
Q=floor(log L/(4 log 3)). By the pigeonhole principle some
1≤q≤Q has δ=||qα||≤1/Q. Distinct positive powers 3^q,4^p
differ by at least one for the nearest integer p to qα, so
δ≥c·3^(−q) by the mean value theorem for log on their bounded
ratio. Hence the multiples j q α for 0≤j≤floor(1/δ) form a
circular δ-mesh (possibly negatively oriented). The required
search range for m=m₀+jq is O(Q·3^Q)=o(L). Select the residue
nearest to L mod 1 and n nearest to L−mα. Both exponents grow
linearly in L and the log-scale error is at most 1/Q. The
nonleading terms contribute only exponentially small relative
error. This gives the stated error. The step from F to F+2
remains exactly 2 throughout.

For exact **membership** in these branches, fixing m makes F
strictly increasing in n. All visit coefficients are positive;
there is a positive coefficient with r=1,s=2, and all
transitive surjection counts are nondecreasing in n. Also
using 3^m=Σ binom(3,r) p_r and similarly for 4^n gives
F≥(2/3)3^m4^n, since each c≥12 and all binomial products ≤18.
So at a prescribed odd target N, check O(log N) choices of m,
binary-search n in its O(log N) range, and accept exactly a
literal F or F+2. This subfamily still has only O((log N)²)
distinct counts up to N. The relative interval theorem does
not imply arbitrary exact placement.

Replay: `python verify_mixed.py`. It runs the seam controls,
an independent two-module compressed enumeration, direct
tournament path-count comparisons, and exact inversion tests.
