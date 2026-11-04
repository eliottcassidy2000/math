# Breaking the mirror is not enough: exact one-feedback synthesis

Session 175. All general assertions below are proved; computed lists are explicitly finite certificates for a class with a proved size bound. The repository already proves that 7 and 21 are globally impossible (THM-1370). This packet does not prove the remaining completeness assertion that every other positive odd integer is attained. No historical novelty claim is made.

## 1. A fully controllable, nonuniform one-feedback constructor

Take the transitive tournament on vertices 1,...,n, with i→j for i<j, and add z. Encode z→i by b_i=1 and i→z by b_i=0. Call the result T(b). This includes **every tournament whose deletion of one vertex is transitive**, with a chosen deletion vertex and transitive order. Unlike the Session 174 mirrored wrapper, the intermodule orientations are not constrained to be uniform.

For a nonempty word b, let I={i:b_i=0}, O={j:b_j=1}. Every Hamilton path consists of an increasing subset S, then z, then the increasing complementary subset. If both parts are nonempty, put i=max S and j=min S^c. There are precisely:

- one path when i<j, which forces j=i+1 and b_i b_j=01;
- 2^(i-j-1) paths when j<i, provided b_j b_i=10: the vertices strictly between j and i can be assigned to either part;
- b_1 paths with z first and 1-b_n paths with z last.

Consequently the exact positive formula is

    H(T(b)) = b_1 + 1-b_n + #{i:b_i b_(i+1)=01}
              + sum_(j<i, b_j=1,b_i=0) 2^(i-j-1).       (1)

This is a literal bijection, not an arbitrary formal path profile.

**Streaming theorem.** Start with (P,L,h)=(0,0,1). Reading a bit c applies

    if c=0: (P,L,h) -> (2P,0,h+P+L);
    if c=1: (P,L,h) -> (2P+1,1,h).

Then h=H(T(b)) and P is the binary value of b. To verify the update, an appended zero creates new inversion pairs of total weight P; the terminal contribution increases by L. An appended one creates a 01 adjacency exactly when the terminal contribution disappears. The initial empty-word convention is the singleton z, with H=1. Thus the algorithm is linear in word length in arithmetic operations, with exact integer registers.

## 2. The finite-core theorem gives a complete inverse solver

Initial zeros are universal sources and terminal ones universal sinks. Deleting them preserves H and leaves either an empty word (H=1), or a core beginning with 1 and ending with 0. Such a core is strongly connected: z reaches vertex 1, vertex 1 reaches every later base vertex, every earlier base vertex reaches the last vertex, and the last vertex reaches z.

**Bound.** For a core of length n≥2, equation (1) counts 2^(n-2) paths from the inversion pair (j,i)=(1,n), plus the two paths with z at an end. Hence

    H(T(b)) >= 2^(n-2)+2.                              (2)

For any bound B≥3, every attainable count at most B therefore has a representative whose core length satisfies

    n <= 2 + floor(log_2(B-2)).

Enumerating those finite words is a **complete decision procedure for the entire deletion-one-transitive class**, even with arbitrarily many vertices. It is not merely a search through a fixed order. The code uses integer bit lengths rather than floating point logarithms.

In particular the class cannot realize 7, 13, 21, 29, 35, 41, 43, 49, 55, or 57. The certificate supplies the entire excluded list through 199. This is an obstruction to this proposed route, not a global exclusion: 13, for example, has a five-vertex witness outside this class, checked independently by subset DP in the certificate.

## 3. Explicit infinite constructions and composition

The four-run core w=1 0^a 1^b 0 (a,b≥1) has

    H(T(w)) = (2^a+1)(2^b+1).                    (3)

After reading 1 followed by a zeros, the registers are P=2^a and h=2^a+1. Reading b ones gives P=2^(a+b)+2^b−1 without changing h. The final zero adds P+1. This proves (3). In particular these products have strong witnesses, whereas ordered sums only immediately give nonstrong witnesses.

More generally cyclically orient nonempty tournaments A→B→{z}→A, dropping the mirror requirement B=A*. Let a=H(A), c=H(B), u=F_A(2)-H(A), v=F_B(2)-H(B). The resulting tournament is strong: a directed walk can leave any module, traverse the other two cyclically, and enter any chosen vertex of the original module. Then

    H = uv+2ac.                                       (4)

Indeed, with q_X(2)=2!p_X(2), collapse consecutive runs in each module. The skeleton types give
3ac + q_A(2)c + a q_B(2) + q_A(2)q_B(2).
Since F_X(2)=2H(X)+q_X(2), this equals (4). A module cannot occur more than twice, since another return would require revisiting z. This also proves no skeleton types are missing.

For transitive A,B of positive orders r,s, F_A(2)=2^r and F_B(2)=2^s, so

    H=(2^r-1)(2^s-1)+2.                               (5)

The mirror-breaking choice r=2,s=3 gives 23≡7 (mod8), whereas a mirrored wrapper always has H≡3 (mod8). The symmetric control r=s=2 gives 11, consistent with the mirrored obstruction. The transitive/transitive specialization itself never has a factor 7: the nonzero Mersenne residues modulo 7 are 1 and 3, whose possible products are 1,2,3, never −2=5. To genuinely escape the norm obstruction, use the four-vertex B with outgoing masks (8,1,3,6), for which H(B)=5 and F_B(2)=32. Taking A transitive of order r≥1 gives the proved strong family

    H_r = 27·2^r−17.                                  (6)

At r=2, H=91=7·13, exhibiting an odd 7-adic valuation forbidden in the mirrored family. For r=3k+2, every count is divisible by 7 and

    H_(3k+2)/7 = (108·8^k−17)/7.

Since 8^k≡1+7k (mod49), this quotient is 6+3k modulo7. Hence v_7(H_(3k+2))=1 for every k≥0 except k≡5 (mod7). This is an infinite family in the difficult 7-times-odd lane, including 7·13, but does not prove that the cofactor is prime infinitely often or that every 7p occurs. The prime lane also has the certified examples H_1=37 and H_3=199; no general primality assertion is made.

Finally H(A⊕B)=H(A)H(B): every path visits all of A before all of B, since no edge returns from B to A. All products of realized counts are realized, with vertex counts additive. This is useful closure, but it does not supply additive closure or realize all odd integers. Formula (5), all one-feedback counts from §§1–2, and any separately certified seeds can be combined this way.

## 4. Scope and source

The falling-factorial profile and uniform-module skeleton method are inherited from [Session 174 spectrum](https://github.com/eliottcassidy2000/math/blob/a55a46e2a9fd3880980b26149c0c6680f3ec8909/05-knowledge/results/session174_palindrome_symmetry_20261103/spectrum/PROOF.md), and the norm obstruction from its sibling synthesis/NOTE.md. The present packet derives the arbitrary binary feedback constructor, its finite-core inverse solver, and the asymmetric formula. It neither claims to solve the global spectrum nor treats arbitrary numerical profiles as realizable.

## 5. Exact insertion accounting, including orphan paths

Let T be arbitrary and prescribe b_x=1 when a new vertex v beats x, and 0 otherwise. For a Hamilton path P of T, let d(P) count adjacent 10 patterns along its b-word. The number of legal positions for v is

    b_first + (1−b_last) + #01 = 1 + #10.

Call a Hamilton path of T+v an orphan if deleting v does not leave a Hamilton path of T. Such a path has internal consecutive x,v,y with x→v→y and y→x. The following formula is exact:

    H(T+v) = H(T) + sum_P d(P) + O(T,b),                (7)

where O counts orphan paths. This partition is disjoint and exhaustive. More explicitly, for nonempty proper subsets S, x∈S, y∉S with b_x=0,b_y=1,y→x, multiply the number of Hamilton paths of T[S] ending at x by the number in T[S^c] beginning at y, and sum. That is O(T,b). Each orphan has a unique S,x,y, so nothing is double counted.

The old-path insertion contribution alone is not the full increment. For a transitive three-vertex base and b=100, H(T)=1, sum d(P)=1, but H(T+v)=5 and O=3. Thus a purported scalar recurrence based only on insertion positions would miss three paths, and even the tempting formula H(T+v)=H(T)+2 sum d(P) fails. Endpoint-sensitive subset data is needed for exact control. This directly identifies what a prime-lane insertion proof must supply: control both the path transition statistic and the orphan term.

## 6. Reproduction

Run `python3 check.py > certificate.json` and `python3 -O check.py > certificate.optimized.json`; compare the receipts. The checker uses only the standard library, explicit exceptions rather than Python assert, literal subset-DP Hamilton counts for all binary words of length at most eight, all 225 pairs of distinct profiles realized at orders at most five, the complete bounded inverse solver, and explicit insertion/orphan controls. No software was installed.
