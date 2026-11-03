# Cyclic three-block tournaments: exact paths and feedback number

## Scope

This note proves exact Hamilton-path counts for actual strong tournaments, not merely residue coverage. No claim of novelty, of a new gap, or that these integers require high feedback vertex number in every tournament is made. The triangular construction is a concrete way to work with the odd-cycle inspiration; it does not assert a connection to Clifford or octonion algebra.

## Construction and proof

Let T_n be a transitive tournament on n vertices. For a,b,c >= 1 construct C(a,b,c) by substituting T_a,T_b,T_c into the three vertices of a directed triangle A -> B -> C -> A. The result is strongly connected. Its minimum feedback vertex set has size exactly min(a,b,c): if every class remains nonempty then picking one vertex from each gives a directed triangle; deleting a smallest class leaves the transitive ordinal sum of the other two classes.

For n,k >= 0 let Q(n,k)=k! S(n,k) denote the number of surjections from an n-set to k ordered nonempty parts, with Q(0,0)=1. Define U=e^x-1, V=e^y-1, W=e^z-1. The number H(a,b,c) of directed Hamilton paths is

    a! b! c! [x^a y^b z^c]
      (U+V+W+UV+VW+WU+3UVW)/(1-UVW).

Indeed break each Hamilton path into maximal consecutive runs of vertices from the same class. Since the three classes form a directed triangle, the sequence of classes of these runs is forced to follow A,B,C,A,... starting from one of the three classes. A run can contain any nonempty subset of the class, but its internal order is uniquely the transitive order. The same class may recur, and assigning its labeled vertices to its successive runs is an ordered set partition, not necessarily an interval partition. Each run thus contributes U, V, or W to the exponential generating function. Grouping by number of complete turns gives the displayed numerator and geometric denominator. This argument includes starts and finishes in every class and counts each path exactly once.

Equivalently, put A_k=Q(a,k), B_k=Q(b,k), C_k=Q(c,k), and sum over k >= 0:

    H = sum_k [3 A_k B_k C_k
      + A_(k+1) B_k C_k + A_k B_(k+1) C_k + A_k B_k C_(k+1)
      + A_(k+1) B_(k+1) C_k + A_k B_(k+1) C_(k+1)
        + A_(k+1) B_k C_(k+1)].

The sum is finite. In the first term k=0 contributes zero because a,b,c are positive. The three classes in each group of terms correspond to the possible start and end positions of the cyclic run word.

## Explicit exact subfamilies

Using Q(n,1)=1, Q(n,2)=2^n-2, Q(n,3)=3^n-3*2^n+3, the two smallest class sizes give:

    H(a,b,1) = 2 + (2^a-1)(2^b-1),

and

    H(a,b,2) = 2*3^(a+b) - 4*(3^a*2^b + 2^a*3^b)
               + 13*2^(a+b) + 2*(3^a+3^b)
               - 11*(2^a+2^b) + 15.

In particular for b >= 2 a strong tournament with feedback vertex number exactly two realizes

    H(2,b,2) = 4*3^b + 5*2^b - 11,

starting 45, 137, 393, 1121, 3225, 9377. The formula is strictly increasing in b. Equal block sizes give a fully explicit witness family with arbitrarily large feedback vertex number:

    (n, H(n,n,n)) = (1,3), (2,45), (3,2721),
                     (4,421425), (5,132518793), ... .

This high-feedback claim concerns the displayed witnesses. It does not say the same count cannot have a different witness with a smaller feedback set.

## Reproduction

Run `python3 check.py` and `python3 -O check.py`. The script independently orients each edge of the literal tournament, counts Hamilton paths by subset/endpoint dynamic programming, and checks the ordered-surjection formula for all 64 triples in {1,2,3,4}^3 and the special formulas when c=1 or c=2. These finite checks complement the proof; they do not establish any infinite claim by themselves.

## A 3-adic consequence of the two-feedback family

For every e >= 1, N_b=H(2,b,2) attains exactly the residues r modulo 3^e with r != 1 modulo 3, each infinitely often as b ranges over integers b >= 2. For b >= e, the 4*3^b term vanishes modulo 3^e. The residue 2 has multiplicative order 2*3^(e-1) modulo 3^e, the full order of the unit group, since v_3(4^(3^j)-1)=j+1 by induction. Thus 5*2^b ranges over all units modulo 3^e as b ranges over a complete period. The possible r are exactly those for which r+11 is a unit, or equivalently r != 1 modulo 3. Replacing b by a large enough congruent integer preserves the residue and allows b >= max(2,e).

In particular, for every t >= 1, taking e=t+1 and r=3^t gives a strongly connected C(2,b,2) with feedback vertex number exactly two and v_3(H)=t. Taking r=2 modulo 3 gives valuation zero. Taking r=0 gives examples divisible by arbitrarily high powers of 3. This is a congruence consequence of the explicit family, not an exact classification of the spectrum.

The companion [carrier packet](../carrier/NOTE.md) gives the same cyclic gluing identity for arbitrary input tournaments; this transitive-block specialization uses C_k(T_n)=k! S(n,k).
