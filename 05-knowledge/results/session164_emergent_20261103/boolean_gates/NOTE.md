# Canonical cut geometry and a sharp strong Boolean escape gate

**Status.** All-order theorems and literal tournament constructions, with independent finite controls. This does not solve the Hamilton-path spectrum problem, supply all three ternary branches at a returned physical source, or assert literature priority. The new levers are an orientation-preserving carrier with a deletion-minor floor, a general linear-forest support rule, and an all-strong family attaining its escape-density bound. All graphs use row masks: bit j of row i means i→j.

## 1. Source credits and scope

The one-apex quadratic response, its two-path-cover coefficients, submodularity, and an integral directed-cut representation are inherited from [S152, integral carrier](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session152_sequences_20261103/cut_induction/README.md) and [S152, reciprocity](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session152_sequences_20261103/recursions/RESULT.md). S152 obtains the carrier by integral max flow and does not require its orientation to follow the tournament. Section 2 supplies a canonical choice supported exactly on the old tournament, with an edgewise deletion bound and a complete equality classification.

The motivating one-exception degree bound and modular-interface questions are in [S163, marked escape](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/nonuniform/MARKED_ESCAPE.md). Sections 3–4 extend the degree argument to arbitrary variable-edge graphs, prove a sharp density theorem, and give new exact Boolean-gate families. The wrapper G(X)=X→a→b→X and multiplicative ordered-slot profiles are inherited from [S162, square transfer](https://github.com/eliottcassidy2000/math/blob/3b532f92e3170b15771d9ad38182688efe1f526e/05-knowledge/results/session162_mrb_square_20261103/square_transfer/NOTE.md). These are targeted source checks, not a comprehensive novelty search.

## 2. A canonical carrier follows every original arc

Let A be a nonempty tournament. For distinct u,v, let Q_uv count ordered covers (P,R) of V(A) by two nonempty disjoint directed paths, with P ending at u and R starting at v. Write K_uv=Q_uv+Q_vu. For D⊆V(A), add an apex z beating exactly D, and call its Hamilton-path count h(D). Let H=H(A).

**Theorem 1.** For every old arc u→v, put c_uv=K_uv/2 and c_vu=0. Then c is an integral carrier and

    h(D) = H + 2 sum_{u→v, u∈D, v∉D} c_uv.                 (1)

Moreover, with H(empty tournament)=1,

    c_uv ≥ H(A−{u,v}) ≥ 1.                                 (2)

Thus the positive-support digraph of this canonical carrier is exactly A. No flow allocation is needed.

**Proof of (1).** Let s_v,e_v count old Hamilton paths starting and ending at v. Splitting a child path at z gives the inherited identity

    h(d)=H+sum_v (s_v−e_v+sum_u Q_uv)d_v
              −sum_{u<v}(Q_uv+Q_vu)d_ud_v.

If u→v, joining a two-path cover at u,v is a bijection with old Hamilton paths containing the consecutive arc u→v. Consequently

    H=s_v+sum_{u→v}Q_uv,
    H=e_v+sum_{v→u}Q_vu.

Subtracting these identities shows

    s_v−e_v+sum_u Q_uv = sum_{v→u}(Q_vu+Q_uv).

Substitution gives (1). Rédei parity makes K_uv even, since every child count is odd and its mixed Boolean difference is −K_uv. Hence c_uv is an integer.

**Proof of (2).** Fix a Hamilton path P=p_1,...,p_m of A−{u,v}. For every slot i=0,...,m and both choices (a,b)=(u,v),(v,u), form two candidate paths

    p_1,...,p_i,a   and   b,p_{i+1},...,p_m.

Let R_i be the subset of {u,v} beaten by p_i, and set R_0={u,v}, R_{m+1}=∅. Write σ for swapping u,v. The number of legal choices at slot i is exactly

    |R_i \ σ(R_{i+1})| ≥ |R_i|−|R_{i+1}|.

The sum over slots is at least 2 by telescoping. Each legal choice gives a cover counted by Q_uv+Q_vu. Covers produced by different P, slots, or choices are distinct: delete their designated endpoints and concatenate the remainders to recover P and the split. Thus K_uv≥2H(A−{u,v}). The empty-path case m=0 gives the two singleton covers directly. ∎

The lower bound is sharp, for example for consecutive vertices of a transitive tournament, where c_uv=H(A−{u,v})=1.

### Equality and strongly connected components

By positivity of all c_uv, h(D)=H precisely when A has no arc from D to its complement. If A has k strongly connected components, ordered C_1→...→C_k, these D are exactly the k+1 terminal unions

    ∅, C_k, C_{k−1}∪C_k, ..., V(A).

So exactly k+1 one-apex orientations preserve H. They simply place the apex between existing strong components; there are no other count-preserving interfaces. More quantitatively,

    h(D)−H ≥ 2 sum_{u→v, u∈D, v∉D} H(A−{u,v}).            (3)

### A constructive minimum-cut consequence

Suppose z must beat a specified set L and lose to a disjoint set R. Minimize the c-weight directed cut with L on the source side and R on the sink side. An ordinary integral s–t minimum cut, enforcing each terminal with capacity M>sum c, gives the exact smallest possible h(D), by (1). This constructs the minimizing tournament interface, not just a lower bound.

The minimum strong extension is also computable from these cuts. A+z is strong exactly when z beats some vertex of C_1 and loses to some vertex of C_k. Enumerate one forced source-side vertex in C_1 and one forced sink-side vertex in C_k (distinct if k=1), and take the smallest resulting minimum cut. This is at most |C_1||C_k| minimum-cut problems; if no distinct forced endpoints exist, report infeasible (in particular, n=1). It optimizes a one-apex extension; arbitrary target counts and full observer returns are not thereby decided by a single minimum cut.

## 3. Linear forests control all orientation polynomials

Fix a vertex set, orient some edges permanently, and let E be a set of freely oriented edges, with one Boolean variable per edge. Its variable-edge graph is the undirected graph (V,E). Let P(d) be the exact Hamilton-path count of the resulting tournament.

**Theorem 2 (linear-forest support).** P has a unique integer multilinear Boolean polynomial. A monomial can have nonzero coefficient only if its selected edges form a linear forest: an acyclic graph of maximum degree at most two. Therefore

    deg P ≤ λ(E),

where λ(E) is the maximum number of edges in a linear forest contained in (V,E). In particular, if every variable edge touches a designated set S of s vertices,

    deg P ≤ min(|V|−1,2s).                                  (4)

If a selected monomial uses r variable edges internal to S, its degree is at most 2s−r.

**Proof.** For each vertex permutation, the indicator that it is a directed Hamilton path is the product of its required consecutive-edge indicators, each equal to d_e, 1−d_e, zero, or one. Every monomial appearing on expansion uses a subset of that permutation's consecutive edges, hence is a linear forest. Summing cannot create a forbidden support. For (4), every selected edge touches S, while the selected degrees at S total at most 2s. An internal edge uses two of these incidences. ∎

This is the intersection of a forest constraint with local degree constraints; it is not asserted that linear forests themselves form a matroid. It gives an exact sparse interpolation rule: compute Möbius coefficients only for linear-forest subsets, since all other coefficients vanish. The Boolean coefficients for a fixed marked observer (H,x,y) obey the same support and degree restrictions when the query-edge orientations are fixed: each query count does, and their always-even differences permit coefficientwise integer division by two.

### An arbitrary-modulus minimum-support lemma

Let f:{0,1}^N→G be a multilinear polynomial of degree at most d, with coefficients in any abelian group G. If f is not identically zero, then

    #{a:f(a)≠0} ≥ 2^{N−d}                                  (5)

when d≤N. This is the Boolean low-degree minimum-distance argument, and it does not require a field.

**Proof.** If d=N the claim is immediate. Otherwise induct on N. If f is independent of a selected variable, both slices are equal, and induction doubles the support. Otherwise write the two slices as f_0,f_1. Their difference g=f_1−f_0 is nonzero and has degree at most d−1. At each point where g≠0, at least one slice is nonzero. Thus support(f)≥support(g)≥2^{(N−1)−(d−1)}. The degree-zero case is immediate. ∎

Apply this to P modulo any integer m>1, or to P−c modulo m. If any allowed interface escapes a specified count congruence, at least a fraction 2^{−λ(E)} escape. With s exceptional vertices this fraction is at least 4^{−s}. The same conclusion holds for a fixed full observer modulo m, taking G=(Z/mZ)^3. It concerns failure of divisibility or of one congruence, not simultaneous coprimality with every prime divisor of a composite modulus. In particular it does not repair S163's separate-prime versus simultaneous-unit obstruction.

A uniform random independent interface search therefore misses an existing escape after r trials with probability at most (1−2^{−λ(E)})^r. This is a success guarantee conditional on existence, not an existence proof.

## 4. A five-vertex strong Boolean gate

Start with the strong four-vertex tournament

    Q=(8,1,3,6).

Add z=4 with z→1, 3→z, z→0 iff x=1, and z→2 iff y=1. Set t=y(1−x). The four literal rows and exact profiles are

    (x,y)=(0,0): (24,1,19,22,2), H=9,  F(2)=84;
    (x,y)=(1,0): (8,1,19,22,3),  H=9,  F(2)=84;
    (x,y)=(0,1): (24,1,3,22,6),  H=11, F(2)=100;
    (x,y)=(1,1): (8,1,3,22,7),   H=9,  F(2)=84.

Here F_B(r) counts covers of B by r ordered directed-path slots, with empty slots allowed. Thus this literal gate satisfies

    H(B)=9+2t,          F_B(2)=84+16t.                       (6)

Every setting is strong, because Q is strong, z beats old vertex1, and old vertex3 beats z. The variable-edge graph is the length-two path 0–z–2. This is a physical two-input Boolean detector, not merely an assignment of integers to inputs.

The gate's full marked observer at its apex z is

    v(B,z)=(9+2y−2xy, 4+3x+5y, 7+2x+4y).                 (7)

Thus the three inactive scalar states are NOT equal marked-source states. In particular (6) does not give a closed physical source bank. Formula (7) explicitly records that limitation.

### Uniform substitution turns the gate into a modular singleton

For every gate, C_1=9+2t and C_2=F(2)−2H=66+12t. Thus modulo3,

    C_1 ≡ 2t,       C_k ≡ 0 for every k≥2.

The latter assertion for k≥3 follows from k!|C_k. Consequently, replace any s vertices of any fixed tournament Q by independent gate modules, orienting all cross-module edges according to Q and leaving other vertices as singletons. Partition paths by their sequence of module runs. Modulo3 only one-run use of each gate survives, so the exact general congruence is

    H(Q[B_1,...,B_s]) ≡ H(Q) 2^s product_i t_i (mod3).     (7a)

This is a compositional substitution rule: the entire gate module acts as a singleton of weight2t modulo3. If Q is strong and H(Q) is not divisible by3, all substituted tournaments are strong and the count detects the AND of the s gate activations. A marked port outside the replaced modules gives the identical scaling rule for the full observer, since the same reasoning applies to both query graphs. No assertion is made for nonuniform cross-module interfaces, which can break the run factorization.

### Exact all-strong AND composition and valuation distribution

Take s≥1 independent gates B_1,...,B_s. Let n be 0 or 1 with n≡s mod2, and form

    Y=B_1⊕...⊕B_s⊕T_n,
    Z=G(Y),     G(Y)=Y→a→b→Y.

T_0 means the empty padding. Every Z is strong, on 5s+n+2 vertices, for every one of the 4^s choices of the 2s interface bits. The same s apex vertices support all variable edges.

The wrapper formula is elementary: paths having a immediately before b correspond to all ordered two-slot covers of Y, including empty slots; paths having b before a must have b first and a last and correspond to Hamilton paths of Y. Hence H(G(Y))=H(Y)+F_Y(2). Both H and F(2) multiply under ordinal sum, while F_{T_n}(2)=2^n. Consequently, if exactly k of the s gates have t_i=1,

    H(Z)=9^(s−k)11^k + 2^n 84^(s−k)100^k.                 (8)

The s+1 counts in (8) strictly increase with k. Exactly binom(s,k)3^(s−k) interfaces yield the kth count, because an active gate has one assignment and an inactive gate has three.

**Theorem 3 (sharp density, all members strong).** For the family (8),

    v_3(H(Z))=s−k,
    sum_interfaces w^{v_3(H(Z))}=(1+3w)^s.                 (9)

In particular, precisely one of the 4^s interfaces has H(Z) nonzero modulo3. Thus the bound 4^(−s) in (5) is sharp even among families in which every resulting tournament is strongly connected.

**Proof.** Put r=s−k. For r≥1, factor (8) as

    H(Z)=3^r[3^r11^k+2^n28^r100^k].

The bracket is nonzero modulo3. If r=0, (8) is congruent to (−1)^s+(−1)^n=2(−1)^s, nonzero by the parity choice. The multiplicities give the generating polynomial in (9). ∎

### The same sharp law for a fixed physical observer

Retain either wrapper vertex as the observer port, outside every gate block. Every inactive gate is 3-run-divisible: its ordered nonempty path-cover counts satisfy C_1=9, C_2=F(2)−2H=66, while k! divides C_k and hence 3 divides C_k for k≥3. In the source and both observer query graphs, all outside vertices see each gate block uniformly. Partition any Hamilton path by its sequence of block runs. Each of the r inactive blocks contributes a factor divisible by3 independently, proving 3^r divides all three query counts. Since 2 is invertible modulo 3^r, 3^r divides each observer coordinate.

The first coordinate H has valuation exactly r by (9), so

    v_3(gcd(H,x,y))=r.

Consequently the full marked observers have the same valuation enumerator (1+3w)^s, and exactly one of the 4^s fixed-port interfaces escapes the zero observer modulo3. This realizes sharpness for vector-valued (5), too. It is a modular zero detector at an unchanged physical port; distinct integer observer states are not identified.

Examples: s=1 gives strong counts 177 and 211, with multiplicities 3 and 1; s=2 gives strong counts 7137, 8499, 10121, with multiplicities 9, 6, 1 and valuations 2,1,0. These are exact count-construction families, not claims of interval coverage or ternary affine source closure.

## 5. Reproduction

Run `python3 check.py` and `python3 -O check.py`. It independently constructs the literal tournaments; computes Hamilton paths with endpoint/subset DP; computes Q_uv by subset convolution; checks the canonical carrier, deletion floors, and SCC equality on every tournament through order4 plus selected larger examples; checks the four gates and their marked observers; tests all sixteen s=2 wrapped settings by literal 12-vertex DP and their fixed-port observer queries; and verifies forest support on a full six-variable order-four cube. Explicit failure checks remain active under Python -O. General proofs, not finite checks, establish all-order assertions.
