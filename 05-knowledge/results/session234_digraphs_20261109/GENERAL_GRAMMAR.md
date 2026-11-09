# A finite grammar for every fixed semicomplete Hamiltonian-path count

This is an elementary structural consequence of the SCC decomposition and the generalized Camion lemma. No novelty claim is made. Digraphs are finite binary adjacency matrices; paths visit every vertex exactly once, so loops are inert. A semicomplete digraph has at least one arc between every pair of distinct vertices. Isomorphisms preserve arc directions and, when included, loops.

## 1. Finite strong-core alphabet, without tournament floor theorems

Let H(D) denote the number of directed Hamiltonian paths. Every strong semicomplete digraph C on m>=2 vertices has a spanning directed cycle, by the self-contained proof DG4 in [PROOFS.md](PROOFS.md). Cutting that cycle at each of its m possible starting vertices gives m distinct Hamiltonian paths: their first vertices differ. This remains valid at m=2, where the two paths are uv and vu. Therefore H(C)>=m. A singleton has H=1, whether looped or unlooped.

The SCC condensation of a semicomplete digraph is a total order, and all cross-component arcs point forward. Hamiltonian paths concatenate component paths uniquely, so H(D)=product_C H(C). Every non-singleton strong component has H(C)>=2.

Fix h>=2. If H(D)=h, every non-singleton strong component C has H(C)=d for an integer divisor d>1 of h, and |C|<=d<=h. Thus the alphabet of possible nontrivial strong cores, up to isomorphism, is finite and can be determined by a finite exhaustive enumeration through order h. This is an effectiveness theorem, not a claim that enumeration through h is computationally practical.

Conversely, any ordered sequence C_1,...,C_k of these cores with product_i H(C_i)=h, together with arbitrary nonnegative singleton gap lengths g_0,...,g_k, gives a semicomplete digraph of count h by ordered sum. Every such digraph arises uniquely this way up to isomorphism. Strong components cannot be permuted by an automorphism: the condensation's total order must be preserved. In particular Aut(D)=product_i Aut(C_i) in the loopless case. Singletons have no nonidentity internal automorphisms.

The core-word length satisfies k<=Omega(h), where Omega counts prime factors with multiplicity. Furthermore, if d_i=H(C_i), then
sum_i |C_i| <= sum_i d_i <= product_i d_i = h.
For k=1 this is immediate. For k>=2, the elementary inequality a+b<=ab for integers a,b>=2 and induction prove the second inequality. Hence the total size of every nontrivial-core word is at most h.

## 2. Exact generating functions without loops

Let W_h be the finite set of ordered isomorphism-class core words with product of their H-values equal to h. For w=(C_1,...,C_k), put m(w)=sum_i |C_i| and a(w)=product_i |Aut(C_i)|.

Let u_h(n) count loopless isomorphism classes of order n and H=h. Its ordinary generating function is

F_h(x) = sum_{w in W_h} x^{m(w)} / (1-x)^{k(w)+1}.

Indeed, the k+1 independent singleton gaps each contribute 1/(1-x). Equivalently,

u_h(n) = sum_{w in W_h, m(w)<=n} binomial(n-m(w)+k(w), k(w)).


Let l_h(n) count digraphs on a fixed labelled n-element vertex set. Its exponential generating function is

L_h(x) = sum_n l_h(n) x^n/n!
       = sum_{w in W_h} x^{m(w)} / [a(w)(1-x)^{k(w)+1}].

Consequently l_h(n)/n! is given by the same binomial sum with each term weighted by 1/a(w). This is the orbit-stabilizer identity: an unlabelled digraph D has n!/|Aut(D)| labelled realizations. The generating function here is an EGF, not the ordinary generating function of labelled counts.

One can group the finite core alphabet by count d. Define
C_d(x)=sum_{strong loopless C up to isomorphism, H(C)=d} x^{|C|},
A_d(x)=sum_{same C} x^{|C|}/|Aut(C)|.
Then F_h is equivalently the sum over all ordered factorizations d_1...d_k=h, d_i>=2, of product_i C_{d_i}(x)/(1-x)^{k+1}; replace C by A for L_h. These polynomials have degree at most d.

## 3. Exact degree, thresholds and asymptotics

If W_h is empty, h is unattained and both generating functions are zero. Otherwise let K(h)=max_{w in W_h} k(w). Then 1<=K(h)<=Omega(h).

For n>=h, all word-size thresholds have been passed. Thus u_h(n) and l_h(n)/n! are exactly polynomials in n, of degree precisely K(h). The degree does not merely have an upper bound: every word of maximum length contributes a strictly positive leading coefficient, so cancellation is impossible.

The leading coefficients are respectively
(number of words w with k(w)=K(h))/K(h)!,
and
(sum_{w: k(w)=K(h)} 1/a(w))/K(h)!.

Equivalently, both nonzero rational generating functions have a pole at x=1 of exact order K(h)+1 and no other finite poles. The bound n>=h is a uniform sufficient threshold, not necessarily the smallest threshold. The sharper sufficient threshold is max_{w in W_h} m(w).

K(h) is the largest length of an ordered multiplicative factorization of h whose factors are attained by strong semicomplete cores. Equality with Omega(h) should not be assumed without an existence argument for the requisite prime-count cores. This grammar alone does not prove that every positive integer is attained.

## 4. Allowing arbitrary loops

Loops do not change H or the SCCs, but they do change isomorphism classes and automorphism groups. A singleton gap is now an ordered word in two singleton types (unlooped or looped), hence contributes 1/(1-2x).

For a loopless strong core C of order m, let G=Aut(C) act on its 2^m binary loop assignments. Its number of decorated isomorphism classes is
q(C) = (1/|G|) sum_{g in G} 2^{c(g)},
where c(g) is the number of cycles of g on the vertices (Burnside's lemma). In general q(C) is not 2^m: distinct loop assignments can be isomorphic. A decorated core's automorphism group is the stabilizer of its loop assignment.

Let W_h^loop be the finite set of ordered decorated-core words with count product h. Writing m(w), k(w), and a(w) for their sizes, lengths and products of decorated automorphism orders, the exact formulas become

F_h^loop(x) = sum_{w in W_h^loop} x^{m(w)}/(1-2x)^{k(w)+1},
L_h^loop(x) = sum_{w in W_h^loop} x^{m(w)}/[a(w)(1-2x)^{k(w)+1}].

Alternatively, group each undecorated word w: its contribution to the ordinary isomorphism OGF is
[product_i q(C_i)] x^{m(w)}/(1-2x)^{k(w)+1}.
For labelled EGFs, orbit-stabilizer gives
sum_{loop-decoration orbits on C} 1/|Aut(decorated C)| = 2^{|C|}/|Aut(C)|.
Therefore L_h^loop(x)=L_h(2x), equivalently l_h^loop(n)=2^n l_h(n), as also follows directly by choosing the n diagonal bits of each labelled matrix independently.

For n>=h, the ordinary isomorphism count u_h^loop(n) is exactly 2^n times a polynomial of degree K(h), and l_h^loop(n)/n! is likewise 2^n times a polynomial of degree K(h). Specifically, a decorated word contributes
2^{n-m(w)} binomial(n-m(w)+k(w), k(w)).
All coefficients at maximum degree are positive. Both nonzero generating functions have their sole finite pole at x=1/2, of order K(h)+1. The degree K(h) is unchanged because adding loops never changes the available core-word lengths.

## 5. Counts one and zero; empty-graph convention

For h=1, every strong component is a singleton, so the loopless graphs are exactly transitive tournaments. There is one isomorphism class of each order and n! labelled graphs. With the convention H(empty graph)=1, both F_1(x) and L_1(x) equal 1/(1-x). Allowing loops gives 2^n isomorphism classes and 2^n n! labelled graphs, and both corresponding generating functions equal 1/(1-2x). If the empty graph is excluded, subtract the constant term one.

No nonempty semicomplete digraph has H=0: each strong core has a Hamiltonian path (indeed a Hamiltonian cycle if nontrivial), and component paths concatenate. With the stated empty-graph convention the count-zero family is empty altogether.

## 6. Relation to the seven/twenty-one classification

The general theorem bounds strong cores for count seven by order seven and for count21 by order21, already proving finite reducibility without Busch or tournament floor results. The stronger spanning-strong-tournament comparison sharply reduces those bounds to order four for H7 and order six for H21. The detailed census then supplies the exact alphabets.

The known H7 alphabet contains three strong cores and gives K(7)=1. The H21 grammar has single strong cores and two-core words using an H3 core and an H7 core, so K(21)=2. There are two H3 core types, not just the directed triangle: its one-edge-doubled version also has H3 because the added reverse arc is inert. Thus there are twelve ordered two-core types for H21 before singleton gaps. These examples illustrate the general rational grammar; they do not establish a universal integer spectrum.

## 7. Small sanity check: h=2

The only strong core with H2 is the bidirected pair: its order is at most two by Section1, and its automorphism group has order two. Thus
F_2(x)=x^2/(1-x)^2, L_2(x)=x^2/[2(1-x)^2].
In particular u_2(n)=n-1 and l_2(n)=n!(n-1)/2 for n>=2.

The bidirected pair has three loop-decoration orbits (zero, one or two looped vertices), not four. Consequently
F_2^loop(x)=3x^2/(1-2x)^2,
L_2^loop(x)=2x^2/(1-2x)^2.
For n>=2 these give u_2^loop(n)=3(n-1)2^{n-2} and l_2^loop(n)=n!(n-1)2^{n-1}. This verifies the distinction between unlabelled decoration orbits and independent labelled diagonal choices.
