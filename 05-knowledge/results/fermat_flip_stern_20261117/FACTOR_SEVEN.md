# Two-slot compensated factor-seven switch

For nonempty A and B with two-cover profiles (u,R) and (v,S),
respectively, the inherited triangle substitution identity is

    H(C3[A,B,1]) = (S-v)R + (3v-S)u
                = (R-u)(S-v) + 2uv.

The SCC ordinal product multiplies both coordinates (each color-class
Hamilton path crosses adjacent SCCs in their forced order).
Consequently the following physically exact profile arithmetic makes
a switch valid for **every** surrounding nonempty A:

    B  = (11,800) ⊕ (53,1640) ⊕ (97,2584)
    B' = (77,1904) ⊕ (53,1520) ⊕ (97,4100) ⊕ (1,2).

The parentheses in this display give each actual block's pair
(H,P(2)); see the standalone check.py for literal outgoing rows
and the independent subset-DP verification of all seven blocks.
The unpadded product increases H by 7 but has half the desired
P(2) ratio; the last singleton exactly restores the missing factor 2.
The reason is an explicit odd-prime factor flow. The successive
lower/upper odd mantissa ratios are

    25/(7*17),  41/19,  (17*19)/(25*41).

They cancel exactly to 1/7. The H factors of the second and third
pair are unchanged throughout this correction. The resulting
profiles are (56551,3390208000) and (395857,23731456000), the latter
exactly seven times the former. Both cyclic tournaments are strong,
since each vertex can reach its own block via a directed path
around the other two nonempty blocks. Their Hamilton counts
therefore obey exact sevenfold scaling under a replacement of
the marked B interface, independently of the choice of A.

This is genuinely prime-useful, not the obvious scaling of a
multiple of 3 or another common factor. The pair (56551,3390208000)
is primitive; in fact its H factors are 11,53,97 and its P(2)
prime factors lie in {2,5,17,19,41}. With A transitive of order
1,4,6 the smaller output is respectively 3390264551, 50852384837,
213579654389, all prime by literal trial division. The larger
outputs are their exact multiples by seven.

## Limits and what changed

A fixed B→C switch within this cyclic architecture scales every
surrounding A by seven **iff** the two-coordinate profiles are
scaled by seven. Necessity already follows from the two choices
A=TT1 and TT2: subtract their equations to obtain scaling of S-v,
then obtain scaling of 3v-S, an invertible two-by-two linear test.
Sufficiency follows by the displayed gluing identity. In particular,
merely inserting extra ordinal summands into B cannot do it:
multiplicativity of the first coordinate would require the product
of inserted tournaments to have count 7, which is forbidden.
The displayed switch changes existing summands and uses two
independent same-H P(2) corrections instead.

This does not contradict the seven-gap atlas's single-module
amplification no-go. That theorem would enlarge a fixed homogeneous
module with H=p and then try to read 7p. Here p is the output of
C3[A,B,1], not the Hamilton count of a p-path module contained in
the larger tournament. The two different host profiles carry 56551
and 395857 paths, and both are read using the one-visit and two-visit
terms. In particular the replacement is not a claim of induced
containment or retaining the original B. The construction realizes
absolute sevenfold lifting for this entire profile family, but
does not show that every prescribed p occurs in the smaller readout.
For a target p the exact inverse still requires a physical A pair
with R=(p+(S-3v)u)/(S-v), with u positive and satisfying the
integrality congruence.

## A density barrier for assemblies of one-apex strong SCCs

An independent restriction applies if both modules' nonsingleton
SCCs are strongly connected one-apex tournaments. Such a strong
core of order n>=3 has a transitive deletion spine whose first
vertex is an apex outcontact and whose last is an incontact
(otherwise the spine extremes could not lie in the same SCC).
The odd cycles using those two vertices and the apex have any
even subset of the n−3 intermediate vertices. By the inherited
odd-cycle packing identity, evaluated at 2, these alone contribute
at least

    2^(n-1) * sum_{k even} binom(n-3,k)2^(-k)
      = 2*(3^(n-3)+1)

to the core P(2). There are at most 2^(n-3) distinct strong one-apex
cores, from pointed binary spine words. For any sigma>log_3(2),
the sum of P(2)^(-sigma) over the allowed strong cores converges.

A whole ordinal profile is determined by a multiset of allowed SCC
types, since its two coordinates are products and are insensitive
to the order of the components. The Euler product
prod_C (1-P_C(2)^(-sigma))^(-1) converges by the previous bound
(and includes the singleton with P(2)=2). This product upper-bounds
the Dirichlet sum over distinct whole profiles; double it as a
product over two independent modules. Since R>=2u and S>=2v, the
cyclic output is at least RS/4. Thus the number of distinct outputs
at most X from this entire two-module grammar is
O_epsilon(X^(log_3(2)+epsilon)). In particular, unrestricted
numbers of one-apex SCC factors on both sides still cannot cover
the prime 7p targets. This is a capacity theorem for this grammar,
not an obstruction to arbitrary tournaments of general structure.
