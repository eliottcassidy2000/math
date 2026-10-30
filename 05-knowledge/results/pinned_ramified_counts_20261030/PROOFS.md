# Pinned cycle indices and the order-p rigidity test

2026-10-30. **PROVED scoped; same-agent mathematical review.** These are
statements about the universal directed-color group H, not an assertion
that every tournament in its bank has full automorphism group H.

## P1. Inherited source and exact question

Use [U1--U2](../universal_ramified_closure_20261030/PROOFS.md) at commit
`c830320f16`: prime p>=5, k>=1, v_p(D)=1,

    R=(Z/p^k)[delta]/(delta^2-D), N=p^(2k), rho(x+y delta)=x mod p,
    g=(1+delta)/(1-delta), eta=g^(p^(k-1))-1,
    H={z -> g^t z+b+eta h(rho(z))}.

This acts on the literal N vertices. Directed pair colors are the orbits
of differences under <g>. Choosing opposite-orbit signs gives actual
tournaments; there are no ties or artificial orientations. H is their
common automorphism group, with order p^(3k+p-2).

Let S be a NONEMPTY set of distinct named vertices, s=|S|, and G=H_(S)
its pointwise stabilizer. Write j=|rho(S)|. If j=1, let tau be the minimum
delta-valuation of pair differences, with tau=2k for a singleton, and
ell=floor((tau+1)/2). Then 1<=ell<=k. Inherited U2 gives

    j>=2: G=C_p^(p-j);
    j=1:  G=C_(p^ell) x C_p^(p-2).

The new task is the full permutation cycle index, not just these orders
or abstract group types. Let x_i mark a cycle of length i, and normalize
Z_G by |G|. Put

    A=N/p, B=N/p^2, J=x_1^A+(p-1)x_p^B.

**Cycle-index theorem.** For j>=2,

    Z_G = x_1^(jA) J^(p-j) / p^(p-j).                       (P1.1)

For j=1,

    Z_G = [x_1^A J^(p-1)
           +p^(p-2) sum_(r=2..ell) (p-1)p^(r-1) M_r]
          /p^(ell+p-2),                                   (P1.2)

where the empty sum is zero and

    M_r = x_1^(p^(2k-2r+1))
          product_(i=1..r-1) x_(p^i)^((p^2-1)p^(2k-2r+i-1))
          x_(p^r)^((p-1)p^(2k-r-1)).                       (P1.3)

Thus occupancy plus LCA determines this cycle index. This does NOT
classify pin configurations up to H; directed pair colors still separate
configurations with the same j and tau, as U2's (0,1)/(0,2) hostile shows.

### Proof

For j>=2 the inherited stabilizer consists of independent socle shifts
on the p-j unoccupied residue fibres. The occupied fibres are fixed
pointwise. On a free fibre, the identity has A singleton cycles; each
of its p-1 nonidentity translations has B cycles of length p. Disjoint
support multiplies the cycle enumerators, giving (P1.1).

For j=1 translate a pin to zero. Every pin now has residue zero. Write
a map fixing zero as g^t z+eta h(rho(z)), h(0)=0. Fixing the pins requires
t a multiple of p^(k-ell). The representation has multiplicity p because
the top clock is a linear fibre shift. The maps whose clock has order
at most p are exactly F_0, all h with h(0)=0; their enumerator is
x_1^A J^(p-1).

If the clock order is p^r, r>=2, adding the socle shift does not change
its cycle partition counts. On the zero fibre, h is zero. On a nonzero
residue fibre, (g^t-1)z has valuation 2k-2r+1<2k-1, so eta h cannot
cancel it. For every power divisible by p the shift vanishes and the
power is the corresponding pure clock. For powers prime to p the same
valuation argument applies. Fixed-point counts of every power, and
therefore cycle counts, agree with the clock.

The fixed-point count of the p^i-th power is p^(2k-2r+2i+1) for
0<=i<r, and N for i=r: it is the size of the annihilator of g^(tp^i)-1.
Successive differences divided by p^i give (P1.3). There are
(p-1)p^(r-1) such clocks and p^(p-2) distinct maps per clock after the
overlap is divided out. This proves (P1.2), including k=1. QED.

This root calculation avoids the unrooted depth-one cycle-type failure
recorded in the preceding ramified packet; it does not extend that false
unrooted formula to k=1.

## P2. Attachment counts and a short LCA recurrence

Set every x_i=Q to count G-orbits of colorings with a named palette of
Q colors. Denote the result B_(p,k,j,ell)(Q), omitting ell when j>=2.
For any prescribed colors on S, the number of extensions is

    B_(p,k,j,ell)(Q) / Q^s.                              (P2.1)

Every group element fixes all s pins, so their Q^s independent choices
factor from every fixed-coloring count. This division is exact, not a
claim that the color assignment can be changed by a palette symmetry.
Repeated pins first require compatible assignments and then count once.

For one occupied residue fibre, fix p,k and write B_ell(Q). Then

    B_1(Q)=Q^A [Q^A+(p-1)Q^B]^(p-1) / p^(p-1),
    p B_ell(Q)=B_(ell-1)(Q)
                   +(p-1)Q^(2p^(2k-ell)-p^(2k-2ell)),   ell>=2. (P2.2)

The exponent is the total number of cycles in M_ell; summing (P1.3)
gives it exactly. This adds one clock layer without enumerating its
permutations. For multiple occupied fibres use (P1.1); meeting one more
fibre replaces J/p by Q^A. These formulas need O(p+ell) polynomial terms,
with binary powering for their evaluations. This is an arithmetic-operation
description, not a polynomial bit-time claim for exponentially long outputs.

For r named binary attachment ports use Q=2^r. A color is the ENTIRE
r-bit answer word at a vertex. Bits within a word may be correlated;
no separate palette gauges or independent-port assumption is introduced.
For a full weight enumerator substitute x_i=sum_a w_a^i instead.

If the full parent tournament group equals H, these are exact counts for
that fixed parent with named pins and named new ports (with any separately
specified internal port tournament). If H is only a subgroup, they count
orbits in the retained arithmetic frame. Forgetting ports or the parent,
imposing primality, and moving by actual complements require additional
transport/realization tests. No A000568 term is added here.

## P3. Higher clock layers are free for the rigidity predicate

Put D0=A-B=N(p-1)/p^2. Color the unpinned vertices independently and
uniformly with Q>=2 colors, with any fixed assignment on S.

**Exact rigidity theorem.**

    Pr(Stab_G(coloring)=1)=(1-Q^(-D0))^(p-j).            (P3.1)

The number of rigid extensions is Q^(N-s) times this expression; the
number of their free G-orbits is that integer divided by |G|. The formula
does not depend on tau or ell. Those coordinates still affect all orbit
counts in P2 and the size of each free orbit.

### Proof: every extra symmetry has an order-p witness

Normalize one pin to zero. The root group is abelian: clocks commute
with fibre shifts. For f(z)=g^t z+eta h(rho(z)), h(0)=0,

    f^p(z)=g^(pt)z.

Thus f^p=1 iff t is a multiple of p^(k-1), which means f lies in F_0.
Intersecting with the pin stabilizer shows that its elements of order
dividing p are EXACTLY F_S, the independent socle shifts on unoccupied
fibres. This includes the identity and works at k=1.

Any nonidentity element of the finite p-group G has a nonidentity power
of order p. If it fixes a coloring, so does that power. A coloring has
trivial G-stabilizer iff it has trivial F_S-stabilizer. On disjoint
fibres, invariance under a tuple of translations is equivalent to
invariance under each of its nonzero component translations. All nonzero
translations in one fibre generate its same C_p subgroup. It therefore
suffices to test ONE generator on EACH unoccupied fibre.

Such a generator fixes Q^B of Q^A fibre colorings, a probability Q^-D0.
The events in different fibres are independent. Their complementary
probabilities multiply, proving (P3.1). QED.

More generally, for IID vertex colors with a named finite distribution
pi, replace Q^-D0 by

    (sum_a pi_a^p)^B.                                   (P3.2)

This allows arbitrary correlations among ports within a joint answer word.
Independence across vertices is essential for this probability formula;
the deterministic order-p test itself needs no probability model.

**Sharp sparse certificate.** Among binary extensions with all pinned
vertices white, exactly p-j black vertices are necessary and sufficient
to kill G. Necessity: an unoccupied fibre without a black vertex is
constant and supports a nonidentity shift. Sufficiency: choose one black
vertex in each such fibre; no nonidentity fibre translation preserves a
singleton, so P3 kills G. If j=p, G is already trivial and zero suffice.

At p=5,k=2, a single named root leaves |G|=3125 and

    Pr(rigid)=(1-2^-100)^4; minimum black support=4.

This is not a claim of unrooted full-tournament rigidity. Without the
root, H contains symmetries moving it, and a particular tournament may
have more automorphisms than H.

## P4. Why the smaller test is not a smaller universal carrier

For p=5,k=2, compare S={0,delta} (tau=1,ell=1) and
S'={0,5delta} (tau=3,ell=2). Both have two pins and j=1. Their rigidity
probabilities and sparse certificate costs agree exactly. Their group
orders are 625 and 3125 and their coloring-orbit counts differ, as the
exact receipt records. Hence dropping ell is valid for the Boolean
rigidity predicate, invalid for enumeration.

This is a precise way to obtain information without an extra measurement:
an order-p witness follows from any proposed higher symmetry. It does not
recover arbitrary joint port tensors or close under complement. The
81/27 hostile remains outside any claim that cyclic triples determine
the group: its 67 identical labelled cyclic triples still give 1620 versus
1728 regular binary masks. We independently replay that inherited control.

## P5. Exact correction for a surviving unsigned phase

This section joins P1--P3 to the incoming [Session112 axis/phase theorem,
S1--S2](../session112_counting_interfaces_20261030/symmetry/PROOFS.md), pinned
at `ecd0185d34`. Assume k>=2 and j>=2. Let W be the full unsigned-color
group, as inherited from U7. Put q=p^k, a=p-j, and let epsilon indicate
whether the actual named pins lie on a surviving reflection axis. Session112
proves that W_(S) is F_S if epsilon=0, and is

    F_S semidirect C_2, inversion action, if epsilon=1.

Each of the p^a elements of that reflection coset fixes exactly q vertices,
q/p in each residue fibre, and has order two. This is actual shared phase
data; j alone does not determine epsilon.

**Unsigned cycle index.** If epsilon=1,

    Z_(W_(S)) = [ Z_(H_(S)) + x_1^q x_2^((N-q)/2) ] / 2.    (P5.1)

If epsilon=0 the directed formula is unchanged. The same substitution and
Q^s division as P2 counts prescribed-color extensions. This follows by
summing the known cycles of each reflection, not by assuming the phase
acts freely on colorings.

**Unsigned rigidity probability.** For uniform independent Q colors,

    Pr(Stab_(W_(S))=1)
      =(1-Q^-D0)^a
       -epsilon p^a Q^(-(N-q)/2)(1-Q^(-D0/2))^a.           (P5.2)

**Proof.** Once a coloring is F_S-rigid, at most one reflection can fix
it: the product of two distinct such reflections would be a nonidentity
element of F_S fixing the same coloring. Thus the remaining bad events
are disjoint on the F_S-rigid locus.

Fix a reflection f. Its color-invariance probability is Q^(-(N-q)/2),
also after prescribing colors at the s fixed pins. In an unoccupied
fibre, f has (A+q/p)/2 cycles. The group generated by f and the fibre's
C_p shift is dihedral of order2p. Its identity fixes A points, nontrivial
rotations fix none, and each of its p reflections fixes q/p. Burnside
therefore gives (A+q)/(2p)=(B+q/p)/2 orbits in that fibre. Imposing the
shift in addition to f loses exactly (A-B)/2=D0/2 color choices in the
exponent. Conditional on f-invariance, the probability of failing the
shift test is Q^(-D0/2), independently in each free fibre. Multiply the
complementary probabilities, then sum over the p^a disjoint reflection
events. Subtract from P3 to obtain (P5.2). QED.

For an independent uniform ordered sample of r pins, repetitions allowed,
condition on j>=2 occupied residue fibres. Session112 gives
Pr(epsilon=1 | j)=p^(2k+j-2-kr). Consequently the joint pin/color failure
term in (P5.2) averages to

    p^(2k+p-2-kr) Q^(-(N-q)/2)(1-Q^(-D0/2))^(p-j).         (P5.3)

For a uniform minimal transversal (r=j=p), directed symmetry is already
gone. The unsigned-rigidity probability after a random Q-color attachment is

    1-p^(-(k-1)(p-2)) Q^(-(N-q)/2).                       (P5.4)

At625 and Q=2 the exact failure probability is 2^-300/125. This new joint
pin/attachment formula credits the incoming pin-only probability1/125;
it is not a fresh claim of that earlier calculation.

**Boundary control.** On two disjoint five-point cyclic fibres with one
COMMON reflection, the directed group has25 elements and900 rigid binary
words among1024. Adjoining that reflection gives50 elements and NO rigid
binary words, despite the orbit count remaining64. This literal abstract
action verifies the correction's necessity; it is not the full unsigned
ramified group at k=1. That shallow group has independent fibre reflections,
so P5 is explicitly restricted to k>=2. Our actual625 controls filter
250 two-pin maps, test all125 reflection cycle types, and independently
compute1000 generated-subgroup orbit counts for the proof's exclusion terms.
