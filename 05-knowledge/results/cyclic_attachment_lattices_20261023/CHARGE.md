# Integer charge after modular gluing: restrict characters before taking totals

**PROVED, scoped.** Let S=r+H be a nonempty solution coset in
G=product_i Z/m_i obtained from the supplied attachment model. The theorem
below decides constancy of a specified integer observation on this SAME
coset. It neither minimizes an arbitrary cost nor proves a Collatz bound.

## 1. A small hostile to testing only group generators

In (Z/3)^3 let

    H={(x,y,z):x-y-z=0 mod3},
    u=(1,1,0), v=(1,0,1),
    P(x,y,z)=x-y-z,

where P uses the actual representatives0,1,2 as integers. H=<u,v>.
The integer cost vanishes on ALL of <u> and ALL of <v>. Nevertheless

    2u+v=(0,2,1),       P(0,2,1)=-3.

Thus even testing every multiple of every chosen generator does not
certify integer constancy. The modular zero becomes a nonzero carry.
The missing implication is treating P as a homomorphism to Z; a nonzero
homomorphism from a finite group to Z cannot exist. Generator tests remain
valid for actual group homomorphisms, not for integer representative maps.

This is a rejected trial shortcut, recorded before any publication.
The strongest survivor is an exact Fourier test, which retains mixed
characters rather than assuming additive integer behavior on the group.

## 2. General restriction theorem

Let a finite abelian group G act as the chart coordinate space. Suppose an
observation is supplied as a sum of local finite tables f_e on named
coordinate sets. Expand each table in the characters of its local group,
with the normalized finite Fourier convention

    f_e(x)=sum_chi fhat_e(chi) chi(x).

On the solution coset r+H, this becomes

    P(r+h)=sum_psi C_psi psi(h),
    C_psi=sum_(e,chi: chi|_H=psi) fhat_e(chi) chi(r).    (1)

**Theorem.** P is constant on r+H if and only if C_psi=0 for EVERY
nontrivial character psi of H. Its constant value, when this holds, is
C_1. For a vector observation, apply the same criterion coordinatewise.

**Proof.** Restriction and the coset shift give (1). Distinct characters
of H are linearly independent: a nontrivial character has zero sum because
translation by an element on which it is nonidentity rescales its sum.
Multiplying by the inverse of each character and summing isolates its
coefficient. This proves both directions, not only a sufficient bound.

The averaging in this proof is finite character orthogonality, not a
randomness hypothesis about actual sources or a density-to-pointwise
inference. Feasibility is checked first; the empty set is not a certificate
that an actual source has a constant value.

At most the sum of the supplied local table sizes contributes Fourier
terms before restriction. Shared-character terms must be ADDED before
testing for zero. No enumeration of all elements of H is required. This
bound is relative to supplied tables: a high-arity table may itself be huge,
and this does not manufacture small tables for Hamiltonian-path observables.

## 3. Explicit compact formula for linear integer charge

Suppose P(a)=sum_i w_i [a_i]_(m_i), using least nonnegative representatives.
Let d_i Z/m_i be the image of H in coordinate i and q_i=m_i/d_i. If B is
the integer lattice basis from the main theorem, then

    d_i=gcd(m_i,B_i1,...,B_it).

Write r_i=c_i+d_i beta_i with0<=c_i<d_i. On H write h_i=d_i t_i, where
t_i is a homomorphism H->Z/q_i. Put zeta_q=exp(2pi i/q). The constant
term from coordinate i is

    w_i [c_i+d_i(q_i-1)/2].                            (2)

For1<=k<q_i its character is psi_ik(h)=zeta_(q_i)^(k t_i(h)), with coefficient

    -w_i d_i zeta_(q_i)^(k beta_i)/(1-zeta_(q_i)^(-k)). (3)

To verify (3), change variables u=beta_i+t modulo q_i in the normalized
Fourier sum; sum_u u z^u=-q_i/(1-z) for z^q_i=1, z!=1. The constant part
c_i contributes only to (2). Therefore there are at most sum_i(q_i-1)
candidate nonconstant terms. This is linear in the total cyclic domain
sizes, even when |H| is their large product.

Characters are compared on the SAME complete generator bank: psi_ik has
the tuple of rational phases

    (k B_i1/m_i, ..., k B_it/m_i) mod1.                (4)

These tuples agree iff the characters agree on H, because the columns
generate H modulo the period lattice. Once equal tuples have been joined,
the coefficients are tested exactly in Q[z]/Phi_L(z), L=lcm_i(q_i).
The checker uses the finite Fourier sum rather than an inverse denominator
to obtain its coefficients, and performs exact rational polynomial reduction.
No approximate complex zero test is involved.

The formula includes q_i=1: that coordinate is fixed and contributes only
its actual value. Offset phases in (3) are essential; omitting r can change
the answer on different cosets of the same subgroup.

## 4. A family with exponentially many lifts and exact charge cancellation

Take2r named consecutive cyclic27-cells, r>=2, joined in a ring. On
alternate edges use the matching block (rotation equality modulo27); on
the other edges use the quotient block (equality modulo9). With cell
indices0,...,2r-1, the solution equations are

    a_(2j)=a_(2j+1) mod27,
    a_(2j+1)=a_(2j+2) mod9,                           (4a)

cyclically. Choose the common residue modulo9 in nine ways. Each of the
r pairs has three independent lifts to modulo27. Hence the full count
is9*3^r. Fixing one actual owner gives3^(r-1); fixing one owner per pair
gives one. Every allowed frame has the exact integer charge

    P=sum_j([a_(2j)]_27-[a_(2j+1)]_27)=0.             (4b)

This has a direct proof because paired coordinates are equal even as
least representatives. It also provides a large positive control of the
general character implementation. At r=32 the checker obtains
16677181699666569 frames from a64-by-64 lattice and tests1664 candidate
terms, combined into584 distinct characters, all with zero coefficient.
The common modulo9 quotient causes additional character coincidences
between different pairs. This is why characters are grouped on the full
shared generator bank rather than paired only by visual resemblance.

The family is an authenticated supplied-cell tournament construction:
orient all remaining cross-blocks uniformly, and retain the named cells
as part of the input. There are54r vertices. The computation operates on
the proved block records; it does not enumerate all full adjacency matrices
or claim that the named cells are recoverable from unmarked degree data.

## 5. The source-specific zero-slot charge is a positive control

The earlier [zero-slot proof](../odd_cell_gluing_20261023/ZERO-SLOT.md)
uses V_0=(1,0), V_1=(0,1), V_k=V_(k-1)+V_(k-2). The first register is
geometric length. Its three units A,B,C have second charges0,1,2. A
preferred proper split loses one z=(0,1); changing C to B loses one more.
For a fixed ordinary Zeckendorf root representing n,

    Y(n)=leaf charge + internal split count + recolour count,
    Y(n)=2n-floor((n+1)/phi^2).                         (5)

These identities and the coherent row family are inherited proved results.
What this session adds is their precise interface with (1).

Give every recolourable C occurrence its actual absolute position j and
bit t_j. Its leaf-charge contribution changes by -t_j; its debt changes
by +t_j. On the SAME choice space these two terms have identical restricted
characters and opposite coefficients. They cancel identically in (1),
at every finite prefix and every refinement depth. The original integer n,
the row's split tree, and its occurrence addresses remain fixed. A new
split adds exactly its declared debt; it does not reset the source register.

If leaf and debt use independent copies t_j and u_j, their nontrivial
characters differ and do not cancel. The finite control with six actual
owners has64 legal joint choices and constant total charge. Erasing their
identifications creates4096 choices and a varying total. This is the
same precise ownership failure as separate cell-attachment choices, now
in a different finite choice group. It does not identify binary colour
choices with odd-order tournament automorphisms.

The map is: typed expansion choices -> leaf/debt coordinate functions ->
restricted characters on their shared choice space. It preserves (5),
forgets a particular expansion only after its contribution is accounted
for, and needs the occurrence-address sidecar. This gives refinement-stable
conservation, not a uniform bound on the accumulating debt or positive
Collatz drift. Those are different inequalities still to be proved.

## 6. A central-kernel version of the same discipline

Let1->K->E->Q->1 be a central extension. Suppose actual edge transports
g_ij in E have a flat quotient connection: their images equal q_j q_i^-1.
Choose one lift u_i of each q_i and retain

    k_ij=u_j^-1 g_ij u_i in K.

Full compatible vertex frames of the form u_i z_i exist exactly when
k_ij=z_j z_i^-1. On a connected graph this is equivalent to trivial K
product on every fundamental cycle; tree propagation proves both
directions. Changing the chosen lifts changes k by a K-coboundary and
does not change the cycle obstruction. Eliminating private path interiors
multiplies the retained K charges; shared interior vertices stay shared.

This statement requires centrality and a flat quotient. It does not replace
general nonabelian lifting by an untwisted additive equation.

For an odd-order hostile take the upper-unitriangular group over F3,

    (a,b,c)(a',b',c')=(a+a',b+b',c+c'+a*b').

Its quotient by the c-coordinate is F3^2. For X=(1,0,0), Y=(0,1,0),
the commutator XYX^-1Y^-1=(0,0,1) is nonidentity although its quotient
is zero. A four-edge return can therefore close in the quotient and fail
in the full group. The analogous quaternion example has central charge-1
over the quotient C2^2. Neither finite group is asserted here to be the
full automorphism group of a particular tournament.

The literature recovery in [NONABELIAN.md](NONABELIAN.md) now gives the
stronger boundary: explicit finite group-coset feasibility is already a
tractable Mal'tsev problem, and full solutions have a compact prefix bank
relative to the listed groups. The remaining challenge is native or implicit
input compression and the actual observation retained through it. Central
cycle charges are one explicit form of that retained information; they
are not a necessary restriction on the general feasibility theorem.
