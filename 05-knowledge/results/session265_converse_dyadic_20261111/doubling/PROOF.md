# A Fano deformation supplies a converse-hidden parity switch

**PROVED conditional only on the already-proved Session264 fusion law; the
seven-vertex input data and literal consequences through17 vertices are
FINITE-EXACT, independently counted here. Chan–Pak remains OPEN.**

This is constructive use of the old octonion/Fano material, not another audit
of it. It produces an arithmetic operation not present in the old150 switch:
the parity channel participates, and one fixed physical inventory already has
count-difference gcd2.

## 1. Inheritance and the actual 7 / 21 / 63 / powers-of-two map

The historical source
`05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md`
constructs octonion multiplication from a cyclic orientation of the seven Fano
lines, with seven nonidentity basis vectors in an eight-dimensional algebra.
There are21 pairs of distinct imaginary basis vectors; each line determines
three arc signs. Reversing any chosen lines gives a128-element binary family.
The inherited classification is16 alternative octonion tables, all Paley with
H189, and112 nonalternative line tables, all with H171. We do not claim a new
classification or identify a nonalternative deformation with standard
Cayley–Dickson octonions.

The current pre-audit `THM-5133-flip-and-spin-kernel-signature-sign-tournaments-of-twisted-algebras-and-the-pg32-extremal-theorems.md`
provides additional inspiration: quaternion lines, zero-divisor parity and
Cayley–Dickson doubling. None of its pending claims is needed for the switch
below. Its higher-dimensional doubling and zero-divisor statements are not
re-proved or promoted here.

Here is a precise additional dyadic/converse representation. Label the seven
line choices by x in F_2^7. Converse is x -> x+1, where1 is the all-one word.
For any real-valued function invariant under converse, Fourier expansion

    f(x) = Σ_s fhat(s) (−1)^(s·x)

has fhat(s)=0 whenever |s| is odd: pair x and x+1 in the Fourier sum. Thus the
128-dimensional function space splits into64 even and64 odd characters. The
converse-folded space has64 points, and **63 nonconstant character directions**.
This is an exact role for63: an augmentation-space dimension, not a forbidden
or attained Hamiltonian count. The roles of7 (line controls),21 (arc locations),
and63 (nonconstant even Fourier modes) are different and explicitly related
by a map, rather than equated because they are multiples of7.

H and every labelled undirected Hamiltonian edge weight W_uv lie in the even
sector. The root response U−V lies in the odd sector, since reversal exchanges
root-start and root-end families. Folding discards this whole sector. The
sidecar needed to compose pointed objects is signed endpoint/parity data.
For this fixed family the executable census gives only7 nonconstant Fourier
terms for H,21 for W_01, and13 for U_0−V_0. These support counts are finite
facts for the stated labels, not invariants of arbitrary tournaments.

The source-target contract is therefore:
- source: seven independent cyclic-line choices on the Fano plane;
- map: retain the tournament, then form W or rooted response polynomials;
- preserved: path reversal, and respectively even or odd Fourier character;
- lost by W: the root-start/root-end chirality;
- consumer: shared-root fusion, whose scalar output can detect that chirality.

## 2. Explicit labelled converse twins invisible to W

Use nonzero binary labels1,...,7. List Fano lines by a<b<a xor b and cyclically
orient (a,b,a xor b). Relabel original binary vertex3 as root0. The resulting
seven-vertex tournament A has out-neighbour bit masks

    A = (26,44,25,112,98,69,7).
    B = converse(A) = (100,81,98,7,13,26,56).

The bit in position j of entry i means i -> j. All input Fano-line
tournaments are strong: every pair belongs to a cyclic Fano triple, which
gives directed paths both ways between that pair. By path reversal A and B have
identical **labelled, rooted, undirected weighted Hamiltonian graphs**, not
merely the same sorted weight multiset. This follows for every converse pair,
independently of the finite numbers below. Their complete unmarked path-cover
polynomials and root-deletion path-cover polynomials also agree by reversal.

For a pointed tournament T define H,h,Q,k,U,V,d,e exactly as in Session264
`fusion/PROOF.md`: h=H(T−root), Q and k are the two-cover sum and alternating
sum on the deletion; U,V are the total weighted rooted path masses; d,e are
their odd-minus-even masses. The exact input table is

    (H,h,Q,k)_A = (H,h,Q,k)_B = (171,43,424,0)
    (U,V,d,e)_A = (260,257,24,25)
    (U,V,d,e)_B = (257,260,25,24).

For transparency, the root-start polynomials, coefficient indexed by length,
are

    P_A = 43w+39w²+37w³+35w⁴+37w⁵+44w⁶+25w⁷,
    P_B = 43w+37w²+37w³+35w⁴+37w⁵+44w⁶+24w⁷.

Consequently P_A−P_B=2w²+w⁷. In particular this trade changes both total mass
(by3) and parity mass (by−1). Its chirality is invisible to W even when the root
and every vertex label remain present. These input polynomials are checked by
a direct subset dynamic program; all5040 permutations independently verify
H and W. This is a finite certificate, not an unsupported symbolic assertion.

The inherited Session264 triangle-wrap endpoint consumer makes this packet
literal. For W(T)=u→T→v→u, old-root-start Hamiltonian paths correspond to a
root-start T-path R followed by v,u and a complementary T-path. The position
of v is |R|+1. Therefore W(A) has260 root-start paths versus257 for W(B),
with exactly2 extra paths at gate position3 and1 extra at gate position8.
Their odd-minus-even difference in the original root-leg length is−1
(and it is+1 in gate position, whose parity is shifted). This uses only the proved
endpoint consumer; it is not an unconditioned scalar projection of a chosen
length channel. The full response polynomials above are the certificate.

## 3. Exact positive scalar consumer

For pointed tournaments L,R, fusion L⋆R identifies their roots and sends every
nonroot L vertex to every nonroot R vertex. Session264 proves that its matrix
carrier multiplies:

    M(T) = [h,U−h,d−h,H−h; 0,Q,0,V−h; 0,0,k,e−h; 0,0,0,h].

It follows directly from the input table that

    M(A)M(B) − M(B)M(A)
      = [0,1143,43,1330; 0,0,0,1143; 0,0,0,43; 0,0,0,0].

Taking the corner after arbitrary physical left and right contexts gives

    Δ(L,R) = H(L⋆A⋆B⋆R)−H(L⋆B⋆A⋆R)
           =1143(h_L V_R+U_L h_R)
            +43(h_L e_R+d_L h_R)−1042h_L h_R.                 (1)

**Theorem.** Δ(L,R)>=1158 h_L h_R>0 for all pointed tournaments L,R.

**Proof.** Root singleton paths give U_L>=h_L and V_R>=h_R. Each parity mass
is a signed sum of the same nonnegative summands, so |d_L|<=U_L and
|e_R|<=V_R. Hence1143U_L+43d_L>=1100U_L>=1100h_L, and the analogous right
inequality holds. Substitute in(1), giving(2200−1042)h_Lh_R. Every tournament
has a Hamiltonian path, so h_L,h_R>=1, including the empty deletion convention.
This also covers one-vertex contexts. QED.

The estimate is deliberately elementary and need not be sharp. The important
point is that the negative parity entries of physical carrier matrices do not
invalidate positivity: the positive total-mass channel dominates them with an
explicit margin. Entrywise positivity of a commutator alone would not suffice.

A word with a copies of A and b copies of B therefore has a bubble-sort ladder
of ab+1 strictly increasing counts on1+6(a+b) vertices. Identical W-shadows of
the factors do not prevent their orders from affecting scalar H.

## 4. One fixed inventory reaches the primitive parity lattice

Let C be the rooted three-cycle. Its data are

    (H,h,Q,k,U,V,d,e)=(3,1,4,0,3,3,1,1).

First H(A⋆B)=60307 and H(B⋆A)=58977, giving1330. With C only on the left the
gap is3616: the counts are266425 and262809. This has gcd2 with the older150
switch, but that comparison alone uses different contexts/inventories.

The stronger statement uses **one common inventory {A,B,C,C}**, with17 literal
vertices after the shared-root identifications:

    H(C⋆A⋆B⋆C)=1191647,  H(C⋆B⋆A⋆C)=1185745,
    H(A⋆B⋆C⋆C)=1087199,  H(B⋆A⋆C⋆C)=1074439.

The two adjacent AB swap gaps are5902 and12760, and gcd(5902,12760)=2.
These follow by multiplication of the displayed small matrices and are also
verified by literal subset Hamiltonian counts of all four17-vertex outputs.

**Corollary.** The gcd of all pairwise H differences among arrangements of
this fixed inventory is exactly2. Every difference is even by Rédei parity,
and two differences already have gcd2. Therefore there is no forced affine
congruence modulus larger than2 on this arrangement set.

This is strictly an arithmetic-lattice result. The fixed inventory has at
most12 words, and the counts are far apart. Bézout combinations use subtraction
and are not physical constructions. No interval, semigroup saturation, or
cofinite odd-H theorem follows. The useful new ingredient for Chan–Pak is a
positive, context-uniform switch whose one fixed inventory is not trapped by
the150 congruence of the earlier two-module inventory.

This is not the first parity-only inventory in the unrestricted fusion grammar.
The rooted two-vertex source and rooted two-vertex sink already fuse in the two
orders to a transitive triangle and a directed triangle, a gap of2. Our specific
new carrier uses only **strong modules**, makes A and B indistinguishable by
the entire labelled rooted W-shadow and all unmarked/root-deleted cover
profiles, and retains a context-uniform positive swap despite parity channels.
Thus the primitive-lattice conclusion is new information about this particular
Fano/triangle family, not a universal arithmetic novelty claim.

## 5. Reproduction and boundary

Run:

    python3 05-knowledge/results/session265_converse_dyadic_20261111/doubling/construct.py

`receipt.json` contains all matrices, polynomials, labelled W and endpoint
matrices, literal counts, the full128-table Fourier supports, and121 ordered
pairs of contexts (all labelled tournaments of orders1,2,3). No claim of a
complete seven-vertex tournament census is made. New exact arithmetic is
checked; old published/repository claims are inherited with scope.

The future construction problem is to exploit the primitive lattice while
controlling the actual size of swaps and their compatibility. A signed Bézout
identity, a converse quotient, a Walsh projection, and a physical tournament
constructor are four distinct objects. This session supplies the last object
for a positive sorting switch, not a scalar projector for arbitrary Fourier
characters and not a Collatz descent proof.
