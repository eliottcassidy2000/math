# Shared-root fusion: a finite matrix carrier and a positive sorting switch

**PROVED from explicit path/cycle decompositions and the inherited odd-cycle
formula. FINITE-EXACT constructor controls accompany the proofs. Chan–Pak OPEN.**
This is a new physical operation and a finite interface for that operation;
no claim of completeness for arbitrary tournament contexts is made.

## 1. From old nonlinear carries to an actual pointed product

The archive sampling first suggested rapidity/group completion: distinguish
physical shared information from formally averaging two marked copies. The
productive continuation was the central-cross-term mechanism in
`01-canon/theorems/THM-3220-root-four-jet-schwarzian-heisenberg-transgression-and-oriented-discriminant-holonomy.md`
and `THM-3449-noncommuting-smooth-hensel-heisenberg-orbit-law.md`. Those concern
formal jets and Hensel automorphisms, not tournaments. Here no theorem is
transported between those objects. Instead we construct a literal operation
whose path count is exactly its own ordered bilinear carry.

The closest current mechanism is Session263's return join. There two distinct
marks are joined by a reversed edge. Here the marks are IDENTIFIED. The least
used sidecar is rooted path-length parity; its role becomes one matrix channel.
The hostile is forgetting the difference between an algebraic group commutator
and an actual positive word. Our matrices are generally singular and no inverse
operation on tournaments is asserted.

A pointed tournament has one distinguished vertex, labelled 0. For pointed
A,B let A⋆B identify their roots, retain each input's induced tournament, and
orient every edge from A−0 to B−0. This is a completely specified tournament
on |A|+|B|−1 vertices. It is associative, with the one-vertex rooted tournament
as identity: for three factors, all edges between nonroot factors go in their
factor order, regardless of brackets. The original distinguished vertex is
retained, unlike the moving-mark sign-flip grammar.

## 2. Eight integers close under every word

For a pointed T with root z write:

- H=H(T), h=H(T−z), with H(empty)=1.
- Q=Σ_(S⊆T−z) H(S)H((T−z)−S), the empty-slot two-path cover count.
- k=Σ_(S⊆T−z) (−1)^|S| H(S)H((T−z)−S).
- P⁺(w)=Σ_(paths R starting z) H(T−V(R))w^|R| and P⁻ analogously
  for paths ending z; U=P⁺(1), V=P⁻(1).
- d=−P⁺(−1), e=−P⁻(−1): odd-minus-even weighted rooted counts.

Define the integer matrix

    M(T) = [ h  U−h  d−h  H−h ]
           [ 0   Q    0   V−h ]
           [ 0   0    k   e−h ]
           [ 0   0    0    h  ].                         (1)

**Theorem F1.** M(A⋆B)=M(A)M(B). Thus this eight-integer interface compiles
arbitrarily long shared-root words with no increasing path-cover-order debt.

**Proof, paths.** Delete the root: (A⋆B)−z=(A−z)⇒(B−z). Hence h,Q,k multiply.
A root-start path either visits a nonroot A-vertex or does not. In the former
case it consists of a nontrivial rooted A-path followed by an arbitrary possibly
empty B−z path, with independent complementary paths. In the latter case it is
a rooted B-path and the whole A−z is in the complement. Therefore

    P⁺_(A⋆B)=(P⁺_A−w h_A)K_(B−z)(w)+h_A P⁺_B,
    P⁻_(A⋆B)=K_(A−z)(w)(P⁻_B−w h_B)+h_B P⁻_A.      (2)

Here K_C(w)=Σ_S H(S)H(C−S)w^|S|. Evaluating at ±1 proves all the noncorner
entries of (1), including the signs in d,e.

**Proof, cycles.** Apply the accepted odd-cycle packing formula
H(X)=Σ_F 2^|F|. A packing has at most one cycle through the common root.
Root-avoiding packings factor. Root cycles wholly in A or B contribute
(H_A−h_A)h_B and h_A(H_B−h_B). Every mixed root cycle consists uniquely of a
root-start A-path and a root-end B-path, both with at least two vertices,
sharing only their root. Its length is |p|+|q|−1, so it is odd exactly when
|p| and |q| have the same parity. All other cycles lie in the two complements.
If the odd/even nontrivial A-path masses are a_o,a_e and the B-path masses
are b_o,b_e, the mixed contribution is 2(a_o b_o+a_e b_e), namely

    (U_A−h_A)(V_B−h_B)+(d_A−h_A)(e_B−h_B).

Consequently

    H_(A⋆B)=H_A h_B+h_A H_B−h_A h_B
              +(U_A−h_A)(V_B−h_B)+(d_A−h_A)(e_B−h_B). (3)

This is exactly the remaining corner of the matrix product. No new cycles
or complementary packings are omitted. ∎

The empty-deletion case is essential: k=1 for the rooted singleton. For every
nonempty C, K_C(−1)=0. One proof expands by the same odd-cycle formula:

    K_C(w)=Σ_F 2^|F|(1+w)^(|C|−|V(F)|) ∏_(γ∈F)(1+w^|γ|).

Every summand for nonempty C has an odd factor vanishing at −1, since either
there is an unused vertex or at least one odd cycle. Thus k=0 for every
nontrivial pointed tournament. The parity channel is a last/first-factor
memory, not a group direction with an available inverse.

A singleton selected set has W=0 in the earlier selected-source carrier, so
its selected offset is D=H>0. This does NOT mean it has unit offset, nor that
h=H(T−z) is that offset. Fusion changes H and pays that change honestly.

## 3. An all-length closed family with a three-term spectral carrier

For a fixed nontrivial pointed T normalize

    ρ=Q/h, b=(U−h)/h, c=(V−h)/h,
    δ=(d−h)/h, η=(e−h)/h, x=(H−h)/h.

Since Q≥2h, ρ>1. Write T^[j] for j fused copies, j≥1. Theorem F1 gives

    H(T^[j])=h^j[1+jx+bc(ρ^j−jρ+j−1)/(ρ−1)^2+(j−1)δη]. (4)

Indeed b_j=b(1+ρ+...+ρ^(j−1)), c_j has the same formula, δ_j=δ, η_j=η,
and x_(j+1)=x_j+x+b_j c+δη. Summing this exact finite geometric progression
proves (4). It follows that the sequence from j=1 satisfies the recurrence
with characteristic polynomial (λ−Q)(λ−h)^2. Do not extend the formula
through j=0 when δη≠0: the identity has k=1 rather than 0.

For a directed triangle with any root,

    (H,h,Q,U,V,d,e)=(3,1,4,3,3,1,1),
    H(T^[j])=(4^(j+1)+6j+5)/9, order 2j+1.

The first five counts are 3,9,31,117,459. This special family also belongs to
old one-apex constructions, so its scalar sequence is an illustration, not
an asserted first discovery. The new object is the general carrier, valid
when root deletion has arbitrary internal cycles.

## 4. Context-uniform sorting from a new full-profile-neutral module trade

The `../returns/` lane constructs two SIX-vertex pointed modules A,B. At
root0 their exact data are

    H=29, h=9, Q=84, k=0, d=4, e=5,
    (U_A,V_A)=(52,41), (U_B,V_B)=(50,43).             (5)

Their literal outgoing-mask rows are A=(50,44,57,17,34,8) and
B=(42,52,57,18,33,8), each rooted at0. They agree in EVERY unmarked empty-slot cover count:
Q_r=r^6+12r^4+16r^2. Their ROOT-DELETED full cover profiles also agree:
Q_r(T−0)=r^5+6r^3+2r. Their labelled odd-cycle selectors differ, so (5) is
not inferred from that scalar profile; the rooted counts are retained.

**Theorem F2 (positive adjacent sorting gate).** For every pointed L,R,
including the rooted singleton identity,

    H(L⋆A⋆B⋆R)−H(L⋆B⋆A⋆R)
       =150[h_L V_R+U_L h_R−h_L h_R] >0.            (6)

**Proof.** Substitute (5) into (1). Direct integer multiplication gives
M_A M_B−M_B M_A=150N, with N_12=N_24=N_14=1 (one-based indices), all other
entries zero. Since the two orders have equal root-deletion h, their H
difference is the (1,4) entry of M_L(150N)M_R. This is exactly (6).
The root-only path contributes h to each of U,V, so U_L≥h_L and V_R≥h_R.
Thus the bracket is at least h_Lh_R>0. ∎

This is stronger than H(A⋆B)>H(B⋆A). The improvement cannot be reversed by
ANY surrounding word in this physical grammar, even with unrelated rooted
modules. No formal subtraction operation is implemented: subtraction only
compares two actual positive constructions.

**Corollary (fixed-inventory monotone count ladder).** With a copies of A and
b copies of B, all pointed words have order 1+5(a+b), and moving from
B^b A^a to A^a B^b by swapping each BA once produces ab+1 strictly increasing
actual Hamiltonian counts. The extreme sorted words minimize/maximize H over
all arrangements, because every word bubble-sorts monotonically. For a=b=t
this supplies t²+1 distinct counts with exactly the same inventory on 10t+1
vertices. This does not assert all binomial(a+b,a) words have distinct counts,
or that any consecutive odd integers occur. The gaps in this ladder are
explicitly positive, but can be large.

The relation between F1 and F2 is constructive: a cycle-address trade found
through the multi-return algebra becomes an order-robust arithmetic switch
only after its rooted response is passed to the shared-root carrier.

## 5. Connection contract, loss ledger, and reproduction

- **Old formal-jet / Hensel carrier → literal root-sharing carrier:** retain
  ordered bilinear cross terms and then the central scalar, rather than
  treating the quotient as abelian. Formal jets themselves are not mapped
  to tournaments. Our (1) is a singular semigroup representation, not a
  Heisenberg-group isomorphism and not an imported root-selector theorem.
- **Old product trichotomy → inspect the actual gluing:** ordinal join,
  shared-root fusion, and disjoint marked converse doubling are different
  operations. Equality of an input profile cannot identify their outputs.
- **Session263 parity filter → one extra signed matrix channel:** d,e,k retain
  exactly the parity needed for odd mixed cycles. Full rooted length data is
  lost, so this carrier is not adequate for arbitrary length-sensitive consumers.
- **New return trade → positive sorting:** the map is the eight-entry rooted
  profile, preserved under ⋆ by F1; discarded cycle addresses are unnecessary
  for THIS context. The decisive test is a common left/right fusion context,
  not an equality of unmarked H or Q_r.

Run `python3 -B 05-knowledge/results/session264_constructive_20261111/fusion/construct.py`.
The standalone counter counts actual directed paths by subset DP; it does not
call the matrix law to produce the literal expected profile. Initial controls:
all121 ordered pairs of pointed labelled tournaments with sizes1–3 (root0),
64 associativity triples from a stated four-object bank, and ten repeated-module
towers through five factors. Sixteen literal common contexts of the two six-vertex trade modules (through
15 vertices) check (6) with direct Hamiltonian subset DP. Four sorting ladders
through 41 vertices are exact matrix compositions, not full-graph recounts.
Larger matrix-only demonstrations are labelled as such. The parity-collision search in sizes1–5 found no hostile pair with the
chosen five equal scalar coordinates; no minimality claim is based on it.

Chan–Pak remains open. The new family gives physical composition with bounded
state and a useful monotone switch. It does not give a requested-integer
inverse cover, unit-offset preservation, or a height-tight all-target grammar.
