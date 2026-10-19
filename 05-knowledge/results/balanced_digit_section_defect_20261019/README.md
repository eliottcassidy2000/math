# A balanced-digit mixed stress raises the 29-adic section-defect bound to 440

Date: 2026-10-19. Status: **PROVED**, with **FINITE-EXACT** controls.

## Result and scope

Let S be the integers of the unramified quadratic extension of Q_29, with its
nontrivial involution, and let G_r be the graph on S/29^r S whose edges have
norm-one differences. Every set-theoretic section of G_2 -> G_1 fails to
preserve at least **440 of the 12,615 edges** of G_1.

The previous straight-cycle argument gave 435. The improvement is strict and
analytic: it does not depend on a search, MILP infeasibility report, or a
large rigidity-matrix calculation. Sharpness at 440 is **not** claimed.
This is a statement about the edge defect of a prescribed section. It gives
no new chromatic bound, no Euclidean realization, and no implication that
chromatic number changes with precision.

The actual cross-family obstruction is a translated twelve-term stress on
three direction families. Its coefficients 1,3,9 cannot cancel for binary
second differences modulo 29. A putative equality support is thereby forced
to be three straight seams; a period rejects all eight possible orientations.
Rotation then amplifies this one-extra-error result to five extra errors.

## 1. Setup and inherited cycle bound

Use the quadratic norm model

    Q(x,y) = x^2 + xy + y^2.

Since -3 is nonsquare modulo 29, this identifies the residue field with
F_29[omega], omega^2-omega+1=0. The same polynomial gives the unramified
quadratic ring modulo 29^2. Write p=29, V=F_p^2, and use the representatives
X_x with both coordinates in {0,...,p-1}. For a direction d with Q(d)=1
modulo p, let

    ell_d(h) = (2d_x+d_y)h_x + (d_x+2d_y)h_y,
    D_d f(x) = f(x+d)-f(x),
    c_d(x) = (Q(X_(x+d)-X_x)-1)/p modulo p.

A section has the unique form X_x+p h_x, and its normalized edge error is

    z_d(x) = c_d(x) + ell_d(D_d h(x)).                 (1)

It preserves the edge exactly when z_d(x)=0. Reverse orientation does not
change the edge error; one representative from each opposite-direction pair
indexes every undirected edge exactly once.

Each direction pair partitions its p^2 edges into p disjoint straight
p-cycles. The sum of z_d around every such cycle is -2. Thus there are at
least p failures per direction family, and equality forces exactly one
failure on every cycle, with its error equal to -2. Consequently, equality
on three direction families has the form

    z_u=-2 I_u,  z_v=-2 I_v,  z_w=-2 I_w,             (2)

where each I_d is binary and has sum one on each d-cycle.

For completeness, the cycle sum follows directly from (1). Write
X_(x+d)-X_x=d-p k_d(x), and Q(d)=1+p K_d. Then

    c_d(x)=K_d-ell_d(k_d(x)) modulo p.

Around a straight cycle, sum k_d=d, so sum c_d=-ell_d(d)=-2 Q(d)=-2.
The motion contribution telescopes. This is the previous checkpoint's
mixed-characteristic carry obstruction, restated to make the proof independent
of any computational witness.

## 2. A compact mixed stress that also annihilates the carry

Choose

    u=(1,0),  v=(0,1),  w=(3,9)=3u+9v.

All three directions are norm one modulo 29, because

    Q(w)=9+27+81=117=1+4*29.

They are pairwise nonparallel. For every section, its errors satisfy

    D_u D_v z_w - 3 D_v D_w z_u - 9 D_u D_w z_v = 0. (3)

This identity is a twelve-term edge stress at each vertex. Here is a proof
which retains and cancels the carry rather than silently discarding it.

First, ell_w=3 ell_u+9 ell_v, and all translation differences commute. This
proves (3) for the rigidity-matrix term ell_d D_d h in (1).

Second, the integer wrap functions obey

    D_e k_d = D_d k_e.

Indeed, both sides equal
(X_(x+d)+X_(x+e)-X_x-X_(x+d+e))/p. Applying one further difference and using
c_d=K_d-ell_d k_d gives

    D_u D_v c_w = -ell_w D_u D_v k_w,
    D_v D_w c_u = -ell_u D_u D_v k_w,
    D_u D_w c_v = -ell_v D_u D_v k_w.

Their combination in (3) is zero. The constants K_d disappear under the
second differences. This proves the full affine identity.

This construction is valid more generally for any three unit directions
w=a u+b v, with the corresponding coefficients a,b. The special arithmetic
below is what makes 3,9 useful at 29.

## 3. Equality forces all three supports to be straight seams

Assume, for a contradiction, that a section has only 3p failures on the
three direction families. Substitute (2) into (3) and cancel -2. At each
vertex put

    A=D_u D_v I_w,  B=D_v D_w I_u,  C=D_u D_w I_v.

These are ordinary integers in {-2,-1,0,1,2}, and

    A-3B-9C = 0 modulo 29.

Its absolute value is at most 2(1+3+9)=26<29, so it is zero as an integer.
Reduction modulo 3 gives A=0; then B+3C=0 gives B=C=0. Thus

    D_u D_v I_w=0,  D_v D_w I_u=0,  D_u D_w I_v=0    (4)

as exact integer identities.

**Binary separability lemma.** If d,e are independent vectors in F_p^2 and
I:V->{0,1} satisfies D_d D_e I=0 over the integers, then I is invariant
under d or under e.

To prove this, use the d,e grid coordinates. The mixed difference identity
integrates to I(s,t)=F(s)+G(t), with integer-valued F,G. If both F and G were
nonconstant, the range of F+G would have width at least two, contrary to the
binary range of I. Thus one is constant.

Apply the lemma to (4). Since each I_d has one selected edge on every d-line,
invariance under an independent direction means its support is exactly one
affine line parallel to that direction. There are only the following choices:

- I_u is supported on a line parallel to v or w;
- I_v is supported on a line parallel to u or w;
- I_w is supported on a line parallel to u or v.

The line offsets remain arbitrary. The next obstruction is independent of
those offsets, so no search over them is needed.

## 4. A horizontal period rejects all eight orientations

It is useful to retain a,b briefly, then substitute a=3,b=9. Write

    eta_u=1 if the I_u seam is w-parallel, and 0 otherwise;
    eta_v=1 if the I_v seam is w-parallel, and 0 otherwise;
    eta_w=1 if the I_w seam is v-parallel, and 0 otherwise.

Let r_d=-2 I_d-c_d. Solving the u- and v-edge equations amounts to finding

    H_u=ell_u(h),  H_v=ell_v(h),
    D_u H_u=r_u,  D_v H_v=r_v.

Their cycle sums vanish, so they always have solutions. The two projections
determine h, since the matrix [[2,1],[1,2]] has determinant 3, nonzero modulo
29. General solutions add arbitrary f(y) to H_u and arbitrary g(x) to H_v.
The w-equation must therefore satisfy the necessary horizontal-period test

    sum_x [r_w-a D_w H_u-b D_w H_v](x,y) = 0.         (5)

Adding f or g changes this sum by zero: the f-term is multiplied by p and
the g-term telescopes in x. Thus any convenient integrated solutions may be
used in (5).

Here are the four sums needed for its evaluation.

1. Each I_u row has its selected edge at x=A_y. Cyclic summation by parts
   in D_u H_u=r_u gives

       sum_x H_u(x,y) = 2(A_y+1) modulo p.

   In fact sum H_u=-sum_x(x+1)r_u(x,y); the carry c_u is -2 only at x=p-1,
   whose weight x+1=p is zero. If eta_u=0, A_y is constant. If eta_u=1,
   its seam equation is b A_y-a y=constant, so
   A_(y+b)-A_y=a. Therefore

       sum_x D_w H_u(x,y) = 2a eta_u.

2. A w-parallel I_v seam has one selected point in each horizontal row;
   a u-parallel I_v seam has either zero or p, whose sum is zero modulo p.
   Summing the b consecutive v-differences gives

       sum_x D_w H_v(x,y) = -2b eta_v.

   The v-carry is -2 at y=p-1 and constant in x, so it contributes zero
   to each horizontal row sum.

3. The I_w row sum is eta_w modulo p. Consequently

       sum_x (-2 I_w)(x,y) = -2 eta_w.

4. For 0<a,b<p, the representative carry for w is

       c_w(x,y)=K_w-(2a+b) 1_(x>=p-a)
                         -(a+2b) 1_(y>=p-b) modulo p.

   Hence its horizontal sum is -a(2a+b). This formula is valid even when
   both coordinates wrap: the quadratic wrap term has an extra factor p.

Together, the left side of (5) is

    P = -2 eta_w+a(2a+b)-2a^2 eta_u+2b^2 eta_v.

For a=3,b=9 and p=29 this is

    P = 16-18 eta_u+17 eta_v-2 eta_w modulo 29.

In lexicographic order (eta_u,eta_v,eta_w) in {0,1}^3, its eight values are

    16, 14, 4, 2, 27, 25, 15, 13.

None is zero. Thus (5) fails for every possible seam orientation and every
choice of its offsets. This contradiction proves:

**Three-direction theorem.** Every section fails on at least 88 edges in
the direction families {u,v,w}.

The number 88 is a proved lower bound; no 88-error section for this subgraph
is asserted.

## 5. Rotational averaging gives 440

The 30 norm-one directions form a group under multiplication. Modulo sign,
there are 15 direction pairs. For each pair g, rotate the triple {u,v,w} by
g. The three-direction theorem applies to each rotated triple. This can be
seen by taking a norm-one lift of g modulo p^2, which exists by the elementary
norm Hensel step, and transporting the section by that rotation. Explicitly,
if N(G)=1+p k modulo p^2, then G(1-p k/2) has norm one modulo p^2. Equivalently,
repeat the preceding proof in the norm-preserving coordinates {g,omega g}.

Each direction pair occurs exactly three times among the 15 rotated triples:
once as gu, once as gv, and once as gw. If M is the total section defect,
summing their bounds yields

    3M >= 15*88 = 1320,
    M >= 440.

This argument does not require five disjoint triples or assume that their
failures occur at compatible positions.

## 6. General coefficient criterion (optional reusable corollary)

The proof identifies a reusable test, rather than a special-purpose rank
certificate. Let p>=5 be prime, and let O be the integers in the unramified
quadratic extension of Q_p, with its norm and conjugation. Form G_r on
O/p^r O using norm-one differences. In the anisotropic residue norm plane
O/pO, let u,v be independent norm-one directions, put
t=Tr(u conjugate(v))/2, and suppose w=a u+b v is another unit direction
with a,b nonzero. Thus

    a^2+2tab+b^2=1 in F_p.

Suppose that the only A,B,C in {-2,-1,0,1,2} satisfying

    A-aB-bC=0 in F_p

are A=B=C=0. Suppose also that none of

    tab,  1-tab,  a^2+tab,  b^2+tab

is 0 or 1. Then the same stress and binary separability argument applies.
The horizontal period requires

    eta_w = a^2(1-eta_u)+tab+b^2 eta_v,

which those four values exclude. Hence every section has at least 3p+1
failures on these three families, and norm-one rotational averaging gives

    total defect >= p(p+1)/2 + ceil((p+1)/6).

Only the p=29, t=1/2, a=3, b=9 instance is needed for the headline theorem.
For the general period derivation, choose lifts U,V of u,v with norm one
modulo p^2. They form a basis of O/p^2 O, and the norm in these coordinates
is X^2+2T XY+Y^2 modulo p^2, where T=t modulo p. Thus the same carry and
period calculation applies with mixed coefficient 2t. Its gradient matrix
is invertible because u,v form a basis in a nondegenerate norm plane.
The anisotropic residue norm has p+1 unit directions, so the rotational
average uses (p+1)/2 direction pairs, as required by the displayed bound.

## 7. Exact controls, exclusions, and attribution

Run

    python3 verify_balanced_stress.py

The script uses only the Python standard library and verifies:

- Q(3,9)=117 and its unit residue;
- all 125 possible small-digit triples, leaving only (0,0,0), with a
  modulus-23 aliasing hostile (2,-1,-2) showing why the range restriction matters;
- the rigidity-row cancellation for all 841 translated twelve-term stresses;
- the carry cancellation for every translate;
- an arbitrary deterministic motion as a positive control for the full
  affine error identity, and all 87 straight-cycle error sums;
- explicit u/v integration and the nonzero horizontal periods for all eight
  seam orientations;
- the offset-independence identities for all 8*29^3 offset/orientation tuples;
- the 15 rotations and exactly three occurrences of every direction pair.

The JSON and saved output record these controls. They supplement the proof;
there is no computer-only infeasibility assertion in the theorem.

Exploratory MILP computations are not used as evidence. In particular, a bounded
p=11 infeasibility experiment is not promoted to a theorem here.

Closest predecessors are
`unit_distance_carry_20261019/README.md` (the carry coset, cycle bound, and
q=3,5 positive witnesses) and
`unit_distance_stress_periods_20261019/README.md` (exact rigidity rank,
triangular-family periods, and the 6,716-dimensional cross-family gap).
The finite-field reduction programme also has the external antecedent
[Madore, The Hadwiger–Nelson problem over certain fields](https://arxiv.org/abs/1509.07023),
and the contemporaneous
[local-colouring research note](https://github.com/decalion89/chromatic-number-of-the-plane/blob/main/notes/local_colourings.md).
The latter's triangular decomposition and lifting context were source-checked;
this proof does not use its reported chromatic computations.

Rigidity linearization, commuting finite differences, and the binary additive
separability lemma are elementary standard ingredients. The contribution here
is their combination with the 1,3,9 modular separation and the offset-independent
period, giving the strict 435-to-440 improvement. No exhaustive literature
priority claim is made.

The global Hadwiger–Nelson status is deliberately outside this theorem. The
prior checkpoint source-checked OpenAI's 2026 reported 6<=chi(R^2)<=7 result
without independently rebuilding its formalization; nothing here depends on
that reported improvement or changes its status.

## Independent reproduction

Read the [independent review](review/REVIEW.md). Run `python3 -B review/independent_checks.py` and compare stdout and its generated JSON with the saved review receipts. Both verifiers use assertions: keep assertions enabled and do not use Python's `-O` flag. Both checkers require only Python 3 and its standard library. `SOURCE_PROVENANCE.json` records source, approved scope-clarification, and published hashes.
