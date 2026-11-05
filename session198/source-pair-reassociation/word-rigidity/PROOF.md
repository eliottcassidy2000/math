# Three singleton-shift tests determine the whole physical source ray

Session 198. **PROVED** all-order equivalence and exact growing defect. This
classifies singleton-word scalar reassociation; it does not settle arbitrary
source replacement or odd-target coverage.

Let v=(D,E,U,d) and w=(F,G,V,f) be actual nontrivial start-one sources. Let N be
the reviewed symmetric C2 bilinear form and let the physical singleton birth be

    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

For k>=0 define

    Delta_k(v,w)=N(v,J^k w)-N(J^k v,w).             (1)

Both terms are Hamilton-path counts of actual strong C2 grafts. They use the
same number of new vertices. Scalar symmetry alone gives Delta_0=0.

## 1. Complete singleton-word equality classification

**Theorem.** The following are equivalent:

1. Delta_k(v,w)=0 for k=1,2,3.
2. The full four-coordinate profiles are proportional:
   w=(F/D)v, with F/D a positive rational number.
3. Delta_k(v,w)=0 for every integer k>=0.

This is equality of numerical profiles up to scale, not graph isomorphism or a
claim that all proportional integer tuples are physical sources.

Proof. Put p=D+d and rho=D-E+U-d. Direct expansion of the reviewed N formula
gives

    Delta_1=rho(F+f)-p(F-G+V-f),                   (2)
    Delta_2-Delta_1=2[U(F+f)-pV],                 (3)
    Delta_3-2Delta_2+Delta_1=2[Df-dF].            (4)

These identities hold for arbitrary four-vectors. On actual sources D>0 and
p>0. If the first three defects vanish, (4) yields f=(d/D)F. Then (3) gives
V=(U/D)F, and substitution in (2) gives G=(E/D)F. Thus (1) implies (2).

Conversely, if w=c v, linearity of J and symmetry/bilinearity of N give

    N(v,J^k w)=c N(v,J^k v)=c N(J^k v,v)=N(J^k v,w)

for every k. This proves (2) implies (3), and (3) trivially implies (1). QED.

For an independent matrix check, the three linear rows of Delta_1,Delta_2,
Delta_3 acting on w are

    [-E+U-2d,  D+d, -D-d,   2D-E+U],
    [-E+3U-2d, D+d, -3D-3d, 2D-E+3U],
    [-E+5U-4d, D+d, -5D-5d, 4D-E+5U].

Their signed maximal minors are

    -4(D+d)^2 (D,E,U,d).

The rank is therefore three at every actual source, and its nullspace is
exactly that source's ray. The direct proof (2)-(4) does not need this rank
calculation.

## 2. The first exceptional fiber fails at every longer singleton word

Let t=(1,0,2,1), the transitive two-vertex source, and let

    w_d=(3d,d,3d,d)

be physical. Such states arise by applying M_2 to any actual balanced source
with H=U=d. In particular d=1 is the triangle, and unbounded examples have
d=3^j by repeated M_2. We do not assert physical existence for arbitrary odd d.

The one-letter defect vanishes: Delta_1(t,w_d)=0. But for every k>=2,

    Delta_k(t,w_d)=4d T_(k-2)>0,                   (5)

where T_0=T_1=1, T_2=2, and T_j=T_(j-1)+T_(j-2)+T_(j-3) for j>=3.

Proof. The transfer matrix A of J satisfies the exact identity

    A^4=A^3+A^2+A.

Thus for a_k=Delta_k(t,w_d), the sequence satisfies

    a_(k+4)=a_(k+3)+a_(k+2)+a_(k+1), k>=0.

Direct evaluation gives (a_0,a_1,a_2,a_3)=(0,0,4d,4d), then a_4=8d. Induction
yields (5). Every T_j is positive. No spectral approximation is used. QED.

At d=1,k=2, the actual eight-vertex C2 counts are

    N(t,J^2 w_1)=77,       N(J^2 t,w_1)=73.

Hence even the nontrivial one-letter equality fiber is not a returning rule
that can be iterated. More generally, the theorem in Section 1 says all-power
transport survives exactly on proportional full source profiles.

## 3. Scope and reproduction

The result excludes a scalar reassociation based only on matching the first
singleton-shift invariant. It leaves target-dependent source changes, mixed
gadget words, changes of gate, and all exact-target coverage questions open.
The classification is specific to N and the singleton source generator J.

Run `python3 verify.py` and `python3 -O verify.py`. The standalone verifier checks
all exact matrix identities on a basis, the rank minors at a bounded integer
universe, the three-test criterion on a complete bounded physical source bank,
the positive-defect recurrence on actual balanced-source rays, and the two
literal eight-vertex C2 witnesses. Finite checks supplement the algebraic proof.
The frozen BILINEAR_C2.md supplies only the physical source grammar, N graft,
and cycle law.
