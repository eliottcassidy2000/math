# The common kernel, and a complete cube-plus-bridge grammar before modular wrap

2026-10-25. **PROVED** linear-algebra and physical-support statements below;
**FINITE-EXACT** independent controls supplement them. This audits and
strengthens the abstract-box theorem in [README.md](README.md). It changes
no norm-graph failure bound and makes no claim that a kernel vector is a
realizable binary error pattern or an actual section.

## 1. A universal kernel lemma

Let V be a vector space over any field, and let A,B,C be arbitrary linear
maps with domain V (their codomains can differ). Put

    W_A=ker B intersect ker C,
    W_B=ker A intersect ker C,
    W_C=ker A intersect ker B,
    J=ker A intersect ker B intersect ker C.

For nonzero s_A,s_B,s_C consider the weighted sum M on the external direct
sum of the three W spaces. Then

    ker M={(u,v,-(s_A*u+s_B*v)/s_C): u,v in J}.       (1)

In particular its nullity is exactly 2 dim J. **Proof:** if a weighted
triple sums to zero, apply A. Its last two terms vanish by their domain
conditions, so nonzero s_A forces A(f_A)=0. Thus f_A belongs to J.
Repeat with B,C. The remaining relation on three copies of J is precisely
(1). The reverse inclusion is immediate. No genericity, direct-sum
splitting or geometric assumption is used.

Applied to line-sum maps in three directions, every coupled relation is
therefore two coefficient channels of a common zero-line-sum array. New
coupled modes are never an extra unexplained interaction between the three
W spaces: they are exactly new modes of J, each counted twice. This does
not mean the old adjacent cubes always span J.

The abstract-box proof in README.md passes independently. Its coordinate
sum maps split using a sum-one basis vector in every characteristic; all
three W formulas are correct. There J=A0 tensor B0 tensor C0, of dimension
L-1, and (1) proves the claimed nullity. Its image dimension is then the
sum of the three domain dimensions minus 2(L-1).

## 2. Exact physical incidence before wrapping

Let p>3 be prime, u=(1,0), v=(0,1), w=(3,9) in F_p^2, and

    T_L={(-a-3b,j-9b): a,b in{0,1}, 0<=j<L},
    2<=L<=p-9.                                      (2)

The range can be empty for small p. In the stated range the four columns
have distinct x-coordinates and each contains L distinct points. The
physical u-line labels are j-9b and range over the integer interval
[-9,L-1], of length at most p. The w-line labels, after dividing their
cross products by 3, are j+3a, in [0,L+2], also without modular aliases.
Thus these integer line equations are exactly the physical F_p equations.
This is a direct incidence calculation; no abstract-box embedding is
asserted. In particular physical extra collinearities remain present.

Write an array on T_L as four length-L column arrays

    A_j=f_(0,0,j), B_j=f_(0,1,j),
    C_j=f_(1,0,j), D_j=f_(1,1,j).

Every array on these finite intervals is extended by zero to all integers.
Define

    h=A, P=A+B, R=A+C, S=B+D.

The physical u-line equations are R_j+S_(j+9)=0. Hence there is a sequence
r supported on 0,...,L-10 such that

    R_j=r_j, S_j=-r_(j-9),                           (3)

with r=0 if L<=9. The w-line equations are

    P_j+(R+S-P)_(j-3)=0.

Because all sequences have finite support, their unique solution is

    P_j=-r_(j-3)-r_(j-6)-r_(j-9).                   (4)

Substitution proves (4), while the difference of two solutions would be a
finitely supported 3-periodic sequence, necessarily zero. Formula (4)
has support inside 0,...,L-1 as required. Finally the v-line equations
are the four column sums. They are equivalent to

    sum h=0, sum r=0.                               (5)

Necessity follows from the A-column sum and the sum of R=A+C; sufficiency
uses sum P=-3 sum r. No division by three is involved in this last step.

Conversely every pair h,r satisfying the stated supports and (5) gives
one common-kernel array, uniquely, by

    A=h, B=P-h, C=R-h, D=S-P+h.                      (6)

These formulas prove both completeness and a constructive parametrization.
They work over any coefficient field for these line partitions; the
physical weighted coupling uses F_p and its nonzero (-3,-9,1).

### Theorem

Throughout (2),

    dim J = L-1 + max(0,L-10),
    dim ker M = 2(L-1) + 2 max(0,L-10).              (7)

The adjacent eight-point cubes supply the h differences. The additional
r differences supply translates of one twelve-point bridge. Together
these two types generate every common-kernel mode in this range, and
their two coefficient channels generate every coupled mode.

For p=29 this proves the complete kernel statement for every L=2,...,20,
extending the original finite test through L=12. It also explains why
individual W dimensions can change at L=5 while the coupled formula stays
unchanged through L=10: the missing coordinate is dim J, not any one
domain dimension. The bridge enters precisely when r has length two and
can support a nonzero zero-sum sequence, at L=11.

## 3. The first twelve-point bridge explicitly

At L=11 set h=0, r_0=1, r_1=-1. Its nonzero column entries are:

| Column (a,b) | Entries +1 at j | Entries -1 at j |
|---|---|---|
| (0,0) | none | none |
| (0,1) | 4,7,10 | 3,6,9 |
| (1,0) | 0 | 1 |
| (1,1) | 3,6 | 4,7 |

Each row here specifies actual points through (2). Every physical line
sum in the three directions is zero, as either direct pairing or (3)--(6)
shows. It is outside the span of adjacent cubes: for every adjacent cube
A+C=0, whereas this bridge has R=r nonzero. Its two coupled channels can
be written (F,0,3F) and (0,F,9F).

For general L, take successive r differences e_i-e_(i+1),
0<=i<L-10. They translate this twelve-point motif by i*v. Along with
the L-1 adjacent cubes, they are independent because (6) is a bijective
parametrization. Thus this is a structural repair of the extrapolation
failure, not merely a larger rank table.

The three shifted copies in (4) come from

    (1-t^9)=(1-t^3)(1+t^3+t^6).

This retains the actual directional relation w=3u+9v. A support count
alone forgets the shifts that create the twelve-point bridge.

## 4. Audit, controls, and remaining boundary

`independent_common.py` independently builds physical lines by testing
point differences for parallelism, forms the common line-sum matrix,
and computes a complete finite-field nullspace. It checks (6), the
explicit twelve-point bridge and its translates, independence and
spanning. The universe is all p in {23,29,31} and 2<=L<=p-9, with no
symmetry quotient or pruning. It does not import `check_box.py`.

The bridge is a hostile control against the earlier cube-only spanning
guess. The free h and r parameters distinguish a relation certificate
from a binary error pattern, which is not claimed. At L>=p-8 the integer
line labels can alias modulo p; the proof above no longer applies.
Updating those periodic line constraints is a precise next target, not
an automatic extension of (7).

Audit verdict for the original root packet: **PASS** for its abstract
theorem and its correctly labelled finite physical controls. The new
parametrization upgrades the scoped physical kernel formula from a
finite observation to a proof in (2). No mathematical correction to the
original packet was needed.
