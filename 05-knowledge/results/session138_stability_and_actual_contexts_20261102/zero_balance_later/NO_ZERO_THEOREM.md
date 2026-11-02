# Exact sign classification of three actual mixed-profile modules

**Status: PROVED, all positive integer side multiplicities; accepted by the
independent Session138 reconciliation proof audit.** Fixed numeric profiles
are independently recounted in `audit.py`.
No bounded parameter scan is used as an all-order premise.

The baseline module is `M=Q4[C3,C3,C3,C3]`, where Q4 has masks `(8,1,3,6)`.
Delete vertex0 and restore at local cut867,995, or996, obtaining M_s. Local
labels are the remaining vertices1,...,11, followed by the new vertex.
The exact masks and ordered path-cover differences are frozen in
`mixed_profile_inputs.json`, copied from the first checkpoint's actual hostile
receipt with its source SHA256. No formal capacity perturbation occurs here.

Put `A_a=TT_a[C3,...,C3]`, for every integer a>=1, and define the actual
one-vertex response

    D_s(a,b) = H(C3[M_s,A_a,A_b])-H(C3[M,A_a,A_b]).

The original context is a C3-substitution of the strong quotient
`C3[Q4,TT_a,TT_b]`. The altered tournament changes only the deleted vertex's
arcs inside the M module; all its outer incidences stay canonical.

## Theorem

For s=867 and s=995,

    D_s(a,b)>0  if min(a,b)=1,
    D_s(a,b)<0  if a>=2 and b>=2.                        (1)

For s=996,

    D_996(a,b)>0 for every a,b>=1.                       (2)

In fact the last statement is stronger: for **any two nonempty tournaments**
A,B, not necessarily ordinal triangle sums,

    H(C3[M_996,A,B])-H(C3[M,A,B])
       >= 135618 H(A)H(B)>0.                            (3)

The coefficient135618 in(3) is sharp: equality holds for A=B=K1. The thin
positive statement in(1) also holds when one side is C3 and the other side
is any nonempty tournament, not merely an ordinal triangle sum.

There is a broader negative chamber. Write `p(m)=m^3+2m`. For s=867,995,
the response is strictly negative for any nonempty side tournaments A,B
whose colored path-cover polynomials satisfy

    P_A(m)=p(m)^2 R_A(m),   P_B(m)=p(m)^2 R_B(m),

where R_A,R_B are nonzero polynomials with nonnegative monomial coefficients.
In particular, each side can be an ordinal sum with at least two directed
triangles as components and arbitrary other tournament components. This is
a sufficient structural chamber, not a necessary classification of all sides.

Thus none of these three mixed-profile seeds gives a strong-quotient H-zero
completion anywhere in this unbounded two-parameter family. The sign reversal
for867/995 jumps across zero between side multiplicities1 and2; it does not
produce an intermediate integer zero. General strong-quotient zero rigidity
remains OPEN.

## 1. Exact tridiagonal response

Let f_W(r) count ordered covers of W by r nonempty directed paths, extended by
zero outside `1<=r<=|W|`. Write

    d_r=f_(M_s)(r)-f_M(r).

For arbitrary side modules A,B, directed run-words in an outer C3 are forced
cyclic sequences. Fixing r occurrences of M, the possible side multiplicities
are `(r,r)` (three possible starting phases), `(r-1,r-1)`, `(r+1,r+1)`, and
`(r-1,r),(r,r-1),(r+1,r),(r,r+1)` (one phase each). Therefore the response is

    sum_r d_r [3 f_A(r)f_B(r)
      +f_A(r-1)f_B(r-1)+f_A(r+1)f_B(r+1)
      +f_A(r-1)f_B(r)+f_A(r)f_B(r-1)
      +f_A(r+1)f_B(r)+f_A(r)f_B(r+1)]
       = f_A^T K_s f_B,                                 (4)

where K_s is the symmetric tridiagonal matrix

    (K_s)_(r,r)   = d_(r-1)+3d_r+d_(r+1),
    (K_s)_(r,r+1) = d_r+d_(r+1),                        (5)

and all other entries vanish. Indices start at1 and d_0=d_13=...=0.
These formulas count actual run words, including singletons, rather than
using an asymptotic coefficient approximation.

For s=867 the entries in rows1,2,3 are nonnegative, row4 has a positive
left neighbor and negative diagonal/right neighbor, and every entry in
rows5 and above is nonpositive. For s=995 the analogous mixed row is5:
rows1–4 are nonnegative and rows6 onward are nonpositive. The full exact
tridiagonals are in `audit.out`. Their positive prefix entries ensure that
multiplying K_s by any positive-support profile yields a vector with at most
one sign change, from positive to negative; the mixed row cannot introduce a
second sign change.

For s=996, K_s is entrywise nonnegative and `(K_s)_(1,1)=135618`. Equation(3)
follows immediately from (4), since `f_A(1)=H(A)>0` and likewise B.

## 2. A binomial-profile likelihood-ratio lemma

Let a polynomial P(m) have binomial expansion

    P(m)=sum_(r=1..N) v_r binom(m,r),

where the finite sequence `v_1,...,v_N` is positive and log-concave, extended
by zero outside its support. If R(m) is any nonzero polynomial with
nonnegative monomial coefficients, write

    P(m)R(m)=sum_r w_r binom(m,r).

Then every w_r is nonnegative, is strictly positive on1<=r<=N, and

    w_r/v_r is nondecreasing for 1<=r<=N.                (6)

There may be additional nonnegative coordinates r>N.

### Proof

Multiplication by m transforms a binomial coefficient sequence by

    (Lv)_r = r(v_r+v_(r-1)).                            (7)

The convolution `v_r+v_(r-1)` is log-concave: expanding its adjacent
log-concavity difference gives

    (v_r^2-v_(r-1)v_(r+1))
    +(v_(r-1)^2-v_(r-2)v_r)
    +(v_r v_(r-1)-v_(r-2)v_(r+1)),

all nonnegative by log-concavity and its monotone adjacent ratios. Boundary
zeros cause no problem. Multiplication by the positive log-concave sequence r
also preserves log-concavity. Thus L preserves log-concavity and expands the
positive support by one coordinate.

Moreover

    (Lv)_r/v_r = r[1+v_(r-1)/v_r]

is nondecreasing on the original support, because `v_(r-1)/v_r` is
nondecreasing. Iteration gives `L^e v / v` nondecreasing: it is the product
of the successive nonnegative nondecreasing ratios. Finally, express
`R(m)=sum_e c_e m^e`, c_e>=0. The ratio w/v is the nonnegative sum of these
nondecreasing ratios. This proves(6). No claim that an arbitrary such sum
remains log-concave is needed or made.

## 3. Apply the lemma to ordinal triangle side modules

The colored path-cover polynomial of A_a is

    P_(A_a)(m)=(m^3+2m)^a.

This follows from the colored odd-cycle formula, or from multiplicativity
under ordinal join. Its binomial coefficients are exactly f_(A_a)(r).
For a=2 the positive profile is

    v=(9,126,684,1656,1800,720),                         (8)

and it is log-concave (the exact adjacent minors are checked in `audit.py`).
For every a>=2,

    P_(A_a)=P_(A_2) (m^3+2m)^(a-2),

whose multiplier has nonnegative monomial coefficients. Hence

    f_(A_a)(r)/v_r is nondecreasing for r=1,...,6,        (9)

and all coordinates beyond6 are nonnegative.

## 4. Single-crossing comparison and the negative quadrant

We use the elementary following consequence of(9). Let h be a finite vector
whose signs are nonnegative then nonpositive, with its last positive
coordinate at some c<=6 and h_r<=0 for r>6. If

    sum_(r=1..6) v_r h_r < 0,

then, for every a>=2,

    sum_r f_(A_a)(r) h_r < 0.                           (10)

Choose `lambda=f_(A_a)(c)/v_c>0`. Below c the profile ratio is <=lambda,
so positive terms can only shrink relative to lambda times the base sum.
Above c it is >=lambda, so negative terms can only become more negative.
The additional coordinates beyond6 are nonpositive contributions. This
proves(10). If there are no positive coordinates, the conclusion follows
directly; in the applications below the first coordinate is positive.

The exact base responses are

    v^T K_867 v = -373209154588134,
    v^T K_995 v =  -93517858646046.                     (11)

The single-crossing row patterns in Section1 imply that `h=K_s v` satisfies
the sign hypothesis, with c at most4 for867 and at most5 for995. Thus(10)
first gives `f_(A_a)^T K_s v<0` for every a>=2.

Now fix any such a. By the same row-sign pattern,
`h=K_s f_(A_a)` again has one positive-to-negative crossing, and every
coordinate beyond6 is nonpositive. Symmetry of K_s and the preceding
inequality give `v^T h<0`. Apply(10) again, now with side parameter b, to
obtain `f_(A_b)^T K_s f_(A_a)<0` for every b>=2. This proves the second
line of(1). Exactly the same argument proves the broader negative chamber:
use the nonnegative multipliers R_A,R_B in the likelihood-ratio lemma in
place of `p^(a-2),p^(b-2)`. If a side is an ordinal sum containing two C3
components, multiplicativity of the colored cycle polynomial supplies this
factorization; every remaining tournament factor has nonnegative monomial
coefficients by the colored odd-cycle formula.

## 5. The thin boundary is strictly positive

If a=1, its side profile is `(3,6,6)` and zero beyond3. Consequently the
run-word coefficient of d_r in(4) vanishes for r>=5. For both867 and995,
all four values d_1,d_2,d_3,d_4 are strictly positive. At least the coefficient
of d_1 is positive, since it includes `3f_(A_1)(1)f_(A_b)(1)>0`. Thus the
response is strictly positive for every b>=1. Symmetry handles b=1.
The same proof uses only the nonnegative ordered profile and positive H of
the other side, so that side may be any nonempty tournament. Together with
Section4 this proves(1) in every positive integer case. For(3), singleton
sides have only profile coordinate1, so the stated lower bound is attained
exactly in that exhibited control; no claim of a full equality classification
is needed.

## 6. Exact audit and scope

`audit.py` independently computes each fixed twelve-vertex profile without
importing any Session137/138 counter: a backward endpoint DP counts every
induced Hamiltonian path; an anchored unordered-set-partition DP counts
path covers; multiplication by r! gives the ordered profiles. All three
vectors match the frozen inputs. It then checks every matrix entry/sign,
log-concavity of(8), the base values(11), the thin positive coefficients,
and selected finite likelihood-ratio controls. Only the exact fixed
identities, not the latter finite controls, are premises of the proof.

The new statement closes one explicit, unbounded construction family. It
neither proves that every mixed-profile seed has this row-sign pattern nor
that a general realizable side context preserves a one-crossing response.
The no-zero assertion for arbitrary strong quotients is still OPEN.
