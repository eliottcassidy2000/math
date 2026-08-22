---
id: THM-3403
title: "Hadamard-core maxdet, Smith, and circuit descent"
status: >
  PROVED (self-contained integer proof) + VERIFIED-EXACT (independent small
  positive, boundary, and hostile controls) + THM-3394 DEPENDENT COROLLARY.
  The classical comparison notes are CITED and carry no priority claim.
source: root-2608-crouzeix-puzzle-2026-08-15
audit: exact block algebra, determinant equality audit, parity-rank audit, Smith pairing audit, p-adic rank audit, and normal/-O byte-identical controls
depends_on:
  - THM-3394-twelve-formerly-missing-hadamard-orders-through-2000
dependency_scope: THM-3394 is used only for the twelve finite existence corollaries in Section 6
related:
  - THM-3393-hadamard-order-668-explicit-certificate
  - THM-447-skew-sylvester-doubling
  - THM-451-skew-tower-hadamard-chirality
  - THM-480-skew-tower-row-code-is-dplus
  - THM-481-paley-gauge-generates-gleason
script: 04-computation/hadamard_core_descent_thm3403.py
output: 05-knowledge/results/hadamard_core_descent_thm3403.out
script_sha256: 4f75f8be55840747790d29f54bbf9d202a445dcbaf09c09bca38f385d01c4894
output_sha256: 0a428ad558077b95f76724cd51f68d0a50a06b14a66400d50a4b37b1240e6322
dependency_thm3394_theorem_sha256: 63e3e841b609f2e9a38bbaa13dbbd665d6c05d7594fdb92351653207c3cadf86
dependency_thm3394_script_sha256: 7ae931b3cf268550287bd0621b9b85b8ea167126fadfb90d57b5106d0f82fb2d
dependency_thm3394_data_sha256: 68f7ceebb67005bf1b968171f7e6897cc33bde68adbd63f14bd45edfeb7b3f06
dependency_thm3394_output_sha256: d8efee90947015a7e6fc28a1685cc3d378357a85e1d4814953b32b17c5cd76a9
hash_basis: working-tree bytes (LF)
---

# THM-3403 -- Hadamard-core maxdet, Smith, and circuit descent

**PROVED (self-contained integer proof) + VERIFIED-EXACT (independent small
controls) + THM-3394 DEPENDENT COROLLARY.**

Let `m >= 1`, put

~~~text
N = 4m,                 v = N-1,
~~~

and let `H` be a Hadamard matrix of order `N`.  Normalize `H` by row and
column sign changes and write

~~~text
    [ 1   1^T ]
H = [           ],      B=(J-K)/2 in {0,1}^{v by v}.          (1)
    [ 1    K  ]
~~~

Here `1` is the all-one column and `J=11^T`.  No symmetry of `K` or `B` is
assumed.  For an integer `a`, the notation `a^r` in a Smith profile means
that `a` occurs with multiplicity `r`.

## 1. Universal core design and polar inverse

For every `m` and every normalized `H` as in (1),

~~~text
K1=K^T1=-1,                 KK^T=K^TK=4mI-J,                  (2)
B1=B^T1=2m1,                BB^T=B^TB=m(I+J),                 (3)
|det B|=2*m^(2m),           B^(-1)=(2B^T-J)/(2m)=-K^T/(2m).  (4)
~~~

Consequently `B` is an incidence matrix of a symmetric
`2-(4m-1,2m,m)` design: there are equally many points and blocks, every
row and column has weight `2m`, and two distinct rows meet in `m` points.
The word *symmetric* describes the design parameters, not the matrix `B`.

Moreover `B` is normal.  Its singular values are

~~~text
2m  once,              sqrt(m) with multiplicity 4m-2,       (5)
~~~

so its Euclidean operator condition number is `kappa_2(B)=2 sqrt(m)`.

### Proof

Orthogonality with the first normalized row and column gives the two sum
identities in (2).  The lower-right blocks of `HH^T=H^TH=4mI` give the two
Gram identities for `K`.  Expanding `(J-K)(J-K)^T/4` and its transpose,
using `J^2=vJ` and (2), gives (3).  Alternatively,

~~~text
BK^T=(JK^T-KK^T)/2=-2mI,
~~~

which proves the inverse in (4).  Finally `I+J` has eigenvalue `v+1=4m`
on `1` and eigenvalue `1` on `1^perp`.  Equations (3)--(5), including
normality, follow.  Taking determinants in (3) gives

~~~text
(det B)^2=m^v det(I+J)=m^(4m-1)(4m)=4*m^(4m),
~~~

which proves the determinant formula in (4).  This proof is entirely
internal; the design and determinant facts are classical context, not a
novelty claim.

## 2. Exact binary maximum determinant and the one-bit tariff

Define

~~~text
D_01(v)=max{|det A| : A in {0,1}^{v by v}}.
~~~

Whenever a Hadamard matrix of order `4m` exists,

~~~text
D_01(4m-1)=2*m^(2m).                                        (6)
~~~

This is a global optimum, not only a construction lower bound.  More
precisely, for arbitrary `A in {0,1}^{v by v}` set

~~~text
       [ 1    1^T ]
S(A) = [            ].
       [ 1   J-2A ]
~~~

Then

~~~text
det S(A)=(-2)^v det A,                                     (7)
|det A| <= 2*m^(2m),                                       (8)
~~~

and equality in (8) holds if and only if `S(A)` is Hadamard.

The inverse in (4) also determines the entire Hamming shell at distance
one.  Put `epsilon=sign(det B)`.  Then

~~~text
adj(B)=epsilon m^(2m-1)(2B^T-J)
      =-epsilon m^(2m-1)K^T.                              (9)
~~~

Thus every cofactor has absolute value `m^(2m-1)`.  If `B^(ij)` is obtained
by toggling any one entry of `B`, then, independently of `i,j` and of the
value toggled,

~~~text
|det B^(ij)|=(2m-1)m^(2m-1)
             =(1-1/(2m))|det B|.                         (10)
~~~

Every one-bit neighbor therefore pays the same exact determinant tariff
`m^(2m-1)`.

### Proof

Subtract the first row of `S(A)` from all lower rows.  The result is block
upper triangular with diagonal blocks `1` and `-2A`, proving (7).  Every
row of `S(A)` is a sign vector of Euclidean norm `sqrt(4m)`, so Hadamard's
determinant inequality gives

~~~text
2^v |det A|=|det S(A)| <= (4m)^(2m)=2^(4m)m^(2m).
~~~

This is (8).  Equality in Hadamard's inequality holds exactly when the sign
rows are pairwise orthogonal, which is exactly the assertion that `S(A)` is
Hadamard.  The matrix `B` from (1) attains the bound, proving (6).

Equation (9) is `adj(B)=det(B)B^(-1)`.  For a toggle let

~~~text
delta=1-2B_ij=K_ij,       B^(ij)=B+delta e_i e_j^T.
~~~

The `(j,i)` entry of `B^(-1)` is `-delta/(2m)`.  The exact rank-one
determinant lemma therefore gives

~~~text
det B^(ij)=det B (1+delta(B^(-1))_ji)
          =det B (1-1/(2m)),
~~~

which proves (10) and the stated absolute loss.

## 3. Odd-quarter parity: one full binary circuit

Assume now that `m` is odd.  Over `F_2`,

~~~text
rank_2(B)=v-1,
ker(B)=ker(B^T)=<1>,
row_2(B)=col_2(B)=E_v={x in F_2^v : sum_i x_i=0}.           (11)
~~~

Equivalently, both the row and column binary matroids of `B` are
`U_(v-1,v)`: there is exactly one nonzero dependence, and its support is
the full set of all `v` rows or all `v` columns.

If in addition `m>1`, then `H` has no closed quadruple of rows and no
closed quadruple of columns.  Here a closed quadruple means four distinct
sign rows or columns whose coordinatewise product is constant.  The case
`m=1` is sharp: in `H_4` the four rows and the four columns themselves are
closed quadruples.

### Proof

Reducing (3) modulo two gives

~~~text
BB^T=I+J.
~~~

Since `v=4m-1` is odd, `I+J` has kernel `<1>` and rank `v-1`.  Hence
`rank(B)>=v-1`.  But `B1=2m1=0` modulo two, so `rank(B)<=v-1`; equality
and both kernel statements follow.  Every row and column has even weight,
so its span lies in `E_v`; dimensions then prove the two span identities.
The kernels say exactly that the only nonzero dependence vector is `1`,
which is the full-circuit statement.

In the normalized matrix a constant row product must be all plus because
the first coordinate is plus.  If a closed row quadruple omits the first
all-plus row, its four core negative-support vectors sum to zero in
`F_2^v`; if it contains the first row, the other three core vectors sum to
zero.  For odd `m>1`, the unique dependency has support `v>=11`, so neither
a three-row nor a four-row dependency exists.  Apply the same argument to
`B^T` for columns.  When `m=1`, the unique core dependency has support
`v=3`, giving precisely the stated boundary case.

## 4. Squarefree quarter-order: exact Smith forms

Assume that `m` is squarefree.  Oddness is not required.  Every Hadamard
matrix of order `4m`, and every descended design from any normalization,
has

~~~text
SNF(H)=(1, 2^(2m-1), (2m)^(2m-1), 4m),                    (12)
SNF(B)=(1^(2m-1), m^(2m-1), 2m).                          (13)
~~~

Thus, for example, the Smith tier applies at order `8` (`m=2`) even though
the parity-circuit tier does not.

### Proof

Let `s_1|...|s_N` be the positive Smith invariant factors of `H`.  Since
the entries of `H` are units, `s_1=1`.  Since `H` modulo two is the all-one
matrix of rank one,

~~~text
s_i is even for every i>=2.                               (14)
~~~

We next prove, rather than cite, the Hadamard Smith pairing

~~~text
s_i s_(N+1-i)=N.                                          (15)
~~~

If `UHV=diag(s_i)` is a Smith reduction, then
`N diag(s_i)^(-1)` is integrally equivalent to
`NH^(-1)=H^T`; in particular each `s_i` divides `N`.  Reordering the
diagonal entries `N/s_i` into divisibility order shows that the invariant
factors of `H^T` are

~~~text
N/s_N | N/s_(N-1) | ... | N/s_1.
~~~

Transposition preserves Smith factors, so uniqueness proves (15).

At the central pair, divisibility and (15) give

~~~text
s_(2m)^2 | 4m.                                            (16)
~~~

Write `m=2^e u`, where `e` is `0` or `1` and `u` is odd squarefree.  No
odd prime can divide `s_(2m)` by (16), and its 2-adic valuation is at most
one.  Equation (14) therefore forces `s_(2m)=2`.  Divisibility gives
`s_2=...=s_(2m)=2`, and pairing (15) gives (12).

Finally, subtract the first row of (1) from every lower row and then the
first column from every other column.  These are unimodular integer
operations and give

~~~text
H ~_Z 1 direct_sum (-2B).                                (17)
~~~

If the Smith factors of `B` are `t_1|...|t_v`, the right side of (17) has
Smith factors `(1,2t_1,...,2t_v)`.  Comparing with (12) proves (13).

Squarefreeness is sufficient, not necessary.  It also cannot simply be
deleted from the theorem for arbitrary Hadamard matrices: the exact
companion's Sylvester matrix of order `16` has

~~~text
(1,2^4,4^6,8^4,16),
~~~

not the standard profile `(1,2^7,8^7,16)`.  Order `36` supplies odd
nonsquarefree branching in the classical literature; see the scope notes.

## 5. Simple prime divisors give self-dual residue codes

Let `p` be an odd prime with `p || m`, meaning that `p` divides `m` but
`p^2` does not; no squarefreeness hypothesis is placed on the other prime
divisors of `m`.  Then the row code

~~~text
C_p(H)=row_(F_p)(H)
~~~

is a Euclidean self-dual code with parameters

~~~text
C_p(H)=C_p(H)^perp,             [4m,2m]_p.                 (18)
~~~

### Proof

Modulo `p`, the identity `HH^T=4mI` makes the row space self-orthogonal;
hence its dimension `r` is at most `2m`.  On the other hand,

~~~text
|det H|=(4m)^(2m),             v_p(det H)=2m.              (19)
~~~

Exactly `4m-r` Smith factors of `H` are divisible by `p`.  Each contributes
at least one to the valuation in (19), so `4m-r<=2m`, or `r>=2m`.
Therefore `r=2m`, and a half-dimensional self-orthogonal subspace of the
nondegenerate Euclidean space `F_p^(4m)` is self-dual.  Notice that this
argument uses only `p || m`, not the complete Smith profile (12).

## 6. The twelve THM-3394 specializations

This section, and only this section, depends on THM-3394.  That theorem
proves Hadamard existence at

~~~text
N = 668,716,892,1132,1244,1388,1436,1676,1772,1916,1948,1964,
m = 167,179,223, 283, 311, 347, 359, 419, 443, 479, 487, 491. (20)
~~~

The second line is `N/4`, and every displayed `m` is an odd prime.  Hence,
for every pair `(N,m)` in (20), Sections 1--5 prove all of the following:

~~~text
D_01(N-1)=2*m^(2m),
SNF(B_N)=(1^(2m-1),m^(2m-1),2m),
SNF(H_N)=(1,2^(2m-1),(2m)^(2m-1),4m),
row_(F_m)(H_N) is a self-dual [4m,2m]_m code,
the core rows and columns each form one full binary circuit,
and H_N has no closed row or column quadruple.                 (21)
~~~

In particular, (6) determines the exact binary maximal determinants at

~~~text
v = 667,715,891,1131,1243,1387,1435,1675,1771,1915,1947,1963. (22)
~~~

The four checked-in THM-3394 dependency hashes are recorded in the
frontmatter.  The companion verifies those bytes before emitting (20)--(22)
but deliberately does not import or execute the THM-3394 renderer.

## 7. Exact companion and independent controls

The standard-library companion constructs Paley-I matrices independently at
orders `4,8,12,20` and a Sylvester matrix at order `16`.  It uses no floats,
randomness, solver, network, subprocess, dynamic import, `assert`, or file
write.  It checks:

- normalization and both Hadamard and design Gram identities;
- Bareiss determinants, the border determinant identity, polar inverse,
  adjugate, and the common absolute cofactor value;
- every one-bit toggle of every Paley control;
- binary ranks and all row and column closed quadruples;
- Smith factors by its own Euclidean integer row/column reducer; and
- the relevant odd-prime ranks and self-orthogonality checks.

The cheapest hypothesis controls are all present:

- `H_4` has one closed row and one closed column quadruple, pinning the
  `m>1` boundary;
- at order `8`, `rank_2(B)=3` rather than `v-1=6`, and there are `14`
  closed quadruples in each orientation, pinning oddness;
- orders `12` and `20` satisfy every applicable tier and every bit toggle
  has the predicted determinant; and
- Sylvester order `16` has the nonstandard Smith profile displayed after
  (17), pinning the squarefree hostile.

Reproduce from the repository root with

~~~bash
python3 04-computation/hadamard_core_descent_thm3403.py
python3 -O 04-computation/hadamard_core_descent_thm3403.py
~~~

Both runs are required to be byte-identical to
`05-knowledge/results/hadamard_core_descent_thm3403.out`.  These are finite
controls for the proof obligations above; they are not a computational
substitute for the universal proof and do not enumerate all binary matrices.

## 8. Classical comparison, THM-451 route, and scope

No priority claim is made for the classical Hadamard-design descent,
maximal-determinant bound, squarefree Smith theorem, or residue-code setting.
The proof package above is self-contained so that none of these comparisons
is a proof dependency.

- The standard Smith profile in the squarefree regime, the pairing (15),
  and the skew-Hadamard extension are surveyed in Peter Sin,
  [*The Smith normal form of Hadamard matrices*](https://people.clas.ufl.edu/sin/files/paley-hadamard.pdf).
- Hacıoğlu and Keman give a short proof of the standard profile for every
  skew-Hadamard matrix in
  [*A shorter proof of the Smith normal form of skew-Hadamard matrices and their designs*](https://dergipark.org.tr/en/pub/hujms/article/525333).
- Orrick's
  [*Switching operations for Hadamard matrices*](https://doi.org/10.1137/050641727)
  supplies the closed-quadruple/Hall-set context and records multiple Smith
  classes at order `36`, showing genuine odd nonsquarefree branching.

THM-451 now follows precisely this separate repair route.  THM-447 proves
that every level of the skew tower is skew-Hadamard, while the cited
skew-Hadamard Smith theorem applies at every level; MISTAKE-388 records why
its formerly conjectural status was wrong.  That cited all-level repair,
including nonsquarefree tower orders, is logically separate from the
self-contained squarefree theorem (12) here and does not become one of its
dependencies.

Finally, (6), (12), (13), (18), and (21) determine values, invariant factors,
and code existence.  They do **not** classify binary maximizers, Hadamard
equivalence classes, or the minimum distances or equivalence classes of the
resulting codes.  Nor does the twelve-order corollary extend THM-3394 beyond
its twelve finite existence certificates or prove the Hadamard conjecture.
