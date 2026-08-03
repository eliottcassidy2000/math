---
id: THM-3248
title: "Q4 paired-owner Stirling compiler"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The Q4 walk resolvent factors through the paired coordinates
  Y=(1+X_1)(1+X_4)-1, Z=(1+X_2)(1+X_3)-1, and C=X_2X_3.  This gives an
  exact finite-difference/Stirling contraction for every path-cover depth of
  Q4[T_r,T_r,T_r,T_r], including a positive one-sum Hamiltonian formula,
  and an O_d(M(r)+r log r) fixed-depth unit-cost arithmetic compiler.
source: root/creative-reframes/2026-08-03
audit: >
  The assertion-independent exact companion derives the denominator and full
  walk numerator from the Q4 adjacency matrix, exhausts all 24 coordinate
  permutations, checks the first ten Hamiltonian values, compares the first
  four with a separate vertex-level Held--Karp engine, and compares eight
  all-depth cells with a separate permutation-and-cut path-cover engine.  It
  also checks 520 direct-versus-convolved finite differences and the
  determinant-only hostile.  Normal, optimized, and stored transcripts agree
  byte-for-byte; the script AST has no assert node or floating literal.  An
  independent hostile audit reconstructed D and N in SymPy, exhausted S4,
  separated the four-element observer group from the trivial tournament
  automorphism group, rederived the all-depth and Hamiltonian contractions,
  independently matched the first ten Hamiltonian values and Held--Karp
  controls, and used an unordered set-partition/path-cover DP to match all 24
  cells r=1,2,3 and 1<=d<=4r.  It also audited the fixed-d unit-cost bound and
  every stated nonimplication.
depends_on:
  - THM-3121-path-cover-walk-content-substitution-kernel
  - THM-3213-tournament-normalized-cyclic-diagonal-and-fast-moving-jet-transform
related:
  - THM-3226-unbalanced-q4-unequal-saddle-and-transcendence-wall
  - THM-3235-transitive-tournament-blowup-saddle-and-hypergeometric-decimation
  - THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence
script: 04-computation/tournament_q4_paired_owner_stirling_compiler_thm3248.py
output: 05-knowledge/results/tournament_q4_paired_owner_stirling_compiler_thm3248.out
script_sha256: 4a47d02a4412832c33b2f11c53deafb60fee4f2a9381cbe8c8f736b7e5513fc8
output_sha256: ca338bb41bf4be6848dd6f803f5d125d31e5cba6415241e965fc9364322bec6c
hash_basis: LF-normalized bytes
---

# THM-3248 -- Q4 paired-owner Stirling compiler

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem turns the unbalanced four-vertex quotient of THM-3226 into an
exact sequence compiler.  Its compression is not the determinant alone: the
full directed-walk numerator survives, but both rational-function pieces use
the same two coordinate pairs.  Transitive equal-block substitution then
turns those pairs into finite-difference rows.

## 1. Quotient, walk observer, and paired owners

Let `Q4` be the tournament on `{1,2,3,4}` with arcs

```text
1->2, 2->3, 3->1, 3->4, 4->1, 4->2.                    (1)
```

For `X=diag(X_1,X_2,X_3,X_4)`, use the row-to-column adjacency convention and
put

\[
 D(X)=\det(I-A_{Q4}X),\qquad
 W(X)={\bf1}^{T}X(I-A_{Q4}X)^{-1}{\bf1}=\frac{N(X)}{D(X)}. \tag{2}
\]

Define the three paired coordinates

\[
 Y=(1+X_1)(1+X_4)-1,qquad
 Z=(1+X_2)(1+X_3)-1,qquad C=X_2X_3.                    \tag{3}
\]

Then the complete quotient resolvent is

\[
 \boxed{D=1-CY,\qquad N=Y+Z+YZ+2CY.}                   \tag{4}
\]

Moreover, the simultaneous coordinate-owner group of `(D,N)` is exactly

\[
 \boxed{\Gamma=\langle(1\ 4),(2\ 3)\rangle\cong C_2\times C_2,} \tag{5}
\]

with owner orbits `{1,4}` and `{2,3}`.

Here `Gamma` is the coordinate symmetry group of the observer pair `(D,N)`.
It is not asserted to be the labelled tournament automorphism group or a
lossless invariant for arbitrary future operations.

### Proof

The matrix in the first determinant of (2) is

\[
 M=\begin{pmatrix}
 1&-X_2&0&0\\
 0&1&-X_3&0\\
 -X_1&0&1&-X_4\\
 -X_1&-X_2&0&1
 \end{pmatrix}.                                        \tag{6}
\]

The matrix determinant lemma, first in the rational-function field and then
as a polynomial identity after clearing `det M`, gives

\[
 N=\det\!\left(M+{\bf1}({\bf1}^{T}X)\right)-\det M.    \tag{7}
\]

The second matrix is

\[
 \begin{pmatrix}
 1+X_1&0&X_3&X_4\\
 X_1&1+X_2&0&X_4\\
 0&X_2&1+X_3&0\\
 0&0&X_3&1+X_4
 \end{pmatrix}.                                        \tag{8}
\]

Direct determinant expansion and regrouping by (3) give

\[
 \det M=1-CY,qquad
 \det\!\left(M+{\bf1}({\bf1}^{T}X)\right)
   =(1+Y)(1+Z)+CY.                                      \tag{9}
\]

Their difference is (4).  The two cubic monomial supports of `D-1` are
`{1,2,3}` and `{2,3,4}`; their intersection is `{2,3}`.  Hence every owner of
`D` preserves the unordered partition `{{1,4},{2,3}}`.  Conversely, the two
independent within-pair swaps visibly preserve both expressions in (4),
which proves (5).  \(\square\)

## 2. Equal transitive blowup and its exact diagonal

Let `T_r` be the transitive tournament of order `r`, and write

\[
 V_r=Q4[T_r,T_r,T_r,T_r].                              \tag{10}
\]

For a tournament `T`, let `pc_T(d)` count spanning covers by `d` unordered
nonempty directed paths, and set

\[
 F_T(d)=d!\,pc_T(d).                                    \tag{11}
\]

Thus `F_T(1)=H(T)`.  Put `u(x)=e^x-1`.  The transitive profile is

\[
 F_{T_r}(k)=k!\,{r\brace k}=r![x^r]u(x)^k.             \tag{12}
\]

THM-3121's walk-content substitution law, in the diagonal form of THM-3213,
therefore gives, for every `r,d>=1`,

\[
 \boxed{
 \frac{F_{V_r}(d)}{(r!)^4}
  =[x_1^rx_2^rx_3^rx_4^r]
    W\bigl(u(x_1),u(x_2),u(x_3),u(x_4)\bigr)^d.}        \tag{13}
\]

The paired coordinates linearize partly under this substitution:

\[
 Y=u(x_1+x_4),\qquad Z=u(x_2+x_3),\qquad
 C=u(x_2)u(x_3).                                       \tag{14}
\]

Equation (14), rather than determinant scalarization, is the mechanism behind
the compiler below.

## 3. The all-depth finite-difference compiler

For `n,k,ell>=0`, define

\[
 E_{n,k}(\ell)=\Delta^k(t^n)|_{t=\ell}
  =\sum_{j=0}^{k}(-1)^{k-j}{k\choose j}(\ell+j)^n.      \tag{15}
\]

Use the convention `E_(n,k)(ell)=0` for `k>n`.  For each fixed `d`, define
the finitely supported nonnegative integers `eta^(d)_(a,b,c)` by

\[
 (Y+Z+YZ+2CY)^d
   =\sum_{a,b,c\ge0}\eta^{(d)}_{a,b,c}Y^aZ^bC^c.       \tag{16}
\]

Then the complete ordered path-cover coordinate has the exact contraction

\[
 \boxed{
 \begin{aligned}
 F_{V_r}(d)
  ={}&\sum_{a,b,c\ge0}\eta^{(d)}_{a,b,c}
      \sum_{k\ge0}{d+k-1\choose k}E_{2r,k+a}(0)\\
    &\quad\cdot\sum_{\ell=0}^{b}(-1)^{b-\ell}{b\choose\ell}
       E_{r,k+c}(\ell)^2.
 \end{aligned}}                                        \tag{17}
\]

Every sum in (17) is finite: the outer support depends only on `d`, and the
finite differences kill `k+a>2r` or `k+c>r`.

### Proof

The elementary exponential identity

\[
 u(x)^K e^{\ell x}
   =\sum_{n\ge0}E_{n,K}(\ell)\frac{x^n}{n!}             \tag{18}
\]

gives the two paired diagonal extractions

\[
 [x_1^rx_4^r]Y^K=\frac{E_{2r,K}(0)}{(r!)^2},            \tag{19}
\]

and

\[
 [x_2^rx_3^r]C^KZ^b
  =\frac{1}{(r!)^2}\sum_{\ell=0}^{b}
       (-1)^{b-\ell}{b\choose\ell}E_{r,K}(\ell)^2.     \tag{20}
\]

Indeed, (20) follows by expanding
`Z^b=(e^(x_2+x_3)-1)^b` and applying (18) separately in the two variables.
Finally,

\[
 (1-CY)^{-d}=\sum_{k\ge0}{d+k-1\choose k}C^kY^k.       \tag{21}
\]

Insert (16) and (21) into (13), then apply (19)--(20).  Their product supplies
`(r!)^-4`, which cancels the normalization in (13), proving (17).  \(\square\)

## 4. Positive one-sum Hamiltonian formula

Write

\[
 A_k=k!{2r\brace k},\qquad B_k=k!{r\brace k},           \tag{22}
\]

with both rows extended by zero outside their natural ranges.  The `d=1`
specialization of (17) reduces to the positive one-dimensional sum

\[
 \boxed{
 H(V_r)=\sum_{k=0}^{r}\left(
   A_{k+1}B_k^2
  +2(A_k+A_{k+1})B_kB_{k+1}
  +(A_k+3A_{k+1})B_{k+1}^2\right).}                    \tag{23}
\]

To see the reduction, use

\[
 E_{n,k}(0)=k!{n\brace k},\qquad
 E_{r,k}(1)=B_k+B_{k+1},                               \tag{24}
\]

and collect the four numerator monomials `Y`, `Z`, `YZ`, and `2CY`.
The first ten values, beginning at `r=1`, are

```text
5
393
187049
310181025
1334332247705
12404417636753793
220149530887653870809
6819797079681288544615425
344577048040258820563492943705
26927377105708317796024318003822593.                   (25)
```

These are theorem consequences of (23), not recurrence guesses.

## 5. Fixed-depth arithmetic complexity

For `0<=K<=n`, the entire shifted finite-difference row in (15) is one
ordinary convolution:

\[
 \frac{E_{n,K}(\ell)}{K!}
  =[z^K]\left(\sum_{j=0}^{n}\frac{(\ell+j)^n}{j!}z^j\right)
          \left(\sum_{m=0}^{n}\frac{(-z)^m}{m!}\right). \tag{26}
\]

At fixed `d`, formula (17) needs the one row `(n,ell)=(2r,0)` and only the
`d+1` rows `(n,ell)=(r,0),...,(r,d)`.  Binary powering evaluates their input
values in `O_d(r log r)` arithmetic; the fixed number of products in (26)
costs `O_d(M(2r))`; and the contraction in (17) costs `O_d(r)`.  Hence, under
the usual constant-factor scaling convention for polynomial multiplication,

\[
 \boxed{F_{V_r}(d)\text{ is computable in }
   O_d(M(r)+r\log r)\text{ arithmetic and }O_d(r)\text{ storage}.} \tag{27}
\]

This is an exact unit-cost arithmetic bound over characteristic zero.  It is
not a bit-complexity bound, and the hidden constants are not uniform when `d`
grows with `r`.

## 6. Information boundary and sharp scope

The numerator in (4) is load-bearing.  At `r=1`, replacing `W=N/D` by the
denominator-only series gives

\[
 [x_1x_2x_3x_4]\frac{1}{1-u(x_2)u(x_3)u(x_1+x_4)}=1,   \tag{28}
\]

whereas the actual quotient has `H(Q4)=5`.  Thus the determinant and its owner
partition do not by themselves compute the sequence; the full walk observer
is the necessary sidecar.

The exact boundary is:

- (17) holds for every `r,d>=1`, but the fast bound (27) fixes `d`;
- the result covers equal transitive factors in this specific Q4 quotient,
  not arbitrary factor profiles or a generic unbalanced quotient;
- the moving `k`-range in (17) supplies no constant-order recurrence;
- no bit-complexity, C-finiteness, or P-recursiveness conclusion follows from
  the compiler; and
- as recorded in THM-3226, unconditional transcendence of the Q4 product
  radius and unconditional non-P-recursiveness remain **OPEN**.

## 7. Exact companion

Run

```bash
python3 04-computation/tournament_q4_paired_owner_stirling_compiler_thm3248.py
python3 -O 04-computation/tournament_q4_paired_owner_stirling_compiler_thm3248.py
```

The standard-library-only companion uses exact integers and rational numbers,
with no floats, randomness, or Python `assert`.  It independently reconstructs
`D` and `N` by permutation determinants and the matrix determinant lemma;
checks all 24 coordinate permutations; compares (23) with vertex-level
Held--Karp counts for `r=1,...,4`; compares (17) with full permutation-and-cut
path-cover enumeration for `r=1,2` and `d=1,...,4`; verifies (26) on 520
cells; and reproduces the determinant-only failure (28).  The normal,
optimized, and stored transcripts are byte-identical.

QED.
