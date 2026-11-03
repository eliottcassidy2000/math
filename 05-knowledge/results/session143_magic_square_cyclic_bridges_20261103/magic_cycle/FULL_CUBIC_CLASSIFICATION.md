# All-real cubic numerical relabelings of nine consecutive states

**Status:** FINITE-EXACT classification with a complete reduction to exact rational
linear algebra. The coefficient domain is all real numbers, not a bounded
integer or rational scan. Both ordinary and `python -O` runs pass, and Gaussian
ranks are independently recomputed from minors. No global originality claim.

## Statement

Let

\[
f(k)=Ak^3+Bk^2+Ck+D,\qquad k=-4,-3,\ldots,4,
\]

with real coefficients. Treat its nine values as an indexed multiset; repeated
values are allowed. They can be arranged as a real 3-by-3 magic square if and
only if

\[
\boxed{B=0\quad\text{and}\quad(A=0\ \text{or}\ C=-15A).}
\]

The constant coefficient D is arbitrary. Thus the only genuine cubic family is
an affine transformation of the values of `k^3-15*k`. This is a classification
of numerical labels, **not** a construction of square-valued labels. It does
not resolve the magic-square-of-squares problem.

## 1. Complete reduction to 945 homogeneous systems

Every real 3-by-3 magic square has its mean as its central entry, and every
opposite pair sums to twice that mean. The mean of the nine polynomial values is

\[
\mu=D+20B/3.
\]

Choose the input index j whose occurrence occupies the center. Then

\[
3Aj^3+(3j^2-20)B+3jC=0.
\]

Pair the eight remaining input occurrences according to their opposite cells.
For each pair {u,v},

\[
(u^3+v^3-2j^3)A+(u^2+v^2-2j^2)B+(u+v-2j)C=0.
\]

There are nine possible indices j and exactly `7*5*3*1=105` perfect matchings
for each: 945 systems, each consisting of five homogeneous linear equations
in (A,B,C). Repeated values cause no gap: their occurrences can still be
assigned distinct input indices, and every magic arrangement determines one
of these indexed matchings.

All matrix coefficients are integers. Rational Gaussian elimination computes
their entire real kernels, since the same elementary row operations are valid
over the reals. It is not an assumption that A, B, or C are rational. The
checker also computes rank independently by testing all 3-by-3 and 2-by-2
minors and nonzero entries.

The exact ranks are:

| Rank | Systems |
|---|---:|
| 1 | 1 |
| 2 | 30 |
| 3 | 914 |

The unique two-dimensional kernel is the plane B=0. All one-dimensional
kernels outside that plane are the following four coefficient rays, where
multiplication by an arbitrary nonzero real scalar is permitted:

\[
(1,3,-18),\ (1,-3,-18),\ (2,3,-39),\ (2,-3,-39).
\]

All other rays are already contained in B=0. Rank-three systems allow only
(A,B,C)=(0,0,0), the constant family. The JSON records every nontrivial kernel
space, its multiplicity, and an example matching, not merely this summary.

**Center distinction:** when B=0 the centered center *value* is zero. Its input
index need not be zero: an odd cubic may have additional zero values. The
reduction enumerates every j and never drops those branches.

## 2. Why the four nonodd rays fail the remaining magic constraints

For a zero-mean magic square, write the two top corners as u and v. Its entire
array is necessarily

\[
\begin{pmatrix}
u&-u-v&v\\
v-u&0&u-v\\
-v&u+v&-u
\end{pmatrix}.
\]

Consequently its four opposite-pair representatives, up to permutation and
individual signs, are

\[
u,\quad v,\quad u+v,\quad u-v.
\]

This condition is also sufficient by the displayed array. It remains valid
with repeated values and zero entries: equality of centrally symmetric
multisets is equivalent to equality of the absolute-value multisets with
multiplicity, including the zero-pair count.

The four nonodd coefficient rays have centered positive representative
multisets

\[
\{20,20,34,36\}\quad\text{or}\quad\{20,54,56,70\}.
\]

Neither can be of the required form. Indeed, if x=|u| and y=|v| with x>=y,
the absolute representatives are `{x,y,x+y,x-y}`. Their largest member is the
sum of two other members, so it is at least the sum of the two smallest.
But the two candidates satisfy `36<20+20` and `70<20+54`. The checker also
independently exhausts all 24 permutations and 16 sign choices for each ray.

## 3. Complete classification of the odd plane

For B=0 and A nonzero, subtract D and divide by A. Put t=C/A. The remaining
multiset is

\[
\{0,\pm p_1,\pm p_2,\pm p_3,\pm p_4\},\qquad p_i=i^3+ti.
\]

The signed representative criterion exhausts 24 permutations (a,b,c,d) of
(1,2,3,4) and 16 sign choices `(s_a,s_b,s_c,s_d)`. Each branch asks

\[
s_ap_a+s_bp_b=s_cp_c,\qquad s_ap_a-s_bp_b=s_dp_d.
\]

These are two linear equations in the one real variable t. The checker
handles inconsistent constant equations and parameter-free equations
explicitly, rather than dividing blindly by their coefficients. Exactly eight
of the 384 branches survive; all have `t=-15`. There is no parameter-free
surviving branch. This finite rational calculation therefore covers every real
t, including parameters producing zero or repeated p_i.

## 4. Sufficiency and actual consequence objects

The linear family uses the centered Lo Shu array

```
 3  -4   1
-2   0   2
-1   4  -3
```

scaled by C and translated by D. The genuine cubic family uses

```
  4  -22   18
 14    0  -14
-18   22   -4
```

scaled by A and translated by D. Its entries are exactly the values of
`k^3-15*k`, each once. Every one of its eight line totals is zero before
translation. The checker verifies these eight totals and both input multisets.

## Reproduction

```
python3 full_cubic_audit.py > full_cubic_audit.json
python3 -O full_cubic_audit.py > full_cubic_audit.optimized.json
cmp full_cubic_audit.json full_cubic_audit.optimized.json
```

`full_cubic_audit.py` uses only the Python standard library and explicit
exception gates. A separate independent parameter-elimination/graph-matching
check of the B!=0 sector is provided in `audit/full_cubic_nonodd.py`.
