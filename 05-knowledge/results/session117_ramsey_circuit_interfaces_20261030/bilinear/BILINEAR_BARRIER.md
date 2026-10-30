# Inner product: a classical Ramsey/Fourier counterbenchmark

Status: proved elementary benchmark, independently checked in small dimensions.
No novelty or improved explicit-NP circuit lower bound is claimed. The purpose is
to prevent a false transfer from strong spectral/rectangle properties to large
unrestricted Boolean circuits.

For m >= 1 let q(x,y) = x dot y over F_2, with x,y in F_2^m. The actual Boolean
input length is n = 2m. Circuit gates may be any binary Boolean function, with
unrestricted sharing and fanout. A rectangle here is A times B in the fixed
x/y partition; it is not a clique or an ordinary graph Ramsey witness.

## 1. Exact unrestricted circuit size

There are m AND gates computing x_i y_i, followed by m-1 XOR gates, so

    C_B2(q) = 2m - 1 = n - 1.

For the lower bound, all 2m inputs are essential: to witness x_i, set y_i=1
and every other coordinate zero, and conversely for y_i. In the ancestor DAG
of an output depending on n inputs and using g gates, the underlying graph is
connected. It has at least n+g-1 edges and at most 2g input-pin edges. Hence
g >= n-1. Constants or repeated pins do not weaken this count.

Thus the function is uniformly computable in linear time, and is in P, hence
NP with an empty witness. At odd lengths one may ignore the final bit to get
a total language; this is only an upper-bound benchmark, not a padding-based
lower-bound claim.

## 2. Perfectly flat Walsh spectrum and full sign-matrix rank

With unnormalized Walsh convention,

    W_q(a,b) = sum_(x,y) (-1)^(x dot y + a dot x + b dot y)
             = 2^m (-1)^(a dot b).

Indeed, summing over y gives zero unless x=b, and gives 2^m when x=b.
Consequently q is bent: every Walsh coefficient has magnitude 2^m. This is
standard bent-function terminology, not a new construction; see for example
[Nyberg, Multidimensional Walsh Transform and a Characterization of Bent
Functions](https://www.tcs.hut.fi/Publications/knyberg/itw07.pdf).
For H_xy = (-1)^(x dot y), character orthogonality gives HH^T = 2^m I,
so H has real rank 2^m. Even exponential rank in the half-input length and a
perfectly flat Walsh spectrum coexist with the minimum essential-input gate
count n-1.

## 3. Sharp monochromatic rectangle areas

For a nonempty zero rectangle A times B, put U=span(A), with dimension r.
Then B is contained in U-perp, so

    |A| |B| <= 2^r 2^(m-r) = 2^m.

Equality holds exactly when A=U and B=U-perp. Such rectangles exist for every
linear U, including U={0}.

For a nonempty one rectangle choose b_0 in B, and again U=span(A), dim U=r.
The functional x -> x dot b_0 is nonzero on U and has value one on A, so
|A| <= 2^(r-1). All points of B solve the same linear equations x dot y=1
for x in A. If consistent, their solution set is an affine coset of U-perp,
with size 2^(m-r). Therefore

    |A| |B| <= 2^(m-1).

This is sharp: take A={e_1} and B={y:y_1=1}. These are exact bounds for
rectangles, not a general circuit lower-bound method.

## 4. Complete maximal affine-flat statement

View q on V=F_2^(2m). Its polar form is

    B((x,y),(x',y')) = x dot y' + x' dot y.

It is nondegenerate alternating. On an affine flat a+U, constancy is equivalent
to

    q(u) = B(a,u) for every u in U.

In particular U is totally isotropic, and dim U <= m. Conversely, for any
isotropic U, q restricted to U is linear, and nondegeneracy lets that linear
functional be represented by B(a,-). Thus some constant cosets exist.

If dim U=m, then U=U-perp. The equations B(a,u)=q(u) have a unique solution
coset a+U. Every other coset has a nonconstant affine-linear restriction of q,
whose signed sum is zero. The global signed sum is W_q(0,0)=2^m. The unique
constant coset has 2^m points, so its value must be zero. Therefore every
Lagrangian direction has exactly one constant coset, always of color zero.

It follows that the largest zero affine flats have dimension m and the largest
one affine flats have dimension m-1. Sharp examples are y=0 (zero) and
x=e_1, y_1=1 (one).

Let L_m be the number of Lagrangian subspaces of V. Count pairs (U,v) with U
Lagrangian and 0 != v in U. Each U supplies 2^m-1 vectors. For each of the
2^(2m)-1 possible nonzero v, the Lagrangians containing v correspond to
Lagrangians of the nondegenerate symplectic quotient v-perp/<v>, of dimension
2m-2. Hence

    L_m (2^m-1) = (2^(2m)-1) L_(m-1),  L_0=1,
    L_m = product_(i=1..m) (2^i+1).

This is also the number of maximal-dimension zero affine flats, since a flat
has a unique direction. It does not count all inclusion-maximal smaller flats.

## 5. Reproducible checks and research implication

Run `python check_bilinear.py` and `python -O check_bilinear.py` from this
folder. Both were run and returned identical JSON, matching `receipt.json`.
The checker uses explicit exceptions, not disableable assertions. It checks
all nonempty rectangles, all linear subspaces and their affine cosets for
m=1,2,3; literally evaluates the displayed circuits and essential-input
witnesses; and checks every Walsh coefficient and sign-matrix orthogonality.
It covers 9, 225, 65025 rectangle pairs and 11, 307, 26387 affine flats,
respectively. Maximum zero-flat counts are 3, 15, 135. The all-m conclusions
come from the proofs above, not extrapolation from these finite checks.

Together with the separately proved point-defect family, this supplies two
opposite spectral controls: highly concentrated and perfectly flat spectra
both admit linear circuits. A proposed Ramsey-to-circuit transfer must use
additional structure that survives unrestricted gate sharing and tracks the
actual input coordinates. None of the statements here improve the known
linear constant for explicit families in NP, nor rule out more refined
Ramsey- or Fourier-based lower-bound methods.
