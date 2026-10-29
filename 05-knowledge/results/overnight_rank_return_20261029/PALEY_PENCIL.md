# The full integral prime-Paley pencil and its ramified correction

2026-10-29 heartbeat. **PROVED elementary / FINITE-EXACT controls.**
Same-agent proof audit and independent modular elimination; no external
priority claim. This extends [ATOM_PRECISION A2](ATOM_PRECISION.md), whose
original coprimality and unramified hypotheses were sufficient and correct.
The full Smith form below removes both restrictions and handles singular
and signed weights. The vertex order must remain prime.

## P4. Complete formula

Let q be a prime congruent to 3 modulo 4, A the q-vertex Paley adjacency
matrix, r=(q-1)/2, and E=aA+bA^T for any integers a,b.
If a=b=0, every Smith factor is zero. Otherwise put

\[
g=\gcd(a,b)>0,\quad x=a/g,\quad y=b/g,\quad
s=x+y,\quad\delta=x-y,\quad N=(s^2+q\delta^2)/4.
\]

If s is nonzero, put R=|rs| and

\[
\boxed{d={\gcd(R,N)\over\gcd(q,s)}}.
\]

Then the q nonnegative Smith factors, in divisibility order, are

\[
\boxed{g\big(1^{[r]},\ d,\ N^{[r-1]},\ RN/d\big).} \tag{P4}
\]

Superscripts in square brackets denote multiplicities. If s=0, the answer is

\[
\boxed{g\big(1^{[r+1]},\ q^{[r-1]},\ 0\big).} \tag{P5}
\]

Thus this is an unconditional algebraic formula for every actual
prime-Paley Hamilton endpoint matrix once its two endpoint counts a,b
are supplied. It is not a formula computing those Hamilton counts.

## Proof, including the integer extension at q

Scaling every matrix entry by g scales every Smith factor by g, so assume
g=1. Write Q=A-A^T and J for the all-ones matrix. The character identity
Q^2=J-qI and QJ=JQ=0 give

\[
E={s\over2}(J-I)+{\delta\over2}Q,\qquad
E^2+sE+NI=cJ,\quad c=(qs^2+\delta^2)/4. \tag{P6}
\]

The fractions are identities over Q; N,c are integers. The characteristic
polynomial is (X-rs)(X^2+sX+N)^r, as inherited from Session54.

For a prime ell different from q, the projector J/q splits the Z_ell
lattice into its constant line and augmentation lattice. The constant
contribution is Z_ell/(rs). If ell divides N, then s is an ell-adic unit:
for odd ell, simultaneous divisibility forces ell to divide both x and y;
for ell=2, even s makes primitive x,y both odd, and then N is odd.
The two roots 0 and -s modulo ell lift to distinct roots of
X^2+sX+N. Their polynomial projectors split the augmentation lattice into
two rank-r free modules. One eigenvalue is a unit and the other has
valuation v_ell(N). This gives (Z_ell/N)^r, exactly as in A2. If q does
not divide s, there is no q-primary contribution.

Now suppose q divides s, and put t=v_q(s)>=1. Primitivity makes delta a
q-adic unit. Therefore v_q(N)=1 and c is a unit in Z_q. Modulo q, E is a
nonzero scalar multiple of Q. The q-rank of Q is r+1: its (i,j) entry is
(j-i)^r in F_q. Expanding this polynomial gives r+1 nonzero binomial
coefficients and two full-rank evaluation matrices for distinct monomials
of degree at most r. Hence E has exactly r nonunit q-primary factors.

In the cokernel, (P6) says N times the module is the image of cJ, a cyclic
submodule. As N is q times a unit, at most one of those r factors can
have q-exponent greater than one. The determinant has valuation r+t.
Thus the q-exponent list is

\[
0^{[r+1]},\quad1^{[r-1]},\quad t+1. \tag{P7}
\]

The uncorrected direct-sum guess would instead give 0^[r],1^[r],t.
Removing precisely one factor q from gcd(R,N), and moving it to the last
factor, produces (P7) while leaving every other primary component intact.
This proves (P4), including d|N and N|RN/d.

If s=0, primitivity gives E=Q or -Q. The cokernel has free rank one.
For ell different from q, the augmentation restriction is invertible
because Q^2=-qI there, so its torsion is q-primary. The sum of coordinates
descends to a surjection from the cokernel to Z. On its torsion subgroup,
the identity Q^2=J-qI implies qz=0. The rank r+1 modulo q leaves r-1
torsion factors besides the free factor, proving (P5).

## P6. The cheapest hostile and its precision law

For q=3 and (a,b)=(1,2), the correct Smith form is (1,1,9), not (1,3,3).
The rational spectral factors and determinant agree with both candidate
groups, but the actual cokernel is Z/9, not (Z/3)^2. The obstruction is
the unavailable integral projector J/3.

For primitive weights in the ramified case, the exact response exponent is

\[
\log_q|\operatorname{im}(E\bmod q^k)|
=(r+1)k+(r-1)(k-1)_+ +(k-t-1)_+. \tag{P8}
\]

The uncorrected splitting undercounts the image by a factor q at precisions
1 through t, and agrees only after that. Equal determinant is not equal
response memory. Common content g shifts the corresponding prime digits
by v_q(g); it must be removed before applying the primitive formula.

## P7. Field degree changes the visible rank

The prime assumption has a sharp mechanism, not just a missing proof.
Let q=p^m=3 modulo 4 and use the actual Paley tournament over F_q.
The character matrix now satisfies

\[
\boxed{\operatorname{rank}_{F_p}(A-A^T)=((p+1)/2)^m.} \tag{P9}
\]

Proof: Q_ij=(j-i)^((q-1)/2) over F_q. The base-p digits of (q-1)/2
are all (p-1)/2. In characteristic p,

    (1+X)^((q-1)/2) = product_(h=0)^(m-1) (1+X^(p^h))^((p-1)/2).

There are exactly ((p+1)/2)^m nonzero coefficients. The full monomial
evaluation matrix in degrees 0,...,q-1 is invertible, so expansion of
(j-i)^((q-1)/2) gives exactly that rank. Its entries lie in F_p, so
extending to F_q does not change rank. If p divides s but not both
primitive weights, E has the same rank modulo p.

At q=27, the actual rank is 8, whereas substituting 27 into the prime-only
rule would predict 14. Over F3[z]/(z^3-z-1), exact Smith forms include

```text
A+2 A^T:  1^[8], 3^[12], 9^[6], 351
A-A^T:    1^[8], 3^[6], 9^[6], 27^[6], 0.
```

The digit factors are the concrete field-degree sidecar. Frobenius rotates
their field coordinates; forgetting those factors destroys the rank law.
This is the binomial-coefficient Lucas mechanism, not a claim about Lucas
prime distribution or an identification of different objects named Lucas.

## Inheritance, literature and verification

The closest mechanism is the local splitting proof in A2. The hostile is
the nonsplit Z/9 above; the corrected near miss is using a rational
isotypic decomposition as an integral direct sum. The underused sidecar
is the augmentation lattice in [Session35](../session35_dimensions_dualities_20261023/recursions/AUGMENTATION-LATTICE.md),
which already warns that rational equivalence loses integral divisors.

For context, Chandler--Sin--Xiang's primary
[Paley graph paper](https://people.clas.ufl.edu/sin/files/paley4.pdf),
Sections 2--4, treats undirected q=1 modulo 4 graphs. It separates
noncharacteristic primes from the characteristic-primary calculation,
using Jacobi-sum valuations and a transfer matrix for digit carries.
That is a useful next method for prime-power tournaments; its graph
theorem is not imported as a theorem about this directed pencil.

The current script checks 556 signed/scaled pencils at q=3,7,11,19,23,31:
86 ramified, 48 singular, 416 unramified and six zero matrices. All have
exact Smith forms; an independent unit-pivot algorithm over Z/q^k checks
every q-primary layer. It also constructs F27, checks every nonzero field
inverse, and certifies three prime-power hostiles.

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/paley_pencil.py --write
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/paley_pencil.py
```

**Subsequent closure.** [PRIME_POWER Q1--Q7](PRIME_POWER.md) now proves
the complete prime-power law using the cited classical Jacobi valuation
and cyclic carry recurrence, with actual F27,F243,F343 checks. The prime
formula here is its m=1 specialization; the failure of substituting a
prime power directly into (P4) remains a valid hostile. These pencil
formulas classify linear responses, not the entire complement carrier.
