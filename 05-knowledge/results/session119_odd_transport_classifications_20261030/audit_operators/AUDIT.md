# Independent audit of odd-function operator classifications

**Verdict: PASS, under the explicitly stated domains and coefficient conventions.**
This is an independent mathematical and direct-matrix audit of the signed-cycle
and reflection-odd multiplier classification. It imports no producer code. It
makes no priority claim and establishes no new general Ramsey or circuit bound.

## 1. Signed cycles and integral obstruction

Use (Qx)_i=s_i x_(i+1), with indices modulo L, s_i in {+1,-1}, and
eta=product s_i. Eliminate x_1,...,x_(L-1) using coefficients that are units in Z.
The remaining equation is

    (1-eta)x_0 = b_0+s_0 b_1+s_0 s_1 b_2+... .

Thus the integer Smith factors of I-Q are L-1 units and 1-eta, including the
zero case and L=1. The rational inverse in the eta=-1 case is the finite
geometric sum divided by two. Over Z, the final signed right-hand sum must
vanish when eta=1 and must be even when eta=-1. Modulo q its only remaining
obstruction is divisibility by gcd(q,1-eta), and a soluble fiber has exactly
that many points. These statements retain inhomogeneous compatibility, not
just kernel dimension. For I+P, put Q=-P; the product becomes (-1)^L times
the original cycle sign.

For a disjoint union of cycles, unit factors precede all factors two and all
zeros. This is already in Smith divisibility order. Fourier diagonalization
alone would not prove this integral statement.

## 2. Literal reflection-odd lattice and divisor strata

Fix an odd N>=3 and gcd(u,N)=1. On the integer lattice f(-x)=-f(x), f(0)=0
is forced. Use values at x=1,...,(N-1)/2 as independent coordinates. The row
for x in P_u sends it to ux modulo N, with a negative sign if the chosen
representative of ux is its negative. This is the matrix constructed by the
independent checker, without using divisor or order formulas to build it.

A residue of additive order d has an ordinary multiplication orbit of length
o_d=ord_d(u). Its reflection-pair orbit closes at the first exponent with
u^j=+1 or -1 modulo d. If -1 occurs, the length is o_d/2 and sign product -1;
otherwise the length is o_d and sign product +1. Since the stratum contains
phi(d)/2 reflection pairs, the cycle multiplicities are respectively
phi(d)/o_d and phi(d)/(2o_d). This proves the stated direct-sum decomposition,
including nonunit residues of the original modulus.

Applying Section 1 gives only Smith factors 1,2,0. In particular there is no
odd-prime cokernel torsion. Over an odd coefficient modulus, factors two become
units but zero factors remain genuine compatibility conditions.

**Coefficient-domain caveat checked.** Reducing this integer lattice modulo an
even q retains f(0)=0. The set of all functions into Z/q satisfying reflection
oddness additionally allows f(0) in the two-torsion subgroup. Both I-P_u and
I+P_u kill that extra coordinate. Its target must therefore be zero, and each
soluble fiber gains gcd(q,2) choices. The producer explicitly states this;
omitting it would have made its even-modulus count false.

## 3. Four-way arithmetic criterion

Put t_p=v_2(ord_p(u)) for p|N. Passing to p^a multiplies the order by an odd
power of p, so this valuation is unchanged. The condition -1 in <u> modulo
an odd divisor d is equivalent to all t_p, p|d, being equal and positive.
For sufficiency write each prime-power order as 2^t s_i with s_i odd and take
j=2^(t-1) lcm_i(s_i). Then j is congruent to half of each local order. For
necessity such an exponent j has v_2(j)=t_p-1 at every prime, forcing equality.
The unique element of order two in the cyclic local subgroup is -1 modulo
an odd prime power when its order is even.

Consequently:

- Common profile zero: every signed cycle is positive of odd length. I+P is
  nonsingular and I-P is singular.
- Common profile one: every signed cycle is negative of odd length. I-P is
  nonsingular and I+P is singular.
- Common profile t>=2: every signed cycle is negative of even length. Both
  operators are nonsingular.
- Nonuniform profile: choose primes p,q with distinct valuations. On the pq
  stratum -1 is not generated, while the order is even. Its positive even
  signed cycle makes both operators singular.

This covers all odd N>=3. It is a rational nonsingularity classification, not
a claim that a determinant power of two is a unit over the integers.

## 4. Full sign vectors require more than a kernel

For a signed cycle of length h and product epsilon, a vector with every
coordinate in {+1,-1} satisfying P sigma=alpha sigma exists exactly when
alpha^h=epsilon. If allowed, the initial sign is free and determines the whole
cycle, giving two choices. Requiring this on every divisor stratum proves:

- Preserving full odd signs exist exactly when every t_p=0.
- Reversing full odd signs exist exactly when every t_p=1.
- In either allowed case the number is 2^C, where C is the number of signed
  reflection-pair cycles. These are labelled functions, not isomorphism classes.

The prime strata establish necessity. The uniform-profile argument establishes
sufficiency on every composite stratum as well.

For N=15,u=2, the literal matrices have Smith factors

    I-P: (1,1,1,1,2,2,0),
    I+P: (1,1,1,1,2,0,0).

Both rational kernels are nonzero, yet neither admits a full sign vector. The
forced zero coordinates on incompatible strata are essential information, so
discarding them would falsely certify a tournament symmetry.

## 5. Independent finite verification

`check_independent.py` uses direct signed matrices in the literal residue-pair
basis. A general Euclidean integer row/column algorithm computes their Smith
forms; it is not a cycle-count implementation. Small matrices are cross-checked
again using gcds of every minor of each size, evaluated by exact Bareiss
elimination. Determinants, modular kernels, and sign vectors provide further
separate controls.

The exact universes are all signed cycles of lengths 1 through 7, all unit
multipliers for odd N=3,...,65, every determinantal minor when N<=11, every
vector modulo 2,3,4 when N<=9, and every full sign vector when N<=15. Normal
and optimized Python runs must agree byte for byte. The receipt records counts.
These finite controls corroborate the uniform proofs above; they are not a
proof-assistant build or independent human review.
