# Rank-88 Gaussian ratio: two minima, three stationary temperatures, and a complete branch-sign diagram

Status: PROVED with independent audit of the accompanying exact interval certificate. No minimum-shell enumeration, angular eigenvalue computation, or Gaussian-energy optimality claim is made. The universal scalar theta and fixed-angular-operator identities are inherited from the Session28 rank88 report and scalar follow-on.

## 1. New quantitative result

For every extremal even unimodular lattice of rank88, write its determinant-one Gaussian Hessian on Sym_0(R^88) as b(t)(r(t)I+A), with b(t)>0 and A independent of t. Then r has exactly three stationary points on (0,infinity):

* the strict global minima t0 and 1/t0;
* the strict local maximum t=1.

The exact rational interval certificate gives

1.098604471215707 < t0 < 1.098604471215708,

11135.8763212 < m := r(t0) < 11135.8763213.

The inherited exact identity and interval computation give M:=r(1) approximately11919.949215863706. The tail limit is L=2723200. On [1,t0], r is strictly decreasing; on [t0,infinity), r is strictly increasing toward L. Reciprocal symmetry determines the other half-axis. In particular, the scalar radial part a(t)=b(t)r(t) is strictly positive at every temperature.

These are computer-assisted analytic theorems, not samples on a finite temperature grid. The script makes every sign decision with integer/rational, outward-rounded interval arithmetic. Its floating-point display fields are not used as proof inputs.

## 2. Stationary equation with the cusp factor removed

Put x=2*pi*t, q=exp(-x), D=q*d/dq, E=E2, and

F=E4^11-2640*Delta*E4^8+1813680*Delta^2*E4^5-244992000*Delta^3*E4^2.

Write U=q^-4 DF, V=q^-4 D^2F, W=q^-4 D^3F and J=q^-4 Delta^4. Then

r=(V-45U/x)/(990J),

P(x,q)=x^2(-W+4EV)+45x(V-4EU)+45U,

and dr/dx=P/(990*x^2*J). Thus sign(dr/dt)=sign(P). Removing q^4 before numerical evaluation is crucial: direct interval evaluation of modular polynomials followed by division by Delta^4 creates artificial dependency/cancellation amplification.

If P=sum_{k>=0}(a_k*x^2+b_k*x+c_k)q^k, all a_k,b_k,c_k are exact integers. The constant coefficients are a0=b0=0, c0=30329640000. Successive ordinary x derivatives along q=exp(-x) are obtained by

(a,b,c) -> (-k*a, 2*a-k*b, b-k*c)

for each exponent k. This removes the need for numerical differentiation.

## 3. Exact tails and interval proof

The script generates F through degree84, hence P through degree80. Let R=1/100. On |q|=R, the elementary bounds sigma3(n)<=n^4<=16^(n-1), sigma5(n)<=n^6<=64^(n-1) give

|E4| <= 1+240R/(1-16R),
|E6| <= 1+504R/(1-64R).

Use |Delta| <= (|E4|^3+|E6|^2)/1728 in the displayed polynomial F to obtain M_F=1468058589186006771/18078415936. Cauchy's coefficient estimate gives |f_n|<=M_F*R^-n. Also sum |[q^n]E2| R^n <=1+24R/(1-4R)=5/4.

For n=k+4, the absolute P coefficient components satisfy

|a_k| <= M_F R^-n(n^3+5n^2),
|b_k| <=45 M_F R^-n(n^2+5n),
|c_k| <=45 M_F R^-n n.

On 0<=x<=16, three ordinary derivatives of the complete coefficient term are bounded by 20000 M_F R^-n n^6 q^k. For t>=1, q<1/500 is certified. Consecutive majorant terms after k80 have ratio less than1/4. Therefore the omitted tails of P and all derivatives used are less than10^-20. This intentionally loose common error is used symmetrically in every sign evaluation. Cauchy errors also dominate the smaller U,V numerator tails.

The interval engine uses signed integer endpoints with denominator2^180. Addition is exact; multiplication and inversion round outwards by integer floor/ceiling. Pi is enclosed by Machin's identity and alternating arctangent remainders. To enclose exp(-x), expand exp(x/16) through degree60, add a positive tail smaller than2^-250, invert, and raise to power16. The certified domain is 0<=x<=16. No library floating-point transcendental is used in proof decisions.

Centered mean-value interval enclosures reduce dependency overestimation: evaluate P (or P') at the midpoint, and bound its variation by the next exact derivative on the entire interval. Adaptive subdivision then certifies:

* P'<0 for 1<=t<=1.02;
* P<0 for 1.02<=t<=1.09;
* P'>0 for 1.09<=t<=1.12;
* P>0 for 1.12<=t<=2.2.

Here primes are ordinary derivatives with respect to x, not t; their sign is unchanged since x=2*pi*t. Reciprocal symmetry implies P=0 exactly at t=1. Exact point signs at1.09 and1.12 give a unique root in that window. Fifty-five exact bisections produce the displayed enclosure for t0.

For the whole remaining tail t>=2.2, let x0=4.4*pi. For every k>=1 and j=0,1,2, x^j exp(-kx) decreases for x>=x0. Evaluate the sum of absolute values of all nonconstant terms at x0 and add the Cauchy tail. The rigorous upper bound is20193127386, less than c0=30329640000. Hence P>0 on the entire infinite tail. This is not extrapolation from the last grid point.

Finally J=product_{n>=1}(1-q^n)^96 is enclosed using its first40 factors. The remaining product lies between1-192*q^41/(1-q) and1, by log(1-y)>=-y/(1-y), q<1/500, and exp(-u)>=1-u. Evaluating r on the certified root interval gives the displayed enclosure for m. Exact rational endpoints are retained in JSON.

## 4. Complete scalar branch-sign classification

For any angular eigenvalue lambda of A, put c=-lambda. The branch has the sign of r(t)-c. Write m=r(t0), M=r(1), L=2723200. Then:

* c<m: positive at every t.
* c=m: nonnegative, vanishing exactly at the two scalar minima.
* m<c<M: four simple zeros and two disjoint negative temperature intervals, exchanged by t ->1/t.
* c=M: two simple outer zeros and a double zero at1; negative between the outer zeros except for the isolated zero at1.
* M<c<L: two simple zeros and one negative central interval containing1.
* c>=L: negative for every finite t.

For an actual perfect lattice, lambda_min(A)>-L, so the last case cannot occur. Also, for any actual integral lattice in this family, the c=M case cannot occur: ARITHMETIC_ADDENDUM.md proves that A has algebraic eigenvalues whereas M is transcendental by Nesterenko's theorem. The abstract c=M case is retained above solely to describe the scalar function completely. The double-zero assertion at1 uses the inherited nonzero negative curvature; the minima are nondegenerate because P'>0 in their certified window.

Consequently all-temperature positive definiteness of the Hessian is equivalent to lambda_min(A)>-m. All-temperature positive semidefiniteness is equivalent to lambda_min(A)>=-m. These equivalences do not establish either spectral inequality for the constructed lattices and do not resolve higher-order stability at a zero Hessian eigenvalue.

The Morse index on the half-axis t>=1 first weakly increases toward t0 and then weakly decreases toward zero for a perfect lattice. Nullity can occur at no more than four temperatures for any one angular eigenvalue. This rules out arbitrary alternating temperature-instability phases despite the still-unknown angular spectrum.

## 5. Procedural exploration and hostile tests

Inherited mechanism: the quartic harmonic theta has only the Delta^4 sidecar, so one angular operator controls every temperature.
Canonical hostile test: the self-dual symmetry does not force a minimum; Session28 rigorously found a local maximum.
Corrected near miss: the old sparse-grid minimum-at1 guess was already superseded; no part of that numerical guess enters this proof.
Least-used sidecar: the scalar stationary numerator after cusp-factor removal, rather than the enormous minimum shell.

Concept board: reciprocal symmetry; fixed angular eigenspaces; exact q-coefficients; interval dependency; tail domination; representation multiplicities.

Candidate1 was direct modular-polynomial interval evaluation. Its natural interval dependencies cause impractical refinement: failure of the representation, not a mathematical counterexample. Candidate2 divided out q^4 and used centered derivative enclosures; the same target became a small finite certificate plus an analytic infinite tail. Candidate3 was automatic splitting using the Hermitian complex structure; it is not promoted: multiplication by sqrt(-23)/sqrt(23) is an orthogonal real complex structure but need not preserve the lattice or its minimum shell. Therefore its Hermitian/anti-Hermitian decomposition need not reduce A merely because the lattice admits a Hermitian tensor description.

## 6. Literature and remaining spectral obligation

Sakamoto's rank88 construction and its source-linked computational claims remain CITED, not independently recertified here. The paper gives a Hermitian tensor construction and automorphism information; these suggest reducing A using the actual finite-group action on Sym_0, with multiplicity spaces retained. They do not justify treating all complex-structure summands as A-invariant. Source: https://arxiv.org/html/2609.15633v2 .

Classical harmonic-theta/design provenance remains Nebe's survey https://arxiv.org/abs/1201.1834 and the Session28 sources. The new contribution in this packet is the exact scalar critical-point certificate and its branch-sign consequences, not a new general harmonic-theta theorem.

A further structural consequence is available for any extremal rank88 lattice, whether or not its minimum shell is perfect. Let T=A+L*I be the positive semidefinite shell frame operator on Sym_0. The minimum shell is perfect exactly when T is positive definite: the scalar identity direction is already detected by the shell, and the two-design identity makes it orthogonal to all traceless directions in the frame operator. If the shell is not perfect, some nonzero H has TH=0 and AH=-L*H. Since r(t)<L at every finite t, this one fixed H has strictly negative Gaussian Hessian in every temperature. Conversely perfection gives positive definite Hessians for all sufficiently large and sufficiently small temperatures, because r approaches L and lambda_min(A)>-L. Thus perfection is equivalent to eventual positive definiteness at either end; nonperfection gives a universal negative deformation direction, not merely an asymptotic failure.

The next decisive input is an exact bound on lambda_min(A), or enough orbit-compressed shell moments to supply one. Full shell enumeration of168498000 vectors was not attempted. Scalar theta alone cannot supply the missing angular operator.

Independent audit: /workspace/shared/session29-gaussian-audit/AUDIT.md accepts the exact scalar certificate and independently checks coefficient generation, derivative arrays, interval primitives, exponential enclosures, and Cauchy tail constants. This audit does not establish an angular spectrum.

## Reproduction

Run python3 certify_scalar.py. It uses the Python standard library only, takes approximately20-30seconds in the current environment, and writes certify_scalar.json. The .out file contains a human-readable transcript. All displayed decimal conclusions are backed by exact rational enclosures in JSON.
