# Prime-power Paley pencils from a cyclic carry recurrence

2026-10-29, 09:00 UTC continuation. **PROVED using a CITED classical
Jacobi-valuation theorem / FINITE-EXACT independent checks.** Same-agent
audit. No literature-priority claim: the carry generating function is
classical; the work here derives its complete directed-pencil application,
retains the exceptional source block, and verifies the actual matrices.

This closes the prime-power Smith obligation left in
[PALEY_PENCIL](PALEY_PENCIL.md). It does not count the Hamilton paths
needed as input weights, classify every tournament, or close the entire
complement/re-marking carrier.

## Q1. Inheritance and the precise scope

Let p be a prime congruent to3 modulo4, m a positive odd integer,
q=p^m, and A the actual Paley tournament adjacency on F_q. Fix integer
weights a,b and the pencil E=aA+bA^T. The vertex relation uses field
subtraction and quadratic residues, not subtraction modulo the integer q.

The closest proved mechanism is the prime-pencil splitting theorem.
The canonical hostile is the F27 character rank8, where the prime-only
formula gives14. The corrected near miss is using rational eigenvalues
without an integral character basis. The least-used sidecar is the
fixed-class lifting filtration in the independently audited
[October19 synthesis](../automata_arithmetic_geometry_session_20261019/SYNTHESIS_REVIEW.md).

The anchor is the full prime-power pencil; the niche is integral lifting;
the wildcard is a finite carry automaton that counts Smith layers without
building a q-by-q matrix. The concept board is character source / cyclic
carry / field phase / integral block / precision / realization guard.

**Classical input.** Chandler--Sin--Xiang,
[The Smith and critical groups of Paley graphs](https://people.clas.ufl.edu/sin/files/paley4.pdf),
Theorem3.4, records the Stickelberger Jacobi-sum valuation. Their Section4
also supplies the carry generating denominator used below. Their graph
theorem has q=1 modulo4; we rederive the directed character blocks for
q=3 modulo4 rather than importing that graph conclusion.

## Q2. The classical carry polynomial, with its exceptional sources

Put h=(p-1)/2 and alpha=(p+1)/2. Define

\[
W_p(z)=\begin{pmatrix}\alpha&hz\\h&\alpha z\end{pmatrix},
\qquad F_m(z)=\operatorname{tr}(W_p(z)^m).
\]

The trace counts closed carry sequences in digit addition by
r=(q-1)/2. The transition from incoming carry c to outgoing carry d
counts digits u with 0<=u<=p-1 satisfying

    u+h+c = v+p*d,       0<=v<=p-1,

and gives that transition weight z^d. Thus c=0 has alpha choices for
d=0 and h for d=1; c=1 has h and alpha choices. The carry is one bit,
but its accumulated weight is retained.

The quadratic denominator gives

\[
\boxed{F_0=2,\quad F_1=\alpha(1+z),\quad
F_m=\alpha(1+z)F_{m-1}-pzF_{m-2}.} \tag{Q1}
\]

This is the classical recurrence, not a new discovery of it. The local
application needs its source correction. Write s_p(i) for the base-p
digit sum of the least residue modulo q-1. For i in
I={1,...,q-2}\{r}, put

\[
c(i)={s_p(i)+s_p(r)-s_p(i+r\bmod(q-1))\over p-1}.
\]

It is the number of cyclic carries. The trace has four extra closed
records: source0 at weight0, sourceq-1 at weightm, and sourcer with
two carry choices of weights0 andm. Therefore

\[
\boxed{C_m(z)=\sum_{i\in I}z^{c(i)}=F_m(z)-2-2z^m.} \tag{Q2}
\]

These are three excluded digit words and four closed records, not four
vertices; the words0 andq-1 represent the same residue moduloq-1.
Forgetting the exceptional clocks gives the wrong Smith
rank. Also c(i+r)=m-c(i), and multiplication of i by p modulo q-1
rotates its padded digits and preserves c(i).

## Q3. Integral character blocks, with the sign appropriate to tournaments

Work over the unramified extension R of Z_p containing the (q-1)st roots
of unity, with Teichmuller character omega. Extending an integer matrix
to this ring preserves its p-primary Smith exponents: p is still a
uniformizer, and every unit remains a unit.

Transpose E, which preserves its Smith form. For Q=A^T-A, the matrix
entry in row y, column x is chi(y-x). In the free R-module on F_q put

    e_i=sum_(x!=0) omega(x)^(-i) [x].

The multiplicative Fourier basis is integral and invertible because
q-1 is a unit in R. With r=(q-1)/2, its invariant blocks are
span(e_i,e_(i+r)), 1<=i<=r-1, and span(1,[0],e_r).

For nonexceptional i, direct convolution gives

\[
Qe_i=J_i e_{i+r},\qquad
J_i=\sum_x\omega(x)^{-i}\chi(1-x).
\]

The classical valuation input gives v_p(J_i)=c(i). Since Q^2=J-qI,
where J is the all-ones operator, the paired coefficients satisfy

    J_i J_(i+r)=-q.

The minus sign matters. It is the tournament value chi(-1)=-1.
On the exceptional block,

    Q(1)=0,       Q([0])=e_r,       Q(e_r)=1-q[0].

These formulas are direct character sums; in particular J(chi,chi)=1.
They prevent an unnoticed import of the undirected sign convention.

Normalize g=gcd(a,b)>0, x=a/g, y=b/g, s=x+y, delta=x-y.
The primitive transpose pencil is

    E0^T=(s/2)(J-I)+(delta/2)Q.

On the nonexceptional two-dimensional block, twice this operator is

\[
\begin{pmatrix}-s&\delta J_{i+r}\\\delta J_i&-s\end{pmatrix}. \tag{Q3}
\]

In the ordered exceptional basis (1,[0],e_r), twice the operator is

\[
\begin{pmatrix}(q-1)s&s&\delta\\0&-s&-q\delta\\0&\delta&-s\end{pmatrix}. \tag{Q4}
\]

## Q4. Complete characteristic-primary factors

For s nonzero and p dividing s, set

    t=v_p(s),       d=min(2t,m).

Primitivity makes delta a p-adic unit. Because m is odd, 2t differs
from m, so N=(s^2+q delta^2)/4 has valuation exactly d: there is no
equal-valuation cancellation boundary in this tournament scope.

If c=c(i), the smaller Smith exponent of (Q3) is
e=min(t,c,m-c), from the valuation of its entries. Its determinant has
valuation d, so the other exponent is d-e. The exceptional block has
two unit factors, visible modulo p in (Q4), and last exponent t+d.

Let C_j be the coefficient of z^j in C_m. Pairing i with i+r gives
exactly C_j blocks of smaller carry j for 0<=j< m/2. Consequently the
full primitive p-exponent counting polynomial is

\[
\boxed{B_{m,t}(z)=2+z^{t+d}+
\sum_{j=0}^{(m-1)/2}C_j
\left(z^{\min(t,j)}+z^{d-\min(t,j)}\right).} \tag{Q5}
\]

If p does not divide s, every p-primary factor is a unit. If s=0,
the primitive pencil is plus or minus Q, and its nonzero Smith factors
are all powers of p, with exponent polynomial

\[
\boxed{S_m(z)=F_m(z)-2z^m,\qquad\text{plus one free factor}.} \tag{Q6}
\]

The factor2 at exponent0 comes from (Q4), while the free factor replaces
its singular last direction. This proves the complete p-primary rule,
including the prime case m=1.

**Weight saturation.** Once t>m/2, the nonexceptional factors equal
those of the skew character matrix; only the exceptional exponent t+m
continues moving. When 2t<m, carries j>=t merge into exponent t, while
the exceptional factor has exponent3t. These are precise truncation laws
for this pencil, not assertions about a Collatz branching rule.

## Q5. Assemble every integer Smith factor

For a=b=0 the matrix is zero. Otherwise first divide out g. If s=0,
use (Q6) and multiply all nonzero factors by g, leaving the free zero.
For s nonzero, put r=(q-1)/2, R0=|rs|, N=(s^2+q delta^2)/4 and
d0=gcd(R0,N). Form the ordered list

    1^[r], d0, N^[r-1], R0*N/d0.

Remove every power of p from each entry of that list. Then multiply
entrywise by the ordered p-powers from (Q5), or q unit factors if
p does not divide s. Finally multiply every entry by g.

This is the complete integer Smith form. Away from p the projector J/q
is integral, so the constant and augmentation proof in PALEY_PENCIL
applies unchanged and gives Z/(rs) plus r copies of Z/N. At p use Q4.
The p-free factors still form a divisibility chain, as do the sorted
p-powers; their entrywise product therefore is the Smith chain.

No factorization of the large norms is needed to construct the list.
The carry histogram requires O(m^2) integer arithmetic operations using
(Q1), rather than constructing a q-by-q matrix. Printing all q factors
still has output cost at least q; compressed multiplicity counts avoid
that output expansion. This is a response classification, not a faster
general Hamiltonian-path counting algorithm.

## Q6. Actual field controls and fixed-class lifting

For the singular skew pencil, the exact nonzero exponents are:

| Field | Multiplicities at exponents0,1,...,m | Free rank |
|---|---|---:|
|F27|8,6,6,6|1|
|F243|32,40,50,50,40,30|1|
|F343|64,108,108,62|1|

For weights (1,p-1), every field above has t=1. Its exponent
multiplicities at0,1,2,3 are respectively

```text
F27:   8,12,6,1
F243: 32,180,30,1
F343: 64,216,62,1.
```

The October19 lifting mechanism now has an exact new consumer. Let K_k
be the subspace of ker(E mod p) whose classes admit some compatible lift
to ker(E mod p^k). If the nonzero Smith exponent multiplicities are b_e
and the free rank is f, then

\[
\dim K_k=f+\sum_{e\ge k}b_e. \tag{Q7}
\]

This concerns existence of a lift of a fixed mod-p class, not extension
of every selected intermediate representative. For the F27 skew pencil,
the dimensions at k=1,2,3,4 are19,13,7,1. For F243 they are
211,171,121,71,31,1. This is the same integral lifting predicate as in
the older padding operator; it is not an identification of those matrices.

## Q7. Phase, provenance and what the small automaton forgets

The carry polynomial aggregates character sources. Frobenius sends
i to pi modulo q-1, rotating the digit and carry words. Each Jacobi
sum is constant on that source orbit, but the basis vectors and named
ports are transported. An aggregate multiplicity polynomial loses the
source orbit and its position. The receipt retains the entire index orbit
for every directly tested Jacobi sum; it does not quotient these phases
away before a later operation.

Thus a two-state automaton computes arbitrarily deep Smith multiplicities,
while an unweighted cycle graph or a single rank cannot encode all those
precision layers. This is a literal finite weighted automaton, not a
reduction from a Turing-machine halting question or a Collatz orbit.
The concrete source map is (character index,digit position,incoming carry)
to the next carry and valuation weight; the digit source is the sidecar
needed for transport.

The incoming [cyclic coefficient guard](../overnight_research_20261029/CYCLIC_COEFFICIENT_GLUE.md)
likewise retains a source lattice congruence before splitting character
sectors. Its prime3 guard and the exceptional block here have different
sources. The common methodological lesson is to retain the integral
embedding and named source, rather than infer a cross-domain isomorphism
from their determinant factors.

## Reproduction and audit boundary

`prime_power_pencil.py` verifies all556 earlier prime cases as exact
specializations; computes full integer Smith forms for90 F27 pencils;
and independently eliminates the actual F243 and F343 matrices for six
pencils each, checking every characteristic-primary layer. The latter
checks do not independently compute their other primary factors.

It constructs the fields with defining polynomials z^3-z-1 over F3,
z^5-z-1 over F3 and z^3-2 over F7, checks all nonzero elements are units,
and verifies the Frobenius action on every adjacency entry. In unramified
quotient rings modulo p^(m+1), it independently constructs Teichmuller
roots and evaluates604 Jacobi sums, checking their carry valuations,
paired product -q and source-Frobenius orbits. The recurrence is checked
through degree101 at p=3,7,11; direct digit enumeration is added whenever
p^m<=20000. All guards use exceptions active under optimized Python.

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/prime_power_pencil.py --write
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/prime_power_pencil.py
```

Remaining targets are integer three-run gluing and retention of the full
complemented relation under later operations. The all-prime-power pencil
itself is now closed in this packet, relative to the classical valuation
input explicitly cited above.
