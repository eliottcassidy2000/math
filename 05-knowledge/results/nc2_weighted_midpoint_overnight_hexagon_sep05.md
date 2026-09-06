# Weighted midpoint deletion: all-root one-factor signs and a sharp abstraction failure

Status: **PROVED / INDEPENDENTLY AUDITED** for the weighted
network and square-pencil deductions (root and observer full written
audits). The explicit abstract joint-sign
implication below is **REFUTED**. The actual binomial joint signed
transport remains **OPEN**; its finite positive controls are not a proof.
No theorem ID is allocated to this sidecar.

## 1. Inheritance and the retained coordinate

The closest proved mechanism is the
[full-support binomial-path factorization](overnight3_20260906_moments_root_geometry.md).
The new consumer is the
[virtual Hadamard doubling and midpoint-defect proof](nc2_hadamard_transport_overnight_hexagon_sep05.md),
which proves the actual-minus-virtual coefficient inequality but leaves
its negative-root evaluation open. Its virtual sign uses
[THM-4440, signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md).
The canonical actual-versus-virtual hostile is support `(-3,1,9)`:
at the first root `-1`, its actual row is `-224`, not the virtual `-20`.
The corrected near miss is truncating the negative-index auxiliary
support; the discarded prefix `21+20t+5t^2` is not real-rooted.

The least-used sidecar is the **selected midpoint vertex set**, not
only its path count. Deleting or weighting that set preserves a planar
Toeplitz network. The live board is: full path support; weighted
deletion; degeneration at a repeated root; factor versus joint
Hadamard operations; endpoint normalization; carry degree.

The source-to-target map in this note retains every path and its
midpoint-hit indicator. Its preserved predicate is total
nonnegativity of the full Toeplitz path matrix. Passing to abstract
real-rooted pencils loses the unit-step endpoint relation; Section 5
proves that this lost information is essential to the desired joint sign.

## 2. Full-support weighted deletion is again a PF sequence

Let `p,q` be positive integers, and `U,V` nonnegative integers with
`U+V>0`. There is no coprimality assumption. Define the entire finite
two-sided sequence

```text
a_j=binom(U+V+(p-q)j, V+p*j)
       if U-q*j>=0 and V+p*j>=0, and zero otherwise;
b_j=binom(2U+2V+(p-q)j, 2V+p*j)
       if 2U-q*j>=0 and 2V+p*j>=0, and zero otherwise;
h_j=sum_l a_l a_(j-l),       d_j=b_j-h_j.                    (1)
```

Write the corresponding Laurent polynomials `a(t),b(t),d(t)`, so
`d=b-a^2`. For every real `w>=0`, the full coefficient sequence of

```text
c_w(t)=d(t)+w*a(t)^2=b(t)+(w-1)*a(t)^2                     (2)
```

is nonnegative and is a Pólya-frequency sequence. Consequently every
nonzero polynomial obtained by clearing its negative powers has only
real nonpositive roots. Zero rows are allowed; added zero roots and
multiple roots must not be suppressed.

Here is the exact network. Put `T=(-q,p)` and `C=pU+qV>0`. The first
row counts northeast unit paths from `(0,0)` to `(U,V)+jT`. The doubled
Toeplitz entry `b_(j-i)` counts paths from

```text
S_i=iT    to    E_j=(2U,2V)+jT.
```

The selected midpoint vertices are all

```text
M_l=(U,V)+lT,       l in Z.                               (3)
```

Sources, selected vertices and sinks lie on the respective parallel
levels `pX+qY=0,C,2C`. Each northeast edge increases this functional
by the strictly positive amount p or q. Hence a path can visit at
most one selected vertex. Paths through `M_l` are counted by
`a_(l-i)*a_(j-l)`; summing l gives exactly `h_(j-i)`.

Give every selected vertex the same weight w and every other vertex
weight one. A path has weight w if it hits the selected set and
weight one otherwise. Thus the weighted entry is exactly (2), not a
polynomial with higher powers of w. At `w=0` this is literal deletion
of the selected vertices. Paths that jump over the midpoint level,
and paths that hit an unselected lattice vertex on it, are retained.

Translation by T preserves the graph and the weighted vertex set,
so the path matrix is Toeplitz. For any finite minor, the usual
first-intersection tail swap cancels intersecting families, preserving
their vertex-weight product. Source and sink order on the parallel
levels forces an increasing pairing for a disjoint family. The
remaining weights are nonnegative. Every finite Toeplitz minor is
therefore nonnegative. Only finitely many paths are involved in any
fixed entry or minor, despite the infinite indexing.

The cited finite-PF characterization is exactly the one already
checked for the incoming factors:
[Brändén, arXiv:math/0403364v2, Section 2, Theorem 2.1 and its finite consequence](https://arxiv.org/pdf/math/0403364#page=2).
No closure theorem for sums of real-rooted polynomials is invoked.
The new proof is the explicit weighted/deleted network, followed by
that inherited characterization.

The original alpha factor is `(p,q,U,V)=(A,A,m-z,z)`; the beta factor
is `(p,q,U,V)=(B-A,B,y,x)`. Thus this theorem applies to both complete
factors in the preceding transport note. It proves real-rootedness
of each **individual** midpoint defect. It does not identify their
joint defect with one of these single networks.

## 3. Square-pencil degeneration retains multiplicity and sign

The following elementary lemma makes the root-position consequence
precise. Let S and D be real polynomials, with S nonzero. Suppose
`D+w*S^2` is real-rooted for an unbounded sequence of positive w.
Let rho be a real root of S of multiplicity `r>=1`. If D is not zero,
then

```text
ord_rho D >= 2r-2.                                       (4)
```

At exact order `2r-2`, the leading coefficient of D in powers of
`t-rho` is strictly negative. In particular

```text
D(rho)<=0;    if r>=2, then D(rho)=0.                     (5)
```

For `D=0` use infinite order; all these conclusions hold. The lemma
requires only an unbounded positive sequence, not the entire pencil.

Proof: write `S=(t-rho)^r(s_0+...)`, with `s_0!=0`. If
`nu=ord_rho D<2r`, put `D=(t-rho)^nu(d_0+...)`, `d_0!=0`, and
`N=2r-nu`. Remove the common real factor `(t-rho)^nu`, substitute

```text
t-rho=w^(-1/N)*z,
```

and take the coefficient limit along the given w sequence. The
resulting nonzero limit is

```text
s_0^2*z^N+d_0.                                          (6)
```

Real-rooted polynomials of bounded degree are closed under nonzero
coefficient limits, allowing degree drop; this follows directly from
continuity of polynomial roots, with any escaping roots sent to
infinity. For `d_0!=0`, (6) can have only real roots only if `N<=2`.
If `N=2`, its two roots are real only when `d_0<0`. This proves
(4)--(5). No claim about the sign of the leading coefficient is made
at the allowed odd order `2r-1`.

The boundaries are exact. With `S=(t+1)^r`, each of

```text
D=0,
D=-(t+1)^(2r-2),
D=(t+1)^(2r-1)
```

gives a real-rooted pencil for all positive w. Changing the minus
sign in the second example creates a nonreal quadratic pair. For
`r>=2`, replacing its exponent by `2r-3` creates a nonreal cubic pair.

For a simple root there is an independent one-line check: the
elementary Laguerre inequality `P'^2-P*P''>=0` for real-rooted P,
evaluated at rho with `P=D+wS^2`, has w-leading term
`-2w D(rho) S'(rho)^2`. It forces `D(rho)<=0`. This derivative proof
alone does not handle multiple roots; the scaling argument above does.

## 4. Every individual first root has a nonpositive true doubled value

For the Laurent path rows (1), every nonzero real root rho of a
satisfies

```text
b(rho)<=0.                                               (7)
```

If its multiplicity in a is `r>=2`, then b vanishes at rho. More
precisely the Laurent defect `b-a^2` has order at least `2r-2` there,
with a strictly negative local leading coefficient at equality.

To keep all signs literal, choose an **even** integer E large enough
that `t^(E/2)*a` and `t^E*b` are polynomials. For example
`E=2*(floor(V/p)+1)` always works. The weighted network theorem gives
the real-rooted pencil

```text
t^E*(b-a^2)+w*(t^(E/2)*a)^2.
```

Apply Section 3. Since `rho^E>0`, the sign is the sign of the raw
Laurent value b, with no lower-carry ambiguity. At an a-root the
defect value equals b. Multiplying by this nonvanishing positive
factor also preserves the local order and leading sign.

Strictness is not automatic. For `p=q=1`, every path hits the full
midpoint level at a selected vertex, so `b=a^2`, the defect is zero,
and every first root has doubled value zero. No individual-factor
coprimality theorem is claimed here.

This is an all-root statement for an individual auxiliary factor.
It is not a constant-term detection theorem for the original
trinomial, whose first row is a Hadamard product of two such factors.

## 5. Compatibility, joint real-rootedness and strict individual signs still fail

The smallest displayed obstruction is entirely explicit. Let

```text
F(t)=G(t)=1+t,
D(t)=2t(2t+3)^2,
A(t)=F(t)^2+D(t)=1+20t+25t^2+8t^3.                       (8)
```

For **every** `w>=0`, `D+wF^2` has nonnegative coefficients and is
real-rooted. Indeed its cubic discriminant is

```text
4w(2w^2+9w+432),                                         (9)
```

strictly positive for `w>0`; all coefficients are then positive, so
its three real roots are negative. At `w=0`, the roots are zero and
the repeated root `-3/2`. Thus both copies of this factor pair satisfy
the entire nonnegative weighted-pencil conclusion, not merely a
finite test. Both even satisfy the strict individual sign
`A(-1)=D(-1)=-2`.

Nevertheless, for ordinary coefficientwise multiplication `star`,

```text
F star G=1+t,                      first root rho=-1,
F^2 star G^2=1+4t+t^2,             value -2,
A star A=1+400t+625t^2+64t^3,      value 162.              (10)
```

The joint actual-minus-virtual analogue is

```text
396t+624t^2+64t^3=4t(99+156t+16t^2),                     (11)
```

which is itself real-rooted with nonnegative coefficients: its
quadratic discriminant is `18000`, both quadratic roots are negative,
and the remaining root is zero. Yet its value at rho is **164>0**.

Therefore coefficientwise inclusion, full nonnegative weighted-pencil
compatibility, strict individual first-root signs, and even joint
real-rootedness do not imply the desired joint sign. The first failed
implication is transporting those factor predicates through joint
Hadamard endpoint doubling. The missing sidecar is the exact unit-path
endpoint/coefficient relation. This is an **abstract preserver hostile**,
not a claim that (8) is an actual binomial endpoint-doubled row and
not a counterexample to the open trinomial statement.

A narrow degree boundary explains where this failure starts. In the
same normalization `F=G=1+t`, suppose A and B have degree at most two,
dominate `(1+t)^2` coefficientwise, and satisfy `A(-1),B(-1)<=0`.
Then `A_1>=A_0+A_2`, `B_1>=B_0+B_2`, and

```text
(A star B)(-1)
 <=-A_0 B_2-A_2 B_0<=-2.                                (12)
```

Thus the same abstract predicates do force the joint sign through
quadratic degree; (8) fails at cubic degree. This is a normalized
explanatory survivor, not a new substitute for the stronger proved
[THM-4432, trinomial two-channel two-rung noncancellation with carries](../../01-canon/theorems/THM-4432-trinomial-two-channel-two-rung-noncancellation-with-carries.md).

## 6. Exact controls and the remaining target

The standalone checker is
[nc2_weighted_midpoint_overnight_hexagon_sep05.py](../../04-computation/nc2_weighted_midpoint_overnight_hexagon_sep05.py).
It imports no repository mathematical producer. Its complete network
universe is `p,q=1..3`, `U,V=0..3`, excluding `(U,V)=(0,0)`, and
`w in {0,1/2,1,2,7}`: 135 parameter cases and 675 weighted rows.
Zero rows, repeated roots, negative auxiliary indices and nonprimitive
directions are explicitly retained. A literal two-state path DP
independently counts paths that hit or miss the selected vertex set.
All Toeplitz minors of sizes two and three from index sets `{0,1,2,3}`
are checked using separate rational determinant formulas.

Every distinct nonzero first root is either an exact shared root of
the actual doubled row or has a rational-isolation negative-value
certificate. Repeated-root defect divisibility is checked separately.
The degeneration controls include both sharp allowed orders and their
two nonreal-root hostiles. The cubic discriminant (9) is checked as a
symbolic identity, not by sampling w. The quadratic boundary has an
independent finite integer bank as a control for its elementary proof.

The open target is still the sign of the **specific binomial joint
missed-midpoint row at a root of its first Hadamard row**. Deletion
proves a stronger individual predicate than bare real-rootedness, but
(8)--(11) prohibit discarding its endpoint relation before attempting
the joint step. Root placement and the actual coefficient map, not
only the count of real roots, remain essential.

Both root and observer independently audited the complete written
argument, including the full-support network, zero/multiple-root
boundaries, coefficient-limit proof, even Laurent clearing and exact
abstract hostile. The full producer passes **38,974 explicit gates**:
639 literal path endpoints, 35,100 Toeplitz minors, 168 strictly
negative individual doubled root values, 15 shared-root zeros, 87
degeneration controls and 729 quadratic-boundary controls. Its 675
weighted rows include 21 genuinely zero rows. No gate is a Python
assertion.

Reproduce from the repository root:

```bash
python3 -B 04-computation/nc2_weighted_midpoint_overnight_hexagon_sep05.py
python3 -B -O 04-computation/nc2_weighted_midpoint_overnight_hexagon_sep05.py
```

Normal, optimized and [stored output](nc2_weighted_midpoint_overnight_hexagon_sep05.out)
agree byte for byte. Raw LF SHA256 values are:

```text
source 86d7bbc9e586f5563dfa408ad4d3b70396c96f2bf4b940e58ea5bbf4a3ac26f4
output 855cb32478debb5c1375dc27ad733f42ad6068ad8424535d6e5210b1fa089e13
trace  67d2a13a236fe2f00399382fa661e28d44c92dcf840c60c9dd060fd1bb093870
```
