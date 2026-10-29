# Exact precision growth for the fifteen kernel response carrier

**PROVED from an explicit unimodular certificate / FINITE-EXACT alternate
image enumeration.** The map is the ten-by-fifteen quadratic-guard incidence
matrix B of [FIFTEEN](FIFTEEN.md), with arbitrary integer weights on the
fifteen Lagrangians. This is a closed weighted response carrier under named
symplectic transport, not a general tournament attachment classification.

## P1 The integer image and its five hidden directions

There are unimodular integer matrices U,V such that

\[
UBV=[D\ 0],\qquad
D=\operatorname{diag}(1,1,1,1,1,2,2,2,2,4). \tag{7}
\]

[precision.json](precision.json) contains U,V and every elementary operation
that produces them. [precision.py](precision.py) checks (7) by dense integer
multiplication, then separately replays the row and column operations.
Every operation is unimodular, and the positive diagonal entries divide
their successors. Thus (7) is a complete Smith certificate, not a numerical
rank estimate or an unverified library output.

Consequently the image lattice B Z^15 has index 64 in Z^10, and

\[
\operatorname{coker}B\cong(\mathbb Z/2)^4\oplus\mathbb Z/4.
\]

The binary rank five and rational rank ten are two endpoints of a precise
integral filtration. Four directions first appear one precision level later,
and the last direction appears two levels later.

## P2 An all precision recurrence

Let C_k be the number of possible complete response vectors modulo 2^k.
The domain consists of arbitrary weight vectors modulo 2^k. Unimodular
changes preserve the image cardinality, so (7) gives

\[
\boxed{\log_2 C_k=5k+4\max(k-1,0)+\max(k-2,0).} \tag{8}
\]

In particular

\[
C_0=1,\quad C_1=32,\quad C_2=16384,\quad
C_{k+1}=1024C_k\quad(k\ge2).
\]

The first successive bit costs are 5,9,10,10,..., and the cumulative costs
are 0,5,14,24,34,.... Independent subgroup enumeration gives exactly 32
images modulo two and 16384 modulo four. Higher precisions follow from
the certified Smith theorem, not from enumerating their exponentially many
weight vectors.

For k>=2 the integer lattice constraints save six bits relative to ten
unconstrained k-bit responses. They do not save five whole coordinates at
all precisions. A bit omitted over F2 can reappear as a carry at the next
level, and one direction in this example needs two levels before appearing.

The image is closed under every named symplectic frame change: that action
permutes the fifteen object weights and the ten quadratic guards, preserving
B and its image lattice. The Smith coordinates themselves depend on the
chosen unimodular bases; transport must be conjugated into those bases.

## P3 A small general theorem for realization and precision

For any integer observation matrix M of rank r with nonzero Smith factors
d1,...,dr and prime l, write ei=v_l(di). Then the number of reachable
responses modulo l^k is

\[
|\operatorname{im}(M\bmod l^k)|
=l^{\sum_i\max(k-e_i,0)}. \tag{9}
\]

This is the elementary Smith decomposition of the source and response
modules. It specifies exactly when a parity-level dependency remains free
and when it returns at higher precision. If a family of matrices has an
unbounded list of ei or unbounded rank, a fixed low-precision summary cannot
be assumed to close under later observations. A genuine compression theorem
must name the coefficient ring, precision budget and target predicates.

For single realized Lagrangians, the stronger nonlinear relation
t=ad+bc saves the fifth storage bit. For mixtures modulo two, the five
linear coordinates are necessary. For integer mixtures, (7) governs the
additional precision directions. These are three different operations on
the same explicit source, so their different information costs have a
mechanism rather than a numerical analogy.

The question for actual strong-atom attachments is now concrete: identify
their integer joint-response maps and compute the invariant factors before
deciding that a binary dependency is safe under repeated complement,
re-marking or weighted summation. The 81/27 control still requires retention
of the transporter data, which scalar response ranks do not certify.

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/precision.py --write
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/precision.py
```
