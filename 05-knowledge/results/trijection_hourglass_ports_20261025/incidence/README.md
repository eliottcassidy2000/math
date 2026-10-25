# What the 6–4–6 actually measures, and how it grows

2026-10-25. **PROVED** abstract box theorem; **FINITE-EXACT** physical
three-direction controls below, upgraded by the
[complete physical kernel theorem](common_kernel.md) in its stated range.
This is a continuation of the linear incidence
sidecar, not a new norm-graph failure bound or a Campedelli identification.

## Inheritance

The relevant recovered 6–4–6 is in
[the sixteen-point hostile, section 6](../../section_defect_two_error_geometry_20261019/README.md).
There u=(1,0), v=(0,1), w=(3,9) over F29 and
B={0,-u,-w,-u-w}. Two disjoint eight-point cubes have union
T=B+{-1,0,1,2}v. The middle cube is also supported in T.
The three line-sum spaces W_u,W_v,W_w have dimensions 6,4,6; their
external direct sum has dimension 16. Their weighted pointwise sum has
rank 10 and kernel dimension 6. These are dimensions, not a partition
of sixteen code quotients. The extra two kernel modes over the four
endpoint-cube modes are the middle cube in two coefficient channels.

The closest mechanism is that explicit connecting cube. The hostile is
the false implication "disjoint supports imply independent incidence".
The corrected near miss is counting supports without their shared lines.
The least-used sidecar is the full three-direction line partition.

## A structural growth theorem

Take an abstract box with coordinate sets of sizes 2,L,2, L>=2, over any
field F. A coordinate line varies one coordinate with the other two fixed.
Let W_A consist of functions summing to zero on every B-line and C-line;
define W_B,W_C cyclically. For nonzero scalars s_A,s_B,s_C define

    M: W_A direct-sum W_B direct-sum W_C -> F^(2L2),
    M(f_A,f_B,f_C)=s_A f_A+s_B f_B+s_C f_C.

**Theorem.**

    (dim W_A,dim W_B,dim W_C)=(2L-2,L,2L-2),
    rank M=3L-2,             dim ker M=2L-2.

Moreover every kernel vector is a linear combination of the two relation
channels of the L-1 adjacent cubes. Thus the intermediate cubes form a
complete constructive sidecar for this abstract box, at every length.

**Proof.** Write A=F^2, B=F^L, C=F^2 and let A0,B0,C0 be the kernels of
coordinate summation. Then

    W_A=A tensor B0 tensor C0,
    W_B=A0 tensor B tensor C0,
    W_C=A0 tensor B0 tensor C.

Choose a vector of sum one in each factor; this splits each factor as its
sum-zero subspace plus a one-dimensional complement, in every
characteristic (no division by L is needed). The three spaces share
J=A0 tensor B0 tensor C0, of dimension L-1. Outside J their respective
pieces occupy three distinct direct summands of dimensions L-1,1,L-1.
The nonzero coefficients make M injective on those three pieces. On the
three copies of J its kernel is the two independent relations among
s_A,s_B,s_C, tensored with J. Hence its nullity is 2(L-1); its image has
dimension (L-1)+(L-1)+1+(L-1)=3L-2. Adjacent differences form a basis of
B0, while each of A0,C0 is one-dimensional. Their tensor products are
the L-1 adjacent cubes, proving the generator assertion. QED.

For L=4 this gives exactly 6–4–6, rank10, kernel6 on16 points. The number
16 is explained by a 2-by-4-by-2 support together with its line operators.

## Embedding needs its incidence data

The theorem transports to an actual support in a plane only when the
embedding identifies precisely the three coordinate-line partitions with
the three physical line partitions. Injectivity on points alone is not
sufficient. Physical lines can combine coordinate lines, changing their
constraints without merging any support points.

`check_box.py` independently constructs physical supports

    T_L={-a u-b w+j v: a,b in{0,1}, 0<=j<L}

over F29, with the inherited u,v,w. Translation by -v recovers the
sixteen-point support at L=4. It forms every physical line-sum equation
and the pointwise coupling (-3,-9,1), and computes ranks by two separate
elimination procedures. It also builds the adjacent cubes, verifies every
equation on them, and checks their independence.

The full universe is L=2,...,12, with no quotient or genericity filter.
All supports have exactly4L distinct points. The resulting dimensions are:

| L | dim(W_u,W_v,W_w) | coupling rank | kernel dimension | extra beyond adjacent cubes |
|---|---|---|---|---|
| 2 | 2,2,2 | 4 | 2 | 0 |
| 3 | 4,3,4 | 7 | 4 | 0 |
| 4 | 6,4,6 | 10 | 6 | 0 |
| 5 | 9,5,8 | 14 | 8 | 0 |
| 6 | 12,6,10 | 18 | 10 | 0 |
| 7 | 15,7,12 | 22 | 12 | 0 |
| 8 | 18,8,14 | 26 | 14 | 0 |
| 9 | 21,9,16 | 30 | 16 | 0 |
| 10 | 24,11,18 | 35 | 18 | 0 |
| 11 | 27,13,21 | 39 | 22 | 2 |
| 12 | 30,15,24 | 43 | 26 | 4 |

Already at L=4 there are extra physical collinearities, even though the
dimensions agree with the abstract box. At L=5 the individual-space
dimensions detect a difference; the coupled kernel still has exactly the
adjacent-cube basis through L=10. At L=11 two genuinely new coupled modes
appear. Thus keeping only point count, or even the early kernel formula,
would incorrectly extrapolate this family. The follow-up
[common-kernel theorem](common_kernel.md) now classifies these modes:
one explicit twelve-point bridge and its translates supplement the
adjacent cubes. It proves nullity2(L-1)+2max(0,L-10) throughout
2<=L<=p-9, including all L<=20 at p29. Modular wrapping beyond that
range remains a separate problem.

Positive controls at p=101,w=(17,37), L=2,...,10 agree with the abstract
dimensions and generators. These controls are finite and do not replace
the proof or certify a generic plane embedding outside their range.

## Connection to the current lanes

The source-to-target map here is the labelled box embedding; it preserves
support points, while the extra sidecar is the physical line partition.
For ordered tournament charts the corresponding obligation is to retain
directed joins and actual positional reversal, not merely the edge-XOR
word. For hourglass codes it is to retain the actual branch triples, not
merely seven labels. This is a shared method of preserving constraints;
no map identifying the three mathematical objects is supplied.

Reproduce from this directory with `python check_box.py` and
`python -O check_box.py`. No external libraries, floating point, or disabled
assertions are used. `normal.json` and `optimized.json` are the receipts.
