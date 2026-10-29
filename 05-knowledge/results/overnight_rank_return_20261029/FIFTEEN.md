# Four storage bits and five incidence queries for fifteen Lagrangians

**PROVED elementary geometry and transport; FINITE-EXACT complete query
census; CITED arithmetic realization. Same-agent audit.** This extends the
incoming [symplectic versus quadratic boundary](../sha_marked_join_f1_20261029/F1_QUADRATIC_BOUNDARY.md)
and [Richelot and Frobenius carrier](../session97_overnight_arithmetic_carriers_20261029/supersingular/RICHELOT_FROBENIUS_AND_FANO.md).
The Plucker geometry is classical; the explicit observation-cost comparison
is the target here. No external novelty claim is made.

## G1 The fifth bit is determined by a quadratic equation

Write V=F2^4 with coordinates (x1,x2,y1,y2), alternating form
omega(u,v)=x.u*y.v+y.u*x.v, and Q0(x,y)=x1*y1+x2*y2. A Lagrangian is a
two-dimensional omega-isotropic subspace. For a basis u,v define
p_ij=u_i*v_j+u_j*v_i, with indices 0,1,2,3. Because GL2(F2) has determinant
one, these Plucker coordinates are independent of the chosen basis.

Put

\[
(a,b,c,d,t)=(p_{01},p_{03},p_{12},p_{23},p_{02}).
\]

Isotropy says p13=p02. The Plucker equation becomes

\[
\boxed{t=ad+bc.} \tag{4}
\]

Thus the four bits (a,b,c,d), other than 0000, determine the Lagrangian.
Conversely each nonzero four-bit vector gives a nonzero bivector satisfying
the Plucker equation and the isotropy equation, hence a unique Lagrangian.
There are exactly fifteen. Four storage bits are optimal by cardinality.
The fifth coordinate is free through a **nonlinear** relation.

This is an explicit map, not an arbitrary numbering of fifteen objects.

## G2 Symplectic transport closes on the four stored bits

Every g in Sp4(F2) acts linearly on exterior-square coordinates and preserves
the five-dimensional isotropy hyperplane. Its invariant bivector is

\[
\Omega=e_0\wedge e_2+e_1\wedge e_3,
\]

which has (a,b,c,d,t)=(0,0,0,0,1). Quotienting by this invariant line leaves
a four-dimensional linear action on (a,b,c,d). Equation (4) then recovers
the unique isotropic decomposable lift. Therefore the four-bit carrier is
closed under every specified symplectic change of frame, including a named
Frobenius action on geometric two-torsion in odd characteristic.

[fifteen.py](fifteen.py) exhausts all 720 symplectic matrices, derives each
induced four-by-four binary matrix and verifies every action on all fifteen
Lagrangians. The induced action is faithful. Its matrices must be transported
through this construction; substituting the original four-by-four matrix
without the exterior-square quotient is not justified.

If F is the induced Frobenius matrix, the rational Lagrangian kernels are
exactly the nonzero vectors of ker(F-I). Their number is therefore

\[
\boxed{2^{\dim\ker(F-I)}-1\in\{0,1,3,7,15\}.} \tag{5}
\]

The closure here is for symplectic transport. It makes no claim about the
path-native complement operation on tournament attachments.

## G3 Natural quadratic guards need five queries

For v in V let Q_v(x)=Q0(x)+omega(v,x). Exactly ten of these refinements are
split, namely those with Q0(v)=0. The allowed binary query is

\[
b_v(L)=\mathbf1[Q_v\text{ vanishes on }L].
\]

Each query is true on six of the fifteen Lagrangians, and each Lagrangian
answers true to four of the ten queries. The full response distinguishes
all fifteen. With v encoded by bits (x1,x2,y1,y2), the explicit formulas are:

| v | b_v in the five Plucker coordinates |
|---:|---|
| 0 | a+b+c+d+t |
| 1 | a+b+t |
| 2 | a+c+t |
| 3 | a+t |
| 4 | c+d+t |
| 6 | c+t |
| 8 | b+d+t |
| 9 | b+t |
| 12 | d+t |
| 15 | t |

All sums in the table are in F2. Substitution on the fifteen subspaces is
a complete 150-case verification; the script also reconstructs each unique
linear coefficient vector independently.

**Five queries are necessary, even adaptively.** The first query partitions
the bank into six and nine objects. After a false answer, three further
binary answers can distinguish at most eight of the nine possibilities.
Thus depth four is impossible. Five suffice: for example the queries
v=0,1,2,3,4 separate every Lagrangian. The complete census finds 162
five-query subsets that do so, and a dynamic program verifies adaptive
worst-case depth five.

Frobenius can lower that cost if its actual fixed-kernel set is known:

| Number of fixed kernels | Number of symplectic matrices | Minimum fixed queries | Minimum adaptive worst case |
|---:|---:|---:|---:|
| 0 | 304 | 0 | 0 |
| 1 | 300 | 0 | 0 |
| 3 | 100 | 2 | 2 |
| 7 | 15 | 3 | 3 |
| 15 | 1 | 5 | 5 |

Zero denotes an empty candidate set, not an existing kernel reconstructed
without data. These are complete finite-action counts, not probabilities
for a distribution of arithmetic curves. The actual phase matrix tells
which queries suffice, not only how many fixed kernels exist.

## G4 Weighted responses have a different minimum

Let B be the ten-by-fifteen zero-one incidence matrix. Over F2 its row rank
is five: the displayed formulas lie in five-coordinate linear space, and
five independent rows are present. Thus arbitrary **parity-weighted**
superpositions require five linear coordinates, although one actual
Lagrangian needs only four bits. The nonlinear equation (4) does not survive
arbitrary addition of objects.

For rational weights the rank is ten. Two distinct guards share exactly
two Lagrangians, so

\[
BB^{\mathsf T}=4I_{10}+2J_{10},\qquad
\det(BB^{\mathsf T})=24\,4^9\ne0.
\]

All ten independent weighted responses are therefore necessary for an
arbitrary rational mixture. They are sufficient and closed under symplectic
transport because that action permutes the ten quadratic forms.

We have three sharply different minima for one geometry:

| Task | Minimum |
|---|---:|
| Store one Lagrangian as a deterministic object | 4 bits |
| Recover it with allowed quadratic guards | 5 binary queries |
| Linear carrier for all parity-weighted responses | dimension 5 over F2 |
| Linear carrier for all rational-weighted responses | dimension 10 over Q |

This is a concrete reason to keep representation cost, query cost and
weighted composition cost separate. A bit can be reconstructed for free on
the realization locus and become necessary after taking superpositions.

## The arithmetic connection and its boundary

For a principally polarized abelian surface over a field of characteristic
different from two, geometric A[2] is symplectic F2^4. Maximal isotropic
subgroups give polarized (2,2)-isogeny kernels. The fifteen-object geometry
therefore has a direct arithmetic realization. Frobenius stability, handled
by (5), is the rationality guard. Incoming Session97 supplies concrete
supersingular examples at p=11 with fifteen versus three rational kernels;
the supersingular j-value alone does not determine this response.

Ogg's fifteen absolute supersingular primes instead classify when every
supersingular j-invariant lies in its prime field. The cited
[Ogg mechanism](https://celebratio.org/Ogg_AP/article/1157/) identifies the
genus of the Fricke quotient with the number of nonrational conjugate pairs
of supersingular j-values. That is a Frobenius-orbit criterion on a different
set. A bijection from the prime list to the fifteen kernels is not supplied
by either construction. Characteristic two also requires a separate
group-scheme treatment and is excluded from the A[2] statement above.

The useful unification is now precise: retain a pairing, a realization
equation, and the transported Frobenius matrix. The equation can recover
a coordinate on individual objects, while weighted composition may restore
the need for that coordinate or a larger linear carrier.

## Reproduction

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/fifteen.py --write
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/fifteen.py
```

The receipt retains all fifteen subspaces, ten guards, 720 transported
matrices, exact ranks and the complete conditional query census. The
original 81/27 tournament pair remains a hostile to replacing this named
action with separate local summaries; no theorem here identifies those
tournaments with the fifteen-kernel geometry.
