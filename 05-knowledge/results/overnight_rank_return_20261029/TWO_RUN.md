# Actual two-run modules expose integral glue and the Frobenius action

2026-10-29 heartbeat. **PROVED identities / FINITE-EXACT module certificates.**
Same-agent alternate-algorithm audit. These are ordered joint endpoint
responses of actual tournaments. They do not retain every complemented
target relation and therefore do not close the full metagraph carrier.

## R1. Recover the old exclusion carrier and specify the new observer

The closest proved mechanisms are the square-zero port resolvent in
[Session35](../session35_dimensions_dualities_20261023/recursions/CONTEXT-AND-PORTS.md)
and the ordered cover convolution in
[Session69 P1--P3](../session69_prime_ports_power_aliquot_20261026/session69-tournament/REPORT.md).
[Session37](../session37_statistics_parity_cuts_20261023/exclusion/RESULT.md)
already separates unsigned vertex exclusion from signed exterior
cancellation. These identities and their exponential exclusion memory are
inherited. The added object is their integer joint-endpoint response module.

For a tournament T on n named vertices and ordered distinct start and end
pairs (a,b),(c,d), define

\[
M_T[(a,b),(c,d)]
=\#\{(P_1,P_2):\text{ordered spanning disjoint directed paths},
\ P_1:a\to c,\ P_2:b\to d\}. \tag{R1}
\]

Singleton paths are allowed, so a=c or b=d is permitted. Equal vertices
in different paths are forbidden. The matrix has n(n-1) rows and columns.
The weights are raw labelled counts, not automorphism-orbit counts.

In the inherited ring Z[x_v]/(x_v^2), put P=X(I-AX)^(-1). Then

\[
M_T[(a,b),(c,d)]=[x_V]P_{ac}P_{bd}. \tag{R2}
\]

Multiplication imposes disjointness, and full support imposes spanning.
More generally, the ordered k-run tensor is [x_V] product_i P_(s_i,t_i).
This is an exact joint-incidence formula; no independence of the runs is
assumed. Complement targets still require Session69's full rooted
relation U_C, all ordered component endpoints and multiplicity.

The anchor is integral joint response. The niche is the recovered signed
versus unsigned sector distinction. The wildcard is the actual order-21
transport action and its quadratic norm. The live board is exclusion /
joint ports / integral gluing / transporter / precision / realization.

## R2. Exact modules on three small strong atoms

The complete counts and Smith factors are:

| Tournament | Ordered two-run covers | Matrix order | Smith factors |
|---|---:|---:|---|
| C3 | 6 | 6 | 1^[6] |
| T5 with steps1,2 | 90 | 20 | 1^[18],82^[2] |
| Paley7 | 1638 | 42 | 1^[24],2^[6],4^[3],12^[6],25351872^[2],177463104 |

Square brackets denote multiplicity. In particular, the Paley7 matrix has
full rational rank42, parity rank24, and complete two-adic exponent list

\[
0^{[24]},\quad1^{[6]},\quad2^{[9]},\quad6^{[3]}. \tag{R3}
\]

Thus, for this one-step weighted endpoint observer,

\[
\log_2|\operatorname{im}(M\bmod2^k)|
=24k+6(k-1)_++9(k-2)_++3(k-6)_+. \tag{R4}
\]

After k=6, each new digit multiplies the response count by 2^42.
The eventual count is 2^(42k-42). Full rational rank means all 42
coordinates of an arbitrary rational mixture of named start pairs can be
reconstructed from the ending responses. It is a linear-rank statement,
not a lower bound on arbitrary encodings of individual tournaments.

## R3. Run symmetry splits over Q but has parity glue over Z

Swapping both runs preserves M. Order the pair coordinates as forward and
reversed copies of the unordered pairs. Then

\[
M=\begin{pmatrix}B&C\\C&B\end{pmatrix},\qquad M_+=B+C,\ M_-=B-C.
\]

Over Q the sum and difference coordinates diagonalize this matrix.
Over Z their exact image is instead the fiber product

\[
\boxed{\{(u,v)\in\mathbb Z^m\oplus\mathbb Z^m:u\equiv v\pmod2\}},
\qquad m=\binom n2. \tag{R5}
\]

Indeed u=x+y and v=x-y, with inverse x=(u+v)/2, y=(u-v)/2.
The congruence is exactly the realization guard; M_+=M_- modulo2
ensures that the split action preserves it. The guard disappears over
Z_ell for odd ell, but cannot be dropped at ell=2.

For Paley7 the exact sector Smith factors are

```text
M_+: 1^[12],6^[3],12^[3],96^[2],1248
M_-: 1^[12],2^[3],4^[3],81256,1056328,7394296.
```

The unglued direct sum would have two-adic exponents
0^[24],1^[6],2^[6],3^[3],5^[3], instead of (R3). It undercounts the
actual image by a factor eight modulo8,16,32. The determinant agrees,
so a determinant-only check would miss this error.

The sectors have a literal determinant/permanent meaning:

    (M_+)[{a,b},{c,d}] = [x_V](P_ac P_bd + P_ad P_bc),
    (M_-)[{a,b},{c,d}] = [x_V](P_ac P_bd - P_ad P_bc).

These are the permanent and determinant of the two-by-two port matrix,
after square-zero exclusion. The nonzero integral extension is why
knowing the two abstract sector cokernels is insufficient to recover
their joint integer module. This refines the inherited statistics
distinction with an actual small strong-atom witness.

## R4. The literal 21 action, 42 coefficients, and two norm fields

For Paley7 use G={x->a x+b:a in {1,2,4}, b in F7}, of order21.
It acts regularly on directed arcs. Every ordered pair is a unique arc
or reversed arc, so the 42 ports form the regular G x C2 set; C2 swaps
the two runs. Simultaneous transporter covariance makes M a convolution
matrix on this set. Its 42 first-row coefficients describe every entry.
This compresses the operator description, not its full-rank input space.

For each sector use the arc (0,1) as source and g(0,1)=(b,b+a) as target.
The coefficient kernel k(g) obeys M_(g,h)=k(g^-1 h). Summing its seven
translation coefficients at a=1,2,4 gives respectively

```text
plus:  (15,12,12)
minus: (-11,2,2).
```

The constant character has eigenvalue h=sum c_i: 39 and -7.
The other two characters of G/C7=C3 contribute

\[
L=c_0^2+c_1^2+c_2^2-c_0c_1-c_0c_2-c_1c_2,
\]

which is an Eisenstein norm, here 9 and169. This is a precise point of
contact with the incoming [ternary-clock norm calculation](../overnight_research_20261029/EISENSTEIN_CLOCK_JOINS.md).
The source modules and allowed operations differ; equality of the norm
form does not identify them.

Fourier characters of the normal C7 fall into two nonzero orbits
{1,2,4} and {3,5,6}. The order-three multiplier rotates each orbit.
Their three-dimensional blocks are conjugate over Q(sqrt(-7)). If
zeta is a primitive seventh root and eta=zeta+zeta^2+zeta^4, their first
block determinants are exactly

\[
D_+=24-24\eta=12(3-\sqrt{-7}),
\qquad D_-=248+608\eta=8\sqrt{-7}(38+\sqrt{-7}). \tag{R6}
\]

Since eta+bar(eta)=-1 and eta bar(eta)=2, their norms are

\[
N_+=2304,\qquad N_-=650048=2^6\cdot7\cdot1451,
\qquad1451=38^2+7. \tag{R7}
\]

The regular representation has the two three-dimensional constituents
with multiplicity three, so

\[
\det M_\pm=h_\pm L_\pm N_\pm^3. \tag{R8}
\]

For a direct elementary derivation, Fourier-transform the C7 coordinate
in the 21-dimensional regular representation. The zero character gives
the three C3 characters; each nonzero orbit gives three copies of its
three-dimensional block. Dimensions 3+9+9 exhaust the space.

The checker computes the three-by-three blocks modulo Phi7, verifies
invariance under zeta->zeta^2, and independently recovers each norm from
the determinant of the seven-point affine permutation representation.
The same {1,2,4} phase is the actual Frobenius multiplier in the earlier
F8 realization of Paley7. F8 and Q(sqrt(-7)) are different fields; their
specified index action is the bridge.

Here 21 is literally the tournament's transporter group, and also the
number of unordered pairs on seven vertices on which it acts regularly.
No relation to a rank-21 elliptic heuristic is inferred. The new prime
1451 has an explicit quadratic-norm source rather than a numerical analogy.

## R5. Two runs detect the 81/27 control at the first ternary digit

Both twelve-vertex controls still have total Hamilton count462915 and
ordered two-run count9011898. Their complete transporter groups are
again enumerated, with orders81 and27. Nevertheless,

\[
\boxed{\operatorname{rank}_{F_3}M_{81}=24,
\qquad\operatorname{rank}_{F_3}M_{27}=22.} \tag{R9}
\]

This is invariant under relabelling. The single-run matrices vanish
modulo3 and first separate modulo9; joint run incidence separates them
already modulo3. This is a stronger observer, not a counterexample to
the earlier sharpness statement in its single-run channel.

All 132 ternary invariant exponents were determined by modular
elimination through digit12. The exact multiplicities are

| Exponent | 81 control | 27 control |
|---:|---:|---:|
|0|24|22|
|1|21|31|
|2|59|49|
|3|19|24|
|4|8|4|
|5|0|1|
|6|1|0|
|9|0|1|

Every omitted exponent has multiplicity zero. Consequently their
precision exponents are

    f81(k)=24k+21(k-1)_++59(k-2)_++19(k-3)_++8(k-4)_++(k-6)_+,
    f27(k)=22k+31(k-1)_++49(k-2)_++24(k-3)_++4(k-4)_++(k-5)_++(k-9)_+.

They eventually become 132k-234 and132k-231. At the first digit the
order81 control has more responses, but at the second it has fewer
(exponents69 versus75). Group order alone does not order these response
sizes, and no mutual simulation is inferred from cardinalities.

## R6. Scope, guards and reproduction

The six-vertex revisit hostile is replayed in full. Child bitcode7 and
covers ((0,1),(3,2)), ((0,2,1),(3,)) have the same complemented child
type with two outer endpoints. In C3[S4,1,1], their literal complemented
targets have canonical codes81 and83. Joint ordered component ports
separate the child responses. Thus none of the present rank formulas
licenses discarding the complemented relation from Session69's carrier.

The checker uses two independent enumeration routes: subset-path dynamic
programming with disjoint support contraction, and complete permutation
plus cut enumeration through order7. Selected entries through order12
also receive an independent walk inclusion-exclusion check: sum two
walks with n total vertex occurrences and impose full support, which
forces simple disjoint spanning paths. Every transporter action is
checked on the full joint matrix.

Small-atom Smith forms are cross-checked by independent modular unit
elimination at every determinant prime and an exact Bareiss determinant.
The larger hostile calculation certifies only its complete 3-primary
layers, not every primary factor. Explicit exceptions remain active
under optimized Python.

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/two_run.py --write
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/two_run.py
```

The resulting research direction is concrete: retain the integral gluing
of symmetry sectors alongside their individual responses. Denominator
primes of the symmetry projectors are the first places to look for lost
alignment. Here they are2 for run exchange andq for the constant/augmentation
split in [PALEY_PENCIL](PALEY_PENCIL.md). The proposed next test was a
three-run attachment, where S3 introduces both2 and3 as gluing primes.
[THREE_RUN T1--T5](THREE_RUN.md) now verifies both primes on actual
Paley7 responses, gives the exact index54 source guard, and exhibits
T5's four hidden directions destroyed by every legal one-step complement.
The full metagraph minimality question remains open.
