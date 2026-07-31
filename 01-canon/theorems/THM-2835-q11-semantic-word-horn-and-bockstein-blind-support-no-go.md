---
id: THM-2835
title: "Q11 semantic word horn and Bockstein-blind support no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the fixed
  THM-2806/2829 atom,
  the 449 whole-cylinder QB(source)-to-QA(q11 target) labels form an exact
  coefficient horn whose beta-selected ancestry borrow lands in QAB at
  q7.  Both zero- and one-borrow QAB fillers cover the same 449 horn
  labels, so support is blind to the Bockstein carry.  The complete
  three-word relation is triangular and noninvertible:
  R11=R4 o R7 but R2 o R11 is empty.  Its canonical QB response has exact
  rational Hankel rank 13, while the full QB/QAB response bank has rank
  26; regular C13 modules attain both bounds with nonnegative
  input/output maps.  This is an abstract response realization.  The
  characteristic-13 carry quotient of THM-2788 is not a rational module;
  only free rational linearization of its 13 columns gives the matching
  formal character.  Canon supplies no typed basepoint identifying the
  physical address carry with this semantic response.  Literal +9U
  transport loses E3 on
  all 63 repaired cells, and no one of all 567 semantic cells repairs the
  full I/J0/J7 packet.  This is a coefficient/support theorem, not a
  current action, row exclusion, or LRC(14) proof.
audit: >
  audit_2809 independently reconstructed the word relations, both 446
  censuses, the 13/26 Hankel-rank invoice and regular realizations, the
  characteristic-13 versus free-rational-linearization boundary, the
  0/567 physical scan, the 35/14/10/4 E3 census, and the empty
  outer-gauge/endpoint-gauge fibre product; normal and optimized replay
  LF-normalize byte-identically to the stored output.  Final verdict:
  ACCEPT.
source: root/lrc-q11-word-horn-2026-07-28
depends_on:
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2829-q11-semantic-reselection-and-fine-ancestry-phase-obstruction
related:
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2792-cyclic-unit-intertwiner-and-positive-naturality-boundary
script: 04-computation/lrc14_q11_semantic_word_horn_thm2835.py
output: 05-knowledge/results/lrc14_q11_semantic_word_horn_thm2835.out
script_sha256: 207dd94f086338ae1e80b7d7196f99bf41e795893d13b6d48e4e7d516af03523
output_sha256: 1ebe0cbaf7d4ef13defed0bdb5b37df1364880acdbfc6139b243ab9df65f6bf6
hash_basis: LF-normalized bytes
---

# THM-2835 -- q11 semantic word horn and Bockstein-blind support no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2829 isolates a positive `q=11` carrier whose first failure is the
outer THM-2584 `QB` ancestry word.  Testing the complete three-word atlas
does more than locate another positive mask.  It exposes a small
semisimplicial object:

```text
QB(source,a) ---- q11 ----> QA(target,a)
       \                         |
        \                        | beta(11)=9, borrow +1
         \                       v
          -------- q7 ------> QAB(target,a+1).             (1)
```

The triangle is exact on `449` whole ancestry cylinders.  It is not,
however, a horn in the physical current category: its word relation is
one-way, its support cannot see whether the ancestry borrow occurred,
and literal transport destroys the macro `E3` factor.

There is also a positive result hidden in the failure.  The canonical
`QB` response is a cyclic regular `C13` module of dimension `13`; the
full two-source response is exactly two regular modules, of dimension
`26`.  Therefore neither lack of representation capacity nor abstract
positivity is the remaining obstruction.  The missing datum is a typed
basepoint identifying the physical carry column with this particular
semantic response orbit, together with physical transport of `E3`.

## 1. Fixed object and exact word bank

Use the fixed THM-2806/2829 data

```text
p=13,                    D=p^5=371293,
T=297836897838480,       U=T/13,
I=[142004992589460,142005019034340),
J_q=I+431933040+qU.                                      (2)
```

The three disjoint outer words have final two Boolean bits

```text
QA =(in,out),       QB =(out,in),       QAB=(in,in),      (3)
```

in slots `7,8`, whose physical speeds are

```text
W_7=13^3=2197,             W_8=2*13^5=742586.             (4)
```

For `X,Y in {QA,QB,QAB}`, ancestry `a in Z/DZ`, and `q in F_13`, let

```text
S_X(a)    = whole-cylinder source membership in X,
T_Y(q,a)  = whole-cylinder J_q target membership in Y,
C_XY(q)   = sum_a S_X(a)T_Y(q,a).                         (5)
```

Every test in the companion compares midpoint membership with full
cylinder containment.  They agree everywhere, so no half-open seam is
being promoted to a positive sheet.

The source counts are

```text
sum S_QA=0,             sum S_QB=66099,
sum S_QAB=10786.                                           (6)
```

The nonzero directed counts are

```text
C_QB,QA  =449(delta_6+delta_8+delta_9+delta_10
                    +delta_11+delta_12),
C_QB,QB  =66099 delta_0+65612 delta_7,
C_QB,QAB =449 delta_7,

C_QAB,QA =(0,10786,10785,10785,10783,10780,10779,
             0,10329,10329,10328,10328,10328),
C_QAB,QAB=10786 delta_0+10779 delta_7.                    (7)
```

Every count with source `QA`, and every `QAB -> QB` count, is zero.

## 2. The 449-sheet horn and the complete carry census

Define

```text
L={a:S_QB(a)=T_QA(11,a)=1}.                              (8)
```

Then

```text
|L|=449,       min L=59306,       max L=311961,
a=156 mod169 for every a in L.                           (9)
```

All `449` hits are whole cylinders.

For each `a in L`, inspect all three target words after ancestry shifts
`delta=-2,-1,0,1,2`.  In word order `(QA,QB,QAB)`, the exact census is

```text
             delta=-2     -1       0        +1       +2

q=0          (0,446,0) (0,446,0) (0,449,0) (0,0,449) (0,0,449)
q=7          (0,446,0) (0,447,0) (0,0,449) (0,0,449) (0,0,449)
q=11         (0,0,0)   (0,0,0)   (449,0,0) (449,0,0) (449,0,0).
                                                                  (10)
```

In particular, if

```text
B_delta(a)=T_QAB(7,a+delta),
A_delta(a)=T_QB (7,a+delta),                              (11)
```

then on `L`

```text
(sum B_-2,...,sum B_2)=(0,0,449,449,449),
(sum A_-2,...,sum A_2)=(446,447,0,0,0).                 (12)
```

Thus both `B_0` and `B_1` fill every horn sheet.  The binary word support
does not distinguish zero borrow from the compulsory `+1` borrow.

The two labels which fail the reverse `delta=-1` return are

```text
a=107978: slot 1, speed 14, becomes dangerous;
a=154622: slot 4, speed 53, becomes dangerous.           (13)
```

Their outer bits still have the requested `QB=(out,in)` state.  These are
typed low-factor failures, not failures of the outer `QA/QB` switch.

Two numerically equal `446` counts must not be conflated:

1. restricted to `L`, the two-step reverse count is
   `sum_L A_-2=446`;
2. globally on the full `QB` source, the `q7` `QAB` support at shifts
   zero and one has sizes `449` and `895`.  The first is a strict subset
   of the second, so there are `446` new one-borrow sheets, but every one
   has `T_QA(11,a)=0` and is disjoint from `L`.

The equality of those two numbers supplies no bijection or current.

THM-2820's selected label `beta(11)=9` has

```text
(a,11) -> (a+floor((11+9)/13),11+9 mod13)
        = (a+1,7).                                        (14)
```

By `(10)`, every `a in L` lands in `QAB` and in neither `QA` nor `QB`.
The target word transition is

```text
QA=(in,out) -> QAB=(in,in),                               (15)
```

so it flips only outer slot `8`.  Equations `(8),(14),(15)` are the exact
coefficient-support horn in `(1)`.

## 3. The relation is triangular, not a group action

Let `R_q` be the binary relation

```text
X R_q Y  iff  C_XY(q)>0.                                 (16)
```

The complete bank is

```text
R_0       : QB->QB, QAB->QAB;
R_1,...,R_5: QAB->QA;
R_6       : QB->QA, QAB->QA;
R_7       : QB->QB, QB->QAB, QAB->QAB;
R_8,...,R_12: QB->QA, QAB->QA.                           (17)
```

Give the words the filtration

```text
rank(QB)=2,          rank(QAB)=1,          rank(QA)=0.    (18)
```

Every nonzero arrow weakly lowers rank.  Direct relational composition
gives the constructive nonwrapping factorization

```text
R_11=R_4 o R_7.                                          (19)
```

But

```text
R_2 o R_11=empty != R_0,            2=-11 mod13.         (20)
```

Hence `(17)` is neither an action of `C13` nor a functor from the
physical address/future groupoid.  It is a one-way cospan/profunctor.
The apparent 2-simplex fills a support horn but fails the inverse horn.

There is no workaround on the three words themselves.  Any set action is
a homomorphism

```text
C13 -> S3,                                                (21)
```

which is trivial.  Likewise, if an invertible nonnegative `3 by 3`
matrix `M` satisfies `M^13=1`, then `M^-1=M^12` is nonnegative.  A
nonnegative matrix with nonnegative inverse is monomial.  Its permutation
part is trivial by `(21)`, and its positive diagonal entries are
finite-order positive reals, hence all one.  Thus every finite `C13`
action on the three-dimensional simplicial word cone is trivial.

Over `Q`, the same conclusion follows from

```text
x^13-1=(x-1)Phi_13(x),              deg Phi_13=12.        (22)
```

A nontrivial rational `C13` representation has dimension at least `12`,
so a three-dimensional rational action is trivial.

## 4. Exact 13- and 26-dimensional response realizations

The preceding no-go concerns an action on the three visible word states.
It does **not** forbid a larger hidden response module.

Tabulate the source-by-target binary matrices

```text
R_q(X,Y)=1_{C_XY(q)>0},                                  (23)
```

and the weighted matrices

```text
K_q(X,Y)=C_XY(q).                                        (24)
```

For a thirteenth root `z`, put

```text
F_R(z)=sum_q R_q z^q.                                    (25)
```

At `z=1`,

```text
F_R(1)=
 [0 0 0
  6 2 1
 11 0 2],                     rank=2.                    (26)
```

For `z!=1`, the two active rows are

```text
QB : [f(z),       1+z^7, z^7],
QAB: [-(1+z^7),  0,     1+z^7],                         (27)

f(z)=z^6+z^8+z^9+z^10+z^11+z^12.
```

Since `-1` is not a thirteenth root, `1+z^7!=0`; the last two
coordinates in `(27)` give rank `2`.  The weighted transform has the
same rank because its two triangular pivots are

```text
66099+65612z^7,             10786+10779z^7,              (28)
```

and neither can vanish on the unit circle: the two coefficient magnitudes
in each binomial are unequal.

Now restrict to the canonical `QB` source row.  Its binary Fourier vector
is nonzero at every character (`z^7` already occurs in its `QAB`
coordinate), and its weighted version has the nonzero term `449z^7`.
Thus its character rank is one at all thirteen characters.

Equivalently, let

```text
H_(q,X),(h,Y)=K_(q+h)(X,Y),                              (29)
```

with indices modulo `13`, and define the analogous binary Hankel matrix.
Exact Gaussian elimination modulo the prime `1000003` gives the lower
certificates

```text
rank_Q H_QB =13              (binary and weighted),
rank_Q H_full=26             (binary and weighted).       (30)
```

Since reduction modulo a prime cannot increase rank over `Q`, these are
rational lower bounds.  The explicit regular realizations below give
matching rational upper bounds, hence the equalities in `(30)`.

These ranks are minimal hidden-state dimensions.  Indeed, any realization

```text
K_q=O rho(g)^q B                                           (31)
```

through a rational `C13`-module `V` factors the block Hankel matrix
through `V`, so `dim V>=rank H`.  Fourier decomposition gives the same
invoice:

```text
QB source:       1 trivial + 1*12 charged =13;
full bank:       2 trivial + 2*12 charged =26.             (32)
```

Both lower bounds are attained constructively.  Regard the displayed
source-by-target row as the map sending a source basis vector to that
target row.  For the `QB` lane take

```text
V_QB=Q[C13],
B(e_QB)=delta_0,
rho(g)delta_q=delta_(q+1),
O(delta_q)=K_q(e_QB).                                    (33)
```

For the full bank let `S=span_Q{QB,QAB}` and take

```text
V_full=Q[C13] tensor S,
B(s)=delta_0 tensor s,
rho(g)(delta_q tensor s)=delta_(q+1) tensor s,
O(delta_q tensor s)=K_q(s).                              (34)
```

Use the binary matrices instead of `K_q` for the support realization.
All displayed input/output matrices are nonnegative and `rho` is a
permutation action.  Therefore `13` and `26` are exact minima, not only
lower bounds.

The `13` shifted `QB` response columns are independent, so

```text
delta_q -> shift^q(QB response)                           (35)
```

is the unique equivariant isomorphism once its value at `delta_0` is
chosen.  Its forward matrix is entrywise nonnegative; its inverse need
not be.  Abstract algebraic naturality is therefore complete after one
basepoint is supplied.

## 5. What THM-2788 supplies, and what it does not

For `p=13`,

```text
Q[C13]=Q*1 direct_sum Aug_Q(C13),          dim Aug=12.    (36)
```

THM-2788 supplies a physical `13`-column permutation carrier on the
address side.  Its quotient `B/<Z>` is an elementary abelian
characteristic-`13` object; it is **not** literally the rational module
`Aug_Q(C13)` (tensoring that `13`-torsion quotient with `Q` gives zero).
If one instead freely `Q`-linearizes the underlying `13` carry columns,
one obtains the regular rational module in `(36)`.  Equation `(32)` then
shows the following exact formal capacity match:

```text
canonical QB response = one freely linearized carry carrier;
full QB/QAB response = two freely linearized carry carriers.          (37)
```

Equivalently, after this free rational linearization, the charged
dimensions are one and two copies of `Aug_Q(C13)`.  The rational
augmentation module alone has no nonzero invariant pointed convex cone.
If `v` lies in such a cone, then all translates do, while
`sum_q g^qv=0`; hence `-v=sum_{q=1}^{12}g^qv` also lies in the cone.  The
positive regular cone appears only after adjoining the trivial mass line.

But THM-2788 lives on address carry columns.  It supplies neither the
choice

```text
physical delta_0  <->  canonical QB response,             (38)
```

nor the physical input/output maps in `(33)`.  By THM-2611, equivariant
identifications of two fixed `C13` torsors form another `C13` torsor; a
basepoint or a zero-holonomy section is load-bearing.  Thus `(37)` is
only a formal rational-linearized dimension/Fourier-support match.  It
does not say that THM-2788 already supplies the rational response module,
its holonomy section, or a physical intertwiner.  The choice `(38)`
remains the typed sidecar.

This is the sharp conceptual boundary:

```text
not missing abstractly: a response realization or positive regular cone;
missing physically:     rationalized carry-to-response basepoint,
                        lawful intertwiner, and macro-factor transport. (39)
```

## 6. Literal physical transport still fails

The selected interval satisfies

```text
J_11+9U=J_7.                                              (40)
```

On the `63` THM-2829 repairing cells, however, the six-factor containment
signature at `J_7`, in factor order

```text
(E3,clock,q1,q2,c2,c3),                                  (41)
```

has exact census

```text
(0,1,1,1,1,1): 35,
(0,1,1,0,1,1): 14,
(0,1,0,1,1,1): 10,
(0,1,0,0,1,1): 4.                                       (42)
```

Every cell loses `E3`; `28` also lose `q1` or `q2`.  None is complete.
The `35` E3-only cells are exactly

```text
s in {0,1,2,3,8,9,12},       t in {5,6,9,10,11},
clock=1.                                                   (43)
```

The larger all-`567` semantic-cell scan finds no cell whose six factors
simultaneously contain `I,J_0,J_7`.  Clock `1` is the unique clock
carrier containing `I,J_0,J_11,J_7`, so a clock reselection cannot repair
the failure.  The missing map must transport or replace `E3`; merely
choosing another semantic cell does not suffice.

There is a marked-gauge copy, but it has exactly the wrong source typing.
Inside the `2197`-point marked orbit, `+4W` (equivalently the inverse
`+9` move) takes the `q11` marked representative to a `q7`
representative with the same mass and endpoint scalars.  The target gauge
simultaneously changes source harmonic `h0 -> h4` and moves the marked
piece with the carrier.  It therefore loses the original fixed source
`I/h0` and the `449`-sheet horn.  The canonical fixed-source `q7`
coefficient remains zero across all `169` E3 right-endpoint
representatives.  Using the gauge copy would require a new
changed-sheet-to-original-source clutch.

This is not repaired by reindexing ancestry.  The outer representative
gauge has `s_4=4*13^4=114244`; translating the horn labels by `-s_4`
preserves all `449` outer `QB/QA/QAB` incidences.  The endpoint gauge,
however, moves the source interval to `I-4U`.  On that moved source all
`449` reindexed labels have `QA=QB=QAB=0`; the first lost predicate is
the required slot-`8` `b`-danger.  Hence the outer-gauge horn and the
endpoint-gauge mass have empty fibre product.

Finally, the outer slots in `(4)` are not address digits.  On the same
physical carrier the low address changes

```text
8 -> 7,                                                   (44)
```

while the semantic `q=0` bank has `66099` literal `QB -> QB` sheets.
The graph label `q`, outer slot numbers `7,8`, and address residues `8,7`
are three distinct typed coordinates.  Identifying them would manufacture
the missing action.

## 7. Exact consequence and exact non-consequence

What is proved by this theorem:

1. an exact `449`-sheet coefficient/support horn;
2. complete zero/one/two-borrow word censuses and two separately typed
   `446` phenomena;
3. an exact triangular relation with a factorization but no inverse;
4. impossibility of a nontrivial direct `C13` action on the three words;
5. exact minimal positive hidden response dimensions `13` and `26`;
6. exact formal capacity matching after free `Q`-linearization of the
   `13` carry columns, with the characteristic-`13` typing boundary;
7. a precise missing sidecar: the physical basepoint `(38)` plus `E3`
   transport.

What is **not** proved:

- no coefficient in `(7)` is promoted to a canonical THM-2334/2365
  current;
- no lawfully co-shifted complete packet realizes `(17)`;
- no ancestry/address basepoint in `(38)` is selected;
- no fixed old `X,m`, root gauge, owner phase, or endpoint word is
  preserved;
- no row is excluded and the LRC(14) ledger remains unchanged.

## 8. Exact companion

Run

```text
python 04-computation/lrc14_q11_semantic_word_horn_thm2835.py
python -O 04-computation/lrc14_q11_semantic_word_horn_thm2835.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_q11_semantic_word_horn_thm2835.out.
```

The independent hostile audit replayed both modes, reconstructed the
relation, Fourier/Hankel, physical-cell, and marked-gauge arguments by a
separate path, and accepted the theorem after the characteristic-`13`
typing repair in Section 5.
