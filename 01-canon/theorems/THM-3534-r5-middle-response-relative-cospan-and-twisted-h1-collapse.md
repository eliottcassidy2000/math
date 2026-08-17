---
id: THM-3534
title: "R5 middle-response relative cospan and twisted H1 collapse"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER INDEPENDENT AUDIT; NOT IN THE
  PROOF GRAPH.  On the frozen r5 diagonal bank, the full middle-orientation
  difference has rank five and splits into a rank-two source plane plus a
  rank-three endpoint plane.  The dual middle block of the ten-dimensional
  quotient contracts to a different rank-two response plane.  The two planes
  become canonically isomorphic only relative to the unique endpoint-supported
  line at r0=6.  The descended chamber involution splits this relative plane
  1+1, so its natural twisted C4 H1 has dimension one; every response row is
  already exact on the formal digit C13, the literal common chamber line is
  exact, and the surviving chamber line requires the endpoint quotient.  This is a
  static finite-exact representation statement, not a physical current, D5
  flux map, row exclusion, or LRC(14) theorem.
source: codex r5 rank-two cospan session, 2026-08-16
depends_on: []
related:
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
  - THM-3354-inequivalent-h1-carriers-and-typed-obstruction-cospan
  - THM-3431-d5-secondary-h1-descent-defects-and-valuation-persistence
  - THM-3450-marked-d5-carrier-isomorphism-and-full-germ-margin-obstruction
  - THM-3496-marked-graph-kummer-degree-square-and-finite-coefficient-frobenius-flux-extinction
script: 04-computation/lrc_r5_rank_two_relative_response_cospan_thm3534_20260816.py
output: 05-knowledge/results/lrc_r5_rank_two_relative_response_cospan_thm3534_20260816.out
script_sha256: 4eab23ff77de9d190f13d5a41945415ba55734e0fcdc6f8eae3cbb47d4fb101d
output_sha256: f6972a28ba83ce03c7bfa45a2e6b6eaa8af1504f778a2412ee536639a53f4b6c
semantic_sha256: 57a44964888fdb0a9ca1c890abbe4950c6fa7130f7ab2800c0faa8a5d6a0212d
audit_script: 04-computation/lrc_r5_rank_two_relative_response_cospan_flint_audit_20260816.py
audit_output: 05-knowledge/results/lrc_r5_rank_two_relative_response_cospan_flint_audit_20260816.out
audit_script_sha256: 13316f9add30bb0d5de23d0d0339cf3bdcf666419045a84a3a6e1b96eb54e996
audit_output_sha256: a449fb52016453596ffc581a4df0ccde2028d8777fec6696c122f41ff24bd8d6
audit_semantic_sha256: f1dfa230027e7fd8f48c7eb7bac2034d0188ee61a4ce6f24d2579860cf532020
hash_basis: LF-normalized bytes
---

# THM-3534 -- R5 middle-response relative cospan and twisted H1 collapse

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER INDEPENDENT AUDIT.  DO NOT USE
THIS FILE AS A PROVED DEPENDENCY UNTIL ITS STATUS IS EXPLICITLY PROMOTED.**

The refuted `r1`-blind scalar transport has a minimal two-dimensional repair,
but only after one endpoint direction is quotiented.  This repairs a static
response representation.  It does not construct a clock edge, a local system
on the physical address carrier, a word-current, or a map to Jacobian flux.

## 1. Frozen universe and variance

Work over the split prime field

```text
k=F_p,       p=755373809845391722745761.                 (1)
```

The hash-pinned third-digit certificate records six parent kernels

```text
P_e=(k_e(r0,r1))_(r0,r1 in F13),       e=0,...,5.         (2)
```

The six pointed arcs, in order, are

```text
0->1, 1->0, 1->3, 3->1, 3->2, 2->3.                     (3)
```

The common-ten-quotient calculation uses the transposed parent blocks

```text
B_e=P_e^T,       Q_e=k^13/Row(B_e),
Q=direct_sum_e Q_e.                                     (4)
```

Its exact block ranks and quotient dimensions are

```text
rank(B_e)=(11,11,12,12,11,11),
dim(Q_e) =( 2, 2, 1, 1, 2, 2),       dim(Q)=10.          (5)
```

Every one of the thirteen fixed-`r2` children surjects onto the profile in
`(5)`.  Pure arc reversal fails by two dimensions, both on the middle pair;
coupled chamber reflection is exact.

The transpose in `(4)` is load-bearing.  The canonical object that contracts
against response rows is the quotient **dual**

```text
Q_e^vee=Ann(Row(B_e))=ker(P_e^T),                         (6)
```

not an unidentified copy of `Q_e` itself.

The other frozen input is the exception-location tent

```text
h=(12,12,9,3,0,0),
h_odd=(0,0,3,-3,0,0).                                   (7)
```

As already proved by the hostile audit, `(7)` is location data, not response
amplitude.  The construction below replaces one scalar coefficient by two
marked response functions; it does not change that semantic type.

## 2. The rank-five middle response

Put

```text
D=P_(1->3)-P_(3->1).                                    (8)
```

Its support consists of exactly the following 26 address pairs:

```text
r0=0:      r1=11,12;
r0=3:      r1=1,2,3,4,5,7,8,9,10,11;
r0=6:      r1=5,7;
r0=9:      r1=1,2,3,4,5,7,8,9,10,11;
r0=12:     r1=0,1.                                      (9)
```

The twenty entries in rows `3,9` are the old source-middle hostile.  The six
entries in rows `0,6,12` are the endpoint residuals.  Define

```text
S=span_k{D_(3,*),D_(9,*)},
E=span_k{D_(0,*),D_(6,*),D_(12,*)},
R=Row(D).                                                (10)
```

Exact row reduction gives

```text
dim S=2,       dim E=3,       rank D=5,
R=S direct_sum E.                                        (11)
```

In particular, the full orientation defect is rank five, not rank two.  The
rank-two source restriction is minimal for the two reflected exceptional
chambers, and its ten nonroot entries have eight distinct values in each row.
This reproves the scalar hostile while identifying the omitted endpoint
sector exactly.

There is also an all-row cohomology hostile.  Every parent block in `(2)` is
row-sum one:

```text
P_e 1=1                    for every e.                  (11a)
```

Subtracting the two middle equations gives

```text
D 1=0.                                                     (11b)
```

Hence every row of `D`, and therefore every vector in `S`, `E`, and `R`, has
zero seam on the formal oriented digit cycle `C13`.  Prefix integration gives
an explicit vertex potential for all thirteen rows.  Thus

```text
image(R -> H^1_graph(C13;k))=0.                           (11c)
```

The quotient tells the same story contravariantly: the constant digit vector
lies in `Row(B_e)` for every arc, so the constant Fourier mode has quotient
rank zero in all six blocks.  Therefore neither the response plane nor `Q10`
contains an ordinary digit-cycle `H^1` class.  This does not make `r1` a
physical clock; it closes that formal interpretation if one tries it.

## 3. Contraction of the dual middle quotient

Let

```text
L_+=Q_(1->3)^vee,       L_-=Q_(3->1)^vee,
L=L_+ direct_sum L_-.                                    (12)
```

Both lines in `(12)` are one-dimensional.  Choose the deterministic RREF
generator `ell_+` of `L_+`, and link the second marking by chamber reflection:

```text
ell_-(r0)=ell_+(12-r0).                                  (13)
```

The two lines are distinct.  Contract them against `(8)`:

```text
delta_D:L -> R,
delta_D(ell_+,0)=ell_+^T D=:i_+,
delta_D(0,ell_-)=ell_-^T D=:i_-.                        (14)
```

The map `(14)` is injective.  Its image

```text
I=span_k{i_+,i_-}                                        (15)
```

has dimension two, but it is not the source plane:

```text
dim(S+I)=3,       dim(S intersect I)=1.                  (16)
```

Thus the two-dimensional defect of the quotient and the two-dimensional
source response are not the same two-plane.  Equal ranks were hiding a
row/column variance error.

## 4. The unique endpoint-relative repair

The mismatch in `(16)` is exactly one endpoint direction.  Put

```text
E_6=k D_(6,*).                                           (17)
```

Then

```text
(S+I) intersect E = E_6,                                (18)
dim E_6=1.
```

Consequently

```text
V_rel=(S+I)/E_6                                         (19)
```

has dimension two, and both quotient maps are isomorphisms:

```text
S --sim--> V_rel <--sim-- I <--delta_D-- L.              (20)
```

If `Loc_(3,9)=k e_3 direct_sum k e_9`, the marked static response map

```text
tau_2:Loc_(3,9) -> S,
tau_2(e_3)=D_(3,*),       tau_2(e_9)=D_(9,*)             (21)
```

is also an isomorphism.  Equations `(20)--(21)` give the precise replacement
for the impossible scalar map:

```text
Loc_(3,9) --tau_2--> S --sim--> V_rel <--sim-- I
                                      <--delta_D-- L.     (22)
```

The word **relative** is load-bearing.  Before quotienting `E_6`, the two
planes are distinct.

With the chamber-linked normalization `(13)`, row reduction gives

```text
a=194124906132980974193833,
b= 56688580149000751581631,
c=404323792951742926212914                              (23)
```

in `k`, and the exact lift equations are

```text
D_(3,*)=a i_+ + b i_- + c D_(6,*),
D_(9,*)=b i_+ + a i_- + c D_(6,*).                      (24)
```

Both `a+b` and `a-b` are nonzero, and

```text
a^2-b^2=58994564309246164003190 !=0 in k.               (25)
```

Thus `(20)` is not merely a dimension count.  Equation `(24)` gives its
explicit marked transition and its unique endpoint gauge.

## 5. Chamber representation and the twisted-H1 collapse

On response functions define the arc-odd chamber involution

```text
(C f)(r1)=-f(12-r1).                                     (26)
```

It acts by

```text
C D_(3,*)=D_(9,*),       C i_+=i_-,
C D_(6,*)=D_(6,*).                                      (27)
```

Thus `(22)` is `C2`-equivariant.  Adding and subtracting `(24)` gives the
failure anatomy:

```text
D_(3,*)-D_(9,*)=(a-b)(i_+-i_-),                         (28)
D_(3,*)+D_(9,*)=(a+b)(i_++i_-)+2cD_(6,*).               (29)
```

The literal intersection in `(16)` is exactly the chamber-anti-invariant
line in `(28)`.  The invariant lines agree only in the relative quotient,
and their lift discrepancy is the nonzero endpoint term `2cD_(6,*)`.

There is still no lawful `2--0` owner/clock edge.  Nevertheless, the natural
formal test is now decisive.  If a hypothetical `C4` closure uses `V_rel` as
a local coefficient system with monodromy `C`, its cellular cochain complex
on the circle is

```text
V_rel --(C-I)--> V_rel.                                 (30)
```

Since `p` is odd and `C` has one invariant and one anti-invariant line,

```text
rank(C-I)=1,
dim H^0(C4;V_rel,C)=dim H^1(C4;V_rel,C)=1.               (31)
```

More sharply, the literal common line `(28)` is the image of `C-I` and is
zero in twisted `H^1`.  The surviving class is the invariant line `(29)`,
whose two natural lifts differ by the endpoint sidecar `2cD_(6,*)`.

Using the trivial local system would give `dim H^1=2`, but then chamber
reflection is an external decoration rather than the only symmetry that
actually descends from the frozen quotient.  Hence the frozen data supply no
two-dimensional chamber-faithful `H^1` class.

There is a second sharp boundary.  Stabilizing the parent rowspace under pure
arc reversal fills both 13-dimensional middle blocks.  Therefore the repaired
eight-dimensional quotient `Q_A` has zero middle dual:

```text
(Q_A,middle)^vee=0.                                      (32)
```

One may retain the rank-two interface or force arc reversal to descend, but
the latter operation kills the former.

## 6. Why this still does not type a D5 flux map

The cospan `(22)` is a static good-reduction representation.  It is not a
physical word-current for four independent reasons.

1. `Loc_(3,9)` records exception locations, and `(21)` is the marked lookup
   of their two response rows.  No THM-2334/2512 current-to-location chain
   map has been constructed.
2. The same-copy `2--0` edge closing the owner path to `C4` is still absent.
3. `Q^vee`, not `Q`, enters `(14)`, and every child section of `Q` remains
   noncanonical.  No chronology or composition law is gained by dualizing.
4. The endpoint line `(17)` is killed, not reconstructed.  A physical target
   that depends on it cannot factor through `(19)`.
5. On the formal digit cycle all response rows are coboundaries by
   `(11a)--(11c)`; there is no hidden nonzero address-seam class to export.

At the exact finite coefficient level there is also an immediate D5 no-go.
The additive group of `V_rel` has exponent `p`, while the marked Kummer line
of THM-3496 has exponent thirteen.  Since `p!=13`,

```text
Hom_Add(V_rel,F13)=Hom_Add(F13,V_rel)=0.                 (33)
```

There is no unital ring map from `k` to `F13` or `F2`.  Likewise, every
additive map from this finite `p`-group to a characteristic-zero Hamiltonian
response module is zero.  A characteristic-zero reconstruction of `(22)`
would be new data not contained in the finite certificate; it would still
face THM-3354's site/predicate obstruction, THM-3450's margin and gauge
obstructions, and THM-3496's additive-flux/Frobenius hostile.

Equations `(31)--(33)` do not exclude a genuinely derived or nonadditive
correspondence carrying the endpoint row, full word packet, and integral
Hamiltonian filtration.  They do exclude the proposed ordinary rank-two
local-system shortcut from the frozen inputs.

## 7. Tournament and XOR audit

The intrinsic binary relation on the two marked middle orientations is the
chamber involution.  It exchanges both ways.  It is not an antisymmetric
dominance relation and therefore is not a tournament edge.  On all four owner
states the frozen relation remains the bidirected path with census

```text
(both-way,one-way,missing)=(3,0,3),                     (34)
```

not a tournament on four vertices.

The `+/-` decomposition in `(28)--(29)` is a lawful use of the involution.
Calling it XOR would lose the coefficient field and, in characteristic two,
would erase the nonzero mismatch `2cD_(6,*)` while merging the invariant and
anti-invariant lines.  Because there is no coefficient map from `k` to `F2`,
that Boolean shadow cannot prove descent, exactness, or current realization.

## 8. Connection and loss ledger

| field | exact answer |
|---|---|
| source | two marked exceptional chambers `Loc_(3,9)`, then their actual middle-response plane `S` |
| target | the relative response representation `V_rel`; independently, the dual middle quotient `L=Q_middle^vee` |
| map | row lookup `tau_2`, contraction `delta_D`, and quotient by `E_6` as in `(22)` |
| preserved | both reflected response functions, middle orientation, coupled chamber involution, relative zero/nonzero, and the exact two-channel transition |
| destroyed | the `r0=6` endpoint amplitude, outer endpoint rows, absolute lifts, child section, digit chronology, closure edge, source/current semantics, and JC target predicate; ordinary digit-cycle H1 is already zero |
| required sidecar | the actual endpoint line `E_6`, a lawful same-copy closure edge, a physical current-to-response chain map, and a coefficient-compatible filtered JC realization |
| cheapest decisive tests | rank/intersection ledger `(11),(16),(18)`; transition `(24)`; twisted boundary `(30)`; arc-stable hostile `(32)` |
| tournament verdict | intrinsic relation is both-way, so no tournament; XOR erases the load-bearing endpoint mismatch |
| strongest survivor | a minimal chamber-equivariant two-dimensional **relative** cospan, not an `H^1` or flux bridge |

## 9. Exact evidence and current status

The dependency-free companion reads the frozen compact bank, verifies its
hashes, reconstructs `(5)` and all thirteen child quotient profiles, checks
the two-dimensional arc-reversal defect and exact chamber action, and then
proves every rank, intersection, transition, and local-system statement above
by finite-field row reduction.  It includes the scalar hostile, the endpoint
decomposition, the trivial-local-system control, and the arc-stable hostile.

A second implementation imports no primary code and uses
`python-flint.fmpz_mod_mat` for arbitrary-modulus RREF, kernels, inversion,
and multiplication.  It independently reproduces row-sum one, zero digit
seams and constant quotient mode, the rank-five split,
rank-two contraction, intersection dimensions, transition coefficients,
determinant, twisted boundary rank, and middle stabilization.  This is an
independent algebra-engine replay, not yet the requested independent agent
audit of the theorem's typing and scope.

Reproduce with

```text
python -B 04-computation/lrc_r5_rank_two_relative_response_cospan_thm3534_20260816.py
python -B -O 04-computation/lrc_r5_rank_two_relative_response_cospan_thm3534_20260816.py
python -B 04-computation/lrc_r5_rank_two_relative_response_cospan_flint_audit_20260816.py
python -B -O 04-computation/lrc_r5_rank_two_relative_response_cospan_flint_audit_20260816.py
```

Normal and optimized transcripts byte-match the stored output.  Script,
output, and semantic SHA-256 are

```text
4eab23ff77de9d190f13d5a41945415ba55734e0fcdc6f8eae3cbb47d4fb101d
f6972a28ba83ce03c7bfa45a2e6b6eaa8af1504f778a2412ee536639a53f4b6c
57a44964888fdb0a9ca1c890abbe4950c6fa7130f7ab2800c0faa8a5d6a0212d.
```

The FLINT script, output, and semantic hashes are

```text
13316f9add30bb0d5de23d0d0339cf3bdcf666419045a84a3a6e1b96eb54e996
a449fb52016453596ffc581a4df0ccde2028d8777fec6696c122f41ff24bd8d6
f1dfa230027e7fd8f48c7eb7bac2034d0188ee61a4ce6f24d2579860cf532020.
```

The proof candidate constructs no physical current, no `C4` clock closure,
no map to THM-2542's chart class, no JC flux, no row exclusion, and no result
on LRC(14).  LRC(14) remains **OPEN**.
