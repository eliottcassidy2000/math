# A 28-coordinate exact observer for the cap-two language

**PROVED BY EXACT FINITE LINEAR CERTIFICATE:** the degree-cap-two adaptive
language admits six control-dependent rational query charts with dimensions
3,3,1,3,11,7, replacing its 460 path-cover coordinates. Each chart dimension
is sharp among rational linear observers at that fixed control. Their direct
sum gives an exact 28-dimensional representation. For the guarded word
series assigning zero to illegal active-rank words, dimension 28 is also
globally sharp. No corresponding minimum is claimed when illegal outputs
are unconstrained, or among nonnegative-matrix representations.

**PROVED:** the stable control (2,2,1) has two return commands with different
scalar behavior. Repeating rank 0 gives an affine power-of-two sequence.
Repeating rank 1 gives an irreducible degree-eight recurrence, minimal for
every actual source at that control. **OPEN:** the cap-three minimum observer,
constructive entry into all prime and 7p targets, and the numerical image of
any proposed mixed return semigroup.

## 1. Inheritance and the board

The closest proved mechanism is the temporal endpoint-pairing carrier in
[ACTIVE_PORTS.md](../adaptive_metagraph_reformulation_20261103/ACTIVE_PORTS.md).
It has 460 reachable (degree-control, endpoint-signature) coordinates for
cap two and 9,672 for cap three. These were explicit positive carrier sizes,
not minimum linear dimensions. The cap-one reduction in
[ADAPTIVE §2](../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md)
has the sharp three-dimensional positive representation

    R(a,q,s)=(2a,q,q),   N(a,q,s)=(a+q,a+s,a).

The canonical hostile is the same note's actual order-six pair with equal
(n;h,s,q,F)=(6;37,5,13,71) but different older-root continuation spectra.
The corrected near miss is to infer an observer dimension from the count of
reachable forest coordinates, or from one modular rank. The least-used
sidecar in this round is the distinction between the **linear span of actual
reachable sources** and the much larger formal forest-coordinate space.

The incoming [inverse65 inverse](../inverse65_zeckendorf_20261103/INVERSE.md)
and [Zeckendorf audit](../inverse65_zeckendorf_20261103/ZECKENDORF.md)
supply two relevant boundaries: positive integral coordinates need not be
reachable, and a return's characteristic polynomial imposes an obligation
before it can implement a numeral shift. Neither packet supplies a cap-two
rank result or the return recurrence below. No literature-priority claim is
made.

| Live concept | Hypothesis and operation | Cheap test / result |
|---|---|---|
| Anchor: smaller exact observer | Quotient the actual reachable carrier by future-H invisibility | A modest cap-two modular probe suggested 28 chart coordinates; exact closure proves the upper bound. |
| Source versus formal forest | Remove relations valid only on reached source vectors | Full-carrier observable ranks are larger; the 28-chart theorem is source-language-relative. |
| Return arithmetic | Compare the two self-loops at control (2,2,1) | One output recurrence has degree two; the other has irreducible degree eight. |
| Integral entry | Treat positive odd query values as a possible source | A mod-75 necessary condition rejects a small positive odd fake state. |
| Cap-three extension | Transfer the same minimization method | Still open; no second large rank or support census was run. |

The perspective generator was `actual source words × future queries ×
reachable/observable quotients × root-rank controls × return spectra ×
integer lattices`. Its three rounds were rank compression, the return
polynomial, and the source-lattice hostile.

## 2. The actual alphabet and the exact certificate

All graphs start from the named C3 with outgoing rows (2,4,1), followed by
one birth at root 0. A further letter j chooses the j-th currently active
vertex, ordered by birth name; active means current indegree at most two.
The newborn loses only to that named root and beats every other old
vertex. Every word and suffix below uses these dynamic **active ranks**,
not immutable vertex names. The literal construction retains the name list.

The six controls are the tuples of active indegrees

    11, 121, 1221, 211, 221, 2221,

where, for example, 221 means (2,2,1). The initial control is 1221 and its
initial graph has H=5. The prior theorem gives an exact integer carrier
for every legal word, with finite matrices M_j between these controls.

For each control c, the certificate chooses actual prefix words whose
carrier rows R_c are rationally independent. It then verifies, over Q,
that every row of every R_c, advanced by every legal letter j, belongs to
the destination span R_d. The initial source is included. This proves by
induction that the listed spaces contain every actual reachable source,
at every order; no finite word-length cutoff is used for that conclusion.

For a legal continuation w define the physical observation

    Q_w(T)=H(T after the actual continuation w).

The empty continuation observes H itself. A selected list of these
observations gives a row vector y_c(T). The certificate supplies a rational
matrix L_(c,j) with

    y_d(T after j)=y_c(T) L_(c,j).                    (2.1)

It verifies (2.1) on every row of R_c and hence on its whole span. Since H
is the first observation in every chart, equal query vectors at the same
control give equal H after every legal future continuation. Each coordinate
is an actual nonnegative path count. The transition matrices may have
negative or fractional entries; they are used only on reached states.

Here are the exact dimensions and one nonzero actual Hankel minor per
control. Each Hankel entry is H after an actual prefix followed by an actual
suffix; all prefixes in a row block end at the displayed control.

| Control | Forest coordinates | Reachable Q-span dimension | Exact query dimension | Hankel determinant |
|---|---:|---:|---:|---:|
| 11 | 15 | 4 | 3 | 4 |
| 121 | 50 | 4 | 3 | -8 |
| 1221 | 23 | 1 | 1 | 5 |
| 211 | 42 | 4 | 3 | -4 |
| 221 | 73 | 19 | 11 | 917504 |
| 2221 | 257 | 9 | 7 | -9600 |
| Total | 460 | 41 | 28 | — |

Every rational linear observer at control c factors its actual Hankel
matrix through its state space. The nonzero determinant therefore gives
the matching lower bound for that chart. This proves local sharpness.

There is also a separate global lower bound, with its guard stated
explicitly. Extend the actual H function to every word in {0,1,2,3} by
assigning value zero as soon as a letter is illegal at its current control.
The block direct sum of the six charts represents this total word series
in dimension 28. Its 28 selected actual prefix rows and 28 suffix words
of length at most three give the exact Hankel determinant

    12470867040665600 ≠ 0.

The square matrix, including its zero entries for illegal concatenations,
is recorded in `global_guarded_hankel` in the exact receipt. Thus dimension
28 is the minimum rational linear dimension of this **guarded** series.
The argument does not apply to a representation permitted to assign
arbitrary values to illegal words. Nor is it a positive-realization theorem.

For comparison, the first modular probe found full **formal-carrier**
observable dimensions 12,27,21,24,20,31. Thus simply dropping unobservable
forest columns before imposing actual reachability would miss much of this
compression. Those modular numbers are retained only as diagnostics; the
exact chart theorem uses the rational closure and exact minors.

## 3. The continuation-count charts

The selected suffix lists are:

| Control | Query suffixes, in coordinate order |
|---|---|
| 11 | empty, 0, 1 |
| 121 | empty, 0, 1 |
| 1221 | empty |
| 211 | empty, 1, 2 |
| 221 | empty, 2, 0, 1, 02, 12, 10, 11, 112, 110, 111 |
| 2221 | empty, 3, 0, 1, 2, 02, 12 |

Concatenation denotes successive active-rank letters. All prefixes,
Hankel matrices, selected coordinate pivots, and all 19 transition matrices
are frozen in `active_rank_exact.json`. The source-carrier hash is included.
The initial query vector is the one-coordinate vector (5).

The small queries already distinguish mathematical information from an
arbitrary state encoding. For example, at control 221 the three counts
H,H_0,H_1 are insufficient as a rational linear observer for all later root
choices: the actual Hankel rank is eleven. Conversely eleven actual continuation counts are sufficient
for every later cap-two word. The statement is exact on this source
language, not on every strong tournament with the same degree tuple.

The observer is a useful input for induction searches because each trial
word has a fixed rational matrix and each coordinate has a literal
construction meaning. It still requires an actual reachable initial query
vector; the next section gives an explicit entry guard.

## 4. A lattice guard: positive odd coordinates do not ensure entry

At control 2221 write its query coordinates as

    y=(H,H_3,H_0,H_1,H_2,H_02,H_12).

The exact root-0 transition, followed by another rank 0, gives

    H_00 = (-236H+299H_0−40H_1−12H_2−29H_02+37H_12)/75.  (4.1)

Thus the numerator is divisible by 75 on every actual cap-two source at
this control. This is a necessary coordinate-lattice condition, not a
complete list of congruences or a classification of the reachable image.

The actual one-letter prefix (1) has

    y=(13,23,29,25,29,49,43),   H_00=59.

Changing only H to 15 gives a positive odd tuple

    (15,23,29,25,29,49,43),

but (4.1) becomes 3953/75. It cannot be an actual source query vector.
The scalar 15 is attained in the adaptive language; its attainability
does not authorize this state replacement. This is the same logical
boundary as the incoming inverse65 cone hostile, now with an exact
cap-two continuation-coordinate obstruction.

## 5. The stable control has two different return spectra

At control 221, letters 0 and 1 both return to 221. Letter 2 exits to 11.
Their 11-by-11 query matrices A_0,A_1 are integral, even though some
inter-control matrices are rational. The certificate checks the following
polynomial identities on the H-output column e_H:

    (A_0²−3A_0+2I)e_H=0,
    p(A_1)e_H=0,

    p(X)=X^8−2X^7+X^6−4X^5+3X^4−6X^3+2X²−4X+1.    (5.1)

Therefore, for any actual source at 221, repeated rank 0 has

    H_k=H_0+(2^k−1)(H_1−H_0),                       (5.2)

where H_1 here means the count after one rank-0 return. This is the
fixed-protected-root ray. In contrast, repeated rank 1 has

    H_(k+8)=2H_(k+7)−H_(k+6)+4H_(k+5)−3H_(k+4)
              +6H_(k+3)−2H_(k+2)+4H_(k+1)−H_k.   (5.3)

Both identities hold from k=0, without an eventual-tail qualification.

The polynomial p is irreducible over Q. The self-contained checker proves
its irreducibility modulo 31 by the degree-eight Rabin criterion:

    X^(31^8)=X mod p,
    gcd(p, X^(31^4)−X)=1 mod 31.

Since the only prime divisor of eight is two, these are the required
conditions. The reduction is monic of unchanged degree, so irreducibility
lifts to Q. For every actual source the H sequence is nonzero, and its
minimal recurrence polynomial divides the irreducible annihilator p.
Consequently (5.3) has **minimum rational recurrence degree eight for
every actual source at this control**. In particular no such complete
unary return sequence is an affine power-of-two or Fibonacci recurrence.
This does not obstruct mixed words, other controls, or selected isolated
counts.

There is also an exact integral return register: retain the eight consecutive
observations (H_k,...,H_(k+7)). Equation (5.3) shifts this register by an
integer companion matrix of determinant one. Its inverse recovers H_k
integrally because the coefficient of H_k in (5.3) is -1. Thus the register
transition is invertible modulo every modulus. This preserves an arithmetic
clock, but does not reconstruct the discarded source graph or guarantee
that an arbitrary integral register is realized. It is the same distinction
between a numeral inverse and legal source entry emphasized by the incoming
inverse65 work.

The actual active-rank prefix (1,0), corresponding to immutable-root
history (0,1,0) from C3, reaches H=29 at control 221. Its rank-0 ray is

    30·2^k−1,

while its rank-1 tail begins

    29,69,155,349,803,1837,4205,9633,22031,50409,115361,263957.

No infinitude-of-primes assertion is made for either sequence. The new
point is the exact return mechanism and its minimal recurrence, not the
individual numerical attainments.

There is a useful structural explanation for the asymmetry. The active
ports at 221 consist of the protected degree-two root, the previous
newborn now of degree two, and the newest degree-one vertex. Rank 0 keeps
the protected root; rank 1 discards it and shifts protection to the other
degree-two port. Under repeated rank 1, each later newborn selects the
vertex born two steps earlier. This moving boundary carries more future
information than the fixed-root ray.

## 6. Connection contract and the remaining positive obligation

| Source → target / map | Preserved predicate | Lost data and required sidecar | Decisive control |
|---|---|---|---|
| Actual cap-two forest vector → query chart | Every future H under legal cap-two root ranks | Graph, automorphisms, arbitrary attachment response; retain source word and degree control | Exact full-carrier reachable closure and actual Hankel minors |
| Query tuple → attempted inverse state | Literal continuation counts when reachable | Positivity/parity do not impose the source lattice | The H=15 fake tuple in §4 |
| Return matrix → scalar recurrence | Whole unary H sequence from the actual source | Other commands and state entry | Irreducible degree-eight output polynomial |
| Cap-two minimization → cap-three problem | Same method applies in principle | No cap-three rank result has been computed here | Cap two still permanently misses 127 |

Cap-two's inherited permanent exclusion of 127 survives observer
compression. The smaller state does not enlarge its language. A useful
prime/7p induction must either prove entry in a stronger actual language
(for example cap three), or specify source-changing operations and their
new observer. A return matrix alone does not establish entry, and a
formal state satisfying (2.1) is not a source. The successful move here is
to replace 460 coordinates by explicit physical queries with exact closure
and exact lower bounds, then test concrete return actions on that state.

## 7. Reproduction and independent audit

`active_rank_probe.py` uses one modest modular closure experiment over
1000003 to select candidate bases. Its standalone receipt
`active_rank_modular.json` explicitly labels those ranks as modular lower
bounds only. The source is the proved prior active-port transfer, loaded
by an explicit relative path and pinned by SHA-256.

`active_rank_exact.py` reconstructs the integer source vectors, chooses
pivots, inverts square matrices with standard-library exact fractions,
and checks every coordinate of every claimed relation over Q. It verifies
129 full-coordinate reachability transitions and 947 scalar query closure
equations, all six local nonzero Hankel minors, the global guarded
28-by-28 minor, the two output annihilators, and the elementary irreducibility
certificate. Thus no rational upper bound is inferred from a modular rank.
Runs take about 0.7 seconds locally.

```sh
python3 -I -B 05-knowledge/results/adaptive_state_induction_20261103/active_rank_probe.py --cap 2
python3 -I -B 05-knowledge/results/adaptive_state_induction_20261103/active_rank_exact.py
python3 -I -B -O 05-knowledge/results/adaptive_state_induction_20261103/active_rank_exact.py
```

The exact normal and optimized outputs agree with `active_rank_exact.json`.
The independent `metagraph_reset_rank_audit.py` imports none of the producer
implementation: it reconstructs actual named graphs, counts Hamilton paths
by recursive completion, checks all 198 Hankel entries and their Bareiss
determinants, and tests every one of the 19 transfer matrices on an actual
source. It also checks all 784 entries of the global guarded minor,
including 176 illegal-word zeros, and independently obtains its determinant
12470867040665600. The combined audit counts 433 distinct actual graphs,
including the first five values of the rank-1 tail and the mod-75 entry
guard. The independent audit validates the physical queries and minor
values; the all-order upper bounds additionally use the exact linear-closure
proof described above.
