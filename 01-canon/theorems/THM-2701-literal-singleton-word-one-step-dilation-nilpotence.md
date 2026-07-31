---
id: THM-2701
title: "Literal singleton-word one-step dilation nilpotence and the guard-debt exit"
status: >
  PROVED + VERIFIED-EXACT; BASE THEOREM INDEPENDENTLY AUDITED.  On the canonical
  typed row, let Q_a and Q_b be the literal THM-2305 pure terminal words,
  including guard safety and the unshifted source-c1 safety factor.  Under
  B(y)={13y}, every mixed six-state word in {Q_a,Q_b} is empty, while the
  five-state word bbbbb has a strict exact witness.  The first obstruction
  is structural: an a tooth makes the guard dangerous three states later;
  if the first three states are instead b, the first two b teeth make speed
  14 dangerous at state five.  Thus the literal semantic-word language has
  sharp nilpotence index six and no recurrent SCC.  The nearby BABA orbit is
  instead a lawful guard/source-clock debt cospan: its b vertices are
  guard-danger and its a vertices are unshifted-c1-danger.  This enriched
  edge-debt grammar has an intrinsic period-four base SCC, but no chain map
  to THM-2305 endpoints, target current, scalar row, or LRC(14) conclusion is
  proved.  THM-2710 separately owns the central-half literal-word theorem.
  The exact addenda here show that the current guard/source debt operator span
  is target-null under the proved endpoint covectors; a lawful escape must add
  a paired blocker-graft dipole on common ancestry.  On the THM-2698 half-cycle
  cylinders every provisional paired pivot has all twelve moving
  `C_13` characters nonzero, but the atomic left relation residue is absent,
  so no full endpoint character is yet defined.
source: lrc-semantic-target-bridge-2026-07-28
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
related:
  - THM-2310-quantitative-beta-shift-gluing-of-positive-handoffs
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2693-odometer-skew-product-three-event-escape-and-uniform-delayed-depth-four-nilpotence
  - THM-2698-central-half-odometer-full-local-cycle-and-semantic-sidecar-boundary
  - THM-2710-central-half-phase-literal-word-nilpotence-and-prescribed-clock-invisibility
script: 04-computation/lrc14_literal_singleton_word_nilpotence_thm2701.py
output: 05-knowledge/results/lrc14_literal_singleton_word_nilpotence_thm2701.out
script_sha256: 523618278df6d5762c06598c56bf5c54a0f441dc61a2a0ccc6dc4a316529175b
output_sha256: c9776c5e4ed87d8f1621e9630d5378012a3337c37cd810f606c16192ba028302
hash_basis: LF-normalized bytes
---

# THM-2701 -- literal semantic words die at six; the survivor is an edge debt

**PROVED + VERIFIED-EXACT; BASE THEOREM INDEPENDENTLY AUDITED.**

The phase-scheduled BABA carrier found after THM-2693 has positive rail,
present, private-root, and primitive-unit data.  Its remaining target label
cannot be repaired by attaching another endpoint coefficient.  The literal
semantic target words themselves are not closed under its one-step base
dynamics.

The obstruction is smaller and more informative than the carrier:

```text
literal target-a endpoint -> guard debt three events later;

two initial literal target-b endpoints -> speed-14 failure at event five.
                                                                    (1)
```

Consequently every literal mixed word dies at six states.  The positive
BABA object survives only because it retains the complementary guard sector
and a shifted source-clock complement.  Those labels define a useful new
edge-debt grammar, but not a THM-2305 terminal-word transition.

## 1. Literal word and endpoint typing

Work on the canonical typed row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5).                  (2)
```

For `x in T=R/Z`, write `centered(x)` for its representative in
`[-1/2,1/2)`, and put

```text
I=[-1/14,1/14),
D_v={x:centered(vx) in I},

G=[-1/7,1/7),
C_H=T\G.                                                (3)
```

Thus `C_H` is the guard-safe set and `D_v^c` is the ordinary safe factor.
The common safe base is

```text
A0=C_H intersection intersection_(i=1)^5 D_(q_i)^c.     (4)
```

With source owner `c1`, call `a=c2` and `b=c3`.  The literal THM-2305
singleton terminal words are

```text
Q_a=Q_(1,{a})
   =A0 intersection D_(c2) intersection D_(c1)^c
       intersection D_(c3)^c
   =E_a,

Q_b=Q_(1,{b})
   =A0 intersection D_(c3) intersection D_(c1)^c
       intersection D_(c2)^c
   =E_b.                                                 (5)
```

Both the guard-safe bit and the **unshifted** `c1`-safe bit in (5) are
load-bearing.  A THM-2623 guard-danger sector is outside `A0`.  Replacing
`D_(c1)^c` by a complement translated to a nonzero owner-clock window is a
physical clock refinement, but it is not the terminal word in (5).

Let

```text
B(y)={13y}.                                               (6)
```

For a word `w=w0...w_(d-1)` in `{a,b}^d`, define its literal cylinder

```text
Q[w]=intersection_(j=0)^(d-1) B^(-j)Q_(w_j).             (7)
```

Membership in `Q_a` or `Q_b` is already only the terminal end of a semantic
handoff.  The full THM-2305 source-to-target witness is the span

```text
E_1 <- E_1 intersection B^(-2)E_t -> E_t,                (8)
```

because the source clock is `k_1=2`.  The theorem first proves the stronger
negative statement that even the terminal cylinders (7) cannot recur under
the one-step map (6).

## 2. An a word emits a three-step guard debt

If `y in Q_a`, then by (5)

```text
centered(13^3 y) in I.                                    (9)
```

But

```text
B^3y={13^3y}.                                             (10)
```

Every literal word in `Q_a union Q_b` lies in `C_H`, while

```text
I subset [-1/7,1/7).                                     (11)
```

Therefore

```text
Q_a intersection B^(-3)(Q_a union Q_b)=empty.            (12)
```

This is a literal factor conflict, including the half-open endpoints.  It
does not use the five ordinary safeties, the source-safe factor, the other
target, a rail, a clock label, or an odometer lift.

Equation (12) is the semantic meaning of the guard-danger vertices in the
BABA carrier: a target-`a` tooth cannot be followed three one-step events
later by another guard-safe terminal word.

## 3. Two b words emit a state-five speed-14 failure

Suppose that

```text
y in Q_b,                   By in Q_b.                    (13)
```

Put

```text
z=B^5y.                                                     (14)
```

The two target-`b` factors in (13) say

```text
centered(2z) in I,
centered(26z) in I.                                      (15)
```

The first condition lets us write

```text
z=c+e mod 1,
c in {0,1/2},
-1/28<=e<1/28.                                          (16)
```

Since `26c` is integral, the second condition concerns `26e`.  On the range
in (16),

```text
-13/14<=26e<13/14.                                      (17)
```

No nonzero wrap in (17) can land in `I`: the negative endpoint maps to the
excluded right endpoint `1/14`, and the positive wrap would require the
excluded value `13/14`.  Hence

```text
-1/14<=26e<1/14,
-1/(28*13)<=e<1/(28*13).                               (18)
```

Now `14c` is integral, and

```text
-1/26<=14e<1/26,
1/26<1/14.                                              (19)
```

Thus `z in D_14`.  Since every literal terminal word keeps speed 14 safe,

```text
Q_b intersection B^(-1)Q_b
    intersection B^(-5)(Q_a union Q_b)=empty.            (20)
```

This sharpens the homogeneous target-`b` mechanism to the exact two-hit
factor packet needed here.  It is again independent of every omitted
carrier sidecar.

## 4. Sharp nilpotence index six

Take any word

```text
w=w0w1w2w3w4w5 in {a,b}^6.                              (21)
```

If one of `w0,w1,w2` is `a`, equation (12), shifted to that index, makes
`Q[w]` empty.  Otherwise

```text
w0=w1=w2=b,                                             (22)
```

and (20), using only `w0,w1` and the speed-14 safety inside the sixth word,
makes `Q[w]` empty.  Hence

```text
Q[w]=empty                     for every w in {a,b}^6.   (23)
```

The exact certificate partition of the `64` words is

```text
a at state 0 -> guard failure at 3:    32 words;
a at state 1 -> guard failure at 4:    16 words;
a at state 2 -> guard failure at 5:     8 words;
b at states 0,1 -> speed-14 failure:    8 words.          (24)
```

The bound is sharp.  Put

```text
y0=132661/(2*13^5)=132661/742586.                         (25)
```

Its first six states are

| `j` | `B^j y0` | minimum strict safe-factor margin for `Q_b` |
|---:|:---|:---|
| 0 | `132661/742586` | `186041/5198102` |
| 1 | `18417/57122` | `3314/199927` |
| 2 | `841/4394` | `1493/30758` |
| 3 | `165/338` | `66/1183` |
| 4 | `9/26` | `15/182` |
| 5 | `1/2` | speed `14` has distance `0` |

At states zero through four, `c3 B^j y0` is integral, so the target-`b`
danger factor is at its centre, while every factor required safe in (5) has
the displayed positive margin.  At state five, speed 14 is dangerous.
Therefore

```text
Q[bbbbb] nonempty,                every length-six Q[w] empty. (26)
```

The literal singleton-word language has exact nilpotence index six and no
recurrent strongly connected component under `B`.

## 5. The period-four carrier pays exactly the forbidden debt

Every period-dividing-four tail has the form

```text
y=m/(13^4-1),                  13^4-1=28560.             (27)
```

The dependency-free companion exhausts this lattice with the exact half-open
conventions.  If only the five ordinary safeties and the selected
target-danger/other-target-safe bits are imposed, the complete results are

```text
BABA: m in {1176,27384},
ABAB: m in {13272,15288}.                                (28)
```

They are the two orientations/phases of the one orbit

```text
7/170 -> 91/170 -> 163/170 -> 79/170 -> 7/170.           (29)
```

Adding the literal guard-safe and unshifted-`c1`-safe conditions makes both
rows of (28) empty.  More precisely, on the displayed BABA orientation,

```text
b vertices:
  ||y||=7/170<1/7,                  guard-danger;

a vertices:
  ||13y||=7/170<1/14,               unshifted-c1-danger. (30)
```

These are not two unrelated accidents.  Since `c1=13`,

```text
{c1*y_j}=B(y_j)=y_(j+1).                                (31)
```

The missing `c1` safety at an `a` vertex is literally transported to the
missing guard safety at the following `b` vertex.  Thus the positive orbit
has one edge debt:

```text
b -> a edge: safe debt state;
a -> b edge: dangerous debt state.                       (32)
```

### The intrinsic debt grammar

The scheduled labels can therefore be made intrinsic without pretending
they are terminal words.  Define two Boolean/clock-refined types:

```text
B_debt:
  target b, five unit safeties, c2-safe,
  guard-danger, unshifted-c1-safe, owner clock ell=0;

A_debt:
  target a, five unit safeties, c3-safe,
  guard-safe, unshifted-c1-danger,
  shifted-c1 complement at owner clock ell=1.             (33)
```

Each condition in (33) is a literal physical interval factor.  The
guard-safe/danger pair is exactly the lawful Boolean cospan operation of
THM-2623, and the shifted `c1` complement is retained with its clock label,
not identified with the unshifted factor.  On (29), membership in (33)
determines the alternating type, and `B` gives an intrinsic recurrent SCC

```text
B_debt -> A_debt -> B_debt -> A_debt -> B_debt.           (34)
```

This is the strongest positive survivor.  The external tooth schedule has
become an edge-debt state, but the semantic category has been enlarged.

## 6. Why no inherited action turns the debt SCC into terminal words

The apparent target-switching operations around this carrier have distinct
types:

| operation | what it really moves | why it does not supply `Q_a <-> Q_b` on (34) |
|---|---|---|
| odometer lift `T_m(x)={13x+m/13^6}` | fibre coordinate; tail always follows `B` | changes membership pointwise but not the categorical word span; (23) closes the literal terminal language |
| THM-2365 target co-shift | phases in the fixed two-dimensional target quotient | keeps the delayed word fixed and never permutes target roles |
| THM-2640/2657 carry-root translation | absolute deep root and predecessor carry | is target-neutral on the delayed prefix and moves the present carrier |
| reflection `x -> -x` | oriented root/rail halves | every centered danger/safe factor is reflection-invariant; no target exchange |
| THM-2623 guard complementation | the guard-safe/danger cospan | lawfully enlarges the object, as in (33), but leaves `A0` and hence the THM-2305 endpoint category |
| nonzero owner-clock shift | a translated `c1` complement | is a real clock-labelled factor but is not the unshifted source-safe bit in (5) |
| THM-2315 abstract target swap | names of the two boundary axes | swaps clocks/types (`k_a=4`, `k_b=6`) and the fixed `PAT_QB` carrier; it is relabelling, not a physical map of one fixed packet |
| re-rooted pure-owner span, with THM-2310 mixing if needed | an actual owner `E_a` or `E_b` to a later owner | is the genuine semantic option, but changes the one-step chronology and does not preserve the displayed rail, root, unit, or endpoint phase |

All positive rails in the nearby BABA carrier are built from the one fixed
`PAT_QB` bank.  Alternating the delayed tooth in (33) never changes that
rail's categorical target.  Hence (34) preserves positive physical support
but has no map

```text
{A_debt,B_debt} -> {Q_a,Q_b}                             (35)
```

whose image preserves membership, much less a chain map between the
corresponding THM-2305 witness spans.

There is also no existing target-current differential which can absorb the
edge debt.  In THM-2461's exact covector table, both the guard and source
owner are target-neutral.  Pulling a polarized THM-2379 repair through a
zero target covector gives the zero target family; it cannot transport the
guard/source debt.  THM-2623 supplies the positive cospan identity

```text
Q_safe + Q_danger = Q_(omit guard),                       (36)
```

not a signed endpoint differential.  A mapping-cone attempt would have to
retain `Q_danger-Q_safe`, prove covariance on one common ancestry, and map
the fixed `PAT_QB` rail to the correct terminal endpoint.  None of those
three steps is in canon.

The concurrently verified four primitive-unit rows add a useful invoice,
not a repair.  Conditional on identifying their four separately selected
rows along one legal chain map, their product in

```text
F_13[z]/(1+z+...+z^6)
```

is

```text
9+5z+5z^2+12z^3+11z^4+9z^5,                            (37)
```

with multiplication determinant `3` and exact order `168`.  Thus a naive
multiplicative transport would acquire nontrivial holonomy.  No identification
with THM-2542's `C_91` chart holonomy or THM-2607's boundary invoice is made:
their transition types have not been connected to (33).

The two honest exits are therefore

```text
edge-debt exit:
  build a signed guard/source mapping cone and a common-ancestry endpoint
  chain map which explains or cancels (37);

chronology exit:
  re-root at the actual terminal owner and use its clocks 4 and 6, with
  new rail/present/root/unit transport.                              (38)
```

Another positive endpoint search on the same one-step literal-word category
cannot work, by (23).

## 7. Exact verification and scope

Run

```bash
python3 04-computation/lrc14_literal_singleton_word_nilpotence_thm2701.py
python3 -O 04-computation/lrc14_literal_singleton_word_nilpotence_thm2701.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_literal_singleton_word_nilpotence_thm2701.out.
```

The companion is dependency-free and uses explicit optimized-mode guards.
It checks the half-open wrap constants in Sections 2--3, classifies all `64`
length-six words by (24), verifies the strict five-state witness and every
margin in Section 4, and exhausts all `28,560` period-four lattice points in
Section 5.

An independent audit rederived the literal-word typing and both factor-sparse
implications, including the half-open wrap boundary.  It replayed normal and
optimized executions, byte-matched the frozen output, and reproduced the
declared hashes

```text
script: 523618278df6d5762c06598c56bf5c54a0f441dc61a2a0ccc6dc4a316529175b
output: c9776c5e4ed87d8f1621e9630d5378012a3337c37cd810f606c16192ba028302
```

The theorem is confined to the canonical typed row and the one-step base map
`B`.  It does not exclude a changed chronology, an enriched signed debt
category, or a different row.  The row is a typed non-cover control; no one
of the `165` scalar obligations is removed, and LRC(14) remains open.

## 8. Central-half route: canonical ownership

THM-2710 proves and audits the central-half facts previously recorded here:
`B_1^k=B^k` at the prescribed clocks `k=2,4,6`, every six-state literal
`B_1` word is empty, and the half and transverse packet cycles fail the
unshifted-`c1`-safe and exclusive-source endpoint predicates.  It also records
the genuine odd-delay predicate repair and the doubled reflection-sheet debt.
Those proofs and their exact companion are not duplicated here.  Sections
9--10 use only the resulting location of the open problem: an enriched debt
carrier, not another literal endpoint word.

## 9. Target-null debt-cone addendum

The vague mapping-cone exit in Section 6 has a sharp boundary.  In role order

```text
(H,j,k_a,k_b,a,c),
```

the proved target-covector matrix and the two chronological debt edges are

```text
Lambda = [0 0 -1  0 1 0],      rho = [-1  1]
         [0 0  0 -1 0 1]              [ 1 -1]
                                       [ 0  0]
                                       [ 0  0]
                                       [ 0  0]
                                       [ 0  0].          (39)
```

Although `rank(Lambda)=2`, exact multiplication gives

```text
Lambda rho=0.                                             (40)
```

Thus every signed combination generated by guard complementation and
source-clock toggling is target-null.  This remains true before ancestry is
marginalized: THM-2379's common-ancestry table satisfies
`K_f^+(r,0)=0`, while the actual covectors of `H` and `j` are zero.  The
pullback therefore vanishes on the full `156`-dimensional anchored table
space.

The smallest formal nonzero columns are `k_a` or `k_b`; the smallest lawful
directions are the paired dipoles

```text
D_a=e_a-e_(k_a),              D_c=e_c-e_(k_b),            (41)
```

which map to `2I_2`.  A graft complement alone does not realize `(41)` on
the debt ancestry.  The next decisive test must insert the paired
blocker-graft table on one exact debt cylinder and retain the left endpoint
residue, not only the moving residue.

The exact secondary companion is

```text
04-computation/lrc14_guard_source_debt_cone_target_span_probe_20260728.py
05-knowledge/results/lrc14_guard_source_debt_cone_target_span_probe_20260728.out
```

with SHA-256 values

```text
172981c6f60fbf062add30aff419d61f4b11f305710a9518d3794f1ef743efff
bcecc1c193a99dba2ba5ac24713bd4cde2fe298563c1d08998c7d2ac8fd63de5
```

This no-go is confined to the operator span generated by the proved debt
operations and endpoint covectors.  It does not forbid a new paired dipole,
an actual common-ancestry chain map, or re-rooting at a different source.
Together with THM-2710 it forbids another literal-endpoint search on the same
`B` or `B_1` categories.

## 10. Paired-cylinder addendum: the moving characters survive

The paired-dipole test in Section 9 has a positive partial answer on each
THM-2698 event cylinder.  Use

```text
a=13^3,              c=2*13^5,
(x0,r0)=(55232507/(24*13^6),7),
(x1,r1)=(58313459/(24*13^6),3),
epsilon=1/301082946198216.                              (42)
```

For every unresolved pivot `k in {27,40,53,66}`, define

```text
K_(i,k)(s)=integral_(|x-x_i|<epsilon)
 d_1(cx-r_i/13) u_1(ax-s/13) u_1(kx+s/13) dx.           (43)
```

All three factors are constant on the interval.  Exact wall arithmetic gives

```text
K_(i,k)(s)=2 epsilon * 1_(s in S_(i,k)),                 (44)
```

where every one of the eight support sets `S_(i,k)` is a nonempty proper
subset of `F_13`.  For every `h!=0`, its Fourier character is nonzero:

```text
sum_(s in S_(i,k)) zeta_13^(h s) !=0.                   (45)
```

Indeed a vanishing sum would make the zero--one polynomial
`sum_(s in S)z^s` divisible by `Phi_13=1+z+...+z^12`; since its degree is at
most twelve, all thirteen coefficients would have to be equal, contrary to
`S` being nonempty and proper.  The exact companion also checks all `96`
cyclotomic remainders directly.

Thus numerical moving-endpoint cancellation is not the next obstruction.
For `eta=e_a-e_(k_a)`, `(45)` supplies the moving term `-eta.v`.  A THM-2334
endpoint character requires `eta.(u-v)`, but the half-cycle packet retains no
atomic left-present relation index `u`.  Its rail/root source labels are not
relation indices.  Hence `eta.u` is undefined, not zero.  The pivot identity
among four candidates is also unresolved, and the old-head `k_a`-danger zero
plane used by THM-2563 is absent.

The exact companion is

```text
04-computation/lrc14_half_cycle_paired_debt_cylinder_probe_20260728.py
05-knowledge/results/lrc14_half_cycle_paired_debt_cylinder_probe_20260728.out
```

with SHA-256 values

```text
87f6fa3b229a888ce05518b498802a6774bd0070835eba97b6ecbd560d4d596f
24c3dc756e1521e47244f0cdc0827b22fdaaab9125afbd9257a7d094c7c93ef1
```

Normal and optimized executions byte-match.  The theorem does not identify
the `B_0` BABA debt SCC with the `B_1` half-cycle ancestry; it proves only the
paired moving-character statement on the latter.  The next object must retain
both atomic relation indices `(u,v)` on one ancestry, then test
`eta.(u-v)` and full-`X` aggregation.

The theorem is confined to the canonical typed row.  No one of the `165`
scalar obligations is removed, and LRC(14) remains open.

QED.
