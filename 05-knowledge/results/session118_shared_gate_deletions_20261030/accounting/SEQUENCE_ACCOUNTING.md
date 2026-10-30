# Exact sequence inventories and the coefficient-four boundary

**Status: PROVED accounting for specified legal transitions.** The statements
below do not certify the external candidate's full gate-elimination theorem or
its explicit Boolean family. They apply to any sequence for which the stated
transition bounds and invariants have separately been justified. In particular,
abstract scalar schedules below are not claimed realizable by a circuit.

The reference is Section 7 of the pinned [ternary-pair candidate](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/4a61f024a4993592d97b49ca7d00979029b7fece/submission/elimination.md),
read in Session 117. We also give the exact consequence of strengthening AP2's
three deletions to four. The latter conclusion is conditional on the separately
audited structural lemma, not proved by this accounting.

## 1. Variables and complete slack bookkeeping

Start with n parameters, no marked pairs, and I_0 used unpaired parameters.
Write J=n-I_0>=0. At the end there are m_T parameters, q_T disjoint marked
pairs, I_T used unpaired parameters, and s_T charged gates/equations. Set
h=a_T=m_T-q_T. No dimension is silently replaced by the number of used inputs.

Use these nonnegative integer row counts:

- G: GP steps, with child counts k_j>=1;
- B: AP2 steps;
- A: AP1 steps;
- P: P2 and PE steps combined;
- Z: P0 steps;
- H: PC steps;
- D: elementary D1--D4 deletions.

Preparations M and AR preserve m,q,s and do not increase I. Let U be the
cumulative additional decrease in I beyond the following upper bounds:
GP: -2; AP1: -1; P2/PE: +1; every other row: 0. Thus U>=0 is defined by
actual changes, including any usage loss in preparations or elementary
simplifications. It counts lost used parameters, not deleted gates.

Put K=sum_GP(k_j-1)>=0. Let E>=0 be the total number of charged-object deletions
beyond each row's guaranteed count, including any further losses in P0 or
preparation. Each actual deletion is assigned once to its reducing row or
preparation. Shared descendants must already have been counted distinctly in
the structural row lemma; this algebra never adds two deletions of one gate.

For the original table the guaranteed charged-object counts are
GP: 2+(k_j-1); AP2: 3; AP1: 3; P2/PE: 2; PC: 1; D: 1; P0/M/AR: 0.
For the strengthened table replace only AP2: 3 by AP2: 4 and redefine E
relative to that stronger guarantee. Write r=3 or 4 for this AP2 count.

## 2. Exact global conservation identities

For either table,

    n-h = G+B+A,
    q_T = G-2B-A-P-Z-H,
    I_T = n-J-2G-A+P-U.

The first is safe-core dimension conservation. The second counts pair births
and deaths. The third is the exact usage equation with its nonnegative slack.
Eliminating G and P gives the particularly useful identity

    h = B+A+q_T+Z+H+U+I_T+J.                       (1)

Every term on the right is nonnegative. In particular

    #AP2 + #AP1 <= a_T,
    G >= n-2a_T,
    P >= n-4a_T.

For the last inequality, substitute G=n-h-B-A and the pair inventory:
P=n-h-3B-2A-q_T-Z-H. The subtracted expression after h is at most 3h by (1).
Negative lower bounds are harmless. If the algorithm terminates at h=d+1,
all AP1/AP2 steps together occur at most d+1 times, whereas all but O(d)
dimension-reducing steps are GP. At least n-4(d+1) pair-cleanup steps P2/PE
are needed. These are restrictions on the full sequence, not independent
local worst-case choices.

**Important consequence:** improving AP2 by one deletion can improve the
boundary/error term by O(d), but cannot alone change the leading coefficient
when d=o(n). The abundant candidate bottleneck is GP followed by P2/PE,
not a positive-density population of AP2 steps.

## 3. Exact gate-slack identities

The actual charged-object inventory is

    s_0 = s_T + 2G+rB+3A+2P+H+D+K+E.             (2)

Substituting (1) into (2) yields the following exact identities.

### Original three-deletion AP2 table

    s_0 = 4n-9h
          + s_T + 5(I_T+J+U) + 3q_T
          + 2A + 3Z + 4H + D + K + E.           (3)

### Strengthened four-deletion AP2 table

    s_0 = 4n-8h
          + s_T + 4(I_T+J+U) + 2q_T
          + A + 2Z + 3H + D + K + E.            (4)

Thus the new table gives 4n-8(d+1), conditional on full structural coverage,
when the terminal dimension is at most d+1. Neither identity by itself proves
that coverage or constructs an eligible disperser. The quantities in these
identities are exact and do not presume fanout additivity.

These formulas also classify near-tight traces: if the surplus over the base
term is S, each displayed nonnegative summand is at most S. In particular the
number of GP children beyond one, extra deletions, and elementary deletions
have total at most S. A linear positive lower bound on any surviving surplus
term would give a genuine leading-coefficient improvement. Merely knowing
that AP2 has bounded frequency does not supply one.

## 4. Complete aggregate equality classification

For either (3) or (4), equality with its displayed base term holds if and only if

    s_T=I_T=J=U=q_T=A=Z=H=D=K=E=0.

Equations (1) and the inventories then force

    B=h, G=n-2h, P=n-4h.                         (5)

Hence n>=4h is necessary. Every GP must have exactly one child, each charged
row must attain its minimum deletion count, every P2/PE must attain its usage
increase of one, all starting inputs must be used, and no extra cleanup can
occur. The mix of P2 versus PE is not determined by scalar accounting.
Conversely, these conditions and (5) make the scalar identity an equality.
This converse concerns only the inventories; it does not assert a legal DAG,
a terminal output, a function, or a source-theorem equality example.

## 5. Prefix-feasible abstract schedules: why inventories alone stop here

For integers h>=1 and n>=4h consider the word

    (GP,P2)^(n-4h), GP^(2h), AP2^h.              (6)

Use k=1, exact usage bounds and exact minimum charged-object losses throughout.
Start with m=I=n, q=0 and s=4n-(12-r)h. Each first two-step block has net change
(m,q,I,a,s)=(-1,0,-1,-1,-4). After these blocks, m=I=4h and s=(4+r)h: this is 7h for r=3 and 8h for r=4.
The 2h GP steps leave m=4h,q=2h,I=0,s=rh,a=2h. The last h AP2 steps leave
m=h,q=I=s=0,a=h. Throughout, m>=2q, 0<=I<=m-2q, s>=0, and every required
pair/unpaired inventory exists. Before the last step a>=h+1.


This shows that pair/usage conservation, nonnegativity, the prefix domain-size
constraints and the loop's numerical dimension condition do not improve
coefficient four or the corresponding scalar boundary penalty. However, (6)
has **not** been realized as a reduced B2 DAG. In particular, absence of a
nonconstant terminal output might refute an equality schedule without refuting
its scalar inventories. A stronger theorem must use retained incidence,
output ancestry or another structural invariant that these scalars omit.

The shared-gate distinction matters: a proof that GP/P2 blocks force losses
must show that proposed extra deletions are distinct gates, cannot be re-used
by later blocks, and survive affine role exchange and pair-domain restrictions.
This inventory calculation identifies exactly which high-frequency blocks
such a proof must target, and how much certified surplus would translate to
an improved lower bound.

## 6. Reproduction

`check_sequence_accounting.py` checks exact integer identities on a bounded
aggregate grid and verifies all prefixes of the schedule (6) for both tables.
It uses explicit exceptions and runs unchanged under Python optimization.
The all-parameter assertions above follow from the displayed identities;
finite tests are regression checks, not universal DAG coverage.
