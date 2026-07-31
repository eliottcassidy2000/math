---
id: THM-2882
title: "Event-twisted all-q coefficient carry lift and selector-parity obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Zero-chart truth plus u3 give a charged C169 lift: reduced/vertical
  omega^3, closed flat.  External frequency acts by diag(E,1), not on the
  raw scalar.  q7/q11 need independent origin parity with no typed QA/QAB.
  No row exclusion or LRC(14) proof follows.
source: root/lrc-event-twisted-coefficient-2026-07-29
depends_on:
  - THM-2835-q11-semantic-word-horn-and-bockstein-blind-support-no-go
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent
  - THM-2874-endpoint-kummer-galois-clutch-and-bockstein-seam-transgression
  - THM-2878-endpoint-factor-exit-carry-transducer-and-flat-vertex-boundary
  - THM-2880-q0-q3-one-fibre-selector-provenance-obstruction
related:
  - THM-2859-horn-collar-q0-hinge-minimal-v4-globalization-and-witt-endpoint-obstruction
script: 04-computation/lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.py
output: 05-knowledge/results/lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.out
script_sha256: 3ed346e0c631b34bd61f0c4d27d7f161e8d35b70decfb95f5207c5f57893d005
output_sha256: 0faa0a24f6ba8b6c88b6bbfc4f225e38667097b1a937d977741453499884901c
hash_basis: LF-normalized bytes
---

# THM-2882 -- Event-twisted all-q coefficient carry lift

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem proves no row exclusion and no LRC(14) conclusion.  It
constructs an exact coefficient-level `C169` lift from a literal endpoint
event, and then identifies the independent selector coordinate still
missing from the physical current.

## 1. Verdict

There is a genuine positive coefficient construction.

1. At the zero endpoint character, truth-polarizing by `PAT_E3` or its
   full complement gives one full translated target atom at every
   `q in C13`.
2. The THM-2868 Prony split on each atom gives a charged line
   `U_(r+1)=omega^3 U_r`.
3. Decorating the base-13 chart `n=13a+q` by section
   `r=a+q-3 mod13` and transporting `a` by THM-2878's literal `u3` exit
   realizes the nonsplit carry law on all 169 states.
4. The reduced q3/q11/q7 seam has character-three curvature `omega^3`,
   exactly the vertical ancestry clutch.

There is an equally exact negative.

- At fixed frequency, raw target translation has ratio one.
- The twist acts on the charged Prony channel by `diag(E,1)`, not scalarly
  on the recombined coefficient.
- Restoring the two endpoint origins makes every common-selector q11/q7
  coefficient cancel.
- A nonzero q7/q11 coefficient requires an origin-odd selector parity
  which is independent of the carry event and has no typed source/QA/QAB
  operation.

Thus the abstract event-decorated rank-one coefficient transition is
complete; the physical clutch is not.

## 2. Thirteen truth-polarized atoms

Let `J_q` be the translate by `qT/13` of the selected target atom.  At the
zero endpoint character, its complete literal word is

```text
q0  SSSSSSSSD       q7  DSSSSDSSD
q1  SSDSSSSSD       q8  DSSSSSSSD
q2  SSDSSSSSD       q9  SSSSDSSSD
q3  SSSSSSSSD       q10 SSSSDSSSD
q4  SDSSSSSSD       q11 SSSSSSSSD
q5  DDSSSSSSD       q12 SSSDSSSSD
q6  DSSSSDSSD.
```

The full `PAT_E3` contains `J_q` exactly at

```text
q in {0,3,11};
```

its set-theoretic complement contains every other `J_q`.  Choosing the
true block therefore produces thirteen full, distinct, translated
weight-one atoms.  This is a zero-chart Boolean coefficient atlas, not yet
a common two-origin selector.

All thirteen atoms have the same endpoint parameters and common Prony
nodes.  Using the inherited formal sections

```text
m_r=1+42r
```

and offsets `e_4=1`, `e_8=2`, all other `e_r=0`, gives 26 adjacent samples
and 26 corresponding endpoint frequencies, all units modulo `91`.  In
both certified fields, every atom/section split satisfies

```text
U_r=U_0 omega^(3r),       V_r=V_0,                    (1)
```

with both channels nonzero.  The common source scalar `P_src` occurs
algebraically in both channels; this reuse does not itself construct an
all-q fixed-source packet operation.

## 3. The nonsplit 169-state lift

Write each `n in C169` uniquely as a **set chart**

```text
n=13a+q,             0<=a,q<13.                       (2)
```

Equation `(2)` is not a product-group identification
`C169=C13 x C13`.  The carry witness

```text
L_1(0,12)=(1,0)
```

is precisely the nonsplitting.

Assign state `(a,q)` the section and decorated vertex

```text
r(a,q)=a+q-3 mod13,
C_(a,q)=(J_q,U_(r(a,q))).                              (3)
```

These 169 vertices are pairwise distinct.  Let

```text
kappa(q,h)=floor((q+h)/13),       0<=q,h<13,           (4)
```

the positive-path count of the address-uniform `u3:D->S` event.  Define

```text
L_h(a,q)=(a+kappa(q,h), q+h mod13).                    (5)
```

The exact charged coefficient ratio is

```text
E(q,h)
 =U_(r(L_h(a,q)))/U_(r(a,q))
 =omega^(3h) omega^(3 kappa(q,h))
 =F(q,h) tau(q,h).                                    (6)
```

The first factor is the old flat frequency transition; the second is the
event twist.  The companion checks all 169 reduced edges and all 2,197
lifted state transports in each field.

For `s=h+k mod13` and `epsilon=floor((h+k)/13)`,

```text
kappa(q,h)+kappa(q+h mod13,k)
 =kappa(q,s)+epsilon,                                  (7)

E(q+h mod13,k)E(q,h)
 =omega^(3 epsilon) E(q,s).                            (8)
```

All 2,197 reduced and 28,561 initial-state versions of `(8)` hold.  Thus
`E` is a projective coefficient system on the reduced residue groupoid and
an honest system on the nonsplit lift, where `omega^(3 epsilon)` is the
vertical `T` transition.  With

```text
phi(q)=omega^(-3(q-3)),
```

all 169 edges obey

```text
tau(q,h) phi(q)=phi(q+h mod13) E(q,h),                 (9)
```

so the event-twisted line is gauge-equivalent to the ancestry character.

## 4. Curvature is not closed-loop holonomy

For the reduced seam

```text
q3 --8--> q11 --9--> q7
q3 --------4-------> q7,
```

the event counts are `(0,1,0)`, and

```text
E(11,9)E(3,8)/E(3,4)=omega^3.                         (10)
```

This is the THM-2851 Bockstein value, but it is not holonomy of a closed
path upstairs.  The via path ends at `(a+1,7)` and the direct path at
`(a,7)`.  Equation `(10)` is their vertical `T`-clutch, equivalently the
projective curvature after forgetting `a`.  Appending `T^-1` closes the
lifted loop and contributes `omega^-3`; the honest closed loop has
holonomy one.

This distinction prevents a reduced projective obstruction from being
misreported as nonflatness of the already faithful lift.

## 5. Charged channel, not raw scalar current

The lift acts on the two Prony channels by

```text
(U,V) -> (E U,V).                                      (11)
```

It does not act by the scalar `E` on `S=U+V`.  Exhaustion of the 169
reduced edges gives

```text
S'=E S fails on 144 edges;
the 25 survivors are exactly E=1.                     (12)
```

Moreover, at one fixed frequency section the raw target-translation ratio
is one on all 169 `(q,h)` edges: both endpoint translation exponents
vanish.  The nontrivial factor `(6)` is created by the external assignment
`(a,q)->r(a,q)` plus the event decoration.  It is not a raw translation
phase of one physical current.

## 6. Exact selector XOR calculus

At the zero and stepped endpoint origins, the E3 truth supports are

```text
origin 0:  {0,3,11},
origin 1:  {0,11}.                                    (13)
```

Encode E3 by one and complement by zero.  Let `t_o(q)` be the true block
containing `J_q`, and `s_o(q)` the selected block.  Origin `o` contributes
exactly when `s_o=t_o`; hence the signed coefficient is nonzero exactly
when

```text
s_0 XOR s_1 XOR t_0 XOR t_1 = 1,                      (14)
t_0 XOR t_1 = delta_(q,3).                            (15)
```

The script exhausts all four selector pairs at every q and both fields.
Equations `(14)`--`(15)` give:

- every origin-independent rule `s_0=s_1` has signed support `{3}` on all
  26 samples;
- q11 E3 and q7 complement atoms are individually nonzero and equal at the
  two origins, hence cancel;
- the everywhere-full choice `s_o=t_o` cancels identically; and
- any q7 or q11 support requires `s_0 XOR s_1=1`.

Within the Boolean E3/complement class, there is a unique selector which
keeps the full zero-origin atom and the empty stepped-origin block at every
q.  On the seam it is

```text
q3:  (E3,E3)
q11: (E3,not-E3)
q7:  (not-E3,E3).                                     (16)
```

Its signed coefficient equals the all-q zero-origin atom, so `(16)` removes
q11 cancellation and supplies q7 coefficient support.  It does so by a
target-dependent projection of the stepped origin.  No source/QA/QAB
operation currently types this origin-odd selector.

## 7. Selector parity and carry are independent coordinates

The parity of selector `(16)` is

```text
p(q)=s_0(q) XOR s_1(q)=1 XOR delta_(q,3).              (17)
```

It toggles only at `q2->q3` and `q3->q4`.  At the unique carry edge
`q12->q0`, `p` remains one.  Conversely, the THM-2880 path `q0->q3`
changes `p` while

```text
kappa(0,3)=0.                                          (18)
```

Thus neither coordinate determines the other.  Multiplying by the `u3`
event phase cannot manufacture or type the origin-odd selector.  The
polarization in `(16)` lies precisely in THM-2880's explicitly unclosed
broader signed/independent-origin correspondence class, not its closed
source-induced fibre route.

## 8. Physical stopping gate

THM-2851/2835 provide the exact support incidence

```text
449 sheets: QA(q11,a) -> QAB(q7,a+1).                  (19)
```

But `a` in `(2)` is an externally assigned first-carry coordinate.  No
coefficient-to-ancestry basepoint or sheetwise alignment is constructed.
The marked-gauge hostile still moves the retained q7 endpoint copy's fixed
source into the all-safe/no-word chamber.  Together with the inherited
q11/q7 common-selector zeros and untyped selector `(16)`, this leaves the
physical `QA/QAB` clutch open.

The exact missing sidecar is no longer “another phase.”  It is an
origin-resolved, source-preserving map from the 449 semantic sheets to the
charged coefficient line which simultaneously justifies the origin-odd
q7/q11 selector.

## 9. Reproduction

```bash
python3 04-computation/lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.py
python3 -O 04-computation/lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.py
```

Normal and optimized modes byte-match the stored output.  The script has no
executable Python `assert`.

```text
script 3ed346e0c631b34bd61f0c4d27d7f161e8d35b70decfb95f5207c5f57893d005
output 0faa0a24f6ba8b6c88b6bbfc4f225e38667097b1a937d977741453499884901c
```

The independent audit reconstructed the all-q truth atlas, carry event,
base-13 nonsplit chart, every transport/composition/gauge identity, reduced
curvature versus closed lift, all Boolean selector pairs, channel scope,
fixed-frequency hostile, and selector/carry independence in both fields.
