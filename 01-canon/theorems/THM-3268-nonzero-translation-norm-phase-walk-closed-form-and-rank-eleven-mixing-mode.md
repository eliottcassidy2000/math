---
id: THM-3268
title: "Nonzero-translation norm-phase walk closed form and rank-eleven mixing mode"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In THM-3267's
  chosen F_169 model, the freely varying nonzero-increment walk on nonzero
  points has norm-phase quotient Q=14*J_12-I_12, with spectrum 167 once and
  -1 eleven times. Its phase-preserving and each fixed phase-shift path
  sequence have explicit two-term closed forms, rational ordinary generating
  functions, and recurrence C_(n+2)=166*C_(n+1)+167*C_n. The same theorem
  holds over F_(q^2)/F_q with quotient (q+1)*J_(q-1)-I_(q-1). Repeating one
  fixed translation does not obey this quotient iteration. This is an
  abstract variable-increment sequence theorem, not a physical ancestry
  walk, current, row exclusion, or LRC(14) decrement.
source: root/creative-synthesis-cont/2026-08-03
audit: >
  The assertion-independent companion hash-pins THM-3255, THM-3267 and its
  exact artifacts; independently rebuilds F_169, the norm phase and every
  translation edge; checks the equitable quotient at all 168 sources; and
  performs a point-level dynamic program through length five without using
  the quotient update. It verifies the closed forms, recurrences, generating
  functions, centered-marker eigenspace, general-q algebraic consequences and
  fixed-translation hostiles. Normal and optimized executions byte-match the
  frozen transcript. An independent implementation rederived the field,
  quotient, formulas, general-q proof and fixed-translation counts; its sole
  wording repair made explicit that an increment is allowed only when the
  next point remains nonzero.
depends_on:
  - THM-3267-norm-phase-factorization-ladder-and-projective-determinant-blindness
related:
  - THM-3255-twelve-balance-multiplicative-singer-rank-defect-and-phase-marker-boundary
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
script: 04-computation/lrc_norm_phase_translation_walk_closed_form_thm3268.py
output: 05-knowledge/results/lrc_norm_phase_translation_walk_closed_form_thm3268.out
script_sha256: fd26af6cba1cad5019556fd2235b319f7d089273b64d7f0dc01d2877ca0a4985
output_sha256: a9d8f74c6f9a101b81fd1c6a652c8f9379d945d2883c502452252fdbfc1fe844
hash_basis: LF-normalized bytes
---

# THM-3268 -- the norm-phase translation census is a two-mode walk

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3267 computes the one-step phase-difference census across every nonzero
translation. This theorem identifies the operation response when the
translation may be chosen afresh at every step. The exponential path family
collapses to two exact modes, while repetition of one fixed translation gives
the canonical hostile.

## 1. Variable-increment walk and equitable quotient

In the chosen model, let

~~~text
X=F_169^*,       phi(alpha^j)=j mod 12.                (1)
~~~

An allowed step has x,t,x+t all nonzero and sends x to x+t. Equivalently,
every ordered pair of distinct points x,y in X is one edge, labelled by the
unique nonzero increment t=y-x. The point graph is the loopless complete
digraph on 168 vertices and has outdegree 167.

For n>=0 and d in C_12 define

~~~text
C_n(d)=#{(x_0,...,x_n) in X^(n+1):
          x_(i+1)!=x_i for 0<=i<n,
          phi(x_n)-phi(x_0)=d}.                        (2)
~~~

Each of the twelve norm-phase fibres has fourteen points. From a fixed point,
there are thirteen targets in the same phase and fourteen in each other
phase. The partition is therefore equitable, with quotient

~~~text
Q=14*J_12-I_12,                                        (3)
~~~

where J_12 is the 12 by 12 all-ones matrix. Its spectrum is

~~~text
167 with multiplicity 1,       -1 with multiplicity 11. (4)
~~~

The first eigenspace is the constant line and the second is the C_12
augmentation space.

## 2. Closed forms, recurrence and generating functions

The two spectral projectors give

~~~text
Q^n=(167^n/12)*J_12
    +(-1)^n*(I_12-J_12/12).                            (5)
~~~

Summing the appropriate entries over the 168 possible starting points gives

~~~text
C_n(0)=14*(167^n+11*(-1)^n),                          (6)
C_n(d)=14*(167^n-(-1)^n),             d!=0.           (7)
~~~

Their sum over d is 168*167^n. For every fixed nonzero d,

~~~text
C_n(0)-C_n(d)=168*(-1)^n.                             (8)
~~~

Both sequences, and more generally every scalar entry or fixed linear
observable of Q^n, satisfy the characteristic recurrence

~~~text
C_(n+2)=166*C_(n+1)+167*C_n.                          (9)
~~~

For the two profiles in (6)--(7), the ordinary generating functions are

~~~text
sum_(n>=0) C_n(0)z^n
  =168*(1-153z)/((1-167z)*(1+z)),                     (10)

sum_(n>=0) C_n(d)z^n
  =2352*z/((1-167z)*(1+z)),          d!=0.            (11)
~~~

Thus length-n multiplicities can be evaluated by one integer exponentiation,
using O(log n) integer multiplications rather than enumerating 168*167^n
paths. At n=1, (6)--(7) recover THM-3267's exact counts 2,184 and 2,352.

## 3. Uniform finite-field theorem

Let q be any prime power. The multiplicative group F_(q^2)^* is cyclic of
order (q-1)(q+1), and the norm is z |-> z^(q+1). Hence its image has q-1
elements and every norm fibre has q+1 elements.

For the analogous loopless complete walk on F_(q^2)^*, quotient by terminal
norm divided by initial norm. The equitable matrix is

~~~text
Q_q=(q+1)*J_(q-1)-I_(q-1),                            (12)
~~~

with eigenvalues q^2-2 once and -1 with multiplicity q-2. If C_(q,n)(eta)
counts paths whose terminal/source norm ratio is eta in F_q^*, then

~~~text
C_(q,n)(1)
 =(q+1)*((q^2-2)^n+(q-2)*(-1)^n),                    (13)

C_(q,n)(eta)
 =(q+1)*((q^2-2)^n-(-1)^n),          eta!=1.         (14)
~~~

For q=2, (14) is vacuous. Each sequence in (13)--(14) obeys

~~~text
C_(n+2)=(q^2-3)*C_(n+1)+(q^2-2)*C_n.                 (15)
~~~

This is a fibre-size proof for all prime powers, not a finite interpolation.

## 4. The eleven-dimensional marker mode

THM-3255's centered phase markers have the form

~~~text
M_r=12*e_r-1.                                          (16)
~~~

They span the eleven-dimensional augmentation space, and (3) gives

~~~text
Q*M_r=-M_r.                                            (17)
~~~

After normalizing Q by its outdegree, the whole marker sector has multiplier
-1/167 per freely varying step. This explains the alternating discrepancy in
(8). It does not realize the marker physically: the quotient has forgotten
the point, increment word, intermediate geometry and ancestry.

## 5. Repeating one translation is the sharp hostile

The fixed translation t=(1,0) has the same one-step relative histogram as a
row of Q,

~~~text
n=1: (13,14,14,...,14).                               (18)
~~~

That coincidence does not authorize iterating Q. Reuse the same t and reject
a path as soon as it hits zero. Exact enumeration gives

~~~text
n=2:  (12,14,14,...,14),      166 valid paths,
n=13: (156,0,0,...,0),        156 valid paths.         (19)
~~~

At length thirteen, characteristic 13 returns every surviving point to
itself. The twelve nonzero points on the line F_13*t have encountered zero;
the other 156 survive. Thus (6)--(15) count freely varying admissible
increments, not powers of one affine translation, a constrained H_13 word,
or a path on a literal ancestry sheet.

## 6. Verification and scope

The companion performs a direct point-level dynamic program through n=5,
retaining all 168 vertices and the initial phase without using Q in its
update. It matches every phase profile from (6)--(7), verifies the recurrence
and generating functions, and freezes both fixed-translation hostiles.

Run

~~~text
python3 04-computation/lrc_norm_phase_translation_walk_closed_form_thm3268.py
python3 -O 04-computation/lrc_norm_phase_translation_walk_closed_form_thm3268.py
~~~

and compare LF-normalized bytes with the declared transcript.

The theorem is an exact abstract sequence law. It supplies no deterministic
physical transition, same-ancestry allocation, endpoint current, cellwise
safety statement, row exclusion, or decrement for LRC(14).

QED.
