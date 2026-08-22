---
id: THM-3262
title: "Joint two-face selector bit, Pareto threshold policy, and signed collision tax"
status: >
  PROVED + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED on the promoted
  support-(1,2), bank-I2 and support-(1,3), bank-I2 faces. After padding the
  small face into eight multiplicity coordinates, exactly 19 identical
  exclusive state vectors demand opposite choices from rows 2 and 10; hence
  no face-blind function of the count vector can be lawful on both faces.
  One persistent face-origin bit is necessary and sufficient. In the declared
  axis-threshold-tree class augmented by that bit, minimum depth is 6 with 23
  leaves, while the absolute minimum 22 leaves is first attained at depth 8.
  Forgetting the face identifies 239 raw state pairs and 80 signed-Parikh
  pairs. The two distance sequences nevertheless have one compact catalytic
  closed form u*D_12(z)+v*D_13(z). This is a two-face selector and
  enumeration theorem, not FC(3), SFC(3), or a global moment functional.
source: root/creative-synthesis-cont/2026-08-03
audit: >
  The repaired exact companion hash-pins and replays THM-3244 as well as its
  six original dependencies; reconstructs the complete 239- and 4,319-state
  banks and all response directions; exhausts the threshold-tree recurrence;
  verifies the recovered trees on every nonreset state; computes all raw,
  signed and unsigned collision quotients; and independently multiplies the
  coordinate distance factors. Normal and optimized runs byte-match the
  frozen transcript, and the source has no assertion node or floating
  literal. An independent reconstruction reproduced both universes, every
  conflicting fibre, the policy-independent Parikh counts and distance
  polynomials, and all policy-specific continuation counts. The audit also
  confines the DP optima to the declared tree model, the face bit to a static
  persistent origin tag, and the ordered-word counts to the displayed
  controller.
depends_on:
  - THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge
  - THM-3256-optimal-two-row-threshold-policy-injective-signed-trace-and-factored-distance-enumerator
related:
  - THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go
  - THM-2000-support-harmonic-abel-dini-figurate-surface
  - THM-2005-support-dirichlet-automatic-tournament-atlas
script: 04-computation/fc3_joint_two_face_policy_sidecar_beach_20260803.py
output: 05-knowledge/results/fc3_joint_two_face_policy_sidecar_beach_20260803.out
script_sha256: 1d69c98d226648d1289c36c72250126e9dbf561c73aa6c3e846ddf3e1a6b47e8
output_sha256: 67aa3d2e7e6bc836d5a1d2e79caf16f9d1d9cfb5059cb1666fe2e8b74da792f0
hash_basis: LF-normalized bytes
---

# THM-3262 -- the face marker is one bit, not decoration

**PROVED + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED on two promoted
physical faces.**

THM-3256 compiles a small row-2/10 selector on the support-(1,3), bank-I2
face, while showing that its signed operation trace retains every physical
start. The same two rows cover the support-(1,2), bank-I2 face. This theorem
asks whether the two selectors can be merged after forgetting which face
supplied the count vector.

They cannot. The missing face coordinate is logically load-bearing. Once it
is restored as one bit, the joint policy has an exact depth/leaf frontier,
and the same marker becomes a catalytic variable in the joint distance
generating function.

## 1. Two exact state boxes

Embed both faces into n=(n1,...,n8). The small and full boxes and their unique
resets are

~~~text
small capacities=(4,3,2,1,1,0,0,0),
q_12=(1,2,1,1,1,0,0,0);

full capacities =(4,3,2,2,2,1,1,1),
q_13=(1,0,2,1,1,1,1,1).                               (1)
~~~

Deleting the empty physical multiset leaves

~~~text
small states=239,       full states=4319,
tagged joint states=4558.                               (2)
~~~

Rows 2 and 10 cover every nonreset state on both faces. Their exclusive
response records have the census

~~~text
small: row 2 only=52,  row 10 only=19,  overlap=167;
full:  row 2 only=304, row 10 only=103, overlap=3911.   (3)
~~~

Only the 71+407=478 exclusive records constrain a deterministic selector.

## 2. Nineteen contradictory fibres

Forget the face tag but retain the padded count vector. Exactly 19 vectors
occur as exclusive records on both faces with opposite required labels:

~~~text
small face requires row 2,
full face requires row 10.                              (4)
~~~

The smallest witness is the literal one-pole state

~~~text
(5),       n=(0,0,0,0,1,0,0,0).                       (5)
~~~

Therefore no function of the untagged count vector, regardless of its circuit
or lookup representation, can choose a lawful member of {2,10} on both
faces. This conclusion is stronger than a threshold-tree lower bound.

For a static origin sidecar, zero bits are impossible by (5), while one bit
encoding the two faces is sufficient. A previous-chart bit initialized
identically on both faces also fails at the first decision (5). This does not
rule out later history-dependent memory, and a face-initialized bit is
equivalent to the face tag only while it is retained.

## 3. Exact threshold-tree Pareto frontier

The declared tree model has binary internal tests

~~~text
n_j<=h                                                     (6)
~~~

for every nontrivial coordinate threshold, together with the optional Boolean
face test. There are 11 small-face thresholds, 16 full/joint thresholds, and
17 joint tests after adding the face bit.

For a classified subset S, let L_d(S) be the least leaves in a tree of depth
at most d. Pure sets cost one, mixed sets at depth zero are impossible, and

~~~text
L_d(S)=min_T [
 L_(d-1)(S intersect T)+L_(d-1)(S minus T)].            (7)
~~~

Exact bit-set dynamic programming gives:

| policy class | minimum depth | leaves at minimum depth | absolute minimum leaves | first depth attaining it |
|---|---:|---:|---:|---:|
| small face, axes | 4 | 8 | 8 | 4 |
| full face, axes | 5 | 15 | 15 | 5 |
| joint, axes only | impossible | - | impossible | - |
| joint, axes plus face bit | 6 | 23 | 22 | 8 |

A gate-first composition of the two separate trees has depth 6 and 23 leaves,
so it attains the joint minimum depth. One recovered minimum-depth tree first
tests n8=0 and asks for the face only on that branch. It consults the bit on
all 238 small nonreset states and 2,159 of the 4,318 full nonreset states.

The leaf optimum has a genuine tradeoff: a depth-8 tree with two conditional
face-test nodes has 22 leaves and consults the bit on only 198 small and 599
full nonreset starts. These are optima only among binary decision trees built
from (6) and the face test; no DAG, arbitrary arithmetic-circuit, or stateful
policy lower bound is asserted.

## 4. Three different quotient taxes

Every small raw count vector also occurs in the full box. Thus forgetting the
face in the disjoint tagged union gives

~~~text
4558 tagged states -> 4319 untagged vectors,
raw face-forgetting tax=239.                            (8)
~~~

For a Q-monotone route, the signed Parikh vector is q_f-n. It reconstructs n
when f is retained. Equality across faces is equivalent to

~~~text
m=n+(q_13-q_12)=n+(0,-2,1,0,0,1,1,1).                 (9)
~~~

The admissible cross-face choices number

~~~text
5*2*2*2*2=80.                                          (10)
~~~

Hence the policy-independent signed quotient is

~~~text
4558 face-tagged signed classes,
4478 face-forgotten signed classes,
signed-Parikh collision tax=80.                        (11)
~~~

Unsigned Parikh vectors are coarser. The small and full faces have 96 and
1,536 values, the entire small set lies in the full set, and therefore

~~~text
unsigned untagged classes=1536,
unsigned face-tagged classes=1632.                     (12)
~~~

For one reproducible controller--the minimum-depth tree followed by the
least-coordinate lawful direction--ordered words retain more data than their
Parikh vectors. Its continuation counts are:

| retained output | face forgotten | face tagged |
|---|---:|---:|
| chart | 313 | 345 |
| unsigned coordinate | 4125 | 4208 |
| chart plus coordinate | 4520 | 4548 |
| signed edit | 4536 | 4558 |

The final gap consists of 22 cross-face pairs with identical complete signed
words. These 22 collisions and the switch counts are controller-dependent;
the 80 in (10)--(11) is intrinsic to the two tagged state boxes.

## 5. Closed forms and the catalytic face variable

Coordinatewise absolute-deviation enumeration gives

~~~text
D_12(z)
 =(1+2z+z^2+z^3)(1+z)^4(1+2z)-z^6,                   (13)

D_13(z)
 =(1+2z+z^2+z^3)(1+z)^4
   (1+z^2)(1+z+z^2)(1+2z)^2-z^8.                     (14)
~~~

The subtracted monomials remove the forbidden empty states. Their sum is

~~~text
D_joint(z)
 =(1+2z+z^2+z^3)(1+z)^4(1+2z)
   [1+(1+z^2)(1+z+z^2)(1+2z)]-z^6-z^8               (15)

 =2+19z+82z^2+220z^3+426z^4+648z^5+803z^6
   +821z^7+687z^8+467z^9+250z^10+101z^11
   +28z^12+4z^13.
~~~

The coefficients sum to 4,558 and independently match the complete route
distance histogram. If a downstream consumer needs the face, the correct
closed form is the catalytic refinement

~~~text
u*D_12(z)+v*D_13(z).                                   (16)
~~~

Specializing u=v=1 is lawful for distance counting and unlawful for response
selection. The face coordinate is therefore a control sidecar before
abelianization and a catalytic sequence variable after it.

## 6. Typed scope and next boundary

The exact information-flow contract is

~~~text
source:      tagged state (face,n);
quotient:    forget face, then optionally abelianize the route;
preserved:   pole multiset, distance after u=v=1;
destroyed:   which response row is available, reset origin, cross-face
             multiplicity of signed traces;
sidecar:     one persistent face-origin bit, or u/v after enumeration;
hostile:     state (5) for policy, 80 doubled signed Parikh values for counts.
                                                                    (17)
~~~

This theorem concerns only the two promoted bank-I2 faces and lawful rows 2
and 10. It supplies no selector for the other maintained faces, no globally
optimal previous-chart memory controller, no proof that adaptive row choice
is a single scalar Gaussian-moment functional, and no conclusion about
FC(3), SFC(3), or the Jacobian conjecture.

## 7. Exact reproduction

Run

~~~text
python3 04-computation/fc3_joint_two_face_policy_sidecar_beach_20260803.py
python3 -O 04-computation/fc3_joint_two_face_policy_sidecar_beach_20260803.py
~~~

and compare LF-normalized bytes with the declared transcript. The companion
pins all eight directly consumed dependency artifacts, replays the inherited
banks, verifies both recovered trees on every state, and checks every count,
hash and polynomial using exact integer arithmetic.

QED.
