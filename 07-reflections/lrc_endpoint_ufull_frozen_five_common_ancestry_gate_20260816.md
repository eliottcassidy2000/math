# The frozen U_full role residues pass the scalar Boolean gate but not the ancestry gate

**Status: FINITE-EXACT typing and realization sidecar.**  The incoming
`13^3` endpoint bank and its graph nonvanishing are inherited and are now
incorporated in promoted THM-3479.  This note supplies the strict realization
boundary: it does not construct a physical current, decide U_clock, or remove
an LRC(14) row.

## 1. Inheritance and the exact next question

The completed U_full guard-refinement package gives five distinct normalized
inverse-DFT images in the certified split-prime chart

```text
p=572252886246508880869:

(0,0,0)  -> 405336876493642499425,
(0,1,0)  -> 518539850465495448196,
(1,0,0)  -> 503604956476841920373,
(1,0,1)  -> 320618948602619577408,
(1,12,0) ->  15703541686881447885.                  (1)
```

The guard bridge is

```text
A(1,0,1)-A(1,0,0)
 =389266878372286537904 mod p !=0,                  (2)
```

and the bridge, left `K4`, right `K4`, and their product are nonzero in
all `72` declared labelled role-contract charts.  Those are exact endpoint-aggregate graph
statements.

The closest proved realization mechanism is THM-2471, whose weights are
marginals of one finite Boolean ancestry stalk.  THM-2538 identifies the
positive common-base source/arrival product that would be needed after a
genuine later field exists.  The corrected near miss is MISTAKE-293:
sharing an ancestry fibre does not mean evaluating every factor at one circle
point.  The canonical hostile is MISTAKE-261 and its descendants: separate
endpoint or graph marginals do not supply the joint target/current.

The precise question is therefore not whether the five numbers in (1) make a
nonzero graph polynomial—they do—but whether they carry enough information to
produce an atom-level U_full coupling before endpoint marginalization.

## 2. A positive result: there is no scalar Boolean obstruction

Let the common base consist of five atoms, one for each class in (1), with
counting measure and nonnegative weight equal to the displayed representative.
For each role label, take the indicator of its class atom.  Then

```text
integral F 1_(class(label)) = displayed role value.             (3)
```

The repeated roles `c2,q3,q4,q5` use the same event.  Thus (3) is a literal
nonnegative weighted Boolean realization of all eight vertex potentials on
one finite base.  Computing the signed edge differences and integer matrix-tree
factors—not reducing first—gives no zero among the `72` products.  Their
absolute values lie between

```text
21899028683147318541426405228605418004021414198897289055920845818968807139208402909708300365999230796601207111498655940614096947873983762889480
```

and

```text
379988029205165189064531036950840125177333602108490073274036822188321102898888717366956823048010693841669688048415379170173005548411429146037650.
```

Reducing those integer chart rows modulo `p` reproduces the incoming graph
digest

```text
b7d8c2c9860e4f1aa542b1c85fdb7b65cf4985aba5a81a84ff3a324834d51c51. (4)
```

So the cheapest scalar gate passes strongly.  The remaining obstruction is
not positivity of some abstract vertex masses, graph incidence, or the tree
polynomial.  It is provenance and coupling.

This construction is deliberately called **synthetic**.  Its atom labels are
the five residue classes themselves.  They have no U_full owner, root, deep,
clock, or endpoint meaning.

## 3. Why reduction cannot certify Boolean positivity

The map from the characteristic-zero endpoint ring to the split-prime image
does not preserve order.  For example, the positive representative of the
guard class and the negative lift

```text
320618948602619577408-p
 =-251633937643889303461                            (5)
```

have the same image modulo `p`.  The first participates in the synthetic
nonnegative realization above; the second cannot be the mass of an indicator
against a nonnegative weight.  Therefore

```text
nonnegative Boolean event mass
```

is not a predicate on the frozen finite-field tuple.  Equation (5) does not
say that the actual cyclotomic U_full lift is negative.  It says the stored
residues alone cannot decide the sign/order gate.

## 4. Why separated endpoint marginals cannot certify ancestry

The endpoint worker makes the loss point explicit.  For every character it
computes

```text
ax = fast_x_sweep(...),
by = fast_endpoint_sum(...),
gamma = phase * ax * by,                            (6)
```

after both endpoint sums have already been taken.  Its output is only

```text
(alpha, interval_count, tuple(gamma scalars)).       (7)
```

No left/right atom identifier, stalk coordinate, or ancestry relation
survives in (7).

The minimal exact hostile uses only two atoms of mass `1/2`.  Put

```text
left          =(1,0),
right_aligned =(1,0),
right_crossed =(0,1).                               (8)
```

Both right fields and the left field have mass `1/2`, so both products of
separately computed marginals equal `1/4`.  But the products formed before
integration are

```text
integral left*right_aligned =1/2,
integral left*right_crossed =0.                     (9)
```

Thus even exact characteristic-zero marginals and their product do not
determine the common-base intersection.  The kernel is the familiar
two-by-two checkerboard coupling direction.  Finite Fourier inversion of the
post-marginalization scalars is linear and cannot restore that lost coupling.

More precisely, if `U` and `V` are the left- and right-atom spaces and
`epsilon_U,epsilon_V` are total mass, the marginal map sits in the exact
sequence

```text
0 -> ker(epsilon_U) tensor ker(epsilon_V)
  -> U tensor V
  -> U x_k V -> 0.                                  (9a)
```

THM-2538 gives the categorical dimension `(m-1)(n-1)` and mixed-Haar basis;
THM-3492 independently realizes the same sequence for Rule 30 slack and phase,
where its kernel is `(q+1)P_(N-1) tensor ker(e_h)`.  The hostile `(8)--(9)` is
the smallest nonzero vector in `(9a)`.  Hence the missing ancestry key is
exactly data selecting a semantically lawful lift through this sequence, not
one more scalar appended to the two marginals.

One can regard (6) as living on an independent product base.  That does not
make it a lawful LRC ancestry: the product base has no support condition
linking the two endpoints through owner-normalization sheets, collision roots,
or chronological nodes.  THM-2471 and THM-2594 retain precisely such linked
coordinates before summing.

## 5. The first unavailable coordinate and API

The first unavailable coordinate is

```text
omega = a shared U_full ancestry key
```

linking one left `fast_x_sweep` atom to one right `fast_endpoint_sum` atom
before either marginal is taken.  Equivalently, the missing API must produce

```text
J_ell(omega_L,omega_R)                              (10)
```

with all of the following properties:

1. `J_ell` is supported on a typed U_full ancestry relation, not the full
   independent Cartesian product;
2. both endpoint factors and the character phase are inserted before the
   atom sum;
3. the atom-level inverse DFT at
   `q_H=(1,0,1)` and `q_q5=(1,0,0)` is defined on that same carrier;
4. its bridge difference reduces to (2);
5. owner/deep/root and chronology labels survive long enough to compare with
   THM-2471/2538.

THM-2594 demonstrates that a linked-node Boolean expansion is possible on the
different canonical row.  It explicitly proves one auxiliary theta-slaved
contraction, not the THM-2334 endpoint bank and not a generic transplant.
Consequently it supplies a model for (10), not the missing U_full map.

The cheapest next computation is only the bridge version of (10).  There is
no reason to build both `K4` factors until the atom-level difference of the two
classes in (2) has been recovered.  A failure there is decisive for this
realization route; a success would inherit the already completed graph audit.

## 6. Connection contract and scope

| field | exact content |
|---|---|
| source | five frozen U_full endpoint residues and the `ax*by` worker |
| attempted target | nonnegative Boolean/common-ancestry role current |
| scalar map | five class atoms with weights from (1) |
| scalar result | PASS; all integer graph products are nonzero |
| first loss | split-prime reduction forgets characteristic-zero sign/order |
| second loss | separate endpoint sums forget left/right atom pairing |
| hostile | (8)--(9), same marginals and separate product, different joint mass |
| needed sidecar | typed atom coupling (10) on one U_full ancestry fibre |
| cheapest test | recover only the `H-q5` bridge (2) before either `K4` |

Nothing here obstructs the existence of a lawful U_full coupling.  Nothing
here proves a grouped exact-address coefficient `C(a;X,m)`, an all-unit
projector `B(q)`, a physical current, a scalar-row exclusion, a U_clock
statement, or LRC(14).

## 7. Exact companion

Run from the repository root:

```text
python -B 04-computation/lrc_endpoint_ufull_frozen_five_common_ancestry_gate_20260816.py
python -B -O 04-computation/lrc_endpoint_ufull_frozen_five_common_ancestry_gate_20260816.py
```

The normal and optimized transcripts agree with the stored output.  LF hashes
are

```text
script:   b897eef476f6dd7f29793b2cc795d6840b66dc2e0e23f6f051e4a4d6b4b704d0
output:   d4043a74c600f5737fb32aab6b7b7536adb77e7af6238145fb76a2cf43ccd652
semantic: 674f52e86cabcc72d97a193528434cc9524ae3111439e53e49a43a811d4bf41a
```
