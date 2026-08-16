# The Fibonacci 24-state packets are one transitive-T4 atlas with two connections

**Status: VERIFIED-EXACT structural sidecar to PROVED THM-3487 +
INDEPENDENTLY AUDITED.**  The finite claims below are reproduced by
`04-computation/fibonacci_t4_atlas_h1_branch_transplant_probe_20260816.py`.
They sharpen the branch-transplant question without supplying a full-tree,
physical-current, LRC, or Jacobian map.

## 1. Inheritance pass

The closest proved general mechanism is THM-2622, not a new orbit-count
argument.  For a cyclic affine torsor, parallel sections are exactly fixed
points of affine holonomy.  A nonzero pure translation of `V4` has no fixed
point; the identity has four.

THM-3339 supplies the specific Fibonacci data:

```text
six orders pi_i of {p,q,r},
closed owners o_i=(0,p,p,r,q,q),
unique linear frame steps L_i:pi_i->pi_(i+1).            (1)
```

The four-state probe supplies the projective Fibonacci cycle

```text
G=(0 1 2 3).                                             (2)
```

THM-3364 already has a different positive `T4`: one Berggren parent and its
three children form a local transitive tournament under hypotenuse order.
The object below is not that chamber.  It is the atlas of all labelled
transitive `T4`s carried by owner-plus-channel-order states along the
Fibonacci path.

## 2. The exact static atlas

For `V4=F_2^2`, an owner `u` and an ordering `(d_1,d_2,d_3)` of the three
nonzero directions define

```text
Theta(u;d_1,d_2,d_3)
  =(u,u+d_1,u+d_2,u+d_3).                                (3)
```

The four entries in (3) are all vertices of `V4`.  Conversely, any total
order `(v_0,v_1,v_2,v_3)` is recovered uniquely by

```text
u=v_0,
(d_1,d_2,d_3)=(v_1-v_0,v_2-v_0,v_3-v_0).                (4)
```

Therefore

```text
V4 x Ord({p,q,r})  <->  all 4!=24 total orders of V4.    (5)
```

Orienting every pair from earlier to later turns (5) into all 24 labelled
transitive tournaments on four vertices.  This is the rigorous sense in
which the six-state order packet and four-state owner packet are “essentially
a tournament of size four.”  Neither factor alone is a `T4`; their product is
the complete transitive-`T4` atlas.

The affine group acts exactly as it should.  If

```text
A_i(x)=L_i x+t_i,                                        (6)
```

then `L_i` carries the ordered directions `pi_i` to `pi_(i+1)`, and hence

```text
A_i Theta(u;pi_i)=Theta(A_i u;pi_(i+1))                 (7)
```

position by position.  Every lawful moving-owner edge is thus an honest
vertex relabelling of the ranked `T4`.  The connection varies with `i`; it is
not one static relabelling.

## 3. The base-preserving conjugacy equation

Let `X` be the frame-line bundle

```text
T_X(i,j)=(i+1,Gj),                                      (8)
```

and let `Y_h` be an affine owner bundle whose edge maps `B_i` have linear
parts `L_i` and total holonomy translation by `h in V4`:

```text
H_B=B_5...B_0=tau_h.                                    (9)
```

A base-preserving fibre map has the form

```text
F(i,j)=(i,f_i(j)),                                      (10)
```

where each `f_i:P1(F_3)->V4` is a point bijection.  Equivariance is exactly

```text
f_(i+1) G=B_i f_i.                                      (11)
```

Once `f_0` is chosen, (11) determines every later `f_i`.  Closing after six
steps gives

```text
f_0 G^6=H_B f_0,
tau_h=f_0 G^2 f_0^(-1).                                 (12)
```

Now `G^2=(0 2)(1 3)` is a double transposition.  On the affine four-point
set, the three double transpositions are exactly the nonzero translations
`tau_p,tau_q,tau_r`.  Equation (12) therefore has:

```text
h=0:       0 point gauges;
h=p,q,r:   8 point gauges for each h.                    (13)
```

The count eight is the centralizer size of a double transposition in `S4`.
The exact probe exhausts all 24 initial bijections and checks the partition

```text
24=8+8+8.                                               (14)
```

For every accepted gauge, recursion (11) closes, the resulting map (10) is a
base-preserving conjugacy, and (7) transports every ranked tournament
positionwise.  This is a genuine `T4`-atlas transplant, not merely an
abstract permutation conjugacy.

## 4. Obstruction and repair are the same datum

The lawful THM-3339 owner connection has `h=0`.  It has four parallel owner
sections by THM-2622, but (12) gives no frame-line transplant.  Every nonzero
repair has

```text
cycles(Y_h)=12^2=cycles(X),                             (15)
```

and exactly eight structured base-preserving transplants.  Yet `tau_h` has no
fixed point, so THM-2622 forbids every closed owner section.

Thus the exact dichotomy is

```text
zero H1 class:
  closed owner sections, no projective/T4 transplant;

nonzero H1 class:
  projective/T4 transplant, no closed owner section.    (16)
```

Provisional THM-3487 counted `288` abstract conjugacies for each repaired
`12^2` permutation.  Equation (13) identifies the structured subset:

```text
8 base-preserving, positionwise-affine T4 conjugacies
inside 288 abstract conjugacies, for each h!=0.          (17)
```

The other 280 forget some combination of base order, affine frame, and
positionwise tournament transport.

An independent audit and the strengthened probe also exhaust all
`4^6=4096` translation cochains.  Each of the four seam classes has exactly
`1024` representatives.  Class zero always has four sections and cycle type
`6^4`; every nonzero class has no section and type `12^2`.

## 5. Why no single K4 relabelling explains the motion

Element orders in `S4` are only

```text
1,2,3,4.                                                (18)
```

The lawful full shift has order six and each repaired shift has order twelve.
Therefore no fixed permutation of the four vertices induces either full
24-state motion.  Equation (7) is instead a time-dependent affine connection:
each edge is a lawful relabelling, while their ordered holonomy is the global
invariant.

This resolves a recurrent false dichotomy.  The packet is not “a tournament”
versus “not a tournament.”  It is one static atlas of tournaments equipped
with inequivalent connections.

## 6. H1 and the branch-tree boundary

For the fixed edge-linear system, after an affine gauge trivializes its
identity return, the translation seam is the cellular class

```text
[c] in H^1(C_6;V4)=V4.                                  (19)
```

Equation (12) says the projective bundle demands `[c]!=0`; THM-2622 says a
parallel affine-owner section demands `[c]=0`.  This is an explicit finite
model of a word/flux compatibility tariff.

The restriction is load-bearing.  If the return linear part is the coordinate
swap rather than the identity, THM-2622 gives two fixed sections and the
24-state skew product has cycle type `(6,6,12)`.  It lies outside the
pure-translation classifier (19); arbitrary affine local systems are governed
by the full equation `(I-A)x=c`, not just a class in constant-coefficient
`H^1(C_6;V4)`.

The full Berggren ancestry graph is a tree before quotients or relations are
imposed.  Ordinary cellular `H^1` of a tree vanishes, so edge gauges can be
chosen recursively there: local set-theoretic transplant is not the hard
part.  The obstruction reappears when one asks that the recursion descend to
a periodic quotient, respect a branch-monoid representation, preserve a
signed current, or close around relations.  Those operations create the
loops or equivariance equations on which holonomy lives.

This suggests a sharper full-branch program:

1. construct the affine/T4 identification recursively on the free ancestry
   tree;
2. compute its defect under every demanded quotient, branch generator, and
   current;
3. package those defects as explicit cocycle classes in the correct local
   systems; and
4. ask whether one sidecar kills all defects simultaneously.

That is closer in grammar to the desired LRC word-current versus JC-flux
`H^1` map than comparing two 24-element sets.

## 7. Harmonic and subset boundary

Once a uniform 24-address origin is declared, any selected collection of the
transitive tournaments is a periodic Boolean word.  Its ordinary harmonic
coefficient sees only its cardinality.  THM-3364, and the polynomial
extension proposed in THM-3490, say that the full character bank recovers the
address subset and cyclic phase.  Neither transform recovers the affine
connection (9) unless the temporal edge maps are retained as a sidecar.

The distinction is exact:

```text
state subset / H0 data:       which T4 atlas entries are marked;
connection / H1 data:         how successive atlas entries are transported.
                                                               (20)
```

A Fourier address word cannot substitute for holonomy, and holonomy cannot
substitute for the marked subset.

## 8. Reproduction and scope

Run

```bash
python -B 04-computation/fibonacci_t4_atlas_h1_branch_transplant_probe_20260816.py
python -B -O 04-computation/fibonacci_t4_atlas_h1_branch_transplant_probe_20260816.py
```

Script/output/semantic LF SHA-256 are

```text
0f2ea46b6df1f58be7299ed6be64aef34911559421422b9c43e11b8a275988e5
76f7dee37b1a527e46f3048ca0600e49f96c0d2a80066944dc4cec46491abe01
3e19be5a6d7679656f4e476ea1cb5a3bd5944a190bab860c43a8ecbd034d8176
```

The probe checks all 24 total orders, every edge/state transport, all four
seam classes, all 24 initial point gauges, all `4096` translation cochains,
all structured conjugacy equations, cycle types, the nontrivial-linear-return
hostile, and the static-`S4` order hostile.  It proves no full-tree
equivariance, no signed-current preservation, no physical owner, no LRC
bispectrum statement, no Jacobian flux, and no case of LRC(14).
