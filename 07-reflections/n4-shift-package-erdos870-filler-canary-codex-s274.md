# n=4 Shift Packages, Erdős-870 Fillers, and Canary Fibers

Anchors: HYP-3146, HYP-3145, HYP-3143, HYP-3142, HYP-3141, HYP-3140,
HYP-3134, HYP-3133, HYP-3124, HYP-3054, HYP-3053, HYP-3049, T1211, LTI-272,
LTT-170, OPEN-Q-108.

The useful realization is that the user's two `n=4` tournament tables are
two different proof moves.  After rebasing over S276/S274, HYP-3143 names the
exact-order subbasis audit, HYP-3145 names the filler-core interface, and this
note, HYP-3146, names the canary/scaffold policy for the remaining fiber mass.

The fixed-Hamiltonian-path table is a cover.  With path `0->1->2->3`, the
free chords

```text
a=(0,2), b=(1,3), c=(0,3)
```

give `T,+,-,S` for `E,a,b,c`, and every mixed pair is `S`.  But the full
fixed-path cube has eight states:

```text
T: E
+: a
-: b
S: c, ab, ac, bc, abc.
```

So `S` has fiber size `5`, with PGF `z+3z^2+z^3`.  It is not a group quotient:
it is a representation-rich cover with a generator-local absorbing behavior.
The triple representative `abc` is the only delete-one-stable `S` state in the
small scout.  That is canary behavior.  The redundancy is exactly what lets a
state survive deletion.

The two-bit table is a finite-filler scaffold.  Fixing

```text
0->2, 0->3, 2->1, 1->3
```

leaves endpoint variables `x=(0,1)` and `y=(2,3)`, with partial outscore
vector `(2,1,1,0)`.  Then `E,x,y,xy` are exactly `T,+,-,S`, and the class table
is the Klein four table.  This model has no hidden `S` fiber: it is quotient
legal, but it has lost the deletion-stable representation cluster.

The hidden bridge is the Boolean compression

```text
x = a OR c,
y = b OR c.
```

The `c` chord is the clustered canary: if it appears, both endpoint variables
turn on.  This is the finite version of the proof pattern in the formalized
Erdős #870 negative answer: sparse core plus finite fillers plus shift
packages plus clustered canaries.  That theorem is additive-basis work, not a
tournament theorem, but the proof architecture is directly useful here.

The improved thought process for LRC packets is:

```text
When a quotient has a dangerous fiber, do not ask "can I ignore it?"
Ask which next operation is coming.

If the next operation is deletion:
  keep the redundant fiber as a canary cluster.

If the next operation is classification or gluing:
  add finite scaffold/filler data until the quotient is a shift package.
```

This is the missing middle between HYP-3053 fixed-path half-tilings and
HYP-3054 observer-extension payloads.  It also sharpens HYP-3141 and HYP-3142:
an edge witness or k=8 bounded-core exit should state whether its local class
fiber is being used as canary redundancy, or whether a finite scaffold has
made the quotient congruent enough to carry a proof.

Concrete packet fields:

```text
shift_package_scaffold_id
fixed_path_cover_fiber_pgf
canary_cluster_fiber
delete_one_stable_representative
monotone_or_compression_word
finite_filler_arc_set
quotient_congruence_status
deletion_stability_status
terminal_exit_or_named_debt
```

Assumption challenge: the live vertices were not runners, arcs, or raw
tournament classes.  They were proof carriers: finite-filler scaffold,
clustered canary fiber, monotone OR compression, fixed-path cover, edge
tip/tail witness, fiber-PGF packet, and raw class-count shadows.  The scout's
tournament has directed 3-cycles among the scaffold/canary/compression region,
which is the right warning: there is no single best quotient.  The next
operation decides which carrier is legal.
