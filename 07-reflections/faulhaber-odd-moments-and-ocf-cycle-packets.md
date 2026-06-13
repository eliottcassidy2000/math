# Faulhaber Odd Moments And OCF Cycle Packets

codex-2026-06-13

The useful new point from the prompt is not only the asymptotic formula.  It is
the sentence "only odd Faulhaber moments matter."  That puts the triangular
tower thread next to OCF in a precise way.

For the power balance, centering at `c=a+n` gives

```text
D_p(c,n)=c^p-2*sum_{r odd<=p} binom(p,r)c^(p-r)S_r(n).
```

The paired terms erase every even moment.  The left endpoint has one unpaired
central term `c^p`, and everything else is an antisymmetric odd moment.  This
is why the p=1 and p=2 anchors are so clean:

```text
p=1: c=n(n+1)
p=2: c=2n(n+1)
```

The higher-p anchors keep the same carrier `u=n(n+1)`, but Bernoulli terms add
a fractional address.

## The Square Pyramid Warning

The p=2 geometry has square pyramids all over it:

```text
6*S_2(n)=n(n+1)(2n+1).
```

That is the cuboid-packing layer.  But the balance equation for the p=2 tower
uses only the odd moment `S_1(n)`.  This is a very good warning for the rest of
the repo: the beautiful visible packing can be true while the proof-driving
imbalance lives one layer below it.

For LRC14, that says not to confuse the pretty shell `C=27`, the product
fibers `d*m`, or the visible endpoint walls with the actual obstruction.  For
code72, it says the formal weight enumerator and the `78` design parameter are
volume shadows until support compatibility is attached.

## Why OCF Is The Right Analogy

OCF does not stop at "count the odd cycles."  It says

```text
H(T)=I(Omega(T),2)=sum_k 2^k alpha_k(T),
```

where `alpha_k` counts compatible packets of vertex-disjoint odd cycles.  The
odd cycles are the atoms; `Omega(T)` is the compatibility object.

The Faulhaber identity gives an atom inventory:

```text
S_1, S_3, S_5, ...
```

What it does not yet give is the compatibility object.  That is the frontier.
The hidden tournament may not have runners as vertices.  It may have odd
moments, boundary events, denominator residues, support packets, or proof
obligations as vertices.

## Merge With HYP-2456

HYP-2456 says the triangular crossover word decomposes into:

```text
Beatty address -> Pell carry -> visible side token.
```

This session adds another layer:

```text
triangular carrier u -> odd Faulhaber moments -> compatibility packets.
```

The first is an endpoint/carry normal form.  The second is a moment normal form.
They are not competitors.  They are two projections of the same triangular
carrier.

## A Productive Next Object

Define a small "odd moment collection function" for a finite shell model.
Possible vertices:

```text
odd moment atoms S_{2k+1}
endpoint wall atoms
Beatty/Pell carry atoms
owner-support atoms in Q27
minimum-word support packets for code72
```

Edges should mean conflict: two atoms cannot be realized by the same hidden
support/lift without spending the same boundary resource.  Then the OCF analogy
would become testable:

```text
single atoms       -> odd Faulhaber moment list
independent sets   -> compatible proof packets
evaluation at 2    -> binary choices / support signs / reversal choices
```

The point is not to force a theorem too early.  The point is to stop throwing
away compatibility as soon as a scalar formula looks pretty.

## What Changed

Before this, the triangular tower thread mostly said:

```text
p=1,2 exact; p>=3 fractional; side overlaps have Pell/Beatty structure.
```

Now it says:

```text
the common carrier is u=n(n+1);
the balance uses only odd moments;
p=1,2 are the vanishing-Bernoulli cases;
HYP-2456 supplies endpoint/carry addresses;
OCF supplies the template for the missing compatibility lift.
```

That is a sharper research object.  It gives LRC14 and code72 a concrete next
question: what are the `alpha_k`-like packet counts after the visible boundary
totals have been fixed?
