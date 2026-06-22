# LRC boundary-state induction from tournament ears

The tournament induction archive says the same thing over and over: deletion
only works after the right boundary state has been chosen.  Strong components
give multiplicative atoms, but the inverse step is not a raw new vertex.  It is
a labelled ear, and its exact contribution is

```text
H(T+x)=start(sig=1)+end(sig=0)+Q(sig=0,sig=1).
```

The boundary state `(start,end,Q)` is the object.  The scalar `H` is the
readout after the state has been propagated.

This is exactly the discipline the LRC proof needs.  Removing a large speed is
not a magic size induction.  It is a finite-comb boundary estimate:

```text
meas(A cap Safe_v) >= (1-2/n)meas(A) - 2 components(A)/v.
```

The boundary state is `(measure, components)`, or for many large speeds,
`(core floor, arc count, resonant-pair graph)`.  If we forget that state, the
induction becomes false under dilation; if we keep it, the unbounded branch
descends to smaller LRC cases.

So the LRC14 proof tree now has a clean shape:

1. missing small prime gives the direct `t=1/p` witness;
2. scale-separated speeds peel by finite-comb induction;
3. up to six large speeds close by the union bound;
4. seven or more large speeds need a second-moment divisibility graph;
5. the bounded covering core is the irreducible finite base.

The bounded core is where tournament analogies should stop being deletion
analogies and become boundary-state analogies.  Use the half-tiling address,
AP/Goddyn-Wong tight locus, Fejer/additive-energy labels, and the
missing-depth parity Newton packets.  Runner deletion is the wrong minor order;
residue/cut/depth boundary ledgers are the right one.
