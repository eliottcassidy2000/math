# A two-arm dyadic family using opposite-direction leaves

A forward leaf at q beats exactly q among the current vertices.
A reverse leaf at p loses exactly to p among the current vertices.
These operations are different and reversing the direction is
essential to this example.

Label the eight-vertex source 0 through 7, with outgoing-neighbor
bitmasks

`(118,124,248,241,224,192,128,3)`.

The source is strong, witnessed by the cyclic path 0→1→2→3→4→5→6→7→0;
both operations preserve strength on a strong source of order at least two.
Let q=7 and p=1. Make a forward leaves at q and b reverse leaves
at p, where a,b≥0. A mixed pair of operations commute up to the
irrelevant order of labels: the reverse leaf beats the forward
leaf regardless of which was attached first. The exact count is

```
H(a,0) = 128·2^a + 1,
H(a,b) = 128·2^(a+b) + 3     (b≥1).
```

Consequently before any reverse move the same two designated
ports realize both unit defects: the forward move has defect
−1 and the reverse move has defect +1, at every height a.
The all-negative route thus keeps a genuine positive alternative
and starts at H129, passing through H257. At a fixed total number
of added vertices n≥1, the actual constructions hit both
128·2^n+1 and 128·2^n+3 with no order penalty. After the
first reverse move, both designated defects become −3 rather
than ±1. This is an all-height fork with a controlled phase
change, not a complete binary tree or a proof of cofinite
coverage.

## Proof

The forward leaves form a transitive module with q its only
outneighbor outside the module. The reverse leaves form a
transitive module with p its only inneighbor outside itself.
A directed Hamiltonian path can therefore visit either module
at most twice: crossing from the first forward visit to a
later forward visit requires exiting to q, and crossing into
a later reverse visit requires entering it from p. Because
q,p are old vertices, each may be used only once.
Compress each maximal visit to a module into a labelled slot,
retaining all old vertices and the order of slots. A skeleton
with j occurrences of a given module has exactly j! S(k,j)
lifts in that module, where k is its size. Indeed partition
the vertices into j ordered nonempty sets, and use the unique
transitive path on each set within the corresponding visit.
The two modules' choices are independent. For j∈{1,2} the
multipliers are 1 and 2^k−2. Thus on the positive grid H(a,b)
is multiaffine in 2^a and 2^b. On either nonempty axis it is
affine in the corresponding power of 2. Subset DP yields the
following 3×3 table (rows a, columns b):

```
129   259   515
257   515  1027
513  1027  2051
```

The four positive-interior corners determine its multiaffine
polynomial, two points on each nonempty axis determine its
affine polynomial, and the origin is checked separately. This
proves the count identities for every a,b.

## Why opposite orientations do not rescue two fixed arms

For any tournament and fixed old p,q, mixed-direction moves
commute physically for the same reason (the reverse leaf
always beats the forward one), independently of this particular
seed. If both defects remained fixed across the other move, the
counts 2(2H−1)+1=4H−1 and 2(2H+1)−1=4H+1 would differ
on isomorphic children. More quantitatively, when the initial
forward defect is −1 and reverse is +1, the new reverse
defect after a forward move exceeds the new forward defect
after a reverse move by exactly 4. The family above realizes
one exact way to absorb that difference: retain +1 down the
forward fan and charge the reverse branch with −3.
