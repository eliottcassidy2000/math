# Arbitrarily wide neutral returns and the next source boundary

**PROVED:** variable-width neutral rewrites, complete prefix classification
in the declared template model, and a mirrored extra +2 exit. **FINITE-EXACT:**
the complete one-apex exclusion at1817 and the stated ready-parent queries.
**OPEN:** unlimited repetition. Definitions and inherited credits are in
[PROOF.md](PROOF.md). The root and Berggren lanes independently found the
generalization below; there is no literature-priority claim.

## 1. The complete varying-width family

Let r>=4,t>=1. Let z start O,end I and satisfy

```
H(z)=c_r=3*2^(r-3)-1.
```

The source word is `O^r I^t z`. Its first base vertex is1-ready. Make the
standard two births x,y. Reverse exactly the r-1 arcs from y to the first
r-1 old base vertices. Every vertex and old arc is retained. The new base
order is `old0,x,old1,...,old(r-2),y,old(r-1),old tail`; its word is

```
O I O^(r-2) I O I^t z.                           (12)
```

For prefixes p=`O^r`,q=`O I O^(r-2) I O`,

```
P(q)=3P(p)=3(2^r-1),
H(p)=1, H(q)=3*2^(r-2)+3, e(p)=e(q)=1.
```

The inherited concatenation law therefore gives

```
H(qv)-3H(pv)=3*2^(r-2)+2-2H(v)-2[v begins I].    (13)
```

For v=`I^t z`, this equals2 exactly when H(z)=c_r. The rewrite preserves
the exact standard target3H+2 at every r,t and every eligible z. The source
count and the indegree-one defect at returned x are

```
Hsource=(2^r-1)(Q(z)+1)2^t+1-5*2^(r-3),
b_x=(2^r-3)(Q(z)+1)2^t+2-3*2^(r-3).             (14)
```

The second formula follows from the deletion prefix `O^(r-2) I O`, whose
state is `(2^(r-2)+1,2^r-3,1,1)`. The retained old-apex future-coordinate
drop is still `(0,2P(z)+2,2P(z),2P(z))`, independent of t.

The family is nonempty at every r. The reverse-complement partners

```
z_1=OO I^(r-3),       z_2=O^(r-3) II
```

both have count `2+3(2^(r-3)-1)=c_r`. For z_1 the source formula is

```
(2^r-1)(2^(r-1)-3)2^t+1-5*2^(r-3).
```

Thus r4 gives75*2^t-9 and r5 gives403*2^t-19. Rewrite width and defect
can grow, genuinely exceeding a fixed finite rewrite bank. The entire
union over r,t,z nevertheless has only O(sqrt X) source values up to X:
see the independently checked [size bound](../fibres/WIDE_REWRITE_SIZE.md).
This family alone cannot cover all admissible odd integers.

## 2. Prefix completeness within the declared template model

Fix a prefix p beginning OO, of length n and binary value P, followed by
`I^t z` for variable t and fixed strong z. The first standard birth inserts
a zero after the first bit, giving `B=P+L`, L=2^(n-1). Reinsert newest y
**after x and within the fixed prefix**, leaving the variable tail unchanged.
If its incidence is epsilon and ell old bits follow it, its prefix value is

```
Pnew=2B-b+epsilon*2^ell,
```

where b is the low-ell-bit value of B. Since Q(I^t z) grows while H(I^t z)
is fixed, a uniform3H+2 identity forces Pnew=3P. Since P<2L and b<2^ell,
epsilon must be0 and b=2L-P>0. Because y follows x, ell<=n-1. Reduction
modulo2^ell gives `2P=0 mod2^ell`. Also `2L-2^ell<P<2L`. Hence

```
P=2L-2^(ell-1).
```

Thus p is a run of O's followed by I's, and the tail after the inserted y
is `O I^(ell-1)`. Absorb the trailing I's into the variable leading-I tail.
This gives exactly p=`O^r` and q=`O I O^(r-2) I O`. At r2 the special
prefix q=OIIO has H5; the count equation would require H(z)=0, impossible.
At r3 it requires H(z)=2, excluded by oddness. For r>=4 it is exactly H(z)=c_r.

This completeness concerns one newest-vertex cut replacement, retained
x->y, a fixed prefix insertion slot, and a leading-I variable tail. It
does not classify slots inside a growing tail, changed old arcs, or
nontransitive deletion charts.

## 3. Mirrored continuation, but not an immediate repeated rewrite

If z ends II, the returned word ends II. Its reverse-complement begins OO,
so its last base vertex is1-ready for the mirrored two-birth operation.
In the original orientation append a vertex beating only that root and
then another beating only that newborn. This preserves old arcs and gives
another exact3H+2 exit. Thus z_2 works for every r>=4 and z_1 for r>=5.

For r4,z_2=OII the actual counts are

```
105*2^t-9 -> 315*2^t-25 -> 945*2^t-73,
```

beginning201->605->1817. Forward no-ready and converse ready are compatible:
the operation must retain its orientation gauge.

The immediate wide-prefix rewrite does not recur on either explicit core.
For z_2 the mirrored initial O-run has length two, below the required four.
For z_1 it has length r-3, also too short when r<=6. For r>=7, removing
that initial run leaves `II O^t reverse-complement(q)`. After stripping
the leading II, its strong core has length t+r+2 and count at least
`2^(t+r)+1` by the inherited bound. The next wide rule would require only
`3*2^(r-6)-1`. So this direct repeated rule fails at every r,t. Other
compensating maps remain open.

## 4. The1817 boundary and its higher-core alternative

A one-apex strong H1817 core has length at most12. Exhausting every strong
word at lengths2..12 finds none; padding changes no H. Thus a uniform
one-apex return of the r4,z2 mirrored family is impossible already at t1.

There are nevertheless two unique-start parents for a1-ready1817 child,
found by the sibling exact six-state compiler:

| Parent order | Grammar rank | Six-state |
|---|---:|---|
|10|344|(135,15,135,9360,1683,4725)|
|10|410|(193,31,241,11200,1625,7403)|

Their child value is h+E(2)-1=1817. This query has a complete order gate:
a1-ready parent has H<=909, whereas the inherited Busch minimum at order14
is1125. All parent classes through13 were included. These are rooted
classes, not an assertion about distinct unrooted classes.

Each nonroot vertex of both parents was tested as a transitive inflation
position using the independent complete run-polynomial engine. None gives
child coefficients `(1817,945,0,...)`, so none supplies the desired
945*2^t-73 family by one such module. This excludes exactly eighteen concrete
inflations, not larger charts, multiple modules, or nonmodule growth.

There is also an actual retained-source bridge at t1. In the literal labels
of `wide_results.json`, the mirrored1817 graph has13 vertices. Its subset

```
(11,0,9,8,1,2,3,5,12)
```

maps in that order to vertices `(0,2,3,4,5,6,7,8,9)` of the rank344
ready1817 child. Every induced arc agrees. This retains the old apex8 and
both mirrored newborns11,12. Delete the other four source vertices and
insert the two missing target vertices with their target cuts. The final
H is1817 and the target's newest vertex10 is1-ready. Thus the actual path
continues `201->605->1817->same-count ready1817->5453`.

An exhaustive induced-embedding search tests every10- and11-vertex subset
of both ready targets and their converses against the actual13-vertex source,
finding none. Nine retained vertices are therefore optimal for return to
a strong1-ready1817 graph in either orientation gauge. This is not an
optimality claim over2-ready graphs or arbitrary non-ready replacements.
The receipt contains the exact matrices, retained map, and the5453 child.

Every newest-vertex transitive insertion slot and either incidence were
also checked in the mirrored family for t1..15, finding no equal-count
replacement. This is FINITE-EXACT only. The all-order1817 exclusion already
pays the stated failure of a uniform one-apex return containing t1.

## 5. Reproduction

Run `python -B wide_return.py` and `python -O -B wide_return.py`.
The receipt separates all-prefix finite controls, literal graph controls,
the complete1817 one-apex universe, the inherited complete ready grammar,
and the eighteen independent module tests. The infinite claims use the
proofs above.
