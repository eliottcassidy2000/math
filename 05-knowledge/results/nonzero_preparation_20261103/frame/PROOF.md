# Terminal replacement forces a new defect coordinate

**PROVED — [independent audit](AUDIT.md) passed.** The result concerns
arbitrary terminal tournaments and arbitrary terminal port sets, at every
size, while the six-vertex core, outer root and prefix interaction remain
fixed. It is not an obstruction to changing that core or the attachment frame.

## Inheritance and the actual object

Use the actual core from
[MPG native six repair](../../multiport_grammar_20261103/component/NATIVE_SIX.md):

```
D=(2,4,57,3,11,27), w=4,
sigma(D,{w})=(17,224,28,44,0,2,4).
```

The seven coordinates are `(h,p,f,g,j,e,s)`, with ordered path slots,
including empty slots unless an endpoint is prescribed. Take a nonempty
transitive prefix A of length k, then `D -> U`. Retain the physical outer
root r, which beats A and every unmarked suffix vertex and loses precisely
to w and a nonempty subset W of U. Require W to meet the terminal SCC of U;
this is the exact strongness condition for the full attachment. Put z=2^k.

The [proved ordered-product law](../../multiport_grammar_20261103/component/PROOF.md)
gives the exact full-graph observations

```
F_r = (17+28z)h + 224z f,
H   = F_r - 44f - 17j - 2s.                       (1)
```

Thus `H` alone loses the constant part of the retained root response. This
is the relevant sidecar, not a choice of a different vertex with equal H.

The starting terminal objects K_m are the actual triangle from
[PREPARED_FOUR](../../multiport_grammar_20261103/component/PREPARED_FOUR.md),
followed by m-1 zero-letter C2 steps at the transported auxiliary, m>=1.
The old marked vertex u is retained. The inherited transport, or direct
application of its rooted identities, gives

```
h=f=N=3^m, j=0, s=N/3,
H_m(k)=3^(m-1)(756z-83),
F_m(k)=3^(m-1)(756z+51).                            (2)
```

The closest mechanism is exact C2 auxiliary transport. The hostile is its
all-auxiliary nonzero-readiness obstruction. The corrected near miss is to
enlarge a terminal port set and assume that its extra coordinates are free.
The least-used coordinate is the count s of paths beginning in that set.
The connection maps an actual terminal replacement to the affine response
(1); it retains all-clock H/F and strongness, but loses the named operation
that constructs U. Positive scalar parameters still require an actual U.

## A cut inequality that charges every port

For every tournament U and every named vertex v, cut each Hamilton path
immediately after v. This injects those paths into the two ordered covers
counted by F_v. Consequently, for q=|W|,

```
f = sum_(v in W) F_v >= q h.                        (3)
```

This holds for weak as well as strong U. The joint coordinate j and the
start count s are nonnegative. With one port, j=0 because two disjoint
nonempty path slots cannot both contain the same named vertex.

## Exact affine resource equation

Fix m and d in {2,4}. Suppose an arbitrary fixed replacement `(U',W')`
in the same frame realizes `H'(k)=3H_m(k)+d` for every legal k (or for
infinitely many k). Comparing the coefficient of z in (1) forces an
integer b with

```
h'=3N-8b,       f'=3N+b.                            (4)
```

Comparing constants then gives the complete remaining equation

```
180b + 17j' + 2s' = 2N-d.                          (5)
```

These equations are necessary and sufficient for the H identity once the
parameters come from an actual terminal object satisfying the strongness
guard. They do not assert existence of that object.

Since f'>=h', (4) implies b>=0. Equation (5) implies b<N/90. If W' had
two or more vertices, (3) would give `3N+b>=2(3N-8b)`, hence
`b>=3N/17`, a contradiction. Therefore **every such return has exactly
one terminal port**, at every m and without an order bound. In particular
j'=0, and

```
s'=N-90b-d/2,
F_r'-3F_m=-136b,
(F_r'-H')-3(F_m-H_m)=-136b-d.                       (6)
```

The changing root defect is therefore compulsory data for any nonzero b.
Adding more terminal ports cannot bypass this equation.

## Parity and the first possible positive digit

Use the inherited exact root-gap identity in
[READY, Section 1](../../two_entry_funnels_20261103/READY.md): for a single
port v in any tournament, with q=H(U'-v),

```
f'-h' = q-s' + delta,       delta is even.           (7)
```

The same cut proof gives q>=s' and delta>=0. Every tournament, including
the empty one, has odd Hamilton-path count, so q is odd. Substituting
`f'-h'=9b` and (6) yields

```
b is odd for d=2;       b is even for d=4.           (8)
```

In particular no +2 return in this frame can triple the old root response:
tripling F forces b=0, whereas (8) forces b odd. This obstruction holds
after arbitrarily many zero-channel steps and for arbitrary terminal size.

For m<=4, N<=81 and (5) forces b=0, so there is no H-only +2 return either.
At m=5, N=243, the sole arithmetically surviving +2 choice is

```
b=1,       (h',f',s')=(721,730,152).                 (9)
```

The apparent b=2 alternative `(713,731,62)` fails (8). Existence in (9)
is **OPEN**. Its root response must be `3F_m-136`, not `3F_m`.
At k=1 the source has `(H,F)=(115749,126603)` and the requested child
would have `(347249,379673)`, with root defect32424.

## The stronger response-preserving statement

More generally take any actual singleton-marked terminal object with
`h=f=N` and start count s. If a replacement in this fixed-D frame satisfies
both `H'=3H+d` and `F_r'=3F_r` at every clock, then comparison in (1)
forces `h'=f'=3N`. Equation (3) forces one port and hence j'=0. The zero-gap
identity forces s and s' to be odd, but the H equation requires
`s'=3s-d/2`. Therefore **d=2 is impossible for every such source**, not
just K_m. This includes the new clone continuation of the actual444 output.

## Next decisive object and scope

The cheapest surviving tail target is the rooted graph in (9). A positive
construction would be the first H-only +2 escape within this fixed core,
and would force precisely the defect change in (6). It would not restore
a two-child class automatically. Other promising operations must change
the older core, its marked vertex, the prefix interaction, or the frame.
The multicolour endpoint forest can evaluate those operations, but the
existence and repeatability of their actual preparations remain open.

The separate literal controls cover all labelled tournaments through five
vertices, all their nonempty port sets, zero-channel terminals through m=5,
and the stated integer resource equations. They corroborate the cut and
parity proofs; no finite census replaces the all-order arguments above.

The standalone [source](../../../../04-computation/nonzero_preparation_20261103/frame.py)
passes 51,259 always-active gates on 1,099 labelled tournaments, 32,767
nonempty port sets and 5,405 rooted cases. Normal and optimized runs have
byte-identical receipts. Reproduce from the repository root:

```text
python3 04-computation/nonzero_preparation_20261103/frame.py
python3 -O 04-computation/nonzero_preparation_20261103/frame.py --output /tmp/nonzero_frame_optimized.json
```

Receipt SHA256:
`b317cdafa85e07160b8f4a558819ac1a9d6cd86b714688e3b4207dc67465088a`.
