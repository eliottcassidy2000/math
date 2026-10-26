# All six low-order ledger failures admit three-black repairs

2026-10-26. **FINITE-EXACT:** the six failed affine source cells of
`port_residual_integral_20261025/LOW_ORDERS.md` are now repaired at their
full source scope. They were failures of a fixed menu, not connectivity
obstructions. Every source in each cell has **three successive black K
moves** to its exact labelled `K_P(T)` target.

The new finite result, combined with the audited residual `I={h,y,c}`
repairs (all valid from n=12), closes every closed marked endpoint-four
path at n=12 and n=14. The previously proved n=13 and n>=15 results then
yield the uniform **n>=12** closed endpoint-four replacement theorem.
This does not settle n=11, open marked paths, any other endpoint layer,
or full black-core connectivity. This finite packet does not claim
that 12 is the optimal threshold.

## 1. Exact universe and inheritance

Use positions `P=(0,...,n-1)=(a,b,x,...,y,c,z)` as literal vertex labels.
A K move retains the path arcs and reverses every other arc. A directed
word is blue iff its actual positional reversal is anti in its source.
The source satisfies: P is directed, blue and closed (`z->a`),
`d+(a)=n-5`, `d+(z)=4`, and

```
N-(a)={z} union I,       N+(z)={a} union {n-1-v:v in I}.
```

No strongness, module, clone, disjoint-port, isomorphism or interior filter
is imposed. The old base cells also fix `T(1,n-2)` (middle open or closed)
and `T(2,n-2)` (the x,c chord). Only closed-middle cells with both y,c in I
are outside that base bank.

| n | All four-branch cells | Base cells | Endpair cells | Old failing cells |
|---|---:|---:|---:|---:|
|12|336|322|14|2|
|14|660|642|18|4|

`check.py` independently repeats this complete base partition using the
pinned, previously audited Gaussian solver and literal old words. All 964
base ledgers are universally legal and have the correct literal target;
exactly the six stated cells can have an old blue word.

## 2. The failing symmetry makes the residual tiny

For each failed cell, impose the anti equations of its **actual offending
old word at its actual intermediate source**. The resulting affine space
has dimension 3 for each n12 cell and dimension 0 for each n14 cell:

| n | I | Middle | x,c chord | Old blue step | Residual sources |
|---|---|---|---|---:|---:|
|12|{2,7,9}|open|x->c|1|8|
|12|{6,7,9}|open|x->c|1|8|
|14|{3,5,10}|closed|x->c|3|1|
|14|{3,7,10}|closed|x->c|3|1|
|14|{5,9,10}|closed|c->x|3|1|
|14|{7,9,10}|closed|c->x|3|1|

At n12 the three independent raw arc bits are exactly
`T(1,3), T(1,4), T(1,8)`; all eight assignments occur. At n14 the
source is uniquely forced in each cell. The other two old words are
universally black in each whole original source cell.

This proves the complete conditional partition: **off** the offending
blue locus, the old ledger already consists of three black moves;
**on** it, there are exactly `8+8+1+1+1+1=20` literal tournaments.
The checker enumerates every solution of the Gaussian systems and
checks equality with the 20 supplied source rows, not merely inclusion.

## 3. Cyclic-shift bridge: an exact operation behind the certificates

This lemma needs only that P is a directed closed Hamiltonian word,
not endpoint scores or blue status. Put `B=K_P(T)` and
`C=(1,2,...,z,0)`. If W is a directed Hamiltonian word of the induced
`tournament B[V minus {0}]`, beginning at1 and ending atz, set

```
Q1=C,           Q2=(W,0),           Q3=(0,W).
```

These three words are directed in their successive sources and reach
B exactly. Indeed, `E(C) xor E(P)={01,0z}`, so `S=K_C(T)` agrees with
B away from those pairs. In S, `z->0` and `1->0`, while the W arcs
are directed. Hence Q2 is legal. Complementing its omitted pair01
makes `0->1` true; retaining all W edges makes Q3 legal. Finally,

```
E(Q2) xor E(Q3)={01,0z},
E(Q1) xor E(Q2) xor E(Q3)=E(P).
```

Three full complements and this XOR identity prove the literal target.
**Blackness is a separate obligation**: this lemma alone does not give a
black detour. Here every required blackness test is checked exactly.

## 4. Six short W recipes

The following table is used only on the offending-blue loci of section2.
Let `u=T(1,3)` and `v=T(1,4)` at n=12. The third free bit does not affect
the choice. Use the bridge words from section3.

| Source cell and actual bits | W |
|---|---|
|Either n12 cell, u=0|(1,3,4,2,5,6,7,8,9,10,11)|
|Either n12 cell, u=1,v=0|(1,4,5,6,7,2,3,8,9,10,11)|
|n12 I={2,7,9}, u=v=1|(1,6,7,8,2,3,4,5,9,10,11)|
|n12 I={6,7,9}, u=v=1|(1,6,2,3,4,5,7,8,9,10,11)|
|n14 I={3,5,10} or {3,7,10}|(1,3,4,5,2,6,7,8,9,10,11,12,13)|
|n14 I={5,9,10} or {7,9,10}|(1,3,4,5,6,2,7,8,9,10,11,12,13)|

There are ten affine recipe subcells: dimensions `(2,1,1)` for each
n12 port cell and four singletons at n14. `check.py` also proves universal
legality, all three blackness conditions and labelled target equality
symbolically in every one of these ten cells. The full 20-row replay is
an additional exact path, not a random spot check.

## 5. Literal ports, flips and clone transport are retained

`found.json` is the compact 20-source certificate. `normal.json` and
`optimized.json` contain, for every source and every successive move:

- the literal Hamiltonian word and complete intermediate tournament rows;
- every flipped unordered pair;
- a specific pair whose two anti-reversal arc bits disagree;
- all labelled 0/z neighborhoods and both actual word-endpoint neighborhoods;
- all global twin pairs, internally ordered maximal transitive modules,
  and complete block-intersection matrices between successive sources.

All 80 states along these 20 detours happen to be clone-free. This is
checked, not assumed. No vertex is contracted and no path is lifted from
a quotient. A positive merger/split control checks
`(6,4,0) <-> (2,4,1)` under the word `(0,1,2)`: a transitive 3-block
merges/splits exactly while a cyclic triangle is never misidentified as
a transitive clone. Maximal transitive modules are recovered as connected
components of the global twin graph and each component is independently
checked to be a transitive module. Completeness follows because adjacent
vertices in the internal order of any transitive module are global twins.

The genuine order-seven hostile is retained:

```
T=(6,4,8,19,103,79,15),
K_P(T)=(122,124,121,112,32,64,16).
```

Its target still has exactly nine Hamiltonian words, all blue. Thus the
checker distinguishes a real isolated-black target from a failed menu.
Reversing a directed candidate is rejected, and one cyclic shift alone
is checked not to equal the desired target.

## 6. Full n12/n14 consequence and inherited boundary

The remaining 14/18 cells have `I={h,y,c}`, with h=2..n-4 and either x,c
chord. The existing `port_residual_integral_20261025/PORT_FRONTIER.md`
section3 exhausts these cells without an extra guard after complementary
cases are combined:

- h=2: `ports_owners_f1_20261025/endpoint.md`, section2, valid n>=11;
  its x->c subcase is `endpoint_four_frontier_20261021/README.md`, section3.
- h=3,3->1: `unprotected_ports_cancellation_20261025/ports/ENDCAP_NEXT.md`.
- h=3,1->3: `port_residual_integral_20261025/h3/README.md`.
- h=4,1->3: `unprotected_ports_cancellation_20261025/ports/README.md`.
- h=4,3->1: `port_residual_integral_20261025/h4/README.md`.
- h>=5: `unprotected_ports_cancellation_20261025/ports/README.md`.

All but the stronger h=2 result start at n=12. They are complete labelled
three-black replacements, not partial ledgers spliced after a different
operation. Thus this packet's repaired base and those inherited residual
families exhaust n12 and n14. Together with the already certified n13
head and n>=15 tail, the closed endpoint-four replacement threshold is
now n>=12. Full black-core connectivity is not a consequence.

## 7. Reproduction and dependency provenance

```
python check.py > normal.json
python -O check.py > optimized.json
cmp normal.json optimized.json
```

Both runs return PASS and are byte-identical. The checker uses explicit
exceptions so optimized Python does not disable its gates. It imports
two byte-identical copies of the inherited independent audit files,
verified by their previously published SHA256 values before execution.
This is pinned reuse of that solver, not a claim of a newly independent
Gaussian implementation. The new literal row replay and solution census
are local. The independent review, when attached, is a separate path.

Exploration scripts and the separate `next11/` diagnostic are not part of
this result's proof or stated threshold.
