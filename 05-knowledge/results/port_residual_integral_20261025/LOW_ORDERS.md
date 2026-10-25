# Low-order audit of the inherited endpoint-four base bank

2026-10-25. **FINITE-EXACT:** at order13, all464 complete affine source
cells in the three inherited base cases admit the unchanged three-black
ledger. At order12, two cells admit a blue first word; at order14, four
cells admit a blue last word. The six second-word budget systems all
exclude three incoming middle ports.

The order13 consequence is submitted for independent review. This packet
does not alter [PORT_FRONTIER.md](PORT_FRONTIER.md), lower an inherited
all-order threshold, or treat a failed ledger as a connectivity obstruction.

## 1. Inheritance, exact question, and pinned sources

The closest mechanism is the audited
[endpoint-four base partition](../endpoint_four_repairs_20261021/README.md),
with its24 affine second-word budget systems at n15..26. The canonical
hostile is a legal word whose actual positional reversal is blue. The
corrected near miss is assuming that the original cyclic-cut bound
n>4k-2 may be lowered without retesting its actual source guards.
The least-used sidecar is the complete affine family at one small order:
it can certify every interior completion without enumerating tournaments.

The board is exact incoming ports / middle orientation / chord x,c /
first-word symmetry / second-word port budget / final-word symmetry.
This audit asks only whether the unchanged inherited base ledgers work
at n12,13,14. It does not search for replacement Hamiltonian paths.

The checker imports two independently audited sources, both checked by
SHA256 before execution:

| Source in endpoint_four_repairs_20261021/audit/ | SHA256 |
|---|---|
| independent_audit.py | 376782e9c9eab3f3dae563a22735271d064fe73f4ef0361de044beb56cdc85b0 |
| audit_core.py | 9e0f86d30773b46dad359c04ead55dc07e606cad54156fb4eb0025b42a7339a4 |

The first supplies the literal old words and the independently reviewed
minimum-port systems; the second supplies Gaussian elimination over F2.
The new checker independently performs each complete-source replay and
implements its own integer-row literal replay for hostile witnesses.
It is a pinned reuse of the earlier independent solver, not a claim to
be a third independently implemented elimination algorithm.

## 2. Full source universe and logical meaning of a cell

At each n in{12,13,14}, use the standard marked word P=(0,...,n-1),
a=0,b=1,x=2,y=n-3,c=n-2,z=n-1. A raw independent bit for every unordered
pair encodes its forward arc. Impose exactly:

- every P path arc and closing arc z->a;
- the full P anti-reversal equations;
- a chosen three-element set I from{2,...,n-2}, with
  T(a,v)=1 exactly when v is outside I, for1<=v<=n-2;
- the chosen middle bit T(b,c) and chosen chord bit T(x,c).

There is no strongness, mixedness, clone, module, port-disjointness,
isomorphism or extra interior filter. Every source system is explicitly
checked for consistency. These data fix endpoint scores(n-5,4) and
partition the full marked source universe: choose I and the two arc bits.

Only cells with a **closed middle and both y,c in I** are deferred:
these are exactly the I={h,y,c} residuals covered by the other packets.
Thus the old three base cases being tested are every open middle,
every closed middle with c outside I, and every closed middle with
c in I but y outside I. The optional fourth guarded case of the old
base theorem is unnecessary and is not used.

For each selected ledger, symbolic joins must all reduce to the
constant1. At each actual intermediate source, adjoining the complete
anti-reversal equations for that word must be inconsistent to certify
blackness for every completion. The path-edge XOR and full symbolic
labelled final matrix are also checked. Hence a passing cell proves a
universal statement about all its completions. A failing blackness
test says only that at least one completion has that word blue.
The checker constructs such a completion and replays it literally.

## 3. Exact results

| n | All source-branch cells | Base cells tested | Closed-endpair cells deferred | Base cells admitting a blue word |
|---|---:|---:|---:|---:|
|12|336|322|14|2|
|13|480|464|16|0|
|14|660|642|18|4|

Every displayed ledger in all1428 tested cells is uniformly directed,
has the correct edge XOR, and reaches the literal target. The only
failures are blackness failures listed below.

For the additional second-word systems, the affine minimization projects
T(a,v),2<=v<=n-2, without fixing I. Each free parity component contributes
min(p,q) incoming ports, plus any forced incoming constants. Independent
root choices attain that minimum; the inherited solver constructs and
checks a minimizing full assignment.

| n | Open penultimate Q2 system | Closed one-chord Q2 system |
|---|---|---|
|12|Inconsistent|Inconsistent|
|13|Exact minimum4; example I={3,6,7,10}|Exact minimum4; example I={4,5,6,7}|
|14|Inconsistent|Inconsistent|

The n13 four-port minimizers are diagnostics for the relaxed systems,
not endpoint-four sources: the latter require exactly three middle
incoming ports. Thus Q2 is black throughout all six relevant source
families. The old proof's obstruction at12 and14 is elsewhere.

### Complete order13 consequence

The464 base cells exhaust every source except the16 closed-endpair
cells I={h,10,11},2<=h<=9, with both values of T(2,11). The proved
second-port theorem handles h=2 from n>=11; the paired new and inherited
h=3/h=4 repairs and the h>=5 theorem handle all remaining h from n>=12.
These inherited repairs impose no additional chord restriction after
their complementary cases are combined. Their exact domains and reviews
are linked in [PORT_FRONTIER.md](PORT_FRONTIER.md), section3.

Consequently the full order13 closed endpoint-four claim follows from
this complete finite certificate plus those proved residual repairs:
every such marked path has three black moves to its exact labelled
K_P(T). This consequence awaits the requested independent review before
the maintained frontier is updated. The n>=15 proof is unchanged.
The isolated order13 conclusion does not imply a threshold n>=13,
because order14 still needs a different repair in the cells below.

## 4. Explicit hostile witnesses to the unchanged ledgers

The following are the entire failed-cell list. Each cell has the stated
I, middle bit and x,c chord. These are failures of the old words,
not examples of disconnected black components.

| n | I | Middle | x,c chord | Possible blue word |
|---|---|---|---|---|
|12|{2,7,9}|open|x->c|Q1 of open_penultimate|
|12|{6,7,9}|open|x->c|Q1 of open_penultimate|
|14|{3,5,10}|closed|x->c|Q3 of closed_one|
|14|{3,7,10}|closed|x->c|Q3 of closed_one|
|14|{5,9,10}|closed|c->x|Q3 of closed_zero|
|14|{7,9,10}|closed|c->x|Q3 of closed_zero|

At n12, I={2,7,9}, a literal source is

    (1402,3748,1945,2674,1698,2372,3478,2857,2586,1121,2472,533).

Its inherited words are

    Q1=(0,1,10,11,2,3,4,5,6,7,8,9),
    Q2=(1,10,0,11,2,3,4,5,6,7,8,9),
    Q3=(1,2,3,4,5,6,7,8,9,10,0,11).

All three are directed, have path XOR E(P), and reach K_P(T).
Their actual blue flags are(true,false,false). Thus merely extending
the old first-word claim to n12 is false.

At n14, I={3,5,10}, a literal source is

    (7126,12116,13352,7507,12452,13643,8852,
     13614,2644,13500,2385,13044,9538,1289).

Its inherited words are

    Q1=(12,1,13,0,2,3,4,5,6,7,8,9,10,11),
    Q2=(13,12,1,0,2,3,4,5,6,7,8,9,10,11),
    Q3=(0,13,1,2,3,4,5,6,7,8,9,10,11,12).

All three are directed, have path XOR E(P), and reach K_P(T).
Their actual blue flags are(false,false,true). The other n14 branch
uses the same final word; the old common-final lemma has a genuine
ledger obstruction here. The receipts give explicit literal rows and
all three words for every one of the six cells.

The first failed implication is therefore the claimed absence of the
displayed anti-reversal, not legality, target equality, or second-word
port budget. The strongest survivor is complete base closure at n13,
and all but those six affine cells at n12/n14. The missing datum for a
new repair is the actual feasible symmetry of the offending word;
a new operation would have to use or break it.

## 5. Reproduction and limits

From the repository root:

    python 05-knowledge/results/port_residual_integral_20261025/low_orders.py
    python -O 05-knowledge/results/port_residual_integral_20261025/low_orders.py

Both runs passed and produced byte-identical
[normal](low_orders.normal.json) and
[optimized](low_orders.optimized.json) receipts. All gates are explicit
exceptions. This is an exact affine-family audit at three fixed orders,
not an inference about infinitely many orders and not a large tournament
or Hamiltonian-path census. No source or status file from an earlier
packet was modified.

