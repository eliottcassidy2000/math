# Independent live-port, rank, and fresh-block audit

**PASS: the all-order two-live-port carrier and width equivalence; exact
rational Hankel rank17; the stated fresh-block formulas. FINITE-EXACT:
the explicitly bounded controls below.** The independent audit did not rerun
the target saturation through999. It imports no producer or repository Python
module. SymPy is used only for exact rational linear algebra.

## 1. The finite carrier proof survives arbitrary port forgetting

Retain exactly two names. A command chooses one as root, creates x with
root->x->every other old vertex, and discards one of the three names. A
discarded vertex never becomes a root again, so every future newborn points
to it. Hence a path-cover component whose end has been forgotten can never
acquire a future outgoing edge. Two such ends cannot occur in the restriction
of a final Hamiltonian path. Discarding these signatures is sound.

This is a necessary pruning condition, not a claim that every surviving
signature extends under every continuation. The transfer still enumerates
actual covers exactly: choose at most one old component entering x and one
leaving x, reject choosing the same component twice, and then forget the
discarded name. Component pairing prevents cycles; multiplicities count
different choices with the same endpoint signature. Uniform adjacency to
forgotten vertices removes their identities from future counts.

With two names, there are at most three surviving components. The independent
component-chain implementation closes the actual C3 seed support to exactly
24 signatures and reconstructs all six producer matrices entry for entry.
It verifies every full cover vector for all43 programs of lengths0,1,2,
using literal permutations and all cut sets, deduplicated as unordered
collections of actual directed paths. This exercises forgetting either old
port and forgetting the newborn, including paired endpoints and multiplicity.

## 2. Future-live width is exactly the two-port condition

For a finite chronological root word `(r_0,...,r_(m-1))`, put

`L_t={r_i : i>=t} intersect {vertices present before birth t}`.

Necessity: each member of L_t is needed as a future root and therefore cannot
already have been forgotten. The two retained names contain L_t.

Sufficiency: initially choose two names containing L_0, filling any empty
place arbitrarily. If the retained set P_t contains L_t, the required root
r_t is available. Moreover `L_(t+1) subset L_t union {x_t}`. Since
`|L_(t+1)|<=2`, discard one member of the three-name pool while retaining all
of L_(t+1). Fillers are harmless. Induction realizes the whole word.

The three directed pairs `(0,1),(1,2),(2,0)` in C3 are related by cyclic
automorphisms. Their independently enumerated initial signature vectors
agree, so one oriented retained pair suffices. Port slots are ordered and the
initial arc goes from slot0 to slot1; replacing an old slot by the newborn
preserves the slot convention without asserting a later arc orientation.

The unique ancestry theorem in
[ANCESTRY A1-A3](../adaptive_metagraph_reformulation_20261103/ANCESTRY.md)
recovers the root word up to the original C3 rotation. This rotation preserves
every |L_t|. Thus maximum future-live width is an ordinary-class invariant
within the triangle-seed adaptive language.

The degree and width restrictions are incomparable:

* `(0,1,2)` has future-live width3 and root indegrees `(1,2,3)`. Its six-vertex
  graph has rows `(10,20,33,6,13,27)` and literal H37. Its class is cap3 but
  outside the two-port grammar, by unique ancestry.
* `(0,3,...,3,0)`, with k repetitions of3, has width2 and final root
  indegree k+1. Keep names0,3 after the first birth. This proves unbounded
  root degree inside the two-port grammar, not just a larger finite example.
* The earlier H127 history `(0,0,3,3,0)` has width2 and root degrees
  `(1,1,2,2,3)`. All8! permutations independently reproduce127 on rows
  `(154,4,1,102,14,23,55,126)`.

## 3. Independent exact rank certificate

The producer lists24 actual prefixes and17 actual continuations. The audit
rebuilds their vectors from its own matrices and recounts **all408 Hankel
entries directly on actual graphs**, using backward path-completion recursion
instead of the producer's forward subset-last recurrence. There are254
distinct labelled graphs in this check, of order at most10.

Over the prime field F_1000003, independent elimination gives reachable
rank24, observable rank17 and Hankel rank17. For the upper bound, let C be
the24-by-17 matrix of those actual continuation columns. The audit selects
an invertible17-row submatrix and verifies, over Q and entry by entry,
`M_i C=C Q_i` for every one of the six letters. Its first column is the
Hamiltonian observer. Thus all future observations lie in this17-dimensional
space. The modular lower bound and exact upper bound prove rational
reachable rank24, observable rank17 and minimal scalar-series dimension17.

A fully actual17-by-17 Hankel minor has determinant

`33554432 = 2^25`.

Its prefix-row indices within the producer list are
`0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,18`; it uses all17 suffixes.
This is stronger evidence than a modular rank alone. The carrier still
contains unbounded integer weights; rank17 is not17 possible states or a
membership algorithm for prescribed H.

**Every-characteristic corollary.** The same scalar series has rank17 over
every field of odd characteristic, and rank1 in characteristic2. Here is the
localization argument, rather than an extrapolation from the one modular
test. Let B be the17 selected actual prefix rows and K=BC the displayed
Hankel minor. Then

`A_i=K^(-1) B M_i C`

has all entries in Z[1/2], since det K=2^25 and B,M_i,C are integral. The
already audited rational observable closure gives `M_i C=C A_i` exactly.
Both sides lie in Z[1/2], so these identities, the integer initial row uC,
and its first-coordinate H observer reduce over every odd characteristic.
They give an upper bound17. The same actual Hankel minor stays nonzero in
each such characteristic, giving the lower bound17. Extension of the base
field changes neither bound. In characteristic2, every program constructs
an actual tournament and its H is odd by
[THM-001, redei](../../../01-canon/theorems/THM-001-redei.md).
Thus the complete Hankel matrix is identically1, and its rank is exactly1.

## 4. What the earlier rank17 does and does not explain

[October30, Exact minimal linear memory is seventeen](../square_sum_spectrum_overnight_20261030/README.md)
uses a different operation. Two permanent marked vertices x->y are removed
to leave a transitive chain. Each of four letters0,1,2,3 specifies which
marked vertices point to a new chain vertex. Its positive gap-assignment
carrier has36 coordinates: two orders of x,y, three states in each of the
first two gaps, and two in the last. Its rational reachable/observable ranks
are28/20, with scalar Hankel rank17. The old actual Hankel minor is
`2^24*3^2`; the new minor above is different.

There is a genuine **proper common sublanguage**. Never discard either of
the two original live ports, so use only live letters2 and5. Reverse the
chronology of all unmarked vertices to make them an old-style transitive
chain. Root0 gives old incidence letter1, root1 gives letter2, and the
initial third C3 vertex supplies the final letter2. Thus a live program's
old word is the reversed sequence of these incidence letters followed by2.
This is a literal graph relabelling, not just a matching count. The audit
checks all15 words of lengths0 through3 as adjacency matrices.

On this common subseries, exact closure over F_1000003 gives Hankel rank13
in both representations. Their reachable/observable ranks are13/15 in the
24-coordinate live model and15/20 in the reversed36-coordinate static model.
These are exact finite-field statements; no rational-rank13 upper bound is
claimed by this audit.

At this first stage, the common subseries did **not** identify the full
operation languages or explain why both full rational ranks are17. Port
replacement leaves the permanent-feedback model, and the newborn-only live
sublanguage has unbounded static feedback cost. A universal common boundary
quotient needed an explicit intertwining map, compatible operations, initial
state and output; equality of dimensions was insufficient.

**Subsequent refinement:** the coordinator's twelve-command two-port envelope
provides that missing response-quotient connection, independently audited in
[COMMON_OBSERVER.md](COMMON_OBSERVER.md). Its25-coordinate carrier covers
arbitrary incoming masks at the two ports; seventeen static continuation
queries close exactly under all twelve commands, including port replacement.
The full old static language embeds after reversing chain chronology, and
the C3 six-command language embeds after a literal triangle entrance. This
explains the shared response dimension17 while preserving the distinction
between the full graph languages. The rank13 subseries calculation above
remains valid, but is no longer the strongest connection available.

## 5. Fresh transitive block: independent formula check and proof

For `(h,s,q,F)` with the conventions of
[ADAPTIVE section2](../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md),
adjoin k transitive twins at the old root v: v points to every twin and every
twin points to T-v. Later twins point to earlier twins. Mark a new vertex w
with A newer and B older twins, so k=A+B+1. Set `a=h-s`, `d=F-q`.

The audited formulas are

```
h' = a+q+(2^k-1)d
s' = 2^B d                         if A>0
     a+(2^B-1)d                    if A=0
q' = h                             if k=1
     a+q+(2^(k-1)-1)d              if k>1
F' = 2^A a +(2^(B+1)3^A-2^A)d +2^(A+1)q +1_{B=0}(s-q).
```

The q' formula follows by deleting w, leaving k-1 identical twins. For s',
an old first path ending at v must be nontrivial when preceded by w; these
old two-slot covers number d. If A>0, the newer twins force a nonempty
post-v twin run, while each older twin chooses the initial or post-v run,
giving2^B d. If A=0, the case with no post-v twin run contributes a, and
the other assignments give `(2^B-1)d`.

For F', count ordered covers `(P,Q)` of the enlarged graph, with P ending
at w. Old vertices have no edges into the twins except from v, so the twin
runs are forced into the following three disjoint cases:

| Case | Contribution |
|---|---|
| v belongs to P | `2^B 3^A d +2^A q` |
| v belongs to Q, with no twin run after v | `2^A a +1_{B=0}s` |
| v belongs to Q, with a nonempty twin run after v | `(2^B 3^A-2^A)d +(2^A-1_{B=0})q` |

In the first case each newer twin can occupy P's initial run, P's final
run before w, or Q's initial run; an older twin has only the first and third
choices. If the old first path is just v, its initial run is forbidden,
giving the q term. In the third case the three choices are P (before w),
Q before v, or Q after v; subtract the assignments with an empty post-v run.
If Q's old first path is just v, its initial run is forbidden; the all-P
assignment must be subtracted exactly when B=0. Summing the table proves
the displayed F' formula. Empty second paths are included throughout.

The h' formula is the inherited fixed-root block count and also follows
from the same forced-run decomposition. At k=1 the four formulas reduce
exactly to the previously proved newborn map N. At A=1,B=0 they imply
`h'-q'=2s'`, so this individual command has singular image; this does not
force the union of all block commands into dimension three.

The independent block script verifies all four coordinates by direct
subset/end-path recursion, with no producer imports. Its240 controls are
every rooted labelled tournament of orders1 through3, regular C5 at every
root, and the adaptive H37 six-vertex source at every root; for each it
tests k=1,2,3 and every new mark. Maximum child order is9. This finite
verification is separate from the all-order combinatorial argument above.

## Reproduction

Run each audit normally and with `python3 -O`. Both pairs agree byte for byte.

* `live_port_independent_audit.json` SHA256:
  `7f1ada92fc09d700db7c397b4b4ec40127e702737fdc642ff01ed1a8b289930d`.
* `block_switch_independent_audit.json` SHA256:
  `9c69b041d1a9f29d22284277facd819fe0643d6c56fb239ea76aac8c50afb559`.

The live-rank audit reads and hashes the frozen `live_port.json`, treating
its listed words as proposed certificates and independently checking them.
It does not verify or rely on the999-support enumeration. The scripts never
use `assert` for a load-bearing check, so optimized replay preserves checks.
