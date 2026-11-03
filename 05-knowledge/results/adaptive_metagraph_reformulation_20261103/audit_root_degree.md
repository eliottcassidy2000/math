# Independent audit of bounded-indegree adaptive roots

**PASS. FINITE-EXACT search and witness checks; PROVED all-size restricted
nonattainment and active-root bound. No universal cap3 or unrestricted
spectrum theorem.** This audit imports no producer implementation. It uses
Boolean adjacency matrices and a memoized recursive Hamilton-word counter,
whereas `root_degree_probe.py` uses outgoing bitsets and forward subset DP.
The literal127 check uses all8! permutations as a third count method.

## Scope and outcome

Start at the named triangle `0->1->2->0`. A birth at an old root v adds x
with `v->x` and `x->u` for every other old vertex u. The cap d restricts the
root's **indegree immediately before the birth** to at most d. It does not
bound its eventual degree or the final tournament's maximum indegree.

The independent search checks the complete cap2 language with H<=199 and
retains every named history. It tests1,699 candidate children. Its layers
agree with the producer exactly:

| Parent order | Surviving parents | Surviving children |
|---:|---:|---:|
|4|1|4|
|5|4|14|
|6|14|42|
|7|42|122|
|8|122|262|
|9|262|128|
|10|128|0|

The attained counts are precisely the odd values from3 through199 excluding
7,21,127. Thus cap2 does not suffice for the ordinary prime lane: **127 is
absent at every order of this restricted language**, not merely through a
fixed vertex order. This follows from complete count-bounded pruning and
strict increase, proved below. It does not assert that127 is globally absent.

Three further searches independently use first roots0,1,2 at cutoff127.
Each tests767 children, has identical layers, and misses127. This directly
checks the hostile without relying on the producer's first-root gauge.
The explicit C3 rotations also prove that gauge at all later depths.

All312 saved witness histories across `root_degree_99.json`,
`root_degree_199.json`, and `root_degree_3_199.json` are reconstructed from
scratch. Recursive counts, root caps, tournament strongness and target bounds
all pass. In particular all97 allowed odd targets in[3,199] have actual cap3
witnesses. This positive assertion needs only the supplied witnesses; this
audit does not independently repeat the producer's complete cap3 search.
The cap1 cutoff99 receipt is checked for its positive witnesses and internal
missing-list consistency, not independently searched here.

## Literal127 witness

The root history is `(0,0,3,3,0)`. Root indegrees before those births are
`(1,1,2,2,3)`. Thus its last step uses precisely the extra degree allowed by
cap3. The final outgoing bitsets are

`(154,4,1,102,14,23,55,126)`.

They represent a strong8-vertex tournament. Literal enumeration of all40,320
vertex permutations counts exactly127 directed Hamiltonian paths, agreeing
with the independent recursive counter. The audit also checks127 is prime
by trial division through11. This witness shows cap3 can repair the specific
cap2 failure, not that cap3 always suffices or that every127 witness has this
ancestry.

## Why the finite rejection is all-order in the restricted language

Every birth from a strong tournament remains strong. Old vertices reach v
inside the old tournament, then reach x through v. The new vertex directly
reaches every old vertex except v, and it reaches v through any old
in-neighbor of v. Such an in-neighbor exists because the parent is strong.

For strict H growth, retain every old odd-cycle collection and note that an
old in-neighbor u of v supplies the new directed triangle
`v->x->u->v`. The odd-cycle-collection identity
[THM-002 / ocf](../../../01-canon/theorems/THM-002-ocf.md) therefore gives
`H(child)>=H(parent)+2`: the old collections contribute H(parent), and the
singleton collection consisting of this new triangle contributes2. This
argument uses the established OCF theorem, not a numerical growth fit.

Consequently an over-cutoff history can never return below the cutoff.
Every later history below the cutoff must descend from one of the retained
children. Once a whole next layer is empty, no further such history exists.
This validates the exhaustive rejection of127 at arbitrary order. The
search does not rely on a guessed vertex cap.

The triangle rotations `i -> (i+r) mod3` send first root0 to root r, fixing
the first newborn vertex. Birth commutes with relabeling, and indegree and H
are relabeling invariant. Extending that permutation by fixing later newborn
names gives a bijection on all later root histories. The script checks the
initial4-vertex isomorphisms explicitly; the equivariance argument gives the
all-depth normalization.

## All-size active-root theorem

For an integer d>=0, write `A_d(T)={v: indeg_T(v)<=d}`. In every strong
post-birth tournament in this language,

`|A_d(T)| <= 2d`.

More generally this holds after a birth from any strong tournament, provided
the resulting order is at least4. The qualification matters: the C3 seed has
`|A_1|=3>2`.

**Proof.** If k=|A_d|>0, the arcs internal to A_d contribute exactly
`k(k-1)/2` to its total indegree, which is at most kd. Thus k<=2d+1. Equality
k=2d+1 forces every member to have indegree d and permits no incoming arc
from outside A_d. If A_d is proper, all arcs go from A_d to its complement,
contradicting strongness. If A_d is the whole tournament, it is d-regular.
But the last newborn vertex has indegree1, forcing d=1 and order3, excluded
by the post-birth hypothesis. For d=0 strongness directly gives A_0 empty.
This proves the bound. The first4-vertex birth has two vertices each of
indegree1 and2, giving equality for d1 and d2. QED.

Old indegrees obey the exact update

`indeg_child(u)=indeg_parent(u)+[u!=v]`,

while `indeg_child(x)=1`. Therefore an old vertex that has left A_d can never
return; its indegree never decreases. Only a new vertex can replenish the
active-root set. The independent script checks both statements at every
visited state/birth and checks the active bound for every d through the
current order. Those tests are controls; the proof is all-size.

The bound controls the number of roots eligible at one step, not the amount
of continuation information needed. The selected root preserves its
indegree, other old roots age by one, and names can matter to future H.
Nothing here identifies a closed2d-state automaton or bounds the memory
needed to update Hamiltonian counts.

## Reproduction and independence

Run from the repository root after the three producer receipts are in the
packet:

```sh
python3 -B 05-knowledge/results/adaptive_metagraph_reformulation_20261103/audit_root_degree.py > /tmp/root-degree-audit.json
python3 -B -O 05-knowledge/results/adaptive_metagraph_reformulation_20261103/audit_root_degree.py > /tmp/root-degree-audit-optimized.json
cmp /tmp/root-degree-audit.json /tmp/root-degree-audit-optimized.json
cmp /tmp/root-degree-audit.json 05-knowledge/results/adaptive_metagraph_reformulation_20261103/audit_root_degree.json
```

For the initial audit the receipts resided in `/tmp`; the equivalent runs
used `--receipt-dir /tmp`. Output records only receipt basenames and content
hashes, so relocating identical input files does not change the receipt.
All mathematical gates use explicit exceptions and survive `python -O`.
The output is deterministic. No large tournament census was performed.

The producer's search logic was also read: it preserves actual adjacency,
checks the cap before each birth, retains all later histories, and uses strict
growth for cutoff pruning. The independent audit shares only the mathematical
operation and saved target words, not transition or counting code. Its
complete search and first-root controls cover the load-bearing negative
claim; literal permutations cover the main positive hostile.
