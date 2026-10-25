# Every closed endpoint-four line is repaired from order fifteen

2026-10-25. **PROVED by audited inheritance and the new residual repairs:**
every directed, closed, blue marked path with endpoint outdegrees
(n-5,4), at every n>=15, has three successive directed black K moves
to its exact labelled K target. The inherited proof includes a
**FINITE-EXACT** 24-system certificate at orders15..26. There is no
strong-target, mixedness, module, disjoint-port, or isomorphism filter.
Full black-core connectivity across all endpoint layers remains **OPEN**.

## 1. Exact statement

Let T be any tournament on n>=15 vertices. Let
P=(a,b,x,...,y,c,z) be a directed Hamiltonian path such that z->a and
its actual positional reversal is an anti-automorphism of T. Assume

    d_T^+(a)=n-5,          d_T^+(z)=4.

For a directed word Q, let K_Q retain its path arcs and reverse every
other arc. Call Q black when its own positional reversal is not anti
in its actual source. There exist directed words Q1,Q2,Q3, each black
in its successive source, such that

    K_Q3(K_Q2(K_Q1(T))) = K_P(T).

The equality concerns every literal vertex and arc. No clone contraction
or quotient lifting is performed. The result includes arbitrary clone
configurations and shared or reflected ports.

Use P positions as labels0,...,n-1, and put
a=0,b=1,x=2,y=n-3,c=n-2,z=n-1. There is exactly one three-element set I
of interior incoming ports with

    N^-_T(a)={z} union I,
    N^+_T(z)={a} union r_P(I).

The path excludes b from I; the possible members of I are2,...,n-2.
The middle is called open if b->c and closed if c->b. This distinction
is independent of the hypothesis that the full path P is closed.

## 2. Inheritance audit: the earlier bank really is unfiltered

The decisive predecessor is
[endpoint_four_repairs_20261021/README.md](../endpoint_four_repairs_20261021/README.md),
not a chronological frontier summary. Its theorem gives three black
moves, for every n>=15, in all three cases:

1. The middle is open, with arbitrary I.
2. The middle is closed and c is not in I.
3. The middle is closed, c is in I, and y is not in I.

The producer explicitly declares its source universe unfiltered.
Its [independent audit](../endpoint_four_repairs_20261021/audit/AUDIT.md)
approves these same quantifiers and checks that the affine universes
contain no strongness, mixedness, module or isomorphism condition.
In particular I and r_P(I) are allowed to overlap.

The inherited operations were checked at their actual scope. The
[arbitrary-k open-middle ledger, R2](../endpoint_three_repairs_20261020h/root.md)
requires n>4k-2; here k=4 gives exactly n>=15. The
[degree-free movable-cut corollary, section9](../endpoint_three_repairs_20261020h/theory.md)
requires c->b,z->b,a->y,d_T(z)>=3. These follow respectively from
closed middle, c in I, y outside I, and endpoint degree4. Its cut is
at the first actual outgoing port after b, so no extra fixed-port
condition is concealed.

The two other middle-word tests in the base are proved by three
disjoint endpoint-row parity chains for n>=27. They would force at
least four incoming middle ports, contradicting |I|=3. At n15..26
the base uses24 complete affine systems and minimizes the exact
incoming-port count, without first fixing I. The independent Gaussian
audit constructs minimizers and verifies the same contradiction.
Thus the finite head is an explicit part of the inherited proof,
rather than an extrapolation from sampled tournaments.

The later
[odd-prefix frontier](../endpoint_four_frontier_20261021/README.md),
[cyclic-prefix theorem](../trijection_hourglass_ports_20261025/endpoint/THEOREM.md)
and [port-bank frontier section5](../ports_regular_actions_hn_20261025/ports/README.md)
successively reduce the residual to a closed middle with
I={h,y,c},3<=h<=n-6,c->x,y->x and failed prefix guards. That statement
is consistent with the earlier theorem. However, none of those last
chord or failed-guard conditions is needed for the final exhaustion
below. In particular the cyclic-prefix theorem's actual guards
3->c and (n-h)->x are not silently dropped.

## 3. Exhaustion of the last possible port sets

If none of the three base cases applies, the middle is closed and
both y,c belong to I. Since |I|=3 and b is excluded,

    I={h,y,c},           2<=h<=n-4.

Every possible h is now covered:

| Remaining port h | Complete repair and actual scope |
|---|---|
| h=2 | [Second-port theorem](../ports_owners_f1_20261025/endpoint.md), section2: every closed-middle I={x,y,c} source at n>=11; it splits the relevant chords and retains both directions of x versus y. |
| h=3, 3->1 | [Two-triangle endcap repair](../unprotected_ports_cancellation_20261025/ports/ENDCAP_NEXT.md): all such endpoint-four sources at n>=12. |
| h=3, 1->3 | [New h3 theorem](h3/README.md): all such endpoint-four sources at n>=12. |
| h=4, 1->3 | [Unprotected-port theorem](../unprotected_ports_cancellation_20261025/ports/README.md), section1a: all such sources at n>=12. |
| h=4, 3->1 | [New h4 theorem](h4/README.md): all such endpoint-four sources at n>=12. |
| 5<=h<=n-4 | [Unprotected-port main theorem](../unprotected_ports_cancellation_20261025/ports/README.md): all such endpoint-four sources at n>=12, with arbitrary middle. |

A tournament has exactly one direction on the pair{1,3}, so the two
h=3 cases and the two h=4 cases are exhaustive. All order thresholds
are satisfied by n>=15. The [unprotected/endcap review](../unprotected_ports_cancellation_20261025/ports/INDEPENDENT_REVIEW.md)
and new [h3 review](h3/INDEPENDENT_REVIEW.md) and
[h4 review](h4/INDEPENDENT_REVIEW.md) explicitly verify their full
source universes, intermediate blackness, edge XOR and labelled target.

Every case supplies a complete three-step replacement starting at the
original T. No certificate is spliced after part of a different
ledger. In every case the unordered path-edge identity is

    E(Q1) XOR E(Q2) XOR E(Q3) = E(P),

and legality and blackness are separately proved. Three full
complements together with this identity give the asserted target.
This proves the theorem.

Consequently a walk may replace any closed blue endpoint-four move
at n>=15 by a black three-step segment, preserving both endpoints.
This is a local replacement theorem. It does not imply that every
other blue move belongs to this endpoint layer, that all blue moves
are closed, or that the whole black core is connected.

## 4. Genuine small-order hostile and scope boundary

The [order-seven audit](../shape_ports_20261025/tournaments.md), section3,
retains the explicit endpoint-four source and target

    T=(6,4,8,19,103,79,15),
    K_P(T)=(122,124,121,112,32,64,16)=C3 -> 1 -> C3.

Here P=(0,...,6) is directed, closed and blue, with endpoint
outdegrees(2,4). The target has exactly nine Hamiltonian paths and
all nine are blue. The blue/black status of a word is unchanged
across its own K edge, so no black move can enter that target.
It has no positive-length black detour.

This is a real obstruction to replacing n>=15 by an unfiltered
n>=7 statement. It is not a claim that15 is the optimal threshold.
The present synthesis does not settle any remaining orders11..14.
The failed fixed-cut hostiles in the earlier packets likewise remain
failures of their ledgers, rather than global connectivity obstructions.

## 5. Bounded reproduction of the inheritance gate

The parent packets retain their full normal/-O reproduction scripts
and independent reviews. For this synthesis, the following additional
control was run under both normal and optimized Python. It reruns
the24 inherited finite boundary systems using the independent
Gaussian solver, then independently verifies the exact order-seven
source, target, and all5040 target permutations. It is not a census
of large-order tournaments.

```python
import importlib.util, json, sys
from pathlib import Path
from itertools import permutations
root=Path.cwd()
p=root/'05-knowledge/results/endpoint_four_repairs_20261021/audit/independent_audit.py'
sys.path.insert(0,str(p.parent))
spec=importlib.util.spec_from_file_location('endpoint_base_review',p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
rows=[m.boundary(n,b) for n in range(15,27) for b in ('open','closed')]
if any(not r['inconsistent'] and r['minimum_middle_inports']<4 for r in rows):raise RuntimeError('boundary escaped')
S=(6,4,8,19,103,79,15);B=(122,124,121,112,32,64,16);n=7;P=tuple(range(n))
def arc(R,u,v):return (R[u]>>v)&1
def directed(R,q):return all(arc(R,u,v) for u,v in zip(q,q[1:]))
def anti(R,q):
 r=dict(zip(q,q[::-1]));return all(arc(R,u,v)==arc(R,r[v],r[u]) for u in range(n) for v in range(n) if u!=v)
if not (directed(S,P) and anti(S,P) and arc(S,6,0) and S[0].bit_count()==2 and S[6].bit_count()==4):raise RuntimeError('hostile source')
E={frozenset((u,v)) for u,v in zip(P,P[1:])}
B1=tuple(sum((arc(S,u,v)^int(frozenset((u,v)) not in E))<<v for v in range(n) if v!=u) for u in range(n))
if B1!=B:raise RuntimeError('hostile literal target')
HP=[q for q in permutations(range(n)) if directed(B,q)]
if len(HP)!=9 or not all(anti(B,q) for q in HP):raise RuntimeError('pure-blue hostile')
print(json.dumps({'status':'PASS','boundary_systems':len(rows),'consistent_minima':[(r['n'],r['branch'],r['minimum_middle_inports']) for r in rows if not r['inconsistent']],'order7_hostile_target_HP':len(HP),'order7_target_black_HP':0},separators=(',',':')))
```

Both runs return PASS. The consistent boundary minima are open6 at15,17,
open8 at19,21, open10 at23,25; closed6 at17, closed8 at21 and closed10
at25. All other systems are inconsistent. The order-seven target has
nine Hamiltonian paths and zero black ones.

**Audit disposition:** the full closed endpoint-four theorem for n>=15
is an unconditional consequence of the audited source proofs. No
additional bank filter or anchor assumption survives the partition.

