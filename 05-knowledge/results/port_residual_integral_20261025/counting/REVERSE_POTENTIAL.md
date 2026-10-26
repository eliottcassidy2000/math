# A reverse valuation potential for normalized binary digits

2026-10-25. **PROVED algebraic construction:** a reverse minimum-cost
pass removes a common binary baseline from the ordinary scalar transfer
without requesting additional residue digits. The ordinary state graph
is unchanged. Integrating this pass into the production ordinary engine
and measuring its cost remain **OPEN implementation work**. The
independent proof review is [REVERSE_REVIEW.md](REVERSE_REVIEW.md).

## The finite acyclic construction

Take the structurally reachable and terminal-reaching graph of the
ordinary recurrence in [README.md](README.md). A vertex includes stage
ell and the exact (S,D) key. Its integer transfer along an edge e:v->w is

    Y_w += T_e*2^J_e*Y_v,        T_e positive odd, J_e>=0.

Here excess is zero. Thus T_e is the odd-factorial ratio in(I5), and
J_e is(I7) with r=S'. One initial vertex has value1. Identify terminal
vertices by zero-cost unit edges if necessary. The stage decreases at
every edge; no termination assumption about an infinite orbit is used.

Define a suffix valuation by a reverse pass:

    V(terminal)=0,
    V(v)=min_(e:v->w) (J_e+V(w)),
    B=V(initial).

Every edge has reduced cost

    J'_e=J_e+V(w)-V(v)>=0.                         (R1)

Use the coordinate change Z_v=2^[V(v)-B]*Y_v. Its formula can contain a
negative exponent; nevertheless it is integral on every reachable
partial sum. Indeed a prefix path p to v has cost c(p)=sum_(e in p)J_e.
Append a minimizing suffix. The resulting complete path has cost
c(p)+V(v)>=B. Its prefix contribution to Z_v is therefore an odd integer
times the nonnegative power 2^[c(p)+V(v)-B]. This is pathwise, before
any merger. Because V is a function of the destination state, addition
commutes with the change of coordinates.

The transformed recurrence is consequently

    Z_initial=1,
    Z_w += T_e*2^J'_e*Z_v.                         (R2)

At the terminal it gives Y_terminal/2^B. It can be run modulo2^K from
the start; every shift is nonnegative and every remaining division is
by an odd unit. No already truncated Y_v is ever divided by a power
of two. This theorem is about a reachable lattice, rather than the full
product lattice of arbitrary integer Y_v.

## Identification of B for the ordinary tournament recurrence

A complete path selects one odd cycle partition lambda of N. The
intermediate factors O(u)*2^[S*(m-b)] in the preceding gauge telescope.
The terminal contribution is O(N) times the original rational Burnside
summand. The product of odd transfer factors is

    product_e T_e=O(N)/product_d[d^k_d O(k_d)],

which is odd. Hence the sum of edge costs is exactly

    e(lambda)-sum_d v2(k_d!),

the valuation of that summand. The
[termwise baseline theorem](../precision/README.md) therefore identifies

    B=floor(N/2), N>=3;       B=0, N=0,1,2.         (R3)

Equality is attained by the single N-cycle at odd N and by the pair
(1,N-1) at even N>=4. B is the minimum summand valuation, not necessarily
the valuation of the final sum. For instance N=50 has an extra factor
of two in the sum.

Finally Y_terminal=O(N)*a(N), so(R2) gives

    a(N)/2^B modulo2^K

after multiplying by the modular inverse of O(N). This is a third
precise observation route, alongside the unchanged absolute-residue
engine (which needs B+K digits for that question) and the bounded-cycle
partition filter. It needs the extra reverse graph pass; no smaller
state graph or complexity advantage has yet been demonstrated.

## Hostile boundaries and finite control

The graph must be generated from structural reachability, not from
nonzero absolute residues. A branch invisible modulo2^K may affect the
normalized answer. Stage must remain part of a vertex. A single
subtraction of B at the initial edge is also invalid: a zero-cost skip
can precede a costly continuation. The state-dependent V makes the
cost inequality hold on every edge.

The sign-free ordinary case has positive arithmetic transfer factors.
That is not a construction of individual unlabelled tournaments in each
state. The signed clone recurrence requires its own graph and valuation
bookkeeping; no positive clone model is inferred.

The following standalone control tests the generic finite-DAG mechanism
on arbitrary nonnegative costs and positive odd factors. It checks every
intermediate scaled value, every reduced edge cost, terminal equality
and K=1..6 residues. It does not claim to be an implementation of the
ordinary tournament state graph.

```python
from fractions import Fraction
from random import Random
import json
rng=Random(2026102501)
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
cases=0;coordinates=0;modular=0;negative_naive=0
for size in range(2,19):
 for trial in range(12):
  E=[[] for _ in range(size)]
  for i in range(size-1):
   for j in range(i+1,size):
    if j==i+1 or rng.randrange(4)==0:
     E[i].append((j,rng.randrange(9),2*rng.randrange(6)+1))
  V=[0]*size
  for i in range(size-2,-1,-1):V[i]=min(J+V[j] for j,J,T in E[i])
  B=V[0];Y=[0]*size;Z=[0]*size;Y[0]=Z[0]=1
  for i in range(size):
   scale=Fraction(2)**(V[i]-B)
   need(scale*Y[i]==Z[i] and Z[i]>=0,'integral prefix identity')
   coordinates+=1
   for j,J,T in E[i]:
    Jp=J+V[j]-V[i];need(Jp>=0,'reduced cost')
    Y[j]+=T*2**J*Y[i];Z[j]+=T*2**Jp*Z[i]
  need(Y[-1]==2**B*Z[-1],'terminal normalized value')
  negative_naive+=any(J<B for j,J,T in E[0])
  for K in range(1,7):
   mod=1<<K;W=[0]*size;W[0]=1
   for i in range(size):
    for j,J,T in E[i]:W[j]=(W[j]+T*pow(2,J+V[j]-V[i],mod)*W[i])%mod
   need(W[-1]==Z[-1]%mod,'fixed normalized precision');modular+=1
  cases+=1
need(negative_naive>0,'hostile naive initial subtraction exercised')
print(json.dumps({'status':'PASS','acyclic_graphs':cases,'prefix_coordinates':coordinates,
 'modular_comparisons':modular,'naive_initial_shift_failures':negative_naive},sort_keys=True))
```
