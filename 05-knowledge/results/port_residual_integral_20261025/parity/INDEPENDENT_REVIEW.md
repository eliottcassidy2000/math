# Independent review of the odd-action parity and clone transform

2026-10-25. **ACCEPTED:** equations (P1)--(P9) in [README.md](README.md).
The argument proves the stated modular identities, including their empty
objects and positive-even-order boundary. It does not provide extra
equivariant operations, an integer binomial inflation transform, or black
detours. This review did not modify the producer proof or checker.

## Proof audit

For any finite odd-order group, every orbit of subsets of a fixed size has
odd cardinality. Adding those orbit cardinalities and reducing modulo two
proves (P1) without assuming a faithful action. Substitution of `-1-u`
then gives exactly `u^n`, including `n=0`. The leading coefficient remains
one, so the characteristic-two depth is exactly `n`, not only at least
`n`. The positivity assertion for the integral shifted polynomial uses the
odd-cycle product formula in the inherited positive-jet proof; parity alone
would not prove positivity.

The correction `B` in (P3) is a nonnegative **integer weighted count**.
Its rank-one coefficient for the regular `C5` action is two. A transitive
five-element set cannot map equivariantly onto a two-element set: any
equivariant image is transitive and its size divides five. This excludes
that prescribed quotient, not every imaginable geometric realization of
the coefficient. The producer replaced its original elementary power
calculation with an actual equivariance-map enumeration during this audit;
the revised gate was read and rerun normally and under `-O`. The independent
verification below also enumerates all two-point `C5` actions and all 32
maps. It finds only the trivial two-point action and two constant
equivariant maps, neither surjective.

The constant-XOR descent criterion is also necessary, not only sufficient.
If `gF != F`, the equivalent source subsets `F` and `gF` are sent to the
empty subset and the nonempty subset `F symmetric-difference gF`.
Those outputs cannot be in the same orbit. Thus a fixed toggle descends
iff `F` is invariant, independently of the coefficient field. This is
stronger justification than merely checking that the XOR difference is
nonzero.

The canonical decomposition input is
[THM-4463, homogeneous tournament contacts, section 2](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md),
with the exact weighted-fibre treatment in
[CLONE_INVOLUTION.md](../../ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md).
For a pair-free quotient `R`, positive ordered block sizes give its exact
inflations modulo `Aut(R)`. Inflating cannot merge two distinct quotient
blocks into one transitive module: such a merger would produce a
nontrivial transitive module, hence a homogeneous pair, in `R`. Conversely,
the canonical modules of the inflated tournament recover the quotient and
the sizes. An isomorphism lifts uniquely because a finite transitive
tournament has no nontrivial automorphism. These facts justify the actual
isomorphism classes underlying (P4), rather than just a formal substitution.

Every weight-vector orbit has odd size, so (P4) is valid grade by grade.
The empty quotient contributes only order zero. The substitutions in (P5)
have zero constant term and are well-defined inverse substitutions in
`F2[[x]]`. The all-tournament parity input is independently explained in
[THM-4750, Camion function, section 6](../../../../01-canon/theorems/THM-4750-the-camion-function-half-speed-interpolation-and-the-camion-lift.md):
each reciprocal odd automorphism order is one modulo two. The groupoid
mass has valuation `binom(N,2)-v2(N!)`, which is positive for `N>=3`
(`N=3` separately; for `N>=4`, use `v2(N!)<=N-1`). It equals one at
orders zero, one and two. Hence `A=1+x+x^2`, and substituting
`y=x/(1+x)` gives `y+y^2=x/(1+x^2)`. This proves (P6) at every order.

The integer transform is correctly rejected. At order four, the three
labelled placements of weight two on `C3` constitute one rotation orbit.
The incorrect integer binomial transform gives six instead of four.
This is compatible with, and illustrates, the modulo-two identity.

For (P7)--(P9), a tournament of order `N` with `q` canonical modules has
exactly `N-q` gaps, so its summand reduces to
`(1+u)^q u^(N-q)`. Grouping by the **pair-free canonical quotient** and
using (P4) is therefore legitimate. Odd positive `q` alone survives.
For `N=2m+1`, (P8) is

    (1+u) sum_j binom(m,j) (1+u)^(2j) u^(2m-2j)=1+u.

For `N=2m>0`, it is

    u(1+u) sum_j binom(m-1,j) (1+u)^(2j) u^(2m-2-2j)=u+u^2.

Both use `(1+u)^2+u^2=1` in characteristic two. Orders one and two
are included without a limiting convention. The empty polynomial is one.
No positivity or integer-size information is recovered from this reduction.

## Reproduction and independent literal route

The producer checker passed under normal Python and `python -O`; the two
outputs were byte-identical when reviewed. It enumerates all tournament
classes through order six and checks the signed interpolation through
order twenty. Its order-six result is

    R_6=[26,175,491,738,628,288,56].

The standalone code below imports no producer code. It constructs every
tournament class through order five from literal adjacency matrices,
determines its full automorphism group, and enumerates actual weight-vector
orbits. Thus its signed interpolation uses integer weighted-quotient
objects, rather than the producer's rational odd-partition expansion.
A second path sums the literal cut-subset orbit contributions. It also
tests positive-weight inflation for every pair-free class through order
five and total size through nine, including the order-four integer hostile.
Explicit exception gates remain active under optimization.

```python
from itertools import permutations, product
from math import comb
import json

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def compositions(total, slots):
    if slots == 0:
        if total == 0:
            yield ()
        return
    for x in range(total + 1):
        for rest in compositions(total - x, slots - 1):
            yield (x,) + rest

def models(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    permutations_n = list(permutations(range(n)))
    seen = set()
    result = []
    for bits in product((0, 1), repeat=len(pairs)):
        if bits in seen:
            continue
        T = [[False] * n for _ in range(n)]
        for bit, (i, j) in zip(bits, pairs):
            T[i][j], T[j][i] = bool(bit), not bit
        codes = [(p, tuple(int(T[p[i]][p[j]]) for i, j in pairs))
                 for p in permutations_n]
        seen.update(code for p, code in codes)
        aut = [p for p, code in codes if code == bits]
        gaps = [frozenset((i, j)) for i, j in pairs
                if all(T[i][k] == T[j][k]
                       for k in range(n) if k != i and k != j)]
        require(len(aut) % 2 == 1, 'even tournament automorphisms')
        result.append((T, aut, gaps))
    return result

def weight_orbits(model, excess):
    T, aut, gaps = model
    representatives = {
        min(tuple(w[p[i]] for i in range(len(T))) for p in aut)
        for w in compositions(excess, len(T))
    }
    return len(representatives)

bank = [models(n) for n in range(6)]
classes = [len(row) for row in bank]
cores = [sum(not g for T, a, g in row) for row in bank]
require(classes == [1, 1, 1, 2, 4, 12], 'class census')
require(cores == [1, 1, 0, 1, 2, 5], 'core census')

weight_cases = 0
for q in range(1, 6):
    for model in bank[q]:
        if model[2]:
            continue
        for N in range(q, 10):
            count = weight_orbits(model, N-q)
            require(count % 2 == comb(N-1, q-1) % 2, 'weight parity')
            weight_cases += 1

polynomials = []
signed_counts = []
for N in range(6):
    b = [(-1)**(N-M) * sum(weight_orbits(t, N-M) for t in bank[M])
         for M in range(N+1)]
    signed_counts.append(b)
    R = [sum(b[M]*comb(M, j) for M in range(j, N+1))
         for j in range(N+1)]
    cuts_R = [0]*(N+1)
    for T, aut, gaps in bank[N]:
        g = len(gaps)
        q = N-g
        cut_orbits = set()
        for choose in product((False, True), repeat=g):
            C = frozenset(gaps[i] for i in range(g) if choose[i])
            images = [frozenset(frozenset(p[v] for v in edge) for edge in C)
                      for p in aut]
            key = min(tuple(sorted(tuple(sorted(e)) for e in image))
                      for image in images)
            cut_orbits.add(key)
        for key in cut_orbits:
            k = len(key)
            for j in range(q+k+1):
                cuts_R[j] += (-1)**(g-k)*comb(q+k, j)
    require(R == cuts_R, ('weight/cut fibre mismatch', N))
    require(all(v >= 0 for v in R), 'positive interpolation')
    require(R[0] == cores[N] and R[-1] == classes[N], 'endpoints')
    expected = [0]*(N+1)
    if N == 0:
        expected[0] = 1
    elif N % 2:
        expected[0] = expected[1] = 1
    else:
        expected[1] = expected[2] = 1
    require([v % 2 for v in R] == expected, 'R parity')
    if N:
        inflation_count = sum(weight_orbits(t, N-q)
                              for q in range(1, N+1)
                              for t in bank[q] if not t[2])
        require(inflation_count == classes[N], 'canonical inflation')
    polynomials.append(R)

require(1 + 3*cores[3] + cores[4] == 6 != classes[4], 'integer hostile')
action_count = map_count = surjection_count = 0
for sigma in permutations(range(2)):
    power = list(range(2))
    for _ in range(5):
        power = [sigma[i] for i in power]
    if power != list(range(2)):
        continue
    action_count += 1
    for f in product(range(2), repeat=5):
        if all(f[(i+1) % 5] == sigma[f[i]] for i in range(5)):
            map_count += 1
            surjection_count += len(set(f)) == 2
require((action_count, map_count, surjection_count) == (1, 2, 0),
        'C5 equivariant quotient hostile')
print(json.dumps({'status': 'PASS', 'classes': classes, 'pair_free': cores,
    'weight_parity_cases': weight_cases, 'signed_orbit_counts': signed_counts,
    'R': polynomials, 'C5_two_point_actions': action_count,
    'C5_equivariant_maps': map_count, 'C5_equivariant_surjections': surjection_count},
    sort_keys=True))
```

The independent check passes normally and under `-O`, with byte-identical
output. There are 53 individual pair-free-quotient/total-size weight-parity
checks. Its exact signed weighted-orbit rows are

```text
N=0: [1]
N=1: [0,1]
N=2: [0,-1,1]
N=3: [0,1,-2,2]
N=4: [0,-1,3,-4,4]
N=5: [0,1,-4,8,-12,12].
```

The two independent constructions produce, in increasing powers of `u`,

```text
R_0=[1]
R_1=[1,1]
R_2=[0,1,1]
R_3=[1,3,4,2]
R_4=[2,9,15,12,4]
R_5=[5,29,68,80,48,12].
```

These match the producer’s literal-tournament results. The `C5` hostile
has one two-point action, two equivariant maps and zero equivariant
surjections. Neither this bounded check nor the all-order parity law
claims a full integer recurrence for the class counts.
