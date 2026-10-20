# Every fixed positive black-deficit centre series is a polynomial

**PROVED, all orders.** For every even integer beta>0, let D_beta(z)
enumerate strong self-converse tournament types T satisfying

```text
beta(T)=h(T)-b(T)=beta,
h(T)=H(T)/|Aut(T)|,
b(T)=number of blue Hamiltonian-path orbits.
```

Let L be the least power of two strictly larger than beta+2. Then

```text
deg D_beta <= (beta+2)(2L-1),                              (1)
```

with the zero polynomial understood to have no support. All odd-deficit
centre series vanish. In particular these are finite polynomials, with
an effective quadratic bound on their support; they are not merely
formal inventories whose finiteness is conjectured.

This does not determine the individual coefficients of every D_beta.
The separately completed
[deficit-two theorem](../deficit_two_centers_20261020/README.md)
gives the much sharper exact answer `D_2(z)=z^4`.

## 1. Inheritance and the generalization

The incoming deficit-two proof supplied two competing mechanisms:

* a clique of t half-supported even-order path reroutings forces
  `beta >= 2(t-1)`;
* large automorphic rotation orbits force vertex connectivity.

These are proved in
[the cutoff manuscript, sections 1 and 5](../deficit_two_centers_20261020/session15-deficit-centers/PREFIX_CONNECTIVITY_CUTOFF.md),
with its [approved independent audit](../deficit_two_centers_20261020/session15-deficit-cutoff-audit/AUDIT.md).
The four-vertex exit lemma used there was tailored to beta2. The new
operation is to choose a power-of-two cycle, so **every** nontrivial
rotation has even order, and use a majority of its exit vertices.

The hostile is a cycle whose length has an odd factor: a nontrivial
rotation can then have odd order, invalidating the clique argument.
Another hostile is overlapping a prefix with its positional reflection;
that is why the proof explicitly separates `n<2L` from `n>=2L`.
The retained sidecars are the concrete cycle, its exit vertices, and the
rotation subgroup. A degree inventory alone supplies none of them.

## 2. A general connectivity upper bound

**Lemma.** If T is strong and self-converse, `n>=2L`, and its black deficit
is less than L-2, where L>=4 is a power of two, then

```text
kappa(T) <= L.                                           (2)
```

Suppose instead that vertex connectivity is at least L+1.
Choose an induced strong L-vertex tournament S, and choose any vertex v
outside S. The existence of S follows by taking the induced tournament
on an L-cycle. Here we use Moon's vertex-pancyclic theorem:
each vertex of a strong tournament belongs to a directed cycle of every
length from three through its order.
The exact statement is Theorem 1, page 297 of
[Moon, *On subtournaments of a tournament* (1966)](https://doi.org/10.4153/CMB-1966-038-7);
it was checked in the primary paper during this session.

Either at least L/2 vertices of S beat v, or at least L/2 are beaten by v.
In the latter case apply any anti-automorphism of T to both S and v.
Thus we may arrange at least L/2 exits from S to v.

Since removing L vertices leaves a strong tournament, T-S has a Hamiltonian
path beginning at v, obtained by cutting a Hamiltonian cycle. Fix that
path. For each exit vertex of a fixed Hamiltonian cycle of S, cut the cycle
immediately after the exit and append the fixed path on T-S. This gives
at least L/2 Hamiltonian words of T with a common suffix.

Between two such words the relative vertex permutation is a nonidentity
power of the L-cycle, supported on the first L positions. Its order is
even because L is a power of two. Its support is disjoint from its
positional reflection because `n>=2L`. The inherited clique lemma gives

```text
beta(T) >= 2(L/2-1)=L-2,
```

a contradiction. This proves (2).

For clarity, the clique lemma uses both kinds of loss: two words in the
clique cannot both be blue, and two black words cannot lie either in the
same path orbit or in its converse-paired orbit. Those exclusions follow
because the corresponding automorphism would have even order, with its
two factors supported on disjoint halves. Tournament automorphism groups
have odd order. Thus t-1 black words force at least 2(t-1) black orbits.

## 3. Cycle cuts give the opposite bound

Choose any directed Hamiltonian cycle, with successor permutation R.
Write

```text
m=|Aut(T) intersect <R>|,       q=n/m.
```

The group has odd order, so m is odd. Its q vertex orbits have size m,
and R^q acts transitively on each. Each orbit therefore induces a regular
tournament. The original cycle gives directed perfect matchings around
the cyclic sequence of those q orbits.

The n cyclic cuts fall into exactly q Hamiltonian-path orbits. Among them
the blue cuts give zero or `gcd(2,q)` path orbits; this is the
[cycle-cut coset lemma](../paired_black_leaf_reduction_20261020/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md).
It follows that

```text
beta >= q-gcd(2,q),       q <= beta+2.                    (3)
```

The rotation-orbit cut argument proves

```text
kappa(T) >= ceil(m/2)       when q>=2,
kappa(T) >= ceil(n/3)       when q=1.                     (4)
```

Briefly, a regular m-vertex tournament split as X->Y after s deletions
has `|X|,|Y|<=s`, by averaging outdegrees, hence `m<=3s`.
If q>=2 and fewer than m/2 vertices are removed globally, there can be
at most one split rotation orbit. With no split orbit a surviving
perfect matching forbids a Y-to-X transition around the orbit cycle.
With one split orbit, its two neighboring matchings and the regular
degree inequalities give `m<=2s`. Both contradict the assumed small cut.
The full proof, including q=2 and m=1 boundaries, is in the cited
cutoff manuscript; no stronger bound is applied at q=1.

## 4. The effective polynomial bound

Take L to be the least power of two strictly greater than beta+2.
Then `beta<L-2`. If `n>=2L`, (2) and (4) imply:

* when q>=2, `m<=2L`; oddness sharpens this to `m<=2L-1`.
  Equation (3) gives `n=qm <= (beta+2)(2L-1)`;
* when q=1, `n<=3L`, which is at most `(beta+2)(2L-1)` for beta>=2.

The omitted case `n<2L` already satisfies (1). Thus every centre in
D_beta has the claimed bounded order.

The converse/reversal involution on path orbits pairs all black orbits,
so beta is even. This proves the vanishing for odd beta.

| beta | Chosen L | General support bound | Additional exact information |
|---:|---:|---:|---|
|2|8|60|Incoming classification improves this to `D_2=z^4`|
|4|8|90|Coefficients not classified here beyond the inherited census|
|6|16|248|L=8 would give equality in the obstruction and is invalid|
|8|16|310| |
|10|16|372| |
|12|16|434| |
|14|32|1,008| |

Since `L<=2(beta+2)`, we also have

```text
n < 4(beta+2)^2.                                        (5)
```

Equivalently, positive black deficit in strong self-converse tournaments
grows at least on a square-root scale with order:
`beta > sqrt(n)/2-2`. This is a proved coarse lower bound, not a sharp
asymptotic or an additive SCC law.

## 5. Exact algorithmic coefficients and the power boundary

Formula (1) supplies a terminating determination procedure for every
fixed beta. Enumerate tournament types only through that explicit bound;
retain the strong self-converse ones; enumerate Hamiltonian words and
their positional anti-automorphisms to compute h and b exactly. This
determines every coefficient and the exact degree when the computation
is carried out. No feasible runtime bound is claimed.

Equivalently, for `n` at most (1), the coefficient is the finite sum

```text
[z^n]D_beta(z)
 = (1/n!) sum_[labelled strong self-converse T,
               (H(T)-BluePaths(T))/|Aut(T)|=beta] |Aut(T)|.
```

It is zero beyond (1). The weight converts each labelled isomorphism orbit
to one class; it is not an unweighted labelled count.

Thus the user's powers now have three distinct mathematical meanings:
ordinary convolution powers count SCC placements, the centre sector in
[centers.md](centers.md) has literal square-minus-linear deficits
`q^2-q`, and this proof uses power-of-two cycle lengths to produce
even-order rerouting cliques. Each meaning has a specified operation.
None assumes that black deficits add under SCC concatenation.

The theorem and its boundaries were independently checked by the centre
agent. The new argument is symbolic; the complete finite residual replay
in [incoming_replay.json](incoming_replay.json) audits the stronger
deficit-two result, not every higher-deficit coefficient.
