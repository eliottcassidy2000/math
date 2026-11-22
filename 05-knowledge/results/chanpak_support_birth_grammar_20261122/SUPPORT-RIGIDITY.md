# The Hamiltonian arc support remembers the whole tournament

**Status: PROVED.** The new reconstruction uses a structural theorem
of Busch–Jacobson–Reid on dead arcs. It yields a grammar over partial
orientations, not just enumeration of tournament completions, and
a parity criterion that recovers all stripped arcs. It does not give
a target-size bound or prove the Chan–Pak spectrum conjecture.

## 1. The stripped carrier

Let S(T) be the spanning oriented digraph obtained by retaining only
arcs traversed by some directed Hamiltonian path of a tournament T.
Call a deleted arc *dead*. The directed Hamiltonian paths of S(T)
are exactly those of T, directly from the definition.

If T has ordered strong components C_1⇒…⇒C_q, its Hamiltonian paths
concatenate one directed Hamiltonian path in each component. In a
nonsingleton strong component, cutting a directed spanning cycle
gives all starting and ending vertex choices. It also shows its
internal support is strongly connected. Consequently S(T) consists
of the internal supports S(C_i), all forward arcs between adjacent
components, and no arcs between nonadjacent components. In particular
the SCCs of S(T) are exactly those of T, and every deleted cross-arc
points forward in their recovered order.

## 2. A dead arc has a mod-four reverse-order count

For two distinct labels a,b, let h_(a<b)(D) count directed Hamiltonian
paths of D in which a precedes b, not necessarily consecutively.

**Theorem.** If a→b is dead in a tournament T, then
h_(b<a)(S(T)) is divisible by four. In particular a precedes b in an
odd number of directed Hamiltonian paths of S(T). The latter criterion
orients all missing unordered pairs of S(T) uniquely.

**Proof.** For a cross-component arc, h_(b<a)=0. If a,b belong
to the same strong component, the other components multiply the
count by their ordinary Hamiltonian counts, so first reduce to T
strong.

Use the structural separator theorem for a dead arc in a strong
tournament: there is a vertex z distinct from a,b such that T−z has
ordered SCCs C_1⇒…⇒C_k, with k≥4, a in C_1, and b in C_k.
See [Busch's dissertation, Theorem 4.12](https://digital.auraria.edu/downloads/ma7yx-pcr25/AA00005831_00001.pdf)
and [Busch, Jacobson and Reid, On arc-traceable tournaments](https://doi.org/10.1002/jgt.20173).
Only the stated separator properties are used; the paper imposes
additional restrictions on arcs incident to z.

Any Hamiltonian path with b before a must split at z. Its part
before z follows component order and includes b, so ends in C_k.
Its part after z includes a, so begins in C_1. The split is forced
into b on the prefix side and a on the suffix side. For each
middle component C_i, partition its vertices arbitrarily between the
prefix and suffix and take one Hamiltonian path on each assigned
subset, using the empty path when appropriate. Each middle component
therefore supplies the independent factor

```text
Z_C(2)=Σ_(A⊆C) H(C[A]) H(C−A).
```

Endpoint factors account for the path in C_1 after z starting at an
outneighbor of z and containing a, and the path in C_k before z
ending at an inneighbor of z and containing b. Each endpoint factor
includes the unconstrained complementary path within its own
component. Conversely, every such data set concatenates to a unique
path with b before a: all cross-component moves in either part are
forward, and the only checks across z are exactly these endpoint
conditions. Thus the ordinary reverse-order count is divisible by the
product of Z_C(2) over middle components. Each Z_C(2) is even, since
modulo two every term is one by Rédei (and H(empty)=1).
There are at least two middle components. The opposite-order count
is therefore divisible by four. As H(T) is odd, the observed
a-before-b count is odd. ∎

In a strong component the reverse-order count is positive even for a
dead arc, by cutting a directed spanning cycle at b. Hence the
middle components in a separator as above have at most log_2 H(T)
vertices in total: Z_C(2)≥2^|C|. This illustrates why "dead" means
unseen as an adjacency, not that its direction agrees with
**every** path order. The parity assertion is what recovers it.

## 3. Canonical completion without a missing-edge search

For a finite spanning oriented graph D with every arc on some
Hamiltonian path, compute H(D). If H(D) is even, reject it as a
candidate support. If it is odd, complete every missing pair
{a,b} by orienting it in the direction of odd relative order
among the Hamiltonian paths of D; retain already present arcs.
Call the resulting tournament T*(D).

**Recognition rule.** D is a tournament support if and only if
H(D) is odd and H(T*(D))=H(D). Necessity follows from rigidity.
For sufficiency, all paths of D survive in its completion, and
total-count equality says there are no additional paths. Every
arc of D was on a surviving path, so its support is exactly D.

In particular a support has **one** completion with the same path
family, not a cube of free missing arcs. The conclusion is on labelled
objects; all constructions respect relabeling. Thus two tournaments
are isomorphic iff their support digraphs are, and
Aut(S(T))=Aut(T). The same completion commutes with taking the
converse because every directed Hamiltonian path reverses as a word.

A separate direct seam fact needs no separator theorem. Holding all
arcs except between a,b fixed, at least one orientation of that pair is
on a Hamiltonian path of its completion. Take a Hamiltonian path of the
remaining vertices. If its first vertex does not dominate both a,b,
prepend a pair whose outgoing endpoint dominates that vertex. Otherwise
insert the pair just before the first later vertex not dominating both:
its predecessor dominates both and the exiting endpoint can be
chosen. If all remaining vertices dominate both, append either pair
orientation. With no remaining vertices either orientation works.
It follows that reversing a dead arc preserves the old Hamiltonian paths
and adds an even number at least two. The [primitive seam note](PRIMITIVE-SEAMS.md#a-primitive-strong-seam-changes-an-arc-active-at-both-ends) strengthens this for strong tournaments: reversing a dead arc adds at least four paths, so a primitive strong unit seam cannot use it.

## 4. A recursive exact grammar and generating-function coefficient

Fix n and a target odd h. Fix the rooted spine word (1,…,n) and start
with its oriented path graph, not a tournament. At each step choose
another permutation word not already a directed Hamiltonian path of
the current graph, with its consecutive arcs compatible with all
arcs present. Adjoin its consecutive arcs. Discard a state if its
Hamiltonian count exceeds h, and deduplicate equal labelled arc graphs.
Each state has every arc on a Hamiltonian path: the paths that
introduced the arcs survive each addition. Adding a genuinely new
word raises the path count by at least one. The grammar therefore
never needs more than h−1 additions after the root and is finite for
each n,h.

At H(D)=h apply the canonical completion test; if rejected, stop
that branch anyway, since any further addition would exceed h.
Every standardized support with h paths is generated: start at its
root spine and add any of its paths still missing from the current
path family. Its final arc union is S(T) by definition, and the
intermediate path count never exceeds h because all interim arcs
belong to S(T). Thus one can generate and filter **non-tournament
partial orientations**, then recover the hidden arcs once, rather
than enumerate all completions.

Let G_(n,h) be the accepted labelled support objects with fixed root
word. Then the exact generating function by order is

```text
F_h(z) = Σ_(n≥1) z^n / h · Σ_(D∈G_(n,h)) |Aut(D)|.
```

For an unlabelled tournament T with h paths, standardizing one
chosen path to the root word gives h/|Aut(T)| distinct labelled
objects. This is because the automorphism action on the full
ordered paths is free. Use reconstruction to replace Aut(T) by
Aut(D); the weights for the objects representing that class sum
to one. The same grammar with strongly connected accepted
supports computes the strong-component series A_d(z).

In terms of those strong series the familiar ordinal-component
generating rule becomes, for h>1,

```text
F_h(z) = Σ_(k≥1) (1−z)^(-k−1)
             Σ_(d_1···d_k=h, d_i>1) A_(d_1)(z)···A_(d_k)(z).
```

Every Hamiltonian path has one piece in each strong component;
the h-factors multiply. Singleton runs can be inserted in the
k+1 gaps, explaining the denominator. The all-transitive h=1
case gives z/(1−z) when the empty tournament is omitted.

Each A_d is a polynomial: a nontrivial strong tournament has a
directed spanning cycle and hence at least as many Hamiltonian
paths as vertices, by its distinct cuts. It is therefore enough
to execute the strong-support grammar on orders at most d.
This proves the entire rational function F_h is effectively
computable, independently of any bound on a first nonstrong
realization. Busch's sharper strong minimum can of course reduce
the finite search. For attained h>1, let K be the largest length
of an ordered decomposition into realized strong factors.
Positivity of the strong coefficients shows the pole of F_h at
z=1 has exactly order K+1. The large-order class counts are
eventually polynomial in the order n, of degree K with leading
coefficient

```text
(1/K!) · Σ_(d_1···d_K=h) A_(d_1)(1)···A_(d_K)(1).
```

Thus birth depends on a *minimum total size*, whereas eventual
growth depends on a *maximum number* of strong factors; these
optimize different decompositions.

This separates **two sources of multiplicity at birth**. The first
order for h is the minimum of the summed first strong orders of
its ordered multiplicative decompositions; singleton padding never
helps. The leading coefficient is the sum of products of the first
strong multiplicities over the minimizing lists. Therefore a
minimal representative with at least two nonisomorphic strong
components immediately gives multiple smallest classes by
permuting the components. A unique decomposable class at birth
must consist of identical copies of one minimal strong class,
and have no competing minimizing list. For prime h, there is
just one possible strong factor, so uniqueness or multiplicity
is entirely the count of smallest accepted strong supports.
This is an exact criterion rather than a classification of
which odd h have particular minimizing lists.

### A small birth contrast with an actual stripped arc

The complete source lists at orders three through five give a unique
minimal class for h=5 (order four, outgoing masks (2,8,3,5)); its
support retains all six arcs. At h=9 the first order is five and
there are two classes, with rows (2,4,17,7,11) and (2,16,3,7,13).
In the first, the 4→1 arc is the only dead arc, and the opposite
relative order occurs in exactly four paths, while the second is
arc-traceable. Thus the two birth mechanisms are already distinct
on the stripped non-tournament carrier: nine paths supported on
nine edges versus nine paths supported on ten edges. This is a
finite explanatory example, not a general claim that the number
of active edges classifies birth multiplicities.

## Scope

The separator/rigidity theorem controls the invisible arcs after
a support has been found. This is not a proof that a support with
prescribed h exists at an efficiently bounded order, and does not
provide a uniform physical straddler. It suggests a sharper seam
search: generate two target support languages, recover their
unique completions, and check one-arc distance rather than assume
that a shared partial pattern leaves free orientations.


## Reproduction

The C++17 companion checker reads the complete stored representative
decks at orders three through eight. For each literal tournament
it independently enumerates all spanning directed paths, strips
the arcs, re-enumerates paths from that graph, checks every missing
pair's opposite-order divisibility by four, and reconstructs the
original rows exactly. Coverage of the discrete decks is finite
control, not the source of the all-order theorem. It tests order by
label masks, not a sorted numerical tournament profile.
