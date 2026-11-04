# Source-relative module separation after dyadic transitive thickening

**Status: integration proof with exact literal strong-tree controls.**
Every fixed nonempty tournament has an explicit uniform transitive inflation
whose source response differs from every distinct target obtained by XOR of
modules of that inflated source. The original tournament is not thereby
reconstructed, and arbitrary pairs of targets are not compared.

## 1. Exact source, action, and constructive size

Let T be a finite nonempty named tournament. Let P(T) be the finite collection
of prime quotients in its canonical strong-module tree. A strong module is a
module overlapping no other module, not necessarily a strongly connected
tournament. The two-child/ordered quotient convention is linear, not prime.
In particular every quotient in P(T) has order at least three.

For a tournament quotient Q of order m>=2, the separately proved Seidel
energy theorem defines the positive rational 2-adic integer

    E(Q)=2^(6-m) sum_(nonempty even I) ||S_I^(-1)1||^2,

and the explicit exponent L_Q=1+v2(E(Q)). Set

    L_T = max( {1} union {L_Q : Q in P(T)} ).

If P(T) is empty, this gives L_T=1. Each energy is computed by finitely many
exact rational linear solves; if E=p/q is reduced, q is odd and L_Q is one
plus the number of factors of two in p. No implicit threshold is used.

Take any positive integer N divisible by 2^L_T. Replace each old vertex v
by the transitive tournament TT_N on names (v,1),...,(v,N), retaining v's
incidences between cells. Denote the resulting named source by

    T[N] = T[TT_N,...,TT_N].

Let F be any family of modules of this ORIGINAL inflated tournament T[N].
Toggle each old physical edge if an odd number of members of F contains
both endpoints; call the resulting tournament U. The family can cross or
nest, and can include new modules inside the transitive cells. No assumption
that its selected modules remain modules in a sequential ordering is made.

**Theorem.** If U differs from T[N], then their complete named two-port
Hamilton responses differ. This holds for every N satisfying the displayed
divisibility condition, not merely along an unspecified subsequence.

The case |T|=1 is allowed: its inflated source is transitive and the proof
uses only the linear-node theorem. It is not an application of a nonexistent
order-one Seidel noncancellation theorem.

## 2. Reviewed inputs

The new input is `seidel_energy/PROOF.md`:

1. If N is a positive multiple of2^L, TT_N has complete profile
   g_h=(-2)^(h-1) modulo2^L.
2. The normalized g-profile class is closed under arbitrary quotient
   substitution. In particular, for every nonempty tournament A,

       q_h(A[TT_N]) = g_h (mod2^L),
       H(A[TT_N]) = 1 (mod2^L).

3. Every quotient Q of order at least two, with normalized-g-profile
   children modulo2^L_Q and corresponding profile-preserving target
   children, has a child-uniform quotient-converse separator in every fixed
   ambient. Its complete augmented contrast is congruent to
   (delta,2delta,0,...) modulo 2^L_Q, and the ambient multiplier is odd.

The integration also uses the following separately reviewed Session188
proofs at commit `8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7`, under
`05-knowledge/results/session188_structural_bridges_20261104/`:

- `module_decomposition/PROOF.md`: the original strong tree remains a common
  recursive modular partition under original-module XOR; a prime quotient
  can only remain unchanged or be conversed; scalar Hamilton equality
  implies equality of every complete ordered-cover profile on every fixed
  physical induced subset.
- `highest_quotient_reduction/PROOF.md`: highest changed quotient nodes are
  disjoint, contain all changed edges, and any highest changed linear node
  is detectable. Equal-profile other highest nodes can be restored without
  changing the selected node's actual lifted query.
- `covered_laminar/ORDERED_BLOCKS.md`: the underlying fixed-ambient ordered
  block separator, permitting equal-profile child replacements.

The classical prime/linear strong-module decomposition is attributed to
Gallai, as in the inherited structural packet. The next lemma spells out the
substitution step used here rather than assuming new prime quotients away.

## 3. Uniform transitive inflation introduces no prime quotient type

**Lemma.** Every prime quotient in the canonical strong-module tree of T[N]
is the same quotient, on correspondingly inflated children, as a prime
quotient in T's strong-module tree. If A is such an old child, its new
physical child set carries exactly A[TT_N]. Linear quotients can expand or
flatten; no new prime quotient is created inside a transitive cell.

**Proof.** Induct using the prime/linear decomposition.

At a prime node Q[A_1,...,A_m], its inflated realization is
Q[A_1[TT_N],...,A_m[TT_N]]. In a substitution with prime Q of order at least
three, these displayed children remain the maximal proper strong modules.
To see the relevant obstruction, suppose a module M overlaps a child X_i.
It meets at least two children and misses part of X_i. Its support J in Q
is a quotient module: any vertex in a child outside J sees the nonempty
intersections with all children of J uniformly. Primality forces J to be
all quotient names. But a vertex of X_i outside M would then have to see all
other child sets in one direction, making i a source or sink in Q. A prime
quotient of order at least three is strongly connected, so this is impossible.
No proper module can unite several full children either, by primality.
Therefore the prime quotient and its physical children persist exactly.

At a linear node, the children are its ordered strongly connected
components, with singleton components permitted. A nontrivial strongly
connected tournament remains strongly connected after substituting nonempty
cells: paths in its quotient connect any two cells, and a directed quotient
cycle permits travel between vertices in the same cell. An old singleton
component becomes TT_N, whose strongly connected components are its N
singleton vertices. Thus the linear quotient only expands by these ordered
singletons. Recurse within the other components. The base case TT_N is
transitive and has no prime quotient. This proves the lemma. QED.

Consequently every prime node of the inflated source has quotient Q in P(T),
and every one of its children has normalized profile g modulo2^L_T, hence
also modulo2^L_Q. The canonical tree need not remain the target's canonical
tree after the module operations; only the original-source tree is used.

## 4. Highest-node proof of the theorem

Suppose, to the contrary, that U differs from T[N] but their complete named
responses agree. The actual query a=c=0 has the ports preceding every old
name, so its Hamilton count is H of the old tournament. Hence H(T[N])=H(U).
Hereditary profile equality gives, on every fixed physical subset A,

    q_h(T[N][A]) = q_h(U[A])   for every h.

There is at least one highest changed quotient in the original strong tree
of T[N]. If any such quotient is linear, the inherited ordered-block
separator and equal-profile restoration contradict response equality.
Otherwise choose one highest changed prime node M. By Section3 its quotient
is an original Q in P(T), and each source child lies in the normalized
geometric profile class modulo2^L_Q. The target quotient is Q^op and all
corresponding child profiles agree by hereditary equality.

Form V by changing just M's source interior to its U interior, leaving all
exterior names and arcs as in T[N]. The Seidel theorem gives an actual
child-uniform query separating T[N] from V, with both ports imitating M at
every exterior vertex. This is the complete ambient-safe query promised by
the all-cover identity, not a generic lift of a scalar gap.

For every other highest changed node K, its source and target profiles agree
and it is disjoint from M. Both ports see K uniformly because they inherit
M's exterior incidences. Thus K remains a module in this same actual port
extension. Replacing its source interior by its target interior preserves
the queried Hamilton count exactly. Performing all these disjoint
replacements turns V into U, since all changed edges lie in the highest
nodes. The separator between T[N] and V is therefore a separator between
T[N] and U, the desired contradiction. QED.

## 5. Boundaries and checks

This is a theorem about the newly inflated source's ORIGINAL-module XOR
family. It does not show that the uninflated T has this property, that
arbitrary same-size tournaments are reconstructed, or that every pair of
family outputs is distinguished. It does not authorize arbitrary sequences
using newly created modules. The blowup can be very large, and the displayed
exponent is sufficient rather than minimal or polynomially computable.

The source-profile closure and exact noncancellation are supplied by the
separately audited all-order Seidel proof. The independently checkable
structural controls here enumerate modules and strong nodes directly on
small literal inflations, verifying that every prime node/child partition
comes from the original source. These finite controls test the substitution
lemma, not an all-size inference from small cases or an exhaustive search
through the potentially enormous module-family space.

## 6. Reproduction and finite-data scope

Run `python3 verify.py > normal.json` and `python3 -O verify.py > optimized.json`,
then compare the receipts. The verifier imports no Seidel producer code. It
enumerates every module and original strong node directly in 155 literal
inflations, including all 75 labeled base tournaments through four vertices
with TT2 and TT3 cells and five larger named controls. It tests 635,370
nonempty subsets, retains all physical cell sets and module lists, and
compares 90 prime-node/child partitions exactly. The largest literal source
has 18 vertices. `structural_controls.json` retains the complete data.

These N=2 and N=3 structural controls are not claimed to satisfy the larger
energy-based sufficient exponents. The modulus and actual normalized-profile
claims are the all-order Seidel theorems, independently reviewed in that
packet. The finite checks here test the strong-tree substitution step only.
Frozen files are PROOF.md, README.md, verify.py, structural_controls.json,
normal.json and optimized.json.
