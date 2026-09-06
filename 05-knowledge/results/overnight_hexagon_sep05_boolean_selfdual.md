# Parity trace, self-complementary two-graphs, and the recovered complement count

Status: **PROVED**, independently audited by root and nc2_seed;
small complete dictionaries are **FINITE-EXACT**. This is a trace and
parity-basis transfer, not a Boolean spectrum transfer. No new ID.

## 1. Inheritance: two correct maps with different scopes

The closest proved mechanism is the full, all-order Fourier orbit basis in
[THM-4078 — signed dual spectrum](../../01-canon/theorems/THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation.md),
read including its historical boundary and current continuation. The current
[Peck/parity note](overnight_hexagon_sep05_boolean_peck.md) supplies the exact
native parity index. The canonical hostile is already at n=4: a self-dual
Fourier mode is a weighted zero mode but not a Boolean zero mode.

The corrected near miss is the representative-level identification:
[THM-1430 — graph switching and E_n](../../01-canon/theorems/THM-1430-graph-switching-is-exactly-E-n.md)
explicitly distinguishes its odd-order unique Eulerian representative from
the even-order count identity. The bicycle explanation is in
[THM-1440, Section C](../../01-canon/theorems/THM-1440-seidel-spectra-are-sine-and-the-odd-n-projection-is-one-parity-fact.md).
The least-used sidecar is the involution on the **dual orbit set**, rather
than only its scalar eigenvalues. The live board is: parity trace;
Fourier orbit sums; complementary orbit pairs; odd-order projection;
complementing permutations; actual Boolean kernel.

An older independent strand is
[THM-516, Section 4 — Burnside enumerator family](../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md)
and the executable
[self-complementary graph engine](../../04-computation/self_complementary_graphs_burnside_s566.py).
Its live table/output correct the old n=13 header near miss to 5600. We use
the exact complement-orbit mechanism below, not an unverified sequence match.

## 2. The all-order trace dictionary

Let Z be the binary Eulerian cycle space on n vertices, W its dual cut
quotient, and G=S_n. Write K for the full edge set and

```text
eta(F)=(-1)^|F|,             kappa([H])=[H+K] on W.
```

Relabelling fixes K, so kappa is an involution of W/G. Let q_+,q_- count
the even- and odd-edge classes of Z/G, and Delta=q_+-q_-. On the invariant
function space C[Z]^G, multiplication by eta has trace Delta in the class-
indicator basis. In THM-4078's Fourier orbit basis,

```text
Psi_Omega(F)=sum_(h in Omega) (-1)^<h,F>,
eta Psi_Omega = Psi_(kappa Omega).                   (1)
```

Thus its trace is the number of fixed dual orbits. For every n>=3,

```text
Delta_n = #{Omega in W/S_n : kappa Omega=Omega}
        = #self-complementary two-graph isomorphism classes.  (2)
```

Here a two-graph is the set of odd triangles of a signing; adding K
complements that entire triangle set. “Self-complementary” in (2) allows
vertex switching followed by relabelling, not just ordinary graph
complementation. No canonical Eulerian representative is required for (2).

There is a fuller dictionary, not just a trace. Every fixed Omega gives
an even-supported function Psi_Omega. Every free complementary pair
`{Omega,kappa Omega}` gives

```text
Psi_Omega+Psi_(kappa Omega), supported on even-edge Eulerian classes;
Psi_Omega-Psi_(kappa Omega), supported on odd-edge Eulerian classes.
```

These are bases of the respective parity subspaces, since they diagonalize
the permutation involution on the complete Fourier basis. Consequently

```text
q_+ = #(W / <S_n,kappa>),
q_- = #unordered pairs of distinct complementary S_n-orbits in W. (3)
```

This is an exact change of parity bases. It does not identify the native
Boolean adjacency with a diagonal operator in those bases.

## 3. Exactly which zero-mode statement follows

For odd cycle length j, adding K reverses every cycle sign. If Omega is
self-complementary, its negative-j-cycle count therefore equals half the
total. THM-4078 gives

```text
M_j Psi_Omega=0 for EVERY odd j,                     (4)
```

where M_j is the multiplicity-weighted cycle operator. These Delta_n
functions are independent joint weighted zero modes. A paired, non-fixed
orbit with half its triangles negative could add weighted zero modes;
(2) does not assert weighted kernel dimension exactly Delta_n.

Booleanization does not preserve (4). At n=4 the unique self-complementary
switching orbit is the single-negative-edge orbit. On the Eulerian classes
empty,C3,C4 its Fourier sum is

```text
Psi=(6,0,-2),
M_3 Psi=0,
B_3 Psi=(0,4,0) !=0.                               (5)
```

The Boolean kernel is instead spanned by `(1,0,-1)`. In addition this
particular labelled cut class has **no Eulerian representative**. This
separately blocks both an unweighted eigenvector transfer and an all-order
representative map. At n=3 there is no self-complementary dual orbit.

The native Boolean nullity bound `nullity(B_3)>=Delta_n` comes from its
rectangular bipartite block, as proved in the Peck/parity note. Equation
(2) interprets that dimension deficit, but does not supply a natural map
from the weighted vectors (4) to native Boolean kernel vectors. Such a map
would need to retain the flattened edge multiplicities or solve a new
native incidence relation. A Laplacian-gap conclusion is still absent.

## 4. Odd-order representatives and their even-order hostile

For odd n, let d(H) be the set of odd-degree vertices of H. It has even
cardinality. The graph

```text
theta([H])=H+cut(d(H))                              (6)
```

is Eulerian: an even-sized cut on odd n has degree-parity vector exactly
its vertex indicator. Conversely no nonzero cut is Eulerian at odd n,
so (6) is the unique representative. It is equivariant under relabelling.
Since K is Eulerian when n is odd, theta also commutes with complement.
Thus (2) then specializes to

```text
Delta_n = #self-complementary EULERIAN graph classes on n vertices. (7)
```

The adjective matters. At n=5 there are two ordinary self-complementary
graph classes but only one Eulerian such class, namely C5. The Fourier
mode corresponding to that class is not its class indicator: (1) still
puts the mode on the even-edge side, although C5 itself has five edges.

At even n, a cut changes every vertex degree parity by the same value
`|S| mod2`. A class whose degree vector is neither zero nor all-ones has
no Eulerian representative; otherwise it has `2^(n-2)` of them. This is
the retained boundary in THM-1430, independently replayed below. It is not
a reason to restrict the Fourier trace identity (2), which remains valid.

## 5. Recovering an older count one order lower

Let SC(m) count ordinary self-complementary undirected simple graph
isomorphism classes on m vertices. Complement commutes with relabelling,
so the same invariant-trace argument gives its familiar twisted Burnside
formula: a permutation contributes zero unless all of its unordered-edge
orbits have even length, and then contributes two choices per edge orbit.

At m=4k, the admissible vertex-cycle types are exactly `4mu`, mu a
partition of k. Odd vertex cycles of length>=3 have odd internal edge
orbits; a cycle of length2 mod4 has an odd antipodal matching orbit. At
most one fixed vertex is possible, and the total multiple of four rules
it out once all other lengths are divisible by four. Conversely all
unordered-edge orbits have even length when every vertex cycle has length
divisible by four.

If mu has r parts and `S=sum_(i<j)gcd(mu_i,mu_j)`, the edge-orbit count
of 4mu is `2k+4S` and `z(4mu)=4^r z(mu)`. Hence

```text
SC(4k)=sum_(mu partitions k) 2^(2k-2r+4S)/z(mu).     (8)
```

This is exactly the weight found by the new Eulerian parity criterion
at cycle type `(1,4mu)`. Therefore, for every k>=1,

```text
Delta_(4k+1)
 = #self-complementary two-graphs on 4k+1 vertices
 = #self-complementary Eulerian graphs on 4k+1 vertices
 = SC(4k)
 = #ordered pairs of loopless digraphs on k vertices, up to joint relabelling. (9)
```

The last equality was derived by ordered-arc orbit counting in the
Peck/parity note. Values 1,10,720,703760,9168331776 are therefore a
recovered self-complementary sequence, not an unidentified new sequence.
The cycle-type map preserves the exact fixed weight. It is not an explicit
objectwise bijection between the different sets in (9), except for the
odd-order representative map between its first two graph interpretations.

## 6. Exact controls and scope

The [checker](../../04-computation/overnight_hexagon_sep05_boolean_selfdual.py)
materializes only n=3,...,6. It builds switching orbits by adjacent
transpositions and independently tests each complement-fixed status by
literal odd-triangle sets under every vertex permutation. It reconstructs
the entire orbit Fourier matrix, verifies (1) entrywise, tests every fixed
mode against both M_3 and B_3, and counts the Eulerian members of each
cut class. Every Fourier matrix is checked invertible exactly.

The recovered SC(4k) identity is separately checked for k=1,...,5 against
the older literal edge-orbit complement engine, without enumerating graph
states at those larger orders. Ordinary SC(5)=2 is checked by a separate
literal graph census. These finite controls test the all-order proofs;
they are not their justification.

```bash
python3 -B 04-computation/overnight_hexagon_sep05_boolean_selfdual.py
python3 -B -O 04-computation/overnight_hexagon_sep05_boolean_selfdual.py
```

The companion [output](overnight_hexagon_sep05_boolean_selfdual.out) retains
the complete small counts and the minimal nontransfer witness. Constructive
native Boolean kernels, exact all-order nullity, and the Boolean Laplacian
gap remain **OPEN**. The scalar/canonical-representative/Fourier/kernel
maps above have deliberately not been conflated.

Root and nc2_seed independently read and passed the complete proof,
including the odd/even representative boundary and the recovered
complement-twisted Burnside weights. Normal and optimized executions
agree on all 398 optimization-live gates. Frozen checker SHA256:
`d3d9d6ec6520dba9fc9bb6a1b916f1780bc8dd2bed7cd4266a11b52e2d1def48`;
output SHA256:
`cfcb4194b3b34d11c2d9df9323431d42649a0fab89a7dced701f4a331bad555c`;
semantic trace:
`18f878ce55bd8686ee112df72dd45f5fefd678a09e2138cee22720c2b367a317`.
