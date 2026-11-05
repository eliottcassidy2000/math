# Transitive cloning: matching rank, binary amplifiers, and a prime-coverage barrier

Session 196, 2026-11-05. **All-order exact physical construction and geometric classification. The full `{7,21}` Hamilton-path spectrum remains open.** No claim of literature priority is made.

## 1. Why change the operation?

The reviewed [S195 guide](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/00-navigation/S195.md) develops rooted one-seam transfers. Its [rooted interface](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/rooted-interfaces/PROOF.md) has a closed three-coordinate state, but a scalar binary or ternary return need not be available from every state. Its [path-cover algebra](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/pathcover-algebra/PROOF.md) supplies a more general run decomposition and explicitly uses `k! S(m,k)` for ordered path covers of a transitive tournament.

Here the physical operation replaces **one vertex** by a transitive module. It is not another attachment from that fixed seam bank. The new conclusions are an exact matching formula for the largest possible number of clone-runs; a graph-theoretic classification of genuinely binary affine amplification; a physical same-sign return; and an all-order prime-coverage obstruction for finite template banks even when arbitrary ordered products are admitted. The underlying maximal-run substitution mechanism is inherited, not claimed as new.

This does not contradict the [S194 whole-input module barrier](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session194/module-barrier/PROOF.md): it is a single vertex that is replaced by a module, not the entire nontransitive input. The original tournament is an induced subtournament of the result after choosing one clone. Consequently induced-core endpoint floors still apply.

All tournaments are finite, and `H(T)` counts directed Hamiltonian vertex words. An empty tournament has count one. Bit j of adjacency row i means `i -> j`.

## 2. A literal clone-run expansion

Let `(T,r)` be a nonempty rooted tournament. For `m>=1`, write `T[r:m]` for the tournament obtained by replacing r by a transitive m-vertex chain. Every clone has the same incidences to `V(T)-r` that r had. Put `h_m=H(T[r:m])`.

Define `a_k(T,r)` as the number of legal words using every vertex of `V(T)-r` once and k indistinguishable r-symbols, with no two r-symbols consecutive. Legality uses the original arcs at every consecutive pair. Thus `a_1=H(T)` and all `a_k` are nonnegative integers. There are at most `|T|` r-symbols.

**Proposition 1 (exact expansion).** For every positive m,

    h_m = sum_(k>=1) a_k(T,r) k! S(m,k),                 (1)

where `S(m,k)` is a Stirling number of the second kind.

**Proof.** Cut a Hamilton word in `T[r:m]` into maximal clone-runs. Compress each run to r. The result is exactly one word counted by `a_k`. Conversely distribute the m distinct clones among the k ordered, nonempty runs, then put the clones in each run in their unique increasing order. There are `k! S(m,k)` ordered set partitions. The cross incidences are legal independently of the partition. These operations are inverse. This is an exact identity on actual graphs, not a fitted recurrence. QED.

**Positive multiplicities only.** None of these formulas interprets m=0 as deletion. For a directed triangle, the positive-m formula is `2^m+1`, whose formal value at zero is two, while deleting the root leaves a transitive edge of count one. The same restriction applies to the synthesis below: Z is nonempty.

## 3. The run degree is a bipartite matching number

Let

    O = {x: r -> x},       I = {y: y -> r}.

Form the bipartite graph `B_r` with sides O,I and edge xy exactly when `x -> y` in T. Its edges are precisely the directed triangles `r -> x -> y -> r` through r. Let `nu_r` be its matching number.

**Theorem 2 (matching rank).**

    max{k: a_k(T,r)>0} = 1 + nu_r.                     (2)

Equivalently, the maximum number of separated clone-runs possible with enough clones is one plus the maximum number of pairwise vertex-disjoint `O -> I` edges.

**Upper bound.** Between any two consecutive r-symbols in a legal quotient word is a nonempty directed path. Its first vertex is in O and its last vertex is in I. Along that path there is an O-to-I transition, and the corresponding arc is an edge of `B_r`. Select one such edge from each gap. The gap vertex sets are disjoint, so the selected edges form a matching. Thus `k-1<=nu_r`.

**Constructive lower bound.** Choose a maximum matching `x_i -> y_i`, for `1<=i<=nu_r`. In the actual tournament with `nu_r+1` clones start with the partial path

    r_0, x_1, y_1, r_1, x_2, y_2, r_2, ..., x_nu, y_nu, r_nu.

Every displayed edge is legal. Insert the remaining named vertices one by one into this directed path. A vertex can always be inserted: put it first if it beats the first vertex, put it last if the last vertex beats it, and otherwise the orientation toward it changes somewhere from an old vertex beating it to it beating the next old vertex. Insert there.

Insertion preserves the relative order of all existing vertices. In particular, each successive pair of clones remains separated by its original named x_i,y_i. There are no additional clones to insert. The resulting Hamilton path therefore has exactly `nu_r+1` nonempty clone-runs. Their endpoints may change during insertion, which is harmless; only their continued separation is required. Compress the runs to obtain `a_(nu_r+1)>0`. This argument includes the hostile case where `B_r` is just a matching of disjoint edges. QED.

**Exponential interpretation.** The standard inclusion-exclusion formula

    k! S(m,k) = sum_(j=0)^k (-1)^(k-j) binom(k,j) j^m

makes `h_m` an integer linear combination of `1^m,...,(nu_r+1)^m`. Its highest-base coefficient is exactly `a_(nu_r+1)>0`. Distinct positive exponential sequences are linearly independent, for example by a Vandermonde determinant. Hence the highest exponential base really is `nu_r+1`; negative lower coefficients cannot cancel it.

The matching number is also retained at **every clone** after further transitive cloning of r. In the cross graph at a chosen clone, earlier and later clones create no new cross edges: an earlier clone beats a later one, every original I vertex beats every clone, and every clone beats every original O vertex. The remaining cross graph is exactly the old `B_r`, with extra isolated vertices. This is a graph-level return property.

## 4. Precisely when cloning is binary affine

A nonempty bipartite graph has matching number one if and only if all its edges share one vertex. Indeed, two edges must meet; if edges `xy` and `xy'` have different y-endpoints, any edge not incident with x is disjoint from one of them. The one-edge case is included.

**Corollary 3 (binary amplifier classification).** The following are equivalent:

1. `nu_r=1`.
2. There is at least one directed triangle through r, and all directed triangles through r share another vertex c.
3. There exist integers `alpha>0` and beta such that for every `m>=1`,

       H(T[r:m]) = alpha 2^m + beta.                  (3)

The coefficients are exactly

    alpha = a_2 = (h_2-h_1)/2,
    beta = a_1-2a_2 = 2h_1-h_2.                      (4)

**Proof.** The equivalence of 1 and 2 is the star criterion. When `nu_r=1`, (1) reduces to `a_1+a_2(2^m-2)`, and `a_2>0`. Conversely (3) and the positive highest exponential coefficient in Theorem 2 force the largest base to equal two. QED.

The case `nu_r=0` is separate: `h_m=H(T)` for all m, so it is a constant sequence and not a positive binary amplifier. There is no assertion that beta is plus or minus one merely because the star condition holds. For example

    rows=[16,1,3,7,14], root=0

has `alpha=7,beta=-5`, hence counts `9,23,51,107,...`. Another five-vertex star input can have `beta=3`. Binary exponential growth and literal `2H±1` are distinct conditions.

Every tournament Hamilton-path count is odd; equivalently here `beta=h_1-2alpha` is odd. For the actual repeated operation, inserting one more clone gives

    h_(m+1) = 2h_m-beta.                              (5)

Thus beta equal to -1 supplies a literal `2H+1` operation, and beta equal to +1 supplies a literal `2H-1` operation. Because cloning a clone simply enlarges the same transitive module, the **same sign** returns at every step. This does not automatically produce the opposite sign from the same input.

### Explicit returning examples

The strong four-vertex tournament

    D=[8,1,3,6], H(D)=5

has two different choices.

- Root 0: `alpha=3,beta=-1`, so `H(D[0:m])=3*2^m-1`. The counts are `5,11,23,47,95,191,383,767,...`; every step is `2H+1`.
- Root 1: `alpha=2,beta=1`, so `H(D[1:m])=2^(m+1)+1`. The counts are `5,9,17,33,65,129,257,513,...`; every step is `2H-1`.

The strong H13 tournament

    E=[2,12,9,17,7]

at root 0 has `alpha=6,beta=1`. It gives

    H(E[0:m])=3*2^(m+1)+1,
    13,25,49,97,193,385,769,1537,... .

For m=2 its actual rows are `[34,12,41,49,39,2]`. These are physical strong graphs at every m: a transitive expansion of a vertex of a strong tournament remains strong, because a directed cycle through that vertex gives every clone a return path through the exterior. The formulas are not numerical-primality claims, nor are they asserted to be new scalar families in the literature.

### Why two fixed clone buttons cannot return both signs

**Proposition 4 (commutation obstruction).** Two fixed distinct clone classes cannot both retain literal `2H+1` and `2H-1` operations after every possible sequence of cloning either class.

Cloning the first class and then the second produces the same tournament up to relabeling as doing them in the opposite order. But the claimed two scalar maps do not commute:

    2(2h-1)+1=4h-1,       2(2h+1)-1=4h+1.

This is a contradiction. The proposition fixes the two clone classes and their signs. It does not forbid changing the chosen root or operation according to the full current graph, and does not turn into a general impossibility theorem for all binary constructions.

## 4A. A complete two-piece synthesis and a unit-determinant target

The star characterization is constructive in both directions. Take any two nonempty rooted tournaments `(X,c)` and `(Y,c)`, identifying only their root. Write `A=X-c`, `B=Y-c`, and add a transitive m-chain Z. Orient

    A -> Z -> B,      A -> B,      Z -> c,             (6)

and keep all X and Y arcs. Empty A or B is allowed. Call this actual tournament `F_m(X,Y)`.

Use the rooted cover counts from the S195 interface: `R_c(X)` counts ordered spanning covers `(P,Q)`, with P possibly empty and Q nonempty starting at c. `L_c(Y)` counts covers `(P,Q)`, with P nonempty ending at c and Q possibly empty. Put

    h=H(X), e=e_c(X), d=H(X-c), R=R_c(X),
    g=H(Y), t=H(Y-c), L=L_c(Y),
    U=R-d,       sigma=R-h-d+e.

**Theorem 4A (complete star synthesis).** For every positive m,

    H(F_m(X,Y)) = U L 2^m + d g - sigma t.            (7)

The vertices of Z form a positive binary amplifier exactly when `U>0`. Every positive binary amplifier in Corollary 3 is obtained from this construction, possibly after taking the converse tournament. Consequently its two coefficients have the exact form

    alpha = U L,       beta = d g - sigma t.          (8)

**Geometric completeness.** Suppose a root r has a nonempty star cross graph, with center c. If c is an in-neighbor of r, take the converse tournament; Hamilton counts and clone-run counts are unchanged. We can now assume `r->c`. Put `A=I(r)` and `B=O(r)-{c}`. Because every O-to-I edge is incident with c, every edge between A and B points `A->B`. Also `A->r->B` and `r->c`. Thus the original graph is precisely `F_1(X,Y)` with `X=T[A+c]`, `Y=T[B+c]`, and its r-cloning family is (6). Conversely (6) has no cross edges at a Z-vertex except those from c to A, so its cross graph is a star. `U` counts ordered two-path covers of A whose second path is nonempty and whose first vertex is beaten by c. It is positive precisely when c beats some A vertex: a singleton second path and any Hamilton path on the remaining A vertices witness positivity. This proves the claimed positive-amplifier condition and completeness.

**Exact coefficient proof.** Every Hamilton word can be written `P,c,Q`. Each of P,Q traverses its A-, Z-, and B-parts in that order, with empty parts permitted.

For both P and Q to contain Z vertices, the A-part of Q must be nonempty and start at an out-neighbor of c, since `Z->c` rules out going directly from c into Z. The ordered A-cover has U choices. The B-part of P is either empty or ends at an in-neighbor of c; the ordered B-cover has L choices. The Z vertices can be assigned to their two nonempty increasing runs in `2^m-2` ways. This gives `UL(2^m-2)`.

When all Z vertices occur in P, there are two disjoint cases. If Q has A vertices, the same U and L cover choices give UL. If Q has no A vertices, all A lies in P and the B+c subword is an arbitrary Y Hamilton path, giving dg.

When all Z vertices occur in Q, its A-part must be nonempty. If P has B vertices, there are `U(L-t)` choices. If P has no B vertices, all B lies in Q and the X subword is a Hamilton path not ending at c, giving `t(h-e)`.

Adding these four contributions and the two-run contribution gives

    UL(2^m-2) + UL + dg + U(L-t) + t(h-e)
      = UL 2^m + dg + t(h-e-U)
      = UL 2^m + dg - sigma t,

as required. This proof uses actual path runs rather than formal subtraction as a construction. QED.

The inherited zero-slack theorem in the [S195 guide](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/00-navigation/S195.md) implies `sigma>=0`, but the counting identity and its proof do not require that inequality as an assumption.

In particular, a unit-offset binary ray is synthesized when two **actual** rooted pieces satisfy

    d H(Y) - sigma H(Y-c) = +1 or -1.                 (9)

This is a determinant/Bézout condition with a simultaneous amplitude constraint `alpha=(R-d)L_c(Y)`. Coprimality of d and sigma is necessary for (9), but not sufficient: neither an arbitrary integer Bézout solution nor an arbitrary desired cover response need be realized by a tournament. Thus (9) is a concrete new inverse-target obligation, not an all-target solution.

For example,

    X=[18,8,3,5,14], root 0;   Y=[0], root 0

gives `alpha=15,beta=-1`, with actual `F_1` rows

    [18,40,35,37,46,1].

Its counts are `29,59,119,239,479,959,...`. The other sign at the same amplitude is obtained from

    X=[10,0,3,6], root 0;      Y=[4,1,2], root 0,

giving `alpha=15,beta=1` and counts `31,61,121,241,481,961,...`. The first seven-vertex rows are `[42,112,115,118,1,16,49]`. These examples demonstrate actual synthesis from distinct rooted pieces. They are separate same-sign rays, not two returning buttons on a common parent.

## 5. Finite cloning banks cannot cover the required prime atoms

A stronger all-order limitation allows arbitrary root changes and any number of clone steps, provided no operation creates a new template.

**Theorem 5 (fixed-template bound).** Fix a tournament Q. Let s be the number of vertices of Q in nontrivial strong components. Among all positive transitive blowups

    Q[Tr_(m_1),...,Tr_(m_n)],        m_i>=1,

the number of distinct Hamilton counts at most X is at most

    max(1, floor(log_2(X-1))^s)                       (6)

for `X>=2`. A finite bank has the sum of these bounds.

**Proof.** Every vertex v in a nontrivial strong component lies in a directed triangle. To see this, divide that component into the out- and in-neighbors of v. If there is no O-to-I edge, O cannot reach v, contradicting strongness. Choose one such triangle.

In a blowup, the m_v clones of v and one clone of each of the other two triangle vertices induce the cyclic three-block tournament `C3[Tr_(m_v),1,1]`. Its Hamilton count is `2^(m_v)+1`: apply (1) to a triangle, whose quotient coefficients are `a_1=3,a_2=1`.

Hamilton paths of an induced subtournament inject into those of the full tournament, by deterministic insertion of each missing named vertex. Deleting the inserted vertices recovers the original path. Therefore `H>=2^(m_v)+1`, so every active multiplicity satisfies `m_v<=floor(log_2(X-1))`.

The remaining vertices of Q are singleton strong components. Blowing up one of them merely inserts a longer transitive chain in the linear order of strong components and does not change H: every Hamilton path is the concatenation of component paths. Those multiplicities can be discarded. There are at most the stated number of active multiplicity vectors. QED.

**Corollary 6 (ordered products do not repair prime coverage).** Start from any fixed finite bank. Permit arbitrary transitive vertex cloning, arbitrary choice of the root at each step, and arbitrary ordered joins of previously built tournaments. Then infinitely many odd primes are still not Hamilton counts of the generated class.

Every nontrivial strong component of any output is a transitive blowup of a strong component in the finite original bank. Cloning preserves that assertion, and an ordered join only concatenates the component lists. A prime total Hamilton count requires exactly one component to have that prime count and all others to have count one. By Theorem 5, the number of available prime counts at most X is `O((log X)^s)` for a fixed s.

For completeness, the number of primes is not bounded by any fixed power of log X for all X. Such a bound would make `sum_p 1/p` converge, by dividing primes into dyadic intervals. Convergence would bound the Euler products `prod_(p<=Y)(1-1/p)^(-1)` uniformly. But expanding each finite product includes `1/n` for every `n<=Y`, contradicting divergence of the harmonic series. This is enough to show infinitely many prime omissions. No prime number theorem is needed.

This limitation applies only to finite **cloning** template banks, even with ordered products. It does not apply to the S195 one-seam grammar, whose new cross arcs can continually change the template. Nor does it preclude using target-dependent, unbounded template banks; doing so would leave the inverse-target problem to those templates.

## 6. Verification and exact scope

Run in this packet:

    python3 verify.py > normal.stdout
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout

The checker uses explicit exceptions, so its checks remain active in optimized Python. It enumerates all 1099 labeled tournaments of orders one through five and every root, without isomorphism or value filtering. It independently computes:

- the literal coefficients a_k by permutations of the exterior vertices and all sets of insertion gaps;
- the bipartite matching number by augmenting paths;
- actual Hamilton counts of inflated adjacency matrices by subset endpoint DP;
- the star condition and preservation of matching number at every clone;
- explicit insertion witnesses from an independently brute-forced maximum matching on a specified sub-universe;
- a five-vertex hostile whose cross graph has two disjoint edges and therefore a strictly positive third-run term;
- the basic affine families through m=8 and simultaneous positive clone expansions of three strong templates;
- all 841 rooted pairs X,Y of orders at most three, at m=1,2,3, checking the two-piece synthesis on 2,523 actual graphs, and both synthesized examples through m=6.

The finite table is a test universe, not a strong-order classification. In particular the observations that the listed labeled five-vertex H13 inputs have only minus amplifiers and that H15 has no unit-offset amplifier at that order do not assert the corresponding all-order impossibilities. The all-order conclusions in this packet are Theorem 2, Corollary 3, Proposition 4, Theorem 4A, and Theorem 5 / Corollary 6, proved above.

The required open step remains an inverse construction for every required prime and 7p atom. This operation supplies exact physical binary rays and a geometric diagnostic, but its finite-template closure cannot supply that theorem.
