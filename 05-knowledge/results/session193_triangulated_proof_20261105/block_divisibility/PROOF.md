# All-face module divisibility, apex-budget growth, and an exact escape

Session 193 continuation, 2026-11-05. **PROVED elementary theorems; FINITE-EXACT controls; independently reviewed.** The unrestricted Hamilton-path spectrum conjecture with exceptions `{7,21}` remains **OPEN**. This packet does not enlarge the 135-offset window or prove an unbounded interval. It proves a scalable obstruction to one proposed amplification mechanism, including *all* attachment cuts rather than a single balanced slice, and constructs an unbounded physical family escaping that obstruction at exactly the next apex budget.

## 1. Inheritance and the changed question

Let H(T) count directed Hamiltonian vertex words, with H(empty)=1. The starting result is [Session 193 power windows](../power_windows/PROOF.md): 135 fixed offsets around every sufficiently large power of two, using finite physical cores and exact response R=2^q. Its arbitrary binary continuation theorem shows that padding a fixed boundary cannot enlarge its pure-dyadic response fibre. The closest complementary construction is [Session 192 cut sumsets](../../session192_surprising_approaches_20261105/cut_sumsets/PROOF.md), which supplies actual exchange cubes and the full ordered-join cut capacities. Its cofactor obstruction applies to one boundary-balanced Cartesian slice; it expressly leaves open whether the union of different slices escapes.

Here that union is allowed. Fix an ordered join of many old cores, permit **every** cut, and then permit several arbitrary new apexes. The invariant is no longer one cut's slope, but a common divisor of ordered path-cover counts of each old core. Untouched modules contribute that divisor even when their local cut boundaries vary.

Anchor: scalable amplification of the exact windows. Niche: all-cut divisibility across many physical cores. Wildcard: the first change in the number of apexes that removes a real modular obstruction. Concept board: dyadic response fibres; physical exchange cubes; ordered multi-path covers; exposed old modules; multiplicative prime atoms; gcd-one survivors.

No literature-priority claim is made. Sources for the inherited cut and padding laws are the packets named above; the exact prime-atom reduction is [Session 193 prime induction](../prime_induction/PROOF.md), sections 2–3. No external order bound is needed here.

## 2. The exact local invariant is a path-cover fixed divisor

For a nonempty tournament V and j>=1 define

    B_j(V) = sum_(A_1 disjoint-union ... disjoint-union A_j = V)
                        product_(t=1)^j H(V[A_t]).                 (1)

The j slots are ordered and may be empty. Equivalently, B_j counts ordered j-tuples of directed paths whose vertex sets partition V. In particular B_1=H(V), and B_2 is exactly the two-path count denoted B in S192's join-capacity formula.

Let A_k(V) count ordered k-tuples of **nonempty** paths partitioning V. If n=|V| then

    B_j(V) = sum_(k=1)^n binom(j,k) A_k(V),                       (2)
    A_k(V) = sum_(j=0)^k (-1)^(k-j) binom(k,j) B_j(V),           (3)

where B_0(V)=0. Choose the occupied slots to prove (2); binomial inversion gives (3). Notice A_n=n!, so these are positive finite, computable invariants, not formal signed realizations.

For an apex budget r>=0 set

    g_r(V)=gcd(B_1(V),...,B_(r+1)(V)),
    g_infinity(V)=gcd(A_1(V),...,A_n(V)).                         (4)

The triangular integer transforms (2)–(3) show

    g_infinity(V)=gcd(B_1(V),...,B_n(V))
                 =gcd{B_j(V): j>=1}.                           (5)

Thus g_r stabilizes to g_infinity once r>=n-1. Also g_r and g_infinity are odd because they divide the odd integer H(V). A nontrivial automorphism group gives a sufficient but unnecessary source of a persistent divisor: its action on the ordered path-cover tuples is free, since fixing all tuple positions fixes every vertex. Consequently |Aut(V)| divides every B_j(V). The examples below do not rely on identifying automorphism groups.

## 3. All cuts obey a common divisor

Take nonempty tournaments V_1,...,V_m and their ordered join

    U=V_1 -> V_2 -> ... -> V_m.

Add one new vertex x with **arbitrary** incidences to U. Put

    g_a=g_1(V_a)=gcd(H(V_a),B_2(V_a)).

For each prime p, sort the numbers v_p(g_a) decreasingly, padding with zeros if necessary. Define D_1 by

    v_p(D_1)=sum_a v_p(g_a)
              -(largest valuation)-(second largest valuation). (6)

### Theorem 1: all-face divisibility

For every attachment subset S, without any restriction on its size or boundary sums,

    D_1 divides H(U+x_S).                                      (7)

Singleton blocks can occur anywhere and do not weaken the value of D_1 contributed by the other blocks.

**Proof directly from the inherited exact capacities.** Write h_a=H(V_a), b_a=B_2(V_a), and P=product h_a. The cut formula is H(U+x_S)=P+Delta_U(S). An internal capacity in block a is multiplied by P/h_a. A cross capacity between blocks a<b is

    (product_(t<a or t>b) h_t)
      [R_i(V_a)L_j(V_b) product_(a<t<b)b_t
         +[b=a+1] e_i(V_a)s_j(V_b)].                            (8)

Every term is divisible by the product of the g_t with at most the two endpoint blocks a,b omitted. The adjacency term in (8) has no intermediate blocks, so the same statement holds. Internal terms omit at most one block, and P omits none. Taking the gcd over all products with at most two omissions gives exactly (6), proving (7). A singleton has h=b_1=1, B_2=2, hence g=1. QED.

In particular, m repetitions of a fixed core V obey

    gcd(H(V),B_2(V))^max(m-2,0) divides every H(U+x_S).           (9)

This applies simultaneously to the union of *all* boundary slices, every exchange cube, and all arbitrary binary attachment choices on any extra singleton padding. It is stronger than the earlier slice-specific cofactor lattice. It still does not apply to arbitrary changes of the old block incidences, or to unrelated parent tournaments.

For mixed cores, (6), not a single common gcd, is the exact general guarantee supplied here. For example if a prime divides the local g_a for t modules, then at least t-2 of those divisibility contributions survive every one-apex cut, with their valuations retained by (6).

## 4. Several arbitrary apexes expose at most twice their number of modules

Now adjoin r new vertices, with **arbitrary** arcs among themselves and arbitrary individual arcs to every old vertex. No uniformity within an old block, no fixed attachment word, no common rooted interface, and no strongness assumptions are imposed. Only the old ordered-join arcs remain unchanged.

Put g_a=g_r(V_a). Let t=min(2r,m). Define

    D_r = gcd { product_(a notin E) g_a : E subset {1,...,m}, |E|=t }.
                                                                    (10)

Equivalently, for every prime p,

    v_p(D_r)=sum_a v_p(g_a) - sum of the largest t valuations.       (11)

### Theorem 2: bounded-apex module divisor

Every resulting tournament T satisfies

    D_r divides H(T).                                             (12)

**Proof.** Partition its Hamilton paths by their following exact boundary data:

1. the order of the r apexes in the Hamilton word;
2. which gaps between successive apexes, and the two exterior gaps, contain old vertices;
3. the actual old vertex immediately preceding or following each apex, whenever it has such a neighbor.

There are at most 2r recorded old vertices, and therefore at most 2r old blocks containing a recorded vertex. Call those blocks exposed. Adjacent apexes must have their prescribed arc; discard any inconsistent boundary data.

Removing the apexes leaves r+1 ordered old-vertex runs. Each run traverses the old blocks in increasing block order, because every interblock arc points forward. Its recorded first and/or last vertex bounds the possible blocks of that run. The first exterior run has only its last vertex constrained; the final exterior run has only its first vertex constrained. Empty runs permit no blocks.

For an **unexposed** old block V_a, these bounds specify some number q_a between 0 and r+1 of available run slots. There are no endpoint restrictions inside that block. Its vertices may be partitioned arbitrarily among those available slots, with one directed path in each chosen subset, allowing empty slots. The number of choices is exactly B_(q_a)(V_a), where B_0(V_a)=0. For an exposed block, all remaining internal endpoint constraints produce some nonnegative integer count; no divisor claim for that count is needed.

The block choices multiply. Concatenating their path pieces in block order and inserting the fixed apexes constructs a unique global Hamilton word; conversely the word recovers all pieces and the boundary data. Thus each boundary class has cardinality divisible by product_(a unexposed)g_a. This is a multiple of D_r because there are at most 2r exposed blocks; extra omissions only remove integer factors. Summing the disjoint boundary classes proves (12). QED.

For r=0 this is ordinary ordered-join multiplicativity: D_0=product H(V_a)=H(U). For r=1 it recovers Theorem 1 with a separate path decomposition proof. This gives an independent explanation of why the two endpoint blocks in the cut-capacity expression are precisely the possible divisor losses.

### The omission budget is sharp

Let C3 have outgoing rows `(2,4,1)`. Two ordered copies, plus one apex beating one chosen vertex in each copy, have H=41. Therefore no factor 3 survives in general when m=2,r=1, although g_r(C3)=3 for every r.

More generally take an ordered join of r such H41 gadgets. Regard their r gadget apexes as the new vertices and their 2r cyclic triangles as the old blocks. The resulting tournament has H=41^r, not divisible by 3. Thus no universal replacement of 2r by 2r-1 is possible, even for one fixed old-core type. The example is not strong when r>1; the theorem itself allows all tournaments.

There is also all-m strong sharpness for r=1. If every C3 block has an attachment of cardinality k_a in {1,2}, then only adjacent-block capacities can contribute at valuation m-2: those capacities equal 10*3^(m-2). All internal capacities are divisible by 3^(m-1), and nonadjacent cross capacities by 3^m. Hence

    H/3^(m-2) = sum_(a=1)^(m-1) k_a(3-k_(a+1))  (mod 3).        (13)

Choose every k_a=1 when m is not 1 mod3; otherwise choose the last k_m=2 and all other k_a=1. For every m>=2, (13) is nonzero modulo 3. These physical apex extensions are strong because the apex beats a vertex in the source block and loses to a vertex in the sink block. Therefore the exponent m-2 in (9) is attained exactly for every m>=2, not only at the smallest example.

## 5. Explicit profile values and their atom consequences

The following exact finite profiles are independently recomputed by the verifier. Rows name actual tournaments; the A-column gives every nonempty-cover coefficient, which by (2) certifies the universal j-quantifier.

| Core | Outgoing rows | H | A_1,...,A_n | g_1 | g_2 | g_infinity |
|---|---|---:|---|---:|---:|---:|
| C3 | `(2,4,1)` | 3 | `3,6,6` | 3 | 3 | 3 |
| C5 | `(6,12,24,17,3)` | 15 | `15,90,210,240,120` | 15 | 15 | 15 |
| V9 | `(20,1,2,7,14)` | 9 | `9,66,186,240,120` | 3 | 3 | 3 |
| V13 | `(18,8,3,5,14)` | 13 | `13,78,198,240,120` | 13 | 1 | 1 |
| V25 | `(2,12,9,49,7,23)` | 25 | see receipt | 1 | 1 | 1 |

The inherited prime-atom induction separates odd primes other than 7 and targets 7p for odd primes p>=5,p!=7, plus the finite atoms 49,63,343. Theorems 1–2 give immediate constraints on proposed modules for those lanes:

- If m>=2 cyclic triangles occur as old blocks and m>2r, every output is divisible by 3 and is at least 3^m, so it lies in neither infinite atom lane. The same conclusion holds for repeated V9 cores. Arbitrary singleton padding does not help.
- Repeated C5 blocks give 15^max(m-2r,0) dividing every output. With m>2r, the factor 15 already excludes both infinite atom lanes.
- With only one apex, at least three repeated V13 blocks force divisibility by 13. The output is at least 13^m, so it is neither a prime nor the only 7p atom divisible by 13, namely 91.

The lower bound H(T)>=product H(V_a) follows by the elementary insertion injection: insert each new vertex into a legal slot of each old Hamilton path, choosing one deterministic legal slot; deleting it recovers the previous path. This is not a use of an asymptotic strong-tournament order bound.

Consequently the apex budget must grow at least linearly with the number of repeated C3, C5, or V9 modules in any construction aimed at new members of the infinite atom lanes. This obstructs a fixed number of apexes over arbitrarily many such cores, even while every cut face and every local exchange is allowed. It does **not** obstruct arbitrary tournaments, unlimited apex complexity, new old-arc patterns, or gcd-one cores.

## 6. Even gcd-one repeated cores have an all-cut sparsity bound

There is a different obstruction that does not require a nontrivial divisor.
Let V be a fixed strong tournament of order n with h=H(V)>1, and let B=B_2(V).
Allow every one-apex attachment to every repeated ordered join V^m, m>=1.
Write F_V for the set of counts of the **strong** outputs.

### Theorem 3: all-cut repeated-core sparsity

One has B>2^n. With

    alpha = n log(2)/log(B) < 1,

    #{N in F_V : N<=X} = O_V(X^alpha),                       (14)
    sum_(N in F_V) N^(-s) < infinity for every s>alpha.      (15)

**Proof.** Every subset of V has at least one Hamilton path. Thus every one
of the 2^n terms in B=sum_A H(A)H(V\A) is at least one. The empty and full
subsets contribute h>1, giving B>=2^n+2(h-1)>2^n.

For m>=2, strongness of the apex extension requires the apex to beat some
i in the first block and lose to some j in the last block. The corresponding
outgoing cut capacity (8) is at least

    R_i(V)L_j(V) B^(m-2) >= h^2 B^(m-2).                  (16)

Here R_i>=h and L_j>=h: split every Hamilton path immediately before i,
respectively after j, to inject it into the appropriate ordered two-path
covers. Consequently H>=h^2 B^(m-2). There are at most 2^(nm) attachment
masks at length m. If H<=X, (16) bounds m by
`2+log_B(X/h^2)`. Summing the geometric presentation bound through that
length proves (14). For (15), sum directly over presentations:

    sum_(m>=2) 2^(nm) [h^2 B^(m-2)]^(-s)

is a geometric series with ratio 2^n/B^s<1. The finite m=1 bank adds no
problem. Counting repeated presentations only increases either upper bound.
QED.

For finitely many fixed strong core types, take the **union of their separate
repeat families** and finitely many extra nontrivial seeds. The resulting
bank still has a convergent reciprocal sum. Its multiplicative closure also
has one, since its sum of reciprocals is bounded by

    product_(b in bank) (1-1/b)^(-1) < infinity.           (17)

The convergence follows from sum_b 1/b<infinity and b>=3. Therefore this
closure misses infinitely many ordinary primes and infinitely many 7p
atoms: each of those infinite lanes has divergent reciprocal sum, by the
prime reciprocal divergence recalled in [Session 193 prime induction](../prime_induction/PROOF.md).
This conclusion includes unlimited ordered joins of the produced witnesses.

This extends the earlier fixed-clock sparsity obstruction to an exponential
number of independently chosen root incidences at each module count. In
particular it applies to the gcd-one V25 core, for which n=6 and B=304.
Coprimality is not an exact-coverage theorem.

### A mixed-core criterion using the actual cut interfaces

For a fixed core V, let q(V) be the number of distinct records

    (Delta_V(S), R(S), L(S), e(S), s(S)),  S subset V.

These are the S192 physical local cut value and its four boundary sums.
The exact join-capacity law (8) shows that the records of the blocks determine
the complete final H. Thus for a specified word V_1,...,V_m of core types,
the number of possible H values is at most product_a q(V_a), not necessarily
the much larger number of masks. Every record is witnessed by an actual mask;
this bound makes no independently realizable-sumset assumption.

Let L be any finite library of nontrivial strong cores. If for some s>0,

    sum_(V in L) q(V) B_2(V)^(-s) < 1,                      (18)

then the strong-count bank from **all arbitrary mixed core sequences** in L,
with **all one-apex cuts**, has a convergent sum of N^(-s). The same holds for
its multiplicative closure, and its counting function is O_L,s(X^s).

To prove this, put c=(min_(V in L) H(V)/B_2(V))^2>0. The first-to-last cut
argument gives H>=c product_a B_2(V_a) for every strong output with m>=2.
Sum the record bound over all core words of length m. The contribution to
the Dirichlet sum is at most

    c^(-s) [sum_(V in L) q(V)B_2(V)^(-s)]^m.

This is geometric under (18). The length-one bank is finite. The Euler-product
argument (17), with 1/b replaced by b^(-s), proves the closure assertion.
Finally N<=X implies 1<=X^s N^(-s), proving the counting bound.

For the four actual cores C3,C5,V13,V25, exact enumeration of their local
records gives respectively

    q = 4,8,32,64,        B_2 = 12,120,104,304,
    sum q/B_2 = 1134/1235 < 1.                            (19)

The inequality can be strengthened to the exact exponent s=99/100 without
relying on a floating-point estimate. Indeed

    304*16^100 < 17^100,
    (17/16)*(1134/1235) = 19278/19760 < 1.

Consequently every B_2^(1/100)<17/16 and (18) holds at s=99/100. Thus the
entire arbitrary mixed-core bank on these four cores, **and its unlimited
ordered-join closure**, has O(X^(99/100)) distinct counts through X and
misses infinitely many members of both infinite atom lanes. This explicitly
includes a gcd-one core and many physical cut faces.

**Scope matters.** An arbitrary finite library need not satisfy (18); adding
core types can destroy the inequality. Nor does this argument allow arbitrary
internal singleton padding with freely chosen incidences. A singleton has
two distinct records and B_2=2, so it alone contributes 1 at s=1 and makes the
strict library bound fail. The all-face divisibility theorems do allow such
padding, but this sparsity theorem does not. Old-arc changes and more apexes
are also outside this assertion. The family E_m below has two apexes and is
not being ruled out by this one-apex bound.

## 7. The second apex really removes the V13 divisor

The profile change g_1(V13)=13 to g_2(V13)=1 has an exact physical consequence at every scale, not only a formal disappearance of the bound.

Take m>=2 ordered copies of V13. Add apexes a->b. Both apexes beat every vertex of the **first** block; every later block beats both apexes. Call this tournament E_m. It is strong: the apexes reach the first strong block, that block reaches all later blocks, and the last block returns to both apexes.

### Theorem 4: all-clock two-apex escape

For every m>=2,

    H(E_m) = 552*471^(m-1) - 675*104^(m-1) - 7*13^(m-1),       (20)
    H(E_m) = 2*3^m (mod 13),                                  (21)

so 13 never divides H(E_m).

**Proof by the same path-cover interface.** More generally let a repeated core V have h=B_1, b=B_2, c=B_3, and a_1=A_1=h, a_2=A_2, a_3=A_3. Put t=m-1>=1. In either apex order a Hamilton word has three old-vertex runs L0, apex, L1, apex, L2. In the first block let J be the nonempty set of occupied runs; its contribution is a_|J|. In the later t blocks:

- L1 and L2 can be used only if they already started in the first block;
- if L0 was occupied in the first block, a later vertex must repair its final incidence into the first apex;
- if L1 was occupied in the first block, a later vertex must repair its final incidence into the second apex.

An empty L1 permits only the order a,b; a nonempty L1 permits both apex orders. Inclusion-exclusion for these required nonempty later slots gives the following seven contributions:

    J={0}:       a_1 h^t
    J={2}:       a_1 b^t
    J={0,2}:     a_2(b^t-h^t)
    J={1}:       2a_1(b^t-h^t)
    J={0,1}:     2a_2(b^t-2h^t)
    J={1,2}:     2a_2(c^t-b^t)
    J={0,1,2}:   2a_3(c^t-2b^t+h^t).

Their sum is

    2(a_2+a_3)c^t +(3a_1+a_2-4a_3)b^t
                         +(-a_1-5a_2+2a_3)h^t.               (22)

For V13, `(a_1,a_2,a_3)=(13,78,198)` and `(h,b,c)=(13,104,471)`, giving (20). Modulo 13, only the first term survives: 552=6 and 471=3 mod13, giving 6*3^(m-1)=2*3^m. QED.

The t>=1 condition is essential: at m=1 the apexes are just a source chain and H=13, whereas (20) does not apply. The first two values are verified on the actual 12- and 17-vertex tournaments by subset DP. The module transfer supplies many further exact controls; the all-order claim rests on the run decomposition (22).

This construction escapes the repeated-V13 one-apex restriction by allowing three old-vertex path runs instead of two. It does not produce a growing consecutive interval, does not claim that infinitely many of its values are prime, and is not closed under the prime-lane inverse-parent induction. It is a real boundary crossing for one modular obstruction.

## 8. What this changes for scalable window growth

The original 135-offset bank has 132 distinct core tournaments. A separately retained finite inventory finds g_1=1 for 85 of those cores. Thus the new theorem does **not** dismiss that whole bank. It says that a scalable search using ordered modules should explicitly track `B_1,...,B_(r+1)` and reject a proposed prime-atom construction when (11) forces a forbidden factor. The V25 example gives a small gcd-one core that survives this test already at one apex; survival is not coverage.

Three losses are now separated:

1. **Fixed physical boundary plus binary padding:** the S193 response-injectivity theorem prevents enlarging a pure-dyadic fibre.
2. **Many physical modules with arbitrary cut faces:** Theorem 1 can force an all-face divisor, stronger than the old balanced-slice cofactor.
3. **More apexes:** Theorem 2 measures how much module divisibility the added endpoints can remove; Theorem 4 shows that changing the path-cover budget really can remove a divisor.

Theorem 3 further shows that merely repeating one gcd-one core with every possible cut is still too sparse for either infinite atom lane. Its record-count criterion also rules out arbitrary mixtures from the explicit four-core library. Larger mixed libraries failing (18), and arbitrary internal singleton words, remain unclassified by this argument.

A completion still needs exact coverage of all remaining prime and 7-prime targets. One viable unproved target is an unbounded cut/intercept construction on gcd-one old cores; another is a target-dependent construction whose apex budget grows fast enough. Neither follows from coprimality, complete residues, a large finite window, or the escape family (20). The full spectrum and a scalable interval-filling mechanism remain OPEN.

## 9. Reproduction and audit boundaries

Run:

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

The verifier has explicit exceptions, not removable assertions. It independently computes ordered cover profiles by subset convolution; reconstructs nonempty-cover coefficients by binomial inversion; checks full one-apex cut universes for six named parent families; compares selected actual tournaments with direct vertex-subset DP; tests arbitrary two-apex incidences; verifies exact valuation sharpness on physical all-m examples; and checks the two-apex module transfer against direct DP and literal permutations. The receipt states every universe and finite check count.

The [power-window seed bank](../power_windows/seeds135.json) is used only for the scoped 132-core inventory, not for any theorem above. The elementary proofs and the standalone verifier are sufficient for the new statements. Independent review checked the unexposed-block factorization in Theorem 2 and the seven-term endpoint table in Theorem 4.
