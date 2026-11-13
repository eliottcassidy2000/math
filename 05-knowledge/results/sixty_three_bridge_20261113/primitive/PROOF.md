# A two-slot explanation of63 and a primitive seven-adic bridge

**PROVED:** the arbitrary-module formulas, sharp63-multiplier fibre, literal
marked continuation, and seven-adic residue theorem below. **FINITE-EXACT:**
the small independent subset-DP controls. The two global63 classes and the
7p reduction are **inherited proved inputs**. **OPEN:** exact prescribed
cofactor coverage and Chan–Pak. No assertion of infinitely many primes in
the displayed exponential sequences is made.

## 1. What63 is already known to mean

The global inventory in
[Session24, §§3–5](../../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md)
has exactly two strong ordinary H63 classes, both order8, a converse pair.
Since7 and21 are impossible and Hamilton counts multiply along the ordered
strong-component decomposition, every tournament with H63 is exactly
`TT_a -> C -> TT_b`, where C is one of those two cores. Consequently its
ordinary multiplicity at ordern is2(n−7), n≥8. This all-order conclusion
uses the inherited complete strong inventory and its strong-order bound;
the present packet does not repeat the census.

Their W5 interpretation is also already known:
[Session230, §4](../../session230_creator_grammar_20261109/primes/REPORT.md)
identifies the two cores as W5 substitutions with block sizes
(1,2,2,2,1) and (2,2,2,1,1). The one-apex signatures1001100 and1100110
occur in
[THM-4748](../../../../01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md).
Thus63 is a numerical multiple of7 without a strong component of count7;
it is not a new substitution-prime object.

The exact prime/7p reduction is inherited from
[THM-4094, §§5–6](../../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
The strong49,63,343 carry atoms are already sufficient to reduce the full
conjecture to the two prime lanes stated there. The question here is what
portable marked operation the63 atom supplies, rather than rediscovering
that reduction.

## 2. An arbitrary-module, two-slot formula

Let M be any **nonempty** tournament. Write h=H(M) and p=P_M(2), where
P_M(2) counts ordered pairs of spanning directed paths, allowing either
slot to be empty. The empty induced tournament has H=1 in this definition.

Form an ordinal chain of four nonempty blocks

    A -> M -> C -> D,       A=TT1, C=TT2, D=TT2,

and add z with z→A,C and M,D→z. Denote the result B(M). Keep z as the
named root. The quotient is the inherited W5, so B(M) is strong even if M
is not strong. There are |M|+6 vertices. The exact identities are

    H(B(M)) = 13p+11h,
    state_z(B(M)) = (13p+11h, h, h, 26p−h).             (1)

Our state convention is (H,s,q,F): s counts Hamilton paths starting at
the root, q=H after deleting it, and F counts ordered spanning covers
(P,Q), where P is nonempty and ends at the root and Q may be empty.

**Proof of the first coordinate.** A Hamilton path visits M in at most
two runs: its only possible outside predecessor is the single vertex A.
A one-run path compresses to a path in B(TT1), which has37 paths. For a
two-run path the first run starts the entire path; between its two runs
lie z,A in that order, possibly preceded by some C,D vertices. That
prefix of C,D is either empty, or has a nonempty D part, since C loses
to z. The choices are exactly

    1 + (2^2)(2^2−1) = 13.

All selected C,D vertices occur in their transitive order; unselected
ones follow the second M run. Therefore, with C_2(M)=p−2h counting
ordered covers in two **nonempty** slots,

    H(B(M))=37h+13C_2(M)=13p+11h.

For completeness,37 follows directly from the inherited one-apex path
partition. A four-block word1^a0^b1^c0^d has

    3+(2^a−1)(2^b−1)
      +2^(b+c)(2^a−1)(2^d−1)+(2^c−1)(2^d−1).       (2)

The3 counts the two endpoint placements and the single adjacent01 cut;
the other terms count earlier1/later0 pairs with weight2^(distance−1).
Substitute (a,b,c,d)=(1,1,2,2).

Deleting z leaves an ordinal chain, giving q=h. A Hamilton path starting
at z must next enter A and then traverse the remaining ordered blocks,
giving s=h. To compute F, the deleted chain has two-slot count32p. Its
first-slot ending covers at vertices beaten by z number h+6p: endpoint A
contributes h; the two C endpoints contribute2p and4p. All nonempty
first-slot ending covers number32p−h. Append z to one ending in M or D,
and add the singleton first-slot case{z}, of weight h. Thus

    F=h+(32p−h)−(h+6p)=26p−h.

These are universal path decompositions, not interpolation from module
samples. Internal paths of M retain their order and multiplicity.

### Exactly which modules give the multiplier63?

The identity H(B(M))=63H(M) holds **if and only if** M is TT2 or C3.
Indeed (1) requires p=4h. Cutting each Hamilton path at its |M|+1
positions injects into its ordered two-slot covers, so p≥(|M|+1)h.
Hence |M|≤3. The complete possibilities at those sizes are

| M | h | p |
|---|---:|---:|
| TT1 |1|2|
| TT2 |1|4|
| TT3 |1|8|
| C3 |3|12|

Only TT2 and C3 satisfy p=4h. The outputs are strong graphs of counts63
and189 respectively, at orders8 and9. This is a sharp fibre of the actual
two-slot operation; it is not an SCC multiplication construction. In
particular63 arises from the positive sum37+13·2, with no graph of count7
or21 as an input factor.

## 3. A physical primitive bridge from63 to7·71

Use the inherited root-child operation: add x with root→x→every other
old vertex. R retains the old root; N instead retains x. Their exact maps
are

    R(H,s,q,F)=(H+F−s,q,q,2F−q),
    N(H,s,q,F)=(H+F−s,H−s,H,H+F).                    (3)

See [ROOTED_RESETS R1](../../two_port_return_induction_20261103/ROOTED_RESETS.md).
Every birth preserves all old arcs and strongness for the nontrivial
strong inputs used here. Apply N to (1):

    (39p+9h, 13p+10h, 13p+11h, 39p+10h).            (4)

The new root has indegree1 and positive readiness defect h, since both
q−s and F−H equal h. Two more births, N then R, give the balanced state

    (b,t,t,b),
    b=117p+29h,  t=39p+9h,  b−3t=2h.               (5)

When M=TT_j, j≥1, we have h=1,p=2^j. The graph (4) is unit-ready and
(5) is primitive: gcd(b,t)=gcd(2,t)=1 because t is odd. The concrete j=2
chain is

    (63,1,1,103)
      --N--> (165,62,63,166)
      --N--> (269,103,165,331)
      --R--> (497,165,165,497).                     (6)

Thus the same physical eight-vertex63 source remains induced throughout
an eleven-vertex strong realization of497=7·71. The primitive claim is
about the full marked state, not merely the integer497. The existing
last root-child pair is primitive too. Immediately before the final R,
for M=TT_j the state is

    (65·2^j+9, 26·2^j−1, 39·2^j+9, 78·2^j+19).

Its first coordinate minus its second and third equals1, so the parent's
full coordinate gcd is1. Apply the inherited root-child sealing theorem
[ROOTED_RESETS R3](../../two_port_return_induction_20261103/ROOTED_RESETS.md)
to that literal last birth. Separately every **fresh seal** of the final
balanced root also has continuation content1, by gcd(b,t)=1. These are
two different uses of the same interface theorem; no forgotten old pair
is substituted for a new seal.

For general M, (5) has gcd(b,t)=gcd(h,39p): h is odd, and p is even
because swapping the two slots has no fixed cover for nonempty M. Thus
t is odd, and gcd(b,t)=gcd(2h,t)=gcd(h,39p). No automatic primitivity is
asserted. The special transitive input deliberately removes that sidecar
obstruction. The already proved finite connector in
[CHAIN_REUSE, §3](../../h17_exact_repairs_20261112/constructors/CHAIN_REUSE.md)
connects every primitive balanced source to the inherited repeating
controller after its stated physical guard preparations. No new claim
that this controller selects all exact targets is made.

## 4. Every seven-adic cofactor phase, with primitive ports

Specialize j=3m+2, m≥0. The initial and final counts in §3 are

    H(B(TT_j))=52·8^m+11 =7p_m,
    b_m=468·8^m+29       =7r_m,
    t_m=156·8^m+9,
    p_m=(52·8^m+11)/7,
    r_m=(468·8^m+29)/7=9p_m−10.                    (7)

Both divisions are integral. We have

    p_0=9,  p_(m+1)=8p_m−11,
    r_0=71, r_(m+1)=8r_m−29,
    gcd(b_m,t_m)=1,  t_m≡4 (mod7).                 (8)

The final graph has3m+11 vertices and contains its original strong
3m+8-vertex source. For any distinct nonnegative i,j,

    v_7(p_i−p_j)=v_7(i−j),
    v_7(r_i−r_j)=v_7(i−j).                          (9)

Proof: for d>0, binomial expansion first gives
`v_7(8^d−1)=1+v_7(d)`. One may prove this without an external theorem:
if7∤d the linear term in(1+7)^d has valuation1 and all higher terms
have valuation≥2; raising a number1+7^e u, 7∤u, to its seventh power
raises its difference-from1 valuation exactly once. Induct on v_7(d).
Now subtract (7). The factors52,468 and8 are units modulo7, and the
division by7 removes one valuation, proving(9).

Consequently, for every e≥1 the m in[0,7^e) give **every cofactor residue
modulo7^e exactly once**. The statement applies both before and after the
physical primitive bridge. In particular, for every a≥1 and every unit
c modulo7^e there is a displayed primitive-port source whose Hamilton
count is congruent to7^a c modulo7^(a+e), and hence has exact valuation a.
Choose r_m≡7^(a−1)c modulo7^(a+e−1).

This is the strongest coverage statement here. The counts grow
exponentially in m, so complete7-adic residue coverage is not exact
integer coverage. Already r_1=539 is composite; no primality hypothesis
is hidden. The bridge leaves the original numerical branch:
`13·2^j+11≡11 mod13`, whereas `117·2^j+29≡3 mod13`. Therefore (5) is
not silently recycled as a new B(TT_j) source. Reuse is through its actual
primitive balanced interface, with new target equations still required.

There is also a direct repeatable growth operation within the final
family, provided the module position is retained as a sidecar. Keep the
entire NNR exterior and append one new singleton at the end of M. It is
dominated by all old M vertices and A, dominates C,D,z, and is dominated
by each of the three NNR newborns. This completely specifies its one
newborn cut; every old arc and both last-pair names remain unchanged.
Ordinal slot multiplication gives (h,p)→(h,2p), and therefore

    (b,t) -> (2b−29h, 2t−9h).                        (9a)

The state stays balanced, with b−3t=2h. It preserves gcd(h,39p), since h
is odd. For h=1, three repetitions give
`(b,t)→(8b−203,8t−63)`, producing precisely the next seven-divisible
member in(7). This is a literal one-vertex extension and can be repeated
indefinitely; it is not an N/R command at the balanced root. Forgetting
the module position would lose the specified newborn cut. The constructor
still covers only its explicit arithmetic family.

### Unbounded nontransitive modules: squares place the clock in the right subgroup

The transitive module is not essential to the primitive bridge. Let
Q=(2,4,9,3), the strong order4 tournament, for which H(Q)=5 and P_Q(2)=32.
For r≥0 let M consist of2r ordinal copies of Q followed by TT_a. Choose
a_r∈{0,1,2} with a_r+2r≡2 modulo3, and put a=a_r+3m, m≥0. The module
is nonempty: when r=0 we have a_r=2. Ordinal composition multiplies both
H and P(2), slot by slot, so

    h=25^r,  p=2^(10r+a_r)8^m.

Equation(5) gives a strong graph of order8r+a_r+3m+9 and primitive state

    b=117·2^(10r+a_r)8^m+29·25^r,
    t=39·2^(10r+a_r)8^m+9·25^r.                     (10)

Primitivity is immediate from gcd(25^r,39·2^(10r+a_r)8^m)=1.
Furthermore p/h≡2^(2r+a_r)≡4 modulo7, so b≡0 modulo7. The same proof
as(9) shows that, for each fixed r, the cofactor b/7 is a seven-adic
isometry in m and covers every residue modulo every7-power.

This is an exact squares/doubles mechanism: (32/5)^2≡4 modulo7 lies in
the subgroup generated by2, while (32/5)≡5 does not. The even number of
Q blocks makes an available transitive-tail phase reach the required
ratio4. In this particular ordinal family an odd number of Q blocks
cannot reach7|b for any tail length, since its ratio is outside that
subgroup. This necessity concerns the specified constructor, not all
possible modules or all realizations of a target. For r>0 the module is
growing and nontransitive; no fixed-core completeness is assumed.

## 5. Reproduction and connection ledger

Run from the worktree root:

    python3 -B 05-knowledge/results/sixty_three_bridge_20261113/primitive/bridge.py

`receipt.json` retains the literal63 rows, all retained-root histories,
six module controls (TT1 throughTT4, C3 and the strong order4 tournament),
the three named induced embeddings for(9a), and the finite clock checks
through7^4. Endpoint subset DP independently
counts the full states; it does not call the formulas being checked.
The universal proofs are the run/cover decompositions and valuation
argument above. The inherited snapshot is pinned in the script.

| map | preserved | changed | necessary sidecar / boundary |
|---|---|---|---|
| M→B(M) | literal module and uniform incidences | SCC structure becomes strong | p=P_M(2), not h alone |
| TT_j→TT_(j+1) inside B | every old vertex and arc | H→2H−11 | a particular growing block, not arbitraryH |
| NNR | old graph induced; physical root names | strong count and primitive readiness | full(H,s,q,F) |
| appended module-tail vertex | old graph, root and last-child names | (b,t)→(2b−29h,2t−9h) | retained module location |
| m→m+1 in(7) | seven-adic phase clock | graph grows by3vertices | residue coverage, not exact target supply |

The least-used sidecar here is the ordered two-slot response: it explains
how a forbidden numerical factor can appear in a strong count without
being represented by any component. The canonical hostile is a module
with the same h but different p (TT1 versusTT2); a count-only compiler
would give the wrong output. The corrected near miss is to infer exact
target coverage, or infinitely many prime cofactors, from(9). Neither
follows.
