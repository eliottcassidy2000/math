# Session132: executable FRACTRAN tournament counters and recursive graph structure

**PROVED program semantics and structural identities; FINITE-EXACT replays and graph certificates.** This session supplies actual fixed fraction lists, not merely an appeal to FRACTRAN universality. Classical counting and the inherited repaired pure-blue classification are credited. No fraction-minimality, fast binary-enumeration or moonshine-causation claim is made.

## 1. Three concrete programs

| Target | Ordered fractions | Input and readout |
|---|---:|---|
| A000568(n), ordinary unlabelled tournaments | 5,543 | Input `157*5^n`, n>=0. Successful halt marker37991; exponent of3 is the answer, error exponent at2 is0. Other work registers may remain. |
| Literal metagraph lines ell_n | 16 | Input `306*13^(n-3)`, n>=3. Halts at `59*2^ell_n*3^(n-1)`. |
| Pure-blue merged-node counts p_n | 22 | Fixed seed817086270. At each23-divisible checkpoint, exponent of19 is n and exponent of2 is p_n. |

[The compiled A000568 program and compact line program](fractran/THEORY.md) include the actual lists, a readable arithmetic source, compiler, control/register metadata and executable interpreter. The large artifact has3,979 control states and36 data registers. Its list is independent of n and contains no lookup table of sequence values.

[The pure-blue stream](pure_blue_fractran/README.md) prints its entire22-fraction list. It reaches its first25 checkpoints, n=0..24, in10,686 actual rational steps, with outputs

`1,1,1,2,1,3,2,4,3,6,4,9,6,13,9,19,13,28,19,41,28,60,41,88,60`.

Direct artifacts: [A000568 fraction list](fractran/a000568.fractions.00.txt), [its arithmetic source](fractran/a000568_source.py), and [the22-fraction list](pure_blue_fractran/pure_blue.fractran).

Output conventions are essential: these numbers are exponents at specified checkpoints or halts, not the physical integers and not arbitrary intermediate states.

From the fractran directory, the saved fraction files can be run directly:

```
python3 run_program.py a000568 3
python3 run_program.py literal16 6
```

They return2 and512 respectively. The implementation parses actual ordered fractions. Independent tests additionally scan the entire large list literally for n=0 and use a proven one-hot eligibility index for larger exact exponent runs.

## 2. What the large program computes, and why it halts

[The counting engine](counting/THEORY.md) writes an odd-part partition of n with multiplicities m_d and defines

`z=product d^m_d m_d!`,
`E=sum m_d(d-1)/2 + sum d*binom(m_d,2)`
`  + sum_(d<e) gcd(d,e)m_d m_e`.

A permutation with an even cycle cannot preserve a tournament. For an odd-cycle type there are E unordered-edge orbits, each with two orientations. Burnside therefore gives

`A_n=sum_(odd partitions) 2^E/z`.

The FRACTRAN source evaluates the integer numerator

`sum_(odd partitions) (n!/z)*2^E`

and divides by n!. Mixed-radix counters enumerate every allowed multiplicity vector once; an exact weight test selects partitions. Every loop has a finite counter or a decreasing Euclidean remainder. The emitted INC/DECJZ translation preserves exactly one control prime and implements first-eligible semantics. Self-transitions are routed through guards so cancelling a control factor cannot accidentally enable a fraction in another state.

This proves output and termination for **every allowed n>=0**. Finite replays supplement the proof: register execution through n=6 gives1,1,1,2,4,12,56; actual fraction-exponent execution runs through n=5, and producer literal-integer execution through n=3. The n=6 register run takes32,725,255 primitive steps. We do not claim to have run the FRACTRAN list through n=25.

The separate efficient arithmetic checker computes the counts through n=25 and compares all26 A000568 values against the corrected OEIS controls, including the historically misprinted n=24,25 values.

## 3. The recursive structure is richer than one scalar sequence

Define

`F_n(q,x)=sum_(odd partitions lambda) q^E(lambda) x^length(lambda)/z_lambda`.

The unified identities are

`A_n=F_n(2,1)`,
`S_(2m)=F_m(4,1)`, `S_(2m+1)=F_m(4,2)`,
`M_n=(A_n+S_n)/2`.

S counts self-converse classes and M converse-merged nodes. The packet supplies an integer cycle-insertion recurrence with retained state, signed vertex-extension formulas and strong-component recurrences. Counts through25 are independently checked by explicit permutation pair actions; external S/M comparisons have smaller, separately stated ranges.

[The fold theorem](converse_folding/THEORY.md) explains the half-order identities by actual inverse maps. Choose an anti-involution, orient each of its pairs by its internal winner, and retain the winner tournament together with the symmetric cross graph. Flipping winner arcs on that graph produces a second tournament. At odd order the fixed center contributes one binary vertex color.

This is an equivalence of **anti-involution-marked** tournament groupoids with ordered pairs of smaller tournaments, or colored pairs. Forgetting the involution gives a bijection of isomorphism classes because all anti-involutions are automorphism-conjugate. But pair automorphisms recover the centralizer, not the full tournament automorphism group: C3 has full group order3 and folded stabilizer order1.

The two smaller tournaments share a vertex set and a simultaneous relabeling. Their count is not A_m squared: S4=2 whereas A2²=1, and S6=12 whereas A3²=4. Relative alignment is the missing information. This is a real recursive carrier, but a scalar count still does not specify metagraph adjacency, colors or multiplicity.

## 4. Short programs follow genuine metagraph constructions

Literal line counts satisfy

`ell_n=2^(binom(n-1,2)-1)`,
`ell_(n+1)=2^(n-1)*ell_n`, for n>=3.

The recurrence comes from endpoint deletion of marked tilings: complement commutes with restriction, and each lower literal line has2^(n-1) lifts. The16-fraction program implements exactly this multiplication. For t=n-3, its exact runtime is

`6*ell_n - 5 + (7t²+29t)/2`.

It counts literal line instances, retaining loops and parallel multiplicity, not distinct endpoint pairs or merged nodes.

For pure-blue nodes, the inherited repaired strong-component classification gives two wrappers: add singleton components at both ends, or C3 components at both ends. The four centers are empty, K1, C3 and the regular cyclic T5. Therefore

`P(z)=(1+z+z³+z⁵)/(1-z²-z⁶)`,
`p_(n+6)=p_(n+4)+p_n`.

The22-fraction stream shifts this six-term window. The even/odd outputs interweave the Narayana counts c_m of words with weights1 and3: p_(2m)=c_m and p_(2m+1)=c_(m+2), with c_m=c_(m-1)+c_(m-3). For m>=2, these words are equivalently subsets of possible trimer starts1..m-2 separated by at least3; fill all unused cells with monomers. Their count is sum_k binom(m-2k,k). This is a precise nonconsecutive representation connection, with a different gap rule from ordinary Zeckendorf. Its relation to tournaments is this component construction, not an arbitrary tournament imposed on its register primes. Empty/one-/two-vertex counting conventions are distinguished from the literal-line geometry, which begins at n=3.

The emitted fractions have an integer grading G whose increment is exactly1 at every step. It gives the exact checkpoint time

`12p_n+10p_(n+1)+8p_(n+2)+6p_(n+3)`
`+18p_(n+4)+16p_(n+5)+8n-108`.

It also rules out nonempty periodic integer orbits from **any** seed. The valid stream has proved output/round termination; arbitrary multi-control seed behavior is not otherwise classified. [The recursion/cost analysis](RECURSION_BOUNDARIES.md) gives the grading and exponential magnitude envelope explicitly.

Ordinary, self-converse and merged counts grow faster than every exponential, so they cannot satisfy any fixed constant-coefficient linear recurrence like this pure-blue stream. A fixed FRACTRAN list can still compute them because its counters are unbounded. Exponent readout carries an unavoidable unary-output cost; the large program is not a fast binary counting algorithm.

## 5. The five prime groups and two different25-graph correspondences

[The graph atlas](prime_graphs/PROOFS.md) preserves all supplied lists. They contain25 entries but only24 distinct primes:71 is repeated and67 is absent from the first25 primes. No repair is silently made. The two rendered sheets, delivered separately, use the canonical first25 primes and mark these supplied-group memberships explicitly. Their edge colors encode square-root parity, not the metagraph blue/black reflection predicate.

**Prime-order reading: p_i labels G_(p_i).** Among the first25 primes:

- Disconnected:2,3,5,7,11,13, exactly the first supplied group.
- Connected but without a Hamiltonian path:19.
- Hamiltonian path but no cycle:17,23,29,31.
- Hamiltonian cycle: every prime from37 through97, including67.

All14 positive prime-cycle cases have explicit certificates. The first match is explained by the connectivity threshold14; it does not establish moonshine as its cause. The later supplied groups mix graph regimes.

**Rank-order reading: p_i labels G_i, i=1..25.** The nontrivial Hamiltonian-path orders are15,16,17,23,25, carrying prime labels47,53,59,83,97 and group labelsB,D,C,E,E. Their path counts modulo reversal are1,1,1,3,10. G1 has its trivial path; none of G1..G25 has a Hamiltonian cycle.

The duplicated/missing pair has the same coarse flags under either reading: G67 and G71 both have cycles, while G19 and G20 both lack a Hamiltonian path. Thus these flags alone cannot detect the indexing defect.

The atlas includes exact graphs and invariants for every order1..97. Its uniform structural theorems are

- G_n is connected iff n=1 or n>=14;
- G_n is2-vertex-connected iff n>=31;
- it is a forest through14, has girth5 at15..29, and girth3 from30;
- the arrival edge count is `floor(sqrt(2n-1))-floor(sqrt(n))`.

G31 supplies a sharp boundary: it is2-connected but has no Hamiltonian cycle. The arrival-root set is also too coarse to identify the supplied groups:41,43,47 all have the same set{7,8,9}, despite belonging to C,D,B.

**Value-indexed variant.** If one instead uses A000568(n) itself as the square-sum graph order, the [externally cited Gerbicz cycle-existence theorem](https://oeis.org/A090461), also retained in the [prior square-sum packet](../session131_reflection_numeration_moonshine_20261031/square_sum/PROOFS.md), already applies from n=6 onward. Indeed A6=56 and the general bound A_n>=2^binom(n,2)/n! exceeds32 at n=6 and then increases by the factor2^n/(n+1)>1. The earlier orders are1,1,1,2,4,12. This is a coarse inherited-existence corollary, not the interpretation used in either25-panel sheet and not a new construction of enormous graphs.

## 6. Verification and remaining boundaries

[Verification](VERIFICATION.md), [independent audit](independent/AUDIT.md) and `verify_session.py` bind the final programs, source, data and exact receipts. The independent fraction audit reads the emitted list rather than using the compiler to execute it. Other methods use signed pair actions, full small tournament censuses, direct graph-pair checks and explicit path/cycle certificates.

Two verified checkpoints preserve the distinction between current proved canon and superseded proof fragments. The historical affine-only argument for pure-blue atoms is not reused; its independently repaired classification is the cited input. Prime labels can be renamed consistently in FRACTRAN, so their numerical identities alone cannot supply additional tournament or moonshine structure.

The deliverable is an executable and recursive description with its lost information identified: counting state, relative pair alignment, control guards, path/line incidence, and graph indexing. It is not a universal bare-node recurrence or an optimized arithmetic-circuit result.
