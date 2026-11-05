# A sharp cycle-budget threshold at H=821 with two feedback vertices

Session197. This is an all-order exclusion in a specified graph class, paired with an explicit positive construction. It is not an unrestricted Hamilton-path obstruction or a general coverage theorem.

## 1. The result

Among strong tournaments with at most two feedback vertices and exactly821 Hamilton paths, the minimum number of Hamilton cycles is **four**.

The negative half is a complete capped-transfer closure at all word lengths: no fixed-two-port word has H=821 and 1<=C<=3. The positive half is the physical word003100001, with H=821 and C=4.

More generally the words

    0031 0^k 1, k>=1,

give actual tournaments of order k+7, with

    H=49*2^k+37, C=4,

and feedback number at most two. This is an infinite constructive family, not an assertion of a constant-response return from each output or of coverage of all odd integers.

## 2. The graph language and exact transfer

Start from two retained vertices a->b. A mask in {0,1,2,3} appends a vertex beating the selected subset of {a,b}, where bit0 selects a and bit1 selects b. Every other old vertex beats the new vertex. Always forget the newborn as a port and retain the same physical a,b.

Deleting those two ports leaves the transitive birth order. Conversely, any tournament with a feedback deletion set of size at most two is represented this way: extend the deletion set to two vertices if needed, order their arc as a->b, and append the transitive complement. All strong graphs relevant here have enough vertices. Thus this language is exactly the two-feedback graph class, up to relabelling.

Use the independently audited27-coordinate nonnegative path-cover transfer. The copied twoport_core.py and automaton.inc are the same physical transition system as the published lowcycle-two-port packet. The Hamilton-path readout is the sum of coordinates0,...,8. The cycle count after each mask is respectively

    0, x1+x2, x3+x5, x2+x5.

A cycle count is measured on the birth edge, not recovered from a possibly merged successor state. All operations and readouts are nonnegative integer sums.

## 3. Why both pruning rules are valid at every depth

Clip every coordinate at822. For any nonnegative integer matrix M,

    cap822(Mx)=cap822(M cap822(x)).

The same identity holds for the path and cycle readouts. Consequently the predicates H<=821 and C<=3 are exact under this cap, and any retained H<=821 or C<=3 readout has its exact value. Large hidden coordinates cannot cancel, because there is no signed readout.

It is safe to discard a prefix with H>821. Every Hamilton path in an induced subtournament extends injectively when a new tournament vertex is added: insert the vertex in a permitted gap, or at an endpoint, using a deterministic choice. A gap or endpoint always exists by the change from an in-neighbor to an out-neighbor along the path. Deleting the new vertex recovers the old path. Hence H cannot decrease under an induced extension.

It is also safe to discard a prefix with C>3 when seeking a final strong graph with C<=3. Such a prefix has a Hamilton cycle and is itself strong. Hamilton-cycle count is monotone under a strong induced extension. For completeness, let S be a proper strong induced subgraph of a strong tournament T. If some outside vertex has both incidence signs with S, insert it in each old Hamilton cycle, deterministically, to obtain a strong one-vertex extension and an injection of cycles. Otherwise each outside vertex uniformly dominates or loses to S. Strongness of T forces vertices b,a with S->b->a->S. Replace a fixed edge of each old cycle by that two-vertex ear. Deleting b,a recovers the cycle. Iteration reaches T and proves C(T)>=C(S).

The comparison is directly between the strong rejected prefix and the eventual strong target. It does not assert monotonicity through arbitrary weak intermediate graphs.

## 4. The complete closure certificate

support.cpp starts at the exact two-vertex cover vector. At each state it tries all four masks, rejects cycle readouts above3, applies the cap822 transfer, rejects path readouts above821, and stores every new capped vector. The breadth-first queue is exhausted; there is no word-length or graph-order cutoff. A separate emergency state limit returns a non-success/incomplete result if hit. The certified run does not hit it.

The frozen SUPPORT.json records:

- 374,692 distinct reached states, all374,692 processed;
- 441,723 admissible transition edges;
- 213,747 edges rejected by the cycle test;
- maximum shortest representative word length23;
- missing odd strong outputs from3 through821 exactly7,21,821.

Every other recorded value comes with an actual word and cycle count. The closure proves in particular that no strong two-feedback graph with H821 has at most three Hamilton cycles, at any order. The maximum representative length is an observed property of the completed closure, not a search cutoff. No published strong-order bound or external graph census is needed for this negative result.

## 5. A physical four-cycle escape and an infinite family

The seed word0031 has rows

    [14,44,56,48,35,1].

At retained root a=0 its exact profile is

    H=37, ending count e_a=5,
    reserved-start cover response R_a=49,
    starting count s_a=4, H(T-a)=5.

These small constants are checked by literal vertex permutations and ordered cuts in the verifier. A mask0 appends a universal sink. It preserves H and s_a, makes e_a zero, and doubles R_a: append the sink to either slot of the reserved-a ordered two-cover, allowing the other slot to be empty.

After k>=1 such sinks the relevant profile is therefore

    (H,e_a,R_a,s_a)=(37,0,49*2^k,4).

A final mask1 beats only a. Cutting a Hamilton path at this new vertex gives exactly H-e_a+R_a paths. Cutting a Hamilton cycle at it gives exactly s_a cycles; the source has more than one vertex, so no singleton exception occurs. This proves H=37+49*2^k and C=4, hence strongness. Both physical ports remain exposed throughout.

The condition k>=1 matters: the word00311 at k=0 has H=81, not the formula's86, because the seed ending response is5 rather than zero.

At k=4 the word003100001 has rows

    [974,2028,2040,2032,2019,1985,1920,1792,1536,1024,1].

Direct subset counting gives H=821,C=4. Deleting vertices0,1 leaves the transitive birth order. Together with Section4 this proves the sharp minimum in Section1.

## 6. Reproduction and exact dependencies

The transfer dependency is the published `session197/lowcycle-two-port/PROOF.md`, frozen SHA c4a16105ebc9f9aa9cb3cd61887ec4c437db58c361b2469af630899d536c5c73. Included twoport_core.py has SHA733155121ac7a5936a38a0ec1f3a602759f868c1e4e547fd4b6f89bced25f9cb, and automaton.inc has SHAde75a17dfb95fc1cbdfd789a969ea0ca03b1f573b9600e51a6bd2c9c4d914b10. These are unchanged audited files.

Run `python verify.py` and `python -O verify.py` after copying this packet. The verifier regenerates the C++ table from the combinatorial transfer, compiles optimized and unoptimized C++ closures, compares both complete results with SUPPORT.json, performs small independent Python closures, replays every stored physical witness without capping, and directly counts the seed and H821 graph. It writes VALIDATION.json beside itself. Only Python's standard library and an installed C++17 compiler are required. No binary, network access, census data, or dirty checkout is used.
