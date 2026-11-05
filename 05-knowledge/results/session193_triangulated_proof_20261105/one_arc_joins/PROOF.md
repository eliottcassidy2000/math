# Session 193: one-apex joins as a complementary atom method

**Status:** PROVED elementary family reduction and bounds; FINITE-EXACT independent coverage through 2^20. The reported 2^32 cover remains SOURCE-REPORTED and UNAUDITED by this packet. Chan–Pak completeness remains OPEN.

## 1. Named provenance and exact scope

The main source is
`05-knowledge/results/old_mark_bound_lorentzian_gadgets_one_apex_join_coverage_20261104.md`
at commit `94805064e921a9649c61b621e784b6cb88fea172`, blob
`2e92d52a0056c8f7bf62b76931268199809beecb`. Its associated code is:

- `04-computation/experiments/chanpak_one_apex_join_cover_20261104.c`, blob `e61b21d1b5468f6fdef64a99ddd0598e9dd58ef6`;
- `04-computation/experiments/rooted_triple_bank_20261104.c`, blob `56221490c8e3e4b8272c04e959e0d843dca31960`;
- `04-computation/experiments/one_apex_triple_check_20261104.py`;
- `04-computation/experiments/chanpak_one_apex_join_cover_20261104.out`.

These named files, not a generic commit diff, were read. The main named summary was also read on main on 2026-11-05; it adds an exploratory residue profile absent from the pinned version. That profile does not alter the formulas below. The pre-audit commit was `61397dd494b4eaf168500bdab47a6acf9ad0aa38`. The reported 2^32 run and residue-growth observations are not independently reproduced here. The source itself explicitly disclaims an all-order proof.

The included `source_bank_dump.c` changes only the pinned bank program's output: it includes each graph's outgoing-neighbor masks and its mark. The independent Python code checks the resulting graphs and triples directly. We do not independently prove the source's isomorphism-class enumeration complete; completeness is unnecessary for a positive coverage certificate built from these actual graphs.

## 2. Definitions and join law

For a rooted tournament `(B,b)`, put `h=H(B)`, let `e` count Hamiltonian paths ending at b, and set

    R = sum_{S subset V(B)-{b}} H(B[S]) s_b(B[V(B)\S]),

where `H(empty)=1` and `s_b` counts Hamiltonian paths starting at b. Thus R counts ordered two-path covers whose second path is nonempty and starts at b; the first path may be empty.

Let T(v) be a transitive chain plus an apex, with the apex's outgoing arcs to the chain encoded by the bits of v from most to least significant. The empty word is used for v=0. The inherited exact rooted triple is

    (g(v), e(v), v+1),   e(v) = 1 if v is even, 0 otherwise,
    g(0)=1,
    g(2t+1)=g(t),
    g(2t)=g(t)+t+(t mod 2).

The even identity includes t=0 consistently.

For A=>B with its cross arc a->b reversed, paths not using b->a contribute `H(A)H(B)-e_a(A)s_b(B)`. Paths using the reversed arc split uniquely into a two-path cover of A starting at a and the converse two-path cover of B ending at b. Consequently the inherited join law is

    H = x_h y_h - x_e y_e + x_R y_R,

where x is the rooted triple of A and y that of B's converse.

The source families are:

    F1(u,v) = (u+1)(v+1)+g(u)g(v)-e(u)e(v),
    F2_y(v) = h g(v)-e e(v)+R(v+1),

with u,v>=0 and y from the order-at-most-seven bank. F0 consists of plain bank h counts and plain one-apex g counts. These formulas and the family are inherited, not claimed as new here.

## 3. Exact simplification: the one-apex part of F0 is redundant

**Proposition (PROVED).** Every value g(v) is already an F1 value.

For every t>=0,

    F1(0,t) = (t+1)+g(t)-e(t)
            = g(t)+t+(t mod 2)
            = g(2t).

Starting from v, repeatedly delete a trailing binary 1. The identity `g(2t+1)=g(t)` preserves its count. The process terminates at an even number 2t, including 0, so `g(v)=F1(0,t)`. This is a terminating explicit arithmetic map, not a density argument.

Only F0's plain **bank** counts need to remain beside F1 and F2. We do not claim the entire F0 is contained in F1. In particular, the source's finite F0 loop `v<=X` cannot create a coverage gap relative to adding all unrestricted one-apex counts, because F1 already supplies them.

## 4. Exact search bounds and atom filters

Since `g(v)>=1` and `h>=e>=0`, every F2 value satisfies `N>=R(v+1)`. Since `g(u)g(v)>=e(u)e(v)`, every F1 value satisfies `N>=(u+1)(v+1)`. Thus the finite loops used here are exhaustive for F1/F2 **values at most X**, relative to the supplied physical bank. The F1 symmetry permits u<=v.

These bounds also give short explicit certificates. Let `L(v)=ceil(log2(v+1))`, so T(v) has L(v)+1 vertices. An F1 count N has a witness of order at most `ceil(log2 N)+3`. An F2 count N with an m-vertex bank graph has order at most `m+1+ceil(log2(N/R))`. Leading zero bits merely add a universal source and need not be retained.

For a bank triple, put d=gcd(h,R). An exact branch filter is:

- v odd: d divides N;
- v even: N is congruent to -e modulo d.

Therefore an odd-v prime witness for p with d<p requires d=1, while an odd-v 7p witness with d<p requires d=1 or 7. These are necessary arithmetic filters, not a coverage theorem. In particular, the exploratory mod-3 bias does not by itself supply the missing all-order representation statement.

## 5. Independent bounded audit and concrete complementarity

`independent_check.py` uses a freshly written terminal-removal subset DP, rather than the source's forward-extension DP. It validates every supplied graph's tournament axioms and all **3,613 rooted rows**, giving **2,732 distinct physical triples on 532 graphs**. A separate permutation counter cross-checks the DP on **all 1,099 labelled tournaments of orders 1 through 5**. It checks 128 canonical one-apex words and 256 directly constructed F1 joins. No source numerical output is imported as an expected answer except the explicit bank cardinality and the advertised exclusion set used as controls.

The arithmetic cover uses the direct recurrence for g, not the source's optimized eight-bit suffix evaluation, and no source coverage code. Its universe is all odd N with `1<=N<=2^20=1,048,576`, all F1 parameters satisfying the exact product bound, and all F2 parameters satisfying `R(v+1)<=2^20` for the validated bank.

Results:

- The union of plain bank h, F1, and F2 misses exactly **7 and 21**.
- F1 alone misses **10,895** odds; adding plain bank h leaves **10,882** misses.
- There are **82,023** odd prime atoms p<=2^20 with p!=7. Of these, **2,553** require F2 relative to F1 plus the plain bank counts.
- There are **13,829** numbers 7p<=2^20 with p an odd prime >=5. Of these, **325** require F2 relative to the same baseline.

This is genuine finite complementarity on the prime and 7p lanes, rather than merely a second method reproducing exactly the first method's coverage. It is not an all-order complementarity theorem.

For example, the prime **163** has no F1 or plain-bank representation but has the physical bank triple `(h,e,R)=(15,2,40)` and v=2:

    163 = 15*g(2)-2*e(2)+40*3 = 45-2+120.

A six-vertex bank graph has outgoing masks `[2,4,33,7,15,27]`, rooted at zero-based vertex 3. Joining its converse to T(2) gives a nine-vertex witness. The composite bank count 15 alone loses the endpoint/cover information needed for this prime output; the rooted triple supplies it.

## 6. H127 / induced-H41 hostile control

The source bank gives **59** graph/mark/parameter rows with F2 value 127, corresponding to **22** distinct `(v,bank triple)` parameter signatures. The search over v is exhaustive by `R(v+1)<=127`. Every resulting graph was constructed, H=127 was independently checked, and all induced-subtournament counts were obtained in the same subset DP.

**None of these 59 outputs contains an induced subtournament with H=41.** Their orders are at most 15. This is a finite negative result for this supplied bank, not an exclusion for all possible F2 banks or all H127 tournaments. The absence remains true upon padding the one-apex word with leading zeroes: these only add universal sources and cannot create a new induced H count.

F1's feedback-vertex-number-at-most-two structural limitation is already recorded in the parent packet. F2 can leave that architecture in principle; this bounded test shows that this bank does not by itself exhibit the particular H127/H41 escape.

## 7. What remains open

The exact join laws and F0 simplification expose a clean arithmetic target, but do not establish that every prime atom or 7p atom belongs to the union. A proof still needs a quantified all-order representation argument, potentially combined with other direct constructions and valid smaller-parent operations. A finite bank, a growing tested interval, a decaying miss density, and an empirical mod-3 preference do not supply that argument.

The parent's binary transport proof identifies a separate obstruction to a naive bank induction: the even branch has a nonphysical triple and a subtractive carry. Nothing in this packet compiles that branch into a physical all-order construction. The finite 2^32 statement remains source-reported and unaudited here.

## Reproduction

From this directory run `gcc -O2 -o bank_dump source_bank_dump.c`, then `./bank_dump 7 > source_bank7.graphs.txt`, followed by `python3 independent_check.py 20` and `python3 probe_h127.py`. The generated graph bank is a source of actual witnesses; its completeness is not needed for the positive coverage certificate.
