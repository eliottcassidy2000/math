# Session55: recovered Parker / square-sum connections

## Scope and source pin

Repository inspected at commit `128996d7b87b8243c9644e20f2bdf2ff3ff481f1` (Session54 final), in `session50-publication/clean-repo`. All repository paths below are relative to that commit. No inherited result is claimed as new. No whole paper is included. Startup packet, AGENTS, START-HERE, CURRENT-FRONTIER, ACTIVE-GUARDRAILS, RESEARCH-PROTOCOL and CORE-PAPERS were consulted. Targeted sparse-file recovery used `git show HEAD:<path>`.

## Inheritance quartet

1. **Closest proved mechanism:** `05-knowledge/results/ap_square_pythagorean_20261019.md`, Theorem 1 (PROVED + audited): every positive increasing AP square-sum triangle is uniquely
   `(x,y,z)=2g²(C²−4AB,C²,C²+4AB)` for a primitive Pythagorean triple with size-ordered legs `A>B>0`, `C²>4AB`, and integer `g≥1`. Pair-sum roots are `2g(A−B),2gC,2g(A+B)`. The inverse passes from three AP square roots to a Pythagorean triple. This is a genuine reversible map, not a match of small constants.
2. **Canonical hostile:** `(A,B,C)=(4,3,5)` gives labels `(-46,50,146)`. All pair sums are squares, but positivity fails. Gaussian norm/factor counts without the angle/positivity sidecar overcount the intended graph objects.
3. **Corrected near miss:** `05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.md`, opening status and §4: the old claim that both edges at every degree-two vertex are forced was RETRACTED. A degree-two vertex may be an endpoint. At Q23 an actual path ending at 22 omits edge 14–22. The executable supplied here independently reproduces this hostile.
4. **Least-used relevant sidecar:** the actual signed AP difference, not merely its common norm or number of Gaussian representations. For path growth, its analogue is the pair of remaining named endpoints after a two-path join. Both are compatibility coordinates lost by scalar counts.

## Connection board and exact contracts

### A. Gaussian norm fibres → square APs → magic compatibility

`01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md`, §5: PROVED + audited fixed primitive hypotenuse fibres retain split-prime conjugate factor choices; the first affine V4 fibre is at C=1105. `THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md` retains content and section dependence under Gaussian multiplication. Do not treat primitive normalization as a free globally composable action.

Independent elementary reformulation: put `c=e²` and
`D(c)={d>0: c−d and c+d are positive integer squares}`.
A distinct positive 3×3 magic square of squares with center c exists iff there are `a>b>0`, `a≠2b`, with `{a,b,a−b,a+b}⊂D(c)`. Indeed the normal form is

    c+a     c−a−b   c+b
    c−a+b   c       c+a−b
    c−b     c+a+b   c−a

All eight line sums are 3c. The four opposite pairs give the four differences. Conversely any magic square has this normal form over the integers; rotations/reflections and sign choices permit the normalization a>b>0. Its nine entries are distinct exactly when the four nonzero magnitudes are distinct, equivalently a≠2b under this normalization. Positivity is included in D.

Map: integer representations `x²+y²=2e²` with `0<x<e<y` map bijectively to D(c) by `d=e²−x²`. Preserved predicate: each opposite pair consists of squares. Lost by taking only representation counts: additive identities among signed differences. Required sidecar: labelled D(c), not its size. Cheapest test: search a,b and exact memberships a±b. The supplied test gives e=1105 as a common-center hostile with many APs and no full compatible quadruple. This is a bounded certificate at that center, not an all-center impossibility result.

The AP-square-sum note §6 proves, using the CITED four-rational-squares-in-AP theorem, that neither endpoint can equal twice a square, whereas its center always does. Thus endpoint-to-center gluing of those graph triangles is impossible even after square scaling. Common-center gluing is possible in arbitrarily large windmills. The bare windmill with at least three triangles has no Hamiltonian path (delete its center); extra cross-edges matter. Neither gluing result solves the magic compatibility condition above.

### B. Finite quotients → integers: the precise obstruction

Elementary deduction: if an integer MSOS has minimum m and maximum M, reduction modulo every n>M−m preserves all distinct entries and all magic equations. Hence unbounded Parker quotient rings Z/nZ would obstruct an integer solution. The converse does not follow: compatible local points need not supply one bounded integral point.

Crucial type warning: infinitely many Parker extension fields is not this statement. Any map Z→F_(p^r) factors through F_p. Increasing r increases field size but does not improve injectivity on integer labels. In particular all characteristic-two fields have unavoidable repeated magic entries, yet this supplies no new integer obstruction beyond reduction modulo two. The missing sidecar is the characteristic/integer-image kernel, not field cardinality. Characteristic-two duplication is checked exhaustively over F2 in the test, and its linear identities hold in every characteristic-two field.

### Local 3/8/9 filters: useful, but not global closure

For a primitive integer MSOS (gcd of its nine roots equal to 1), all roots are odd and none is divisible by 3. Proof: modulo 4, each row sum is 3c. If c is even-square, every row has zero odd entries, so all roots are even; if c is odd-square, every row has three odd entries. Modulo 3, each row has either zero or three nonzero square residues; columns force the same choice for every row, and the all-zero choice contradicts primitivity. Consequently all entries are 1 modulo 24, and the normal-form differences a,b are divisible by 24. Modulo 9 their permitted unit-square residues are 1,4,7. This is a necessary sieve, not an obstruction: unit-square residues modulo 9 admit nonconstant local magic patterns, although repeated residues are unavoidable with nine positions. Keeping modular collisions is essential; integer distinctness does not require distinctness at every small prime.

### C. 15..25 square-sum transition → compatible endpoint surgery

Sources: `01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md`; `05-knowledge/results/square_sum_ports_20261023/README.md` (PROVED recurrences and guards; FINITE-EXACT census).

For Q_n on 1..n, i~j iff i≠j and i+j is square:
- Exact unoriented path counts for n=15..25: `1,1,1,0,0,0,0,0,3,0,10`.
- The first path, at 15, runs from 8 to 9. At 16 and 17 new vertices append at existing endpoints. At 18 the new leaf attaches to an interior position; path feasibility disappears.
- 15,23,25 are births from a previously nontraceable graph, hence the newcomer must be interior. New-neighbor count alone does not distinguish these from nonbirths.
- Birth degree is `floor(sqrt(2n−1))−floor(sqrt(n))`.
- Exact arrival formula: `h(H+v)=Σ_(a∈N) E_H(a)+Σ_({a,b}⊂N) C_H(a,b)`, where C counts spanning two-path covers with a,b exposed in different components. Retain the other endpoints for the next step.
- At Q25, edge 1–15 is unused in all ten paths, although it was present in the unique Q15 path. Edge roles are not permanent birth labels.

The same root-level source proves the cut budget `b_P(S)+e_P(S)≤|S|+1−c(G−S)`. At Q25 cuts alone leave three endpoint possibilities for a hypothetical selected 1–15; named-degree propagation and proper-cycle exclusion finish the proof. The map preserves actual spanning-path feasibility, not just degrees. This explains exactly why graph growth does not automatically interweave magic-square line obligations.

### D. All-order square-sum theorem: status is not conjectural

`05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.md` §§6–7,9 separates:
- **PROVED and audited in the repository:** Gerbicz's 25-fold construction sends an odd-length path from 1 to 3 to one on 25n+12 with the same endpoints, hence to a cycle; it gives lengths `(71·25^k−1)/2`.
- **CITED, not independently replayed in the repository:** all n≥25 have a path and all n≥32 a cycle. General construction uses 49-fold shifted blocks, a same-positional-parity pair of paths of lengths n,n+1, and a finite base 41..2032. The old audit read the method but did not fetch/replay its glue tables.

Fresh primary-author manuscript confirmation: Berlekamp–Guy, [Fibonacci Plays Billiards](https://arxiv.org/html/2002.03705), Appendix on square necklaces, explicitly records Gerbicz's all-order result and links code. This supports CITED status, not a claim that this session audited the full constructive proof.

MISTAKE-668 (`01-canon/MISTAKES.md`) additionally repairs the claim that only the welding square 25 explains the 3–4–5 connection: THM-4510 already handles the full tight square triple (9,16,25). Component sizes 3,4,5 alone remain numerology; doubling-square graphs also contain those sizes among more components.

## Primary-source boundary checks and wildcard

Cain [1908.03236v2](https://arxiv.org/html/1908.03236v2), §§3–9, is not blanket authority for every displayed guess. Its same-norm AP construction still needs imaginary fourth-power cancellation to close an hourglass, and a second cancellation for a full square. Conjecture 7.2 and the subsequent finite-field total omit the characteristic-two family proved Parker by Lemma 5.1; read any finite exceptional-field assertion with an explicit characteristic restriction. Its final implication concerns infinitely many Parker Z/nZ, not arbitrary fields. These checks are source reading plus the elementary typing argument above, not an audit of every theorem in the paper.

Auel–Singer [2609.09351v1](https://arxiv.org/html/2609.09351), September 2026: **PREPRINT, CITED / computation unaudited here**. Appendix A Proposition A.1 asserts that the MSOS surface contains no geometric projective lines. The proof uses 36 Schubert charts and Magma ideal-unit checks; the text was inspected, code not replayed. This is a concrete hostile to seeking a nonconstant affine-linear family of root coordinates. It does not exclude nonlinear rational curves or isolated rational points. Next decisive test: replay the chart certificates before importing the no-lines conclusion as a computational dependency. General integer existence is not resolved.

## Reproduction and limits

Run `python3 check_connections.py` and `python3 -O check_connections.py` from this directory. The outputs are identical. Plain unpruned vertex-order DFS independently recovers every path on Q15..Q25; no incorrect degree-two forcing is used. Other tests check normal-form line equations, collision guards, the common-center hostile, the positive AP map, and all 512 binary arrays. No exhaustive search outside these declared universes, no full field census, and no general integer nonexistence claim is made. No files outside this lane were intentionally edited and no Git publication was performed.
