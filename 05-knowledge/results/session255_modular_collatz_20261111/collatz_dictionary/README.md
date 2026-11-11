# Central modular tables are an exact signed Collatz edge atlas

Status: PROVED elementary dictionary and finite compiler theorem; FINITE-EXACT executable controls. Session 255, 2026-11-11. No convergence claim or external novelty claim.

## Inheritance and scope

The actual convention is the central 2-by-2 block of the multiplication table indexed by **1,...,m−1**, with entries reduced to **0,...,m−1**. This is not the entire multiplication table, nor an unspecified modular graph.

Closest proved mechanism: `05-knowledge/results/tiling_modular_atoms_20261003.md`, §3, gives the central blocks for m=4N±1. `01-canon/theorems/THM-4792-two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification.md` merges their rational ±1/4 interpretation with signed-sheet symmetries. The latter's modular sheet-alternation classification is a different question: it is not needed here. The old `05-knowledge/results/collatz_spectrum_equality_and_stitched_sheets_20261005.md`, §2, preserves the owner's mod9/mod11 examples; only the elementary central-block identity is inherited, not its numerical spectrum conclusions.

A direct earlier Collatz linkage is `05-knowledge/results/paired_prime_collatz_20261003.md`, §3 (present at commit `124b7ef29692d39b00b65079cd2bbf9468567da1`). Its shortcut map, denoted here H_s to avoid conflicting with our accelerated T_s below, is H_s(n)=n/2 for even n and (3n+s)/2 for odd n. It proves

    (H_+(2N−1),H_+(2N))=(3N−1,N),
    (H_−(2N),H_−(2N+1))=(N,3N+1).

These output pairs sum to the table moduli 4N−1 and 4N+1. Thus central blocks already encoded signed Collatz **shortcut outputs of adjacent source labels**. That note explicitly retains family/sign and full exponent sidecars, and warns against inferring the sign solely from a displayed 3N±1 expression. We credit that central Collatz linkage, not claim it anew. Our atlas instead reads an entry N as the **raw-edge source**, restricts N by parity, and decodes the even-source edge with the owner's division by six. The two conventions therefore differ in their source labels and map, not in their arithmetic: the block (N,3N+1) is a pair of minus-shortcut outputs at (2N,2N+1), yet for odd N it is also the raw plus edge N→3N+1. Similarly (N,3N−1) is a reordered pair of plus-shortcut outputs, yet for odd N a raw minus edge. The user-requested mod8 split in §2 makes that change of convention explicit.

The inverse-fibre mechanism x→4x+1 and its 3-adic odometer are inherited from `05-knowledge/results/arithmetic_braids_20260917_collatz.md`, §§1–3. New here is an explicit graph-type bridge for the owner's four mod8 classes and an arbitrary odd-modulus compiler for the signed inverse ladders. These are local constructions, not claims to have invented inverse Collatz formulas.

Canonical hostile: center 2→7 in the mod9 table is not a positive Collatz edge. Corrected near miss: treating that center pair as an ordinary edge without the owner's subtraction, division, and orientation changes. Least-used sidecar: the integer quotient and edge orientation, which a residue-class-only picture discards.

Concept board: table center / integer lift / orientation / signed sibling ladder / finite residue compiler. Anchor is the table dictionary; niche is signed edge coding; wildcard is the all-scale clock compression. The compiler links all five rather than asserting that the repeated picture implies descent.

## 1. Exact center calculation

For odd m=2k+1≥3, let d be the diagonal central entry and c the off-diagonal entry. Then

    [[d,c],[c,d]],   d+c=m,   4d≡1 (mod m).

Indeed k+1≡−k and 2k≡−1, so d≡k²≡1/4 and c≡−1/4. Least nonnegative representatives give, for every N≥1,

    m=4N+1: (d,c)=(3N+1,N),
    m=4N−1: (d,c)=(N,3N−1).

Composite moduli work unchanged. Squaring and multiplying the center labels are exactly the diagonal and off-diagonal operations; no prime or quadratic-residue hypothesis enters.

## 2. Four congruence classes, with the missing coordinate restored

Let s∈{+1,−1}, and define the ordinary signed Collatz map on positive integers by

    C_s(x)=3x+s if x odd;  C_s(x)=x/2 if x even.

The graph contains loops/cycles as the formula dictates. In particular the minus graph includes 1→2→1.

For all odd x≥1, the table of modulus

    m=4x+s

has the two central values x and 3x+s. Select them in that order to obtain exactly its odd-source edge. Therefore:

- m=8n+5 (n≥0): x=2n+1, edge x→3x+1.
- m=8n+3 (n≥0): x=2n+1, edge x→3x−1.

For every h≥1, the table m=8h+s has central values 2h and 6h+s. The owner's “reverse, remove/add one, divide by six” has the exact type

    (2h,6h+s) → (6h+s,2h) → ((6h+s−s)/6,2h) = (h,2h).

This is a **reverse-halving edge**, oriented h→2h. Reversing its arrow yields the forward even-source edge 2h→h. Hence:

- m=8h+1: subtract 1 before dividing by 6.
- m=8h−1: add 1 before dividing by 6.

The mod1 and mod7 families encode the *same* even-source graph on the two signed sheets. The former combines with mod5 for C_+; the latter combines with mod3 for C_−. The modulus 1 is excluded because its nonzero table is empty; this is why h starts at 1.

**Dictionary theorem (PROVED).** For either sign, the disjoint union of its odd-edge table family and its even-edge table family is in bijection with all edges of C_s, when the above decoding and orientation flags are retained. Each source has exactly one code: m=4x+s for odd x, m=4x+s for even x as well, but the even decoder uses h=x/2 and discards the center entry 3x+s only after recovering h. Thus a particularly compact encoding is simply

    code_s(x)=4x+s,   x≥1.

The parity of x, equivalently m mod8, selects the decoder. Distinct x have distinct codes. The inverse is x=(m−s)/4. This proves both surjectivity and injectivity; no counting or density heuristic is used.

“Half the edges” means the partition by source parity (equal natural density among source labels), not that each finite component or orbit contains equally many of both kinds.

Examples: m=5 gives 1→4; m=9 gives the decoded forward edge 2→1, not 2→7; m=7 gives minus even edge 2→1; m=3 gives minus odd edge 1→2. At m=13 the center 3,10 is the plus edge 3→10.

Connection contract: source = labeled central integer residues plus sign and mod8 family; target = directed Collatz edges; preserved = exact source/target incidence after the stated decoder; lost by an unmarked picture = orientation, distinction between a raw center entry and a decoded vertex, and integer height. The sign/family/quotient sidecar restores everything. This is not the assertion that all table cells are edges.

## 3. Infinite inverse edges compiled by one affine generator

Let T_s(x)=oddpart(3x+s) for odd positive x. Fix a positive odd target y with 3∤y. Let k0∈{1,2} be the unique exponent satisfying 2^k0 y≡s (mod3). Then every positive odd inverse source, and no others, is

    x_j=(2^k0 4^j y−s)/3,   k_j=k0+2j,   j≥0.       (1)

All x_j are positive, including y=1 on either sheet. The target numerator is exactly 2^k_j y, so the valuation is exactly k_j. Conversely reducing 2^k y=3x+s modulo3 forces the stated parity of k. Consequently

    x_(j+1)=4x_j+s,     m_j=4x_j+s,
    m_(j+1)=4m_j+s.                                  (2)

Both source and **its table modulus** run under the same signed affine rule. Under this rule the unhalved target multiplies by 4, so all these tables encode odd arrows landing on the same dyadic tower and the same accelerated target y. This is a commuting construction, not merely similarity of diagrams.

Multiples of 3 are excluded only as accelerated targets: they remain legitimate starting sources. There is no formula (1) for them because 3x+s is a 3-adic unit.

## 4. Finite compiler theorem (PROVED)

Fix any odd L≥1. For a ladder (1), set

    D=gcd(y,L),   Q=3L/D,   P=ord_Q(4).

(The modulus Q is at least 3, so the order is well-defined.) Then both sequences x_j mod L and m_j mod L have **exact least period P**, independent of the sign except through their initial phase. Their entire infinite integer lifts are represented by the finite list r=0,...,P−1 and an unbounded nonnegative quotient q:

    j=qP+r,
    x_(qP+r)=4^(qP)x_r+s(4^(qP)−1)/3,
    m_(qP+r)=4^(qP)m_r+s(4^(qP)−1)/3.              (3)

Thus a finite residue table plus one height coordinate encodes infinitely many actual edges without losing the valuation k0+2(qP+r).

*Proof of exact period.* Write A=2^k0 y=3x0+s. Then

    3(x_(j+t)−x_j)=A4^j(4^t−1).

As L is odd, 4 is invertible modulo3L; as 3∤y, gcd(A,3L)=gcd(y,L)=D. Therefore L divides the difference iff Q divides 4^t−1. This is independent of j and proves the least-period assertion. Since m_j=4x_j+s and 4 is invertible modL, m has the same least period. Formula (3) is geometric summation of (2). QED.

Useful specialization: if L=3^a H, gcd(H,6)=1, then D=gcd(y,H) and

    P=lcm(3^a, ord_(H/D)(4)),

where ord_1(4)=1. Indeed ord_(3^(a+1))(4)=3^a by the elementary valuation v3(4^t−1)=1+v3(t), and CRT takes the lcm. This splits the compiler into a ternary odometer and a coprime multiplication clock. The finite engine is exact at every chosen scale; it is not one universal bounded-state automaton for all L simultaneously.

Example: y=5, L=45 gives D=5,Q=27,P=9 on either sheet. For y coprime to45, Q=135 and P=18. A target's common factors with L collapse the coprime clock but never erase the ternary odometer.

## 5. What this construction enables, and what remains open

A finite residue predicate on an inverse ladder can now be decided by scanning P records, and every successful record yields the complete infinite family (3). Multiple predicates at odd moduli combine by their lcm. The valuation and height sidecars permit exact integer reconstruction and forbid confusing local coverage with actual descent.

A next constructive use is to attach a finite word of **actual** C_s edges to a compiler record and prove an affine inequality uniformly in q. This would turn a residue observation into an infinite certificate family. It still requires both path validity and the height inequality; no such universal cover is claimed here. The known minus cycles and arbitrary long plus growth prefixes remain necessary controls. Collatz convergence is OPEN.

## Reproduction

From the repository root:

    python3 05-knowledge/results/session255_modular_collatz_20261111/collatz_dictionary/check_dictionary.py

The script directly squares the central labels, independently decodes actual edges for all sources 1..2000 on both signs, checks exact inverse valuations, and checks the proved minimal periods for odd L≤99 and odd targets y≤31 not divisible by3. Controls include the mod9 non-edge and the identical even-source decoders. These are bounded checks supporting the algebraic proofs, not a large trajectory verification campaign.
