# Independent self-converse spectrum audit

**PASS.** The frozen packet proves the exact self-converse Hamilton-count spectrum through 41 at arbitrary order, and its four extra omissions are genuine restrictions of the self-converse class. The SCC square-factor reduction and the use of the strong-order bound are both valid. No correction to the frozen source is required.

## Frozen source snapshot

All source files are in `session197-orbit-wildcard/`:

- `CONVERSE_SPECTRUM.md`: `9c5e01897617341ce6efd2441275f8c7deaeb106339e4bc59ea56503d255abd8`.
- `verify_converse_spectrum.py`: `41a539df5bdc324067b6795b89f2274e3cc74219c0313fe943e571361c6d0de5`.
- `converse-spectrum-normal.json`: `d34792589e6a8e508460a5a836f81bc61a969c53791fa788a1b8aa4dba9eea22`.
- `converse-spectrum-witness.json`: `8e535d5bfcf76f411fcf928b7b736c53803681073c9eae9ee6d8ecc3da95a41a`.

The exact all-order support through 41 is

    {1,3,5,9,11,13,15,17,25,27,29,31,33,37,41}.

The unrestricted witnesses of counts 19, 23, 35 and 39 therefore cannot be self-converse at any presentation or order. These are not new unrestricted spectrum gaps.

## Structural proof review

The ordered strong-component list of a tournament is unique. An anti-isomorphism reverses this list and sends each component to the converse of its paired component. Thus paired components have equal Hamilton counts, and the sole possible middle component is self-converse. The first half is itself an actual tournament of some Hamilton count `d`, so the full count is `d^2 h`. Conversely the actual ordered join `U -> S -> converse(U)` realizes every such product, with the middle part omitted or replaced by a singleton when `h=1`. This proves both inclusions in the claimed exact spectrum law.

In particular a squarefree self-converse count, including a prime or `7p` with `p!=7`, must already be the count of a strong self-converse component. Singleton padding and paired weak components cannot evade this restriction.

The external strong-order input was checked directly in [Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, The Electronic Journal of Combinatorics 13 (2006), N3, Theorem 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf), DOI `10.37236/1141`. It implies `5^(n-1)<=H^3` for a nontrivial strong tournament. The integer comparison `78125>68921` therefore bounds every strong component of count at most 41 by order seven. The proof applies this only to the middle strong component, never to an arbitrarily padded weak graph.

For `d^2 h<=41`, positive odd `d` can only be 1, 3 or 5. These values are all attainable. Squaring gives no further counts below the cutoff beyond the already listed 9, 27 and 25. The finite support and the exact SCC law thus establish the all-order claim without extrapolation.

## Independent computations

`check_converse_spectrum.py` uses no author checker. It generates the complete uncapped tournament isomorphism inventories through order seven, with counts `1,1,2,4,12,56,456`, and tests self-converseness by exact canonicalization of the converse. Hamilton counts use vertex-cover inclusion-exclusion over directed walks, independent of the source subset path DP and source literal enumeration.

For every one of the **114 nonempty self-converse classes** through order seven, the audit checks the complete ordered SCC pairing, component product, central self-converseness and exact `d^2 h` factorization. Its low-count strong support agrees with the frozen set. It also independently checks all emitted self-converse witness graphs and anti-involutions.

The four unrestricted escape graphs pass both walk inclusion-exclusion and literal Hamilton-word enumeration, and direct converse-isomorphism tests show that each is non-self-converse. The first three are strong; the count-39 example has component orders three and five, as claimed.

The independent checker agrees byte for byte in normal and optimized modes. A copy of the source verifier and its helper was also replayed in an audit-only directory, normally and under `-O`; its receipts and emitted witness bytes match the frozen source exactly. The source directory was not modified.

## Additional checked sharpening, separate from the frozen claim

The same official Busch theorem gives the sharper order-eight minimum 45. At order at least nine, already `5^8>43^3`. Hence strong count at most 43 still forces order at most seven. The independent inventory contains no self-converse count 43, and the same square closure therefore proves that the displayed support also remains exact through 43. A directly checked unrestricted witness is

    [124,121,74,32,12,20,56].

This observation was sent to the author for a separately frozen sharpening; it is not silently substituted into the through-41 source snapshot. Any later mod-four or larger-cutoff claim requires its own review.

## Reproduction

Run `python check_converse_spectrum.py` and its `-O` version in this audit directory. Dependencies are the independent `check_graph.py` and `check_converse_mass.py` utilities.

- Independent checker SHA-256: `119355edd5bbfd59d5ba7b11124cb61e6359233eef3e4af74f658927750b47e7`.
- Independent receipt SHA-256: `4473d2358315bcaf9a35cfef599f51f07933cb7845a02dcb5eb09a9189a77168`.

No literature-priority or unrestricted coverage claim is certified by this audit.
