# Multiplicative atoms and the 7p inverse slice

Session 197. **PROVED sharper finite-bank witness bound and FINITE-EXACT 7p slice.** This is a consequence of the target-bounded inverse theorem, not a theorem excluding any of these counts from arbitrary tournaments. The full `{7,21}` spectrum remains open.

## 1. Atomic targets need only one active factor

Use the grammar in `INVERSE_BANK.md`: K1, ordered joins, and a fixed finite bank of whole-module rooted templates, each adding at most M vertices. Let S be its set of attainable Hamilton counts. Call P>1 an atom of S when P cannot be written ab with a,b>1 and a,b in S. This is multiplicative indecomposability in S, not ordinary primality.

**Proposition.** Every attainable atom P has a witness in the same grammar with order at most floor(log_2 P)+M.

**Proof.** Apply the finished-outer-shell decomposition from `INVERSE_BANK.md`. All retained nonunit scalar factors are themselves attainable: active subexpressions are already graph expressions, while each passive factor c is represented by its actual template Q=C[K1]. Because P is an atom, there is exactly one retained nonunit factor. The active-input order bound gives at most floor(log_2 P)+M vertices, and a passive factor has a representative of at most M+1 vertices. No padding inside an active input is removed. QED.

If some prime ell is absent from S, then every target ell*p with p prime is an atom of S whenever it is attainable. Its only possible nontrivial factorization is ell times p, also when p=ell. The proposition therefore gives the prime-style cutoff for these composite targets too. This observation does not assume that ell*p itself is absent or present.

## 2. Seven is absent from the D/E/F/G grammar

Use the physical four-operator bank from the preceding packets; M=5. Its last-generator scalar forms, for positive odd h and integer u>=0, are

    D: 3h+6u,     5h+12u;
    E: 5h+18u,    11h+36u;
    F: 13h+36u,   25h+72u;
    G: 19h+60u,   39h+120u.

None equals seven. For an expression with H=7, primality forces every final ordered join to have a count-one side. Discarding only these finished join factors eventually reaches a last unary generator; all four bank generators are active. The scalar forms exclude it. Thus 7 is absent from this grammar at every order. This argument does not need any global theorem about tournament count seven.

Consequently every 7p, for p an odd prime, has complete witness cutoff

    floor(log_2(7p))+5.

For p<=127, the uniform cutoff is 14, since 7*127=889. Full-profile dynamic programming through degree 14, with H<=889, contains exactly 1,745 profiles. It gives the following complete result for all 30 such targets.

Attained:

    49, 91, 133, 329, 749, 763.

Omitted by the D/E/F/G grammar at every order:

    21, 35, 77, 119, 161, 203, 217, 259, 287, 301, 371,
    413, 427, 469, 497, 511, 553, 581, 623, 679, 707, 721, 791, 889.

The six positive targets have explicit constructions:

    49  = H(F(Tr_3)),
    91  = H(E(D(Tr_2))),
    133 = H(D(X^2 D(Tr_2))),
    329 = H(D(F(Tr_2))),
    749 = H(E(G(K1))),
    763 = H(E(X^3 D(Tr_2))).

Here X means an ordered-join singleton. These are actual graphs, not merely formal positive profiles. The verifier replays them into adjacency rows and independently counts their Hamilton paths by vertex-subset endpoint DP. In particular, missing the prime factor seven does not prevent realizing multiples such as 49 and 91.

## 3. What the slice rules out

The bank's 13 seed is not enough to attain 7*13=91 by multiplication, because seven is unavailable. It nevertheless attains 91 through the mixed active construction E(D(Tr_2)). Conversely, the other 24 failures cannot be remedied merely by running the same bank at a larger order: the atom cutoff already proves completeness. A new physical operator, a nonuniform port operation, or another grammar extension is needed to reach any missing value with this method.

This is deliberately a statement about a specified constructive language. It neither promotes arbitrary positive matrices to graph moves nor treats the table as evidence that the full tournament spectrum has additional exceptions.

## 4. Reproduction and dependencies

Run `python3 verify_seven_prime_bank.py` and `python3 -O verify_seven_prime_bank.py`. Both use explicit exception gates. The verifier imports the released full-profile DP and physical graph operations, checks the seven scalar obstruction, checks the 1,745-profile census, and independently replays all six displayed witnesses. The inherited DP correctness, order theorem, and H-monotonicity are proved in `INVERSE_BANK.md`; the operator identities are proved in `PROOF.md` and `UNIVERSAL_CONTEXTS.md`.
