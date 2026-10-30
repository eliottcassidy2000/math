# External source receipt and precise scope

**CITED paper / VERIFIED local Lean build / FINITE-EXACT identity replay.**
Read on 2026-10-30. No issue, message, or submission was sent externally.

## Identity and provenance

- Isaac Rudich and Louis-Martin Rousseau, *Deterministic Factoring in Time
  N^{o(1)}*, Hexagon 2610.00600v1, version dated October 29, 2026.
  [Landing page](https://hexagonmath.org/2610.00600v1),
  [PDF](https://hexagonmath.org/pdf/2610.00600v1),
  [Palomar record](https://palomar-registry.org/entry.html?id=PALOMAR-2026-10-28-000015&version=1).
- License: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), as
  stated on the landing page. The preserved PDF is unmodified. The
  commentary and deductions in this packet are separate from the source.
- Exact downloaded PDF: 427078 bytes, SHA256
  `e57f24481d2474a6f2c342d86312572f024710835ca583beb31e2bb1f9ceae40`.
  Six pages; text extracted, pages 3 and 4 also visually inspected.
  Web PDF retrieval failed; direct retrieval with a browser user-agent
  succeeded. The external Git repository's PDF has 429302 bytes and is
  not asserted to be identical. Hexagon and arXiv identifiers are distinct.
- [External repository](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring),
  pinned commit `c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.

The authors explicitly distinguish their statement audit from reading the
entire generated proof. We preserve that distinction rather than treating
the short statement as a short mathematical proof.

## Build and statement audit

In the clean external clone at the pinned commit:

```text
lake exe cache get
lake build
```

Both completed successfully; the build reports 1949 jobs. Lean is
`leanprover/lean4:v4.35.0`; the Mathlib manifest pins
`56a005855ea17e0b0958f73ec801ae140b0a1679`.
The Solution target printed:

```text
'Factoring.subexponentialRate' depends on axioms: [propext, Classical.choice, Quot.sound]
```

The theorem has the quantifier order

    exists M, for every k>=1, exists C,
       FactorsWithin M (n -> C*2^(n/k)).

The machine is a Mathlib `TM2ComputableInTime` machine, computing
`Nat.minFac` on `encodeNat` input and output. `FactorsWithin` requires
every stack alphabet to be finite and bounds `M.time` at every input
length. We read these definitions, compared their code with Challenge,
and matched the theorem signatures. The placeholder in Challenge is
intentional and lies outside the 145-file local Solution import closure.
The comment-stripped scan of that closure found no `sorry`, `admit`,
`axiom`, `native_decide`, `unsafe`, `implemented_by`, or `partial def`.
This scan is a syntactic cross-check; the Lean axiom report is the
relevant dependency evidence. The proof uses ordinary `decide` for the
finite crossing identity.

All 150 Lean/configuration source Git blobs are hashed in
[source_receipt.json](source_receipt.json). The full original build log's
hash and a concise [output receipt](lean_build_receipt.txt) are retained.
Reproduce collection after rebuilding with:

```text
python audit_external.py --repo <pinned-clone> --build-log <lake-build-output> --pdf <exact-Hexagon-v1.pdf>
```

No practical running-time measurement, independent review of all proof
steps, or proof of `SqrtLogLogConjecture` is supplied. The latter is
explicitly defined as a proposition rather than established as a theorem.

## Mechanisms actually reused

| Source location | Exact role here | Guard retained |
|---|---|---|
| PDF pp. 3-4; `Definitions.lean`, `Solution.lean` | Statement and witness-search target | One uniform machine; constants depend on k |
| PDF pp. 4-5; `MulCount.lean`, `Search.lean` | Prefix counts of divisors with their unique partners | Witness multiplicity; current prefix has a completion |
| `Transform.lean` | Dyadic roots in rings `Z/(2^m+1)` | Transform length must be invertible for inversion |
| `Cross.lean`, `crossing` | Four named ports and the exact signed `R+C` identity | XOR indices and negative branch values |
| `Moduli.lean` | Central-binomial residue bound | Actual nonnegative count height |
| `EngineGrid.lean`, `EngineGroup.lean` | Full recursive multiplier | Odd as well as dyadic content |
| `EngineMain.lean`, `deepMult_le`, `eval_mainE` | Integer bound and modular zero test | Multiplier bound enters the modulus budget |

The paper's §4 prose describes a power-of-two-scaled count. The generic
Lean implementation instead has

    mult(G,lambda)=(tB(lambda)+1)^nonRects(G)
                    *2^(tB(lambda)*card(rects(G))),
    tB(lambda)=(lambda+2)/3,
    card(rects(G))=binom(G,2)^2,
    nonRects(G)=G^4-card(rects(G)).

`deepMult` multiplies these factors over recursion levels. The engine
correctly bounds it by `2^mlogN` and chooses `M=n+mlogN+4`.
The generic parameters G=2, lambda=4 produce `3^15*4`; they disprove a
power-of-two claim about this generic function. They are not asserted to
occur in the top-level chosen schedule, and do not refute the built theorem.
Their role is to stop an unsupported power-two-only optimization.

C1 in [THEORY](THEORY.md) gives the general content-aware replacement
criterion. C2 gives a smaller bank for pure dyadic content. To claim an
implementation speedup one must bound the effective period after every
actual scheduled odd factor and account for its evaluation cost. That
obligation is open in this packet.
