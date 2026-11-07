# Session214 one-arc escape controls audit

**Verdict: PASS (FINITE-EXACT constructions plus the stated minimum-edit consequence).** Each new restricted-class omission is attained by reversing one old chain-to-apex arc of an actual Session213 comb, retaining its original exit and adding a second.

Audited author bytes:

- ONE_FLIP_CONTROLS.md: `59d998fa03f04552ff90408e3059788152b0ca0a2c0ffa0506d4804fd2a2350d`
- one_flip_verify.py: `434b4921524b176809dd48709e0acae3d77bb04f009913389a2fc962ec734967`
- ONE_FLIP_RESULT.json: `30cac455d955772bd88daffc249c649ca489147ceb123234e9a8a70c325c67d5`
- ONE_FLIP_SOURCES.json: `812c42f41d59e50ea975dcfcec0b623d0ab579c599321a95b1b3546336dcba66`

`one_flip_audit.py` independently validates every parent edge against the comb definition, the transitive spine, the sole old exit, the binary word, tau and apex orientation. It checks that precisely the two row bits of one directed edge change. It then enumerates literal directed Hamiltonian words, importing neither the author's subset DP nor the run transfer:

- B^-(8,14,0), order12: 3187 paths; reverse x_0->a; child7477.
- B^+(9,179,1), order13: 24693 paths; reverse x_5->a; child27371.

Both child restricted apexes have exactly two chain exits. Since every unmodified original comb belongs to the independently audited class F, neither target is attainable with zero flips. Thus one is the exact minimum when one may choose any original comb and reverse old arcs on that same vertex set. This does not claim a minimum order or cover additions, deletions or arbitrary rebuilding.

Reproduce from the surgery-audit release directory:

    python3 one_flip_audit.py ../surgery/ONE_FLIP_RESULT.json > /tmp/session214-one-flip-audit.json
    cmp one_flip_audit.json /tmp/session214-one-flip-audit.json

These two constructions are positive controls, not an all-target recurrence. The independent literal counts also confirm the distinction between the12-vertex one-flip7477 witness and the separate11-vertex two-exit witness.
