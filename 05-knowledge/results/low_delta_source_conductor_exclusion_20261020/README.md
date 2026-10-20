# Low-delta plane polynomial algebra packet

**Status: proved; independently reviewed.**

Main result: a two-generator polynomial subalgebra of k[y] of finite codimension δ≤2 and fraction field k(y) contains a monic quadratic, even without immersion. Trace duality plus the residue at infinity gives a direct proof of symmetric degree gaps. For actual primitive monic Keller source algebras, specialization saturation globalizes a generic quadratic without assuming constant fiber delta.

The sharp δ=3 frontier is quadratic or cubic/quartic. Two exact immersed nonquadratic examples are included, one with an ordinary triple point and one with three ordinary nodes.

- Proof: `LOW_DELTA_QUADRATIC_AND_RELATIVE_GLOBALIZATION.md`
- Primary-source pins: `PRIMARY_SOURCES.md`
- Finite exact controls: `verify_low_delta.py`, `results.json`

Reproduce from this directory:

    python3 verify_low_delta.py > normal.json
    python3 -O verify_low_delta.py > optimized.json
    cmp normal.json optimized.json

The finite controls illustrate and audit formulas; they do not replace the proofs. Combined with the independently reviewed cylinder-generator obstruction at commit `9605e62f324416ddb1a448c2e923d12e7806166f`, a noninvertible actual monic/primitive Keller source has generic conductor degree at least six. At degree six its generic algebra is necessarily nonquadratic and admits cubic/quartic generators over C(x). No global cubic/quartic presentation or Jacobian-conjecture conclusion is claimed.

Producer controls: 33 exact checks, with normal and optimized output identical. The independent review and its exact computations are packaged under `review/`.
