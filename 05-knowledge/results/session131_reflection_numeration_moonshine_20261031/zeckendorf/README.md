# Reflection–numeration bridge

This packet proves actual maps and operation laws connecting the integral
reflection cycle lattice to CG and nonconsecutive Zeckendorf arithmetic.

## Main results

1. **Integral gluing is a binary carry extension.** For each twin pair,
   `s=2u+epsilon,d=2v+epsilon`; addition carries `epsilon*eta` into **both**
   integer coordinates. Reflection sends `v` to `−v−epsilon`.
2. **The actual order-five gluing code has 10 bits.** In the pinned 12-black-
   orbit order, `epsilon_0=0` and `epsilon_1+...+epsilon_5=0 mod 2`. All 1024
   codewords occur as integral cycles; boundary constraints remain explicit.
3. **Two exact minimal 17-state automata** recognize the paired CG and paired
   Zeckendorf words whose values can be integral sum/difference coordinates.
   The receipts include all reachable-state and distinguishing-suffix witnesses.
4. **A reversible two-state interweaving map** turns nonconsecutive even/odd
   Fibonacci strands into signed even-Fibonacci digits via
   `c_i=e_i+o_(i−1)−o_i`. It is value-preserving, not a count coincidence.
5. **The missing shift correction is exactly a boundary carry.** In positive
   CG addition, the number of rewrites touching weight 1 equals
   `P(a)+P(b)−P(a+b)`, hence is zero or one for two canonical summands,
   independently of rewrite strategy. This gives the sharp ±1 correction
   needed to shift reflection eigen coordinates integrally.
6. **Exactly three zero-appends restore universal compatibility.** In both
   CG and genuine Zeckendorf, equal-parity pairs either lift at every iterate,
   or exactly at shifts divisible by three. This follows from the exact
   order-three parity transfer matrix, with explicit sharp examples.

These are metagraph cycle and numeration theorems. No map here claims
square-sum Hamiltonian growth. The prior rooted-forest extension is retained
as a separately typed operation. Full proofs and scope boundaries are in
`PROOFS.md`; source commitments are in `SOURCE_PINS.json`.

## Reproduce

Python standard library only; no repository checkout or Internet required:

```sh
python3 check_bridge.py > replay.normal.json
python3 -O check_bridge.py > replay.optimized.json
cmp replay.normal.json replay.optimized.json
cmp replay.normal.json normal.json
```

The `--test-failing-gate` flag must exit 1 in both modes, demonstrating that
optimization does not remove checks. `normal.json` and `optimized.json`
are byte-identical verified receipts with full finite-state certificates.
`VERIFICATION.json` records final replay and portability checks.

Universe: integers 0–50,000; 512 random 512-bit integers; all signed-even
words through length 8; every CG digit-pair word through length 4 and every
Zeckendorf pair word through length 7; 1024 large pairs per DFA; all
nonnegative summand pairs in `[0,128]^2` under two rewrite strategies;
256 large multi-sums; signed carry pairs in `[-100,100]^2`; 1023 integral
cycle samples; and the complete 1024-word atlas gluing code. The receipts
state exact shift-period test boxes and counts.
