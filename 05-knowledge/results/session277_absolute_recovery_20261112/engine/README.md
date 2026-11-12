# Exact completion of an approximation, not optimistic rounding

**PROVED algorithm; finite runtime/control measurements below.**

Given order n and retained real-support cutoff s, `complete(n,s)`:

1. obtains the exact rational lower/tail interval from `../absolute`;
2. rounds its endpoints inward to an integer interval;
3. requests only b=bit_length(interval width) low bits of A000568(n);
4. computes those bits with an exact modular Burnside response recurrence;
5. returns the unique congruent integer in the interval.

No floating-point value chooses the answer or the precision. The test cache is
not used by the counting function. [Proof and provenance](PROOF.md).

The compatibility rule is important: the raw nonnegative-shift response can
be reduced modulo2^b at every step. The signed projective gauge from the exact
engine cannot be reused that way without carrying additional precision.
The note proves a general backward precision-debt recurrence for that sidecar.

## Results

For support cutoff3 (identity and one3-cycle real corrections):

| Order | Exact answer bits | Low residue bits requested | Bits not requested |
|---:|---:|---:|---:|
|40|621|521|100|
|64|1721|1544|177|
|80|2766|2538|228|
|100|4426|4132|294|
|120|6480|6122|358|
|150|10303|9846|457|
|200|18655|18034|621|
|250|29489|28703|786|

This answers the integer-accuracy demand without pretending the residual is
small in absolute terms. It also shows the cost: excellent relative forecasts
still leave most of the count's bits unresolved. Increasing s reduces b at the
cost of more real-sector work; the table in `results.json` compares0,3,6,9.

The raw Python modular recurrence reconstructed the exact count through200 in
the selected large-order controls. Its250 query was killed with exit137 before
producing a result. This is retained as a resource boundary, not hidden as a
successful computation. The C++/GMP backend is tested separately under an800MiB
address-space limit and90-second timeout per large query; `cpp_results.json`
records success at100,150,200,250. The250 reconstruction completed in14.42s
including launch, with441,859 peak response states. Its code is explicitly adapted
from Session270's raw counter, not claimed as a new exact recurrence.

At200 the bounded C++ residue query took1.72s including subprocess launch in
one run. Session270's historical full exact GMP query took1.096s median; these
are not matched workloads or machines and **no speed win is claimed**. The
hybrid's achievement is a rigorous certificate and a modular implementation,
not a claim that hybrid counting always beats direct counting.

Replacing a modular inverse per partition sector by one inverse of odd(n!)
reduces inverse overhead. A same-process one-run direct-partition comparison at
80 took0.836s termwise versus0.180s with one inverse. These timing observations
are descriptive, not statistically robust performance guarantees.

## Reproduce

From this directory, using Python3; C++ tests additionally need g++ and GMP:

```sh
python3 check.py
python3 response_check.py
python3 precision_check.py
python3 cpp_check.py
python3 dyadic_completion.py 100 --support 3
```

- `check.py`:987 controls of two partition modular implementations and exact
  interval reconstructions through40; selected reconstructions through100 at
  four support cutoffs; hostile intervals rejected.
- `response_check.py`:410 raw-response residue controls through40 and selected
  exact reconstructions through200. Its reproducible command intentionally no
  longer launches the previously killed unbounded Python250 case.
- `precision_check.py`:79,560 exact-divisible signed-shift chains, plus12
  examples where naïve fixed precision loses a bit on division.
- `cpp_check.py`:186 compiled residue controls through30, plus bounded
  large-order reconstructions and explicit resource receipts.

`results.json`, `response_results.json`, `precision_results.json` and
`cpp_results.json` record the scoped outputs. These are not an all-order runtime
benchmark. The theorem supplies exactness for each finite input; full output
requiresTheta(n²) bits and the response state count is not proved polynomial.

No Chan–Pak coverage statement follows from modular tournament-class counts.
The transferable principle is narrower and useful: preserve the information a
future operation needs, even when that information is invisible to the current
readout. Here that sidecar is precision debt, not a new Hamiltonian seed.
