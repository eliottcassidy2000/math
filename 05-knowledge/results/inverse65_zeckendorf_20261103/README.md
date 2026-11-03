# A source-preserving inverse and quadratic lanes for the 65-symbol bank

Status: **PROVED** terminating inversion, append-only branching obstructions,
signed-charge recurrences, and the stated modular-clock laws; **FINITE-EXACT**
source/ear/automaton certificates. Full prime and 7p coverage is **OPEN**.
This packet does not prove that 7 and 21 are the only unrestricted gaps.

The bank is precisely the 65 marked atoms of
[the incoming seam packet](../golden_seam_recurrence_20261103/COVERAGE.md).
All five Q4 presentations are retained. Words here run in **processing order**,
the reverse of source-component order. A source word, physical port, and vertex
clock accompany every construction; equality of Hamilton-path counts never
authorizes replacing the supplied source.

Concurrent [golden-inverse work](../golden_inverse_defects_20261103/README.md)
independently establishes the realization inverse and strengthens the attachment
classification. Our [integration](INCOMING.md) credits that overlap and derives
the dual error/charge channels, with exact integral gluing indices2,5,13 and
two further actual source replacements giving3H-10 and3H-20. Their return
closure remains open.

The main results are:

- [Exact inverse](INVERSE.md): the positive norm `6h+9x+5y` grows by at least
  9/5 per symbol. Inversion is well founded and returns a literal source word,
  or certifies that a proposed state is absent. Scalar inversion is a complete
  finite decision procedure for this bank, with no claim that every required
  scalar succeeds. The certificate has logarithmically many symbols.
- [Quadratic lanes](CHARGES.md): exactly seven individual symbols have a
  nonzero rational linear charge preserved up to sign. Two explicit source
  families have counts `7p_k`, where `p_k` satisfies a forced quadratic
  recurrence. They include `7*19861` and `7*21001`; primality is separately
  certified. Their clocks modulo `13^e` have tail e and period `13^e`.
- [Actual ear lift and generating functions](EAR_GF.md): 28 of the 30 code-8
  ears lift into the bank on their actual six vertices. The two H=41 ears do
  not; adding one omitted marked atom repairs both. This is a source/order
  obstruction despite H=41 already being attainable elsewhere in the bank.
- [Zeckendorf audit](ZECKENDORF.md): a prefix-correct append-only decoder is
  impossible. A more permissive finite-state decoder correct only at the end
  would need a zero-input cycle of at least 12 states. This lower bound does
  not rule out such a decoder with more states or a different operation.

All six bank presentations of H=55 fail all 128 one-vertex attachments to 111;
an inherited positive 55-to-111 source outside this family still succeeds.
The 81/27 pair retains its identical cyclic triples and common H, and separates
under the same named-port attachment. See [audit](AUDIT.md).

Run from this directory:

```text
python -X utf8 inverse.py 133
python -X utf8 inverse.py 37 --x 18 --y 19
python -X utf8 verify.py
python -X utf8 branch_probe.py
python -X utf8 charges.py
python -X utf8 ear_bridge.py
python -X utf8 zero_cycle_prepare.py
```

The first command constructs an actual tournament. The second rejects a
positive cone state even though the scalar 37 is attainable. See
[recovered work](RECOVERY.md) and [next obligations](HANDOFF.md).
