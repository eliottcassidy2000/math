---
id: THM-2069
title: "The k-deletion CRT wheel is the low-weight evaluation code and matroid cogirth"
status: >
  RESERVED WITH A COMPLETE PROOF UNDER AUDIT. For integer rows spanning Z^r,
  reduction modulo p gives an injective evaluation code. A primitive covector
  fails some k-deletion gcd exactly when its codeword has weight at most k.
  Hence the exact local wheel factor is the initial code weight distribution,
  and its first failure radius is the column-matroid cogirth. The k=1 layer is
  THM-2062's coloop wheel; k=2 is exact cosimplicity. CRT multiplication, the
  Paley-e8 design application, and the [72,36,16] support-gate formulation are
  being independently checked.
source: codex-2026-07-21-code-cogirth-wheel
depends_on:
  - THM-2062
related:
  - THM-211
  - THM-480
  - THM-481
  - THM-487
  - HYP-2430
  - HYP-2764
---

# THM-2069 -- The deletion wheel is a low-weight code

This ID reserves the exact higher-deletion generalization of THM-2062. The
proof is under audit. No existence claim for an extremal `[72,36,16]` code is
made; its use here is to identify the first missing support layer, not to
realize it.
