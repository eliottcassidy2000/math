# Every Mersenne target from a returning three-port construction

Every integer2^k-1 with k>=4 is realized by an explicit strong tournament with one Hamilton cycle and feedback vertex number one, on exactly k+2 vertices. This order is optimal within that feedback-one class. Together with the triangle for3, the construction realizes every Mersenne prime except7, without asserting that infinitely many such primes exist.

Two rays,16*4^t-1 and32*4^t-1, use the same literal fixed-three-port left/right period and return by H->4H+3. The full two-block formula proves that cancellation occurs exactly for a palindromic middle word with the stated balanced zero lengths. This is simultaneous growth of two base blocks, not one-root cloning.

- [All-order proof and scope](PROOF.md)
- [Standalone verifier](verify.py), with normal and optimized receipts
- [Actual graph rows and full prime certificates](witness.json) for31,127,511,8191 and the initial15 graph
- [Independent audit](audit/REPORT.md), checker and both receipts
- [Inherited physical compiler](dependencies/ROOTED.md)

The only7p member of these Mersenne targets is511=7*73. No unrestricted minimum-order or full-spectrum closure is claimed.

Run `python3 verify.py --witness generated-witness.json` and its `-O` form, then run `audit/check.py` in both modes. All scripts use only the Python standard library. The checker writes a witness only when explicitly given its path.
