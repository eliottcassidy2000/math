# Independent bounded review: one-swap trinomial gate

Date: 20 October 2026.

Verdict: **PASS under the zero-intercept convention used by the verifier.**
The note should explicitly define

    c_k = floor((k+1)p/q) - floor(kp/q).

Without that convention, “lower mechanical word” can also mean a word with
a nonzero intercept, for which the printed cut/intercept correspondence
needs an offset. This is a specification clarification, not a failed
arithmetic mechanism. It does not block the separately approved main
local-ideal theorem.

## Proof review

- The swapped cyclic pair occupies adjacent linear positions for all cuts
  except r=i+1 modulo q. At every other cut, the changed carry is exactly
  epsilon*2^j*3^s. Both bits' weights are fixed and the outside context
  contributes the displayed factor; the sign 10→01 is positive.
- For the displayed zero-intercept c, cutting at r changes the mechanical
  intercept to rp modulo q. Therefore r_t=t p^(-1) and
  t_star=(i+1)p modulo q are correct.
- The THM-4716 intercept-shift formula gives the displayed monomial M_t.
  Subtracting the two perturbation formulas gives M_t+E_r-E_r'.
- The intercept graph is a path with q vertices and q-1 edges. Omitting the
  exceptional vertex removes one edge when it is an endpoint and two
  otherwise. The retained count is exactly q-2 or q-3, respectively.
- Since gcd(q,p)=1, every binary word of that shape is primitive: a repetition
  count greater than 1 would divide both q and p. Thus distinct cuts give
  distinct words, and carry injectivity makes every retained trinomial
  nonzero.
- The main local-ideal theorem implies g divides each full difference.
  Hence g divides gcd(|D|, all retained trinomials). This is the stated
  divisibility upper bound. The deleted vertex and possible disconnectedness
  prevent claiming equality from these edges alone. Completing a connected
  graph restores the main theorem. No all-q proper-divisor result has been
  proved, and the note correctly marks that target OPEN.

The minimal witness itself is exact: c=0001001 is the zero-intercept
mechanical word of shape (7,2); swapping positions 3 and 4 gives w=0000101.
Its carry is 3*2^4+2^6=112, its clock is 2^7-3^2=119, and their gcd is 7.
This refutes one-swap freeness while still excluding the word as an integer
3x+1 cycle.

For clarity about the convention: if the initial lower mechanical word has
intercept t0, the correct correspondence is

    r_t = (t-t0)p^(-1) modulo q,
    t_star = t0+(i+1)p modulo q.

The unshifted version can fail at q=4,p=1,t0=1: c=0010, so
B(c)-B(rotation_1(c))=4-2=2, whereas the unshifted M_0 formula gives 4.

## Independent exact controls

`independent_one_swap.py` imports no producer/repository mathematics and
uses the positional carry sum. Normal and optimized runs give identical
receipts.

- All allowed swaps of zero-intercept coprime mechanical words for
  4<=q<=14: 308 cases, 2,548 retained edges, and 18,200 explicit checks.
- A separate inexpensive complete occurrence probe for 3<=q<=80 reproduces
  the producer's counts: 52,792 occurrences and 3,821 with nonunit content.
- It independently reproduces the first hostile above and the maximum
  observed g=11,243 at (q,p,i)=(66,35,11), including the exact word and carry.
- There are no observed positive integer candidates with D>1 in that finite
  universe. This is not an all-q family exclusion.

The producer's normal and optimized receipts were compared byte for byte;
the larger q<=40 trinomial run was inspected but not repeated.

Reproduction:

    python3 session11-collatz-review/independent_one_swap.py
    python3 -O session11-collatz-review/independent_one_swap.py

Both outputs must match `independent-one-swap.json`. SHA-256 source and
artifact pins are recorded separately in `ONE_SWAP_AUDIT_PINS.json`.
