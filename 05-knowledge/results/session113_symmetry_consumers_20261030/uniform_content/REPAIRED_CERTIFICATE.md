# One-digit phase repair: certificate protocol

**Status: PROVED at E=9 by an independently checked exact certificate.**
The auditor checked525,000,000 transition inequalities,126 entry inequalities,
valid-domain closure, finite values, maximum sign, and both infinite-letter
tails. The certificate hash is recorded in
`../audit_collatz/repaired9_independent.json`. All-E feasibility remains OPEN.

Use q=5^E,D=8q,N=5q−5. Store u=m modulo5 and z moduloq. Four ordinary modes
have u=1,2,3,4; a fifth mode means u=4 with a pending high-layer bonus.
The pending mode is allowed only when h=min(v5(z),E)>0. There are therefore
4q+q/5=(21/5)q valid states. Initial all-one prefixes are separately paid by
λ_E times their original clock. Their terminal u is either1 or4,z=0; allowing
both possibilities with zero remaining excess is a safe overestimate.

For each new epoch (letter a>1), let u',z' be the actual closed quotient update.
Charge h(z') if u'=1. If u'=4 and h>0, charge0 and enter the pending mode.
All other new-epoch charges are0. All charges below are in units of content;
subtract λ_E*a from each transition.

On a letter1 from pending mode, u'=1; charge h and clear pending. There was
no high-layer density charge before this first possible return.

On any other letter1, if the epoch is active and h>0, charge

    [u'=1] + λ_h − 1/2.

If h=0 or the epoch is inactive, charge0. At each epoch's first possible
return, the first term is its exact level1 return and h−1 pays the high-layer
envelope. Later ones receive exact level1 returns and the tail density.
Thus the total is exactly F_1 from THEORY.md §5, with truncated h≤E.
It remains a one-sided envelope, not an identity-return observer.

Integer form of the density is D(λ_h−1/2)=q−5q/5^h for h≥1. Entry and
transition inequalities use a vector V whose maximum is≤0, exactly as in
the older quotient. Each gross reward is≤E, including pending transitions,
so the all-letter tail is still DE−Na. Checking all letters through64 and

    maxV+DE−65N≤minV,  DE−65N≤minV

covers every omitted positive letter and every omitted entry. Excluded
pending-height0 states are outside the state domain, and the valid domain
must be checked to be closed by the independent auditor.

Producer command:

    g++ -O2 repaired.cpp -o /tmp/repaired
    /tmp/repaired E /tmp/repairedVE.bin

`repaired.cpp` is a producer only; its output and its own diagnostic checks
are not a substitute for an independent certificate auditor. The word
obstructions to the old envelope do not obstruct this repaired consumer.
An all-E certificate remains OPEN even if a specific E is verified.

## Verified ninth-layer theorem

For every finite positive valuation word w with original binary clock A,

    S_9(w) ≤ (488281/781250) A.

This is sharp:781250 consecutive ones have reward488281. The finite word
is realized by the positive source2^781251−1; no infinite all-ones positive
orbit is asserted. The content convention and all source/carry data remain
those of the inherited packet.

The checked certificate has minimum−199232020 and maximum−19531240=−2N.
Consequently every word containing at least one letter>1 satisfies the
stronger bound

    S_9(w) ≤ (488281/781250)(A−2).

The initial run of ones is paid by its original clock; the negative
certificate maximum supplies this two-clock gap for the rest of the word.
This is a fixed-precision, all-word theorem, not an all-layer theorem and
not a statement of Collatz convergence.

The valid state domain has8,203,125 states. The binary file stores9,765,625
slots including1,562,500 excluded pending-height0 slots. Its SHA256 is
`83e016db40ab86a1ae8d7e6aef605a4710c0d1279f23f1741e11f617a187670f`.
The auditor, unlike the producer, independently verifies every edge and
entry, state closure, and all-letter tail; corrupted-vector and `python -O`
controls ensure validation is not silently skipped.

Reproduce the independent audit after producing the vector:

    python ../audit_collatz/audit_repaired.py --help

The full command is:

    python ../audit_collatz/audit_repaired.py 9 /tmp/repairedV9.bin --output /tmp/repaired9_audit.json

The verified receipt is `../audit_collatz/repaired9_independent.json`.

The producer accepts only1≤E≤9 and explicitly serializes signed64 values
as eight little-endian two’s-complement bytes per slot, including the excluded
slots. The auditor reads exactly that format; binary files need not be published.

The two-clock gap is also sharp: append one letter2 to781250 consecutive
ones. The terminal multiplier becomes nonsquare, so that last letter adds
no five-adic content. The word has A=781252 and S9=488281=λ9(A−2).
The same positive source2^781251−1 realizes it: after the initial k=781250
ones the value is2*3^k−1, and the next numerator2(3^(k+1)−1) has exact
2-adic valuation2 because k+1 is odd. E10 remains unattempted and OPEN.
