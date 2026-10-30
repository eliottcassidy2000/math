# Historical correction: a missing lift-delay hypothesis

The historical THM-491, `01-canon/theorems/THM-491-lrc-ramification-tower-pisano-signature.md`, displayed pi(p^k)=p^(k-1)pi(p) generically. The concrete p=3 calculation is correct. The unrestricted prime statement needs a hypothesis.

For odd p put P=ord_p(Q), h=min entry valuations of Q^P-I. Then Q^P=I+p^h C with at least one entry of C a unit, and the odd-prime binomial argument yields ord_(p^k)(Q)=P p^max(0,k-h). A return modulo p^k must first be a return modulo p, so only multiples of P need be considered. This proves the formula without assuming that h=1. When h>=2 the usual first lift is delayed.

This is the classical order-lifting mechanism already used by the repository's October15 U1 theorem. The historical statement is unproved as an unrestricted assertion, not refuted by a claimed Wall–Sun–Sun example. No such example is produced. The repaired bound and prime3 tower are sufficient for the numerical example recorded there; no LRC inference follows from the equality of clock sizes.

Mechanism: a finite example with immediate lifting was promoted to a universal non-exceptionality assumption. Remedy: retain the matrix return depth h in the state, just as a carry or a source-relative guard cannot be dropped because it vanished in initial examples. The source file receives a narrow warning below its title, preserving all original text and unrelated status claims.
