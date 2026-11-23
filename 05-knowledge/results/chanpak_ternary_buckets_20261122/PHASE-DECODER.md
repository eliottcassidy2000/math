# An odd modulus chooses the source and the digit

**PROVED with an 18-phase certificate.** Take the five source families \(a2^s+b\), \(s\ge1\), with
\[
 (a,b)=(21,-13),(33,7),(55,7),(280,7),(1019,219).
\]
The first, fourth and fifth families have common centers of height \(3h\) with both physical corrections. On the second and third, tripling and the positive correction always exist; the negative correction on the specified centers occurs only at 139 and 117, respectively. [Single-reversal family](DIRECT-FORK.md) · [Bounded alignments](PHYSICAL-FORKS.md) · [Four-vertex families](ELEVEN-RAYS.md)

Modulo \(3591=27\cdot7\cdot19\), the 90 source residues (five rays times eighteen phases) are distinct. All 270 arithmetic child residues, including the digit choice \(-2,0,2\), are also distinct. The two sets are disjoint. Indeed the order of two modulo the three factors is respectively 18, 3 and 18. Thus testing eighteen consecutive exponents exhausts the residues. The [integer certificate](phase_decoder.py) builds the table, checks the three injectivity/disjointness counts and implements the target inverse.

Given an odd target \(T\), look up \(T\bmod3591\). Its entry gives a unique candidate \((a,b,s\bmod18,d)\). Demand that \((T-d-3b)/(3a)\) be an integer power of two at least two. If it equals \(2^s\), draw the marked tournament at precisely that address, respecting the two finite exceptions for the negative branch on the four-vertex construction. This chooses the source from the absolute target on the union of the five families. For example \(1703\) selects \((a,b,s,d)=(280,7,1,2)\).

**No arithmetic child from these five families starts any of the five families again.** In particular their finite and infinite physical successes cannot be concatenated into a ternary grammar merely by matching the next height. This conclusion uses the odd modulus. For *every* \(k\ge1\), the child \(3(33\cdot 2+7)=219\) equals the source \(1019\cdot2^k+219\) modulo \(2^k\). Arbitrarily deep binary agreement can therefore suggest an impossible continuation. The missing odd phase supplies a finite certificate of failure.

**Perfect-number boundary.** For all odd \(p\ge3\), none of the five designated central unit edges crosses \(N_p=2^{p-1}(2^p-1)\). This includes every perfect number with odd Mersenne exponent. Since \(N_p\equiv1\pmod3\), a center \(3h\) must equal \(N_p-1\). Look up this lower endpoint for each odd \(p\bmod18\). Only \(p\equiv7\pmod{18}\) survives, and its sole entry requires \(a=280,b=7,s\equiv13\pmod{18}\) and digit zero. In particular exponent 61 selects this specific candidate. But that center would make \(N_p=3(280\cdot2^s+7)+1\equiv2\pmod5\). Put \(A=2^{p-1}\). For odd \(p\), \(A\equiv1\) or \(-1\pmod5\), so \(N_p=2A^2-A\equiv1\) or \(3\pmod5\). This settles every exponent in these five families. Other arc choices, source families and the remaining Farey spines are outside this obstruction.

Run `python3 -B phase_decoder.py 1703` in this directory.
