# Narrow independent audit of incoming conditional-energy bridge

Pinned source: `f128982d824452a919602687554ffea9ef33155c`, Eliott Cassidy,
`05-knowledge/results/odd_eighths_uniform_families_20261030/`.
Audit scope: `COLLATZ.md` E1–E3 and `PHASE_BRIDGE.md` scalar decoder.
Verdict: **PASS within the stated conditional-Haar scope; no changes required to Session115's claims.** No novelty or convergence certification is implied.

## Proof audit

- Each positive word of cost A is indexed by an odd source residue modulo 2^A; its length equals the number of ones in the A-step shortcut parity string. The final odd-endpoint bit is not part of the definition.
- On the original cylinder R=c mod 2^h, unequal costs above h are orthogonal, and they are orthogonal to costs at most h. At a common cost A=L+h, the nonzero Gram blocks have the same prefix residue modulo 2^L. Prefix affine offsets cancel in the squared modulus. The suffix multiplier is q^{-j}, with j prefix ones, and the number of such prefixes is binomial(L-1,j-1). Summation gives 3^{-j}. This proves E1 for every fixed h>=1, odd q>=3, odd c, and finitely supported complex coefficient sequence. Uniform convergence of the finite collection of V_m justifies removing the cost cutoff.
- E2's outermost pair j,j+h has exactly one window; its coefficient is 4^{-h}3^{-j} times a unit phase. The all-one suffix address -(q-2)^{-1} is an odd unit. Hence bandwidth h is sharp and the ordered complex Gram observer distinguishes every odd c modulo 2^h. This is not integer recovery from finite precision.
- For M>=h, the short-cost part of S_M is exactly -1/2. Full suffix windows vanish because the complete nontrivial character sum is zero; only the h boundary windows remain. This proves E3 with the stated M>=h and fixed-h quantifiers. The coefficient period divides the multiplicative order of q modulo 2^h. No estimate uniform in growing h or assertion about a specified integer follows.
- For q=3,h=3, the three partial-sum norms are 1, 2-chi_8(u)*sqrt(2), 1. This gives K_M(c)=1/4+(-1)^M(3sqrt(2)/64)chi_8(c), M>=3. Multiplying its centered value by (1024sqrt(2)/315)(-1)^M yields (32/105)chi_8(c), exactly recovering w_8(c)=1/c+1/(8-c).
- The LRC decoder averages K over actual maximizing points, with repetitions of phase labels retained, and multiplies by N_m/(7m). This agrees with the Session115 physical-point mass formula. It is not an unweighted average over four abstract phases. Clock parity is necessary, and the scalar decoder distinguishes reflection orbits, not all four oriented phases.

## Independent finite controls

`check_independent.py` constructs parity counts by forward shortcut iteration, independently of the source's printed polynomial. Literal truncated character Gram matrices are grouped by frequency modulo 1/2^h and compared entrywise with E1 using its finite geometric coefficient, not the infinite coefficient.

Normal execution: 93 universes, 7,914 matrix entries; h=1,...,5, q=3,5,7, every odd c, cost cutoff T=h+5. Maximum floating-point error 1.12e-16. An additional 48 controls cover all odd c modulo 8 and M=3,...,14 for the E3 coefficient and decoder. These numerical controls supplement the proof audit; they are not exact arithmetic or proof certificates. No claim to have replayed the incoming packet's full checker is made.

## Suggested synthesis note

A concurrent result, credited to the incoming odd-eighths packet at `f128982d824452a919602687554ffea9ef33155c`, supplies a genuine scalar link: the parity-corrected, centered conditional source energy recovers the LRC reflection-even weight w_8. This is compatible with Session115's obstruction to positive character-unitary transport: conditional-energy decoding is a different map, not that unitary transport. The bridge must retain clock parity and actual maximizer multiplicity; it does not identify the two dynamics, prove fixed-seed H1, or imply positive-integer Collatz convergence. Its conditional-Haar proofs and decoder were independently reviewed, with finite Gram controls recorded here.
