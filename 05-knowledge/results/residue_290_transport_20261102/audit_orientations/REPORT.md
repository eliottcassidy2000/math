# Independent audit of the signed290 wheel

**Verdict: PASS O1–O4 in their stated scopes.** The reviewed candidate is
`orientations/NOTE.md`, SHA-256
`729baac8ea4a5d7435d11e5e23b7cb20aceb3c4fea9ca17b5820c2a18d84380e`.
No mathematical or wording correction was requested. This is an independent
agent's proof review and independent implementation, with no producer imports.

## Proof review

- O1 correctly separates a half-orbit negation from negation exchanging two
  whole orbits. Unit cancellation gives equal lengths; the signed wrap has
  characteristic polynomial X^(d/2)+1, while the even wrap has minus1.
  The pullback convention for T does not alter these polynomials.
- O2 uses the odd lift147 of2 modulo290. Literal doubling is not a unit
  operation there. The order4 and order28 calculations prove the displayed
  four28/four14 orbit counts; these say nothing about increasing primes.
- O3's stronger error term follows by a summable O(log²n/n²) difference and
  Abel summation of a bounded zero-mean prefix. For every tail-start phase,
  the signed coefficients have partial sums in one of {0,1} or {0,-1}.
  Since f decreases beyond e, the claimed f(M+1) bound follows for M>=3.
  The positive paired-tail claim explicitly excludes the first n1 pair.
  Antisymmetry is sufficient, not necessary, for zero mean.
- O4's trace collision criterion follows from the quadratic equation for z.
  Differentiating the actual power-map identity is legitimate because
  z−z^(-1)!=0 for primitive odd N>1. It gives +2^k or−2^k exactly, rather
  than inferring a multiplier from the graph alone. The affine q7 and q9
  reversals are restricted three-point transports, not global conjugacies.
  The fixed eta condition gives−29/16 and multiplier35/8; period14 does not
  turn into period3 merely because29 occurs in both notations.

The current arithmetic-braids note was read at the explicit reversal and
q7/q9 sections; PRIME-MRB at its periodic-mean section. Source snapshots
are hash-bound, not presented as a re-audit of every claim in those files.
The auditor separately opened [DLMF27.8, equation6](https://dlmf.nist.gov/27.8#E6)
and [Conrad, Theorem2.5 and Example3.2](https://kconrad.math.uconn.edu/blurbs/galoistheory/cyclotomic.pdf)
on2026-11-02. They support the cited character and degree inputs. No external
priority or new cyclotomic irreducibility claim is made.

## Fresh independent finite controls

`check_independent.py` gives **174,568 active gates**, with identical normal
and optimized receipts. Its universe is explicit:

1. Every odd conductor3..201: all unit doubling orbits, all quotient residues,
   signed edge products, and both characteristic polynomials. The latter are
   assembled from actual quotient edge signs, not copied from the O1 formula.
2. Every residue modulo145: unique odd lift; all112 units have the literal
  147 multiplication transport modulo290.
3. At N7,9,145: cyclotomic polynomials obtained by exact monic division,
   trace polynomials reconstructed independently from Laurent coefficients,
   every proper pre-return iterate, and exact polynomial derivative products.
   This certifies the finite equations; the rational degree theorem stays cited.
4. The odd-character projector on all112 units, computed exactly in Q(zeta28)
   using CRT logs modulo5 and29, plus all three periodic coefficient means.
5. All290 tail-start phases, each with580 following coefficients, check the
   one-signed prefix condition; two independent rational decreasing-weight
   controls check the finite Abel bound. No floating approximation is used
   to certify the actual infinite MRB sum or its constant.
6. Exact q7/q9 affine polynomial identities, image multipliers1 and19, the
   rational−29/16 three-cycle, and its fixed-trace equation.

The producer's51,902-gate checker was also replayed normally and with-O;
that replay is not counted as additional independent checks. All four
historical source hashes were compared to their declared Git revision.

## Reproduction and exclusions

Run `python3 -I -B check_independent.py` and
`python3 -I -B -O check_independent.py` in this directory. The output must
match normal.json and optimized.json. The audit establishes neither a
prime-waiting distribution, a new closed MRB constant, a global affine
conjugacy, nor any famous open conjecture.
