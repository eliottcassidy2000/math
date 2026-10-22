# A positive survivor: density-one finite descent on the oddpart family

Session29 phase2, 2026-10-22. **PROVED scoped:** exact fixed-h endpoint-descent exponent densities, an explicit exponential tail bound, and density-one finite descent, including each of the six mod15/mod35 source registers. **Not proved:** convergence, a uniform first-descent bound, universal entry, or153N coverage. This is a transported Terras–Everett-type finite-prefix argument, not a new general Collatz density theorem or a literature-priority claim.

## 1. Input and counting measure

Use U(n)=(3n+1)/2^v2(3n+1) on odd positive n, with U(1)=1 for endpoint bookkeeping. Let u_k=(3^k+1)/4 for odd k>=1. All densities below are natural densities among odd positive exponents k, not among integer sources n.

The preceding `../REPORT.md`, §5 proves that any finite valuation word(a1,...,ah), A=sum ai, has exponent density2^(−A) on this family. It realizes one odd k class modulo2^(A+1). That exact cylinder theorem, including the terminal odd bit, is the only distributional input. The argument is deterministic. It makes no assumption of independence along one fixed integer orbit.

The parent proposed the endpoint-density argument; this packet proves and tests it, then adds a finite-register corollary.

## 2. A uniform carry bound makes the infinite word alphabet harmless

For a prescribed h-letter actual word write

U^h(n)=3^h n/2^A + beta_w.

The nonnegative carry satisfies the recurrence beta_new=(3 beta_old+1)/2^a, starting from0. Since a>=1,

0<beta_w <= B_h=(3/2)^h−1.

This follows by induction from B_(h+1)=(3B_h+1)/2. Equality is attained for the all-ones word. The bound is uniform over every h-letter word, including arbitrarily large valuations.

Define the exact integer clock threshold

s_h=min{s integer:2^s>3^h}.

In code this is `(3**h).bit_length()`, with no floating logarithm. The equality2^s=3^h is impossible for h>=1. Put

T_h = B_h/(1−3^h/2^s_h),

a positive rational. For every positive odd source n>T_h,

U^h(n)<n  iff  A>=s_h.

Indeed A<s_h makes the slope3^h/2^A>1 and the positive carry forces strict growth. If A>=s_h, the uniform bound gives U^h(n)<=3^h n/2^s_h+B_h<n. Thus possible failures of slope-only classification lie in the finite source interval n<=T_h. The increasing family u_k contributes only finitely many exceptional exponents. We do not need to enumerate or assume convergence for them.

This uniform threshold avoids an unjustified exchange over infinitely many possible valuation words.

## 3. Exact endpoint densities

There are binom(s−1,h−1) positive h-letter words with total clock s. Therefore the non-descent endpoint density is exactly

q_h = sum_(s=h)^(s_h−1) binom(s−1,h−1) 2^(−s),

and the endpoint-descent density is d_h=1−q_h. The sum is finite. Equivalently,

q_h = 2^(−(s_h−1)) sum_(j=h)^(s_h−1) binom(s_h−1,j).

The second identity views the h-th success time in fair independent bits and is independently tested as a rational binomial identity; it is not an additional stochastic hypothesis for Collatz.

For h=1,2,3,4, the exact descent densities are1/2,1/2,11/16,21/32. They are NOT the cumulative first-descent distribution: d4<d3. A concrete family witness is u9=4921: its first odd step is3691<4921, but its second is5537>4921. Endpoint descent and descent sometime before that endpoint are different events.

## 4. An explicit exponential tail, using only integer comparisons

For a formal independent geometric letter a with P(a=j)=2^(−j),

E[t^a]=t/(2−t),  0<t<1.

The preceding finite-cylinder law gives exactly the same distribution for each fixed finite word on our exponent family. Set t=3/4. If the total A lies in the bad range2^A<3^h, then3^5<2^8 gives5A<8h. Hence t^A>t^(8h/5), and the elementary exponential Markov estimate yields

q_h <= (3/5)^h (3/4)^(−8h/5)
     = (65536/84375)^(h/5).

In particular q_h→0 exponentially. Equivalently the exact rational certificate is

q_h^5 <= (65536/84375)^h.

The verification script checks this equivalent inequality, never a floating fifth root. This rational constant is convenient, not optimized.

## 5. Density-one finite descent without a uniform stopping depth

Let D be the odd exponents k for which some finite iterate U^j(u_k) is strictly below u_k. For every fixed h, the endpoint-descent set is contained in D, so the lower density of D is at least d_h. Since d_h→1, D has natural density1 among odd exponents.

There is no countable-union density interchange here. A lower-density inequality for each fixed h suffices. Nor does the theorem contradict the explicit all-ones cylinders from the first packet: for each fixed h a positive-density class fails to descend during those h steps, yet that density can tend to zero as h increases.

This is only FINITE DESCENT, not reaching1. The descended value generally leaves the family u_k, so the result cannot be reapplied inductively as if family membership were preserved. A density-zero exceptional exponent set is not empty merely because its density is zero. The k=1 source u1=1 has no strict descent and is one harmless exceptional exponent.

## 6. The result survives each mod15/mod35 register

The first packet proves that k=1,3,5,7,9,11 modulo12 gives the six source register pairs

(u_k mod15,u_k mod35)=(1,1),(7,7),(1,26),(7,22),(1,21),(7,12).

For every fixed word, the exponent modulus is a power of2. Fixing k mod12 imposes its compatible mod4 condition and a mod3 condition independent by CRT. In particular, a compatible word has density2^(1−A) within each fixed odd class modulo12, twice its unconditional density; an incompatible word has density0.

If k≡3 modulo4 then u_k≡3 modulo4 and the first valuation is1. If k≡1 modulo4 then u_k≡1 modulo4 and its first valuation is at least2. Thus for either type of register the bad endpoint density is at most2q_h. This tends to zero, proving density-one finite descent separately in each of the six registers. The finite exceptional source interval from §2 remains finite after restriction.

The actual source registers are preserved as information. Their presence does not supply native parent identities, terminals, clock signs, pair ownership, or the unused153N interaction budget. The first packet's at-most-two-sources-per-annulus obstruction still applies.

## 7. Literature inheritance and source access

Riho Terras, *A stopping time problem on the positive integers*, Acta Arithmetica30 (1976),241–252, DOI10.4064/aa-30-3-241-252. Publisher record checked: https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/30/3/101028/a-stopping-time-problem-on-the-positive-integers . C. J. Everett's1977 paper is cited in the first packet. Direct access to the old Terras PDF and Everett publisher article returned403 in this run; no unavailable text was treated as read.

Jeffrey C. Lagarias, *The3x+1 Problem: An Overview*, §3.3, explicitly attributes the uniform finite-parity mechanism to Terras and Everett and distinguishes it from long-orbit probabilistic modeling: https://arxiv.org/pdf/2111.02635 . The PDF was accessible and checked. The result here transports this classical mechanism through the proved oddpart exponent-residue bijection, with an explicit finite-exception threshold. It is not inferred by restricting a density-one set of integers to an arbitrary sparse family, which would be invalid.

## 8. Tests and remaining obligations

`python3 verify_density.py` and `python3 -O verify_density.py` produce byte-identical JSON. Checks include both exact tail formulas and the rational Chernoff inequality for h1..200;14970 literal endpoint/carry tests over odd k3..999 and h1..30; the source-dependent uniform threshold; the endpoint-versus-first-descent witness; and exact guarded-tail domination through h24. The sample saw no slope-classification exception, but the proof retains the finite-exception threshold rather than deleting it on that evidence.

No unbounded orbit is computed. The all-order statements come from the proofs above, not the finite tests. Next obligations remain pointwise coverage or a source-preserving descent/ownership mechanism. Density-one descent, even with all six modular guards, does not establish either.
