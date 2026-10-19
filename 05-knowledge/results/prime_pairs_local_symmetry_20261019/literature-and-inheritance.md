# Literature, inheritance, and claim boundaries

Checked 19 October 2026. Sources are primary papers or author-maintained
material. Search results advertising solutions of major open problems were
not treated as established mathematics.

## Recent literature that affects the framing

1. Julia Stadlmann, *Bounded gaps between primes*, submitted 31 August 2026.
   https://arxiv.org/abs/2608.31126
   Full text: https://arxiv.org/html/2608.31126v1

   The preprint states H_1 <= 240, improving the familiar Polymath8b bound
   246, by combining Bombieri-Vinogradov with smooth-modulus distribution.
   Status in this packet: CITED PREPRINT RESULT. The abstract and paper were
   consulted; this session did not independently audit the proof. Do not
   describe 246 as an unqualified current record in the session summary.

2. D. H. J. Polymath, *Variants of the Selberg sieve, and bounded intervals
   containing many primes*, Research in the Mathematical Sciences 1 (2014).
   https://doi.org/10.1186/s40687-014-0012-7
   https://arxiv.org/abs/1407.4897

   The published paper gives the unconditional 246 bound and the conditional
   bound 6 under generalized Elliott-Halberstam. It also discusses the parity
   obstruction for that sieve setting. Its results concern bounded gaps,
   not a proof that the gap 2 recurs infinitely often.

3. Kevin Ford and James Maynard, *On the theory of prime producing sieves*.
   https://arxiv.org/abs/2407.14368
   Author version: https://www.ford126.web.illinois.edu/wwwpapers/prime-producing-sieves.pdf

   The framework formulates upper and lower prime bounds from specified
   Type I and Type II data, and constructions showing optimality in many
   parameter regions. It is the appropriate reference for carefully typed
   sieve-information claims. The elementary Boolean model in this packet
   is only an illustration of a particular information loss; it does not
   replace this theory or establish a barrier for all algorithms.

4. Lasse Grimmelt and Joni Teräväinen, *The Exceptional Set in Goldbach's
   Problem with two Chen Primes*, version 2, 29 July 2026.
   https://arxiv.org/abs/2508.16400

   The abstract gives representations by two Chen primes for numbers
   n=4 mod 6 outside a power-saving exceptional set. A Chen prime is a prime
   p such that p+2 has at most two prime factors. This is a relevant current
   additive result, with a deliberately weaker object and an exceptional
   set. It is not a theorem about two twin primes at every large target.

5. Harald A. Helfgott, *The ternary Goldbach problem*.
   https://arxiv.org/abs/1501.05438
   https://webusers.imj-prg.fr/~harald.helfgott/anglais/book.html

   Every odd integer greater than 5 is a sum of three primes. This solved
   ternary problem should be distinguished from the still-open binary
   problem and the stronger twin-center sum problem in this packet.

## Classical surrounding mechanisms

6. P. X. Gallagher, *On the distribution of primes in short intervals*,
   Mathematika 23 (1976).
   https://doi.org/10.1112/S0025579300016442

   Average singular series are classical. Our Proposition 4 supplies a
   self-contained specialization for the inherited one-parameter twin-center
   weight. It is not a claim to have introduced singular-series moments.

7. D. A. Goldston, Julian Ziegler Hunts, and Timothy Ngotiaoco,
   *The Tail of the Singular Series for the Prime Pair and Goldbach Problems*.
   https://arxiv.org/abs/1409.2151

   This studies weighted second moments of singular-series tails for the
   ordinary prime-pair and Goldbach problems. The shared local nature of
   sum and difference problems is classical; our finite-group formulation
   makes the exact reflection and its positivity loss explicit for twin
   centers and more general symmetric forbidden sets.

8. Michael Rubinstein and Peter Sarnak, *Chebyshev's Bias*, Experimental
   Mathematics 3 (1994).
   https://doi.org/10.1080/10586458.1994.10504289

   The classical analytic setting for one-point prime-race bias. This
   packet's elementary prime-power correction theorem does not independently
   reproduce the conditional distribution theory of those races.

A targeted search also found Tushar Pandey's author-maintained
*An intermediate conjecture between Goldbach and Dubner*:
https://github.com/pandey-tushar/twin-goldbach-note/blob/main/twin_goldbach_note_2.tex
The inspected portions formulate conditional twin-density assumptions,
ordinary pair singular-series estimates over twin members, and the familiar
11 twin-center holes. They are related context, not an imported dependency
or an independently audited research claim here. No global novelty is claimed
against this or other potentially overlapping work.

## Closest repository mechanisms, read without modifying the checkout

- `01-canon/theorems/THM-2447-twin-gap-local-prime-harmonics-and-detrended-spectrum.md`
  supplies the exact c_p(g) in {2,3,4} formula. It is explicitly inherited.
- `05-knowledge/hypotheses/HYP-9025-twin-gap-singular-series-partner-law.md`
  already defines w(g), distinguishes it from availability decay, and records
  a refuted geometric-depth model. This packet proves arithmetic moments of
  that same weight and extends its local role to sum targets.
- `05-knowledge/hypotheses/HYP-1994-twin-goldbach-necklace.md`
  supplies the K+K problem and finite holes. Its all-target conjecture remains
  open. The current exact 20,000-target check is a hostile control, not a
  claimed extension of the much larger existing census.
- `01-canon/theorems/THM-2443-twin-rank-parent-parity-margins-and-boundary-crossing-bridge.md`
  supplies ordered-parent parity, existing census margins, and a conditional
  implication. This packet does not promote those conditional conclusions.
- `05-knowledge/results/prime_power_bias_20261014d.md`
  gives the exact residue-vector prime-power inversion, its scale sidecar,
  and the one-point square correction. Proposition 8 extends the bookkeeping
  to additive convolutions and makes the two summand exponents explicit.
- `05-knowledge/results/prime_power_bias_collatz_collector_20261014.md`
  provides earlier context for the prime-power emphasis. No Collatz claim
  is needed in the new proofs.

The repository guardrails and research protocol were read. The assigned
focus was prime-pair and Goldbach-adjacent research. No LRC result is claimed,
and no theorem number or publication state was changed.

## Status and reproducibility ledger

PROVED: Propositions 1–9 in the main note, including Proposition 7 as an
implication with an open analytic premise; the F_13 hostile; the certified
product interval; the arithmetic moment formulas.

CITED: Prime number and Mertens estimates used only to express the growing
roughness-depth corollary in asymptotic notation; established and preprint
literature as individually labelled above.

FINITE-EXACT: Supporting checks across 365 symmetric subsets of cyclic groups
of orders 2 through 12; all 4,222 residue cases for primes 5 through 199;
all targets of the 385-wheel; 434 interval targets; exact rational prime-power
checks through 500; the r=4 Boolean-poset test; twin-center holes through 20,000.

NUMERICAL: Ordinary floating-point products and the 100,000-target arithmetic
moment experiment. The rigorous internal-surplus constant uses a separate
45-digit directed-rounding computation and a proved tail bound.

OPEN: All actual prime-pair/four-prime asymptotics, the uniform discrepancy
input in Proposition 7, and a transfer of the finite-sieve conditional
surplus to actual prime targets.
