# A generating function for every one-letter reclosure

Status: **PROVED** elementary identities below, **FINITE-EXACT** checker universe. No priority claim; classical divisor identities and cyclic subgroup arithmetic are applied to the inherited cycle-word formula.

## Objects and the distinction that matters

Fix an odd multiplier b>=3 and a nonempty marked positive valuation word u of length k and sum H. Let C=C_b(u). Append a new positive valuation v and *reclose the entire word as a periodic orbit*. This is not taking one step along the original orbit.

The new fixed point and its reduced denominator are

    x_v = N / Delta_v,
    N = b C + 2^H,
    Delta_v = 2^(H+v) - b^(k+1),
    q_v = |Delta_v| / gcd(N, Delta_v).

N is odd and coprime to b. Delta_v is never zero, since a positive power of 2 is even and cannot equal an odd power of b. Positive cycles are precisely the v with Delta_v>0. Their least such v is explicit and finite. The numerator N is independent of the appended valuation: the clock changes the closing equation, not that numerator.

The true odd-step map on a rational p/q with gcd(q,2b)=1 keeps reduced denominator q, since gcd(bp+q,q)=1 and removing powers of 2 cannot cancel an odd q. Thus the transitions 101 to35/175 below are surgeries on words, not transitions between rational denominator sectors of the native map.

## The all-clock cancellation theorem

Write g(v)=gcd(N,Delta_v). For each divisor d of N, let t_d=ord_d(2). For d=1 set t_d=1. Call d *accessible* if there is a solution

    2^v = b^(k+1) 2^(-H) (mod d).

If accessible, the solution is a unique residue class modulo t_d; let r_d be its least **positive** representative, 1<=r_d<=t_d. For d=1 take r_d=1. Then

    g(v) = sum_(accessible d|N) phi(d) 1_(v = r_d mod t_d),
    sum_(v>=1) g(v) z^v = sum_(accessible d|N) phi(d) z^r_d/(1-z^t_d).

Proof: the elementary identity gcd(N,m)=sum_(d|N,d|m) phi(d) is valid for any integer m, including negative m. Since gcd(2,d)=gcd(b,d)=1, divisibility of Delta_v is the stated cyclic-subgroup equation. A point of a cyclic orbit has one residue class of hitting times or no hitting times. Summing geometric series proves the generating function. All prime powers dividing N must be retained, not merely its prime support. QED.

This turns an infinite set of continuation tests into a finite set of phase-labelled cyclic subgroup tests. In particular g is periodic with a period L=ord_N(2), with the convention L=1 if N=1. The theorem supplies a period, not a claim that it is always minimal.

## Local possibility does not imply joint possibility

For coprime d1,d2, suppose both local equations are accessible with phases r1,r2 and periods t1,t2. Their product is accessible **if and only if**

    r1 = r2 (mod gcd(t1,t2)).

This is the generalized Chinese remainder theorem on the *exponent clocks*. The moduli d1,d2 being coprime is not enough: their multiplicative orders can share factors. Thus storing each prime's accessible/not-accessible bit loses the correlation required for cancellation of their product. This is a precise common-clock obstruction, not a claim about simultaneous independent random events.

## The 101 pair has two different exact response series

The inherited words u=(1,1,5) and u'=(1,3,3) have k=3,H=7, old denominator101, and C=19 and31. Their appended numerators are respectively185=5*37 and221=13*17. In both cases Delta_v=2^(7+v)-81, positive for every v>=1.

For N=185:

- d=5: t=4, r=1;
- d=37: t=36, r=25;
- d=185: t=36, r=25, since the two phases are compatible.

Consequently

    g_115(v) = 1 + 4 1_(v=1 mod4) + 180 1_(v=25 mod36),
    G_115(z) = z/(1-z) + 4z/(1-z^4) + 180z^25/(1-z^36).

In one36-step period, g is1 on27 positions,5 on8, and185 on1. Its exact mean is7.

For N=221:

- d=13: t=12, r=9;
- d=17: t=8, r=7;
- d=221: inaccessible, because9 and7 disagree modulo gcd(12,8)=4.

Consequently

    g_133(v) = 1 + 12 1_(v=9 mod12) + 16 1_(v=7 mod8),
    G_133(z) = z/(1-z) + 12z^9/(1-z^12) + 16z^7/(1-z^8).

In one24-step period, g is1 on19 positions,13 on2, and17 on3. Its exact mean is4. These means are outputs of these specific response systems, not a universal significance assigned to7 or4.

At v=1 the common closing denominator is175=7*5^2. The first numerator185 cancels one factor5, giving35; the second221 cancels none, giving175. At v=25 the first cancels all185. The second can never cancel its full221 under a one-letter extension, even though each of13 and17 is individually cancellable for infinitely many extensions. This is a concrete higher-order separation beyond the original one-bit example.

## Exact recurrence for the reduced denominators

The same construction gives a rational generating function for q_v. Let L be any common period of g, let g_r=g(r) for 1<=r<=L, and suppose Delta_v>0 for all v>=1 (as in the anchor pair). Then

    Q(z)=sum_(v>=1) q_v z^v
        =sum_(r=1)^L (z^r/g_r)
          [ 2^(H+r)/(1-2^L z^L) - b^(k+1)/(1-z^L) ].

In particular, for every v>=1,

    q_(v+2L) - (2^L+1) q_(v+L) + 2^L q_v = 0.

Proof: split v=r+jL and sum the exponential and constant terms. The recurrence follows by eliminating the constant forcing from successive L-step updates. If some initial Delta_v are negative, the signed-denominator sequence satisfies the displayed formula; the positive reduced-denominator sequence differs at finitely many positions, so its generating function is changed by a polynomial and its recurrence holds after that finite boundary. No minimal recurrence order is claimed.

The radius of Q is1/2: N is fixed and1<=g(v)<=N, so q_v has exponential rate2. This rationality applies to varying the last valuation of a **fixed marked prefix**. It is not a rational generating function for all Collatz cycles, an efficient unrestricted census, or a descent proof.

## Relation to tournament and shared-port carriers

The source is a marked word; the continuation is an appended valuation; the observable is cancellation of a divisor of the new closing denominator. The finite carrier retains divisor phase and order. Forgetting the phases retains marginal possibilities but loses whether two tests can fire at the same clock. For221, both prime tests are attainable but their conjunction is impossible.

This is the same abstract loss witnessed by tournament port responses that agree separately but disagree under a joint attachment. There is no asserted bijection from these integer divisors to tournament vertices. The transferable theorem is the fiber-product criterion: joint realizability is a compatible pullback over the shared clock, rather than the Cartesian product of marginal possibilities. The exact arithmetic instance is the generalized CRT condition above.

## Verification and provenance

Run `python -B spectra/check_spectra.py` and `python -B -O spectra/check_spectra.py` from the packet root. Outputs must agree. The standard-library checker independently computes direct gcds, divisor sums, reduced rational denominators, anchor histograms, prime-power cases, and the denominator recurrence. Its complete bounded word universe is multipliers3,5,7; lengths1..3; every letter1..4; appended valuations1..40. This is a finite replay of the proof, not evidence for an unbounded Collatz conclusion.

Inherited source: `05-knowledge/results/collatz_fermat_reframes_20261029/README.md`, sections C1-C3, at commit e871dd4513c7dcaa75cb29bb9847aafac35a795c. Its C3 supplies the101/35/175 pair. The present extension supplies the all-clock spectra and recurrence. Classical background on the2-adic conjugacy is Bernstein and Lagarias, [The3x+1 Conjugacy Map](https://websites.umich.edu/~lagarias/doc/bernstein.pdf); no new claim about the existence of that conjugacy is made here.
