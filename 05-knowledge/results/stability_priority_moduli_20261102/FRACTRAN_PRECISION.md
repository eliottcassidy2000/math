# Priority, realization guards, and exact precision at every modulus

Status: **PROVED**, elementary arithmetic; executable finite controls in
`verify.py`. The ordered-rule semantics and universality are classical
[Conway, FRACTRAN (1987), especially sections 2 and 11](https://gwern.net/doc/cs/computable/1987-conway.pdf).
These are explicit precision results for that semantics, not a new universal
machine, a polynomial-time simulation, or a claim of literature priority.

## Correction and scope

A program is an ordered finite list of positive reduced fractions a_j/b_j.
On a positive integer N it uses the *first* j for which b_j divides N,
replaces N by a_j N/b_j, and halts if none applies. Labels j are retained.

THM-4852, file
`01-canon/theorems/THM-4852-primegame-phase-law-both-published-lists-agree-selector-cellular-automaton-in-base-29-primorial.md`,
overextended its tested squarefree setting: prime support <=P does not make
the base-P# units digit sufficient. For (1/4,3/2), integers 2 and 8 have
the same residue modulo 6, but select rule 2 and rule 1, respectively.
The still smaller program (1/4), in base 2, distinguishes 2 (halt) and 4.
PRIMEGAME's actual squarefree denominators and its cited phase law survive.

### P1. Minimal one-step selector modulus

Rule j is reachable exactly when no earlier denominator b_i divides b_j.
Indeed N=b_j realizes it in that case; otherwise an earlier rule always wins.
Let J be these reachable labels. The least positive period of the selector,
including halt, is

    L = lcm_(j in J) b_j,                                (P1)

with lcm(empty)=1. Sufficiency follows after removing unreachable rules.
If d is any period, the selector at b_j and b_j+d is j, so b_j divides d.
This proves minimality, not just a sufficient CRT chart.

For a base B whose prime support contains that of L, the least number of
low digits sufficient to select the rule is

    k = max_(p|L) ceil(v_p(L)/v_p(B)).

This is zero for L=1. Thus one base-P# digit is enough exactly when the
reachable denominator lcm is squarefree. In a general base it is enough
exactly when L divides B. No bound on prime support alone bounds k.

### P2. Minimal one-step successor precision

Fix m>=1. Observe (selected label, successor modulo m), returning only HALT
on a halted state. Its least period is

    Phi(m) = lcm_(j in J) [ b_j m / gcd(a_j,m) ].          (P2)

For necessity, compare b_j and b_j+d: selection forces b_j|d and equal
successor residues force m|a_j d/b_j. These are exactly the displayed
divisibilities. They also suffice, independently on every reachable branch.
If J is empty Phi(m)=1. We do not retain the old N modulo m on a halt.
Adding that observable would change the answer.

For nonempty J put beta_p=max_J v_p(b_j), alpha_p=min_J v_p(a_j).
Because every fraction is reduced, P2 has the explicit valuation rule

    v_p(Phi(m)) = v_p(m)+beta_p,             if beta_p>0;
                  max(0,v_p(m)-alpha_p),    if beta_p=0. (P3)

Thus iterating the sufficient one-step compiler r times adds r beta_p
or subtracts r alpha_p, truncated at zero. This is a statement about
integer precision, not field Frobenius. Iterated Phi need not be the minimal
precision of the actual r-step operation: intermediate states may be unreachable.

**Concurrent result recovered before publication:** incoming `e9dd62edb6`,
`../stability_priority_spectrum_20261102/FRACTRAN_CLOCKS.md`, independently
proves these one-step formulas. Its compact equivalent is
Phi(m)=D*m/gcd(m,G), with D=lcm_J b_j, G=gcd_J a_j and gcd(D,G)=1.
This is immediate from P3. That source also proves the contrasting terminating
family (1/2,2/3), whose exact rank v_2(N)+2v_3(N) decreases by one per step.
P4 below supplies the general guarded-horizon minimum; P5-P6 treat a returning
reciprocal family rather than claiming the one-step formula as unique here.

## P4. Exact finite-clock precision, with all priority guards

Fix a horizon r>=1. The observation O_(r,m)(N) is the actual rule-label word
of length at most r, followed by HALT if it halts before r legal steps,
or by the residue of the state after r steps modulo m otherwise. In
particular a state halting immediately after its r-th step is observed by
its residue, not by an extra unrequested halt test.

For a proposed legal prefix w, write its cumulative multiplier after t steps
in lowest terms as A_t/B_t, with A_0=B_0=1. At step t, applicability of
rule i at the preceding state is the divisibility condition

    d_(t,i) | N,  where d_(t,i)=B_(t-1)b_i/gcd(A_(t-1),b_i).

Let D_w be the lcm of the chosen d_(t,j_t). Let F_w contain every d_(t,i)
for i<j_t. Then the exact realization cell is

    C_w = { N>0 : D_w|N and f does not divide N for all f in F_w }. (P4a)

It is nonempty iff no f in F_w divides D_w; if nonempty, N=D_w is a witness.
This handles *all* positive integers, including arbitrarily large ones.
The formula follows by induction; the accumulated chosen divisibilities
ensure every intermediate state is an integer. In particular B_t|D_w.

The least positive period of O_(r,m) is exactly

    M_r(m) = lcm( { D_w : C_w nonempty, 1<=|w|<=r }
                   union
                   { B_r m/gcd(A_r,m) : C_w nonempty, |w|=r } ). (P4b)

**Necessity.** If d is a period, compare the witness D_w with D_w+d.
Their observed prefixes are the same, so D_w|d. For full-length w, equal
terminal residues additionally require m|A_r d/B_r. These give every factor.

**Sufficiency.** Induct on the common actual prefix of N and N+d. If one
selects j next, that child is realizable, and its D_(wj) divides d. Its chosen
rule is therefore legal at the other state too. If an earlier rule becomes
legal there, its own realizable child divisibility, also preserved by d,
would make it legal at the first state, a contradiction. The same argument
preserves halting. Terminal congruences are the second family of factors.

Consequently M_r(m) divides Phi^r(m), but equality is not automatic. The
compiler enumerates potentially exponentially many guarded words. Exactness
is not an efficiency claim and does not circumvent FRACTRAN universality.
It minimizes the residue modulus for the specified observation; arbitrary
non-residue encodings, unnamed paths, and different output languages are
different minimization problems.

## P5. An infinite reciprocal family and its parity law

Let a,b>1 be coprime and consider the ordered program (a/b,b/a). For every
r>=1 and m>=1,

    M_r(m) = a b^r m,             r odd;
             b^r lcm(a,m),        r even.                (P5)

Its one-step compiler gives Phi^r(m)=(ab)^r m. The discrepancy can therefore
grow exponentially; it is not merely one missing residue at a small size.

To prove P5, a state with maximal b-power t follows rule 1 exactly t times.
If t>=1 it then alternates 2,1, returning to the previous integer after
each pair. If t=0 and a|N it alternates 2,1 from the start; otherwise it
halts. Thus there are exactly r+1 realizable full-length words: 1^r,
1^t followed by the alternating suffix beginning with 2 for 1<=t<r,
and the alternating word beginning with 2. Their prefix divisibilities have
lcm a b^r. After r steps the multiplier is (a/b)^s with

    s in {0,2,...,r}       for r even,
    s in {-1,1,3,...,r}    for r odd.

All these exponents are realized. The largest exponent r forces b^r
times the b-primary part of m. For a-primary precision, the smallest
exponent forces the a-primary part of m for even r, and an additional
factor a for odd r; the selector already requires a. Primes dividing neither
a nor b retain exactly their m precision. CRT yields P5.

For (2/3,3/2), M_2(1)=18, whereas Phi^2(1)=36. If m=2, the even clock
requires 2*3^r, and the odd clock requires 4*3^r: one binary digit is
saved at the even clock by the *realized two-step return*. No bit is removed
without its return and phase certificate. For a single rule 1/p,
M_r(m)=m p^r; hence no fixed finite residue chart preserves every horizon.

## P6. Minimal residue modulus is not minimal information

For the same reciprocal family, suppose the terminal residue is unobserved
(m=1). At horizon h>=1 there are exactly **h+2** observations, even though
the minimal residue modulus is a*b^h. A smallest observation carrier is

    HALT, ALT2, 1,2,...,h,

where ALT2 means b does not divide N but a does, and t>=1 means the maximal
power b^t dividing N, capped at h. HALT means neither a nor b divides N.
These states give different label words: HALT, alternating from 2, or t
initial 1s followed by alternation from 2 (all 1s when t=h).
Thus every exact deterministic observer needs at least h+2 states, and
this one attains the bound. This is minimal among arbitrary finite encodings,
not only congruence classes, for this stated output language.

It is closed under **one legal execution step plus decrement of the clock**:
for h>=2, ALT2 goes to state 1, state 1 goes to ALT2, intermediate t>=2
goes to t-1, and the capped state h goes to the capped state h-1. HALT stops.
At horizon zero every nonhalted terminal state has the same m=1 observation.
It is not closed under indefinitely executing steps while keeping h fixed:
b^h and b^(h+1) share the capped state but their successors can be separated
by a fresh h-step query. The clock is part of the carrier's type.

So the required number of bits at a fixed horizon is ceil(log2(h+2)), not
ceil(log2(a*b^h)). The saving is justified by the restricted operation and
output language. It does not compress arbitrary FRACTRAN machines or supply
a finite all-time quotient. This gives a concrete analogue of the repository's
distinction between a sufficient modular chart and a minimal realized carrier.

## What transfers to tournaments, and what does not

The transfer is a precise method: retain legal-word guards, ordered ports or
rule labels, a common source, and the observation clock; minimize only after
restricting to realizable histories. P4 is the arithmetic analogue of keeping
joint path-component pairings rather than their endpoint marginals. P5 is a
literal return mechanism, not a numerical identification of a or b with an
automorphism order. A general execution can revisit rules and states, so it
is not a Hamilton path on the set of fractions. Completing its trace to a
tournament would discard priority and integrality unless they remain guards.
