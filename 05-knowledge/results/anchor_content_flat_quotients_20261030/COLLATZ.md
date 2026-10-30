# Changing-word content is a return cost, not free anchor reduction

**PROVED:** the algebraic renewal law and uniform fixed-prime logarithmic debt
bound below. **FINITE-EXACT certificate, with a proved telescoping implication:**
the sharp linear bound for the first four five-adic layers, for all word lengths
and all positive valuation letters. Universal endpoint control remains **OPEN**.

## Inheritance and the retained state

The closest proved mechanism is the matrix continuation rule in
[Session101, continuation state](../session101_denominator_connections_20261029/continuations/CONTINUATION_STATE_AND_SYMMETRY.md).
Its hostile pair is `(1,1,5)` versus `(1,3,3)`: the same old denominator does
not determine the reclosed denominator after appending a letter. The
[local carry gates](../collatz_local_prime_power_gates_20261020/LOCAL_CARRY_AND_PRIME_POWER_GATES.md)
already give cyclic-rotation content laws and the two-letter prime-power gate.
Neither is replaced by the prefix/suffix law here. The least-used sidecar in
this continuation is the *whole affine return*, including its translation.

For a nonempty positive valuation word `w=(a_1,...,a_d)`, write

    A=sum a_i,  P=3^d,  Q=2^A,
    F_w(x)=(P x+B)/Q,
    B=sum_(i=1)^d 3^(d-i) 2^(a_1+...+a_(i-1)),
    L=P-Q,  g(w)=gcd(|L|,B).

The initial state is `(P,Q,B)=(1,1,0)`; appending `a` sends it to
`(3P,2^a Q,3B+Q)`. For an expanding word the reduced negative fixed anchor is
`-b/D`, where `(b,D)=(B/g,L/g)`. We study the cancellation `g`; it is not the
whole anchor debt. Every `g(w)` is coprime to six. No word probability is
substituted for the weight of an actual integer source.

## 1. Content renewal, with exact quantifiers

For every two nonempty positive valuation words `u,v`, concatenated in that
order,

    gcd(g(u),g(uv)) = gcd(g(u),g(v)).                 (1)

Indeed

    L_uv=P_v L_u+Q_u L_v,
    B_uv=P_v B_u+Q_u B_v.

Modulo `g(u)` these are `Q_u L_v,Q_u B_v`, and `Q_u` is a unit. Taking their
common gcd with `g(u)` proves (1), including prime-power multiplicities.
In particular the contents of consecutive nonempty prefixes are coprime,
because a one-letter word has carry one. This does not say nonconsecutive
prefix contents are coprime: `(1,1)` and `(1,1,1,1)` have contents 5 and 65.

For any modulus `q>1` coprime to six,

    q | g(w)  iff  F_w is the identity affine map modulo q.             (2)

Consequently two visits of a prefix process to this condition force the
intervening word to be another full affine return. The clock is the actual
sum of its letters, not this sum reduced modulo a period.

If the word is realized by an actual odd source `n` with endpoint `y`, then

    g(w)=gcd(|P-Q|,|y-n|).                            (3)

This follows from `Q(y-n)=(P-Q)n+B` and `gcd(Q,P-Q)=1`. An endpoint residue
alone is insufficient: it must be paired with the multiplier clock. Equation
(3) preserves the marked source, not an independently rotated residue row.

## 2. A uniform bound for any fixed prime

Let `w_i` be the length-i prefix of an arbitrary positive valuation word,
let `A` be the total binary clock, and let `p>=5` be prime. Put

    S_p(w)=sum_i v_p(g(w_i)).

Let `c_e` be any proved lower bound for the binary-clock cost of a nonempty
affine return modulo `p^e`. Equation (1), including the identity starting
state, gives the return-packing inequality

    #{i:p^e divides g(w_i)} <= floor(A/c_e).           (4)

There is a universal elementary choice

    c_e=floor(log_3(p^e))+1.                          (5)

To see this, `B<=3^A-2^A<3^A`: after a letter the carry is `3B+2^A`, and
splitting that letter into ones can only increase this upper envelope.
Thus a positive carry divisible by `p^e` requires `3^A>p^e`.
Writing `L=floor(A log_p 3)` and `H_0=0`, summing (4) gives

    S_p(w) <= sum_(e=1)^L floor(A/c_e)
           <= (log 3/log p) A H_L.                   (6)

This is uniform in word length, letters, and source: a fixed-prime debt of
`O(A log(A+1))`. It is not an estimate summed over all primes.

For a positive Collatz orbit after `t` odd steps, `2^A y<=4^t n`, since
each step satisfies `3n_i+1<=4n_i`. Hence

    A <= 2t+floor(log_2 n).                           (7)

Equations (6)-(7) give an independent source-specific bound at every depth.
They survive averaging against any nonnegative *unchanged* source weights.
They do not establish a residual ROOT tail, a bound for the endpoint `y`,
or control of cancellation moving to previously untracked primes.

## 3. Exact five-adic return costs and the clocks 25, 50, 100

For two letters `(a,b)`, the inherited carry gate becomes

    5^e | g(a,b) iff 2^a=2^b=-3 mod 5^e.

Since 2 has order `4*5^(e-1)`, the phases for `e=1,2,3,4` are respectively
`1 mod4`, `17 mod20`, `57 mod100`, `357 mod500`. For example a two-letter
return modulo25 pays at least34 binary steps. Longer words can be cheaper.

The exact minimum costs over **all** nonempty positive valuation words are:

| modulus | minimum binary clock | one realizing positive source → endpoint |
|---|---:|---|
| 5 | 2 | 7 → 17 |
| 25 | 8 | 147 → 47 |
| 125 | 26 | 43447535 → 250823285 |
| 625 | 33 | 7354388327 → 995085827 |

The complete words and original source weights `4/((n+1)(n+3))` are in
[content_receipt.json](content_receipt.json). These are finite word segments,
not new claims about total stopping times.

Minimality has a small independent certificate. A cheaper word of length d
and clock A must have `3^d=2^A mod q`, with `1<=d<=A`. Below the stated
minima, q=5 and625 have no such pair; q=25 has only `(A,d)=(7,1)`; q=125 has
`(7,1),(14,2),(19,17),(21,3)`. Their respective composition counts are
`1,13,153,190`. None of these 357 words has zero carry modulo125. This
independent enumeration audits a Dijkstra discovery in the affine group;
it does not assume the discovery alphabet covers every possible letter.

For a *constant* letter a=1,2,3 the full affine orders modulo `5^e` are

    2*5^(e-1),  4*5^max(e-2,0),  5^e, respectively.  (8)

These follow from the fixed points -1,1,1/5 and the elementary lifting
identity `v_5(U^r-1)=v_5(U-1)+v_5(r)` for `U=1 mod5`. For a=2 the seed is
`4^4-3^4=175`, with valuation **two**. For a=3 the nonintegral fixed point
costs an extra factor5 beyond the linear-part order. Thus at moduli25,125,625
the three full orders are `(10,4,25)`, `(50,20,125)`, `(250,100,625)`.
An initial candidate incorrectly used valuation one for175; the exact control
caught and repaired it before publication. Equal-looking periods must retain
which affine observable they describe.

## 4. Sharp linear debt through four layers, for every changing word

For `E=1,2,3,4`, define

    S_E(w)=sum_i min(v_5(g(w_i)),E),
    lambda_E=(5/8)(1-5^(-E)).

The exact certificate proves, for all finite positive valuation words,

    S_E(w) <= lambda_E A.                             (9)

The coefficients are `1/2,3/5,31/50,78/125`. They are optimal over this word
universe, including words realized by positive integers. This is stronger
than adding the four separate return-packing bounds.

Here is the complete finite-to-infinite implication. Work in the actual
affine group modulo `q=5^E`, with state `(m,b)=(P/Q,B/Q)` and transition

    T_a(m,b)=(3m/2^a,(3b+1)/2^a) mod q.

Its reward is

    R_E(m,b)=sum_(e=1)^E [m=1 and b=0 mod5^e].

The supplied integer potential H, with `lambda_E=N/D` reduced, satisfies

    0<=H(s)<=H(identity),
    H(T_a s)>=H(s)+D R_E(T_a s)-N a                 (10)

for every state and every integer `a>=1`. Sum (10) along the actual word;
the potential difference is nonpositive because the initial state maximizes
H. This yields (9) with no additive constant and no reset of A.

| E | affine states | checked letters | max H=H(identity) | checked inequalities |
|---:|---:|---:|---:|---:|
| 1 | 20 | 1–6 | 1 | 120 |
| 2 | 500 | 1–8 | 7 | 4000 |
| 3 | 12500 | 1–10 | 119 | 125000 |
| 4 | 312500 | 1–12 | 422 | 3750000 |

The unbounded-letter tail is elementary: if checked letters end at c,
`N(c+1)>=max H+D E`. Every later letter has `D R_E-N a<=-max H`, so (10)
holds automatically. This is a clock-cost estimate, not a modular shortcut.

[five_potential.cpp](five_potential.cpp) constructs H by monotone integer
relaxation and verifies all inequalities. [audit_potential.py](audit_potential.py)
independently reconstructs the modular transitions in Python and verifies the
exported vectors; their byte hashes are in [potential_audit.json](potential_audit.json).
There is no floating-point optimization or presumed convergence in the proof:
only the final inequalities are used.

Optimality follows from `w=(1)^t` with `t=2*5^(E-1)`. Its prefix content is
`3^i-2^i`, whose valuation is zero at odd i and `1+v_5(i/2)` at even i.
Summing gives equality in (9). The finite word is realized, for example, by
the positive source `2^(t+1)-1`; no positive infinite all-ones orbit is claimed.

## 5. Precise next target and boundaries

**OPEN hypothesis:** (9) holds for every E. This would imply

    sum_i v_5(g(w_i)) <= (5/8) A

for every finite positive word. No converse from this weaker untruncated
bound to the sharp bounds(9) has been proved. The four audited levels do not
prove either universal assertion. The [next packet](../layer_lifts_cyclic_intersections_20261030/COLLATZ.md)
adds the fifth exact level, a phase-aware discounted epoch target, and a
marked-entry divisor theorem; the all-layer obligation remains open.
Reduction of affine states from `5^(E+1)` to `5^E` has fibres of size25.
Both coefficient lifts belong to the same word: independently mixing their
marginals is forbidden. Seek an inductive potential on this 25-point fibre,
including the carry in each transition. In unreduced units, the potential
equations scale by5 while the cost numerator gains4; this isolates a concrete
lift obligation rather than an extrapolation of a table.

More precisely take `D_E=8*5^(E-1)`, `N_E=5^E-1`, rescaling the checked
potentials to these units, and let pi be reduction modulo `5^E`. Retain the
nonnegative *old edge slack*

    sigma_E(s,a)=H_E(T_a s)-H_E(s)-D_E R_E(T_a s)+N_E a.

An ansatz `H_(E+1)=5 H_E o pi+K_E` satisfies the next Bellman inequalities
if and only if, on every lifted transition,

    K_E(T_a s)-K_E(s)
      >= 5D_E [T_a s=identity mod5^(E+1)]
         -4a-5 sigma_E(pi(s),a).                    (11)

The combined potential must still be maximized at identity. Formula(11)
retains paid slack rather than throwing it away: old inequalities that were
strict can fund a new high-valuation return. Each source fibre has25 possible
coefficient lifts, but the same letter transports their two coordinates
together. This is the proposed local induction problem; feasibility at
arbitrary E has not been proved.

Even this hypothesis controls only the prime5 part. A universal Collatz
argument still needs the other primes, actual endpoints, and the existing
fixed-source anchor-height bill. Constant words `(1)^t` and `(2)^t` have
unbounded raw content while their reduced fixed anchors stay -1 and1. Thus
raw content, reduced anchors, and positive-source descent are distinct targets.
