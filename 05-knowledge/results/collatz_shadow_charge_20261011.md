# Fixed-source terminal switching: a charge staircase and a carry floor

2026-10-11. **PROVED, scoped and self-audited:** the charge-only bounds,
exact cylinder supremum, and certificate interface below. **FINITE-EXACT:**
the stated controls. **OPEN:** bounding the staircase at every fixed source,
independent universal ROOT positivity, and Collatz. No independent peer audit
or formalization is claimed.

## Inheritance and the changed target

The nearest mechanism is [Q1--Q3 of the quarter-price packet](collatz_quarter_price_source_20261011.md).
The [previous completed family](collatz_terminal_completion_20261011.md)
selects its sources; it does not cover a prescribed integer by changing its
terminal parameter. The [terminal lifts](collatz_terminal_lifts_20261007.md)
already complete arbitrary supplied ROOT cylinders, and the
[fixed-head S-unit obstruction](collatz_terminal_head_finiteness_20261007f.md)
already excludes a finite fixed-middle library on all Mersenne inputs.
We do not re-prove weaker versions of these statements.

The hostile is a succession of genuine rooted neighbours whose expense grows
with the observed precision. The corrected near miss is replacing the marked
source with one of those neighbours. The least-used coordinate is the number
of **odd valuations**, rather than total valuation, word length, or endpoint
height. The session board is marked source / odd charge / terminal phase /
coordinate ring / boundary relation / arithmetic scale bill. The orthogonal
lanes are [rhombus carriers](rhombus_overlap_elliptic_20261011.md) and
[Pell and odd-zeta scales](pell_zeta_scale_20261011.md).

Use the odd map U(n)=oddpart(3n+1), with strict ROOT at the first 1. In the
inherited notation, f(n)=4^-K(3/4)^(tau-1) on rooted n>1 and f(1)=1;
unrooted scores are zero. Put Z(n)=(3n+1)f(n)/4. For a strict ROOT word
a_0,...,a_(tau-1), let o be the number of odd letters. The terminal letter
is even. Q1 already proves

    Z(n)=2^-o product_(i=1)^(tau-1) 3n_i/(3n_i+1),
    0<=Z(n)<=1.

## C1. The reciprocal carry cost is controlled by odd charge alone

For every rooted positive odd n, including 1,

    (3/4)(3/8)^o < Z(n) <= 2^-o.                         (1)

This supplies the carry estimate that was separate from odd-charge control
in Q1. It is not a claim that o is bounded uniformly over all sources.

**Proof.** Every internal nonroot state n_i, i>=1, is coprime to 3 and
therefore at least 5. Partition these states into blocks whose consecutive
connecting letters are even. There are at most o+1 blocks. Write z_i=n_i-1
and H(z)=3z/(3z+4). An even letter has a>=2 and gives

    z_(i+1) <= (3/4)z_i,
    3n_i/(3n_i+1) = H(z_i)/H(4z_i/3).

Since H is increasing, the product in each block telescopes from below to
strictly more than H(z_last)>=H(4)=3/4. Strictness follows because the
first denominator H(4z_first/3) is less than 1. Multiplying the blocks proves
the lower bound. An empty internal list has product 1. The upper bound is
the inherited product with all carry factors at most 1. ∎

Consequently, independently establishing at most B odd letters for the
marked source supplies

    f(n) > f_min(n,B) := 3/(3n+1) * (3/8)^B,
    W(n) >= (9/16) f(n)^2 > (9/16) f_min(n,B)^2.          (2)

The second inequality uses the existing beta-mixture transfer, not a new
weight definition. For n=1 the displayed conservative bounds also hold.

## C2. A finite deadline from the same budget

For n>1 let t be the least nonnegative integer such that

    (n-1) 2^B 3^t < 2 * 4^t,
    T(n,B)=B+t.                                         (3)

An actual source reaches ROOT or spends its (B+1)st odd letter by time T.
Indeed an even letter contracts n-1 by at least 3/4. An odd letter expands
n-1 by at most 2; this follows already from the worst case a=1 and n>=3.
If all T endpoints were nonroot and at most B letters were odd, then

    2 <= n_T-1 <= (n-1)2^B(3/4)^t < 2,

a contradiction. This also proves that a positive infinite orbit must
spend unbounded odd charge; no assumption excluding positive cycles is used.
At n=1 take T=0.

This bound counts expensive events while allowing arbitrary terminal
exponents. It is deliberately conservative. It does not make the budget B
available for an unproved source.

## C3. Exact optimal shadow score at each binary precision

For this paragraph only, extend U by U(1)=1 with valuation 2. This is a
finite-bit observation convention, **not** permission to pad a ROOT proof.

Fix a positive odd n and d>=1. Read the longest initial word w whose total
valuation A is strictly less than d. Stop just before the next letter would
bring the total to at least d. Denote its endpoint by y, its length by t,
its odd charge by o_d(n), and its affine numerator by

    U_w(x)=(P x+B_w)/2^A,       P=3^t.

Every letter read is determined by n mod 2^d. The integer o_d(n) is
nondecreasing in d. Define a **supremum**, not an average or actual source
score, over rooted shadows of this same marked residue:

    E_d(n)=sup {Z(m): m>=n, m=n mod 2^d, m is rooted}.

Then exactly

    E_d(n)=2^-o_d(n).                                   (4)

**Upper bound.** Such a shadow shares w. For m>n, its intermediate states
are strictly larger than those of n, so it cannot terminate earlier along
w, even if the observation convention has reached 1 at n. It has at least
o_d odd letters, and Q1 gives the bound. The case m=n is immediate.

**Attainable supremum by changing the terminal controller.** The unfinished
letter says 3y+1=0 mod 2^(d-A). Since 4 generates the subgroup 1+3Z modulo
3^(t+1), there is a unique phase

    k=k_0 mod 3^t,       4^k=3y+1 mod 3^(t+1).           (5)

This phase follows also by elementary ternary digit lifting: at each level
exactly one of k, k+3^j, k+2*3^j lifts. Choose k in this phase arbitrarily
large, with 2k>=d-A, and set

    h_k=(4^k-1)/3,
    m_k=(2^A h_k-B_w)/3^t.                              (6)

The coprime binary and ternary congruences give m_k=n mod 2^d; for large k,
m_k>n. The strict ROOT word at m_k is w followed by 2k. In particular,
formal observations at 1 in w become genuine nonroot states at m_k; no
completed certificate is padded. All intermediate shadow states tend to
infinity, so their finite carry product tends to 1. Their odd charge is
exactly o_d(n), proving (4). The supremum need not be achieved.

This construction does not hold the *approximating integer* fixed; it holds
the **requested source n and its residue obligation** fixed. Equation (4)
quantifies exactly what these changing controllers buy and what they omit.

## C4. The new coverage obligation, and a finite receipt interface

For each fixed n the following are equivalent:

1. n reaches ROOT.
2. The integer staircase o_d(n) is bounded as d tends to infinity.
3. inf_d E_d(n)>0.

If n is rooted, the staircase stabilizes at its actual odd charge; all later
observations at 1 are even. If not, C2 makes its odd charge unbounded, and
each finite prefix is eventually visible. Thus the equivalence uses actual
integer descent, not an exchange of typical and pointwise quantifiers.

Together with (1), a bound o_d(n)<=B for **all d** pays the original source
weight, not merely its neighbours. But that all-d bound is still the open
input. Equations (4) and (6) alone do not prove it: every source has rooted
shadows, and their optimal scores may tend to zero.

There is also a finite, independently checkable sufficient receipt. Given
n>1 and B, set

    H=1+2^B(n-1),
    b=bit_length(3H+1),       D=1+b*T(n,B).               (7)

Supply **one authenticated strict ROOT shadow** m>=n with
m=n mod 2^D and odd charge at most B. Then n is rooted and has odd charge
at most B. Before ROOT or the (B+1)st odd event, the source is at most H;
every next valuation is at most b. C2 puts that event within T letters,
hence within fewer than D bits. The two sources share that prefix, and the
(B+1)st event would contradict the supplied shadow's charge budget. This
proves the claim without pretending that an arbitrary low-precision shadow
is enough. The shadow may use a changing head and terminal controller; its
identity, binary precision and total odd charge are all retained.

This is a certificate interface, not an algorithm known to produce the
receipt for every source. A source-dependent B established by a separate
argument would make it useful. A B measured retrospectively from a completed
orbit merely repackages that already known proof.

## C5. Retaining the resulting floor through measurement refinement

Retain the original n=2m+1 and the inherited native moduli

    q_d=((3n+1)/2)(4^d-1)/3,
    A_d=nu(m)+4^-d nu(m+q_d),
    nu(j)=(j+1)^-6-(j+2)^-6,
    E_d^alias=9/(448 q_d^8).

Once a budget receipt in C4 has independently grounded n, the lower readout
L_d=A_d f_min(n,B) is justified. Choose any d with

    E_d^alias <= nu(m) f_min(n,B)/2.

Such a d is found by rational integer arithmetic. For this d and every
larger d, the inherited marked-source readout inequality returns at least
f_min(n,B)/2. Choosing d=2^h gives nested moduli, since q_d divides q_(2d).
This transports a paid floor through refinement; it does not supply the
universal budget that would pay it initially.

## Validation and remaining helpers

The [checker](../../04-computation/experiments/changing_controllers_rhombus_20261011.py)
and [exact output](changing_controllers_rhombus_20261011.json) replay all odd
n<=10001, with cap expiry treated as failure to complete the test. They
separately verify the product identity, bound (1), deadline (3), and 224
changing-terminal shadows at odd n<=63 and d=1..7. The proof, not this finite
population, supplies (1)--(7).

The combined packet completes **194,628 always-active exact checks** with
byte-identical ordinary and optimized outputs:

```text
python3 -B 04-computation/experiments/changing_controllers_rhombus_20261011.py
python3 -B -O 04-computation/experiments/changing_controllers_rhombus_20261011.py
```

Source SHA256: `3092653c9a1c39ecaa01eaacb1799701926eb7cf4b109c7f7801af626da59851`.
Output SHA256: `506e4530bc68ac86bf843010f8224f1a83cc1e541cca4269ae5e979b3a210542`.
The tests include separate source-budget stopping events and retained
readout floors; the finite rooted population does not prove a universal
charge budget. The geometric controls include positive colour lifts, not
unverified solver declarations. The session uses the existing research
moves of retaining marked sources and complete incidence; no new universal
meta-pattern is promoted from this one packet.

The next useful helpers are: find a source-specific bound on odd charge
through a recursive controller grammar; produce high-precision low-charge
shadows without first running the entire requested orbit; and ask whether
the retained-credit mechanisms in [the zeta/credit bridge](collatz_credit_zeta_bridge_20261005.md)
can bound this integer staircase. The fixed-head obstruction rules out one
finite library, not these variable-head interfaces. Universal coverage
remains OPEN.
