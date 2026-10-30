# Exact gcd-fibre compression and parity-dual attachment

**PROVED:** an all-core, all-dilation, all-attached-speed transfer formula,
with closed endpoints, singleton components, and an explicit native lift.
**FINITE-EXACT:** independent literal-comb verification in the accompanying checker.
**OPEN:** arbitrary-core LRC(14). Dilation closure adds no new primitive cores.

## Inheritance and precise improvement

Pinned source: `16a4cdb898c71d54dcb2c6cf58784c8e13692d22` (Session119).
The packet `05-knowledge/results/terminal_fees_phase_covers_20261030/LRC.md`
proves exact arbitrary-speed attachment from the complete safe interval set,
and its `lrc14_banks.json` certifies every positive thirteenth speed for 28
specified twelve-speed cores at radius 1/14. Its next compression question
asks how to retain all future attachment verdicts without expanding the full
cover. The present result gives exact compression for the entire dilation
family: use the original core's intervals, never the d-fold lifted cover.

`01-canon/theorems/THM-4446-lrc14-primitive-ten-pack-descent-and-dilation-rays.md`
already uses the open-arc grid cap and a three-tail fibre union bound. We do
not claim that fibre lifting or dilation closure is new. The contribution
here is the exact single-tail response, its odd/even phase duality, the
all-radius interval decision rule, and its integration with the full phase
cover compiler. THM-4025 concerns owner-residue odd dilations, not this exact
safe-set quotient. No external LRC-up-to-13 theorem is needed here.

## 1. Exact fibre response

Let C be a finite nonempty set of positive integers; let 0<r<=1/2, d,v>=1
be integers except for the rational radius. Set

    G_C(r)={s in [0,1]: ||cs||>=r for every c in C},
    g=gcd(d,v), h=d/g, u=v/g.

For each s, all d points t_j=(s+j)/d (0<=j<d) preserve the core:
`||(dc)t_j||=||cs||`. The attached phases form, with multiplicity g,
the translated order-h grid `(u*s+m)/h`, m modulo h. Consequently

    max_(0<=j<d) ||v*t_j||
      = (h-1)/(2h) + ||u*s-(h-1)/2||/h.                 (F1)

**Proof.** The closest grid point to 1/2 has distance
`dist(h/2-u*s,Z)/h` from 1/2. Therefore the largest norm is
`1/2-||h/2-u*s||/h`. Use `||x-1/2||=1/2-||x||` to obtain (F1).
The grid permutation is justified by gcd(u,h)=1. This also covers h=1.

A maximizing native point needs no enumeration. Choose an integer m nearest
to `h/2-u*s`, solve `u*j=m (mod h)` with 0<=j<h, and take t=(s+j)/d.
For h=1 use j=0. A modular inverse and exact rational rounding suffice.

## 2. Exact all-radius quotient, including the parity switch

Write a=h*r-(h-1)/2. Then

    G_(dC union {v})(r) is nonempty
    iff some s in G_C(r) satisfies
    ||u*s-(h-1)/2|| >= a.                              (F2)

This is an equivalence: every original time is a lift of a base time modulo
one, and (F1) optimizes over every lift, not just a chosen owner window.
If a<=0, every base-safe time has a successful lift. If a>0:

* odd h: the test is `||u*s||>=a`;
* even h: the test is `||u*s||<=b`, where `b=h*(1/2-r)`.

The even case is a **near-integer** condition, not an ordinary lower-clearance
attachment. Erasing this parity can reverse the verdict. For example h=2,
u=1,s=1/2 gives true fibre maximum 1/4, not 1/2.

Suppose the original complete safe set has k closed components [L_i,R_i],
including possible points. In the odd case set (lo,hi)=(a,1-a); in the even
case set (lo,hi)=(-b,b). For a>0 the exact criterion is

    exists i: ceil(u*L_i-hi) <= floor(u*R_i-lo).        (F3)

Indeed this says the interval [uL_i,uR_i] meets some closed translate
[n+lo,n+hi]. A successful integer n yields
`s=max(L_i,(n+lo)/u)`, then the lift above yields a native witness.
Thus the decision uses O(k) rational arithmetic operations **independent of
d's numerical size**, plus integer gcd; witness extraction additionally uses
one modular inverse. Bit complexity still depends on input bit lengths.
No residue table, d-fold interval expansion, or loss of singleton/endpoints
is permitted or required. Rational r is needed only for exact implementation.

## 3. Sharp fibre threshold and complete rejection-set transfer

For fixed h, all base-safe times lift automatically if

    r <= (h-1)/(2h).                                   (F4)

This threshold is sharp **per fibre**: choose u*s congruent to (h-1)/2
modulo one and (F1) is exactly its right side. These are realizable safe
fibres, not merely formal phases: take C={1}, s=1/2, d=h, and v=1
for even h or v=2 for odd h. The core clearance is 1/2 and gcd(v,h)=1.
This is not a claim that the
rejection-set law below has an optimal threshold globally over all cores.

For fixed d>=2, let p be its smallest prime divisor. Every proper orbit
h=d/gcd(d,v)>1 is a divisor of d and hence h>=p. Thus every nonmultiple v
of d is automatically accepted when `G_C(r)` is nonempty and

    r <= (p-1)/(2p).                                   (F5)

If d divides v, h=1 and acceptance is exactly that of v/d attached to C.
Define R(C,r) as the set of positive integer rejected attachment speeds.
Under the nonemptiness hypothesis and (F5),

    R(dC,r) = {d*w: w in R(C,r)}.                      (F6)

In particular (F6) holds for every d>=1 whenever r<=1/4; for every odd d
it holds whenever r<=1/3. The d=1 case is identity. Empty base-safe sets are
excluded: their rejection set is all positive integers, so (F6) would fail
for d>1. Beyond the guaranteed range, even a nonempty core need not suffice:
C={1}, d=2, v=1, r=2/5 has G_C nonempty but G_{1,2} empty (maximum 1/3).
This disproves an unrestricted all-radius dilation law, not optimality of
(F5) as a global threshold.

### Scoped LRC(14) corollary

For each of the 28 explicitly recorded cores C in the inherited bank,
every d>=1 and every v>=1 has

    G_(dC union {v})(1/14) nonempty.

This extends the stored finite list to its complete dilation closure without
any new finite-height search. It is a corollary of the inherited base-bank
certificate plus (F6), not an independent certification of arbitrary C.
After normalizing the full speed row, a nonmultiple attachment may remain
primitive, but the **core shapes** add no new gcd-normalized families.
No arbitrary twelve-speed core, no universal owner-window bridge, and no
unrestricted LRC(14) theorem follows.

## 4. Verification and remaining frontier

`fibre_transfer.py` uses only the Python standard library. Every check raises
explicitly on failure and remains active under `python3 -O`. It compares
(F1) and the constructed native witness with every literal lift on a declared
finite universe; compares (F2)-(F3) to an independently organized complete
closed-comb intersection; checks all low-radius rejection laws in that
universe; retains a high-radius counterexample, parity hostile, sharp fibre
points, and a dilation beyond 10^100 / attached speed beyond 10^150 witness.
Normal and optimized output must match byte-for-byte.

`recheck_bank.py` additionally re-derives every stored maximal-width tail
cutoff and repeats all 1,458 finite queries for the 28 inherited cores,
using both the quotient test and complete literal combs. The checked bank
is retained as `source_bank.json`; its scope is unchanged.

The actual unresolved bridge is now isolated: for **primitive core shapes**
not covered by the bank, prove that the base safe cover has no rejecting
attachment. Gcd fibre compression removes redundant dilations; it does not
supply the missing nonemptiness/acceptance assertion. For multiple new speeds,
one cannot maximize each on a different lift: the same label j must work
for all tails. THM-4446's joint-grid reasoning, not composition of (F1), is
needed there. The existing separate-versus-joint-arrivals hostile remains
binding.

Reproduce from this directory:

```sh
python3 fibre_transfer.py > check.normal.json
python3 -O fibre_transfer.py > check.optimized.json
cmp check.normal.json check.optimized.json
python3 recheck_bank.py > bank.normal.json
python3 -O recheck_bank.py > bank.optimized.json
cmp bank.normal.json bank.optimized.json
```
