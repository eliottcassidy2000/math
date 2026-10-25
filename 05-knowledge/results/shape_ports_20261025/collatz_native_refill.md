# Native re-entry immediately after the 27-word fuel is exhausted

**2026-10-25. PROVED explicit specialization / FINITE-EXACT controls.**
For every integer `h>=1`, infinitely many positive integer sources have
exact 27-word fuel `(1,0,h)` at shortcut times `(0,59,60)`. Thus exhaustion
of this particular word can be followed by arbitrarily deep re-entry in
one actual even step. This is a precise state-retention result, not a new
coverage theorem or a globally decreasing Collatz rank.

## Inheritance and the native boundary

Use the shortcut map `T(n)=n/2` for even `n`, `(3n+1)/2` for odd `n`.
The closest proved mechanism is Session53's
[repeated 27 block](../session53_g27_collatz_cooper_20261024/collatz27-repeated-block/REPORT.md),
independently checked in its
[audit](../session53_g27_collatz_cooper_20261024/audit-repeated-block/REPORT.md).
It retains the first 59 parity bits of 27, 37 odd operations and 96
ordinary operations, and gives the exact full-copy fuel. The general
fuel and refill mechanisms are inherited from
[guarded pumping, section 5](../collatz_guarded_pumping_memory_20261004.md)
and [paid controllers, section 7](../collatz_paid_portrait_controllers_20261004.md).
No priority claim is made for those mechanisms or parity cylinders.

The canonical hostile is the difference between a residue involution and
actual integer return; the corrected near miss is interpreting finite-word
fuel as a global rank. The least-used sidecar is the exact affine cofactor
at the first failed repeat. The concept board is **fixed word / integer
guard / immutable source / departure parity / refill cofactor / odd
residue quotient**. This note addresses the native transition immediately
after exhaustion; the companion bit-phase note identifies its exact time.

Write

    q=2^59, a=3^37, b=23q-27a, d=q-a,
    F(n)=(an+b)/q, ell(n)=v2(dn-b), K(n)=floor(ell(n)/59).

Here `b,d` are odd and `b/d` is strictly between 8 and 9, so the valuation
is finite for every positive integer `n`. On the exact source cylinder
`n=27+qt`, `t>=0`, the first 59 shortcut operations are native, their
endpoint is `23+at`, and they first descend below the original source at
shortcut time 59 / ordinary time 96. These are inherited statements.

For `t` odd, that endpoint is even. Therefore one more actual step gives

    n=27+qt  --59 shortcut steps-->  m=23+at  --one even step-->  z=(23+at)/2.

The exact deviations are

    dn-b = q(4+dt),
    dm-b is odd,
    dz-b = (c+ad t)/2,     c=23d-2b=31a-23q.                 (1)

Consequently `ell(n)=59`, `ell(m)=0`, and the next exact fuel is

    K(z)=floor((v2(c+ad t)-1)/59).                            (2)

Equation (2) is the complete refill formula on this specified odd-`t`
boundary. Its cofactor is not supplied by the preceding value `K(m)=0`.
All three displayed numbers are positive integers, with `n>m>z`:
`n-m=4+dt>0`. The original-source comparison has not been reset.

## An arbitrary prescribed exact refill

Fix `h>=1`. Set

    M_h=2^(59h+2),
    t_h=(2^(59h+1)-c)*(ad)^(-1) mod M_h,  0<=t_h<M_h.          (3)

The inverse exists because `ad` is odd. The residue `t_h` is positive and
odd: its defining numerator is odd, so it cannot be zero modulo `M_h`.
Every `t=t_h+M_h u`, `u>=0`, is therefore an admissible positive lift, and

    c+ad t = 2^(59h+1) mod 2^(59h+2).

Thus its valuation is exactly `59h+1`. Equations (1)--(2) give

    (ell(n),ell(T^59(n)),ell(T^60(n)))=(59,0,59h),
    (K(n),K(T^59(n)),K(T^60(n)))=(1,0,h).                    (4)

These events occur at ordinary times `0,96,97`. The new integer `z`
supports exactly `h` complete copies of the 27 word. They consume another
`59h` shortcut / `96h` ordinary steps and end with zero full-copy fuel.
This is a statement about that finite actual prefix; no outcome after it
is inferred. Each fixed `h` has an explicit infinite arithmetic
progression of positive sources. No single positive integer is asserted
to realize every depth or infinitely many prescribed refills.

The mechanism is elementary but operational: a failed repeat reaches an
even integer, division by two changes the relevant anchor comparison,
and the retained odd-coefficient cofactor controls the newly available
dyadic digits. The universal quantifier follows from (3), not from a
finite scan.

## The mod-seven sidecar still cannot predict the refill

Impose the extra condition `t=0 mod 7`. Since `M_h` is a power of two,
there is one suitable lift parameter `u mod 7`, and infinitely many
positive lifts remain for every `h`. Then all the examples have the same
three residue states:

    n=6 mod7,  m=2 mod7,  z=1 mod7.

They also share the exact initial bit depths `ell(n)=59`, `ell(m)=0`.
Nevertheless the next depth is `ell(z)=59h`, with `h` arbitrary. Hence
even `(current bit fuel, current residue mod7)` does not determine this
next native refill. This is a positive-integer hostile to that specified
state compression, not a claim that every finite or arithmetic state
description is impossible.

Connection contract: source is the actual guarded `59+1` integer path;
target is the fuel/residue record. Projection preserves the recorded
clocks and endpoint residues but loses the cofactor in (1). Retaining
`c+ad t`, or enough of its exact binary valuation for the intended
consumer, restores the next-fuel calculation. Two examples with different
`h` and the same residue record provide the decisive test.

## Independent checks and the remaining target

The standard-library checker imports no Session53 producer. It reconstructs
the native word, checks its rational periodic orbit with a fixed odd
denominator, and tests all `h=1,...,64`, each with lift parameters
`u=0,1,2,2^80+1`. All 256 sources are replayed through the first word, the
actual even step, and the exact `h` subsequent words. Original-source
descent, parity guards, clocks and exact exhaustion are checked. Separate
ordinary-map runs verify the first four baseline fixtures at times 96/97.
Sixty-four further fixtures retain the identical mod-seven record.

The checker also independently verifies the companion full bit-phase
claim on all integers `1..512` and on exact valuation boundaries
`ell=0..128` (641 tests): the native orbit agrees with the rational
periodic word for exactly `ell` bits, then disagrees. Its controls give
`n=27`: departure at shortcut time 61, integer 53, ordinary time 100;
and `n=27+q`: departure at time 59. This is a specialized replay of the
inherited odd-slope dyadic isometry, not a new general isometry theorem.

Reproduction from the repository root:

    python3 -B 04-computation/shape_anchors_20261025/collatz27_native_refill.py
    python3 -B -O 04-computation/shape_anchors_20261025/collatz27_native_refill.py

[The receipt](collatz_native_refill.json) is identical in normal and
optimized modes, SHA256
`4a5192a3fbc029698dc22493213c369f74d9510748f63bf78a8ba9eeb5c0b1ef`.
All checks raise explicit exceptions; optimization cannot disable them.

This family was already known to descend at its first 59-step block, so
the refill theorem adds no new global descent coverage. It identifies a
precise native transition and the information its quotient forgets.
The remaining target is a lawful rule for the *other* exhaustion exits,
retaining their exact cofactor and the unchanged original source, together
with a coverage argument. No unrestricted convergence or total-stopping
bound follows from the present family.
