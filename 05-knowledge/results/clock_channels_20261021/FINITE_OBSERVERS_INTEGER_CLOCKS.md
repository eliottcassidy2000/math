# Finite observers do not retain an unbounded integer clock

2026-10-21. **PROVED elementary, with exact literal controls.** This is a bounded information-flow companion, not a universal Collatz proof or an impossibility theorem for finite alphabets, finite rule schemes, or counter machines. The source family is an explicit specialization of the inherited braid map B(x)=4x+1; no historical novelty is claimed.

## 1. The exact factor criterion

For a deterministic map F:S->S and observer pi:S->A, a map f on pi(S) satisfying pi F=f pi exists iff

 pi(x)=pi(y) implies pi(F(x))=pi(F(y)).

Necessity follows by applying f. Sufficiency defines f(pi(x))=pi(F(x)), with the displayed implication making the choice independent of x. If the physical transition also has integer cost ell(x), retaining that cost requires the stronger condition that (pi(F(x)),ell(x)) be constant on every pi-fibre.

This is the deterministic form of a coarse-graining/lumpability test, not a new theorem of information theory. The exact sources63 and273 from the idoneal-information packet have the same compatible15/35 register, but primitive successors95 and51 and braid clocks1 and4. Thus their residue observer fails both coordinates. Idoneality alone also fails to factor the fully odd map:3 and7 share the flag but their successors5 and11 do not.

The next example isolates the clock obstruction even when the terminal is identical and actually grounded.

## 2. One actual grounded family and three different clocks

For j>=0 define

 s_j=(32*64^j-5)/9,   r_j=3j+1.

These are positive odd integers, increasing from s_0=3 and congruent to3 modulo8. Integrality follows from64=1 modulo9. Put

 x_j=(3s_j+1)/2=(4^(r_j+1)-1)/3=B^(r_j)(1).

Every B-iterate before1 is5 modulo8. Hence the **unstopped primitive** transition P takes its a=1 branch, strips exactly r_j B-ancestors, and gives

 P(s_j)=1,  primitive braid clock=1+2r_j=6j+3.          (1)

The actual fully odd Collatz word is instead

 s_j --valuation1--> x_j --valuation(6j+4)-->1,         (2)

because3x_j+1=4^(r_j+1). It reaches ROOT at shortcut clock6j+5 and ordinary unaccelerated clock6j+7. Both odd steps are literal and there is no substituted source. The primitive clock in(1) is not either physical clock in(2).

The first literal descent has shortcut clock4 for every j. Indeed the first four shortcut values are

 (3s+1)/2,  (9s+5)/4,  (9s+5)/8,  (9s+5)/16.

The first three exceed s, while the fourth is below s for s>=3. These are integers along(2), since its final valuation is at least4. Its first odd descent occurs only at ROOT, at clock6j+5. A bounded first literal exit therefore coexists with unbounded odd-exit and ROOT clocks in a genuinely grounded family.

### The actual first-arrival frontier is not forgotten

Let V be any fixed finite frontier containing1. For every sufficiently large j, s_j exceeds max(V), hence is not already in V. Its first P-arrival to V is then the single edge to1, with depth1 and clock6j+3. Small s_j that already lie in V instead have stopped first-arrival depth and clock zero. In particular(1) must not be called the stopped frontier clock at j=0 without checking membership. The eventual statement proves unbounded actual first-arrival clocks for the inherited finite627-point frontier without guessing which small sources it contains.

## 3. Exact same-register pairs with different clocks

For any positive modulus M, write M=2^a M_o with M_o odd. Let L be a positive multiple of the multiplicative order of64 modulo9M_o. Then

 s_(j+L)-s_j=32*64^j*(64^L-1)/9.

Its odd part is divisible by M_o, and its2-adic valuation is exactly6j+5. Consequently if6j+5>=a, the two sources are congruent modulo M. One explicit permissible choice is L=phi(9M_o), by Euler's theorem. Their primitive and both physical total clocks differ by6L, not zero.

For the combined modulus840=lcm(15,35,8), one can take L=6:64 has orders3,2,1 modulo27,5,7 respectively, so64^6=1 modulo945. The2-adic condition holds for every j>=0. Thus

 s_(j+6)=s_j modulo840,
 primitive clocks differ by36,
 and those clocks are equal modulo12.                 (3)

For example s_0=3 and s_6=244335917283 share the registers, with unstopped P clocks3 and39 and physical shortcut ROOT clocks5 and41. For sufficiently large j both sources in(3) are outside any chosen finite V, so the same phenomenon concerns actual first arrivals, not just unstopped transitions.

More generally a multiple of L can also make6L divisible by any prescribed finite clock modulus. Keeping finitely many residue clocks never identifies the full integer clock on this family.

## 4. A precise static-state obstruction, and its limitation

On the family {s_j}, the full integer clock6j+3 is injective in j. No finite-valued static observer can reconstruct it exactly. For any q-valued observer and a uniform choice of j among0,...,J-1,

 H(j | observer) >= log_2 J-log_2 q.

This follows from H(j)=log_2 J and H(observer)<=log_2 q. It is an information-loss bound for that explicitly stated finite uniform experiment, not a density claim about all integers.

A finite palette, the65 idoneal labels, or a compatible CRT register by itself therefore cannot serve as the entire source-and-clock record for this family. But a finite alphabet with an **unbounded word**, or a finite controller with an **unbounded accumulated integer counter**, can retain arbitrarily large clocks. Nothing here rules those out. It also does not show that a sufficient convergence proof must store every exact clock; it addresses claims that a proposed map preserves the actual clock.

The idoneal genus channel in the companion has a different consumer: it retains multiplication inside one fixed quadratic order when the square subgroup vanishes. Its losslessness for that operation is not losslessness for a timed integer dynamical system whose source changes the order. The requested creative analogy becomes an exact rule: a compressed state is adequate only for the operation and output it is proved to determine.

## 5. Reproduction and scope

Run python3 verify_clock_channels.py and python3 -O verify_clock_channels.py. The finite checker independently performs primitive stripping, fully odd iteration, shortcut and ordinary iteration, first-descent detection, source congruence and exact modulus-period controls. It tests j=0..30 and declared modulus pairs, including(15,35,8,840). The all-j statements and finite-observer obstruction are proved above; no infinite trajectory is inferred from a finite replay.

Dependencies: native_cancellation_105_20261018g.md for the primitive map/clock convention; guarded_children_105_20261018g.md for the distinction between shortcut and ordinary cost; the idoneal-information packet only for its separately proved63/273 and flag examples. Neither universal idoneal completeness nor a uniform annulus inequality is used.
