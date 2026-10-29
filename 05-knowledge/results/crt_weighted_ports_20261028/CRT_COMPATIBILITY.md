# Proper clocks, reflected probes, and compatible CRT quotients

2026-10-28. **PROVED scoped / FINITE-EXACT; same-agent audit.**
No general composite-Lucas classification or prime-value distribution claim.
The source is the full set of rooted legal words, including zero and all
nonunit strata. The operation menu is stated separately for every result.

## 1. Inherited reader and an all-size proper-clock obstruction

Put O=Z[phi], phi²=phi+1, R_n=O/(phi^n-1), for odd n. Let W_n be the
binary words indexed modulo n with no adjacent ones, and use

    q(w)=sum_i w_i phi^i,       R(w)_i=w_(-i).

The [Session74 all-period reader](../session74_eleven_mechanisms_20261027/session74-clocks/RESULT.md)
proves q:W_n→R_n is bijective. Our exponent convention reverses its written
word convention; the words and retained rotation action are the same.
Transport R through this bijection, without assuming it is additive.

**Theorem C1 (proper-clock obstruction).** Let A be a nonzero unital ring
quotient of R_n, and let a be the image of phi. If a has multiplicative
order m<n and word reflection descends to A, then A=F4. Conversely this
four-state quotient exists when 3 divides n, and reflection on it is Frobenius.

Proof. Since a is a unit and a^n=1, m divides n. The case m=1 would imply
1=0 by the golden quadratic, so m>=3. Since n is odd and m is proper,
n>=3m. The two legal words with occupied positions

    u={2},       v={0,m+1}

have the same image a²=1+a^(m+1)=1+a. Their reflected images differ by

    a^-2-(1+a^(-m-1)) = -2a^-1.                         (C1)

Thus descent forces 2=0. A is then a nonzero quotient of
O/(2)=F2[phi]/(phi²+phi+1)=F4, and hence is F4 itself. In F4, a^-1=a²;
reflection of any word evaluates to the square of its original value.
Its root has order three, so it is a quotient of R_n exactly when 3|n. QED.

This excludes *nonlinear* descended reflection maps on these proper clocks,
not just a ring-homomorphism shortcut. It leaves full-period factors open.
The weaker additive obstruction in incoming Session94 is consistent with
(C1); it is not being replaced by an assertion that all quotients fail.

## 2. Tripling at a Lucas prime costs more than two field probes

For odd m, write p=L_m and Q=1+phi^m+phi^(2m). The trace identity gives

    L_(3m)=p(p²+3).

When p is prime, m>=5 and p is odd. The ideals (phi^m-1) and (Q) are
comaximal: Q reduces to 3 in R_m=F_p. Consequently

    R_(3m) ≅ F_p × O/(Q),       |O/(Q)|=p²+3.            (C2)

**Corollary C2.** With the first F_p coordinate already known, two additional
p-valued probes cannot identify every state of R_(3m), even if the second
probe depends on the first answer. Each base fibre has p²+3 states, while
a depth-two p-ary decision tree has at most p² leaves. At least three
probes are necessary; this does not assert they suffice for a specified menu.

Examples: p=11,29,199,521 give hidden fibres124,844,39604,271444. Also
(p²+3)/4 is the Eisenstein-field norm of (p+sqrt(-3))/2. This is an exact
norm identity; it does not identify the quadratic fields or their operations.

The source-preserving repeated-word map is instead J(x)=Qx. Its projection
back to R_m is 3x, not x. Repeating the presentation triples its executed
integer clock and odd charge. The full fibre in (C2) must not be confused
with the p-element image of J.

## 3. Complete finite CRT closure beyond 76

The script `crt_probe.py` enumerates *every* legal word at the listed lengths.
It evaluates original and reflected words in the indicated CRT factors,
checks the bijection and involution, and tests every nonempty factor subset.
These rings have squarefree odd factors; no repeated-prime theorem is inferred.
The incoming [all-odd normal form](../session94_recursive_ports_crt_clocks_20261029/crt/COMPATIBLE_CRT_CARRIERS.md)
handles the general prime-power charts separately.

| n | full state count | CRT factors | proper nonzero reflection-closed factor projections |
|---:|---:|---|---|
|9|76|F4 × F19|F4|
|15|1364|F4 × F11 × F31|F4|
|21|24476|F4 × F29 × F211|F4|
|23|64079|F139 × F461|none|
|25|167761|F11 × F101 × F151|none|
|27|439204|F4 × F19 × F5779|F4|

At n=27 the F19 output takes all19 values over *every* fixed F19 input
fibre. Conversely the number of possible F5779 outputs over a fixed F5779
input ranges from1 to75. Thus the period-nine rule “all hidden tags give
distinct reflected outputs off zero” is not a general CRT law.

There is a precise behavioral consequence. Allow all additive translations,
multiplication by phi, and R, and observe any one odd field factor in this
table. Every compatible equivalence is a coset partition by an additive
subgroup I: translation invariance gives x~y iff x-y~0. Multiplication by
phi makes I an O-ideal. In these products of fields, ideals are exactly
factor-deletion ideals. The table therefore proves that the smallest exact
carrier for this menu is the **entire ring**, at each listed n.

This consequence is FINITE-EXACT. It does not assume every subgroup is an
ideal without the phi operation. Session94's directly computed {+1,R}
closures at9,15,21,25 are stronger for that smaller menu and are credited
there. Its sparse period25 and55 certificates also go beyond our census;
we do not recast those results as new here.
The final incoming [weighted CRT census](../weighted_crt_clocks_20261028/CRT.md)
also closes the smaller {+1,R} menu at23. Our larger-menu result at27 is
kept distinct; the extra generators are not silently removed.

## 4. Constructive reflected probes

Given the observed factor pi, retain pi(x) and query pi(R(x+t)). Mixed CRT
translations are allowed here. F4 is encoded by0,1,omega,1+omega as0,1,2,3.
The following complete decoders were selected from a finite seeded menu and
then checked on every source state. These are upper bounds, not minima.

| n; observed factor | translations t in the displayed CRT order | number of reflected probes |
|---|---|---:|
|15; F11|(0,0,0), (1,5,5), (2,5,20), (0,5,12), (0,2,2)|5|
|15; F31|(0,0,0), (3,2,23), (0,2,2)|3|
|21; F29|(0,0,0), (2,14,176), (2,2,157), (1,1,176), (3,2,68)|5|

For example, at15/F11 the successive numbers of distinguishable cells,
starting after the first reflected probe, are111,754,1270,1358,1364.
At21/F29 they are785,11827,23011,24393,24476. Each query acts on the same
retained source; the translations are not independently resampled states.

For rotation and reflection alone, the pair (pi(x),pi(Rx)) already gives
the minimal closed dihedral carrier, as proved in incoming Session94.
The extra probes repair a larger menu. A compressed carrier's validity is
therefore a statement about future operations, not just present cardinality.

## 5. Boundary and next decisive tests

Source: rooted legal words and their full arithmetic label. Target: selected
CRT coordinates or query signatures. Preserved: named observations under
the declared operations. Lost: proper clock phase, odd charge, and sometimes
the nilpotent layer. Restore the actual missing coordinate or a proved
operation-compatible query signature; never infer it from cardinality alone.

**OPEN:** classify the full-period odd-primary factors which nevertheless
fail reflection descent; find all-size sparse witnesses for the resulting
translation closure; determine sharp probe costs for the permitted translation
menu. C1 supplies universal witnesses for proper clocks, not the remaining
primitive-clock factors. Prime-valued Lucas infinitude is untouched. The
incoming fixed-bank sieve theorem already explains why a fixed finite
divisor menu eventually becomes silent on Lucas prime-index candidates.
