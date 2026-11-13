# Which arithmetic geometry is actually visible to the grammar?

**PROVED elementary deductions; CITED odd-cycle polynomial; FINITE-EXACT
controls. No Collatz convergence or Chan–Pak coverage is asserted.**

## 1. The sequence in the question is a residue coordinate

For n≥1, a_n=floor(3n/2)+1 enumerates

    2,4,5,7,8,10,11,13,14,16,17,19,... .

These are exactly the integers at least2 outside the zero residue modulo3.
They include composite numbers. Including1 gives the positive units modulo3.
This is the precise useful predicate: eligibility for an odd Collatz inverse,
not primality and not irrationality.

For ε∈{+1,-1}, consider the accelerated positive odd map

    T_ε(x)=(3x+ε)/2^ν2(3x+ε).

At a positive odd target y, its inverse candidates are

    x=(2^k y-ε)/3,        k≥1.                      (1)

They are integral iff 2^k y≡ε mod3, so y must be a unit modulo3.
The exponent k is even if y≡ε, odd otherwise. All integral candidates
are positive and odd, and have the stated exact halving depth, because y
is odd. Advancing two exponents on the same sign sheet gives

    x(k+2,ε)=4x(k,ε)+ε.                             (2)

Advancing three exponents and switching sign gives instead

    x(k+3,-ε)=8x(k,ε)+3ε.                           (3)

Thus 4±1 is a literal predecessor recurrence; 8±1 occurs in the inverse
numerator, whereas the cross-sheet predecessor recurrence has offset±3.
Keeping the division by3 is essential. These are identities, not a proof
that every positive source eventually reaches1.

The mod8 split determines the first halving depth:

| x mod8 | ν2(3x+1) | ν2(3x-1) |
|---|---:|---:|
|1|2|1|
|3|1|at least3|
|5|at least3|1|
|7|1|2|

This sharpens the earlier [signed modular-table work](../odd_tables_signed_carries_20261111/SYNTHESIS.md)
and [radix transport](../collatz_radix_families_20261018p.md): residue, sign,
exponent and integer carry are separate coordinates.

The concurrent [clone-clock dictionary](../clone_clocks_negative_collatz_20261113/COLLATZ.md)
adds the exact normalization r=(H+1)/2: a permitted tournament move
H→3H+2σ becomes r→3r+(σ-1), for σ=0,1,2. The nonzero offset is σ-1,
not σ. On odd H, with T(x)=(3x+1)/2 at odd x and T(x)=x/2 at even x,
one has3H+2=T(2H+1)=2T(H)+1. This is an arithmetic identity; its tournament
realization still requires the appropriate ready source.

In particular the formal nodes r7=4=3·1+1 and r21=11=3·4-1 locate the
exceptions in that arithmetic tree, but do not prove their exclusion or
realize the other branches. The [rooted atomic grammar](../ternary_core_grammar_20261113/GRAMMAR.md)
has the actual surplus obstruction. A scalar recursion loses precisely
the source and cross-matrix information the selector must recover.

For the one-division signed maps, x≡∓1 mod2^s characterizes s consecutive
odd-branch applications. On positive sources these are s strict ascents,
apart from the minus-map fixed point1. Thus4±1 and8±1 also encode two
and three such applications. This branch clock is distinct from the
single accelerated halving-depth table and the inverse recurrences above.
The incoming wording's off-by-one and fixed-point boundary were corrected
and checked in its verifier's `--dictionary-only` mode.

## 2. A genuine mod3 law for tournament path covers

Let P_T(z) count ordered z-bin spanning path covers, permitting empty bins.
The principal specialization of the tournament odd-cycle expansion gives

    P_T(z)=Σ_C 2^|C| z^(n-Σ_{cycle in C}(|cycle|-1)), (4)

where C ranges over vertex-disjoint directed odd cycles of length≥3,
including the empty packing. This is the specialization of Irving–Omar,
Corollary19 (crediting Grinberg–Stanley Theorem1.39), with reversal of
within-bin paths to reconcile the descent convention.
[Primary source](https://arxiv.org/html/2412.10572v3#S3.SS4).

All exponents in(4) have parity n. Since P_T(1)=H(T), it follows that

    P_T(2)≡(-1)^n H(T) (mod3).                      (5)

If H(T) is a unit modulo3, its two-bin response divided by H(T) records
the vertex parity modulo3. This is an exact map from a graph invariant to
the two unit classes. It forgets the graph, integer size, endpoint counts
and higher cycle information. It does not identify tournament growth with
Collatz iteration.

Adding an ordinal singleton preserves H and doubles P(2), flipping this
mod3 phase. In [the compensated23-core repair](REPAIR23.md), the same doubling
is balanced by a core coefficient60→30, preserving the exact next target.
At modulus7 a different clock, ord_7(2)=3, produces the obstructed a≡0
mod3 phase. The shared mechanism is finite residue phase plus unbounded
integer amplitude; the underlying graph operations remain different.

## 3. Why 17 and224 form an unbounded exact family

For R_n, the transitive tournament with only its extreme arc reversed,
put m=n-2≥1. Every cycle contains that arc. A length k+2 cycle chooses
k middle vertices; for odd cycles k is odd. Two cycles cannot be disjoint.
Equation(4) gives

    P_R(z)=z^(m+2)+2Σ_(k odd) binom(m,k)z^(m+1-k),
    H(R_n)=2^m+1,
    P_R(2)=4·2^m+2·3^m-2.                          (6)

At n6 the pair is(17,224). Its ordinal powers multiply these two
coordinates, giving the family23·17^b+60·224^b2^a. Here powers of2,
powers of3 and the odd core count arise from counting subsets and cycles,
with an explicit reconstruction. No numerological identification is needed.

## 4. Rational, irrational and transcendental observers

The [prior geometry audit](../changing_core_returns_20261113/geometry/HISTORY.md)
explains the finite strict-sign observer boundary: a nonempty open chamber
cannot distinguish dense arithmetic types merely by its finitely many strict
polynomial signs. Equality constraints and arithmetic sidecars can.

For the present problem, actual rational ratios of integer responses are
useful: [neutral context slopes](TRANSPORT.md) classify which repair vectors
survive a specified gluing operation. A real or formal point satisfying the
same linear equations still need not be an actual tournament. There is no
observer-independent optimal arrangement: the23-core preference reverses
between immediate entry and the prefixed later chain operation.

Reproduce the independent permutation/cycle checks and signed inverse
identities with `python3 -B 05-knowledge/results/surplus_selector_20261113/arithmetic_probe.py`.
