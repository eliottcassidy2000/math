# The incoming paired carry grammar and the fixed +4 ray

**PROVED scoped boundary / FINITE-EXACT intersections.** The incoming
[paired-source construction](../../h13_paired_descent_20261111/construction/PROOF.md)
and [repeated-digit actuator](../../h13_paired_descent_20261111/square_graceful/ADDITIVE_CARRY.md)
are used with their actual sources and endpoint sidecars. They are not
reinterpreted as arbitrary-baseline additions.

## 1. Four exact address inverses

The uniform reset gives the target102*2^k+17. A pure-H pair of width2
whose lower value is this target would provide a literal witness for the
remaining102*2^k+19. A width4 pair can instead start at102*2^k+15,
which is three times the old source count34*2^k+5.

The incoming H13 and H25 formulas yield exactly the following equations
for their upper member to be the desired +4 target:

| Carrier | Left address | Equation in B | Gap below upper value |
|---|---:|---|---|
|H13|0|13g(B)+32B+6(B mod2)=102*2^k-1|2|
|H13|1|13g(B)+136B+6(B mod2)=102*2^k-59|4|
|H25|0|25g(B)+122B+16(B mod2)=102*2^k-31, B odd|2|
|H25|1|25g(B)+426B+16(B mod2)=102*2^k-151, B>0|2 if B odd;4 otherwise|

These are instances of the proved complete `linear_inverse`, so each
fixed k is decided over ALL addresses, without a word-length cutoff.
For k1,...,30 the first, second and fourth rows have no solution. The
third row has exactly one: k7,B97. Its actual paired states(H,E,U,d) are

    (13073,25,17141,25),     (13075,25,17141,25).

Both14-vertex graph rows are saved and checked by independent first-vertex
and reverse Hamilton-path tables. Thus this is a real paired construction
at one desired baseline, not merely a numerical difference of two counts.
It remains necessary to connect the uniform-reset source at H13073 to the
lower H25 witness by an allowed source operation. This finite meeting does
not supply that operation, and the complete address inverses do not give
an all-k absence statement.

The incoming C3 unit pair has lower value3g(B)+36B+30, divisible by3.
It therefore cannot be the lower member at102*2^k+17, which is2 modulo3.
This exclusion is universal. The other orientation realizes a lower target
two below the reset value and does not supply the desired +4 branch.

## 2. Drift cancellation between sources is not cancellation at the target

The incoming periodic criterion writes every fixed-prefix, fixed-word
linear-g response as

    C*2^(ell*m)+d*D*m+E.

For a mixed repeated word, D>0. Two different carriers can have
d1*D1=d2*D2, as the incoming13/25 period14 example proves. Their positive
linear terms then agree, which is useful for a repeated source equality.
The common value still retains that nonzero linear term.

Consequently such a fixed-seed mixed-period family cannot equal
102*2^(ell*m+j)+19 for every m, for any fixed j. Subtraction gives
A*2^(ell*m)+d*D*m+C0. If A=0 this nonconstant linear function has at most
one zero. Otherwise its real second derivative has a fixed strict sign,
so it has at most two real zeros. Thus there are at most two meetings
with that matching-rate target subsequence. This is a direct consequence
of the inherited exact drift identity, not a refutation of paired returns.
Changing the seed or the offset with m is outside this statement.

The repeated hexadecimal9 actuator instead changes a defect by4 at each
copy and solves a genuine simultaneous source equality after preparing an
arbitrary prefix. It also changes the absolute target. Applying it to an
address chosen for the fixed +4 ray therefore introduces a separate target
matching equation. That equation is not implied by annihilating the rebase
defect. The incoming future-response comparison adds a further constraint:
the H13 source is cheaper under the common continuations than its H3 rebase.

## 3. Connection contract and next useful freedom

The source is a marked pair with a fixed pure-H gap; the target is the
desired ternary branch. The exact map is the complete address inverse in
section1. It preserves the common(E,U,d) sidecar within each pair, but
forgets how the old reset source reaches the chosen pair. The cheapest
positive test is the explicit k7 lower graph; the missing certificate is
a compatible source change at H13073.

Together with [BIFLANK](BIFLANK.md), the result separates three obligations:
an exact target equation, a repeated-source equation, and an actual source
operation. A changing carry register might solve the first two jointly;
the present work does not claim that a finite family label does so.

Run `python -B carry_ray_probe.py` for the120 complete inverse decisions
and literal paired control in `carry_ray_probe.json`.
