# The 191-term row collapses to one primitive-element gate

**Status: independent audit reflection for proved THM-3517.**  The theorem
and exact companions are the proof sources; this note records why the large
`z` resultant is structurally simpler than its coefficient count suggests.

## Verdict

The new `m=3` all-coordinate result is sound.  A disjoint python-flint
implementation reconstructs, without importing Sympy or the candidate:

- the `28/25/3`-term determinant-one source map;
- the inverse quintic `w^5-w^4+Pw-Q` and its reduced irreducible `D5`;
- the `17`, `29`, and `191`-term `x,y,z` resultants;
- the three exact discriminant/index-square decompositions;
- the `17`, `26`, and `268`-term branch-coprime index cores; and
- the pullback `D5(BC,AC^2)=C^4L5`.

Every eliminant, index-core, and discriminant hash agrees with the incoming
Sympy ledger.  This closes the only named implementation-audit gap in
THM-3517.

## A tiny fibre sees what 191 terms encode

The independent audit also searches the inverse quintic directly over
`F_31`.  The first simultaneous separator is

```text
(A,B,C)=(26,23,1),             (P,Q)=(23,26),
w=(8,9,12,16,18).
```

The reconstructed source rows `(x,y,z,w,gamma)` are

```text
(28,29,19, 8,10),
(10,12,15, 9,28),
( 7, 3,13,12, 9),
(11,30,11,16,17),
( 6,23, 0,18,26).
```

Their three coordinate Vandermonde determinants are `(28,14,1)`.  So each
named source coordinate separates all five sheets.  The giant `z` resultant
is not hiding a second mechanism: it is the coordinate polynomial of another
power basis in the same quintic trace algebra.

## What the trace class sees and misses

Once `x,y,z` are primitive, trace-form congruence forces their discriminants
to differ by basis-change squares.  The inverse quintic supplies `[D5]`, and
target pullback gives `[L5]` because `C^4` is a square.

This common sign class still misses the genuine component `V(C)`.  THM-3448
shows that three sheets escape there in one `C3` orbit.  A 3-cycle is even,
so its discriminant order four disappears modulo squares.  The primitive
coordinate theorem and the Jelonek effectivity theorem therefore coexist
without contradiction:

```text
all three coordinate views primitive
        does not imply
their common sign resolvent detects every nonproper component.
```

The sidecar lost by the square-class quotient is the actual local inertia
cycle, not another coordinate resultant.

## K4, tournament, and XOR interpretation

With three primitive coordinate views and the inverse invariant `w`, one may
place four power bases at the vertices of a labelled `K4`.  The six edge
labels are ratios of basis volumes and obey a multiplicative coboundary law.
This is a useful six-edge carrier, but not intrinsically a tournament:
orientation only chooses whether an edge stores `g_ij` or its inverse.

After squaring, every discriminant ratio is square-trivial.  Before squaring,
a chosen square-class character gives an XOR coboundary, while the bare
predicate “is nonsquare” need not.  Thus the size-four tournament intuition
is accurate only after the orientation observable and character are supplied.

## The all-cyclic-family question closed during integration

The audit's proposed coefficient-level route succeeded for every cyclic
THM-3448 seed.  If `n=ell+2` and `H` is the numerator in

```text
z=H/(bC^2)=gamma[gamma(gamma-1+a)-aw]/(bC^2),
```

then the coefficient of `w^(n-1)` in `H mod T_n` is

```text
((ell+1)^3+1)P^2+(4ell^2+ell-2)P/2+4(ell+1)Q.        (1)
```

It is not the zero polynomial.  Therefore `H mod T_n` has degree `n-1`, so
`z` cannot lie in the base field; otherwise it would give an annihilator of
`w` below the irreducible degree `n`.  The `S_n` point-stabilizer maximality
identified by this audit then proves `z` primitive.  Thus THM-3517 now gives
all three primitive coordinates in every cyclic grade.  Direct exact
remainders through `ell=30`, including the hostile `P=Q=0` where the whole
remainder specializes to zero, accompany the symbolic proof.

This closes the question for the cyclic subfamily, not for every lawful
weighted seed.  A later clean-room audit supplied a second all-degree proof:
at `P=2^(1-n),Q=0`, the roots `0,1/2` have unequal `H` values by a uniform
rational inequality.  That route imports neither Sympy nor the recurrence
companion, and replays two finite same-target/distinct-`z` points through
`n=256` (with the sharp `gamma=0` cubic boundary retained honestly).  The
finite `F_31` separator above remains the independently coded degree-five
check of the full resultant atlas.

## Exact artifacts

```text
04-computation/jc_weighted_odd_family_m3_three_coordinate_independent_audit_20260816.py
05-knowledge/results/jc_weighted_odd_family_m3_three_coordinate_independent_audit_20260816.out
```

Their LF-normalized hashes are

```text
e0eaed0b1f034c799cd95a00ede7fdd60ff4c9258330b39845db45ea159e97e2
b458c072d8220403162bd972ab32dc3c4aaacd217cfb2ebf7a0525fc0e14f96d
```

and the independent semantic digest is

```text
bd208dad9732439dfa14a794ef54dbfe57d66360f5f5a39b26353ffb82b6bba3.
```

The FLINT artifact's scope remains the explicit cyclic weighted `m=3` member;
the separate symbolic companion proves the cyclic all-grade statement.
Third-coordinate primitivity for arbitrary weighted seeds, identification
with the unstored historical family, arbitrary Keller classification,
`JC(2)`, and LRC remain open.
