# Independent audit of the all-parity physical and selected-network ceiling

**Status: ANALYTIC AUDIT PASS after the empty-defect correction;
FINITE-EXACT independent physical head and coefficient-box PASS.** This
audits the corrected [candidate](overnight4_20260906_lrc_parityfree_probe.md):
for every primitive sorted distinct positive ternary-unit triple,
`mu(F_w)<=min_i E_i(w)<=6/55`, with either equality only at `(1,10,11)`.
There is no oddness hypothesis. The selected-network finite columns are
outside this verifier's independent numerical path; the separate root
native interval replay has also passed, with 181,385 gates. Together these
audits accept both the physical and selected-network theorems relative to
the incoming proved additive result and their complete exact finite bases.

## 1. Inheritance and the repaired corner

The relevant inherited result is
[THM-4434 / universal scale-three network projection bound](../../01-canon/theorems/THM-4434-lrc14-universal-scale-three-network-projection-bound.md).
Its address identity, zonotope sections, quadrature, interval counts and
Minkowski argument are parity free. Its even-norm coefficient box and
even-norm cutoff are not. The incoming
[parity refutation](lrc14_parity_empty_core_sep06.md) correctly preserves
`(2,5,7)` and `(2,11,20)` as distinct omitted obstructions. The incoming
[additive theorem](lrc14_additive_parity_empty_core_sep06.md) is a proved,
independently audited dependency only for the selected additive projection.
THM-4437 / all-parity-network-reduction-to-three-low-circuits was inspected
as a **RESERVED / UNPROVED EMPTY STUB** and is not used in this proof graph.

Dropping parity initially lost one premise of THM-4434's strict count:
the allowed defect list need not be nonempty. Pattern `(0,1,2)` has no
admissible integer defect, so `N=F=B=0` and `N<F+B` is false. The primitive
unit-speed example `(1,2,4)` has relation `(0,-2,1)` and physical mass zero.
The corrected candidate now discharges the empty-list case first and uses
the strict sum only for a nonempty list. The independent complete box
confirms that `(0,1,2)` is its unique empty-list pattern. This is a repaired
proof boundary, not a counterexample to the proposed ceiling.

The live comparison board is: literal physical sheets; saturated integer
addresses; owner residues; section width versus roof value; coefficient
convexity; and shortest-relation scale. Removing parity changes the allowed
small circuits and can remove every admissible affine slice. Both changes
must be retained before applying the old count consumer.

## 2. Full-lattice and analytic reduction audit

For primitive `w`, the map `n mod Zw -> w cross n` is onto the full integral
orthogonal lattice: choose integral `u` with `u.w=1` and lift a carrier by
`n=-u cross C`. Its kernel is exactly `Zw`. Thus no unsaturated sublattice
or omitted affine class enters the proof. From the literal sheet predicate,
`y=3x` gives nearest-integer errors of radius `3/14` and owner
`-w_i^(-1)n_i mod3`. Translation `y -> y+1` shifts all owners equally;
the three inverse branches cancel the Jacobian `1/3` in physical measure.

The verifier enumerates nearest-integer owner words, rather than solving
the producer's carrier equations. All **192** affine residue fibers agree:
all-unit relations retain two longitudinal classes on defects divisible by
three; a relation with one zero residue retains one class on each of the
other defects. Primitivity and ternary-unit speeds exclude two zero residues.

The zonotope area `9(a+b+c)/49`, central width bound `6c/(7M)`, and monotone
rectangle discrepancy imply `F/c<=6/49+4/(7M)`, independently of parity.
Resetting the width to zero at strict support endpoints is necessary when
an actual zero coefficient gives a vertical endpoint edge. Open interval
counts give the stated `4/3` or `1` intercept per nonempty admissible slice;
the explicit empty-list split above completes their scope.

The `M<=18` universe contains every admissible magnitude pattern after
permuting coordinates and globally reversing signs. A positive speed
solution has one isolated sign, including when a coefficient is actually
zero. The verifier enumerates sorted triples directly, with no norm-parity
filter. It finds **747** patterns: **698** full support and **49** actual-zero.

For each defect, the producer optimizes by fractional knapsack. This audit
instead intersects the error cube with its affine plane by fixing two cube
facets and solving the third coordinate. Linear objectives attain their
extrema among those vertices, including actual-zero sections. All normalized
speed-polytope vertices are included. Width is a maximum of linear functions
minus a minimum, hence convex in normalized speeds; its vertex maximum
therefore bounds the full polygon. Every resulting exact slope is at most
`15/98`, with the sole boundary pattern `(1,7,8)`. The complete semantic
digest matches the producer, not merely its maximum.

The projected relation lattice has index `c` since `gcd(a,b,c)=1`. The
projected l1-ball area coefficient is strictly greater than `3/4` for
distinct positive normalized speeds; the displayed positive polynomial
identity proves this without parity. Planar Minkowski therefore supplies
`S<4sqrt(c/3)`. The count gap is `547/5390`, giving
`A=1540/547`, `B=21560/1641`. For `S<=18`, `AS+B<64`; for `S>=19`,
`3S^2/16-AS-B` is positive at 19 and increasing. These two integer ranges
are exhaustive. The general sector is strictly closed at every `c>=64`.

The `(1,2,2)` circuit has exactly defect zero, slope `3/14`, and intercept
`4/3`, giving a strict count tail from `c=34`. The omitted `(0,1,1)` forces
equal speeds. The remaining omitted circuits are exactly norm three and
norm four; their count slopes `2/7` are genuine hostile controls to using
the count bound alone.

## 3. Physical exceptional lines and the fixed selected roof

For norms three and four, `|(v.n)|<(3/14)||v||_1<1` forces the entire
dictionary onto `C=kv`. Every coefficient is a ternary unit, so the live
multipliers are precisely the strict roof interval with `3` not dividing k.
The physical profile is even and nonincreasing. Comparing its full integer
sum and its three-step sum gives error at most `2f(0)=6/(7c)`. No residue
class or strict endpoint is discarded.

The additive profile has integral `27/196`; this follows independently from
the incoming additive theorem's three explicit pieces. Thus the physical
tail is strict from `c=50`. For norm four, the doubled coefficient's roof
and common cap have total integral at most `9/98`, because
`c(p+q)-pq<=c^2` for its complementary speeds `p,q<=c`. Retaining additional
roofs can only decrease the physical profile. This gives the strict physical
tail from `c=18` and is sufficient without an exact coarea identity.

The selected-network upgrade requires a pointwise statement. The three
sorted norm-four identities are exhaustive: `c=2a+b`, `c=a+2b`, or
`2b=a+c`; the doubled coefficient cannot occur at c. For its complementary
speeds `p,c`, put `H=3/(7c)`, `k0=3(c-p)/28`, `k1=3(c+p)/28`.
I checked every value in the candidate's endpoint table independently.
At k0 the two other roofs are both at least H; at k1 both are nonnegative.
They decrease on the preceding plateau, and their differences from the
doubled roof are affine on `[k0,k1]`. Hence the same fixed capped roof is
minimal on the entire live line. It follows that `E_j=mu(F_w)=min_i E_i`.
This argument does not interchange a minimum with a sum.

The incoming additive selected-network dependency is correctly typed.
Its two relevant continuum profiles cross at `a/c=1/4`, with maximum
minimum `39/392`. The strict multiplier-count layer-cake estimate remains
valid for the downward cutoff jumps and gives error `4/(7c)`. Its network
tail is strict at `c>=60`; its independently audited 136-row finite base
supplies unique equality. The countercontrol `(4,7,11)` has different
physical and selected masses, so the fixed norm-four roof cannot silently
replace this separate additive proof.

## 4. Independent complete physical head

The audit uses the literal predicate
`all j in {0,1,2}: some w has ||w(x+j/3)||<1/14`.
It places all breakpoints on the integer denominator `42*lcm(a,b,c)`,
tests the exact midpoint predicate on each interval of `[0,1/3]`, and
multiplies the resulting length by three. It neither imports the producer
nor enumerates carriers or the six assignment intersections. Endpoints
have measure zero; strict predicate comparisons are nevertheless retained.

Every one of the **10,074** primitive ternary-unit triples through height 63
agrees with the independently produced raw-table physical mass. The parity
counts are 1325 all-odd, 4375 with one even speed, and 4374 with two even
speeds. There are no violations of `6/55`; its only equality is `(1,10,11)`.
There are 151 violations of the old `6/77` ceiling, first `(2,5,7)` in height
order. The next physical level is `383/3640`, uniquely at `(5,8,13)`.
The exact controls `(1,2,4)`, `(2,5,7)`, `(1,7,8)`, `(1,5,11)`, `(1,10,11)`
and `(2,11,20)` all pass. All analytic tails are strict, so this physical
head supplies the global physical equality boundary.

The separate [native interval verifier](../../04-computation/overnight4_20260906_lrc_parityfree_native.cpp)
independently compares every selected projection and physical mass on this
same head, including the contact-to-carrier factor `3N`. Its 181,385 gates
pass, with selected equality only at `(1,10,11)` and projection vector
`(6/55,12/77,23/154)`. Thus the selected-head independence condition is also
discharged; the two proofs do not rely on the same carrier enumeration alone.

## 5. Reproduction and limits

[Audit source](../../04-computation/overnight4_20260906_lrc_parityfree_audit.py)
and [frozen output](overnight4_20260906_lrc_parityfree_audit.out) retain
**101,824** optimization-live gates. Normal and optimized outputs are
byte-identical. The raw TSV option adds rowwise comparison; without it the
same source still independently computes and verifies the entire literal
head, together with the owner fibers, coefficient box and cutoff gates.

```text
python -B 04-computation/overnight4_20260906_lrc_parityfree_audit.py --head-tsv 05-knowledge/results/overnight4_20260906_lrc_parityfree_probe_head63.tsv
python -B -O 04-computation/overnight4_20260906_lrc_parityfree_audit.py --head-tsv 05-knowledge/results/overnight4_20260906_lrc_parityfree_probe_head63.tsv

SHA-256, LF bytes:
audit source b5f52380e28bd3883f95090ff0c06a89e667d7aafd0e4c241b8fca68c76c7c00
audit output 17b0980955e7b2118eae558152f6f36730f17a08d7678fcaf6d031dd677284d2
producer source 1405ce106b74d533c3f8e330fc95855996594913e031e5f80ea323092445d866
raw table c3d33fdd136245aafe512b04963a6eb6f1b5db6f1a572a3e8535ef59d01a09fa
literal head semantic 10dc0483d65ccec0c4999932983c387bcbb5c081a4c47ed3d76d625df2a24c0e
coefficient semantic cf808062354debbefc1d8ead8ad0d10e9da5427cb42b8f083b6af24d0059c87c
```

Common dilation by a ternary unit preserves the literal physical measure
and permutes the sheet labels, so nonprimitive physical equality is exactly
a positive such dilation of `(1,10,11)`. This does not justify applying a
primitive full-lattice address formula directly to an unsaturated scaled
carrier lattice. No body Haar floor, arbitrary entry, synchronization,
or LRC(14) closure is inferred from this ceiling.
