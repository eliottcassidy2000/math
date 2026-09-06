# Independent literal-interval verification of the complete network head

**Status: FINITE-EXACT + INDEPENDENTLY AUDITED by root, 2026-09-06.**
Every eligible triple through height 601 was evaluated directly from the
native sheet intervals. The root independently reviewed the cursor
construction, strict contact test, simultaneous advancement at tied ends,
and the three native translated contacts per raw carrier. This note reserves
no theorem ID and proves no all-height or LRC(14) statement by itself.

## Verified universe and consequence

The universe is exactly

```text
1<=a<b<c<=601,
a,b,c odd and nonzero modulo three,
gcd(a,b,c)=1.
```

There are **1,317,935** triples. Every one is evaluated before any result
classification; no support, direction, coefficient, or count filter prunes
the native computation. All comparisons are exact integers on denominator
`D=42abc`.

The complete result is

```text
min(E_1,E_2,E_3)<=6/77 on every row,
equality iff (a,b,c)=(1,5,11).
```

A signed norm-four relation with absolute coefficients `(1,1,2)` exists
exactly when one of `c=2a+b`, `c=a+2b`, `2b=a+c` holds. This arithmetic
classification has 9,201 rows. On all **1,308,734** remaining rows,

```text
E_i<6/77 for each i=1,2,3.
```

There are 36 empty rows and 147,048,282 native positive contacts in total.
The non-norm-four maximum of any projection is `12/161` at `(1,19,23)`;
the maximum minimum projection in that class is `12/301` at `(1,11,43)`.
The latter row has all three projections equal to `12/301`.

The non-norm-four class includes one-direction rows. It differs from the
root's raw-congruence H601 census restricted to 1,255,029 multi-direction
rows, whose maxima therefore need not be the same. This verifier does not
classify directions or import the raw carrier formula. The separate
short-relation argument supplies the implication that a multi-direction
row cannot have a norm-four relation. The computation independently checks
the stronger every-projection bound on the entire non-norm-four class.

## Independent native construction

The definitions come from
[the earlier literal sheet implementation](../../04-computation/lrc14_one_ray_overnight_hexagon_sep05.py).
For speed `w_i`, sheet label `s in {0,1,2}`, and multiplier `m_i=abc/w_i`,
the interval centers in `D`-scaled coordinates are

```text
14*((3*owner-w_i*s) mod(3*w_i))*m_i,
0<=owner<w_i,
```

with radius `3m_i`, split at the endpoints of `[0,D]` when necessary.
The sorted center residues are simply `r0,r0+3,...,r0+3(w_i-1)`, where
`r0=(-w_i*s) mod3`. This produces a direct cursor without building or sorting
owner lists. Label zero has two half intervals at the ends; the other labels
have ordinary intervals throughout.

All **six permutations** of the three distinct sheet labels are scanned.
For three current intervals `I_0,I_1,I_2`, a contact is present iff
`max(left_i)<min(right_i)`. At each positive contact the contribution to
projection `i` is

```text
min(length(I_j intersect I_k), length(I_i)), {i,j,k}={0,1,2}.
```

The physical mass contribution is the length of the triple intersection.
The cursor or cursors with the smallest right endpoint then advance
simultaneously. This is equivalent to intersecting each pair first and then
scanning the omitted intervals, but computes all three projections together.
No floating-point, approximate endpoint, or symmetry reduction is used.

Each raw carrier has three translated native contacts; the output therefore
labels its count `native_contacts`. These counts are not primitive direction
counts, and they are not substituted into a raw-cardinality theorem.

## Controls and independent replays

Eight named controls check all three exact projections, physical mass, and
the factor-three contact correspondence. They include the empty-core hostile
`(19,23,29)`, the sharp equality `(1,5,11)`, and three wide rows beyond the
proof head, reaching `(1,1201,1205)`. The norm-four control `(1,5,7)` has
`E_2=6/49>6/77`; it confirms why the universal assertion concerns the minimum
while the every-projection assertion excludes norm-four relations.

All **2,910** eligible rows through height 79 were separately replayed from
the optional TSV dump against both the older Python pair-then-third native
interval implementation and its raw carrier/projection implementation.
All three projections and physical mass matched exactly, with 48,804 Python
explicit checks. Native contact counts were three times the raw counts on
every row. The height-79 optimized C++ output also byte-matched an
AddressSanitizer plus UndefinedBehaviorSanitizer build, with no sanitizer
failure. This supplements, rather than replaces, the complete native H601 run.

## Reproduction and integrity

[C++ source](../../04-computation/lrc14_universal_literal_empty_core_sep06.cpp)
and [complete saved H601 output](lrc14_universal_literal_empty_core_sep06.out).

```bash
c++ -std=c++17 -O3 -DNDEBUG 04-computation/lrc14_universal_literal_empty_core_sep06.cpp -o /tmp/lrc14_literal_sep06
/tmp/lrc14_literal_sep06 601
/tmp/lrc14_literal_sep06 79 /tmp/lrc14_literal_sep06_h79.tsv
c++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer 04-computation/lrc14_universal_literal_empty_core_sep06.cpp -o /tmp/lrc14_literal_sep06_sanitized
/tmp/lrc14_literal_sep06_sanitized 79
```

The optional TSV contains speeds, denominator, each projection numerator,
physical mass numerator, and native contact count for every row. The
reported semantic FNV-1a value `6777675351403439901` is a deterministic
**noncryptographic** checksum of those H601 integers in traversal order.
All proof-relevant checks use explicit failures and remain active under
`-DNDEBUG`. SHA-256 of raw LF artifact bytes:

```text
source 64f2d209f5d93cb8f92250b8565c51a2a5437a61d7294f74bec9a712a502e5d0
output 49e334421cc50bc8768bbcc7d6457e63e86b89ddb106b27be4689690a99033b9
```

The finite output can discharge the stated network head when combined with
the separately proved analytic tail and norm-four arguments. It does not
by itself establish those arguments, chart entry, synchronization, or LRC(14).
