# A residue selector that retains the physical return

**PROVED IN NOTES:** an exact, digit-by-digit residue compiler inside the
primitive module-repair family. **FINITE-EXACT:** small complete residue
tables and joint 7-adic controls. **OPEN:** selecting an arbitrary absolute
Hamiltonian path count, or making the resulting cofactor prime.

This is a response to the concurrent S313 count charts and length bypass,
not a claim that modular completeness implies count completeness. The
source is the physically composable family in [selection](selection/PROOF.md)
and [HANDOFF](HANDOFF.md):

    M(n,m) = D^(ordinal n) → TT_m,
    h = 11^n, p = 100^n 2^m,
    Y(n,m) = 15·100^n 2^m − 11^n,

where D has (H,Q(2))=(11,100). Fix any n≥1. Both strong repair witnesses
exist, at orders 5n+m+5 and 5n+m+6, with the common primitive rooted state

    (Y,5h,5h,Y+6h).

The neutral repair, retained root, R-growth, and common fresh seal are the
ones proved in the selection lane. The residue algorithm therefore
selects actual framed tournaments, not just arithmetic labels.

## 1. Complete unit coverage at every 11-adic precision

For k≥1 let L_k=10·11^(k−1). As m runs through 0≤m<L_k,
Y(n,m) runs through every unit modulo 11^k exactly once.

Proof: 2 has order 10 modulo 11, and

    2^10 = 1+93·11,       93≢0 (mod 11).

The binomial theorem inductively gives

    2^(10·11^(j−1)) ≡ 1+93·11^j (mod 11^(j+1)),  j≥1.       (1)

The order therefore grows by exactly 11 at each precision and is L_k.
Multiplication by 15·100^n permutes units; subtracting 11^n also permutes
units because n≥1. This proves the assertion, with its exact quantifiers.
For n=0 the image instead omits the residue −1 modulo 11; that boundary
must not be silently included.

For a prescribed unit a modulo 11^k, find m in {0,…,9} at precision one.
If m is correct modulo 11^j, replace it by

    m + t·10·11^(j−1),       0≤t≤10.

By (1), exactly one t gives the next target digit. The coefficient of t
after dividing the change in Y by 11^j is

    93·15·100^n·2^m (mod 11),

which is a unit. This gives a constructive inverse using k−1 small
digit lifts, not a search over all L_k exponents. Arithmetic can be done
modulo 11^(j+1); constructing the tournament still requires its actual
number of vertices, potentially exponential in k.

## 2. Simultaneous exact factor 7

The selection lane proves

    7 | Y(n,m)  iff  m≡n (mod 3).

Writing m=n+3b, its exact valuation-one condition is

    v_7(Y(n,m))=1  iff  2+5n+b≢0 (mod 7).                  (2)

Adding L_k to m preserves the entire chosen 11-adic residue. Since
3∤L_k, one of m,m+L_k,m+2L_k satisfies m≡n modulo 3.
From that representative, adding 3L_k preserves both conditions and
changes b by the unit L_k modulo 7. At most one of the first two choices
fails (2). Thus every unit a modulo 11^k has a constructive representative

    Y(n,m)≡a (mod 11^k),   v_7(Y(n,m))=1,   0≤m<6L_k.      (3)

Equivalently, for every unit c modulo 11^k we can construct a count
Y=7u with 7∤u and u≡c modulo 11^k: apply (3) to a=7c.
This does not imply that u is prime. The response state and repair
interface remain available throughout this selection.

## 3. An actual-count valuation chart, and its limitation

For distinct m,m′, cancellation of the common offset gives

    v_11(Y(n,m)−Y(n,m′)) = 0,                  if 10∤(m−m′),
                            1+v_11((m−m′)/10), otherwise.   (4)

Equation (1), followed by factoring a unit exponent, proves (4).
Within one exponent class modulo 10, base-11 common ancestry of the
quotient exponent is therefore recovered from actual count differences.
The next leading unit is multiplied by the explicit unit in (1); it is
not lost by passing to the valuation. With increasing m, the columns
(1,Y(n,m)) have strictly positive ordered 2×2 minors. This is a physical
count version of the address/valuation distinction, although it supplies
an 11-ary hierarchy and does not by itself encode a C3 sibling action.

For fixed n, only O(log X) members satisfy Y≤X. Even allowing both n,m,
all values obey Y≡4 modulo 5. The return family thus gives complete
11-adic unit coverage, joint exact factor-7 control, and reusable ports,
while failing global absolute coverage. It makes the remaining question
more precise: find a new physical digit or module replacement that
changes this absolute range without losing the retained response data.

## Reproduction and provenance

Run `python3 -B 05-knowledge/results/zeckendorf_target_selection_20261113/residue_selector.py`.
It checks all unit targets through 11^3 for n=1,2,3, both direct modular
permutation tables and digit-lift reconstruction, the joint condition
(3), and valuation/leading-unit controls. It does not enumerate large
tournaments; their realization is the previously proved module identity.

The motivating incoming files are S313
`TERNARY_COUNT_CHART.md`, `LENGTH_BYPASS.md`, and `MODULAR_BYPASS.md`
in `session313_zeckendorf_target_selection_20261113` (use the exact paths
in the root incoming-source table). Those results already distinguish
modular coverage, symbolic rank, and absolute count coverage. Our added
point is an explicit simultaneous residue compiler in this retained-port
repair family. No new general theorem about primitive roots is claimed.
