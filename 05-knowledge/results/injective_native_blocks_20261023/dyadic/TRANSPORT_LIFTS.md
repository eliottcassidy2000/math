# Two-sided transport lifts, with the port budget kept separate

23 October 2026. **PROVED:** the finite-ring transport and lifting statements
below. **FINITE-EXACT:** the exhaustive controls in `check.py`. This is an
algebraic completion tool for authenticated native menus; it does not assert
that an arbitrary integer correction has available Hamiltonian-path ports.

## 1. Inheritance and the exact carrier

The inherited [rectangle-bank theorem](../../session44_guarded_native_moves_20261023/rectangle-bank/THEOREM.md)
tests a single source difference against a single target difference. The
least-used sidecar is the *integer* divided residual before reduction: see
[integral lifting connections, section 5](../../automata_arithmetic_geometry_session_20261019/INTEGRAL_LIFTING_CONNECTIONS.md).
The native hostile is a doubled rectangle whose new amplitude needs the
same unique path twice. Linear solvability cannot repair this resource loss.

Let S and T be **disjoint** sets of actual source and target owner nodes,
of sizes s,t >= 2. A signed selected-port change is an integer matrix M,
with positive entries meaning added conductance before multiplication by
the paired-black factor 2. Exact preservation of both margins means

    M 1 = 0,       M^t 1 = 0.                                  (1)

Put R=(e_1-e_s,...,e_(s-1)-e_s) and
Q=(e_1-e_t,...,e_(t-1)-e_t). Every such M has the **unique** form

    M = R C Q^t,       C integer of size (s-1) by (t-1).         (2)

Indeed C is the upper-left block of M, and (1) determines its last row and
column. Thus C retains every entry, not just a total count. The associated
Laplacian change is zero off S union T and has block matrix

    Delta L = -2 [ 0  M ; M^t  0 ].                            (3)

The zero diagonal blocks use the exact zero margins. Sources and targets
that overlap need the original incidence formula, not this block shortcut.

For a marked potential x, write u=Q^t x_T and v=R^t x_S. Since R and Q
contain identity submatrices, preserving x modulo 2^(h+1) is exactly

    C u = 0 modulo 2^h,       C^t v = 0 modulo 2^h.             (4)

Integer current preservation is the same pair of equations over Z. For a
full additive charge signature psi, replace x by each coordinate with its
own modulus. This is the two-sided matrix version of the rectangle test:
every target charge relation must also be valid when viewed from sources.

The source is an authenticated signed port table; the target is its exact
integer matrix C and marked residuals. Preserved predicates are margins
and the stated character equations. Lost information is path existence,
port ownership, shared clock equality, eligibility and capacity. All of
those remain explicit sidecars of a proposed physical move.

## 2. One binary seam is the complete unrestricted correction test

For vectors u in F2^c and v in F2^r, consider

    D -> (D u, D^t v),       D in Mat_(r,c)(F2).                (5)

Given a in F2^r and b in F2^c, a solution exists precisely when

    v.a = u.b,
    u=0 implies a=0,
    v=0 implies b=0.                                          (6)

If u,v are both nonzero, there are exactly 2^((r-1)(c-1)) solutions.
If only u is nonzero there are 2^(r*c-r); if only v is nonzero there are
2^(r*c-c); if both vanish and a=b=0, all 2^(r*c) matrices work.

Proof. Change bases so u and v are their first basis vectors. Then D u
specifies the first column and D^t v specifies the first row. Their common
corner is one shared bit, namely (6); the remaining (r-1)(c-1) entries are
free. Invertible changes of bases preserve the conclusion. The one-zero
and two-zero cases follow by the same column or row argument.

The cheapest hostile is r=c=1, u=v=1, a=0,b=1: the two individual equations
are each solvable, but their shared entry is inconsistent. This is a literal
instance of the repo's shared-coordinate obstruction, rather than an
analogy based on matching dimensions.

## 3. Exact dyadic completion, before capacity restrictions

Suppose C,u,v are integral and (4) holds at level h>=1. Every lift of the
same coefficient residue to the next level is C'=C+2^h D, with binary D.
Put

    a = -(C u)/2^h modulo 2,
    b = -(C^t v)/2^h modulo 2.                                (7)

Division in (7) is performed over the integers. The complete lift test is
(6), with u,v reduced modulo two. In this particular application the first
condition is automatic, because

    v^t(Cu)/2^h = u^t(C^t v)/2^h                              (8)

is an integer identity. Therefore, if both u and v have an odd coordinate,
**every** solution at level h lifts to level h+1, with exactly
2^((r-1)(c-1)) lifts. This iterates through any finite depth. It also gives
compatible 2-adic matrix solutions; it does not imply a bounded integer
matrix or finite port assignment.

If all target differences are even, a must vanish. If all source
differences are even, b must vanish. For example u=v=(2), C=(1) passes
modulo two but has no lift modulo four: a=b=1 while both correction maps
vanish. This is the first genuine valuation obstruction, and it must not
be replaced by a field-rank test on the original integer matrix.

Extra binary character equations on C are unchanged by a correction
2^h D for h>=1. Thus, for the paired-black native order-eight quotient,
one may first impose its four binary character coordinates and then lift
the remaining cyclic coordinate using (7). With multiple independent
cyclic coordinates the correction maps must be stacked; (6) alone is no
longer sufficient. The native full dyadic quotient has only one such
coordinate; its authentication is in [the electrical packet](../electrical/).

## 4. Why a successful lift still needs new actual ports

Take r=1,c=2, u=(1,1), v=(0), and C=(1,1). This preserves (4) modulo two.
There are two lifts modulo four, namely (3,1) and (1,3). An unconstrained
integer lift exists, but there is no lift in the fixed coordinate box
0<=C_ij<=1. The surviving solutions can be represented by (-1,1) and
(1,-1), but that changes orientation and is legal only if the actual
reverse transitions and their ports are available. A residue class is
not an authorization to consume those reverse resources.

This gives a procedural split that is useful for the new native blocks:

1. Test the exact matrix and its next divided residual, not a sum of
   separately rounded rectangle debts.
2. Solve the shared source/target correction equations jointly.
3. Intersect the resulting affine solution set with the literal path,
   shared-owner, capacity and converse constraints.
4. If empty, add genuinely new eligible ports or a new source relation;
   repeating a scarce rectangle does not enlarge that feasible set.

The six-source witness in the sibling search packet takes a stronger route:
its relevant source potential relation and target potential difference are
zero as integers. Its higher scalar debts vanish before lifting. Its full
dyadic character certificate is a separate, authenticated charge calculation.

## 5. Reproduction and boundary

Run `python -B check.py` and `python -O -B check.py` from this directory.
The checker exhausts all vector data and binary correction matrices for
1<=r,c<=3, verifies counts, compares the lifting test with direct lifts,
and retains both the nonunit and resource-capacity hostiles. These are
finite controls for the general proofs above, not native enumeration.
No all-size tournament reconstruction follows from this transport lemma.
