# One added vertex gives a uniform sevenfold lift of a two-count read

**Proved:** an actual nine-vertex quotient and its vertex deletion
give exact factor 7 for *every* tournament substituted at a common
marked port, not just transitive module inflations. Neither coverage
of all primes nor an iterated factor-7 map on an arbitrary witness
is asserted.

## Literal host and proof

Use outgoing row masks on 0,...,8:

    U = [16,37,9,179,38,69,31,375,127].

Mark 8, and let V be U with vertex 4 deleted. For a nonempty
tournament M, substitute M uniformly for vertex 8 in both; all
vertices of M have the same contacts to old outsiders. Deleting
the old vertex 4 from U[M] recovers precisely V[M], so the move
really adds just one vertex and changes none of the other arcs.

Let Q_r(M) count ordered covers by r nonempty paths of M.
In a path of a substitution, compress each maximal M-run to a
visit to the mark. For a fixed sequence of outside vertices and
r occupied slots, the uniform contacts let the r ordered paths
be filled in exactly Q_r(M) ways. Thus the count is the sum of
c_r Q_r(M), where c_r are literal valid compressed word counts.
This is the inherited marked-substitution kernel (THM-3121; see
EXACT_READER_INTERPOLATION.md). The following specific match is
the constructive step:

    V: (c1,c2,c3,c4)     = (49,18,0,0)
    U: (c1,c2,c3,c4,c5)  = (343,126,0,0,0).

The replay independently computes these by a DP on the subset
of outsiders already used, the current tail, and the number
of visits. Starting from an empty state, append a new outsider
if the arc permits, or append a mark if the last symbol can
point to it. Each path of states is a unique legal word. Between
successive visits to the same mark at least two outsiders are
needed in a tournament, so c_r=0 beyond floor((n+1)/2).
This finite computation proves the full coefficient match;
there is no extrapolation in the module size.

Since Q_1=H and Q_2=P_2−2H, with P_2 allowing empty slots,

    H(V[M]) = 13H(M)+18P_2(M) = N,
    H(U[M]) = 7N.

U is strong; substituting any nonempty module into a strong
quotient preserves strong connectivity (to connect two module
vertices, use an outside directed cycle through their quotient
vertex). V has a singleton sink, which is retained in all
parents; this is not a strong-parent claim. The new lift is not
ordinal multiplication by a non-existent seven-path tournament
or a 5h→35h lift from the previous note.

## Exact examples and a restricted target inverse

With a transitive m-vertex module, P_2=2^m and H=1. This gives
18·2^m+13 → 126·2^m+91 for every m>=1, e.g. 157→1099.
Nontransitive modules select further primes:

| Module H | Module P_2 | N | 7N |
|---:|---:|---:|---:|
| 5 | 32 | 641 | 4487 |
| 17 | 224 | 4253 | 29771 |
| 43 | 424 | 8191 | 57337 |

All rows are literal physical extensions independently counted
by subset DP, and primality of each N is tested by trial division.
8191 is the Mersenne prime 2^13−1. The four displayed prime
examples occupy four different residue classes modulo 24;
no statement about prime infinitude follows.

For a prescribed parent target N other than the singleton-module
case 49, find an actual module profile (h,p) with

    h <= N/49,  h ≡ 61 N (mod 72),  p = (N−13h)/18.

Indeed p>=2h; for a module of at least two vertices p≡0 mod4
by complement pairing of induced odd path counts. The inverse
13^(-1) mod72 is 61. Also p>=2^n for a module on n vertices:
every subset gives at least one ordered two-slot cover counted
with its complement. Thus 18·2^n <= N−13. These constraints
give an exact finite target-first test for this reader: enumerate
modules only up to that logarithmic bound and check the physical
profile, not just its congruence. How to cover the remaining prime targets with a bank is still
open and small N already need other templates.

The actionable next step is a bank of proportional marked
coefficient vectors with different linear readouts of (H,P_2),
or of higher cover counts, chosen from N. Matching the first
count alone does not certify that a sevenfold lift survives
module replacement.

Run python3 arbitrary_module_replay.py for the marked coefficients,
all example profiles, literal adjacency substitution, strong-child
checks, primality tests and deletion checks.
