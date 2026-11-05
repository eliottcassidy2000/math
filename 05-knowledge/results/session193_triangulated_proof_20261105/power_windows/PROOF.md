# Exact cut responses, binary padding, and 135-value dyadic windows

Session 193, 2026-11-05. **PROVED elementary all-order construction from FINITE-EXACT seeds; independently reviewed.** The full Hamilton-path spectrum conjecture with exceptions `{7,21}` remains **OPEN**. The width below is a fixed 135 values, not a width tending to infinity or an arbitrary-target algorithm.

## 1. Inheritance and the complementary roles

Write H(V) for the number of directed Hamiltonian vertex words of a tournament V, with H(empty)=1. The inherited S192 directed-cut theorem says that adjoining an apex x beating exactly S changes H by the weighted outgoing-cut capacity. Its weights are positive even integers and retain all longer-cycle contributions. The inherited S192 sink-padding identity transports an exact local cut spacing without its cofactor being multiplied by a nontrivial core count.

The arithmetic packet's residue/size selector does something different: it chooses digits with a prescribed congruence and a magnitude error. Its error bound never reaches the modulus, and one-apex tournaments miss a positive density of odd values. Here exact cut seeds replace the approximation step on a specified finite band of offsets. A binary continuation formula shows exactly why the same digit operation cannot make the boundary fibres grow indefinitely.

Closest mechanism: [Session 192 cut-sumset proof](../../session192_surprising_approaches_20261105/cut_sumsets/PROOF.md), section 5.2, gives four values `2^k+1,3,5,7`. Canonical hostile: the H13 core has a holey one-birth support despite gcd one. Corrected near miss: the sum of abstract cut supports is not a physical tournament sum, and residue equality plus a large size error is not equality. Least-used sidecar retained here: the exact ordered two-path response R(S), together with the end count e(S).

Concept board: physical cut support; two-path boundary responses; binary carry; exact same-slope fibres; finite seed coverage versus unbounded induction. No literature-priority claim is made.

## 2. The exact continuation, with the missing carry retained

Let V be any nonempty tournament, S a nonempty subset, h=H(V), and F=H(V+x_S). Define

    R_i(V) = sum_(A subset V, i in A) s_i(V[A]) H(V\A),
    R(S) = sum_(i in S) R_i(V),
    e(S) = sum_(i in S) e_i(V).

Here s_i and e_i count paths starting and ending at i. R_i counts ordered spanning two-path covers whose second path is nonempty and starts at i; the first may be empty. Put R=R(S), e=e(S).

Append a transitive chain y_1 -> ... -> y_m after V, with every V vertex beating every y_j. The apex beats y_j precisely when the j-th letter b_j of a binary word w is 1. Denote its count by F(w). The old apex incidences on V are unchanged.

**Theorem 1 (all binary words).** Reading a letter changes the three boundary integers by

    0: (F,R,e) -> (F+R+e, 2R,   0),
    1: (F,R,e) -> (F,     2R+h, h).                 (1)

The fixed parameter h remains H(V), because appending a sink to the base tournament does not change its Hamilton-path count.

**Proof.** In the base U -> y, every old Hamilton path ends at y, so old end counts vanish and y has end count h. Each two-path cover contributing to old R_i extends in exactly two ways, putting y at the end of its first or second path; an empty first path becomes the singleton y. This gives R_i'=2R_i. A two-path cover whose second path starts at the sink y has that second path equal to y, and the first path is a Hamilton path of U, so R_y'=h. Thus R'=2R+b h and e'=b h.

If b=1, y is a sink of the full tournament, so F'=F. If b=0, the inherited cut-capacity identity gives F'-F=R+e: the capacity from i in U to the new sink is R_i+e_i. For completeness, its exposed backward words y,i are exactly the two-path covers counted by R_i, whereas its valid forward words i,y are Hamilton paths of U ending at i, counted by e_i. Summing over i in S proves (1). This argument applies at every step, including when the base is weak. QED.

If A(w) is the ordinary binary integer read left-to-right from w, then

    R(w) = 2^m R + h A(w).                            (2)

For another useful closed form, let C(w)=sum_(j:b_j=0)2^(j-1), and let H_apex(w) be the one-apex Hamilton count on the transitive chain with attachment word w. For m>=1,

    F(w) = F + R C(w) + e(1-b_1)
                  + h[H_apex(w)-1].                 (3)

This follows either from (1) or from the S192 ordered-join capacities: core-to-chain capacities are `R_i 2^(j-1)+e_i [j=1]`, and chain-to-chain capacities are h times their transitive-base capacities. Formula (3) retains the exact integer carry rather than reducing it modulo a power of two. At m=0 simply F(empty)=F.

## 3. Why binary digit programming cannot enlarge one boundary fibre

Every Hamilton path of V, split just before a fixed i, gives a distinct two-path cover counted by R_i. Therefore R_i>=h, and in particular R(S)>=h>=1. Also 0<=e(S)<=h.

**Theorem 2 (response separation).** For a fixed initial nonempty S, distinct finite binary words have distinct R(w). Moreover, R(w) is a power of two if and only if R is a power of two and every letter of w is zero.

**Proof.** For fixed length m, (2) is injective in A(w). Its entire possible range is

    [2^m R, 2^m(R+h)-h].                              (4)

Since R>=h, the upper endpoint is strictly below `2^(m+1)R`, the minimum at length m+1. Thus all length ranges are disjoint.

Suppose R(w)=2^k. Because R>=1, k>=m. Reducing (2) modulo 2^m gives h A(w)=0 mod 2^m. Rédei parity makes h odd, hence invertible modulo 2^m; as 0<=A(w)<2^m, A(w)=0. Then w is all zeros and 2^m R is a power of two, so R is a power of two. The converse is immediate. QED.

Consequently, branching over binary attachment digits above one fixed core never produces two distinct words with the same continuation response, even if their final scalar H values collide. In particular it cannot create new pure-dyadic slopes or enlarge one pure-dyadic intercept bank. This is an exact boundary of this operation, not a no-go for other tournament constructions, changing the core, using two independently variable apexes, or solving target-dependent equations without preserving the slope.

## 4. Exact dyadic transport and the finite seed theorem

Taking w=0^t, t>=1, in (1) gives

    F(0^t)=F+(2^t-1)R+e.                             (5)

If R=2^q and b=F-R+e, then

    F(0^(k-q))=2^k+b,   k>=q+1.                     (6)

If V is strong and S is nonempty, the padded child is strong: x reaches V through S; V reaches every padding vertex; the last padding vertex reaches x. These paths bring all vertices into one strong component. The subset S may equal all of V. The order is exactly

    k + |V| + 1 - q.                                 (7)

**Theorem 3 (135 consecutive odd values at every dyadic scale).** For every integer k>=13 and every odd b with -115<=b<=153, there is an explicitly constructed strong tournament with

    H=2^k+b.

The chosen tournaments have orders between k-2 and k+2. There is no claim that they have one common order or that these are minimum realizing orders.

**Proof.** `seeds135.json` contains one explicit core V and mask S for each of the 135 offsets. The independently reproducible integer computations in `verify.py` prove for every record that V is strong, S is nonempty, R(S)=2^q, q<=12, and F-R+e=b. Equations (6)-(7) therefore prove the theorem for every k>=13. The selected core orders are 6 through 9. The exact values of `|V|+1-q` in this bank are {-2,-1,0,1,2}. The finite part consists only of these listed small witnesses and their counts, not an exhaustive classification assumption. QED.

Four targeted new seeds which joined the original smaller intervals are:

| b | core outgoing rows | mask S | h | R | e | F |
|---|---|---:|---:|---:|---:|---:|
| -25 | (2,260,33,7,15,347,287,127,157) | 268 | 359 | 2048 | 92 | 1931 |
| 27 | (256,257,11,51,295,455,31,95,204) | 130 | 711 | 4096 | 226 | 3897 |
| 79 | (2,40,259,21,263,285,319,127,139) | 5 | 371 | 1024 | 206 | 897 |
| 99 | (2,4,257,7,15,287,63,127,219) | 8 | 99 | 256 | 14 | 341 |

The mask convention is that bit i means x beats core vertex i. The base chain is appended after the listed core, and x loses to every chain vertex.

## 5. What was actually triangulated, and what remains

The geometric/cut method supplies exact discrete local supports and a physical attachment mask. The two-path interface supplies the correct slope, which scalar H alone loses. The dyadic digit method provides the expanding scale. We combine them by matching actual slopes R=2^q and exact intercepts b, not by putting a residue estimate beside an unrelated size estimate.

Different cores and boundary slices are essential: their union escapes the S192 cofactor obstruction for one repeated boundary-balanced join slice. But a fixed bank still has fixed width. Theorem 2 shows why simply appending the residue-programming digits cannot turn this theorem into all offsets. A completion would require a new core-changing operation that supplies an unbounded intercept interval (with controlled q), or a different target-dependent exact construction. No such operation is proved here.

The dyadic windows alone have density zero as a set of positive integers. They do not prove a positive-density spectrum theorem, cover all primes, imply eventual interval coverage, or settle `{7,21}`. The one-apex positive missing-density theorem and H13 obstruction remain valid; the cores used here need not be one-apex.

## 6. Reproduction and scope

Run from this directory:

    python3 verify.py > receipt135.json
    python3 -O verify.py > receipt135.optimized.json
    cmp receipt135.json receipt135.optimized.json

All checks use explicit exceptions. The verifier computes directed Hamilton paths by subset dynamic programming, the two-path response by the literal subset convolution, and every padded witness at k=13 as an actual tournament with 11 through 15 vertices. It also checks 90 physically built binary-word continuations and nine independent literal-permutation controls. The universal k-quantifier rests on (5), not finite extrapolation.

Exploration source: `generate_all.cpp` is adapted with attribution from S191's class generator; `probe_bank.cpp`, `scan_seeds.cpp`, and `target_nine.cpp` search physical attachments. The final theorem does not use the completeness of any generated bank. `seeds135.json`, this proof, and the standalone verifier suffice. Source and result hashes are retained separately.
