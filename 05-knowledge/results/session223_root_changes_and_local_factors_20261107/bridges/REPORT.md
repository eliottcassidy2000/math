# Session223: physical cancellation, the product diagonal, and square reflections

Status: elementary PROVED deductions below, inherited identities explicitly credited; FINITE-EXACT checks in `verify.py`. Independent audit PASS after a wording correction; see `../audit/BRIDGE_AUDIT.md`. No literature-priority claim or Chan–Pak completion.

## Inheritance and scope

Source checkout: `59ff4cd82c66e89cc7abf47483956ae0a7e5fc51`.
Read `05-knowledge/results/root_change_carmichael_chamber_20261107/PROOF.md` and its `audit/AUDIT.md`. The corrected chamber theorem is THM-4904 (`01-canon/theorems/THM-4904-root-change-word-fibers-carmichael-chamber-parity-law-sqrt-density.md`). In particular M_7(1897)=1 prevents extrapolating the two-root-change no-singleton conjecture to all K. Its section 6 already gives the diagonal `(2^a+1)^2`; that is not new here.

Relevant other sources:
- `05-knowledge/results/session221_computation_to_construction_20261107/physical/PROOF.md`: every even arc-flip difference is constructible, but the absolute H sidecar remains essential.
- The same packet `automata/PROOF.md`: fixed-apex compression preserves completed H, not continuation equivalence.
- `05-knowledge/results/session222_tiling_recursions_20261107/README.md`: non-H-preserving L/diagonal identification and separate doubling invariants.
- THM-4851, `01-canon/theorems/THM-4851-the-prehistory-of-the-square-sum-numbers-unique-cycles-at-32-and-33-terminal-spectrum-with-the-hole-at-39.md`: ordinary Q_N has vertices 1,...,N and square-sum edges; Hamilton paths exactly for N in {1,15,16,17,23} union [25,infinity), cycles for N>=32, cited/reproduced. This is not an unresolved conjecture in this repository.
- THM-4860, `01-canon/theorems/THM-4860-the-prime-square-sum-graph-is-bipartite-by-residue-mod-4-two-primes-1-mod-3-are-never-adjacent-chebyshev-bias-obstruction-to-prime-square-sum-paths-none-below-615731-and-the-25-primes-below-100-by-ogg-layer.md`: a DIFFERENT graph on primes, including the open question whether every prime has a prime square-sum partner.
- THM-4861, `01-canon/theorems/THM-4861-the-five-fermat-primes-are-hamiltonian-path-counts-of-tournaments-with-least-orders-3-4-6-8-11-a-vertex-transitive-tournament-has-a-prime-count-only-for-c3-and-fermat-primes-are-free-vertices-of-the-prime-square-sum-graph.md`: known Fermat count constructions and prime square-sum residue obstructions.

Concept board: physical vertex birth, rooted state versus scalar readout, multiplicative diagonal, affine reflection composition, prime predicate, Hamiltonian occupancy. Closest mechanism: JZ=OJ and dagger symmetry. Hostile: equal H need not be equal continuation. Least-used sidecar: square-root edge labels plus all visited vertices.

## 1. Physical redistribution, not scalar cancellation

Let s=(D,E,U,d), t=(1,0,2,1), F(s)=D-E+U, with physical moves

O(s)=(D,0,2U,d), Z(s)=(D-E+U,d,2U-d,d),
J(s)=(D-E+U,D-E,D+U,D).

As inherited from Session205 and quoted in THM-4904 section 3,

    s.JZ = s.OJ = (D+2U,D,D+2U,D).

The intermediate D histories are respectively

    D -> D-E+U -> D+2U,
    D -> D -> D+2U.

Thus the increments U-E and U+E in the first program replace 0 and 2U in the second. There is no negative birth increment; two nonidentical histories redistribute growth and preserve the entire final rooted state. Consequently every common suffix agrees, not just the final H. Iterating gives J Z^r = O^r J for every r>=0. This is an inherited identity's elementary block consequence, not a novel relation claim.

## 2. A complete two-run multiplication rectangle

**Proposition.** For integers a,b>=1,

    t.Z^a O^b = (2^a+1,0,2^b(2^a+1),1),
    F(t.Z^a O^b)=(2^a+1)(2^b+1).

Proof: induction gives t.Z^a=(2^a+1,1,2^a+1,1) for a>=1. The first O clears E and each O doubles U. Substitution proves the formulas.

This extends the already recorded chamber diagonal to its full product rectangle. The final Z or J gives a literal tournament with that product as H, by the physical grammar. Swapping a,b preserves F, as also follows from the inherited dagger. But it does NOT preserve rooted state except on a=b.

A decisive quantitative boundary is

    F(t.Z^a O^b O)-F(t.Z^b O^a O)=2^b-2^a.

For example (a,b)=(1,2) gives equal F=15, states (3,0,12,1) and (5,0,10,1), followed by F=27 versus 25 after O. Thus the scalar equality cannot be promoted to an arbitrary-suffix rewrite. The product map loses the two factor labels and the U coordinate; retaining the rooted state repairs continuation.

The square on a=b arises from equal factors, not from an additive square-sum edge. No primality is used in any of these formulas.

## 3. Exact square-sum graph on the generalized Fermat factors

Put r_a=2^a+1 for a>=1. On the infinite set of these distinct integers, the entire square-sum graph has exactly ONE edge:

    r_1=3 -- r_5=33, because 3+33=36.

Proof: if a,b>=2 then r_a+r_b=2 (mod 4), not a square. If one exponent is 1, the equation becomes s^2=2^b+4. For b=1 it gives 6, not square. For b>=2, s=2u and (u-1)(u+1)=2^(b-2). The case b=2 gives u^2=2, impossible. For b>=3 both positive factors must be powers of two, differ by 2, and hence are 2 and 4. Thus u=3 and b=5. (The boundary u=1 gives zero, also impossible.)

If r_a is prime, a must be a power of two: write a=2^m q with q odd; q>1 factors x^q+1 at x=2^(2^m). The converse is false: 2^32+1=641*6700417. Thus the Fermat numbers F_m=2^(2^m)+1 are only a subfamily, and Fermat primes an even smaller subfamily. Since 33 is composite, the prime members of {r_a} have NO edges among themselves. This sharpens the map boundary, not the known theorem that no two Fermat numbers are square-sum adjacent.

Connection contract: source = two-run physical root programs; target = integer factors r_a and their square-sum graph; map = read the two run lengths; preserved = factorization of the Hamiltonian count; lost = graph incidences and rooted continuation; sidecars = a,b and the full state. The cheapest hostile is a=b>=2, where H=r_a^2 but r_a+r_a is not square. The only off-diagonal square-sum edge uses a composite factor and does not enter the prime graph.

## 4. The real square-sum cancellation law

An edge x+y=s^2 is exactly the reflection y=R_s(x)=s^2-x. For a sequence of positive integer roots s_1,...,s_k, direct induction gives

    R_(s_k)...R_(s_1)(x)
       =(-1)^k x + sum_(i=1)^k (-1)^(k-i) s_i^2.

For even k, this composition fixes a starting value iff the alternating square sum is zero, in which case it is the identity on every integer. For odd k, it has exactly one possible fixed starting value, half the alternating sum, which must be integral and in range. These are precise cancellation constraints on local moves, analogous in role to equality of rooted transfer products, not an identification of carriers.

The missing sidecars are indispensable: every intermediate label must lie in 1..N (or be prime in the prime graph), consecutive labels must be distinct, and a Hamiltonian cycle must visit every required vertex exactly once. Cancellation alone gives none of these. Roots (2,2,2,2) cancel but starting at 1 gives 1,3,1,3,1, not a simple cycle. This is the same logical obstruction as response coverage versus absolute H coverage: an exact algebra of transformations is not an occupancy theorem.

For four distinct vertices, replacing edges (a,b),(c,d) by (a,c),(b,d) preserves every vertex degree. A legal simple square-sum switch additionally requires removed edges to exist and inserted edges to be permitted and absent from the remaining graph. In root labels this forces a+b=s^2,c+d=t^2,a+c=u^2,b+d=v^2 and therefore s^2+t^2=u^2+v^2. The equality is necessary, not sufficient without the endpoint equations and distinctness. Even a legal switch may split one cycle into two; connectivity is the sidecar beyond degrees. This explains why the existing square-sum constructions retain terminal interfaces rather than only square identities.

## 5. What this contributes toward Chan–Pak

The robust useful operation is a relation preserving full transfer state (JZ=OJ), which can be inserted inside any construction. Equal scalar outputs from paired non-H-preserving programs are weaker: the explicit product rectangle supplies a family where the exact next-step defect is nonzero. The Fermat-shaped chamber is a multiplicative coincidence of equal factors with a precise derivation, not evidence for an additive square-sum correspondence. A viable next proof must control the retained rooted/endpoints sidecars and prove occupancy of target H values. Neither finite census coverage, the product rectangle, nor square-reflection closure does that.

## Reproduction

Run `python3 verify.py` and `python3 -O verify.py`; outputs are identical to VERIFY.json. Universes: 3,280 reachable states through word length 7, testing JZ=OJ and block exponents 0..5; 1,600 product pairs 1<=a,b<=40 including next-step defects; all distinct exponent pairs through 200 for the exact one-edge graph; 16,380 reflection compositions at root lengths 1..6, roots 1..4, starting labels -3,0,7; explicit closure-not-simple hostile. Every check uses explicit exceptions and remains active under optimization. These finite checks support, not replace, the proofs.
