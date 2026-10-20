# A source-adjunction automaton with exact binary cycles and a faithful 2-adic limit

2026-10-20. **PROVED elementary operation-preserving bridge; FINITE-EXACT controls.** No external-priority claim. This extends the repository's source-free/Hasse module by constructing and classifying a specific autonomous binary automaton on formal tournament combinations.

## 1. The actual carrier and operation

Let T_n be the finite set of isomorphism classes of n-vertex tournaments, including the empty tournament at n=0. Put t_n=|T_n|, t_(−1)=0, and

    V_N = ⊕_(0≤n≤N) F2^(T_n),       D_N=Σ_(n=0)^N t_n.

An element is an F2-valued statistic on each size up to N. Let s(A) be the tournament obtained by adding one new universal source to A. Define source pullback h by

    (h f)_n(A) = f_(n+1)(s(A)),       0≤n<N,
    (h f)_N = 0.

Thus h lowers degree. Set U_N=I+h. Since h^(N+1)=0, U_N is invertible and defines a functional graph on **2^(D_N)** binary statistic states.

For the dual W_N=V_N*, identify its basis with tournament classes of size at most N. The transpose S_N=h* is truncated source adjunction:

    S_N[A]=[s(A)]  if |A|<N,       S_N[A]=0  if |A|=N.

Consequently U_N*=I+S_N updates a formal mod-two tournament combination by toggling each present tournament and its source extension. This dual description is an explicit operation; it is not inferred from matching graph counts.

This is not a mutation rule on individual tournaments and not an orientation tournament formed from pairwise comparisons of LFSR states. The states are binary sums/statistics, and the operation is source adjunction. That distinction is the connection's required carrier sidecar.

### Inheritance

`05-knowledge/results/tournament_hasse_frobenius_20261019/PADDING_FROBENIUS_PROOF.md`, especially its source-free core count and integral Taylor decomposition, proves the native relation h D_j a=D_(j−1)a for a source-free primitive. It also proves that source adjunction is a split injection on classes. We use that actual operation; the Jordan decomposition below can also be read directly from source stripping, independently of the induced-subtournament basis.

The Hasse/Frobenius mechanism is the same one in `01-canon/theorems/THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier.md`, now realized on source-adjunction strings. `01-canon/theorems/THM-2043-period14-parity-hasse-jet-completeness.md` is the warning that a faithful local module does not automatically retain the original problem's physical state.

## 2. Complete Jordan and finite-cycle classification

A nonempty tournament has at most one universal source: two sources would each have to beat the other. Source adjunction bijects T_(r−1) with the r-vertex tournaments having a source. Hence the number of source-free cores of order r is

    k_r=t_r−t_(r−1),    including k_0=1.                      (1)

Repeatedly stripping the unique source gives a unique source-free core A and a unique nonnegative source-prefix length j. The corresponding source tower is

    A, s(A), s²(A), … .

At cutoff N, a core of order r gives a string of length N−r+1. The strings partition the tournament basis; there is no claim that arbitrary induced subtournament counts reconstruct a tournament.

### Theorem 1: exact elementary divisors and cycle census

On V_N, h is the direct sum, over source-free cores A of order r≤N, of one nilpotent Jordan block of length N−r+1. Therefore U_N has elementary divisors

    (x+1)^(N−r+1), repeated k_r times, 0≤r≤N.              (2)

Its characteristic polynomial is (x+1)^(D_N), its minimal polynomial is (x+1)^(N+1), and its exact order is

    ord(U_N)=2^ceil(log_2(N+1)).                              (3)

The convention at N=0 gives order 1. For a≥0 put

    R_N(a)=Σ_(r=0)^N k_r min(2^a,N−r+1)
          =Σ_(r=max(0,N−2^a+1))^N t_r.                    (4)

Then

    #Fix(U_N^(2^a))=2^(R_N(a)),                             (5)
    C_N(1)=2^(t_N),
    C_N(2^a)=[2^(R_N(a))−2^(R_N(a−1))]/2^a,   a≥1,         (6)

and there are no cycles of any other length. The same counts hold for U_N*.

**Proof.** On basis indicators along each source tower, h sends the jth indicator to the (j−1)st and kills the core indicator. This proves (2). The empty core supplies a string of length N+1, proving the minimal polynomial. In characteristic two,

    (I+h)^(2^a)=I+h^(2^a).

On a string of length e, ker h^(2^a) has dimension min(e,2^a), proving (4)–(5). The least dyadic power killing the largest string gives (3), with sharpness supplied by that string. The second expression in (4) follows by telescoping k_r=t_r−t_(r−1): only the last min(2^a,N+1) tournament class counts remain. Exact-period states are the difference of successive fixed-point counts; division by the period gives (6). At a=0, R_N(0)=Σ k_r=t_N. Transpose preserves the ranks of every polynomial in U_N, giving the same cycle census. ∎

This is a noncyclic family in general. Its characteristic polynomial alone would lose all string multiplicities in (2); the full elementary divisors are essential. This is precisely the general-linear extension of the companion formulas in `BINARY_LFSR_PRIMARY_CLOCKS.md`.

### Small exact examples

The tournament class counts through order 5 are

    t_0,…,t_5 = 1,1,1,2,4,12.

For N=3, D_N=5, the 32-state graph has 4 fixed points, 2 two-cycles, and 6 four-cycles.

For N=4, D_N=9, the 512-state graph has

    16 fixed points, 24 two-cycles, 48 four-cycles, 32 eight-cycles.

For N=5, D_N=21, the 2,097,152-state graph has

    4,096 fixed points,
    30,720 two-cycles,
    114,688 four-cycles,
    196,608 eight-cycles.

These counts sum with cycle-length weights to 2^(D_N). The N=5 statement is a proved formula with independent exact matrix-rank verification, not a claim that the checker traversed all N=5 states.

## 3. The inverse-limit distinction matters

The inclusions V_N→V_(N+1) commute with h and U. Natural degree truncations V_(N+1)→V_N do **not** commute with h: a discarded top-degree component can contribute to degree N under h. Therefore taking an inverse limit of V_N with these truncations would not define the advertised action.

The dual system repairs this exactly. Degree truncations

    W_(N+1) → W_N

commute with the degree-raising source adjunction S. Let

    W=lim← W_N = ∏_(n≥0) F2[T_n].                           (7)

Its elements are formal tournament series, with finitely many possible basis terms in every degree. Source adjunction S is continuous, injective, and raises degree by one. Put U=I+S. Its inverse is the degree-convergent series

    U^(−1)=I+S+S²+… .                                      (8)

This is a compact characteristic-two module, not the direct sum of finite statistics. In the direct sum V=⋃V_N every state has a finite dyadic orbit. In W, the next theorem shows that no nonzero state has a finite orbit.

## 4. A faithful Z₂ action on every nonzero state

For α∈Z₂ with binary digits α_a∈{0,1}, define

    U^α = ∏_(a≥0) (I+S^(2^a))^(α_a).                     (9)

The product converges in every finite-degree quotient, since sufficiently high factors are the identity there. It agrees with the usual integer powers for nonnegative α, and density gives a continuous action of the additive group Z₂.

### Theorem 2: every nonzero orbit closure is the dyadic odometer

For every nonzero w∈W, the map

    Z₂ → W,     α↦U^α w                                  (10)

is injective and is a homeomorphism onto the closure of the integer orbit of w. Under this homeomorphism, U is addition of 1. Thus every nonzero orbit closure has supernatural clock 2^∞. Zero is the unique finite-period state.

**Proof.** Let γ∈Z₂ be nonzero, and let a=v_2(γ). In the formal series ring F2[[z]], the binary product for (1+z)^γ has first nonconstant term z^(2^a), with coefficient 1. Hence

    (1+z)^γ−1 = z^(2^a) B(z),     B(0)=1.

Substituting the degree-raising S is valid and gives

    U^γ−I = S^(2^a) B(S).                                (11)

The operator B(S) is invertible by its nonzero constant term, and S is injective. Thus U^γ−I is injective, so it cannot fix any nonzero w. The map (10) is injective. It is continuous from compact Z₂ into Hausdorff W, hence a homeomorphism onto its image. The integers are dense in Z₂, so that image is exactly the orbit closure. Any nonzero integer is a nonzero element of Z₂; therefore a finite orbit would contradict (11). ∎

The orbit identification is a pointed torsor identification and an adding-one conjugacy; the orbit is not asserted to be an additive subgroup of the F2-state module. The global phase space W has many orbit closures. It is not itself asserted to be a single odometer, nor are all these states isomorphic as module vectors. The theorem is a pointed-orbit classification with an explicit faithful map.

### The empty tournament gives a literal binary-digit encoder

Let e be the basis vector of the empty tournament. Its source tower consists of the transitive tournaments Tr_j. For every nonnegative integer t,

    U^t e = Σ_(j=0)^t (binom(t,j) mod2) [Tr_j].             (12)

For α∈Z₂ the same formula is interpreted coefficientwise, by finite quotient limits. Lucas's digit rule gives

    coefficient of [Tr_(2^a)] in U^α e = α_a.              (13)

Thus the infinite formal orbit explicitly stores the binary digits of its 2-adic time. This is stronger than matching a sequence of cycle lengths: the observation coordinates invert the clock parameter.

The same Frobenius identity that makes repeated factors produce dyadic periods in an LFSR makes (12) a Pascal/Sierpiński coefficient pattern. No topological or physical identification is inferred from the visual pattern; the source-adjunction operator supplies the actual map.

## 5. Exact tests, hostiles, and next usable consequences

Run `python3 verify_lfsr_tournament_clocks.py` in this packet.

Independent controls:

1. Enumerate all labelled tournaments through order 5 and canonicalize under every vertex permutation. This recovers the class counts without importing their table.
2. Construct the source-adjunction matrix directly on canonical tournament classes, without assuming the Jordan decomposition.
3. For N≤5, square the matrix over F2 and compare all relevant fixed-space dimensions with (4).
4. For N≤4, traverse every binary state and compare the complete cycle census with (6).
5. For each tested N and all time residues through twice the clock length, compare the actual orbit of the empty class with the Pascal/Lucas prediction (12).

The main hostile is built into §3: truncating the wrong graded direction does not commute with the action. Another is the module/graph ambiguity in the companion note: identical cycle counts need not recover source-free coordinates or a linear conjugacy.

A concrete continuation is now available without inventing new graph orientations: any chosen tournament statistic is an output functional on these source towers. Its exact observed period can be smaller than (3), and its minimal output recurrence can be extracted from the corresponding cyclic observable module. To use that in an original tournament problem one must specify which statistic is observed and why its source-adjunction evolution preserves the target predicate. The present theorem supplies the automaton and its clock, not that additional application.
