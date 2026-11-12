# A lattice-exact observation compiler for arbitrarily switched physical blocks

**Session280, 2026-11-12. PROVED conditional on the inherited physical block
law; FINITE-EXACT implementation controls. Chan–Pak remains OPEN.**

This lane removes a spectral assumption by changing what is observed. It does
**not** recover an arbitrary source packet from its one scalar count. Instead,
it gives a lattice-unimodular three-count chart, and an exact scalar recurrence
for every nonperiodic switching schedule. The distinction is load-bearing.

## 1. Inheritance and a distant-source connection

The physical moving-root block transfer is inherited from
[Session279, physical Sections1,8](../../session279_interval_tightening_20261112/physical/PROOF.md),
which credits the O/Z/J grammar and its corrected provenance in S25. Its source
packet is `x=(A,W,d)`, where A is the Hamiltonian-path count, W=o+e and d=o-e
for the root path/complement sums defined there. Actual source packets from
nonempty block chains satisfy `0<=d<=A<=W` and `W=d mod2`.

For an active word w (starts1, ends0), let P be its binary value,
Q=sum(1-w_j)2^(j-1), t=2^|w|, a=P+1, b=Q+1, and h=g(P), with
`g(0)=1; g(2v+1)=g(v); g(2v)=g(v)+v+(v mod2)` for v>=1. Then

    M_w(A,W,d)=(hA+bW-d, aA+tW, A).                   (1)

These are actual tournament operations, not signed formal graph sums.

The seeded random source sampler used seed2807 and sixteen root-level result
files (recorded in `SOURCES.json`). Two useful non-tournament readings were:

* [finite-seed bootstrap boundary](../../finite_seed_bootstrap_boundary_20261005.md),
  Sections2,4: a locally valid recurrence needs its actual terminal owner and
  closure, not merely a contracting formula. This motivates the physical
  initialization gate in Section5.
* [mixed Smith factors](../../third_20260906_smith.md): rationally invertible
  coordinates can retain substantial integral index. This motivates computing
  an exact integer observation chart rather than only a nonsingular matrix.
* [procedural Collatz approach foundry](../../collatz_procgen_20260922_foundry.md):
  change the observation rather than declaring a scalar-observation barrier a
  system barrier. This is methodological inspiration, not an imported theorem.

Concept board: **switching cocycle / physical sensor / integral index / scalar
recurrence / owned initialization**. The nearest hostile is Session279's
non-Pisot active word10101010. The corrected near miss is treating contraction
of ratios as shrinking absolute uncertainty. The least-used sidecar here is
the *two nearby physical experiment counts*, retained as part of the object.
No claim of literature novelty for elementary observability theory is made;
the result is this explicit physical chart and its all-word guarded compiler.

## 2. Equal-order probes give a unimodular integer chart

Starting from the SAME rooted source, append separately the active words100
and110. Both add four vertices by the inherited one-apex operation. Write

    B=H(source followed by100)=5A+7W-d,
    C=H(source followed by110)=5A+5W-d.

Set `Y=(B-A)/2`, `Z=(C-A)/2`. In native integral coordinates `(A,o,e)`,

    (A,Y,Z)^t = T(A,o,e)^t,
    T=[[1,0,0],[2,3,4],[2,2,3]],    det(T)=1.          (2)

Consequently the *exact* reconstruction is

    o=2A+3Y-4Z,   e=-2A-2Y+3Z,
    W=Y-Z,       d=4A+5Y-7Z.                         (3)

**Theorem 1.** Equation(2) is an integral bijection Z^3 -> Z^3. On actual
root packets it is implemented by three Hamiltonian-path counts: the original
count and two equal-order physical probes. No extra congruence is hidden in
this chart beyond making Y,Z integral. It works for every source size, including
sources with arbitrarily large W, and has no spectral or suffix-depth premise.

Proof: the displayed inverse(3) has integer coefficients and multiplication
checks both inverse identities. The first rows of the inherited block matrices
for100 and110 give B,C. Since W=d mod2, both B and C have parity A. ∎

The cone becomes the explicit integer polyhedral guard

    A>=1,  0<=4A+5Y-7Z<=A<=Y-Z,                      (4)

with A odd for tournament counts. Conversely(4) only constructs a **formal
cone packet**. It does not supply a rooted tournament. For example A7,W7,d1
passes the cone/parity tests and produces B83,C69, but cannot have a physical
source because the inherited exclusion of H7 applies. The chart does not
circumvent that exclusion. Nor are B,C themselves promised attainable for an
arbitrary formal input.

The factor-two normalization is essential: the map `(A,o,e)->(A,B,C)` has
absolute determinant4. Treating the unnormalized counts as a unimodular chart
would lose their two parity constraints.

## 3. Universal nonautonomous three-count observation

Let block j have parameters `(h_j,a_j,b_j,t_j)`. Write A_j for the count after
j blocks, and let the initial source packet be `(A_0,W_0,d_0)`. For ANY two
successive active blocks,

    W_0 = [A_2-h_2 A_1+(1-b_2 a_1)A_0]/(b_2 t_1),
    d_0 = h_1 A_0+b_1 W_0-A_1.                       (5)

**Theorem 2.** The map from `(A_0,W_0,d_0)` to the three counts `(A_0,A_1,A_2)`
has determinant `b_2 t_1>0`. Thus these three counts determine the full initial
packet over Q, for arbitrary active block choices, including non-Pisot blocks.
The exact lattice gate is divisibility of the numerator in(5) by `b_2 t_1`,
followed by the cone and parity checks on the recovered packet.

Proof: the observation rows are `e_A, e_A M_1, e_A M_2 M_1`. The last two
columns of the last two rows are

    [[b_1,-1],[h_2 b_1+b_2 t_1,-h_2]],

whose determinant is `b_2 t_1`. Alternatively substitute
`W_1=a_1 A_0+t_1 W_0` and `d_1=A_0` into A_2. ∎

This is a different integral geometry from the parallel probes: chronological
observations have an ambient lattice index `b_2 t_1`, rather than a unimodular
chart after parity normalization. A rational solution without divisibility is
not a valid root packet.

**A physical two-count hostile.** The two block chains

    (100,0,110): (A,W,d)=(145,225,15), order10;
    (100,1,110): (A,W,d)=(145,223,9),  order10

have the same original count145 and the same triangle child count1095. Their
second triangle child counts are7145 and7121, respectively. So the actual
observation pair `(A, H(after triangle))` is insufficient even among equal-order
physical sources, whereas every triple in Theorem2 is sufficient. This is not
a claim that every imaginable pair of nonlinear sensors fails on every
restricted physical subfamily. On the full-dimensional cone lattice any two
linear sensors have a rational kernel and cannot be injective globally.

## 4. An exact scalar recurrence for arbitrary switching

For j>=2, eliminate W_(j-1) from(1). The result is the integer identity

    b_j A_(j+1)
      = (b_j h_(j+1)+b_(j+1)t_j) A_j
        +(b_j b_(j+1)a_j-b_j-b_(j+1)t_j h_j) A_(j-1)
        +b_(j+1)t_j A_(j-2).                         (6)

Unlike a fitted constant-coefficient recurrence, this retains the full switching
schedule. For j>=1 the local reconstruction is

    W_j=(A_(j+1)-h_(j+1) A_j+A_(j-1))/b_(j+1),
    d_j=A_(j-1).                                    (7)

**Theorem 3 (exact all-word realization criterion relative to a seed bank).**
Fix any finite or infinite supplied bank S of actual rooted sources, with each
source's packet and graph construction. Fix any finite active switching word
`w_1,...,w_m`, m>=2, and proposed integers `A_0,...,A_m`. These are exactly the
prefix counts of this switched physical construction from a source in S iff:

1. the initial packet reconstructed by(5) is integral and equals a packet of
   an actual source in S; and
2. recurrence(6) holds at every j=2,...,m-1.

For m=2 condition2 is empty. The same statement applies to an infinite block
schedule and infinite count sequence, with every finite recurrence required.

Proof: necessity is(1),(5),(6). For sufficiency choose the supplied source
with that packet and apply the actual blocks. Its first three counts agree by
(5). Since b_j>0, recurrence(6) uniquely determines each subsequent scalar,
so induction gives equality of every prefix count. This constructs actual
finite tournaments at each prefix; no infinite tournament with a finite
Hamiltonian count is asserted. ∎

Cone, parity and local divisions are useful *rejection filters* during target
search; with a genuine initial seed and exact recurrence they are automatic.
They must not replace that seed gate. Conversely, checking only finitely many
recurrence instances cannot certify an arbitrary infinite proposed sequence;
an all-length symbolic proof or a recurrence-defined sequence is needed.

All denominators in (5) and (7) are known from the schedule. This gives a
nonautonomous scalar compiler with explicit carries, including schedules
containing any number of10101010 blocks. It does not require a common Perron
ray, a joint spectral gap, periodicity, or bounded source scale.

## 5. What changes, and what does not, for Chan–Pak

The arbitrary-switching problem separates into two questions:

* **Observation:** can lost root state be recovered from retained experiments?
  Theorems1–3 answer yes, uniformly, after augmenting the scalar object by two
  physical counts. This survives even a failure of one-count Pisot squeezing.
* **Coverage:** for a prescribed odd N other than7,21, can one choose an actual
  source and schedule with final A_m=N? None of the theorems supplies this.
  It remains a coupled integer reachability problem, not merely recurrence
  fitting. The other two observation counts cannot be assigned arbitrarily.

The constructive benefit is an exact interface for combining methods: a
candidate target corridor can be represented by integer triples `(A,Y,Z)` and
never lose its reusable endpoint state. Different proof lanes can constrain
those coordinates, while the seed and switching compiler give a decisive
physical acceptance test. A global coverage theorem would still need to show
that every target fiber meets this physically initialized reachable set.

## 6. Connection contract and controls

| Source -> target | Preserved | Lost / sidecar | Test |
|---|---|---|---|
| Root packet -> equal-order probe chart | full integer packet, exact counts | graph isomorphism lost; actual source retained | determinant1 and literal graph recount |
| arbitrary active schedule -> scalar recurrence | all prefix counts and reconstructible packet | schedule and first three counts essential | non-Pisot schedules and random schedules |
| local arithmetic triple -> physical chain | only after seed match | cone/parity alone do not imply physical source | formal H7 guard passes but source forbidden |
| finite-seed Collatz caution -> tournament compiler | terminal ownership obligation | no Collatz theorem transported | absent seed must reject |

`verify.py` uses standard-library exact arithmetic, explicit checks that remain
active under `python -O`, and independent subset-DP recounts on literal
constructed tournaments. `results.json` records the finite test universe.
Run it from any working directory. Infinite claims follow from the displayed
matrix identities and induction, not from the finite checks.

## 7. A common dyadic lattice survives every active switch

There is also a genuine arbitrary-switching contraction, but on an **integral
exterior lattice**, not an archimedean cone. It is useful for exact inverse
carry guards and must not be mistaken for real scalar-interval shrinking.

Use exterior coordinates `(A wedge W,A wedge d,W wedge d)`. Direct minors of(1)
give

    E_w=wedgesquare(M_w)
       =[[ht-ab,a,t],[-b,1,0],[-t,0,0]].              (8)

Every active word has h,a,b odd and t divisible by4. Define the common lattice

    L0={(u,v,z) in Z^3: u+v is even and z is even}.

**Theorem 4 (sharp common dyadic carry depth).** For every active block,

    E_w Z^3 subset L0,     E_w L0 subset 2 L0.        (9)

For an arbitrary sequence of k>=1 active blocks of total spine length
`ell=sum |w_j|`, let P=M_(w_k)...M_(w_1). Its Smith invariant factors are

    (1, 2^r, 2^(ell-r)),  k-1 <= r <= floor(ell/2).  (10)

The lower bound is sharp at every k: repeating the block100 gives

    Smith(M_100^k)=(1,2^(k-1),2^(2k+1)).              (11)

Proof of(9): modulo2, E_w has first two rows(1,1,0) and third row0, so its
image lies in L0. If u+v and z are even, each output coordinate is even.
Moreover the sum of its first two coordinates is

    (a+1)(v-bu)+ht*u+t*z,

which is divisible by4; its third coordinate is `-t*u`, also divisible by4.
Thus the output divided by2 lies in L0. Inductively every entry of the product
of k exterior matrices is divisible by2^(k-1), proving the minor-divisibility
part of(10).

All M_w have determinant t and the same matrix modulo2,

    F=[[1,1,1],[1,0,0],[1,0,0]].

F^2 is the all-ones matrix J and FJ=J, so every product has an odd entry. Since
its determinant is a power of2, the gcd of all entries is1. The first Smith
factor is therefore1. The gcd of all two-by-two minors is a power2^r; it is
the second Smith factor because the first is1. Determinant and the divisibility
ordering of Smith factors prove(10). These statements use the standard
integer Smith-minor characterization, not a finite computed pattern.

For sharpness take100: `(h,a,b,t)=(5,5,7,8)`. A basis matrix of L0 is

    U=[[1,0,0],[1,2,0],[0,0,2]].

In this basis the map E_100/2 on L0 is the integral matrix

    D=U^(-1) E_100 U/2
     =[[5,5,8],[-4,-2,-4],[-2,0,0]].                 (12)

Modulo2 it is idempotent with first row(1,1,0) and other rows zero.
Also `U^(-1)E_100(1,0,0)^t=(5,-6,-4)^t`, whose image under every power of D
has odd first coordinate. Hence `E_100^k(1,0,0)/2^(k-1)` has odd first
coordinate for every k>=1. Its entries are not all divisible by2^k, proving
r=k-1. Since ell=3k, the last invariant is2^(2k+1). ∎

**Meaning and boundary.** Every active-switch image has purely dyadic ambient
lattice obstructions, and after k switches two-dimensional minors have at least
k-1 compulsory binary carry levels. Modulo any odd modulus P is invertible.
This does **not** imply that a finite physical seed set occupies all residue
classes, that the cone image is the lattice image, or that count7 is possible.
The lattice here is ambient Z^3 in `(A,W,d)`, rather than the stricter parity
sublattice of physical packets. Its Smith factors are not advertised as the
Smith factors in a different integral basis. Equation(2) retains that parity
sublattice explicitly. Dyadic thinning, real contraction, and physical seed
supply are three separate mechanisms; the first survives arbitrary switching
without assuming either of the latter two.
