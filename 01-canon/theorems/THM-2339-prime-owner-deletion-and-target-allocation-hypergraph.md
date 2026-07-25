---
id: THM-2339
title: "Prime-owner deletion law and composite-target allocation hypergraph"
status: >
  CLAIMED + VERIFIED-EXACT, PROOF CANDIDATE UNDER INDEPENDENT AUDIT. The
  candidate derives every composite-target lift as a finite min-plus
  allocation of its labelled prime-summand occurrences among the current
  partition blocks. Boolean Möbius inversion of each block's target-subset
  score gives canonical, blockwise-unique coefficients relative to the
  supplied score tables and a signed coloured hypergraph energy; no global
  uniqueness after summing over colourings is claimed. Its full tied minimizer
  set, rather than a tournament, is the lossless owner object. A subset
  dynamic program evaluates the lift from a supplied Gordian-distance table.
  If a partition block is itself a nontrivial prime target, deleting that
  owner identifies the whole owner-preserving target face, including every
  merge drop, with the unknot face of the complementary packet. The first
  cover which merges the owner has an exact collision formula in terms of
  translation gain, bypass, and the competing blocks' descent excesses. A
  four-dimensional conical word-metric pair has identical root packet data,
  identical unary prime-target costs, and identical unpartitioned target
  distance, but different pair Möbius coupling and composite lift
  obstruction. No new Gordian distance, unknotting number, catalyst, or knot
  classification is asserted.
source: codex-2026-07-25-prime-owner-deletion
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2330-partition-lattice-gordian-lift-spectrum
  - THM-2336-prime-target-gordian-owner-diagram-and-bypass-split
related:
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
script: 04-computation/prime_owner_allocation_hypergraph_thm2339.py
output: 05-knowledge/results/prime_owner_allocation_hypergraph_thm2339.out
script_sha256: 5e1a6726c1b4976f3af7fabdd136260ca7591f59d3bde44b540ef964e794dbf3
output_sha256: 58dbac059b0062375c1153cfdcebcf5fddcd26314a1147605139100d17570b9e
hash_basis: working-tree bytes (LF)
---

# THM-2339 -- prime-owner deletion and target allocation

**CLAIMED + VERIFIED-EXACT; PROOF CANDIDATE UNDER INDEPENDENT AUDIT.**

This file contains a complete proof candidate and an exact companion, but it
is not yet a proved dependency.

THM-2336 evaluates THM-2330's lift spectrum when the target is prime. There
is then one prime token, so every endpoint assigns it to a single block;
several minimizing blocks may tie. For a composite target, the missing
object is not another binary relation on the blocks: it is the finite
allocation of all prime occurrences, with genuine within-block couplings
between target tokens. Boolean Möbius inversion exposes those couplings
exactly.

The prime boundary has an additional rigidity. If a block already equals
the prime target, deleting that owner turns the whole face on which it
survives into the root partition face of the complementary packet. The
first merge which destroys the distinguished owner has a separate exact
law. It distinguishes translation contraction from bypass and records when
ownership stays with the merged block, ties, or migrates to an untouched
block.

## 1. Setup and the finite prime-token fibre

Let `I` be a nonempty finite label set, let

```text
x=(K_i)_(i in I)
```

be a packet of oriented knots, and let `pi` be a partition of `I`. For every
block `B in pi`, put

```text
K_B=#_(i in B)K_i.
```

Retain THM-2330's target fibre, lift cost, obstruction, and merge drop:

```text
F_pi(J)
 ={(L_B)_(B in pi):#_(B in pi)L_B=J},

Lambda_x(pi;J)
 =min_(L in F_pi(J))
    sum_(B in pi)d_G(K_B,L_B),

Omega_x(pi;J)
 =Lambda_x(pi;J)-d_G(K_I,J),

c_x^J(pi,rho)
 =Lambda_x(pi;J)-Lambda_x(rho;J)       for pi<=rho. (1)
```

Choose a labelled list of the nontrivial oriented prime occurrences in the
Schubert decomposition of `J`:

```text
J=#_(s in S)P_s.                                    (2)
```

The occurrence set `S` is allowed to be empty, in which case `J=U`. For
`A subset S`, write

```text
J_A=#_(s in A)P_s,             J_emptyset=U.         (3)
```

For each block define its normalized target-subset score

```text
w_B(A)=d_G(K_B,J_A)-u(K_B).                         (4)
```

Thus

```text
w_B(emptyset)=0.                                    (5)
```

An allocation is a map

```text
f:S->pi,

A_B(f)=f^(-1)(B).                                   (6)
```

Its endpoint and normalized energy are

```text
L_B(f)=J_(A_B(f)),

E_pi(f)=sum_(B in pi)w_B(A_B(f)).                   (7)
```

> **Prime-token allocation theorem.**
>
> ```text
> Lambda_x(pi;J)
>  =sum_(B in pi)u(K_B)+min_(f:S->pi)E_pi(f).       (8)
> ```

### Proof

Schubert uniqueness says that every factorization of `J` distributes each
prime occurrence in (2) to one of the labelled blocks. Therefore every
endpoint in `F_pi(J)` is `L(f)` for some allocation `f`. Conversely every
allocation gives an endpoint in the fibre. For that endpoint,

```text
sum_(B in pi)d_G(K_B,L_B(f))
 =sum_(B in pi)[u(K_B)+w_B(A_B(f))]

 =sum_(B in pi)u(K_B)+E_pi(f).                      (9)
```

Minimizing proves (8). QED.

If several occurrences `P_s` are isomorphic as oriented primes, different
labelled allocations can give the same endpoint. Let

```text
G_J=product_P Sym({s:P_s isomorphic to P})          (10)
```

act by permuting equal-prime occurrences. The actual endpoint
factorizations are exactly the `G_J`-orbits of allocations. Formula (8) is
unchanged by this overcounting.

Define the full tied optimum

```text
Opt_pi(J)=argmin_(f:S->pi)E_pi(f).                  (11)
```

The lossless endpoint-owner object is `Opt_pi(J)/G_J`, together with the
energy values if more than the optimum is to be retained. No arbitrary
tie-breaking is made.

The boundary cases recover the preceding canon.

- If `S` is empty, (8) is THM-2330's root formula
  `Lambda_x(pi;U)=sum_B u(K_B)`.
- If `|S|=1`, an allocation chooses one block for the unique prime token,
  and (8) is exactly THM-2336's prime-target owner formula.
- If `|S|>=2`, different prime tokens may have different owners. In
  general there is no distinguished single owner block.

## 2. Exact subset dynamic program

Enumerate the blocks as

```text
pi={B_1,...,B_m}.
```

For `A subset S`, define

```text
F_0(emptyset)=0,
F_0(A)=+infinity                     for A!=emptyset,

F_j(A)
 =min_(T subset A)
    [F_(j-1)(A\T)+w_(B_j)(T)]        for 1<=j<=m.   (12)
```

> **Allocation recursion.**
>
> ```text
> min_(f:S->pi)E_pi(f)=F_m(S),
>
> Lambda_x(pi;J)=sum_(B in pi)u(K_B)+F_m(S).        (13)
> ```

Indeed `F_j(A)` is the minimum energy for distributing exactly the tokens
of `A` among the first `j` blocks. The last block receives some subset
`T`; the preceding blocks receive `A\T`. This gives (12) by induction and
(13) at `j=m`.

Given the table

```text
{d_G(K_B,J_A):B in pi, A subset S},
```

the recurrence uses `O(m 3^|S|)` min-plus comparisons and `O(2^|S|)`
working memory. Storing all minimizing subsets `T`, rather than one, recovers
every labelled optimum in (11). The recurrence evaluates supplied distances;
it does not compute a new Gordian distance oracle.

For repeated prime types one may instead use a multiset state recording, for
each oriented prime type `P`, how many of its occurrences have been assigned.
The labelled recurrence is simpler and has exactly the same minimum.

## 3. Boolean Möbius inversion gives the allocation hypergraph

For every block `B` and nonempty `T subset S`, define the Boolean Möbius
coefficient

```text
mu_B(T)
 =sum_(A subset T)(-1)^(|T|-|A|)w_B(A).             (14)
```

Because `w_B(emptyset)=0`, Möbius inversion gives

```text
w_B(A)=sum_(nonempty T subset A)mu_B(T).            (15)
```

Substituting (15) in (7) yields the exact colouring energy

```text
E_pi(f)
 =sum_(nonempty T subset S)
    sum_(B in pi)
      1_[f(s)=B for every s in T] mu_B(T).          (16)
```

Equivalently, if `f` is constant on `T`, let `f(T)` denote its common value.
Then

```text
E_pi(f)
 =sum_(nonempty T subset S, f|_T constant)
    mu_(f(T))(T).                                   (17)
```

Equations (14)--(17) give a canonical signed coloured hypergraph
representation relative to the supplied block tables `{w_B}`:

- the vertices are the labelled prime occurrences `S`;
- the colours are the current blocks `B in pi`;
- the hyperedge `T` has colour-dependent weight `mu_B(T)`;
- it contributes exactly when all tokens in `T` are assigned to `B`;
- the owner patterns are all minimizing colourings, modulo `G_J`.

For each fixed block, (14) is the unique Boolean Möbius transform of its
specified table `w_B`. No uniqueness is claimed for an arbitrary
hypergraph presentation of the resulting **global** colouring function. For
example, adding the same unary constant to every colour at one token and
subtracting it from every colour at another token changes the block tables
but leaves every global colouring energy unchanged. The canonical gauge here
is fixed by retaining the geometrically supplied tables (4) block by block.

For a pair of tokens,

```text
mu_B({s,t})
 =w_B({s,t})-w_B({s})-w_B({t}).                    (18)
```

It is the exact additional cost, positive or negative, of housing those two
target primes in the same endpoint block after the unary costs have been
removed. Higher `mu_B(T)` are the corresponding irreducible higher
couplings.

If

```text
mu_B(T)=0       for every B and every |T|>=2,        (19)
```

then the canonical energy (17) separates token by token:

```text
min_f E_pi(f)
 =sum_(s in S)min_(B in pi)w_B({s}),                (20)
```

and every optimum is obtained by independently choosing a tied unary owner
for each token. Thus (19) is an exact sufficient decoupling condition. A
nonzero `mu_B(T)` with `|T|>=2` is exactly the failure of the **individual
block table** `w_B` to be unary-additive. It is a coupling in the canonical
blockwise gauge, but cancellation or reparameterization across colours may
still make the global colouring function unary. Dominance may also prevent a
genuine block coupling from changing the optimum for a particular packet.

This is the promised Möbius owner hypergraph. It is not a tournament:
weights can have either sign, ties are genuine, the observable depends on
the external target factors, and hyperedges can have arbitrary arity.
Orienting pairwise comparisons discards the magnitude, colour, and every
cohabitation term of arity at least three.

## 4. Deleting a prime owner identifies an entire partition face

Now let `J` be a nontrivial prime knot. Suppose `D subset I` is a block
whose packet sum is exactly the target:

```text
K_D=J.                                              (21)
```

Put

```text
I^o=I\D,
x^o=(K_i)_(i in I^o),
R=K_(I^o).                                         (22)
```

If `I^o` is empty, use the conventions `R=U`, an empty partition, and empty
sums equal to zero. For a partition `eta` of `I^o`, define the fixed-owner
embedding

```text
i_D(eta)=eta union {D}.                             (23)
```

> **Prime-owner deletion law.** For every `eta`,
>
> ```text
> Lambda_x(i_D(eta);J)
>  =sum_(B in eta)u(K_B)
>  =Lambda_(x^o)(eta;U).                            (24)
> ```
>
> Consequently, for every `eta<=theta`,
>
> ```text
> c_x^J(i_D(eta),i_D(theta))
>  =c_x^U(i_D(eta),i_D(theta))
>  =c_(x^o)^U(eta,theta).                           (25)
> ```

### Proof

THM-2336 shows that a block equal to the prime target owns itself. Its
endpoint assigns `J` to `D` and `U` to all complementary blocks, so

```text
Lambda_x(i_D(eta);J)=sum_(B in eta)u(K_B).          (26)
```

The right side is exactly the root lift cost of the complementary packet,
proving (24). Subtract (24) at `eta` and `theta` to obtain the first
equality in (25). At the root,

```text
Lambda_x(i_D(eta);U)
 =u(J)+sum_(B in eta)u(K_B),                        (27)
```

so the constant `u(J)` cancels under coarsening and gives the other
equalities. QED.

The proof does not require the owner to be unique. Another block may attain
the same prime-target score; (24) chooses one canonical minimizing endpoint
without discarding the tied ones in `Opt`.

The obstruction itself splits into an internal complementary defect and a
two-body interaction with its total. Put

```text
delta_eta
 =sum_(B in eta)u(K_B)-u(R)
 =Omega_(x^o)(eta;U).                              (28)
```

Recall THM-2176's directional split

```text
sigma(R,J)
 =u(R)+u(J)-u(R#J),

C_J(R)
 =u(R)-d_G(R#J,J),

B_J(R)
 =d_G(R#J,J)+u(J)-u(R#J),

sigma(R,J)=C_J(R)+B_J(R).                          (29)
```

Substitution of (24), (27), and (28) gives

```text
Omega_x(i_D(eta);J)
 =delta_eta+C_J(R),                                (30)

Omega_x(i_D(eta);U)
 =delta_eta+sigma(R,J),                            (31)

Omega_x(i_D(eta);U)-Omega_x(i_D(eta);J)
 =B_J(R).                                          (32)
```

The right side of (32) is independent of `eta`. Thus the entire fixed-owner
face is a translated copy of the complementary root face: all internal
coarsening variation is `delta_eta`; translation contributes the constant
`C_J(R)` at target `J`; bypass is exactly the constant gap from the root
slice.

## 5. The first owner collision

The first cover which merges `D` with another block is not in the face
(23), but it too has an exact formula. Let `A!=D` be a block of `pi`, put

```text
K=K_A,
Rcal=pi\{D,A},
```

and let `rho` be the cover obtained by replacing `D,A` with

```text
M=D union A,             K_M=J#K.                   (33)
```

For every untouched block `C in Rcal`, use THM-2176's descent excess

```text
e_C(J)
 =e(K_C;J)
 =d_G(K_C,J)+u(J)-u(K_C)>=0.                       (34)
```

Set

```text
q=min({B_J(K)} union {e_C(J):C in Rcal}).           (35)
```

The convention in (35) gives `q=B_J(K)` when `Rcal` is empty.

> **First-collision law.**
>
> ```text
> c_x^J(pi,rho)
>  =sigma(K,J)-q
>  =C_J(K)+[B_J(K)-q],                              (36)
>
> c_x^U(pi,rho)=sigma(K,J),
>
> c_x^U(pi,rho)-c_x^J(pi,rho)=q.                   (37)
> ```
>
> In particular,
>
> ```text
> C_J(K)<=c_x^J(pi,rho)<=sigma(K,J).                (38)
> ```

The complete post-collision owner set is

```text
M belongs to Own_rho(J)       iff B_J(K)=q,

C belongs to Own_rho(J)       iff e_C(J)=q.         (39)
```

Here `Own_rho(J)` is THM-2336's tied prime-target owner set.

### Proof

At `pi`, the self-owner law gives

```text
Lambda_x(pi;J)
 =u(K)+sum_(C in Rcal)u(K_C).                       (40)
```

At `rho`, compute the prime-target scores. For the merged block,

```text
s_M(J)
 =d_G(J#K,J)-u(J#K)
 =-u(J)+B_J(K).                                    (41)
```

For an untouched block,

```text
s_C(J)
 =d_G(K_C,J)-u(K_C)
 =-u(J)+e_C(J).                                    (42)
```

Equations (35), (41), and (42) prove (39). THM-2336's owner formula now
gives

```text
Lambda_x(rho;J)
 =u(J#K)+sum_(C in Rcal)u(K_C)-u(J)+q

 =u(K)+sum_(C in Rcal)u(K_C)-sigma(K,J)+q.          (43)
```

Subtracting (43) from (40) proves (36). At the root, this cover simply
merges `J` and `K`, so THM-2330 gives

```text
c_x^U(pi,rho)
 =u(J)+u(K)-u(J#K)
 =sigma(K,J),                                      (44)
```

which proves (37). Both `B_J(K)` and all `e_C(J)` are nonnegative, while
`q<=B_J(K)`; (38) follows from
`sigma(K,J)=C_J(K)+B_J(K)`. QED.

The equality and migration boundaries are now explicit.

- The merged block remains an owner exactly when
  `B_J(K)<=e_C(J)` for every untouched `C`. Then `q=B_J(K)` and the
  target cover drop is exactly `C_J(K)`. Equalities produce genuine ties.
- If some untouched block has `e_C(J)<B_J(K)`, ownership migrates to the
  blocks of least descent excess and

  ```text
  c_x^J(pi,rho)
   =C_J(K)+B_J(K)-min_(C in Rcal)e_C(J)>C_J(K).     (45)
  ```

- If another untouched block also equals `J`, its descent excess is zero.
  Then `q=0`, it remains an owner, and the target and root cover drops are
  both `sigma(K,J)`.

The collision observable is therefore not an orientation of the pair
`(D,A)`. The outcome also depends on the external competitors through the
minimum of their descent excesses.

## 6. Brittenham--Hermiller calibration

Let

```text
K=T(2,7),                 J=mirror(K).
```

THM-2176 proves, without assuming the unknown exact value of `u(K#J)`,

```text
C_J(K)=0,

B_J(K)=sigma(K,J)>=1.                               (46)
```

For the two-block packet `(K,J)`, the first-collision law has no untouched
competitor, so `q=B_J(K)`. Therefore

```text
c_(K,J)^J(0hat,1hat)=0,

c_(K,J)^U(0hat,1hat)
 =sigma(K,J)
 =B_J(K)>=1.                                       (47)
```

The target lift cost is unchanged when the target owner is merged:

```text
Lambda_(K,J)(0hat;J)
 =Lambda_(K,J)(1hat;J)
 =3.                                               (48)
```

The root lift cost does drop. Thus the canonical nonadditive pair is
calibrated exactly at the collision boundary: its saving is pure bypass,
not translation contraction. Swapping `K` and `J` gives the same statement
in the other direction. Equations (46)--(48) assert no exact value for
`u(K#J)`.

## 7. Unary data do not determine pair coupling: an exact hostile control

The following example is an abstract conical metric monoid, not a knot
example. It shows that even root packet data, every one-token target cost,
and the unpartitioned composite-target distance can agree while the
composite allocation and obstruction differ.

Let

```text
M=N^4,

p=e_1,       q=e_2,       a=e_3,       b=e_4,

g=p+q-a-b,
v=p+q-a.                                         (49)
```

Define two translation-invariant integer word metrics on `Z^4` and restrict
them to `M`.

- `d_0` has generators `+/-e_i` and `+/-g`, all of cost one.
- `d_1` has the same generators and also `+/-v`, of cost one.

Their word lengths have the exact forms

```text
ell_0(z)
 =min_(k in Z)[|k|+||z-kg||_1],

ell_1(z)
 =min_(k,l in Z)
    [|k|+|l|+||z-kg-lv||_1].                       (50)
```

Indeed, after fixing the net coefficients of the non-coordinate generators,
the cheapest remaining coordinate correction has its `l_1` cost. Conversely
the displayed words realize every candidate in (50).

Use the definitions in (1) for this additive monoid, with connected sum
replaced by addition and `u` replaced by its word length. Since `M` is free
commutative on `e_1,...,e_4`, the factorizations of `p+q` are exactly the
allocations of the two atom occurrences `p,q`.

Direct minimization gives the following complete distance table needed
below:

```text
quantity                         d_0             d_1
-------------------------------------------------------
u(a), u(b)                       1,1             1,1
u(a+b)                           2               2
u(p), u(q), u(p+q)               1,1,2           1,1,2
d(a,p), d(a,q)                   2,2             2,2
d(b,p), d(b,q)                   2,2             2,2
d(a+b,p), d(a+b,q)               2,2             2,2
d(a,p+q)                         2               1
d(b,p+q)                         2               2
d(a+b,p+q)                       1               1. (51)
```

For example, the exceptional vector is

```text
p+q-a=v=g+b,
```

so it has `d_0`-length two and `d_1`-length one. The last vector is
`p+q-a-b=g`, of length one in both metrics. Formula (50) supplies the
matching lower bounds and all other entries in (51). More explicitly, every
listed value is at most two. Hence `|k|>=3` for `ell_0`, or
`|k|+|l|>=3` for `ell_1`, cannot improve it; the remaining finite
coefficient cases in (50) give the table.

Take the two-block packet `x=(a,b)`, its discrete partition, and the
two-token target `J=p+q`. The two metrics have identical:

```text
root subset data on (a,b);
unary scores w_a({p})=w_a({q})=w_b({p})=w_b({q})=1;
prime-target values Lambda(0hat;p)=Lambda(0hat;q)=3;
prime-target values Omega(0hat;p)=Omega(0hat;q)=1;
the tied unary owner set {a,b} at both p and q;
unpartitioned distance d(a+b,p+q)=1.                (52)
```

But their pair Möbius coefficients at the first block differ:

```text
mu_a^0({p,q})=-1,
mu_a^1({p,q})=-2.                                  (53)
```

At the second block the coefficient is `-1` in both metrics. Enumerating the
four token allocations gives

```text
                         both at a   split   both at b
-------------------------------------------------------
d_0 lift cost                 3         4         3
d_1 lift cost                 2         4         3. (54)
```

Consequently

```text
Lambda_x^0(0hat;p+q)=3,       Omega_x^0(0hat;p+q)=2,

Lambda_x^1(0hat;p+q)=2,       Omega_x^1(0hat;p+q)=1. (55)
```

For `d_0`, the two co-located endpoints tie; for `d_1`, the endpoint placing
both tokens at `a` is the unique optimum. Thus the missing coordinate is
exactly the pair hyperedge weight. Unary prime-target owner diagrams, root
defects, and the total unpartitioned target distance do not reconstruct it,
even under the formal metric-monoid and unique-atom-factorization laws.

## 8. Knot, binary-relation, and tournament interpretation

The structural ladder is now exact:

```text
prime target:
  one token -> target-dependent tied weak order on blocks;

composite target with unary Möbius support:
  independent tied owner choice for each prime occurrence;

composite target with pair support:
  coloured monochromatic-edge allocation energy;

general composite target:
  signed coloured hypergraph energy of all arities;

prime block equal to the target:
  fixed-owner face = complementary root partition face;

first merge into that block:
  translation + bypass + competitor descent-excess collision law. (56)
```

The source-to-shadow loss ledger is:

| source | shadow | preserved | destroyed information | exact repair |
|---|---|---|---|---|
| target factorization fibre | one winning block | a prime target's cheapest owner | split ownership, ties, endpoint multiplicity | `Opt_pi(J)/G_J` |
| subset-score tables `w_B` | unary prime-target scores | one-token costs | same-block target coupling | Möbius coefficients `mu_B(T)` |
| allocation hypergraph | pair graph or tournament | selected binary comparisons | signs and magnitudes, target colours, ties, all higher hyperedges | full energy (17) |
| fixed prime-owner face | a single two-knot defect | terminal translation/bypass split | internal complementary coarsening | deletion law (24)--(32) |
| owner collision | orientation of merged versus untouched block | a chosen winner | winning margin and dependence on all competitors | `q` and tied owner set (35),(39) |

There is no intrinsic tournament in (56). At a fixed target, score
comparison is a total preorder with meaningful ties. At a composite target,
the natural finite problem is a colouring with higher interactions. At the
first owner collision, even the binary-looking merge is controlled by the
minimum over every untouched block. Any lawful tournament reduction would
need a proved target-independent asymmetric observable, a tie rule which
preserves all minimizers, and the complete hypergraph/competitor data as a
sidecar.

## 9. Scope and failure boundary

1. Formula (8) is an exact finite reformulation of THM-2330 using Schubert
   prime decomposition. Its distance table can contain unknown Gordian
   distances; the dynamic program does not make them known.
2. Repeated isomorphic prime occurrences are labelled only to run the
   allocation. Actual endpoints are the quotient by (10).
3. The Möbius coefficients are signed. No positivity, submodularity,
   pairwise completeness, or bounded interaction arity is asserted.
4. The deletion law requires a block exactly equal to a nontrivial prime
   target and applies only while that block remains intact.
5. The first-collision formula applies to a cover merging that owner with
   one current block. Later coarsenings are governed by THM-2330's general
   cocycle and allocation law.
6. Owner ties are retained in (11) and (39); no generic uniqueness is
   assumed.
7. The Brittenham--Hermiller calibration uses only the already audited
   signature equalities and their cited shortcut. It does not determine
   `u(T(2,7)#mirror(T(2,7)))`.
8. The word-metric pair is a hostile axiomatic control, not an asserted
   Gordian realization.
9. None of these objects classifies knots or supplies a positive example of
   translation catalysis. They refine exactly where connected-sum
   nonadditivity lives: target allocation, internal target-factor coupling,
   owner-preserving partition loss, and first-owner collision.

## 10. Exact companion

The optimization-safe companion independently checks:

```text
the fifteen word lengths used in the N^4 hostile table;
all four allocations in each metric;
the normalized subset scores and every Boolean Möbius inversion;
the allocation DP and hypergraph energy on every colouring;
375 integer instances of the fixed-owner deletion/cocycle identities;
3900 integer instances of the collision split, bounds, and owner tests.
```

Reproduce with

```bash
python3 04-computation/prime_owner_allocation_hypergraph_thm2339.py
python3 -O 04-computation/prime_owner_allocation_hypergraph_thm2339.py
```

Both commands must equal

```text
05-knowledge/results/prime_owner_allocation_hypergraph_thm2339.out
```

byte for byte. The script uses explicit runtime checks rather than `assert`,
so optimized execution retains every verification. It verifies the finite
abstract control and the displayed algebraic identities; it is not a
Gordian-distance computation.
