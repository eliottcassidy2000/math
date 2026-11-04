# Prime-seven quotient separation except one even carousel branch

## 1. Statement, physical names, and exact remaining branch

Let Q be a prime tournament on seven fixed quotient names. Let A_i be arbitrary nonempty tournaments on fixed disjoint physical child sets, and compare T=Q[A_i] with U=Q^op[B_i], where each B_i uses precisely the same physical set as A_i and has the same complete ordered path-cover profile. Every other name and edge in a fixed ambient tournament is unchanged.

This packet proves that a child-uniform actual two-port query separates T and U in every fixed ambient, except for the following single remaining class: Q is isomorphic to the carousel

    (14,28,56,112,97,67,7)

and all seven children have even order at least four. That branch is **not** claimed here. A separate modulo-sixteen certificate can supply it. In all other cases this packet supplies one query with complete augmented-profile contrast either (2,0,0,...) modulo four or (4,0,0,...) modulo eight.

The actual ports have the fixed edge x->y. Bits a_i,c_i mean that every physical vertex of child i points to x,y, respectively. Both ports inherit the core's original incidence to every ambient exterior vertex. Their augmented core is thus a common module. Source and target are never independently relabeled: an isomorphism used to choose a quotient representative transports every child set and both query masks together.

The result is finite in quotient order and unbounded in child size. It is not unrestricted named-response injectivity. The finite response banks below are complete named polynomial banks for the explicitly specified small physical models, not full response banks for arbitrary large children.

## 2. Reviewed inherited facts

The proof uses these separately reviewed Session188 results at immutable commit `8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7`, under `05-knowledge/results/session188_structural_bridges_20261104/`:

- `parity_profile_transfer/PROOF.md`: a nonzero quotient endpoint-parity bit gives a physical gap two modulo four for arbitrary nonempty children, and the gap survives every fixed ambient.
- `uniform_block_starts/PROOF.md`, Section1: for every union I of whole children, the exact identity is

      Delta(e_I,0)-2 Delta(0,e_I)=2[S_I(T)-S_I(U)].

- `coarse_mod8_affine/PROOF.md` and `five_profile_types/PROOF.md`: the five normalized complete cover types and the exact small-run formula, including the fixed-point-free return-segment involution.
- `mod8_transfer/PROOF.md` and `quotient_cycle_parity/PROOF.md`: the complete higher-cover residues and the collapse of a child-uniform Hamilton-cycle parity to the bare quotient extension.

Here q_h counts ordered spanning covers by h nonempty directed paths, q_1=H is odd, and q_h is divisible by h!. Every normalization by H is taken modulo eight; its inverse exists. Corresponding profiles are equal, so the odd product L=product_i H(A_i) is common to both sides.

## 3. A singleton-sensitive endpoint-anchor theorem

This section applies to **any quotient order**, not just seven. Let V be exactly the set of non-singleton children. Write S_i,E_i for the bare quotient Hamilton-start and Hamilton-end counts and p_i=S_i+E_i modulo two. Let R_ij(Q) count directed run words which start at i, visit every quotient name once, and visit j one additional time. The two copies of j are indistinguishable.

Define

    t_i(V)=S_i-E_i+2 sum_(j in V)[R_ij(Q)-R_ij(Q^op)] modulo4.

Let Q_i be Q extended by a single vertex y, with j->y exactly when j=i. Let K_ij(Q_i) count its run words which visit each symbol once and child j twice. Define

    d_i(V)=H(Q_i)-H((Q^op)_i)
            +2 sum_(j in V)[K_ij(Q_i)-K_ij((Q^op)_i)] modulo4.

Then, for the actual arbitrary child profiles,

    Delta(0,e_i)=d_i(V) modulo4.

If p_i=0, then also

    Delta(e_i,0)=2t_i(V)+2d_i(V) modulo8.                 (1)

Thus these anchor residues depend only on which children are singletons, not on the size or orientation of any other child.

**Proof.** Modulo four, normalized profiles have q_1/H=1, q_2/H=2 for a non-singleton and zero for a singleton, q_3/H even, and q_h/H=0 for h>=4. Two repeated children vanish modulo four. A three-visit child has two nonempty disjoint return segments; swapping them is fixed-point-free and preserves the endpoints. Its coefficient is even, so its even q_3 contribution vanishes. This holds both for full run words with a port and for start-constrained old run words.

Consequently the aggregate start skew is L t_i(V) modulo four. For query (0,e_i), x is a universal source, so it must start the Hamilton path. Delete x; the remaining count is exactly the one-port extension above. The same run truncation gives its contrast as L d_i(V) modulo four. The contrast is even because both physical Hamilton counts are odd; equivalently the normalized contrast is even. Multiplying an even residue modulo four by odd L fixes it, proving the first assertion. If p_i=0, t_i(V) is even, so 2L t_i(V)=2t_i(V) modulo eight. The exact uniform-start identity now gives (1).

If d_i(V)=2, query (0,e_i) has complete augmented-profile contrast (2,0,...) modulo four and survives every ambient. If p_i=0 and the right side of (1) is4, query (e_i,0) has complete contrast (4,0,...) modulo eight. For even augmented order this follows from the inherited higher-cover formulas. For odd order, every Hamilton cycle contains x->y because no old vertex enters y. Deleting x,y identifies it with an old Hamilton path ending in child i. Its converse parity difference is p_i=0. Therefore the odd-order cycle correction vanishes as well. Every h>=4 coordinate is separately zero modulo eight.

The ambient substitution is sum_h c_h q_h, where c_1 is the odd Hamilton count of the contracted fixed ambient. Hence each stated pure first-coordinate gap survives. No nonzero scalar gap is lifted without this complete-profile argument.

The same return truncation proves a useful broader fact: **for every fixed V, every complete child-uniform converse contrast modulo four equals the contrast obtained by replacing each child in V by a transitive pair and every other child by a singleton**. This is a congruence of the entire actual response, including all cubic and quartic channels. It follows exactly as above; the common odd factor fixes the even contrast modulo four.

## 4. Complete seven-name census and singleton-position reduction

The independent census enumerates all 2^21 labeled tournaments and tests every nontrivial proper subset for the module condition. Exactly 973,680 labeled tournaments are prime. They are exactly the disjoint union of the relabeling orbits of the 197 explicitly retained representatives in `quotients.json`. The owner array assigns every labeled prime to its representative and every composite tournament to zero. Every representative's 14 integer start/end counts is independently recomputed by all 7! permutations.

Of the 197 prime representatives, 159 have p nonzero and are covered for every child assignment by the inherited endpoint theorem. The remaining38 representatives have p=0. For each, all128 possible V are evaluated using Section3. Of these4,864 patterns,4,854 have a successful anchor. The ten remaining patterns are retained completely in `anchor_exploration.json`.

Six of the ten contain at least one singleton. The following physical block queries have gap two modulo four, by the broader transfer at the end of Section3. Here V is the non-singleton bitmask, with quotient names0,...,6.

    quotient rows                         V     (a,c)
    (6,12,24,49,99,71,15)                 65     (1,4)
    (6,12,48,21,99,75,15)                101     (1,2)
    (6,20,24,35,105,71,15)                83     (1,4)
    (8,5,17,70,43,15,55)                  54     (1,6)
    (14,28,56,112,97,67,7)                 0     (1,2)
    (14,52,88,98,41,69,19)                 0     (1,2)

The complete841-coordinate block-value contrast is retained for every row, not merely the listed witness. Consequently **every singleton-containing prime-seven substitution is already separated in every ambient**. The remaining four patterns have V=127, on the following four quotients:

    Q1=(8,5,17,70,43,15,55),
    Q2=(14,28,56,112,97,67,7),
    Q3=(14,52,88,98,41,69,19),
    Q4=(66,76,25,49,35,7,60).

These exact rows, rather than unverified class names, determine all subsequent claims.

## 5. The full nontrivial profile domain modulo eight

For a nontrivial child, its normalized profile is one of

    P=(1,2,0), E=(1,6,4), O=(1,6,6), C=(1,2,2),

listed as coordinates q_1/H,q_2/H,q_3/H, with higher coordinates zero modulo eight. Put b=1 for E,O and c=1 for O,C. Here c is the child's order parity. An actual C3 has H=3 and is not replaced by a literal H=1 tournament.

For each block query, let K0 be the all-single-run word count; K_i the count with child i twice; L_i the count with child i three times; and K_ij the count with children i,j twice. The inherited complete formula gives

    Phi_Q=K0+sum_i(2+4b_i)K_i+2sum_i c_i L_i
                   +4sum_(i<j)K_ij modulo8.             (2)

All physical endpoints, both possible port orders, and every cubic/quartic term are included. The exact extraction from physical transitive small models is

    K0=F0,
    K_i=(F_i2-F0)/2,
    L_i=(F_i3-3F_i2+2F0)/6,
    K_ij=(F_ij22-F_i2-F_j2+F0)/4.

The L_i numerator is divisible by12, reflecting the even return coefficient. Small models have at most nine old physical names. Their source and target interiors retain the same physical order and names.

For each Q1,...,Q4, the converse contrast from(2) is an affine function of the14 bits(b_0,...,b_6,c_0,...,c_6), with values divisible by four. Its14 derivative columns are independent over F_2 in the complete841-coordinate block polynomial bank. Solving the full affine system gives:

    Q1: rank14, no zero contrast type;
    Q2: rank14, unique zero contrast type EEEEEEE;
    Q3: rank14, no zero contrast type;
    Q4: rank14, no zero contrast type.

These are complete finite14-bit domain statements, not random tests. The Boolean query basis consists of all masks a,c of individual weight at most two. Because each port has at most two incident edges in a Hamilton path, its full physical response has bidegree at most(2,2). Thus this basis determines every coefficient and every block query exactly.

For even total child order, every nonzero contrast in this branch is four modulo eight and automatically has complete profile gap(4,0,...) modulo eight. It therefore survives every fixed ambient.

## 6. Odd total order and safe actual queries

Local nonvanishing alone would not justify an odd-total ambient conclusion. We independently count Hamilton cycles, rooted at the unique port x, in the two bare quotient extensions at every basis query. Define a query to be safe when those two integer cycle counts have equal parity. The inherited cycle-collapse theorem shows that this exact test also applies after substituting arbitrary actual child tournaments.

Restrict the entire affine response system to the safe **actual query values**, rather than merely selecting polynomial coefficient coordinates. Its complete results are:

    Q1: 445 safe queries, rank14, no zero type;
    Q2: 589 safe queries, rank14, only EEEEEEE zero;
    Q3: 421 safe queries, rank13, only PPPPPPP and EEEEEEE zero;
    Q4: 517 safe queries, rank14, no zero type.

Both Q3 zero types have c=0 and even total order. The Q2 zero type is also even. Thus **every odd-total actual type assignment has a safe query with Hamilton contrast four modulo eight**. Its higher-cover contrast vanishes: for odd augmented order,

    Delta q2=Delta q3=2 Delta H-4 Delta Z modulo8,

and here Delta H=4 and Delta Z=0 modulo two. All later cover coordinates vanish by factorial divisibility. This proves the ambient assertion for every odd-total branch without assuming that an unsafe local gap survives.

Sections4–6 leave exactly the even carousel EEEEEEE branch excluded in Section1.

## 7. Data retention, independent checks, and reproduction

`named_banks.json` and `extra_named_banks.json` retain every integer value at the full named bidegree-(2,2) basis and every integer coefficient, including every zero, cubic, and quartic, for every specified physical small model on both orientations. Each record includes its complete physical outgoing rows, physical child sets, quotient rows, child sizes, port convention, masks, and full block banks. The block banks are obtained by idempotently aggregating the named monomials, not by forgetting physical names before extraction.

`affine_banks.json` retains both orientations' complete36-channel run coefficient banks (841 coordinates per channel), not only ranks or hashes. An independent mixed-radix run-word engine verifies all242,208 integer run-channel values against the physical model extraction. The physical engine counts distinct physical vertices by subset DP; the run-word engine counts indistinguishable repeated child symbols by a different state recurrence.

`cycle_banks.json` retains both exact rooted cycle banks and their full named coefficient arrays. `ambient_analysis.json` retains every cycle-parity contrast and every safe query index. `anchor_exploration.json` retains all3,990 exact anchor primitive coordinates and all4,864 pattern banks. An independent mixed-radix program recomputes all3,990 primitives. `residual_mod4.json` retains every coordinate of the six residual singleton banks. `quotients.json`, `census_values.txt`, and `prime_owners.u16le` retain the complete finite census, exact endpoint data, and disjoint orbit assignment.

The standard-library Python and installed C++17 compiler suffice. No software installation, startup helper, or unrelated repository scan is required. All verification checks use explicit exceptions and remain active with optimized Python. Exact physical orders in the counted small models are at most11 including ports, so integer Hamilton counts are below11!; every retained computation fits safely in64 bits.

Run `python3 verify.py > normal.json`, `python3 -O verify.py > optimized.json`, and compare the two receipts. The source pipeline regenerates every mathematical data artifact and compares it byte-for-byte with the retained packet. Independent audit is still required before integration with the separate modulo-sixteen branch.
