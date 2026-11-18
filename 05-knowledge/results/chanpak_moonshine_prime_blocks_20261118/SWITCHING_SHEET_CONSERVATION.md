# Switching sheets: factorial conservation inside the tournament metagraph

**All-order:** an exact action, orbit-weight identities, and balanced complement flow. **Finite-exact:** n<=5 sheets and color-resolved flows. **Not proved:** any new universal spectrum or A000568 recurrence. Elementary tiling fibre and blue/black laws are inherited from THM-643 and THM-796; the point here is their coupling to switching-isomorphism sheets.

## 1. The transversal and induced action

On labelled vertices 1,...,n let X_n consist of tournaments containing the fixed directed path 1->2->...->n. Its free-coordinate dimension is d=(n-1)(n-2)/2 and |X_n|=2^d. A cut switch multiplies each skew edge sign b_ij by s_i s_j, with s_i in {+1,-1}; simultaneous negation has no effect.

Every labelled switching class has exactly one representative in X_n. Indeed the path requirements recursively determine s_(i+1)=s_i b_(i,i+1), uniquely up to a common sign. Call this representative nu(T).

The symmetric group acts on switching classes by relabelling, because relabelling carries a cut into a cut. Therefore it acts on X_n by

    g * x = nu(g x).

This is genuinely an action: nu(g nu(hx))=nu(ghx), since the intermediate normalization differs from hx by a cut and g carries that cut to a cut. Identity acts identically. Its orbits are the **switching-isomorphism sheets**. They are not the ordinary tournament isomorphism fibres, nor generally connected components of the complement metagraph.

Let G_S be the permutation stabilizer of one labelled switching class belonging to sheet S. Its order is independent of the chosen representative. Orbit-stabilizer gives

    |S| = n! / |G_S|.

G_S consists of permutations preserving the class up to a cut. It is not the ordinary automorphism group of an oriented tournament and can have even order.

## 2. Two conservation laws in every sheet

For an ordinary tournament isomorphism class D, let a_D be the number of vertex permutations preserving every directed edge, and h_D its Hamiltonian path count. Its fixed-path tiling fibre has size h_D/a_D. Proof: write the tournament along an ordered Hamiltonian path; two paths give the same labelled fixed-path tournament exactly when an automorphism maps one ordered path to the other. That action is free, so a_D divides h_D.

Fix a labelled switching class in S. The group G_S acts on its 2^(n-1) cut outputs. Its orbits are exactly the ordinary isomorphism classes D in S: any isomorphism between two cut outputs necessarily preserves the underlying switching class. The stabilizer of a cut output is precisely its ordinary automorphism group. Thus D occurs among cuts exactly |G_S|/a_D times. Consequently

    sum_(D in S) 1/a_D = 2^(n-1)/|G_S|,
    sum_(D in S) h_D/a_D = n!/|G_S|.

The second identity also follows from the tiling orbit size. Equivalently, with W_S(z)=sum_D z^(h_D)/a_D,

    |G_S| W_S(1)=2^(n-1),
    |G_S| W'_S(1)=n!.

The polynomial |G_S| W_S has nonnegative integer coefficients: it is exactly the Hamiltonian-count histogram of the cut class. These equations refine the global tiling identity sum_D h_D/a_D=2^d, but do not determine the support of W_S.

An elementary direct proof of the second cut identity is that every ordered permutation is a Hamiltonian path in exactly one cut output: the consecutive edges form a spanning tree and uniquely solve the sign equations.

## 3. Exact complement flow, with merged-node cautions

Let kappa complement every free coordinate while leaving the fixed path intact. For n>=3 it is fixed-point-free. This is the full complement-partner line system, not the hypercube graph of single-coordinate flips. Define

    M_ST = #{x in S : kappa(x) in T}.

Then M_ST=M_TS, sum_T M_ST=|S|, and M_SS is even. Thus every sheet has a conserved supply of line endpoints and a symmetric transfer matrix. Kappa need not preserve sheets and is not asserted to arise from the S_n action above.

Let rho take the tournament converse and then reverse all vertex labels. It preserves X_n and on free bits is reflection of their positions. Rho commutes with kappa. The usual blue lines are kappa pairs whose two tilings are rho-fixed; other pairs are black. Replacing M by its rho-fixed subset gives a symmetric blue flow B, with even diagonal. Its row sum is the number of rho-fixed tilings in that sheet. A sheet may be exchanged with a different sheet by rho; no assumption that each sheet is self-converse is needed for these identities.

Ordinary converse merging combines the D and D* tiling fibres. A chiral merged node has mass 2h_D/a_D; a self-converse node has mass h_D/a_D, not twice that. Each blue/black self-line counts twice toward degree. A node is pure blue exactly when all tilings in its merged fibre are rho-fixed. A self-converse node need not be pure blue. No new pure-blue classification is asserted.

Globally rho fixes 2^((d+f)/2) tilings, with f=floor((n-1)/2) reflection-fixed free positions. Hence for n>=3 there are 2^((d+f)/2-1) blue lines and 2^(d-1) total lines, inherited elementary laws. These power-of-two identities must not be replaced with a claimed universal mod-24 law.

## 4. The exact 24 and the missing 3,5,13

At n=5 there are two sheets, of sizes 40 and 24, with switching stabilizer orders 3 and 5. Their ordinary class data (h,a) are:

    S40: (9,1),(13,1),(5,1),(15,3),(3,3),(5,1),(3,3),(3,3).
    S24: (9,1),(11,1),(15,5),(1,1).

The repeated entries are distinct ordinary isomorphism classes. Their tiling masses sum to 40 and 24. Their cut histograms are respectively

    3:3, 5:6, 9:3, 13:3, 15:1;
    1:5, 9:5, 11:5, 15:1.

Thus the 24-sheet misses counts 3,5,13 at this order. But the complement metagraph crosses between sheets with exact directed flow, in the order (S40,S24),

    M = [[24,16],[16,8]],     B = [[4,4],[4,4]].

There are 16 cross-sheet lines, four blue. Within S40 there are 12 total lines, two blue; within S24 there are four total lines, two blue. This is a physical route between source sheets, not a theorem that the route hits an arbitrary prescribed count.

Here 24 is rigorously 5!/5, not a mysterious common residue. At n=4 the two sheet sizes are instead 6 and 2, with switching stabilizers 4 and 12. At n=3 the single sheet has size2 and stabilizer3. These initial cases rule out interpreting the 24 as a universal sheet size or modulus.

## 5. Historical connections and exact boundaries

- THM-643 (`blue-black-line-parity-structure`) and THM-796 (`three-sorted-recursive-tiling-line-node-incidence`) supply the earlier fibre and colored-line framework.
- `metagraph_20261020_closed_open_tilings_pure_blue_alphabet.md` explicitly inventories those prior laws and records a scope correction for the later pure-blue alphabet proof. We use none of that provisional classification.
- THM-4602 supplies the earlier projective/Legendre switching interpretation. The Paley pole chart in `paley_target_bridge/PROJECTIVE_SWITCHING_SELECTOR.md` physically realizes every cut but retains its sheet. Complement flow now describes a different physical operation that can leave the sheet.
- `planar_quintic_s5_petersen_resolvent_20261019.md` recovers K5/Petersen incidence and an A5 Galois reduction. Its lesson that the permutation structure and its oriented realization are different objects applies here. The switching stabilizers need not be tournament automorphism groups.

The Leech lattice, Monster, and supersingular-prime collection {4,2,5,3,7,13} can motivate looking for group actions and stabilizer quotients. The present calculation establishes no lattice embedding, moonshine relation, or structural identification with those objects. In particular the appearance of 24 and 13 in this finite table is not such an identification. The precise reusable principle is **refine a count by its stabilizer, then conserve weighted mass across physical transitions**.

## 6. Replay

`python switching_sheet_replay.py` independently enumerates all fixed-path tilings for 2<=n<=5, all label permutations, the normalized action orbits, ordinary automorphisms, switching stabilizers, Hamiltonian counts by subset DP, the two weighted identities, and full/blue complement flows. It also checks the multiplicity |G_S|/a_D directly in every sheet. Normal and optimized runs agree. JSON retains all ordinary class data and canonical encodings. This is not a census at larger orders and supplies no new universal spectrum result.
