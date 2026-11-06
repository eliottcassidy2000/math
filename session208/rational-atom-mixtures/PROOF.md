# Rational atom mixtures for inverse Hamiltonian-path order

Status: PROVED finite-dictionary theorem; FINITE-EXACT dictionary through strong order eight; global construction bound only. No global minimum-order formula or new proof of the attached paper's main theorem is claimed.

## 1. The precise paper bridge

Hajiaghayi, *Logarithmic Embeddings into Probabilistic Spanning Trees for Distances and Capacities*, arXiv:2611.05189v1, Corollary 8.2 (pp.31–32), compresses rational mixtures by eliminating linear dependencies while preserving constraints. The transferable mechanism is rational nonnegative combinations, linear prices, and support elimination. Tree-distance geometry itself does not enter this result. The paper's probability normalization is replaced here by prime-exponent equalities; the feasible set is a cone slice, not a probability simplex. Everything below has its own elementary proof and does not assume the paper's logarithmic embedding claim.

The inherited Session207 inverse recurrence states that Hamiltonian-path count multiplies under dominance joins of strong components and orders add. For any h>1 with a witness, m(h) is the minimum sum of strong-component orders over multiplicative factorizations of h. Session207 also proves existence of lim_k m(h^k)/k by subadditivity.

## 2. Rational support compression becomes an integer construction

Fix distinct primes p_1,...,p_r. A finite dictionary D consists of actual nontrivial strong tournaments C_j with counts H(C_j)=product_i p_i^(a_ij)>1 and costs c_j=|C_j|. Write A=(a_ij); all columns are nonzero nonnegative integer vectors. Duplicate exponent columns may be replaced by their least cost. Let b be a nonzero nonnegative integer vector and h=product_i p_i^(b_i). Assume Ax=b has a nonnegative integer solution (the dictionary realizes h). Define

 M_D(k)=min {c·z: Az=kb, z nonnegative integer},
 λ_D(b)=min {c·x: Ax=b, x nonnegative real}.

M_D(k) is the minimum order among tournaments whose nontrivial SCCs are chosen from D and whose count is h^k. Singleton SCCs cannot help minimize order.

Theorem. The rational LP has an optimal rational solution supported on at most rank(A)≤r atom types. There are integers D≥1 and C≥0 such that, for every k≥1,

 k λ_D(b) ≤ M_D(k) ≤ k λ_D(b)+C.

For every multiple of D equality holds. In particular lim M_D(k)/k=λ_D(b), and the global inverse satisfies lim m(h^k)/k≤λ_D(b).

Proof. The feasible set is nonempty and bounded: for each column j choose i with a_ij≥1, so x_j≤b_i/a_ij≤max_i b_i. A linear objective attains its minimum. Among optimal points choose one with smallest positive support. If its supported columns are dependent, choose a nonzero rational dependence v with Av=0. Because the columns are nonzero and nonnegative, v has both signs. Small perturbations x±tv remain nonnegative. Optimality forces c·v=0, or one sign would decrease cost. Move until a coordinate vanishes, contradicting minimal support. The supported columns are independent; their solution to Ax=b is rational. Thus the support has size at most rank(A).

Clear denominators with D so w=Dx is a nonnegative integer vector. Then Aw=Db and c·w=Dλ. Fix an integer anchor z0 with Az0=b and set C0=c·z0. If k=qD+s with 0≤s<D, concatenate q copies of the w-block and s copies of the anchor. This costs qDλ+sC0 = kλ+s(C0−λ), giving C=ceil((D−1)(C0−λ)). The LP relaxation gives the lower bound. At s=0 they coincide. QED.

A price vector y in R^r with A^T y≤c certifies λ≥b·y. A feasible x with c·x=b·y certifies equality, without a numerical optimizer. Prices need not be nonnegative in the general statement. A dictionary certificate is not a global lower bound: that would require the inequality for EVERY strong atom whose primes lie in the chosen set.

## 3. A checked 3–5 mixture

Use all strong tournaments of order at most eight whose counts have no prime divisors other than 3 and 5, retaining the smallest order for each count. The complete Session207 census gives these fifteen entries, written (H; exponent of 3, exponent of 5; minimum strong order):

 (3;1,0;3), (5;0,1;4), (9;2,0;5), (15;1,1;5),
 (25;0,2;6), (27;3,0;6), (45;2,1;6), (75;1,2;7),
 (81;4,0;7), (125;0,3;7), (135;3,1;7), (225;2,2;8),
 (243;5,0;8), (375;1,3;8), (405;4,1;8).

For b=(1,1), prices y=(16/11,24/11) give integer slacks 11c−16a−24b respectively

 17,20,23,15,18,18,10,13,13,5,5,8,8,0,0.

Thus every dictionary realization of 15^k has order at least 40k/11. On the other hand, three 375-atoms and two 405-atoms give

 375^3 · 405^2 = 3^11 · 5^11 = 15^11,

on forty vertices. The rational optimal weights are 3/11 and 2/11, and λ_D(1,1)=40/11. Equality holds for all k divisible by 11. Moreover any dictionary LP optimizer uses only the two zero-slack types, so these weights are uniquely determined. A one-type optimizer is impossible. This is support size two, not two SCCs: the denominator-cleared block has five SCCs.

For k=11q+s, 0≤s≤10, a simple uniform global construction bound is

 m(15^k) ≤ 40q+5s,

using s copies of an order-five H=15 atom. The independently computed smaller residue costs are

 d_s = [0,5,8,13,15,20,23,27,30,35,38],

so the stronger bound m(15^k)≤40q+d_s also holds. Consequently

 lim_k m(15^k)/k ≤ 40/11 < 4.

The comparison 4 comes from repeatedly using an order-eight H=225 atom. None of these global upper bounds assert that a large strong atom cannot do better; the LP lower bound applies only to the fixed dictionary. The forty-vertex construction is not asserted to be a globally minimum witness for 15^11.

## 4. The exact all-dictionary envelope

Fix an attainable h>1, with its prime set. Over all finite strong-atom dictionaries realizing h,

 lim_k m(h^k)/k = inf_D λ_D(b).

Proof. The previous theorem gives the ≤ direction. Conversely, for any k choose a minimum-order witness for h^k, collect its finitely many nontrivial strong SCCs, and add the SCCs of any witness for h. This dictionary realizes h. Dividing the multiplicity vector of the h^k witness by k gives a feasible rational point of cost m(h^k)/k. Hence inf_D λ_D≤inf_k m(h^k)/k, which equals the limit by Session207's subadditivity argument. QED.

This is an exact variational characterization, not an effective finite certificate of the global infimum. No attainment, positivity, convergence rate, or global duality theorem for the infinite dictionary is asserted.

## 5. Reproduction and scope

`verify.py` reads the attached copy of the Session207 census, independently counts paths of the 375 and 405 witnesses by all 8! permutations, verifies their strong connectivity, reconstructs the complete smooth dictionary, checks every exact rational price inequality, checks the exponent identity, and computes the residue construction costs by integer DP. Exhaustiveness of the dictionary inherits the audited Session207 census; the new verifier does not rerun that full census. It also verifies that each retained witness has the advertised count and is strong using subset DP and reachability.
