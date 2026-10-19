# Ideal-frontier positivity, sharp contraction, and singular bridgeless examples

Date: 2026-10-19 UTC. This note proves statements for ordinary rank-layer ideal transfers. It does not assert a solution of the 1/3–2/3 conjecture, a finite invariant-cone certificate, or literature novelty.

## 1. Definitions

Let P be a finite poset whose incomparability graph has maximum degree at most D, with D>=1. Let I_k be the set of all order ideals of cardinality k. All such ideals are reachable prefixes of linear extensions.

For 0<=k<=k+L<=|P|, define the full rectangular rank-layer transfer matrix K_(k,L), with rows indexed by I in I_k and columns by J in I_(k+L), by the number of ways to append L admissible elements, starting at I and ending at J. Equivalently,

    K_(k,L)(I,J)=0                    if I is not contained in J,
    K_(k,L)(I,J)=e(P[J minus I])      if I is contained in J.

The latter equality holds because I and J are ideals. Every predecessor outside J is irrelevant, and every predecessor of a point of J is already in J. Once I is selected, the only remaining ordering constraints are those of the induced poset on J minus I.

These are the actual reachable ideal states, not an enlarged mask alphabet containing impossible states. Artificial padding states can of course produce permanent zero rows or columns; none of the conclusions below refer to those.

## 2. Strict positivity after D steps

**Lemma 1.** If I,J are ideals and |J|-|I|>=D, then I is contained in J.

**Proof.** Suppose x belongs to I minus J. Every y in J minus I is incomparable to x. Indeed, y<x would force y into I, and x<y would force x into J. But

    |J minus I|=|J|-|I|+|I minus J|>=D+1,

contradicting the incomparability degree bound on x. QED.

**Corollary 2.** Every full transfer product of L>=D consecutive ranks is entrywise strictly positive. For L=D,

    1 <= K_(k,D)(I,J) <= D!.

No ordinal indecomposability, bridge condition, reference labeling, or stationary transfer rule is needed.

The positivity delay D is best possible. In a (D+1)-element antichain, choose a singleton ideal I and its complementary D-element ideal J. Their rank difference is D-1 and the transfer entry is zero. For D>=2 this incomparability graph is connected and bridgeless.

There is also a width-two sharp-delay example: the disjoint union of chains C_2 and C_D. Its incomparability graph is K_(2,D), connected and bridgeless for D>=2. Take I to be the first element of C_2 and J to be all of C_D. Again the rank difference is D-1 and the transfer entry is zero.

## 3. Binomial diagonal scaling

The factorial entry bound is substantially weaker than the following structured bound.

**Lemma 3.** Fix a D-step kernel K=K_(k,D), and let U be the union of all its source ideals. Put q=|U|-k. Then 0<=q<=D, and there are positive row factors a_I and column factors b_J such that

    1 <= K(I,J)/(a_I b_J) <= binom(D,q).

One may take

    a_I=e(P[U minus I]),
    b_J=e(P[J minus U]).

**Proof.** A union of ideals is an ideal. By Lemma 1, every source ideal is contained in every target ideal, so U is contained in every J. Therefore U minus I has q elements, while J minus U has D-q elements.

Any extension of U minus I followed by any extension of J minus U is an extension of J minus I, because U is an ideal. This gives the lower bound a_I b_J.

Conversely, any extension of J minus I determines its restrictions to these two parts and a shuffle of lengths q and D-q. There are at most binom(D,q) such shuffles. This gives the upper bound. QED.

A dual improvement is available for an individual kernel. Let W be the intersection of all target ideals, and put b=k+D-|W|. The same argument, splitting at W instead of U, gives an upper scaling factor binom(D,b). Thus either splitting may be used, and the smaller of the two factors is a valid bound.

In particular, define

    C_D=binom(D,floor(D/2)),
    tau_D=(C_D-1)/(C_D+1).

Every D-step kernel is diagonally equivalent to a positive matrix whose entries lie in [1,C_D]. Here “diagonally equivalent” means multiplication on the left and right by positive diagonal matrices; it is not a claim about similarity or eigenvalues.

## 4. Sharp Hilbert contraction bound

For positive row vectors x,y of the same length, define Hilbert distance by

    d_H(x,y)=max_i log(x_i/y_i)-min_i log(x_i/y_i).

Positive diagonal coordinate scaling and normalization by positive scalars leave this distance unchanged.

**Theorem 4.** Every D-step ideal kernel satisfies

    d_H(xK,yK) <= tau_D d_H(x,y).

Its image has projective diameter at most 2 log C_D. Both constants are best possible uniformly over posets with maximum incomparability degree D.

### Proof of the upper bound

After the scaling from Lemma 3, it is enough to handle a positive matrix H with entries in [1,C], where C<=C_D. Its cross ratios satisfy

    H_ij H_rl/(H_il H_rj) <= C^2.

A direct elementary argument proves contraction by (C-1)/(C+1).

First, if probability distributions mu,nu have likelihood ratio values in [a,b], then a<=1<=b and

    TV(mu,nu) <= (b-1)(1-a)/(b-a).

This follows by bounding the convex function (z-1)_+ on [a,b] by the chord joining its endpoint values and averaging with respect to nu. If b/a<=R, maximizing the right side gives

    TV(mu,nu) <= (sqrt(R)-1)/(sqrt(R)+1).

Indeed, with b=Ra the expression is

    (R+1-Ra-1/a)/(R-1),

maximized at a=1/sqrt(R). Degenerate cases follow by continuity.

Now put h_i=log(x_i/y_i) and x_i(t)=y_i exp(t h_i), for 0<=t<=1. For output coordinate j, set f_j(t)=(x(t)H)_j. Its logarithmic derivative is the mean of h under the probability distribution

    mu_i^j(t)=x_i(t) H_ij/f_j(t).

For any two output coordinates j,l, the likelihood ratios between mu^j and mu^l have maximum/minimum ratio at most C^2. Therefore

    |d/dt log(f_j(t)/f_l(t))|
       <= ((C-1)/(C+1)) (max h-min h).

Integrating from 0 to 1 and taking the maximum over j,l gives the asserted contraction. The diameter bound follows directly from the cross-ratio bound, or by comparing coordinate ratios of arbitrary nonzero nonnegative input images. QED.

This is the standard Birkhoff–Hopf mechanism, supplied here with an elementary proof. A published reference for its exact general matrix formulation is Theorems 2–3 in Gautier, Hein and Tudisco, “The Global Convergence of the Nonlinear Power Method for Mixed-Subordinate Matrix Norms,” Journal of Scientific Computing 88 (2021): https://link.springer.com/article/10.1007/s10915-021-01524-w . The poset-specific contribution within this note is the ideal nesting and binomial diagonal-scaling argument, not the general contraction theorem.

### Proof of universal sharpness

For D>=2, let P be the disjoint union of two chains of length D, and set r=floor(D/2). Represent ideals by the two prefix lengths. Consider the D-step kernel from rank r to rank D+r. Its two extreme source ideals are

    I1=(r,0), I2=(0,r),

and two extreme target ideals are

    J1=(r,D), J2=(D,r).

The corresponding submatrix is exactly

    [[1,C_D],[C_D,1]].

The off-diagonal entries count shuffles of two chains of lengths r and D-r; the diagonal entries count a single D-element chain. Its cross ratio is C_D^2, so the image diameter is at least 2 log C_D and hence equals the universal bound.

For the contraction constant, concentrate positive input mass on these two source states, with all other source weights tending to zero. Perturb their two weights from (1,1) to (exp(epsilon),exp(-epsilon)). For the two selected output coordinates, the output log-ratio changes to first order by (C_D-1)/(C_D+1) times the input Hilbert distance. Taking the other weights to zero and epsilon to zero proves sharpness as a supremum. D=1 has C_D=1 and contraction 0.

The sharp examples already have connected, bridgeless incomparability graphs K_(D,D), so a bridge exclusion cannot improve this general one-block constant.

For D=7, the exact sharp uniform bound is

    C_7=35, tau_7=17/18.

This improves the cruder factorial-based coefficient 5039/5041.

## 5. Consequence for boundary forgetting

Let K1,...,Kt be consecutive D-step kernels, with their possibly different actual state dimensions. For any nonzero nonnegative starting vectors x,y,

    d_H(x K1...Kt, y K1...Kt)
       <= 2 log(C_D) tau_D^(t-1),       t>=1.

The first positive block places both vectors in a cone of bounded projective diameter; each later block contracts by tau_D. The statement also holds in the reverse direction using transposed kernels.

For normalized positive vectors, the likelihood-ratio lemma gives

    TV(p,q) <= tanh(d_H(p,q)/4).

Consequently the corresponding normalized boundary-state distributions differ in total variation by at most

    tanh((log C_D)/2 * tau_D^(t-1)).

Multiplying both vectors coordinatewise by any fixed positive continuation weights does not change their Hilbert distance, so this bound remains valid after conditioning on the same right-hand continuation. Any bounded [0,1]-valued observable of that continuation has expectation difference bounded by the same total-variation estimate.

This gives a rigorous approximation and stability mechanism. It does not, by itself, prove an exact finite rational cone certificate or force any comparison probability across the threshold 1/3. A witness might approach a threshold without a uniform margin, and the useful witness might change as boundary conditions vary. Those remain separate tasks.

## 6. Positive products can remain singular without an incomparability bridge

For every m>=2, define a poset Q_m on two chains

    a1<...<am, b1<...<bm,

with cross relations

    ai<bj and bi<aj exactly when j>=i+2.

Its only incomparable pairs are ai,bj with |i-j|<=1. Thus its maximum incomparability degree is 3 once m>=3. The incomparability graph is connected. Every edge lies in a K_(2,2) square on two consecutive levels (using the previous square for the final diagonal edge), so it has no bridges. Connectedness also makes the poset ordinal-indecomposable: any nontrivial ordinal sum would disconnect the incomparability graph.

An ideal is specified by the two chosen prefix lengths i,j. It is an ideal exactly when |i-j|<=2, subject to 0<=i,j<=m. At an interior even rank 2h the three states, in increasing i order, are

    (h-1,h+1), (h,h), (h+1,h-1).

At rank 2h+1 the two states are

    (h,h+1), (h+1,h).

Away from endpoints, the one-step matrices are therefore

    U=[[1,0],[1,1],[0,1]],
    V=[[1,1,0],[0,1,1]].

The full even-to-even two-step transfer is

    T=UV=[[1,1,0],[1,2,1],[0,1,1]].

It has kernel vector (1,-1,1), and eigenvalues 0,1,3. Hence every power T^q, q>=1, has rank exactly 2. In contrast,

    T^2=[[2,3,1],[3,6,3],[1,3,2]]

is entrywise strictly positive. Taking m arbitrarily large yields arbitrarily long positive, singular, full even-rank transfer products in ordinal-indecomposable 3-thin posets with no incomparability bridges.

Important scope: this disproves an implication from full-frontier singularity to a bridge decomposition. It does not prohibit compression to an effective two-dimensional space or a better choice of cut phase. Indeed the odd-to-odd transfer is VU=[[2,1],[1,2]], which is nonsingular. No claim is made that every possible phase or state compression must remain singular.

## 7. Verification scope

From the repository root run `python3 05-knowledge/experiments/verify_ideal_transfer_20261019.py`. Its published exact output is `05-knowledge/results/ideal_transfer_positivity_20261019.json`.

It checks:

- the zero entry at delay D-1 in the antichain sharp-delay construction, for D=2,...,7;
- positivity and exact maximum cross ratio C_D^2 for the two-chain sharp-contraction construction, for D=2,...,7;
- the degree, graph connectedness, and absence of every possible incomparability bridge for Q_6;
- the exact matrices T and T^2, their ranks, and positivity of both three-step and four-step kernels at a selected interior layer of Q_6.

These are small fixed controls, not a census. The all-size conclusions follow from the proofs above.

## 8. Short primary-source comparison and novelty boundary

Existence and exhaustion-independence of limiting probabilities for locally finite thin infinite posets are already established work: Graham Brightwell, “Linear extensions of infinite posets,” Discrete Mathematics 70 (1988), 113–136, DOI https://doi.org/10.1016/0012-365X(88)90087-8 . Brightwell–Felsner–Trotter, “Balancing pairs and the cross product conjecture,” Order 12 (1995), Sections 1 and 7, explicitly state and use that convergence theorem; the author-hosted primary text is https://trotter.math.gatech.edu/papers/97.pdf . The same paper identifies the infinite Fibonacci thin order as the sharp bulk obstruction. Thus existence of a limiting law and the boundary/bulk distinction are not claimed as new here.

Exact counting and probability computation on the ideal lattice are also established: De Loof, De Meyer and De Baets, “Exploiting the Lattice of Ideals Representation of a Poset,” Fundamenta Informaticae 71 (2006), 309–321, https://doi.org/10.3233/FUN-2006-712-309 . Its primary abstract explicitly describes extension, rank, and comparison counting on the ideal lattice.

The short search checked these sources, the current Aires–Kahn and Aires–Chan–Pak–Panova preprints, and targeted combinations of bounded incomparability, ideal layers, transfer matrices, and projective contraction. It did not identify the precise D-step nesting/binomial-scaling/sharp-contraction proposition above. That absence is not a novelty proof: the full 1988 paper was not accessible in this check, and a comprehensive priority search was not conducted. The publishable claim is the self-contained scoped proposition and sharp constant, with historical priority left unasserted.
