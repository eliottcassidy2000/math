# Independent audit: odd-order blue-spectrum pairing and hollow-frame cubes

Audit date: 2026-10-20. Verdict: APPROVED, subject to the definitions and scope cautions below. No fatal mathematical gap was found. This is a proof audit, not merely a finite verification.

## 1. Definitions that must remain distinct

Let G = Aut(T), and let I(T) be the set of involutory anti-automorphisms. For rho in I(T), let P_rho(T) consist of ordered directed Hamiltonian paths p such that rho(p_j)=p_{n+1-j}. Write h_rho=|P_rho|. The blue count B(T) is the number of blue tilings, equivalently the number of G-orbits of Hamiltonian paths whose positional reversal is an anti-automorphism.

The correct identities are

    B(T) = h_rho / |C_G(rho)|,
    H_anti(T) = sum_{rho in I(T)} h_rho = B(T)|G|.

THM-644 defines H_anti as the aggregate, not a fixed-rho count. THM-647's first proof of oddness of h_rho is sound independently of any ambiguous use of H_anti in its second proof. Never replace B by h_rho without proving the centralizer trivial.

Example: T=C3 -> 1 -> C3 has |G|=9, |I(T)|=3, |C_G(rho)|=3, h_rho=3, H_anti=9, and B=1.

## 2. Anti-involution existence and conjugacy

For n>=2, the automorphisms and anti-automorphisms form a signed permutation group G^pm of order 2|G| whenever T is self-converse. Tournament automorphism groups have odd order: an involutory automorphism would swap the endpoints of an arc, which is impossible; apply Cauchy's theorem.

A Sylow 2-subgroup of G^pm has order two, so an anti-involution exists. All anti-involutions are conjugate under G^pm. They are in fact conjugate under G: if g rho g^{-1}=rho' and g is anti, the automorphism g rho has the same conjugation effect on rho.

An anti-involution has at most one fixed vertex, since the arc between two fixed vertices could not reverse. At odd order n=2k+1 it consequently has exactly one fixed vertex f and k transposed pairs.

Consequence essential to the stronger theorem: forgetting rho gives a bijection from isomorphism classes of pairs (T,rho) to self-converse tournament isomorphism classes. This is why a construction on marked anti-involution frames genuinely descends to unmarked classes.

## 3. Hollow frame and odd automorphism group

For each pair choose its unique labeling (x_i,y_i) with x_i -> f -> y_i and rho(x_i)=y_i. Delete only the k pair edges {x_i,y_i}. Call the resulting directed partial graph D.

The x-subtournament A and the off-diagonal matrix b_ij=[x_i -> y_j] determine D. Anti-symmetry says b_ij=b_ji and determines the y-subtournament as the rho-reversed copy of A. A and the symmetric off-diagonal b are otherwise unrestricted.

The directed graph D intrinsically determines f and rho: f is the unique vertex incident to no missing edge, and the missing edges are the rho pairs. Any automorphism of D fixes f, permutes the x_i among themselves, and carries y_i correspondingly. Thus

    H = Aut(D) = Aut(A,b) <= Aut(A).

This action on pair indices is faithful. Therefore |H| is odd, and every H-orbit O on the pair indices has odd cardinality. It is harmless to define H using the marked frame (D,rho,f), but the markings are redundant.

## 4. Completions, stabilizers, and symmetric paths

Encode the deleted matching orientations by v in {0,1}^k, taking v_i=[x_i -> y_i]. Every v is a valid tournament completion T_v with rho anti-involutory. Isomorphism classes of marked completions of D are exactly {0,1}^k/H. By the conjugacy lemma, these are also the unmarked self-converse classes with hollow-frame type D; distinct hollow-frame types cannot overlap on an unmarked class.

The stabilizer H_v is exactly C_Aut(T_v)(rho). A permutation commuting with rho fixes f and preserves D; preserving the matching edges says precisely that it fixes v. Conversely, every H_v element is such a centralizing tournament automorphism.

In an odd-length rho-symmetric Hamiltonian path, x_i and y_i occupy positions j and n+1-j. These positions cannot be adjacent because n is odd. Consequently no such path uses a deleted matching edge. The set

    P(D) = P_rho(T_v)

is independent of v.

The action of H on P(D) is free, since a permutation fixing an ordered Hamiltonian path fixes all vertices. Hence h(D)=|P(D)| is divisible by |H|. Anti-Rédei, or its direct twisted-involution proof, gives h(D) positive and odd. Thus b0(D)=h(D)/|H| is a positive odd integer, and

    B(T_v) = h(D)/|H_v| = b0(D) * |H.v|.

This optional formula is stronger than required for the pairing.

## 5. Why the centralizer formula is exact

A blue HP orbit consists of words p whose positional reversal alpha_p is an anti-involution. By conjugacy, some representative has alpha_p=rho, so every blue orbit meets P_rho. If p,q in P_rho lie in one G-orbit, q=g p, then

    rho = alpha_q = g alpha_p g^{-1} = g rho g^{-1}.

Thus g lies in C_G(rho). The intersection is exactly one centralizer orbit, and that orbit is free. This proves B=h_rho/|C_G(rho)| without confusing fixed points upstairs with fixed quotient orbits.

## 6. Canonical single-orbit switch

For a precise canonical selection, fix a total binary encoding of rooted pairs (A,b;i). Minimize this encoding over all relabelings that send root i to index 1. Equality of the rooted canonical codes holds exactly when the roots belong to the same H-orbit. The indices whose rooted code is lexicographically least therefore form exactly one H-orbit O(D). Its size is odd, and its selection is equivariant under frame isomorphism.

Set v'=v+1_O over F2. Since 1_O is H-invariant,

    H_{v'}=H_v.

The set P(D) is unchanged, so B(T_v')=B(T_v). The same frame D and same orbit O are used after the switch, so this is an involution on marked completions. The preceding conjugacy and equivariance arguments make it independent of the chosen rho and representative T, hence an involution on unmarked SC classes.

## 7. Cyclic-triangle parity

The action of rho on unordered triples pairs all nonfixed triples. Because rho reverses arcs, it preserves whether a triple is cyclic, so each nonfixed triple orbit contributes an even number to c3.

The only fixed triples are {x_i,f,y_i}. Since x_i -> f -> y_i, this triple is cyclic precisely when y_i -> x_i, that is, v_i=0. Therefore

    c3(T_v) = k - sum_i v_i (mod 2).

In particular switching an H-orbit O changes c3 parity by |O|=1 mod 2. The resulting class involution has no fixed points and pairs the two c3 shores with identical B-values. Every odd order n>=3 is covered. n=1 is correctly excluded because there is no pair orbit to switch.

## 8. Free cube strengthening

Let O_1,...,O_r be all H-orbits on the pair indices. Each indicator 1_{O_j} is H-invariant, so the translations by these indicators descend to a commuting (C2)^r-action on {0,1}^k/H, preserving stabilizers and B.

Define the parity vector

    p([v])_j = sum_{i in O_j} v_i (mod 2).

It is well-defined on H-orbits. Toggling O_j translates p by the j-th standard basis vector, because |O_j| is odd. Therefore no nonidentity cube translation can fix [v]. The action is free.

Each cube contains 2^r distinct unmarked SC classes, all with one B-value, with precisely 2^{r-1} classes on each c3 shore. Since k is the sum of r odd orbit sizes, r=k mod 2. If n=1 mod 4 and n>=5, k is positive even, so r is positive even and hence at least two. Every blue-value multiplicity on each c3 shore is therefore even. Equivalently each blue-value multiplicity among all SC classes is divisible by four at these orders.

A fixed hollow stratum also has total blue mass b0(D)2^k, with b0(D)2^{k-1} on each shore, because the H-orbit sizes sum to 2^k.

## 9. Scope cautions and concrete counterexample

These auxiliary cubes need not be subgraphs of the actual antipodal blue metagraph. The switch is not the all-tile antipode. It preserves the chosen centralizer, not necessarily the full automorphism group, total HP count, total tiling mass, or pure-blue status.

The canonical switch on C3 -> 1 -> C3 at n=7 must toggle all three matching edges because H=C3 acts transitively on the three pair indices. In the current census:

    source node271: B=1, c3=2, H=9, |Aut|=9, tiling mass=1, pure blue;
    switched node267: B=1, c3=11, H=135, |Aut|=3, tiling mass=45, mixed.

The source's sole actual antipodal blue neighbor is node21, which has B=5 and c3=5. Thus its switched partner is not even an adjacent blue node. This is an explicit warning against interpreting the class pairing as an edge-respecting perfect matching.

## 10. Independent finite source checks

Reading census_n3/5/7.json independently gave equal shore spectra:

    n=3: each shore {1:1}.
    n=5: each shore {1:2, 3:2}.
    n=7: each shore {1:5, 3:6, 5:8, 7:16, 9:9}.

The C3 -> 1 -> C3 example above was reconstructed explicitly, its three matching edges reversed, and both source and target classified using hp_fibre from census.py. These checks support but are not substitutes for the all-order proof.

Files read: session14-metagraph-census/census.py and census_n3/5/7.json; math/01-canon/theorems/THM-643-blue-black-line-parity-structure.md; THM-644-gridsym-fiber-law-and-anti-redei.md; THM-647-anti-redei-theorem.md. No repository or remote files were modified.
