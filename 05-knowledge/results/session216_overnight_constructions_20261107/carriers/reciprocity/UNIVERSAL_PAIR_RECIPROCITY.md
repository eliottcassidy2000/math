# Rooted reciprocity obstructs universal paired sources

Status: PROVED, independently audited. This restricts the universal all-O/Z/J paired-source route; it does not obstruct arbitrary paired families or settle the Hamiltonian-path spectrum.

## Theorem

Suppose physical rooted tournaments U1,U2 supply a universal synchronized O/Z/J pair within the Session203 criterion: all eight nonstart coordinates (h,e,L,R,d,C,b,T) agree and s1=s2+2. Then:

1. The orders of U1,U2 have OPPOSITE parity. In particular such a pair cannot be two roots of the same tournament, two root cuts on the same deletion tournament, or any pair of tournaments of equal order parity.
2. Both starts are even. Under width8 the only possible ordered start pairs are (4,2), (10,8), (16,14), (22,20).
3. Every one of their eight shared nonstart coordinates is divisible by3. More sharply, L is congruent to6 modulo12, C and b are divisible by12, and T is divisible by6. In particular h,d are positive odd multiples of3.

These are necessary conditions, not a nonexistence theorem for all universal pairs. They explain why same-graph deletion-hypomorphic constructions cannot realize the desired start gap2, even when their deletion and global profiles agree automatically. A common O/Z/J prefix preserves the relative order parity and both starts, so it cannot remove these particular obstructions from a fixed pair.

## 1. Integer-coefficient rooted cover polynomials

For a rooted n-vertex tournament U define F_U(k) as the number of ordered spanning k-path covers with the first bin a nonempty path STARTING at the root. The other k-1 bins may be empty. Define G_U(k) similarly for the first path ENDING at the root. For positive integers k,

F(1)=s, F(2)=R, F(3)=T,
G(1)=e, G(2)=L.

Both extend to polynomials in Z[k]. Indeed, group covers by the number j of nonempty other bins. Choosing and ordering these bins gives a falling factorial (k-1)_j times an integer: the root path is distinguished, while the other nonempty paths are an unordered set. Equivalently the ordered-cover coefficient is divisible by j!, because permuting distinct nonempty paths acts freely.

## 2. Rooted reciprocity, with a self-contained proof

For every rooted tournament of order n,

F(1-k)=(-1)^(n-1) G(k).                         (1)

Work in the commutative squarefree ring Z[x_1,...,x_n]/(x_i^2). Let A be the zero-diagonal adjacency matrix, X=diag(x_i), x the variable column, and 1 the all-ones column. Put M=I-A X. It is invertible because AX is nilpotent in the squarefree ring. Define

P=1+x^T M^(-1)1,
S_r=x_r [M^(-1)1]_r,
E_r=x_r [(I-A^T X)^(-1)1]_r.

These are respectively the generating series of all directed paths including the empty path, paths starting at r, and paths ending at r. Repeated vertices vanish in the squarefree ring, so matrix walks count exactly simple paths.

The tournament identity A+A^T=11^T-I gives

I+A^T X = M-X+1 x^T = M(I-X)+1 x^T.

Since X^2=0, the inverse of I-X is I+X. Applying the rank-one inverse formula to B=M(I-X), one obtains

(I+A^T X)^(-1)1 = (I+X) M^(-1)1 / P.

Here the denominator is P because x^T(I+X)=x^T. Multiplication by x_r annihilates the extra diagonal x_r term. Therefore

E_r(-x)=-S_r(x)/P(x).                         (2)

Reversing every directed path leaves its squarefree monomial unchanged, so the unrooted path series for A equals that for A^T. The same inverse identity gives

P(-x)=1-x^T(I+A^T X)^(-1)1=1/P(x).            (3)

All divisions are by units with constant term1, so these are identities in the indicated ring. Cover enumeration gives

F(k)=[x_1...x_n] S_r P^(k-1),
G(k)=[x_1...x_n] E_r P^(k-1).

Replace x by -x in G(k), use (2),(3), and note that top-degree coefficients acquire factor(-1)^n. This gives (-1)^n G(k)=-F(1-k), proving (1). Polynomial continuation extends the identity from positive integers to all k.

## 3. The same-parity obstruction

Let two sources have equal e,L,R,T and the same order parity, with common sign sigma=(-1)^(n-1). Set D(k)=F_1(k)-F_2(k), an integer-coefficient polynomial. Reciprocity and equal coordinates give

D(-1)=D(0)=D(2)=D(3)=0.

Thus D is divisible IN Z[k] by the monic polynomial k(k+1)(k-2)(k-3). At k=1 this product equals4. Consequently

s1-s2 is divisible by4.                       (4)

In particular it cannot equal2. This already rules out every same-graph or same-order candidate, independently of h,C,d,b.

## 4. Opposite parity and local divisibility

Now suppose the orders have opposite parity and equal e,L,R,T. Take sigma=(-1)^(n1-1). The same D has D(2)=D(3)=0, so D(k)=(k-2)(k-3)Q(k) with Q in Z[k]. Reciprocity yields

D(0)=2 sigma e=6Q(0),
D(-1)=2 sigma L=12Q(-1).

Hence e is divisible by3 and L by6. Also F(1) congruent F(-1) modulo2 for every integer-coefficient polynomial. Since F(-1)=sigma L, both starts are even. This proves conclusions1 and2, and the initial e,L restrictions.

## 5. Universal propagation forces all eight coordinates divisible by3

A universal pair remains equal in all eight nonstart coordinates after every common O/Z/J word. Its two orders still have opposite parity. Therefore Section4 applies to every such descendant: e_w and L_w are divisible by3 for every word w.

The inherited transfer formulas now give the remaining coordinates, without any realizability assumptions on differences:

- e_Z=d, so d is divisible by3.
- e_J=h-e, so h is divisible by3.
- L_J=C-L, so C is divisible by3.
- L_Z=L+b, so b is divisible by3.
- e_(ZJ)=h+R-e-d, so R is divisible by3.
- L_(JJ)=C+2T-L, so T is divisible by3.

This exhausts the eight coordinates. For the sharper congruences, D(1)=s1-s2=2 gives Q(1)=1. Since Q has integer coefficients, Q(-1) is odd. But Q(-1)=sigma L/6, so L is congruent to6 modulo12. Apply this to the universal descendants Z,J,JJ: L+b, C-L, and C+2T-L must each be congruent to6 modulo12. Hence b and C are divisible by12 and T is divisible by6. The usual oddness of Hamilton-path counts makes h,d odd multiples of3. QED.

## Controls and scope

`verify_reciprocity.py` builds the complete integer polynomials by explicit subset path counting, checks reciprocity and integrality on every labeled tournament through order5 with root0, and checks the same-parity mod4 collision consequence among those literal profiles. All roots are represented up to relabeling in this labeled set. These finite controls audit the formulas; the theorem is proved in the squarefree ring above.

No assertion is made that opposite-parity pairs satisfying the modular conditions exist, or that such a pair would establish eventual spectrum coverage. Their existence would still be only one constructive ingredient.

## Targeted finite search after the theorem

The root-incidence bank contains809 physical bounded-start sources with even s and all eight nonstart coordinates divisible by3. Of these425 pass the stronger cover congruences. A targeted search from all809 sources through every O/Z/J word of length at most7 explores2,653,520 presentations and939,817 full-profile-and-label keys, finding no universal pair. This is only a finite search in this specified menu. In particular it does not exclude a word that first enters the necessary modular class after starting from another source, or sources not represented here. The source and exact receipts are modular_deep_search.cpp and MODULAR_DEEP_COUNTS.txt.

Earlier same-graph shift-chain probes are retained as provenance, but Section3 now rules out that entire universal-pair strategy at every size. The C3 external family tested3060 graphs through order12; the circulant5 external family tested15120 graphs through order12. There is no reason to continue analogous same-graph universal searches after the theorem.
