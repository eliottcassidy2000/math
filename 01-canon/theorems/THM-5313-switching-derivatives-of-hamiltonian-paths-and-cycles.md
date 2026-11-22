---
id: THM-5313
title: Divided powers and switching derivatives of Hamiltonian paths and cycles
status: PROVED
source: codex-2026-11-22
related: [HYP-9669, THM-5139, THM-4524]
---

# Switching derivatives of Hamiltonian paths and cycles

Let \(T\) be a finite labelled tournament. Write \(H\) for its number of directed Hamiltonian paths and \(C\) for its number of directed Hamiltonian cycles counted up to rotation. Switching \(T\) at \(S\) reverses precisely the edges crossing \(S\). For either count put

\[
 \Delta_U f(T)=\sum_{S\subseteq U}(-1)^{|U|-|S|}f(T^S).
\]

**Theorem.** At every tournament and every nonempty set \(U\) of \(d\) vertices,

\[
 v_2(\Delta_U H)\ \ge\ \lceil d/2\rceil,\qquad
 v_2(\Delta_U C)\ \ge\ \lceil d/2\rceil-1 .
\]

Take \(v_2(0)=\infty\). This proves both clauses of [HYP-9669, the switching derivative law](../../05-knowledge/hypotheses/HYP-9669-two-adic-derivative-law-of-the-hamiltonian-path-count-on-a-switching-class.md). It extends the mod-four path and cycle-parity laws of [THM-5139, endpoint graph of a switching class](THM-5139-the-endpoint-graph-of-a-switching-class-hamiltonian-cycle-parity-and-h-mod-four-are-quadratic-on-switching-classes.md). The inequalities are universal lower bounds. Individual switching classes can have higher divisibility.

There is a further factorial bound. Set
\[
 L(d)=v_2((d+1)!)-\lfloor\log_2(d+1)\rfloor \quad(d\ge0).
\]
Then \(v_2(\Delta_U H)\ge L(d)\), and \(v_2(\Delta_U C)\ge L(d-2)\) for \(d\ge2\). For example, every order-fifteen mixed switching difference of \(H\) is divisible by \(2^{11}\). These statements hold for the path and cycle counts on any induced vertex set, using its own switching variables.

At **all primes** one also has
\[
        \lceil d/2\rceil!\mid\Delta_U H,\qquad
        (\lceil d/2\rceil-1)!\mid\Delta_U C.
\]
Consequently, knowing \(H\) on cuts of at most four vertices determines \(H\bmod24\) everywhere on the switching class.

## 1. A rank-two formula for switching

Work over the commutative ring \(\mathcal R=\mathbb Z[z_v:v\in V]/(z_v^2:v\in V)\), with polynomial marking variables \(x_v\). Write \(A\) for the adjacency matrix, \(Z=\operatorname{diag}(z_v)\), \(\mathbf 1\) for the all-ones column, and

\[
 R=Z(I-AZ)^{-1},\quad
 q=\mathbf1^tR\mathbf1,\quad a=\mathbf1^tRx,\quad b=x^tR\mathbf1,\quad c=x^tRx.
\]

The inverse is a finite walk expansion. Repeating a vertex kills its monomial. The coefficient of \(z_W=\prod_{v\in W}z_v\) in \(q\) counts the paths on \(W\), including the singleton path for \(|W|=1\). The other three series mark the last endpoint, the first, and both respectively. They all have zero constant term in \(z\).

For a Boolean specialization \(x=1_S\) put \(s=\mathbf1-2x\), \(D=\operatorname{diag}(s)\), and let \(A_x\) be the switched adjacency matrix. Entrywise on and off the cut,

\[
 DA_xD=A-x\mathbf1^t-\mathbf1x^t+2xx^t=A+PKP^t,\quad
 P=[\,\mathbf1\ \ x\,],\quad K=\begin{pmatrix}0&-1\\-1&2\end{pmatrix}.
\]

Put \(M=P^tRP=\left(\begin{smallmatrix}q&a\\ b&c\end{smallmatrix}\right)\), and \(h=(1,-2)^t\). The determinant lemma and the geometric series give, with \(q_x=\mathbf1^t Z(I-A_xZ)^{-1}\mathbf1\),

\[
 E=\frac{\det(I-A_xZ)}{\det(I-AZ)}
    =\det(I-KM)=(1+a)(1+b)-(2+q)c,\qquad
 q_x=h^tM(I-KM)^{-1}h.
\]

Multiplying the last expression by the two-by-two adjugate gives \(E q_x=q+2-2E\). In particular

\[
            2+q_x=\frac{2+q}{(1+a)(1+b)-(2+q)c}.                 \tag{1}
\]

For cycles use \(-\log\det(I-AZ)=\sum_{j\ge1}\operatorname{tr}((AZ)^j)/j\). Every simple directed \(j\)-cycle has exactly \(j\) starting points. Thus \(-\log E\) generates the change in the cycle count on each support. The logarithm is taken in the rationalized nilpotent ring, and these coefficients are integers at Boolean specializations.

One identity specific to tournaments will supply a missing power of two. Modulo two, \(A^t=A+I+\mathbf1\mathbf1^t\). A diagonal step in the walk series vanishes, and the remaining rank-one update gives \(R^t=R+R\mathbf1\mathbf1^tR/(1+q)\). In particular \(b=a/(1+q)\). Every element of positive \(z\)-order squares to zero in characteristic two, whence

\[
                     ab=0,\qquad q(a+b)=0\quad\bmod 2.          \tag{2}
\]

These identities also hold with the \(x_v\) as ordinary indeterminates. The analogous divisibility question for an independent cube of individual arc flips has a different answer. The rank-two update above uses the common endpoint signs of a cut.

## 2. The divisibility calculation

For \(y\) of positive \(z\)-order, the divided power \(\gamma_j(y)=y^j/j!\) is integral. In its expansion every nonvanishing product uses \(j\) nonempty disjoint supports, and their \(j!\) orderings give identical terms. In particular \(y\gamma_{2i}(y)\equiv\gamma_{2i+1}(y)\pmod2\).

Use \(\mathbb Z_{(2)}\) for divisions by odd integers. For \(s\ge1\) define
\[
 B_s=\frac{(2+q)^s}{s!}
     =\sum_{r=0}^s \frac{2^r}{r!}\gamma_{s-r}(q).
\]
Since \(v_2(r!)<r\) for \(r>0\), \(B_s\) is 2-integral and \(B_s\equiv\gamma_s(q)\pmod2\). For odd \(s\) this reduction has a factor \(q\).

**Coefficient lemma.** Give \(a,b\) degree one in the marks. For every \(t\ge0\), the total degree-\(t\) part of
\[
 B_s(1+a)^{-s}(1+b)^{-s}
\]
is divisible by \(2^{\lceil t/2\rceil}\).

Indeed, omitting the common sign \((-1)^t\), the \(i,j\) summand for \(i+j=t\) is
\[
 B_s\,s^{\overline i}s^{\overline j}\gamma_i(a)\gamma_j(b),\qquad
 s^{\overline i}=s(s+1)\cdots(s+i-1).
\]
If \(s\) is even, the rising products already supply at least \(\lceil i/2\rceil+\lceil j/2\rceil\) factors of two. If \(s\) is odd, each supplies at least the corresponding floor.

For even \(t\), only two odd indices can be short by one factor. After extracting the floor powers, the two odd divided powers contribute \(ab\) modulo two, which vanishes by (2). For odd \(t\), group the indices as \((2u+1,2v)\) and \((2u,2v+1)\). Their common rising product supplies \(2^{u+v}\). The remaining rising factors \(s+2u\) and \(s+2v\) are odd, so after that division the sum contains \(B_s(a+b)\) modulo two. This vanishes by (2) and the factor \(q\) in \(B_s\). The proof includes the empty rising products and divided powers at index zero.

Expand (1) in \(c\). Its \(c^k\) term is
\[
 (k+1)!\,k!\ \gamma_k(c)\,B_{k+1}(1+a)^{-(k+1)}(1+b)^{-(k+1)}.
\]
The scalar factorials contain at least \(k\) factors of two, by counting even factors. As \(c\) has mark degree two, the degree \(d=2k+t\) piece has valuation at least \(k+\lceil t/2\rceil=\lceil d/2\rceil\).

For cycles write
\[
 -\log E=-\log(1+a)-\log(1+b)
 +\sum_{k\ge1}\frac{(2+q)^k c^k}{k(1+a)^k(1+b)^k}.
\]
The degree-\(t\) pure logarithms have factor \((t-1)!\), giving at least \(\lceil t/2\rceil-1\) factors of two. The mixed term has scalar \(k!(k-1)!\) multiplying \(\gamma_k(c)B_k(1+a)^{-k}(1+b)^{-k}\). This contributes at least \(k-1\) powers of two before the coefficient lemma. This proves the cycle bound at degree \(2k+t\).

Finally specialize \(x_v^m=x_v\) for \(m>0\). A coefficient supported on \(U\) only receives terms of ordinary degree at least \(|U|\). The two bounds are monotone, and the resulting multilinear coefficient is exactly \(\Delta_U\). This completes the path and cycle proofs.

## 3. The factorial refinement

First consider arbitrary primes. The polynomial \(E-1\) has positive \(z\)-order and mark degree at most two. In the geometric series for \(E^{-1}\), its \(j\)-th power is divisible by \(j!\) by the divided-power argument. In \(-\log E\) the same term is divisible by \((j-1)!\). A monomial of mark degree \(d\) requires \(j\ge\lceil d/2\rceil\). This proves the two integer factorial divisibilities above. In particular on every switching class \(H\bmod p\) has degree at most \(2p-2\) for every prime \(p\). Equivalently in characteristic \(p\), \(E^p=1\) and (1) gives \(2+q_x=(2+q)E^{p-1}\).

Keep the full rising products before using (2). The integral divided powers and \(B_s\) leave the explicit scalar factors

\[
 (k+1)(k+i)!(k+j)!\quad\text{for paths},\qquad
 k(k+i-1)!(k+j-1)!\quad\text{for mixed cycle terms}.
\]

For \(d=2k+i+j\), the factorial arguments sum to \(d\) and \(d-2\), respectively. A pure logarithm contributes \((d-1)!\). It therefore suffices to identify \(\min_{u+v=n}v_2(u!v!)\). This minimum is nondecreasing, since deleting a factorial factor reduces \(n\) by one. Thus it also survives Boolean specialization.

Writing \(m=\lfloor\log_2(n+1)\rfloor\), the maximum of \(v_2\binom n u\) is \(m-v_2(n+1)\). In the factorial floor formula for the valuation each power \(2^j\le n+1\) contributes at most one, and contributes zero whenever it divides \(n+1\). Equality holds at \(u=2^m-1\), from
\(\binom n {2^m-1}=2^m\binom{n+1}{2^m}/(n+1)\).
The binomial on the right is odd since the top binary digit of \(n+1\) is \(2^m\) (equivalently expand \((1+X)^{n+1}\) modulo two). Subtracting the maximum from \(v_2(n!)\) gives exactly \(L(n)\).

## 4. Recovering counts across a switching class

Fix one vertex's switch bit at zero to remove the global double representation. Modulo \(2^r\), every value of \(H\) on the resulting cube is determined by its values on subsets of size at most \(2r-2\), by Boolean Möbius inversion. The radius improves to the largest \(d\) with \(\max(\lceil d/2\rceil,L(d))<r\), allowing \(d=0\). For example, modulo 2048, switches of at most 14 vertices suffice on tournaments of **any order**. The corresponding radius \(2r\) works for \(C\). If a radius exceeds the dimension of the cube one just uses the whole cube.

Also, on every set of \(2k\) vertices the parity of \(2^{-k}\Delta_U H\) is a switching-class invariant. Switching an unmarked vertex changes the difference by a derivative of order \(2k+1\); switching a marked one negates it. The original binary-digit consequence of HYP-9669 follows as well: the \(k\)-th bit of \(H\), starting with the units bit, has Boolean degree at most \(2^k\). Indeed, choose nonnegative representatives of multilinear coefficients modulo \(2^{k+1}\), and split them into binary summands \(2^j x^U\). For every nonconstant summand \(|U|\le2j\le2^j\). Over \(\mathbb F_2\), extract the coefficient of \(t^{2^k}\) from the product of factors \(1+x^U t^{2^j}\); it equals the desired bit by the binomial identity. Every contributing monomial has degree at most \(2^k\).

These reconstructions predict residues of the ordinary total in a specified switching class. The choice of a class and a switch attaining each prescribed odd integer other than 7 and 21 remains the open Chan–Pak problem. The earlier absolute Walsh degree \(2\lfloor |V|/3\rfloor\) provides a separate order-dependent cutoff. [THM-4524, loop gauge and Walsh expansion](THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md)

## 5. Multiplicity of a prescribed residue

Let \(m=|V|-1\), and let \(N_h\) count the distinct labelled cuts in this class with \(H\equiv h\bmod 2^r\). For \(r\ge2\), set \(D=2^{r-1}\) and \(W=2^r-2\). Then

\[
  2^{\max(0,\lceil(m-W)/D\rceil)}\mid N_h,\qquad
  N_h>0\ \Longrightarrow\ N_h\ge 2^{\max(0,m-W)}.
\]

In particular every attained odd residue modulo eight occurs on at least \(1/64\) of the labelled cuts in its class. The class can still miss a residue entirely. Multiplicities here count cuts before identifying isomorphic tournaments.

**Proof.** An even target gives an empty fiber by Rédei parity. For an odd target the fiber is the common zero set of \(r-1\) Boolean polynomials, one for each nonconstant digit condition, with degrees bounded by \(2,4,\ldots,D\). Its indicator therefore has Boolean degree at most \(W\). A nonzero Boolean polynomial of degree at most \(w\le m\) has support at least \(2^{m-w}\): choose a maximum-degree monomial, fix the variables outside it arbitrarily, and observe that its coefficient stays one. Each such restriction has a nonzero evaluation. This proves the size bound.

For divisibility, write each digit condition's zero indicator as an integer-valued product expansion. If a polynomial is the parity of its monomials \(M_i\), its zero indicator is \((1+\prod_i(1-2M_i))/2\). Selecting \(\ell>0\) monomials from this expansion brings a factor \(2^{\ell-1}\). In a term of the product of all the indicators put \(s=\sum_{\ell_j>0}(\ell_j-1)\), and let \(u\) coordinates be absent from its monomial. Summing over the cube contributes another \(2^u\). The covered coordinates number at most \(\sum_{\ell_j>0} 2^j\ell_j\le W+Ds\), where digit indices begin at one. Hence \(m\le W+Ds+u\), giving \(s+u\ge\lceil(m-W)/D\rceil\) for every term. This elementary Boolean counting argument applies independently of any estimate for the frequency of individual integer counts.

**Finite witness.** The nine-vertex source in the accompanying script has \(H=1465,\ C=73\). The face on \(\{1,2,3,5\}\) gives mixed differences \(396\) and \(42\), respectively, attaining both 2-adic lower bounds in degree four. All 256 cut residues for \(H\bmod24,\ C\bmod4\) are reconstructed from the 163 cuts of radius four. Every odd residue modulo 24 occurs in this one class. The least actual counts in its residues 7 and 21 are 511 and 189. The program counts ordinary paths and cycles directly by subset dynamic programming. [Reproduction script](../../05-knowledge/results/chanpak_switching_divided_powers_20261122/check_switches.py)

## 6. Exact odd-prime divisibility at fixed tournament order

There is a sharp answer away from two. Fix \(n\ge3\), an odd prime \(p\), and a specified \(d\)-set \(U\) of vertices, \(d\ge1\). Put \(r=\lceil d/2\rceil\). If \(3r>n\), both mixed differences vanish for every \(T\). Otherwise let \(b\) be the least integer at least \(r\) having the parity of \(n\). Then, over **all labelled tournaments of this order**,

\[
 \min_T v_p(\Delta_U H(T))=v_p(b!),\qquad
 \min_T v_p(\Delta_U C(T))=v_p((b-1)!).
\]

In particular the maximum possible switching-cube degree of \(H\bmod p\) is exactly \(2\min(\lfloor n/3\rfloor,b_p)\), where \(b_p=p-1\) for even \(n\) and \(b_p=p-2\) for odd \(n\). For modulus three the degree ceiling at every odd order is two; at every even order at least six it is four. Both ceilings are attained by some class at each stated order. Thus at odd order the residues modulo **12** of all switchings are determined by the radius-two ball. The minimum statement allows the minimizing source to depend on \(p\) and on the observable.

**Proof.** Fix reference directions and write \(\chi_F(T)\) for the product of edge signs on \(F\). The Walsh expansion of the path count has support on linear forests all of whose paths have even edge length. For such a forest with \(j\) nontrivial paths and \(c\) components including isolated vertices, the coefficient has magnitude \(2^{j-(n-1)}c!\). Order the \(c\) blocks and independently orient each nontrivial block. An odd-length block cancels by reversal. This is the inherited forest expansion. [THM-4761, Walsh chaos of the Hamilton path count](THM-4761-random-tournaments-the-forest-identity-is-the-walsh-chaos-decomposition-of-kappa.md)

For cycles the same forest has coefficient of magnitude \(2^{j-n}(c-1)!\), by cyclically ordering the blocks. The signs of these two coefficients agree for a common reference orientation. The only additional possible cycle coefficient is a whole spanning cycle, whose vertex degrees are all even. Indeed any proper subset of the edges of a spanning cycle is a linear forest. These counts also apply when a single nontrivial block spans the vertices.

Under vertex switching the edge character on \(F\) acquires a factor equal to the product of vertex spins on its odd-degree set \(\partial F\). Taking the mixed difference multiplies its coefficient by \((-2)^d\) if \(U\subseteq\partial F\), and kills it otherwise. In particular a spanning cycle contributes nothing. On an even linear forest the odd-degree vertices are its \(2j\) endpoints, so \(j\ge r\), \(c\ge j\), and \(n\ge 3j\). Every component has odd order, giving \(c\equiv n\bmod2\).

Conversely, when \(n\ge3r\) there is such a forest with \(j=r,c=b\). Take \(r\) two-edge paths, \(b-r\) isolated vertices, and lengthen one path in pairs using the remaining \(n-(2r+b)\) vertices. The remainder is nonnegative and even. Relabel the vertices so that the endpoints include \(U\). All contributing coefficients are thus divisible over \(\mathbb Z_{(p)}\) by the claimed factorials, and this forest attains their exact valuations. The change of basis between edge characters and function values over all labelled tournaments is invertible over \(\mathbb Z_{(p)}\): Walsh inversion divides only by a power of two. The two lists generate the same ideal. This proves both minima and the modular-degree formula. Fixing one switch bit at zero preserves the upper bound; at the claimed even top degree the endpoint set can be chosen away from that fixed vertex, so it preserves attainability as well.

**Recovering the path and cycle relation.** More generally, for integral \(a\), the exact minimum for \(\Delta_U(H-2aC)\) at odd \(p\) is the minimum of \(v_p((c-a)(c-1)!)\) over
\[
                  b\le c\le n-2r,\qquad c\equiv n\bmod2.
\]
The same construction realizes every listed \(c\) with exactly \(r\) nontrivial paths, and path and cycle coefficients on that forest differ by the factor \(2c\). For example, at the first possible order \(n=3r\) there is only \(c=r\), so \(\Delta_U H=2r\,\Delta_U C\) identically over the integers. The cycle records the last-to-first edge that disappears when a block order is read linearly.

This gives a three-way boundary for integer interpolation. Write \(n=3r+t\), \(0\le t\le2\), \(r\ge1\). At \(t=0,1\), the observable \(H-2(r+t)C\) has switching degree at most \(2r-2\) on **every** source: for \(d=2r-1,2r\) its only admissible component count is \(r+t\); above that range both derivatives vanish. At \(t=2\) both \(r\) and \(r+2\) occur at degree \(2r\), so no constant multiple of \(C\) cancels the top degree of \(H\) on all sources. This distinction comes from fitting one extra pair of isolated vertices alongside the same marked endpoints.
