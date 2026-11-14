# Busch floors, Fermat-target rigidity, and marked inverse fibers

**Proved below:** rigidity for entire arithmetic families, including composite Fermat-number targets, and explicit inverse fibers with five neutral parameters. **Inherited:** Busch's sharp strong floor, orbital path freeness, and the nearly transitive realization. These results narrow a constructive stratum; they do not prove a global Chan–Pak threshold.

## 1. Primary and repository sources

- Arthur H. Busch, [*A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3](https://www.kurims.kyoto-u.ac.jp/EMIS/journals/EJC/Volume_13/PDF/v13i1n3.pdf). Use Theorem 1 and Corollary 1, not only the abstract's shorthand. For strong order \(n\ge3\), the exact minimum is \(3\beta^{n-3},\beta^{n-1},9\beta^{n-5}\) in residue classes \(0,1,2\) modulo three, where \(\beta=5^{1/3}\). In particular it is at least \(5^{(n-1)/3}\). Lemmas 1–2 distinguish structural overlapping-core lower-bound injections from equalities of counts or ready-interface gluing.
- [Session 133's finite core budget](../session133_fractran_spectrum_refinement_20261101/core_reduction/README.md) records the sharper exact inverse budget. [Fermat count packet](../fermat_monster_hamiltonian_counts_prime_classes_20261102/README.md) gives least orders \(3,4,6,8,11\) for the five known Fermat primes and the vertex-transitive prime-count obstruction; its explicit examples are not needed for the rigidity proof.
- The requested **Neményi** is P. F. Neményi's inverse elasticity problem for stress-trajectory nets. His primary survey is [*Recent Developments in Inverse and Semi-Inverse Methods in the Mechanics of Continua*](https://www.sciencedirect.com/science/article/pii/S0065215608703004); the precise modern source is A. W. Sáenz, [*Determination of stresses from their stress trajectories in plane elastic systems: the five constant theorem*](https://link.springer.com/article/10.1007/s00707-008-0118-3), *Acta Mechanica* 208 (2009), 215–225. The modern abstract gives **two through five** real parameters under its hypotheses, not five freely specifiable parameters for arbitrary observed fields and boundary conditions. It cites Neményi's 1933 *Strömlinien und Hauptspannungstrajektorien* and the original Neményi–Sáenz plane stress paper.

## 2. A rough-target rigidity certificate

Write \(H(T)\) for the number of directed Hamiltonian paths. In a tournament, the strong components are linearly ordered, and every spanning path concatenates a spanning path of each component in that order. Thus \(H(T)=\prod_C H(C)\). The automorphism group is likewise the product of those of the components; their positions are fixed by every automorphism.

If a prime \(p\) divides \(|\operatorname{Aut}(C)|\), then \(p\le|C|\) by the faithful vertex action, and \(p\mid H(C)\) because the action on ordered Hamiltonian words is free. Busch now gives
\[
  5^{(p-1)/3}\le 5^{(|C|-1)/3}\le H(C)\le H(T).
\]
Consequently, if \(h>1\) has least prime factor \(P\) and \(h<5^{(P-1)/3}\), **every tournament with \(H=h\) is rigid**. More sharply one can test each individual \(p\mid h\) against the exact Busch floor at order \(p\). This does not assume the tournament is strong or that the original tournament order is small: adding singleton padding does not create automorphisms.

For prime \(h=p\ge5\), one has \(p<5^{(p-1)/3}\). For example \(5^4>5^3\), and \(5^{n-1}/n^3\) strictly increases for all integer \(n\ge5\), since the ratio at consecutive \(n\) is \(5(n/(n+1))^3>1\). Thus every prime-count realizer at \(p\ge5\) is rigid; there is exactly one nonsingleton SCC. At \(h=3\) the single core is \(C_3\) and is not rigid, so the exceptional boundary is real.

**Fermat-number targets, even if composite.** Put \(m=2^r\) with \(r\ge1\). A prime factor \(p\) of \(2^m+1\) has \(\operatorname{ord}_p(2)=2m\): its order divides \(2m\) but not \(m\), and \(2m\) is a power of two. Hence \(p\ge2m+1\). On the other hand, \(m\) is even and
\[
2^m+1=4^{m/2}+1\le5^{m/2},\qquad
|C|\le1+3\log_5(2^m+1)\le1+3m/2<2m+1.
\]
No prime divisor of a putative automorphism group can occur. This establishes rigidity for **all** realizers of every \(F_r\), \(r\ge1\), without deciding Fermat primality.

There is always a strong realizer: take the transitive tournament of order \(m+2\) and reverse only its extreme arc. Exactly one spanning path avoids the reversed arc; a path using it is determined by the subset of interior vertices before its tail, traversed in increasing order there and after its head. The count is \(2^m+1\). This construction is inherited from Busch/Moon; the point is that **the entire exact-count fiber** is rigid, including unrelated realizers, not just this one construction. Compare [THM-4644](../../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md): rigidity in its Fermat-*order* statement concerns the multiplier action on circulant presentations, not absence of translation automorphisms of a circulant graph. The two claims have different inputs and different groups.

**Signed cyclotomic companions.** For odd prime \(\ell\), every prime dividing \(2^\ell-1\) has order \(\ell\), so is \(1\) modulo \(2\ell\). For odd prime \(\ell\ge5\), every prime dividing \((2^\ell+1)/3\) has order \(2\ell\): the only possible other order is 2, which would mean prime 3, but \(9\nmid 2^\ell+1\) since \(\ell\equiv1,5\pmod6\). Thus their least prime factors are likewise at least \(2\ell+1\). Both targets are less than \(2^\ell\), and \(2^3<5^2\) gives \(1+3\log_5 h<1+2\ell\). Every realizer of either family is rigid, regardless of target primality. This says nothing about existence: the first Mersenne example \(2^3-1=7\) is an excluded count.

## 3. Five stress coordinates, and what “same observation” forgets

Here is a self-contained local check in the simple Cartesian isothermal net. In body-force-free homogeneous isotropic plane elasticity, use the standard Airy formulation: a biharmonic function \(\phi\) gives stresses \(\sigma_{xx}=\phi_{yy}\), \(\sigma_{yy}=\phi_{xx}\), \(\sigma_{xy}=-\phi_{xy}\). Requiring the Cartesian axes as principal directions forces \(\sigma_{xy}=0\). On a rectangle the equilibrium equations therefore give \(\sigma_{xx}=A(y)\), \(\sigma_{yy}=B(x)\), and biharmonicity gives \(A''(y)+B''(x)=0\). Each term is constant. Conversely the resulting stress fields
\[
 \sigma_{xx}=a y^2+b y+c,\qquad
 \sigma_{yy}=-a x^2+d x+e,\qquad \sigma_{xy}=0
\]
come from integrating twice in \(y\) and \(x\) separately to obtain a biharmonic Airy function. This is an exact five-real-parameter inverse fiber for the fixed **local model**, without prescribed tractions. To keep distinct principal values near the origin, take \(c\ne e\) and restrict the neighborhood; where they coincide the direction observer degenerates. Three distinct \(y\)-samples of \(\sigma_{xx}\) recover \(a,b,c\); two distinct \(x\)-samples of \(\sigma_{yy}\) recover \(d,e\) after \(a\) is known. Thus the observer loses five parameters, but they can be recovered by five compatible scalar marks. This does not reprove the general curved-net theorem.

The tournament analogue can be made exact without claiming a PDE reduction. Fix an **ordered** chain \(C_1,\ldots,C_k\) of nonsingleton strong components. The isomorphism classes obtained by allowing any number of singleton SCCs between and around them are exactly
\[
\mathbb Z_{\ge0}^{k+1}\ni(g_0,\ldots,g_k)
 \longmapsto TT_{g_0}\to C_1\to TT_{g_1}\to\cdots\to C_k\to TT_{g_k}.
\]
Here \(TT_0\) is empty. This is a bijection because the SCC order is a native, unique decoder. It preserves \(H=\prod_i H(C_i)\) and the full automorphism group, not just a scalar approximation. For fixed size \(N\), letting \(N_0=\sum_i|C_i|\), the ordinary count is \(\binom{N-N_0+k}{k}\) for \(N\ge N_0\). This is the precise difference between five unconstrained coordinates and a fixed-size slice with one linear constraint.

Take the four nearly transitive cores of counts \(5,17,257,65537\) and orders \(4,6,10,18\). Every core is rigid by §2, so the whole family is rigid. The common path count is
\[
5\cdot17\cdot257\cdot65537=(2^{32}-1)/3.
\]
The five gap lengths give an explicit discrete five-parameter family of nonisomorphic tournaments at that one count. At order \(N\ge38\) this subfiber has exactly \(\binom{N-34}{4}\) members. The underdetermination is not caused by graph symmetries. Conversely at a prime count there is only one nonsingleton SCC, so **only two** such padding coordinates per fixed core; one cannot import five available neutral directions to a prime-target proof by analogy.

The common principle is a **marked inverse fiber with a compatibility condition**: one must first show that the invisible parameters actually span lawful solutions, then specify the sidecar that decodes them. Degree counts, generic dimensions, or numeric equalities alone do not transport native endpoint readiness. Here the stress compatibility is harmonic trace and the tournament compatibility is a fixed ordered SCC chain. The constructions preserve the declared observation, not arbitrary prescribed tractions or arbitrary path endpoint constraints.
