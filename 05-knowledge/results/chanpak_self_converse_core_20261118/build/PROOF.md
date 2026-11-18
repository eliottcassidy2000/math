# A strong self-converse embedding with an exact septimal dimension register

Let H(T) be the directed Hamiltonian path count. Let (P_j(T)) count ordered j-slot path covers allowing empty paths, including the unique cover on the empty tournament. The nonempty-slot versions are (Q_j(T)). In particular (h=H(T), a=P_2(T), b=P_3(T)), and (Q_2=a-2h, Q_3=b-3a+3h). These identities count the choices of empty slots. All profiles (P_j) multiply under ordinal sums: each slot visits its left part and then its right part, and the two covers are independently chosen, empty parts included. For a transitive m-block (P_j=j^m), assigning each vertex to one of j ordered slots.

Write (C_3[X,Y,Z]) for cyclic substitution, with all cross arcs (X→Y→Z→X). It is strong if the three blocks are nonempty, regardless of their internal arc orientations: the cycle of modules provides directed paths to any required vertex. If (Y=X^*), the converse of X, and Z is transitive, this tournament is self-converse by the explicit anti-automorphism swapping corresponding X/Y vertices and reversing the Z order. In particular it contains the exact induced source X and its converse.

## One center vertex: norm and exact extraction

For any nonempty T, use (X=T→TT_m), (Y=X^*), and one center vertex. Let (R=Q_2(X)). On contracting maximal same-block runs, the path follows the cyclic module order. The center is visited once; the possible multiplicities ((#X,#Y)) are (1,1), (2,1), (1,2), (2,2), with respectively 3,1,1,1 choices of the run skeleton's starting module. Each ordered list of runs is filled independently by nonempty ordered covers of its block. Consequently
[
M_m=3h^2+2hR+R^2=(2^m a-h)^2+2h^2.
]

This gives exact recovery of the two source scalars from three physical reads:
[
h=\sqrt{M_2-6M_1+8M_0}/3,quad
a=\sqrt{(M_2-3M_1+2M_0)/6}.
]
Their orders are (2|T|+1,2|T|+3,2|T|+5). In fact two reads suffice when (|T|≥4): putting (r=a/h), the ratio (M_1/M_0=(4r^2-4r+3)/(r^2-2r+3)) strictly decreases for (r≥5) (derivative numerator (-4r^2+18r-6)). The bound (r≥|T|+1) comes from splitting every full path at each of its (|T|+1) boundaries into two slots. Those two paths and their order remember exactly which full path and cut were used. The unique root in (r≥5) gives (h^2=M_0/(r^2-2r+3)).

One sufficiently large read also suffices: if (B=2^m≥h), then (x=Ba-h≥h^2) by (a≥2h), and (x^2+2h^2<(x+1)^2). Thus (x=\lfloor\sqrt{M_m}\rfloor), and (h=\sqrt{(M_m-x^2)/2}), (a=(x+h)/B). If only source order is known, take (2^m≥|T|!).

This norm does **not** directly choose arbitrary prime or 7p counts. Both x and h are odd, so every output is 3 mod8. If it is divisible by7, (x^2+2h^2≡0\pmod7) forces both x and h divisible7 since −2 is a nonresidue: (7)-adic valuation is even by repeatedly removing common factors. For prime output p, at most one positive pair x,h represents it: the two ratios x/h of any pair of representations mod p are equal or negatives, so p divides (xh'−x'h) or (xh'+x'h), whose absolute size is at most (p/\sqrt2<p) by Cauchy–Schwarz. The sum is positive and cannot vanish, so the difference is zero, recovering identical positive coordinates. This is a restricted inverse, not a universality claim.

## Two center vertices: an exact physical 7-adic register

Now take center (TT_2), (X=T→D^j→TT_m) where (D^j) denotes j ordinal copies of (D=(2,4,9,3)) in outgoing masks. Write (h',a',b') for X's profile and (R=a'-2h'), (S=b'-3a'+3h'). The same cyclic run counting gives

    F(h',a',b') = 5h'^2 + 6h'R + 7R² + 4RS + 2S².

Indeed, center multiplicity one gives the same expression as before. With center multiplicity two, the (X,Y) multiplicity list is (1,1),(2,1),(1,2),(2,2),(3,2),(2,3),(3,3), with counts (1,1,1,3,1,1,1). Fill its two nonempty ordered slots in (Q_2(TT_2)=2) ways. This adds (2(h'^2+2h'R+3R^2+2RS+S^2)).

The profile of D is exactly (5,32,117), by the literal slot subset DP. **Criterion for the source-entire ordinal padding model:** arbitrary prescribed multiples of seven modulo any (7^k) can be realized by this specific padding family, from any source with
[
b\not\equiv0\pmod7,qquad (h,a)\not\equiv(0,0)\pmod7.
]
Conversely if the criterion fails, no ordinal padding at all (even allowing arbitrary tournaments in place of D/TT) achieves exact valuation one. To see the converse, ordinal padding cannot revive a zero profile coordinate mod7. If (b'=0\pmod7), the form reduces mod7 to (-(a'^2+h'a'-h'^2)), which by discriminant 5 is zero only when both coordinates vanish. If (h'=a'=0\pmod7), the form reduces to (2b'^2), so a zero again requires all coordinates vanish. In either situation a zero is automatically a multiple of49.

For sufficiency find (0≤j,m_0≤5) by the following cases. For a unit h', use ratios (alpha=a'/h'), (eta=b'/h') mod7.

* h,a,b all units: normalize ((alpha,eta)=(2,4)). The multiplicative controls from D and singleton are ((5,1)) and ((2,3)) mod7. In discrete logarithms base3 they are columns (5,0) and (2,1), invertible mod6. So j and m0 always suffice.
* h,b units, a zero: normalize (eta=1) using only TT padding, since 3 generates the units.
* h zero, a,b units: normalize (b'/a'=1) using only TT padding, since (3/2=5) generates the units.

In the first case the normalized form (F(1,alpha,eta)) vanishes mod7 because R=0,S=1 gives 7. In the second it is likewise zero because R=−2,S=4 gives 21. In the last it reduces to (2(b'-3a')(b'-a')), so the specified phase vanishes mod7.

At any such phase, increasing m by (6·7^e u) keeps h' fixed and multiplies a',b' respectively by
[
1+2u7^{e+1},qquad 1+6u7^{e+1}pmod{7^{e+2}}.
]
This follows by binomial expansion from (2^6≡1+2·7), (3^6≡1+6·7\pmod{49}). Substitution into the quadratic gives a change of (u7^{e+1}) times a unit, namely (4h'^2, 6h'^2, 5a'^2\pmod7) in the three cases respectively. For the first evaluate partials in R,S as 10,4 and note (dR=dalpha, dS=deta-3dalpha); the coefficient is (10·4+4(24-12)≡4). For the second (dR=0) and (F_S=8), giving (8·6≡6). For the third the derivative of the factor (b'-a') is (4a'), and the other factor has value (-2a'), giving (-16a'^2≡5a'^2).

Hence the map (t\mapsto F(m_0+6t)) satisfies the exact distance law
[
v_7(F(m_0+6t)-F(m_0+6u))=1+v_7(t-u)
quad(t,u≥0, t\ne u).
]
In particular for each k the first (7^{k-1}) values of t biject onto the multiples of7 modulo (7^k). This is not just a mod7 existence result: each next digit independently selects one of seven dimension offsets, as a Hensel lift. To achieve v7=1 use k=2, giving (j≤5, m≤41), so total host order at most (2|T|+124). The count grows exponentially in the padding dimension, so this does not approximate or hit an arbitrary absolute target.

## Literal 7p control beyond the one-center obstruction

Take (T=(8,5,9,50,7,23,63)), (j=m=0). Its (h,R,S) profile is ((23,546,3666)). The two-center literal strong self-converse tournament has order16 and count (37,050,461=7·5,292,923) with the quotient prime by trial division. This is a positive physical fix to the inert-prime obstruction of the one-center family, not a proof of infinitely many prime quotients. The shape is explicitly cyclic and anti-isomorphic to itself; no blue-edge incidence is asserted here.

Run `python build.py` for literal full Hamiltonian and anti-automorphism controls, exact source profile, quotient primality, and phase lifting on all residues. Use `python build.py --profile 1 2 3 --target 140 -k 4` for a sample prescribed residue. The integer output pair is (number of D copies, width of the final transitive padding).

## All four odd residues modulo eight from four physical centers

For any nonempty X and any nonempty self-converse center Z, the cyclic mirror (C_3[X,X^*,Z]) satisfies
[
H(C_3[X,X^*,Z])\equiv H(Z)+P_2(Z)\pmod8.
]
In the ordered run sum, if all three blocks have multiplicity at least 2 then the product of three even ordered-cover counts is divisible8 (one can freely swap the first two nonempty slots in each). If any block has multiplicity one, no block can have more than two runs. The full remaining contribution mod8 is (H(Z)(3h^2+2hR+R^2)+h^2Q_2(Z)), where (h=H(X),R=Q_2(X)): terms with Z and just one reflected half repeated have total (2hR Q_2(Z)), which vanishes mod8. If |X|≥2, (P_2(X)) is divisible4 because complement-paired subsets give twice a sum of an even number of odd products by Rédei, so (R≡-2h\pmod4). Thus (3h^2+2hR+R^2≡3\pmod8). The singleton also gives 3 directly from h=1,R=0. Finally (3H(Z)+Q_2(Z)=H(Z)+P_2(Z)). The four centers (TT_1,TT_2,TT_3,C_3) yield residues (3,5,1,7) respectively, all with exact induced X and X* and explicit converse symmetry. Only (TT_2) was shown above to have the stated septimal phase criterion.

