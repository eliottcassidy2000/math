# Root O/Z parameter capacity in every fixed residue class

Status: proof for independent review. This is a root-image obstruction, not a universal finite-bank impossibility theorem.

Write H(m)=1+2m-Q(m), with the established Q recurrence. For any fixed odd positive integer p and a in Z/pZ, define

A_a(X)=#{m>=0:H(m)<=X and H(m)=a mod p}.

Then

A_a(X)=(log(2)/p)X+o_p(X).

As all H(m) are odd, this is equivalently the parameter count in the unique odd class modulo 2p lying above a. The root converse involution preserves H exactly and has only O(sqrt X) fixed points on parameters counted at X. Consequently the number I_a(X) of DISTINCT root scalar values in that class satisfies

I_a(X) <= (log(2)/(2p))X+o_p(X).

Each fixed odd residue class modulo 2p therefore contains at least

((1-log(2))/(2p))X-o_p(X)

missing root scalar values at most X. No existence of an image-density limit is asserted.

## 1. A finite primitive digit chain

For p>1 odd use the states (a,q) in (Z/2pZ) x (Z/pZ). Reading binary digits from most significant to least significant, appending digit b in {0,1} sends

(a,q) -> (2a+b mod 2p, q+b(a+(a mod 2)) mod p).

This is exactly the recurrence for (m mod 2p,Q(m) mod p). Give both digits probability 1/2. Every state has two incoming labelled edges: the destination a fixes b by parity, the equation 2v=a-b modulo 2p has exactly two solutions v, and for each v the old q is uniquely determined. Hence the uniform distribution on all 2p^2 states is stationary.

The chain is strongly connected. First, from any state, append sufficiently many digits representing r such that 2^k a+r=0 modulo 2p; choosing 2^k>=2p guarantees such an r in [0,2^k). This reaches a state with first coordinate zero. Next choose b>=1 with 2^b=1 modulo p, and append the binary block for r=2(2^b-1), of length b+1. This preserves the zero first coordinate and changes q by -1 modulo p. Indeed the block identity

Q(2^k v+r)=Q(v)+Q(r)+v rev_k(r)+topbit_k(r)(v mod 2)

has vanishing cross terms modulo p when v=0 modulo 2p, while Q(r)=Q(2^b-1)=2^b-2=-1 modulo p. Thus these loops reach every q with first coordinate zero. Finally append a fixed-length binary representation of any desired first coordinate a_target in [0,2p). The final q equals the initial q plus a fixed constant modulo p, so choosing the prior zero-coordinate q reaches any desired target state. This proves strong connectivity.

There is a self-loop labelled 0 at each (0,q). Therefore the chain is irreducible and aperiodic. The elementary finite-state Markov convergence theorem implies convergence to the uniform stationary law from every initial state.

For p=1 all claims concerning residues are immediate and no chain is needed.

## 2. Uniformity after fixing endpoint digits

Fix L>=1. In a sufficiently long binary word representing m, prescribe its first L digits (including its leading 1) and its last L digits, and let the middle digits vary freely. After the first block the chain starts in a fixed state. As the middle length tends to infinity, its terminal law tends to uniform on the 2p^2 states. Appending the prescribed final block keeps q uniform conditional on the resulting first coordinate: its q update is q plus a function of the old first coordinate. In particular H(m)=1+2a-q modulo p is uniform.

Thus for each such endpoint cylinder C,

#{m in C:H(m)=a mod p}=|C|/p+o_p,L(|C|).

For fixed L this convergence is uniform over all endpoint choices, since there are finitely many of them. The first block need not have any prescribed residue: convergence holds from every state.

## 3. Passing from cylinders to the scalar cutoff

The exact digit identity proved in the parameter-capacity packet gives

Q(m)=m u(m)+O(log(2+m)),

where u(m)=sum_(j>=0) bit_j(m)/2^(j+1). If u_L retains its first L terms, then 0<=u-u_L<2^(-L), so uniformly for m<=2X,

H(m)=m(2-u_L(m))+O(X 2^(-L)+log(2+X)).

The inverse bound H(m)>=m+1 ensures that only m<=X matters; its intersecting dyadic blocks lie in m<=2X. Within an endpoint cylinder in a dyadic block, m varies by at most 2X 2^(-L), while u_L is constant. Thus the variation of H across a cylinder is O(X 2^(-L)+log(2+X)), with an absolute constant. Every cylinder that contains parameters on both sides of H=X is therefore contained in

|m(2-u_L(m))-X| <= C(X 2^(-L)+log(2+X))

for an absolute constant C.

The number of parameters in this boundary region is O(X 2^(-L)+log(2+X)+2^L). To see this without any independence assumption, fix each of the 2^L low-bit residues in turn. Its u_L and hence slope 2-u_L>=1 are fixed. The allowed m lie in an interval of length O(X 2^(-L)+log(2+X)); counting integers in that residue class divides this length by 2^L and adds O(1). Sum over the residues.

For completeness discard m<epsilon X, at a cost at most epsilon X+1 in either count. The remaining parameters lie in finitely many dyadic blocks, with their number bounded in terms of epsilon alone. For fixed epsilon and L, as X tends to infinity all these blocks have middle lengths tending to infinity. Apply Section 2 to all their full endpoint cylinders that lie entirely inside H<=X. Their sum has residue discrepancy o_p,L,epsilon(X). The partially cut cylinders contribute at most the boundary-region bound; cylinders cut by the artificial lower cutoff epsilon X contribute O(X 2^(-L)) parameters as well. It follows that

limsup_(X->infinity) |A_a(X)-A(X)/p|/X <= C' 2^(-L)+2epsilon,

where A(X)=#{m>=0:H(m)<=X}. First let epsilon tend to zero and then L tend to infinity. This proves A_a(X)=A(X)/p+o_p(X). The established exact parameter-capacity theorem A(X)=log(2)X+O(sqrt X) yields the claimed coefficient.

This proof only concerns each fixed p. It provides no uniform conclusion for moduli growing with X.

## 4. A bank excluded beyond the unweighted capacity test

The physical JJ source s=(5,6,2,3) has

H_s(m)=5+6m-3Q(m)-2(m mod 2)
      =3H(m)+2-2(m mod 2).

Its outputs are 2 or 0 modulo 3, never 1. Therefore the union of the root image and this source's tail image still misses at least

((1-log(2))/6)X-o(X)

odd integers at most X in the class 1 modulo 3 (equivalently 1 modulo 6). This proves that the root-plus-JJ bank is not cofinite.

The previously available unweighted bound for this bank is (log(2))/2+(log(2))/3=(5/6)log(2)>1/2, and so does not exclude it. The strict inequality follows from log(2)>3/5, for instance by the midpoint lower bound for the convex function 1/x on the two intervals [1,3/2] and [3/2,2]: log(2)> (1/2)(4/5+4/7)=24/35>3/5. The residue argument supplies genuinely additional information.

More generally, any bank whose non-root members all omit a common odd residue class modulo a fixed 2p is noncofinite. Root O/Z descendants may be substituted for the root without weakening that conclusion. This criterion does not apply when the other members collectively fill every modular class; it is not a general finite-bank obstruction.
