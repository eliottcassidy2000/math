# A precise bracket-to-Paley interface needs a marked residue ruler

Status: PROVED elementary interface and hostile example; the bracket enumeration formulas are supplied and separately checked in the arithmetic lane. No unlabelled tournament isomorphism or prime-count theorem follows.

Fix an odd prime p, let chi be its quadratic character extended by chi(0)=0, and give F_p the explicit integer representatives 0,...,p-1. The square map pushes counting mass on F_p to mass 1 at zero and 1+chi(a) at each nonzero a. Therefore its least-residue first moment is

    R_p = p(p-1)/2 + S_p,    S_p = sum_{a=1}^{p-1} a chi(a).

This is a pushforward with fibre multiplicities, not merely the support of the square map. If p=3 mod4, the row chi(a) at the distinguished origin is the signed adjacency row of the Paley tournament under the convention 0->a iff chi(a)=1. The weighted moment S_p is consequently an exact marked Paley-row readout.

In the relaxed odd-square bracket problem r=p/(p-1), the arithmetic lane's formula is

    N_p = (16p^2-27p+11)/3 + 2 S_p/p.

N_p counts integer sources n>=2 whose r-multiple remains in the same bracket, not primes, and not tournaments. At p=11, QR={1,3,4,5,9}, S_11=-11 and N_11=548. This is a real map between the two computations, but it does not explain the separate original doubling exception set {2,3,11}.

## What cannot be forgotten

The unweighted row sum is zero and loses S_p. The unlabelled Paley tournament also loses the integer cut and the distinguished origin. For fixed ruler weights a, moving the origin to v gives

    S_p(v)=sum_{a=0}^{p-1} a chi(a-v)
          = S_p(0)-p sum_{y=p-v}^{p-1} chi(y),  0<=v<p.

The equality follows by substituting a=y+v modulo p and retaining the wrap carry p. All these rooted vertices lie in the translation orbit of the same tournament, but the ruler is held fixed. For p=11 the moments are

    -11,0,-11,0,11,22,11,0,-11,0,-11.

Thus isomorphic rooted tournaments WITHOUT the transported ruler can give different moments. Transporting both origin and ruler restores covariance. This is the cheapest hostile example against a graph-only transfer; it is not a failure of the exact labelled identity.

For arbitrary original source weights mu_x, the correct moment is sum_x mu_x [x^2]_p and the square-fibre weight is sum_{x:x^2=a} mu_x. The special 1+chi(a) law uses counting mass. Replacing arbitrary weights by uniform residues is not authorized by the square-map representation.

## Composite conductors

At q=341=11*31, the unit square image is QR_11 x QR_31, identified by CRT. Each image point has four roots; the two nonunit axes have two roots per nonzero square, and zero has one. The inherited joint clock subgroup <4,81> equals this 75-point unit image, but it is NOT the entire source carrier. The arithmetic lane's exact floor correction includes the zero axes. The independent root check confirms that dropping them changes the bracket count by36.

This connects the earlier joint-clock object to a new weighted observation while keeping conductor strata, multiplicities, the chosen least-residue lift and original mass. It makes no claim that clock time equals tournament complement time.
