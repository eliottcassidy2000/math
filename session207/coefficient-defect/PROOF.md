# Finite-defect polynomial universality for mixed-prime multiplicities

Status: independently audited PROVED consequence of inherited finite strong inventories and the fixed-count SCC generating formula. This is a separate extension of the mixed-prime release.

For each fixed j>=0, there is a rational polynomial Q_j(b,c) of total degree at most2j such that whenever b,c>=0 and k=b+c>=max(1,j), the coefficient of n^(k-j) in the eventual ordinary-class polynomial N_n(3^b5^c) equals Q_j(b,c)/(b!c!). The polynomial depends only on the finite complete strong inventories with H=3^u5^v and1<=u+v<=j+1. In particular a fixed finite computation determines the first j+1 coefficients for the whole infinite family simultaneously.

Proof. Distinguish every strong isomorphism class A of composite prime weight w_A=u_A+v_A>=2. A core with m_A copies of such atoms has deficit d=sum_A(w_A-1)m_A. It has ℓ=k-d components, consumes U=sum u_A m_A and V=sum v_A m_A prime factors in the composite atoms, and has total size

    s=3b+4c+sum_A(|A|-3u_A-4v_A)m_A.

The remaining components are b-U copies of C3 and c-V copies of S4. The exact contribution of these words is

    (n-s+ℓ)_(ℓ) / [(b-U)!(c-V)! product_A m_A!],

where (x)_(ℓ) denotes the falling factorial. The coefficient n^(k-j) is zero unless d<=j. That condition implies every used w_A<=j+1, at mostj composite atoms, and U+V=d+sum m_A<=2d. Hence only finitely many composite profiles independent of b,c matter.

Multiplying the contribution by b!c! gives the prime falling factorials (b)_U(c)_V, of total degree U+V<=2d, divided by constant product m_A!, times the coefficient at deficit r=j-d of (n-s+ℓ)_(ℓ). That last coefficient is (-1)^r times the rth elementary symmetric polynomial of s-1,s-2,...,s-ℓ. By Newton identities and Faulhaber power sums it is polynomial in s and ℓ of total degree at most2r. Since s,ℓ are affine in b,c, each profile contributes a rational polynomial of degree at most2d+2r=2j.

Sum over the finite profiles with d<=j. Profiles requiring U>b or V>c contribute zero through the prime falling factorial, so the same polynomial formula is valid on every stated nonnegative b,c boundary. Since k>=j, ℓ=k-d>=0 and r=j-d<=ℓ; no negative factorial or negative-order coefficient is needed. This proves the statement.

Q_0=1 and the audited mixed-prime theorem gives

    Q_1=-(b²+4bc+c²+3b+5c)/2.

The finite-input claim uses Busch to bound the orders of the specified finite set of atom counts; it does not presume global spectrum completeness. No claim that Q_j has degree exactly2j or that every coefficient is nonzero is made.

## Provenance and verification

The fixed-count eventual polynomial and SCC-word formula are inherited from `05-knowledge/results/spectrum_basins_20261020e/spectrum.md` section4. The prime atoms C3 and S4 and the explicitly audited Q1 formula are recorded in `session207/mixed-prime-multiplicities/`. Busch supplies finiteness of each fixed-count strong inventory. No novelty is claimed for those ingredients. The additional result is the finite-defect locality and universal polynomial degree bound for the entire infinite mixed-prime family.

`check.py` supplies standard-library exact rational controls. `independent_symbolic_audit.py` supplies separate optional symbolic controls requiring SymPy; it is not a theorem dependency. Both check the pure-three specialization Q2(b,0)=b(b-1)(3b²+17b+2)/24 using the inherited inventories P9=2z5 and P27=2z6+3z7. The known N9 constant2 and N27 linear coefficient10/3 are boundary controls. These finite tests supplement the all-order proof and do not establish it by interpolation.
