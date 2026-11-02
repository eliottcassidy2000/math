# What is closed, and what is not

## Established or independently rechecked

1. **PROVED, recovered exclusion:** no tournament has ordinary Hamiltonian
   path count7 or21. Cited Busch strong-minimum/overlap results reduce the
   exclusion to small orders; the selected272 overlap presentations use two
   exact counters. Independent full labelled enumeration through order6
   corroborates the small spectra. The all-size argument is the written
   proof, not an extrapolation from that enumeration.
2. **PROVED arithmetic refinement:** for the candidate monoid A of positive
   odd integers excluding7,21, the multiplicative atoms are odd primes
   p!=7, the numbers7p for odd primes p!=3, and63,343. This is a theorem
   about A, not a claim that A is already the attained spectrum.
3. **PROVED factor-length bound:** writing a=v7(h), the maximum number of
   nonunit A-factors is Omega(h)-ceil(a/2), except h=7*3^b, b>=2, where it
   is b-1; L_A(1)=0. Every actual tournament with H=h has at most L_A(h)
   non-singleton strong components. Uniform attainment is not asserted.
4. **FINITE-EXACT carry certificates:** strong tournament adjacency rows
   attain49,63,343, with independent DP and permutation counts. Combined
   with SCC multiplication these recover the exact two-prime-lane reduction.
5. **PROVED bridge obstruction:** H_n=1+2Q_n coefficientwise over the integers.
   A square integer-polynomial response map with k raw-H coordinates has
   Jacobian determinant divisible by2^k. Normalizing the count coordinate
   removes that parity factor but not the source-domain constraints.
6. **PROVED legal-domain obstruction:** weighted orientation interpolation
   fills the complete interval between minimum and maximum counts. Actual
   one-edge templates also give normalized responses2-t and9-t, whose
   derivatives are-1, yet forbidden targets3,10 have inverse t=-1 outside
   the permitted{0,1}. Even an integral proper global inverse need not be
   a legal tournament constructor.
7. **PROVED sufficient algebraic inverse, limited scope:** an integer
   polynomial map homogeneous for matching positive input/output weights
   and with Jacobian determinant+/-1 has an integer polynomial inverse,
   by block recursion in increasing weight. This supplies no Boolean or
   nonnegative inverse guarantee; the shear(x+y^2,y) is a hostile control.

## Explicit correction, preserved lineage

The main positive-weight theorem in the inherited Jacobian
`THM-1370-elliptic-jacobian-conjecture-all-dimensions.md` survives. Its broad
sentence excluding every identity-plus-ordinary-homogeneous map from positive
weighted homogeneity does not: the shear(x+y^2,y) with weights(2,1) is a
counterexample. Claims that failure requires an indefinite grading must
also be restricted to maps that admit the asserted grading. The new bridge
proof records these corrections; the canon file now contains the independently
reviewed narrow repair and explicit positive-target/origin argument. Its
historical source snapshot is unchanged.
This is separate from the other theorem sharing the number1370, whose full
slug is `THM-1370-h-spectrum-omits-7-21-all-n.md`.

## The remaining theorem

Construct a strong tournament with H=p for every odd prime p!=7 and one
with H=7p for every odd prime p!=3. Equivalently, supply a legal gadget
atlas with exact normalized response laws and a well-founded inverse descent
covering every q>=0 except3,10, terminating in certified base tournaments.
Lattice membership, rank, a continuous image, a determinant, local inverse,
finite count bank, or arithmetic factorization alone does not prove that.

Royle's reported bounded prefix is safely inclusive80555; the selected local
THM-4104 prefix is80405. These inherited prefixes are cited, not rerun here.
Later repository constructive prefixes, when used, retain their own source
certificates and do not promote the universal claim. No unrelated current
order-eleven claim, LRC(14), planar-Jacobian conjecture or earlier session's
approval is altered or promoted by this release.
