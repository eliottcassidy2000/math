# Independent audit of the explicit (4,6) family and exceptional-point ledger

**Status: PASS / complete analytic read and exact replay.** Root audit,
September 6, 2026. The certificate agent also independently read the full
classification, the cusp-separation argument and the exceptional-point
Euler deduction and reported PASS. The accepted scope is the displayed
family as the whole nonproperness support, with three cusp parameters
remaining OPEN. There is no claim of general JC(2) progress by chart entry.

## Algebraic and local audit

1. The pair-sum substitution is invertible in the collision quotient:
   `p=s+t` is a unit from `p(p^2-2q)=-1`, and `t=p-s` is the inverse.
   The monic degree-six H and monic quadratic in s make its length exactly
   twelve, including every exceptional parameter. This retains the ordered
   pairs, diagonal points and their scheme lengths, not only a resultant.
2. Finiteness follows from the monic equation for t. A generic partner
   would give a positive-dimensional off-diagonal collision locus, contrary
   to the finite algebra. Thus the parametrization is birational and is
   the normalization; the literal pole degrees are 4 and6.
3. The three critical values are distinct. The second/third derivative
   determinant is nonzero at each critical source point, giving an ordinary
   cusp. At each critical parameter the pair-sum polynomial is squarefree;
   exactly one quadratic factor has a double root. The diagonal length is
   two and the other five unordered collisions are reduced.
4. The actual triple-fibre calculation is a coefficient comparison, not
   a projection-multiplicity inference. It forces the displayed quartic
   in A and lambda=-A^2, and its inverse A=(lambda^2+1)/2 makes the triple
   target unique. The cubic discriminant and fourth-root gcds exclude
   repeated preimages and an additional fourth preimage.
5. A cusp cannot share its image with another branch. Both target fibre
   polynomials would contain `(t-r)^2(t-u)`. Their cubic remainder would
   therefore divide the quartic with these multiplicities, forcing the
   triple coefficient relations. The proved squarefree cubic at those
   relations contradicts the double root. This checks multiplicities,
   not just distinct root counts.
6. At every tangency parameter the nonzero linear subresultant has
   invertible leading coefficient. The specialized gcd of H and H' has
   degree exactly one: precisely one double root and no higher root.
   Since diagonal and triple loci are disjoint, this gives one ordinary
   tacnode and four nodes. For immersed branches, the length-two collision
   algebra equals contact order two by local graph elimination.
7. All other singularities are exhausted by critical points and multiple
   fibres of the finite normalization. A single immersed branch is smooth.
   The lambda=-1 control has a reduced collision scheme but a repeated
   projected resultant, which correctly tests the omitted partner data.

## Euler deduction and the actual boundary

The inherited nodal note was independently audited before this argument.
For one exceptional point with r local branches and N ordinary nodes,
normalization gives `chi_c(S_reg)=1-2N-r` and
`chi_c(C^2 minus S)=N+r-1`. Keep the actual fibre count n_z at the
exceptional point. Substituting the proved ordinary-node counts
`d-2delta+omega_p` gives exactly

```text
1=(r-1)delta+n_z+sum_nodes omega_p.
```

No local monodromy or specialization identity is assumed at the
exceptional point. For r>=3 the right side is at least two. For r=2 it
forces delta=1. The unique nonproperness component then has total weighted
generic deleted length one, so every deleted divisor has ramification
index one. Actual source points are already etale. The finite normal
envelope is therefore unramified in codimension one over the smooth target;
purity and simple connectedness contradict any connected degree d>1.
This last step works for all d, not only d=2.

The resulting exclusion applies to the generic nodal, triple and tangency
strata of this particular family. The remaining three cusp parameters obey
`n_cusp+sum omega_p=1`. With five nodes, the integer bound
`omega_p>=max(0,d-2a)` implies d<=2a. These are necessary conditions only.
The inherited Euler framework is credited, and no external priority claim
is attached to its present specialization.

## Exact replay and pins

The final source has **214 always-active gates**. Root ran normal and
optimized Python independently and compared both full output streams with
the frozen output. All three are byte-identical. The source uses universal
identities over Q, direct Groebner ideals at lambda=1,2,-1, and the declared
bounded integer controls; the latter corroborate the analytic inequalities
and are not substituted for their all-degree proofs.

```sh
python3 -B 04-computation/planar_jc48_sep06_curve_probe.py
python3 -B -O 04-computation/planar_jc48_sep06_curve_probe.py
```

```text
source SHA256 97c28adbf35451a8e521486680a2b79267073430fdc8fe69a96e9435ffa923a1
output SHA256 9d227fc32a01f24d13105ea4938cd0485c00ffe95225c723df4e6e293c4f786f
```

No repair of the final mathematical statements was required. The root
promoted the note only after the complete exceptional-point deduction and
the independently checked exact classification were both present.
