# Independent Singer / line-flip audit

Status: PASS, with scope qualifications below. No Git operations performed. Reproduce: `python session35-singer-audit/check.py`.

## 1. Setup and exact scope

Let a finite projective plane have order q, v=q²+q+1 points and lines. Each unordered pair of distinct points lies on exactly one line. Fix any baseline tournament T, write S(u,w)=+1 for u→w and −1 for w→u, and define h(u,L)=sum_{w in L minus u} S(u,w) if u belongs to L, zero otherwise. Flip the entire edge block on each selected line. The blocks partition the edges, so masks x in F₂^v parameterize 2^v distinct labelled tournaments. Let z_L=(-1)^{x_L}. The score is

 d_x(u)=(v−1)/2 + (1/2) sum_L h(u,L)z_L.

This formula holds without requiring a regular baseline. Regularity means sum_L h(u,L)=0 for every u at mask zero. All statements here concern labelled states and labelled line-flip continuations, not a quotient by arbitrary tournament isomorphism.

## 2. Odd-q parity and regular members

For odd q each nonzero local h(u,L) is an odd integer, being a sum of q signs. Reversing line L changes d(u) by −h(u,L), so modulo two the score difference from the baseline is Bx, where B is point-by-line incidence over F₂.

Each line has q+1 even points, two different lines meet once, and v is odd. Thus BᵀB=I+J and B1=0. The kernel of I+J is precisely span(1): (I+J)y=0 implies y=(sum y_i)1, and 1 satisfies this since v is odd. Therefore ker B is contained in span(1), and B1=0 proves equality. In particular rank B=v−1. This argument uses only plane axioms, not Desarguesian coordinates.

If the baseline is regular, any regular family member has zero score difference, hence lies in ker B, so is mask zero or all ones. Both are indeed regular: they are T and its global converse. This does NOT claim there are only two regular tournaments on v vertices, only two in this particular line-flip family.

## 3. Triangle Walsh expansion and positive constant

Let c3 count cyclic unordered triples. For every tournament,

 c3(T_x)=binom(v,3)−sum_u binom(d_x(u),2).

This follows because every noncyclic triple has a unique vertex with two outgoing arcs. Since sum_u(d_x(u)−(v−1)/2)=0, expanding gives a constant plus pair characters z_L z_M, with coefficient

 a_{LM}=−h(w,L)h(w,M)/4, where {w}=L intersection M.

There are no degree-one characters. The constant equals

 A=binom(v,3)−v*binom((v−1)/2,2)−(1/8)sum_{u,L}h(u,L)².

Every a_{LM} is nonzero for odd q. This assertion needs no regular-baseline hypothesis.

The sometimes necessary abstract condition A≠0 is automatic for projective planes. A uniformly random independent line mask makes the three edges of any noncollinear triple independent fair orientations, so its expected cyclic indicator is 1/4. A collinear triple has all its edges simultaneously reversed and keeps its baseline cyclic indicator. Hence

 A=[binom(v,3)−v*binom(q+1,3)]/4 + (number of cyclic collinear baseline triples)>0.

The strict inequality holds for every projective plane of order q≥2, which has noncollinear triples.

## 4. Continuation quotient and translate-linear rank

The alphabet consists of all named line flips, with arbitrary finite words permitted. Their net masks are arbitrary t in F₂^v. Two masks x,y are continuation-equivalent for c3 iff c3(x+t)=c3(y+t) for every t. Fourier uniqueness says this holds iff every nonzero Fourier character takes the same value at x and y. Full pair support says x_L+y_L=x_M+y_M for every pair L,M. Thus x+y is zero or all ones. Conversely global reversal preserves c3, so these are exactly the equivalence classes.

Every mask is reachable. The minimal deterministic output automaton therefore has 2^(v−1) states, by the ordinary distinguishable-futures argument, and the quotient by the global-converse mask realizes the bound.

Over R, Q, or any appropriate characteristic-zero field, the span of all translates of c3 has dimension equal to its number of nonzero Walsh coefficients: 1+binom(v,2). One proof projects each character by the group-average idempotent; the nonzero coefficient recovers that character from translates, while the Fourier expansion proves the opposite inclusion. The same number is the translation Hankel rank and the minimal unconstrained linear realization dimension, if that is the meaning of “linear memory.” It is not a deterministic state count. For q=7 these numbers are 2^56=72057594037927936 and 1597, respectively.

## 5. Singer57, CRT, and symmetry obstruction

The supplied D={8,15,19,38,48,50,51,56} in Z/57 is independently verified to have the 56 ordered differences exactly 1,...,56, each once. Its 57 translates are a (57,8,1) symmetric design and projective plane: every pair of points occurs once, every pair of lines meets once, and the nondegenerate plane axioms follow from these parameters. This certificate independently verifies the claimed cyclic design; identifying this particular model with PG(2,7) via a Singer trace construction additionally uses the producer's stated finite-field construction, not merely our difference check.

CRT sends x to (x mod3,x mod19), with inverse (a,b)↦19a+39b mod57. It transports additive translations to componentwise additive translations. It does not turn incidence, tournaments, or line flips into squaring dynamics. A separate ring-square fibre model carries no such extra structure without a separately proved map.

A tournament invariant under the full projective linear group on PG(2,q) cannot exist. Given two distinct projective points, choose independent representatives u,v and extend them to a basis u,v,w. The invertible linear map exchanging u,v and fixing w induces a projectivity swapping the pair; invariance would require both opposite orientations. The conclusion also holds for PSL: the map u↦v, v↦u, w↦−w has determinant one. It holds for any larger group containing one of these. It does not rule out Singer-cyclic invariant tournaments: v is odd and a cyclic regular tournament is an explicit example.

## 6. Independent finite controls

The checker builds PG(2,2) and PG(2,3) directly from normalized homogeneous coordinates and dot-product incidence. It does not import the producer's code or use extension fields. It enumerates every mask and constructs its actual edge orientations, counts actual cyclic triples, independently checks the score and Fourier formulas, and computes all Walsh support/stabilizers.

- q=2 cyclic-line baseline: all 128 masks regular, c3 always 14, incidence rank 4, one Walsh term, one continuation class. Odd-q claims therefore cannot be silently extended to even q.
- q=3: all 8192 masks checked; exactly two regular, rank 12, 79 nonzero Walsh terms, 4096 continuation classes. Its c3 range is [55,91], with mean 306/4. The producer's Singer-labelled circulant baseline has a different range/mean because its baseline gauge differs; this is not an invariant of the underlying plane.
- q=7: exact difference-set and rank checks only; no exhaustive mask census was attempted or claimed. General proofs establish the full family assertions.

## 7. Nonclaims and audit qualifications

Here c3 is NOT a Hamiltonian-path/cycle count H. Neither the support/rank nor continuation conclusions transfer to H without independent evidence. The regular-baseline premise is essential for identifying the two regular states as the starting tournament and its converse. Pair Walsh nonvanishing and c3 continuation results hold for every baseline when q is odd. The observed range in a sampled q=7 receipt is only a sampled range, not a global extremum. No unresolved mathematical blocker remains for the stated claims.
