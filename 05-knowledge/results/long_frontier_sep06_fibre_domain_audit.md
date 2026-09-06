# Independent audit of the four-anchor admissible domains

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS; FINITE-EXACT replay.**
This audits [the domain theorem](long_frontier_sep06_fibre_domain.md), its
[source](../../04-computation/long_frontier_sep06_fibre_domain.py), and
[frozen output](long_frontier_sep06_fibre_domain.out). All four iff statements,
the endpoint multiplicities, and the moment-resultant bridge are accepted.
No mathematical correction was required.

## 1. Analytic audit

The four disjoint sign-changing C boxes exhaust its quartic degree. Their
F-image intervals order the four resultant roots and identify the correct
boundary rho, rather than merely choosing a positive resultant root.
Alternating values of B_f at 0 and the C roots construct all five B roots
for 0<f<rho. At f=0, B has one simple zero because F'(0)=35.

At f=rho, continuity retains five real nonnegative roots. Adjacent roots
can merge only at a separating C root; only c1 is shared. The exact bound
F'(c1)<-6 excludes a multiple B root there and identifies c1 as the second
B root. Other intervals remain disjoint. This justifies both the weak
endpoint and the asserted simplicity. Conversely, weak interlacing puts
c1 between the first two B roots, where a monic quintic is nonnegative;
therefore f<=F(c1). Necessity and sufficiency are separate and correct.

The exact D factorization and its radical F values give the analogous
first-collision threshold L=14sqrt(2)-16. The source's rational radical
bounds preserve the stated ordering. Since L<rho, the C construction gives
simple B roots even at f=L. The simultaneous domain is precisely [0,L].

For the larger B-only domain, the four derivative sign boxes exhaust the
quartic derivative. Independently recomputed rational Horner enclosures
place their F values in (43/10,22/5), (-10,-9), (26,28), and (-63,-62).
The five monotone branches construct all five roots for 0<f<kappa. At
kappa, the first critical point is a simple derivative root, so the first
two B roots merge into exactly one double root; the other three crossings
remain strict. Conversely, Rolle interlacing, including multiplicities,
gives B_f(eta1)>=0 and hence f<=kappa. The discriminant's four critical
values identify kappa uniquely in the declared interval. Thus the B-only
iff is valid on the repeated-root boundary as well as in the interior.

For degree five, `product_i B'(a_i)=Vandermonde(a)^2` has positive sign.
The moment factorization therefore gives `det H_D=Res(B,D)` with the sign
printed in the theorem. Extension across repeated roots uses a polynomial
identity, not an assertion that an uncancelled pole decomposition persists.
The determinant alone is never used as a sufficient interlacing test.

## 2. Source and exact replay

The producer's **42 explicit gates** remain enabled under Python -O.
Normal and optimized output agree byte-for-byte with the frozen output.
The source uses exact SymPy algebra for polynomial identities and
standard-library rational interval Horner bounds for the eight root boxes;
it imports no research producer and uses no floating-point roots.

The auditor separately checked the four critical sign changes and image
enclosures with a fresh Fraction implementation and recomputed the
discriminant as

```text
3125f^4+125758f^3-4825473f^2-30354226f+211485225.
```

This discriminant recomputation uses the same exact CAS as the primary
producer; it is not advertised as a second algebra engine. The independent
standard-library determinant construction in the
[response producer](../../04-computation/long_frontier_sep06_residue_tail.py)
also verifies the full D fifth-moment factorization by direct permutation
expansion. The primary source and proof text were read in full.

Reproduction from the worktree root:

```bash
python3 04-computation/long_frontier_sep06_fibre_domain.py > /tmp/domain-normal.out
python3 -O 04-computation/long_frontier_sep06_fibre_domain.py > /tmp/domain-optimized.out
cmp /tmp/domain-normal.out /tmp/domain-optimized.out
cmp /tmp/domain-normal.out 05-knowledge/results/long_frontier_sep06_fibre_domain.out
```

SHA-256 manifest checked by the auditor:

```text
source  dacdf6a3323cce90e14002a03bcbd87b105790717fdf243db754c06ef1b93f29
output  41d5f6db13170eac6926051272ad429326e0b3e0ce8fb9c70ac5ec62d828e6bf
semantic 37451b9f35fac399612c2e4b2a475c3e08e31747591e1906621e690118cc012b
```

## 3. Exact consumer and remaining scope

The B-only domain lies below 22/5<5. The separate, independently audited
[response theorem](long_frontier_sep06_residue_tail.md) is valid on the
entire coefficient interval [0,5]. Combining them therefore gives strict
negative full response at every original first-row root for every
nonnegative-root B on the four-anchor fibre, without either interlacer
assumption. The domain theorem alone does not prove this response sign.
Only f=1 is identified with the inherited actual factorial row; other
f values remain model coefficients. The general two-anchor model is open.
