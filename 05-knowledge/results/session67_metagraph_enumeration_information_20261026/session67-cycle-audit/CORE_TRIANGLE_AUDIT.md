# Independent audit: core-size triangle and uniform large-core asymptotic

**Verdict: PASS. No mathematical correction required.** The core-deficit tail,
explicit Burnside support bound, ordinary-class full-core limit, and uniform
large-core column asymptotic in `session67-metagraph/REPORT.md` are valid with
the qualifications stated there. The first audit `REPORT.md` and its original
verifier/receipts are preserved unchanged.

Audited files and SHA-256:

- Peer report: `41aa33c42378135566e15c9f3025146470141ff0f0fd3179cc9f8d03ad7138ef`
- Peer verifier: `731d4322c537ff57ed4361fcddfd50ae51014fe4d65f636616e993c261d4b805`
- Peer receipt: `b34005c58804886147298ccb2e152b1522f82ba8ee9eb5e7e8c5732f4a5d8c5e`

All seven source snapshots were independently compared byte-for-byte with
commit `b4234f728c0516b5fc521c001d69e26f3d2508d5`, and their SHA-256 and Git
blob identities were checked. This is not merely agreement between two
unverified local copies.

## 1. Exact refinement and deficit tail

The formula

    A[n,q] = sum_(lambda odd partition q) c_lambda/z_lambda
                  [x^(n-q)] product_(d in lambda) (1-x^d)^(-1)

correctly records positive chain weights fixed by a core automorphism.
The marker is `u^d`, not `u`, on a quotient d-cycle. Thus its claimed bivariate
form is also correct:

    A(x,u)=Z_C(u^d x^d/(1-x^d))
          =Z_T(u^d x^d/(1-x^d+u^d x^d)).

The second implementation expands
`(1-(1-u^d)x^d)^(-m)` directly; its indices and binomial signs agree with the
formula, rather than accidentally marking the number of automorphism cycles.

For `1<=d<n`, the inequality

    sum_(q<=n-d) A[n,q] <= binom(n-1,d) a_(n-d)

has a genuine surjective proof. From a child with q chains and at least d
surplus vertices, delete d vertices while leaving each chain nonempty.
Reinflate each of the n-d surviving vertices by a positive chain. There are
`binom(n-1,d)` vectors on a fixed representative of the parent. A choice of
isomorphism to that representative transports any reconstruction vector, so
using unlabelled parent representatives creates no omission. Different vectors
may yield isomorphic children, which is harmless for this upper bound.
Conversely, such inflation cannot increase the canonical core order above
n-d. At d=n-1 equality gives the unique transitive class, a useful endpoint
control. The statement is an inequality, not duplicate-free enumeration.

## 2. Burnside support bound

Write `L_n=2^binom(n,2)/n!` and `a_n=L_n(1+epsilon_n)`. Burnside gives

    epsilon_n = sum_(nonidentity permutations g) Fix_T(g)/2^binom(n,2).

Permutations containing an even cycle have zero fixed tournaments and can be
omitted. For an all-odd permutation with m moving vertices:

1. The affected unordered pairs are exactly those having a moving endpoint,
   numbering `m(2n-m-1)/2`.
2. No affected pair is fixed as an unordered pair: that would require either
   two fixed endpoints or a transposition.
3. Every affected pair orbit has odd length and therefore length at least3.
4. Its binary edge-orientation freedom consequently loses at least
   `ceil(m(2n-m-1)/3)` bits relative to the identity.
5. At most `(n)_m=n!/(n-m)!` permutations have this support size. This bound
   deliberately counts many forbidden/extra permutations; it remains valid.

Hence the reported exact finite bound

    0<=epsilon_n<=E_n
      :=sum_(m=3..n) (n)_m 2^(-ceil(m(2n-m-1)/3))

is proved. Since `2n-m-1>=n-1` and `(n)_m<=n^m`, setting
`rho_n=n*2^(-(n-1)/3)` gives

    E_n <= rho_n^3/(1-rho_n),  when rho_n<1.

This step is not justified merely by checking small n; the termwise bound
and geometric summation prove it for every permitted n. In particular
`epsilon_n->0`. A convenient rational majorant is

    E_n <= 2 n^3 / 2^(n-1),  n>=16,

because rho_16=1/2 and rho_n decreases thereafter.

## 3. Ordinary full-core mass and uniformity

The deficit inequality at d=1 gives

    0 <= 1-r_n/a_n
      <= n(n-1)2^(-(n-1)) (1+epsilon_(n-1)),

where `r_n=A[n,n]`. The denominator `1+epsilon_n>=1` was only discarded
upward. Therefore `r_n/a_n->1` for ordinary classes, not just labelled mass.
The new proof correctly supplies what the labelled union bound alone lacks.

Let N_q count nonrigid tournament classes, and R_q count rigid pair-free
classes. The orbit-mass identity is

    a_q-L_q = sum_[T] (1-1/|Aut(T)|).

Every nonrigid tournament has odd automorphism order at least3, so

    N_q <= (3/2)(a_q-L_q),
    R_q >= a_q-(q-1)a_(q-1)-(3/2)(a_q-L_q).

Subtracting both exceptional populations is valid even though they overlap.
For a rigid q-core, all positive weight vectors of sum n are inequivalent.
For any q-core there are at most this many orbits. Therefore

    R_q binom(n-1,q-1) <= A[n,q] <= a_q binom(n-1,q-1),  n>=q>=1.

Dividing gives a squeeze whose two endpoints depend on q alone. This proves
**uniformity over every n>=q**, including arbitrarily fast-growing surplus
n-q. No error estimate that depends on n is being silently exchanged with a
limit. Since `R_q/a_q->1`, both asymptotic equivalences in the peer report hold.

An explicit corollary useful for the final statement is, for `n>=q>=2`,

    0 <= 1 - A[n,q]/(a_q binom(n-1,q-1))
      <= q(q-1)2^(-(q-1))(1+E_(q-1)) + (3/2)E_q.

It implies uniform relative error `O(q^3*2^(-q))` from the coarse support
estimate alone. This is an asymptotic of counts, not a bound on structural
classification complexity or on metagraph-edge transport.

## 4. Hostile boundaries

- **q must grow.** For fixed q=3,

      A[n,3]=(binom(n-1,2)+2*[3 divides n])/3,
      a_3=2,

  so `A[n,3]/(a_3 binom(n-1,2))->1/6`, not1. The claimed theorem correctly
  excludes this limit. At q=2 the entire column vanishes.
- **n=q is allowed.** The claim then reduces to the proved `r_q/a_q->1`;
  it is not a reason to exclude the diagonal.
- **Individual cores need not be rigid.** Their periodic weight corrections
  remain; the squeeze works by bounding the exceptional core population.
- **A finite bank still misses almost all classes.** Its polynomial weight
  growth does not contradict the uniform law at growing q.
- **Labelled and ordinary measures remain distinct.** Neither the original
  labelled pair-free union bound nor generic labelled asymmetry is substituted
  for the new ordinary Burnside estimate.

## 5. Independent checks and replay

The peer command was replayed to an owned file:

    python3 ../session67-metagraph/verify_core_families.py \
      --max-order 40 --class-order 7 --out metagraph-replayed.json
    cmp metagraph-replayed.json ../session67-metagraph/receipt.json

It passed and reproduced the peer receipt byte-for-byte. Separately, the new
`audit_core_triangle.py` uses the original auditor's literal orbit routines,
not the peer's core/automorphism implementation. Its normal and `-O` runs
are byte-identical and execute **10,664 checks**:

- all8,656 nonidentity odd cycle types through order40 obey the claimed
  support loss;
- exact epsilon_n is between0 and E_n for n=0..40;
- all780 nontrivial deficit tails and820 uniform relative-error inequalities
  are checked against the exact triangle;
- every fixed-q3 coefficient through40 agrees with the hostile formula;
- independent exhaustive labelled-orbit classification through order6 gives
  actual rigid pair-free counts `1,0,0,0,2,18`, and verifies the nonrigid
  mass bound and225 rigid-core lower bounds across the triangle;
- every source snapshot's SHA-256, Git blob and commit bytes agree;
- the peer normal and optimized receipts agree.

Reproduce from this directory:

    python3 audit_core_triangle.py > CORE_TRIANGLE_RECEIPT.json
    python3 -O audit_core_triangle.py > CORE_TRIANGLE_RECEIPT.optimized.json
    cmp CORE_TRIANGLE_RECEIPT.json CORE_TRIANGLE_RECEIPT.optimized.json

Hashes:

- `audit_core_triangle.py`:
  `676b704e028f97a05fb1641a4550564206b9cffa37af783a2b3cfc84be52550d`
- Both `CORE_TRIANGLE_RECEIPT` JSON files:
  `1ebd7d541a69c67fd74772ed4add78aefa89f9af929f5a9fc41549adfed8a757`

The asymptotic conclusions are supported by the all-order proofs above;
finite checks audit formulas and boundaries but do not prove the limits.
No priority claim was audited or established for these elementary corollaries.
No shared repository files were edited or pushed.
