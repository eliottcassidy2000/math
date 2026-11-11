# The central square as a reversible core and a prime-power information scanner

Date: 2026-11-11. Status: **PROVED** elementary theorems below; **FINITE-EXACT**
small controls in `checks.json`. This packet does not claim a new proof of
Collatz, a fast factoring algorithm, or literature priority for elementary
ring/channel facts. Its contribution is their exact assembly around the user's
central-table picture.

## Inheritance, convention, and concept board

Closest mechanism: THM-4792,
`01-canon/theorems/THM-4792-two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification.md`,
part (a), identifies the centre with the coherent rational pair `+-1/4`.
The explicit integer representatives come from
`05-knowledge/results/tiling_modular_atoms_20261003.md`, section 3.
The least-used sidecar is the address/phase retained in
`05-knowledge/results/modular_fourier_zeckendorf_routes_20261003.md`.
The corrected near miss is confusing a perfect centre, which exists for every
odd modulus, with a globally reversible table. The hostile pair is `n=9,11`:
both have the same central pattern but only one has every nonzero row invertible.

Portfolio: anchor = central law and prime/composite split; niche = optimal
information sidecars; wildcard = finite monoids compiling infinite multiplier
words. Live board: central complement, kernel size, divisor strata, prime-power
digits, coherent inverse limits. All five meet in the theorems below; the final
section states where the Collatz bridge needs more information.

Use `n>=3` odd, rows and columns `1,...,n-1`, least nonnegative cell residues,
and `k=(n-1)/2`. Four central cells are

    C_n = [[d,c],[c,d]],    d=4^-1 mod n,    c=n-d.

The user's two displayed formulas are precisely:

    n=4N+1: C_n=[[3N+1,N],[N,3N+1]],
    n=4N-1: C_n=[[N,3N-1],[3N-1,N]].

A four-cell square has only four values. There is nevertheless a natural
six-number augmented display preserving the two copies of `n` as well:

    [d  c | n]
    [c  d | n].

Here the last column records the two INTEGER row sums, not extra modular
products. Equivalently the two column sums are `n`. Thus the literal six-number
statement has an exact realization without deleting its `n,n` terms. This is
a proposed interpretation, not a claim to know which display the user had in
mind. If zero is represented by `n`, other augmented drawings are possible.
The standard nonzero-index four-cell block itself contains no zero residue.

The four-cell law is inherited, not newly claimed. Its mechanism is
`2k=-1`, `2(k+1)=1`, so products are `+-4^-1`; complementary positive
representatives sum to `n`.

## CT1. A sharp central radius detects the smallest prime factor

**PROVED.** Let `p` be the smallest prime divisor of `n`. Write the symmetric
row pair at radius `r`, `0<=r<=k-1`, as

    a_r=k-r,    b_r=k+1+r=n-a_r.

Then

    gcd(a_r,n)=gcd(b_r,n)=gcd(2r+1,n).                 (1)

Consequently the largest contiguous symmetric interval of invertible rows
containing the central two rows has exactly `p-1` rows. Its first singular
pair occurs at `r=(p-1)/2` when `n` is composite; for prime `n`, that pair would
have labels `0,n`, just outside the nonzero table. The whole nonzero table
has reversible rows if and only if `n` is prime.

*Proof.* Multiplication by 2 does not change a gcd with odd `n`, and
`2a_r=n-(2r+1)`, `2b_r=n+(2r+1)`. Every positive odd number less than `p` is
coprime to `n`; `p` is not. The admissible invertible radii before this first
failure are `0,...,(p-3)/2`, hence twice `(p-1)/2` rows. Multiplication by a
row label permutes `Z/n` exactly when that label is a unit. QED.

**Interpretation.** The centre is a universally lossless seed. Prime tables
are exactly those where this seed expands without encountering a defective
row before reaching the boundary. Composite tables contain the same seed,
but their first information defect has a precise geometric distance.
The central `(p-1)`-by-`(p-1)` subarray has unit labels in both directions;
this does NOT say each restricted row permutes the restricted interval.
The reversibility assertion is about its full row channel on `Z/n`.

The central pair itself sees no prime/composite distinction. It is the
*extension radius*, not a changed central formula, that carries the split.
For `n=9`, the invertible central interval is `{4,5}`; for `n=25` it is
`{11,12,13,14}`; for `n=15` it is `{7,8}`.

## CT2. Exact channel loss and the smallest sufficient sidecar

**PROVED.** Fix a row `a` of the table and let `g=gcd(a,n)`, `m=n/g`.
For `X` uniform on `Z/n`, the observed product `Y=aX mod n` has

    H(Y)=I(X;Y)=log_2(n/g),    H(X|Y)=log_2 g.        (2)

There are `m` outputs, each with exactly `g` inputs. The product plus a
sidecar taking precisely `g` possible values reconstructs `X`; no smaller
sidecar alphabet works for every input. A canonical decoder is

    u=(a/g)^-1 (Y/g) mod m,    X=u+t*m,  0<=t<g.     (3)

Here `Y` is its representative in `{0,...,n-1}`, necessarily divisible by `g`.
When `m=1` (zero row) interpret `u=0` and use the entire input as `t`.

*Proof.* `aX=aX' mod n` iff `X=X' mod m`. This gives the fibers and uniform
output. Each fiber needs `g` distinguishable tags; (3) provides exactly that.
QED.

This makes "maximal efficiency in information transfer" precise: among
moduli `n>=3`, zero loss for EVERY nonzero multiplier is equivalent to
primality. Individual composite-modulus unit rows are equally efficient.
This is a conditional, known-multiplier communication model; an observer
who does not know the row label loses additional information even over a field.

A tempting smaller-format sidecar `X mod g` works iff `gcd(g,m)=1`, not in
general. Indeed the pair `(X mod m, X mod g)` identifies `X mod lcm(m,g)`;
it is injective on `Z/n` iff `lcm(m,g)=n`. For `n=9,a=3`, inputs `0,3,6`
have both the same product and the same value modulo `g=3`. The missing
coordinate is a HIGH base-3 digit, not another copy of the low digit.

## CT3. CRT turns lost information into erased prime-power digits

**PROVED.** Factor `n=product p^e`. Set `v_p=min(e,ord_p(a))`, assigning
`v_p=e` if the local multiplier is zero. In the CRT component `Z/p^e`,
multiplication by `a` reveals `X mod p^(e-v_p)` after removing its known unit
part, and loses exactly the upper `v_p` base-p digits. Hence

    H(X|aX)=sum_(p^e || n) v_p log_2 p.               (4)

CRT reconstruction of the componentwise digit sidecars is an optimal
reconstruction of the whole input. This is recursive at prime powers: one
additional factor of `p` in the multiplier erases one more digit, until all
`e` digits are lost. Coprime factors run independently for uniform `X`.

*Proof.* In the p-component write `a=p^v u`, `u` a unit when `v<e`.
Dividing the observed product by `p^v` reveals `uX mod p^(e-v)`; invert `u`.
Every choice of the remaining `v` digits is possible. CRT makes a uniform
input into independent uniform components, and the fiber sizes multiply.
This proves (4) and the explicit sidecar construction. QED.

There is a short exact average law as well. If the known multiplier `A` is
uniform on the nonzero labels `1,...,n-1`, then

    E H(X|AX,A=a)
      = sum_(p^e || n) log_2(p) * sum_(j=1)^e (n/p^j-1)/(n-1).   (5)

The notation denotes the average of the conditional loss for each `a`.
It is zero iff `n` is prime, positive for every composite `n`. If `A` is
instead uniform on all of `Z/n`, the loss is the additive CRT quantity

    sum_(p^e || n) [(1-p^-e)/(p-1)] log_2 p.          (6)

*Proof.* For `v=min(e,ord_p(A))`, `E v=sum_(j=1)^e P(p^j|A)`.
The nonzero sampling probabilities are `(n/p^j-1)/(n-1)` and the full-ring
probabilities are `p^-j`. Substitute in (4). QED.

Equations (5)-(6) compare precisely specified channels, rather than asserting
that prime tables have lower raw description length than composite tables.

## CT4. Minimal defect alphabet and an exact finite compiler

**PROVED.** The action of the unit group `(Z/n)^*` on the residue labels has
one orbit for each divisor `d|n`, namely

    S_d={x:gcd(x,n)=d},    |S_d|=phi(n/d).            (7)

The zero orbit is `d=n`. Thus the nonzero labels have exactly `tau(n)-1`
unit-scaling types. This number is as small as possible, namely one, iff
`n` is prime. It is a precise "prime set as small as possible" statement:
the *defect alphabet* is smallest, not the number of distinct residues.

*Proof.* Unit scaling preserves the gcd. Locally in `Z/p^e`, two elements
of the same valuation are unit multiples: cancel their common `p^v`, then
lift the resulting unit ratio from `Z/p^(e-v)` to `Z/p^e`. Choose these units
independently and combine by CRT. Counting the orbit gives `phi(n/d)`.
QED.

The map from multiplication to its defect alphabet is the finite monoid
homomorphism

    q(a)=(min(e_p,ord_p(a)))_p,
    q(ab)_p=min(e_p,q(a)_p+q(b)_p).                 (8)

Its state space is the grid `product_p {0,...,e_p}`, with `tau(n)` states.
For an arbitrary finite or infinite multiplier word `a_1,a_2,...`, one update
per letter gives the exact kernel size, output cardinality, and information
loss of every prefix product. An infinite input is handled by one finite
machine for each fixed modulus; no enumeration of its infinitely many
prefix products is needed to define the monitor.

This quotient is MINIMAL among deterministic summaries from which one can
recover the kernel size of every current product: two distinct valuation
vectors give distinct integer gcd values by unique factorization, hence
distinct kernel sizes. Any such exact summary must have at least `tau(n)`
states. The compiler realizes the lower bound.

The entire odd multiplication family also has a coherent infinite object:
the inverse limit over odd moduli with divisibility transition maps. The
values `d_n=4^-1 mod n` are compatible since reduction preserves inverses.
Thus the centre is the single element `1/4` in that inverse limit and its
negative, not a new unrelated object at each modulus. Integer representative
formulas depend on the modulus and do not themselves preserve an order in
the inverse limit.

## CT5. Two affine sheets force full residue memory

**PROVED.** Assume `gcd(n,6)=1`, `n>=5`. Let `Q:Z/n -> S` be a
deterministic quotient admitting maps `f_+,f_-:S -> S` with

    Q(3x+1)=f_+(Q(x)),    Q(3x-1)=f_-(Q(x)).

Suppose a readout of `Q(x)` determines whether `x=0 mod n`. Then `Q`
is injective: at least `n` states are necessary. The full residue machine
attains this bound. Equivalently, recognizing words in the two affine sheets
that send the initial residue to zero requires all reachable residue states;
all `n` states are reachable from any initial residue.

*Proof.* The maps `F_+(x)=3x+1`, `F_-(x)=3x-1` are permutations of the
finite residue set. Their inverses are positive powers of themselves, so
every composition involving their inverses is represented by a positive
word in the two sheets. In particular

    F_+ composed with F_-^-1 (x) = x+2.

Since `n` is odd, this translation generates all translations `x -> x+t`.
If `Q(x)=Q(y)`, applying any common sheet word preserves equality. Choose
a word representing translation by `-x`; the resulting states have underlying
residues `0` and `y-x`. The zero readout forces `y-x=0`, proving injectivity.
Translations also give reachability from any initial residue. QED.

This places an exact information cost on a natural attempted Collatz bridge:
the `tau(n)`-state multiplicative defect compiler cannot simply be extended
to the two signed affine sheets while still recognizing zero. Its state
budget must rise to `n`. The statement is about arbitrary sheet words, not
parity-guarded integer Collatz trajectories. It does not exclude stronger
compression when the allowed words are restricted by a genuine guard, or
when a weaker observable replaces zero recognition. The condition `3` not
dividing `n` is essential to this proof, since it uses inverse permutations.

## Connection ledger and actual boundary

1. Central table -> reversible channel: map a row label to `X -> aX`;
   preserved predicate = injectivity; lost information = kernel coset;
   sidecar = high prime-power digits; decisive test = `n=9,a=3`.
2. Table -> finite defect monoid: map `a -> q(a)`; preserved = multiplicative
   prefix kernel size and loss; destroyed = unit phase and exact product;
   sidecar = unit residues; decisive test = equal-gcd inputs with different
   products, e.g. `1,2` modulo 9.
3. Centre -> signed affine cycle: inherited `3d-1=-d`, `3(-d)+1=d`;
   preserved = this two-step orbit; lost = general affine dynamics. The
   valuation quotient is NOT an affine quotient: modulo 5, unit inputs 1
   and 3 have the same valuation state, but `3*1+1=4` is a unit and
   `3*3+1=0` is zero. Exact residue/phase is therefore indispensable for
   a Collatz-compatible extension of (8). No forward convergence claim follows.
4. Finite tables -> inverse-limit centre: reduction maps preserve the central
   pair and multiplication; destroyed = the preferred integer representative
   and archimedean size. A sidecar measuring size would be needed for descent.

The strongest next construction is not a broad speculation list: augment
(8) with just enough local unit phase to process an assigned finite set of
affine branches, then minimize the resulting deterministic machine while
retaining the exact branch guard. This is an OPEN construction target,
not a claim that a fixed finite quotient certifies unbounded Collatz descent.

## Reproduction

From the repository root:

    python3 05-knowledge/results/session255_modular_collatz_20261111/central_table/check.py

The checks use exact integer fiber cardinalities and rational probabilities,
not floating entropy comparisons. Their finite universe and positive/hostile
controls are recorded in `checks.json`. Proofs above supply the unbounded
statements; the checks are probes of the new construction, not a broad audit
of inherited research.
