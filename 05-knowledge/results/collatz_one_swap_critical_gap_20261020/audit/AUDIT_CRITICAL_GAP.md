# Independent audit: critical-shape content bound below sqrt(7)^p

Date: 20 October 2026.

**Verdict: the proposed refinement is valid.** No missing Fourier alias or
incorrect inequality was found. The proof below establishes the proposed
constant exactly, not merely numerically. The only presentational qualification
is that the critical shape at `p=2` has `q=4` and is not coprime. Thus `k=1`
follows for an existing coprime critical shape, not from the false assertion
`ceil(2 log_2 3)<4`.

**Central strengthened result:** retaining the full `k=1` legality inequalities
gives `g<7^(p/2)` at every coprime critical shape, and hence the necessary
positive-cycle gap `0<q log 2-p log 3<exp(-p/9)`. This requires only
orthogonality and AM-GM, not the variance argument. Its proof is in Section 0.

**Additional simplification:** the two large-norm endpoints are both forbidden
for an actual swap when `k=1`. Consequently the analytic interior proof and
the two free endpoints already establish the stronger statement
`g<c^p` for every coprime shape `p<q<2p`, not just the critical clock.
The proposed large-endpoint and small-weight estimates remain correct, but
are unnecessary if this admissibility observation is included.

This audit inherits only the proved content/norm and endpoint results in
`../../collatz_one_swap_sparse_norms_20261020/SHARP_ONE_SWAP_CONTENT_BOUND.md`. It does not claim all-weight exclusion,
an effective maximum weight, or an extension to more than one adjacent swap.
No frozen source file was changed.

## 0. Central result: admissibility gives the sqrt(7) bound

For an existing coprime critical shape, `p>=3` and `k=1`, as justified in
Section 1. Write `q=p+r`, `1<=r<p`. If the moved one's residue is `a`, then
the preceding residue is `A=(a+r mod p)` and the following residue is
`a'=(a-r mod p)`. Thus

    delta_before=(r+a-A)/p = 1 iff a>=p-r,
    delta_after =(r+a'-a)/p = 1 iff a<r.

At `k=1`, the minus move is allowed exactly when the preceding gap is at least
`2`, and the plus move exactly when the following gap is at least `2`.
Consequently these are the complete legality criteria, including equality
boundaries:

    plus is legal iff a<r;
    minus is legal iff a>=p-r.

At the critical clock, `r/p<=3/4`. For `p>=7`,

    r<(3/5)p+1<=(3/4)p,

using `log_2 3<8/5`; for `p=3,4,5,6`, the exact values of `r` are
`2,3,3,4`, respectively. The `p=6` shape is not coprime, but the ratio bound
still holds and its inclusion in this check is harmless.

For a legal plus-interior type, `a<r` and the integral inequality `a+1<=r`
give

    mean|Q_+|^2
      =1+4^(a/p)+4^((a+1)/p)
      <1+2*4^(r/p)
      <=1+2*4^(3/4)
      =1+4 sqrt(2)
      <7.

The first strict inequality follows because `a<r`. The final strict
comparison is equivalent to `32<36` after isolating and squaring positive
quantities. It remains valid even in the boundary shape `p=4,r=3`.

For a legal minus-interior type, `p-a<=r` and `p>=3` give

    mean|Q_-|^2
      =1+4^(1/p)+4^((p-a)/p)
      <=1+4^(1/3)+4^(3/4)
      <1+2+3
      =6<7.

The strict rational-power comparisons are certified by `4<2^3` and
`4^3=64<81=3^4`. In each case the exponents in the polynomial are distinct
modulo `p`, so the displayed mean is the exact mean over all `p` conjugates.
AM-GM now gives `|Norm(Q)|<7^(p/2)`.

The only legal endpoint types are plus `a=0` and minus `a=p-1`, whose
prime-to-6 norm parts and contents are `1`. The two large-norm endpoints
violate the exact legality criteria. Hence every actual coprime critical
one-swap word satisfies

    g<7^(p/2).

For a hypothetical positive integer cycle, `D=g>0`, so

    0<q log 2-p log 3
      =log(1+D/3^p)
      <D/3^p
      <(7/9)^(p/2)
      <exp(-p/9).

For the last step, `log(9/7)>2/9`: integrate `1/t>1/9` over `7<t<9`,
or apply `log(1+x)>x/(1+x)` at `x=2/7`. No numerical approximation is
needed. This condition alone is not an all-weight exclusion; an effective
lower bound and any ensuing finite coverage must be proved separately.

The earlier, weaker exponential-content estimate is audited below as an
independent open-band (`p<q<2p`) result, because its proof does not need the
critical restriction `r/p<=3/4`.

## 1. Precise statement and scope

Let `p>=2`, `gcd(p,q)=1`, and let `w` be obtained from the zero-intercept
mechanical word of shape `(p,q)` by one allowed cyclic unequal-adjacent swap.
Let `g` be its carry content. Assume the critical clock

    q=ceil(p log_2 3).

Then, with `c=3 exp(-1/2500)`,

    g<c^p.

For `p>=3`, `log_2 3<8/5` implies `q<8p/5+1<2p` when `p>=3`.
Also `q>p`. Consequently `q=p+r`, with `1<=r<p`, and `k=1`.
For `p=2`, the critical `q=4` violates coprimality, so there is no such shape.
The retained small-weight estimate below is nevertheless valid for `p=2`.

For a hypothetical positive integer cycle in this class, the inherited theorem
first forces the critical clock and `D=2^q-3^p=g>0`. It then follows that

    0<q log 2-p log 3<exp(-p/2500).

The logarithms in this last formula and the analytic argument are natural.

## 2. Plus-interior norm: exact alias audit

Put `rho=2^(1/p)`, `zeta_j=exp(2 pi i j/p)`, and

    Q(z)=1+z^(a+1)-z^a,       1<=a<=p-2,
    Y_j=|Q(rho zeta_j)|^2,
    mu=1+rho^(2a)+rho^(2a+2).

The three exponents are distinct modulo `p`. Hence `mean(Y)=mu<9`, and direct
expansion gives exactly

    Y_j=mu+2 Re(
        rho^(a+1) zeta_j^(a+1)
        -rho^a zeta_j^a
        -rho^(2a+1) zeta_j).

The nonvanishing of `Q` at these roots follows from irreducibility of `Z^p-2`
and `deg(Q)<p`. The triangle inequality gives `0<Y_j<25`, since both
`rho^a` and `rho^(a+1)` are less than `2`.

Write `Y_j-mu=sum_(m=1)^(p-1) C_m zeta_j^m`, combining all residues modulo `p`.
The coefficients are real and `C_(p-m)=C_m`. Parseval gives

    Var(Y)=mean((Y-mu)^2)=sum_(m=1)^(p-1) |C_m|^2.

The exhaustive cases are:

1. If `a!=1,p-2`, the only contribution at frequency `1` is
   `C_1=-rho^(2a+1)`. Frequencies `1` and `p-1` are distinct, so
   `Var(Y)>=2 rho^(4a+2)>2`.
2. If `a=1` and `p>=4`, then `C_1=-rho-rho^3`, again with the distinct
   conjugate frequency. Thus `Var(Y)>2`. At `p=4`, frequency `2` is
   self-conjugate, but it does not affect this coefficient or lower bound.
3. If `a=p-2` and `p>=5`, then
   `C_1=rho^(p-1)-rho^(2p-3)`. Its absolute value is
   `rho^(p-1)(rho^(p-2)-1)`. Since

       rho^(p-1)>=2^(4/5)>sqrt(3),
       rho^(p-2)>=2^(3/5)>3/2,

   it exceeds `sqrt(3)/2`, and the conjugate pair gives `Var(Y)>3/2`.
   Both strict rational-power comparisons follow from `2^8=256>243=3^5`.
4. At `p=3,a=1`, the colliding terms give
   `C_1=C_2=rho^2-rho-2`. Here `1<rho<3/2`, and the function
   `x^2-x` is increasing on this interval. Therefore
   `rho^2-rho<3/4`, so `C_1<-5/4`, in particular `|C_1|>1`.
5. At `p=4,a=2`, the self-conjugate coefficient is
   `C_2=-2 rho^2`. Its single Parseval contribution is
   `4 rho^4=8`. One must not double this self-conjugate contribution; the
   proposed argument correctly does not.

These cases cover every plus-interior type. In particular `Var(Y)>=1`
uniformly. The proof actually has a larger margin, but no larger constant is
needed.

## 3. Strong concavity gives the stated constant

On `(0,25)`, `d^2(log y)/dy^2=-1/y^2<=-1/625`. Taylor's upper bound around
`mu` is therefore

    log y<=log mu+(y-mu)/mu-(y-mu)^2/1250.

Averaging gives

    mean(log Y)<=log mu-Var(Y)/1250<log 9-1/1250.

There is no problem near `Y=0`: all sample values are positive, and the
second derivative only becomes more negative there. Since

    log|Norm(Q)|=(p/2) mean(log Y),

the conclusion is exactly

    |Norm(Q)|<exp(p log 3-p/2500)=c^p.

The factor `1/2500` thus includes the required factor `1/2` when converting
the squared-modulus average to the absolute norm. Since `g` divides the
positive prime-to-6 part of this norm, this bounds `g` as claimed.

## 4. Minus interiors and all endpoints

For a minus interior, `p>=3`, `1<=a<=p-2`, and
`Q(z)=z^(p-a)-z+1`. Its three exponent residues are distinct, so

    mean|Q|^2=1+rho^2+rho^(2(p-a))
             <5+2^(2/3)<7<c^2.

The final comparison follows from the exact lower bound

    c>3(1-1/2500)=7497/2500>299/100,
    (299/100)^2=89401/10000>7.

The elementary strict inequality `exp(-x)>1-x` for `x>0` was used.
AM-GM then gives `|Norm(Q)|<c^p`.

The plus endpoint `a=0` has norm `2` and hence prime-to-6 content `1`.
The minus endpoint `a=p-1` has norm `1`. Their contents are therefore `1`.

### 4.1 The two large endpoints are not admissible at k=1

The moved one has residue `a`, its predecessor has residue
`A=(a+r mod p)`, and its successor has residue `a'=(a-r mod p)`.
The adjacent gap increments are

    delta_before=(r+a-A)/p,
    delta_after =(r+a'-a)/p.

At `a=0`, one has `A=r`, hence `delta_before=0`. When `k=1`, the
preceding gap is `1` and the minus move would collide with the preceding one.
It is forbidden.

At `a=p-1`, one has `a'=p-1-r`, hence `delta_after=0`. When `k=1`,
the following gap is `1` and the plus move would collide with the next one.
It is forbidden.

Thus at `k=1` every actual swap has an interior type or a free endpoint type.
Sections 2–3 and the first part of this section already prove the stronger
all-`k=1` theorem. At `p=2` there are no interiors and every admissible type
is a free endpoint, so no small-weight exception is needed.

### 4.2 The proposed alternate endpoint estimates also check

Although unnecessary for admissible `k=1` words, the proposed estimates of
the formal endpoint gcd formulas are valid. The inherited exact formulas give:

    a=0, minus:   g<=2*3^(r-1)-1,
    a=p-1, plus:  g<=3^(r+1)-2^(r+1).

Both right sides are positive for `r>=1`. Because the critical logarithm is
not integral,

    r<(log_2 3-1)p+1<(3/5)p+1.

Thus either endpoint satisfies

    g<9*3^(3p/5)<9*2^p,

where `3^3=27<32=2^5`. For `p>=6`,

    (c/2)^p>(299/200)^p>=(299/200)^6>9.

The final exact integer certificate is

    299^6-9*200^6=138540961348201>0.

Consequently `g<c^p` for these endpoints when `p>=6`.

For every `2<=p<=5`, the inherited sharp general bound suffices:

    g<=3^p-2<c^p.

Indeed `c^p=3^p exp(-p/2500)>3^p-p*3^p/2500`, while
`p*3^p` is respectively `18,81,324,1215`, all strictly below `5000`.
This also closes any small-weight overlap among the cases.

## 5. Logarithmic consequence and limits

For a hypothetical positive integer cycle, `D=g>0`; therefore

    q log 2-p log 3
      =log(1+D/3^p)
      <D/3^p
      <(c/3)^p
      =exp(-p/2500).

Every inequality is strict. The first follows from `D>0`, and the logarithm
bound uses `log(1+x)<x` for positive `x`.

This is a necessary exponential approximation condition, not a contradiction
by itself. A separately justified effective lower bound for this particular
linear form in two logarithms could turn it into a finite cutoff. No such
lower bound or cutoff is asserted by this audit.

## 6. Reproducible exact checks

Run from the repository workspace:

    python3 session13-collatz-resultants/refinement/audit/verify_critical_gap.py
    python3 -O session13-collatz-resultants/refinement/audit/verify_critical_gap.py

The script uses Python's exact integers and fractions, plus SymPy's exact
resultant implementation. It imports no repository mathematical code and
contains explicit checks that remain active under `python -O`.

Default scopes:

- Directly expand `Q times conjugate(Q)` and audit every plus-interior Fourier
  alias through `p=400`, reducing coefficients exactly modulo `rho^p=2`.
- Evaluate all sparse norms through `p=20`, verify every endpoint norm, and
  compare every interior norm with the slightly stronger rational test
  `|Norm(Q)|<(7497/2500)^p<c^p`.
- Construct each actual cyclic unequal-adjacent swap at every coprime critical
  shape through `p=200`; recompute carry by affine composition; check the
  exact bound `g^2<7^p`, both legality equivalences against the actual adjacent
  bits, free endpoint cases, and absence of large endpoint types.
- Check all displayed numerical certificates using integers or fractions.
- Check inadmissibility of both large endpoints for every coprime `k=1`
  shape through `p=400`.
- Check `|Norm(Q)|^2<7^p` for every admissible type at each coprime critical
  shape through `p=20`.

The saved `receipt.json` records the resulting counts. Finite checks
corroborate the symbolic all-`p` argument; they do not establish its quantifiers
or imply a census of arbitrary Collatz parity words.

The default run passed, with byte-identical output under `python3 -O`:

- 79,401 exact plus-interior Fourier alias checks.
- 418 exact sparse norms, including 150 admissible critical norms satisfying
  `|Norm(Q)|^2<7^p`.
- 14,602 actual critical swap occurrences over 122 coprime critical shapes,
  including 122 cut-crossing occurrences.
- 24,758 direction-legality equivalence checks against the actual word bits.
- 244 free-endpoint occurrences and zero large-endpoint occurrences.
- 97,354 large-endpoint inadmissibility checks over all 48,677 coprime `k=1`
  shapes through `p=400`.

The comparison was reproduced by writing the normal and optimized outputs to
`receipt.json` and `receipt.optimized.json`, then running `cmp` on those files.
