---
id: THM-3974
title: "Height-tower few-weight Darboux support obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For every
  height n>=2, the exact-volume completion algebra B_n has no homogeneous,
  one-by-arbitrary, or two-by-two polynomial Darboux pair in its natural
  weight grading. At height n=2, the exact exponent-two color at u=0 and the
  second compulsory color u=-1 transfer the complete THM-3569, THM-3579,
  THM-3583, and repaired THM-3592 support obstructions: the 2x3, 3x2, 2x4,
  4x2, and 3x3 cells are empty. Hence every height-two Darboux pair would
  need at least seven retained nonconstant weight pieces; the first live cells are
  exactly 2x5, 3x4, 4x3, and 5x2. No Darboux pair or JC(2) counterexample is
  constructed or excluded in unrestricted support.
source: jc-degree6-one-place / post-THM-3973 exact-volume completion support analysis, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (all-frontiers, 2026-08-24). The audit
  rederived the two-ceiling weight pieces, Wronskian bracket, homogeneous,
  one-by-arbitrary, and uniform two-by-two gates. It traced every imported
  height-two ladder and verified that only the distinguished exponent-two
  arm is used for support selection; the second compulsory color supplies
  exactly the terminal degree-two divisor. The missing direct THM-3579
  equal-step dependency and transfer-ledger row were restored before
  promotion. Normal and optimized 1,040-gate runs byte-match the frozen
  output, and all hashes pass.
depends_on:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3579-equal-step-three-by-three-danielewski-darboux-nonentry
  - THM-3583-universal-exponent-two-two-by-four-weight-darboux-nonentry
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
related:
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3971-canonical-debt-determinantal-affine-plane-completion
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
script: 04-computation/jc2_height_tower_weight_support_thm3974.py
output: 05-knowledge/results/jc2_height_tower_weight_support_thm3974.out
script_sha256: 1d3786ac96570159659264ec84229147119b85d53c014cd0e286194294ef0084
output_sha256: 283ff98b6ab888c6d2efac073601ebdbbda053059f87b870673474256625dbf4
semantic_sha256: 954a2edb1f0aac4079d9a7bb2c94671121bb58198b120965ed77d4655a151c68
hash_basis: raw LF bytes
---

# THM-3974 -- exact volume reaches two brackets, but not a few-weight pair

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over
an algebraically closed field `k` of characteristic zero. For `n>=2`, put

```text
u=x^n t,                 z=1+u,
p=x^(-n)u(u+1)=zt,
y=x^(-n-1)u^2(u+1)=x^(n-1)zt^2,
B_n=k[x,u,p,y] subset k[x,t].                           (1)
```

This is the height tower underlying THM-3973. Give it the grading

```text
wt(x)=1,       wt(t)=-n,       wt(u)=wt(z)=0,
wt(p)=-n,      wt(y)=-(n+1).                            (2)
```

The source Jacobian bracket has degree `n-1`. The theorem proves

```text
all n>=2:  no 1 x arbitrary pair and no 2 x 2 pair;

n=2:       no 2 x 3, 3 x 2, 2 x 4, 4 x 2, or 3 x 3 pair;
            every possible pair has at least seven retained pieces.      (3)
```

A *retained piece* means a nonzero homogeneous coefficient after deleting a
scalar weight-zero summand, which brackets to zero. The result is a support
obstruction, not a finite-degree search and not a proof of `JC(2)`.

## 1. Exact homogeneous pieces

Every monomial of weight `r` has the form

```text
x^c p^a y^b
 =x^r u^(a+2b)(u+1)^(a+b),
r=c-na-(n+1)b,                 a,b,c>=0.                (4)
```

Since `u` has weight zero, the weight pieces are `k[u]`-modules. For
`r>=0`, the monomial `x^r` shows immediately that

```text
(B_n)_r=x^r k[u].                                      (5)
```

For `r=-q<0`, the condition in `(4)` is `na+(n+1)b>=q`. The least possible
orders at `u=0` and `u=-1` are respectively

```text
min(a+2b)=ceil(q/n),          min(a+b)=ceil(q/(n+1)).  (6)
```

The first minimum is attained using powers of `p`, and the second using
powers of `y`. The coefficient ideal is an ideal of the PID `k[u]`; its gcd
therefore has exactly the two minima in `(6)`. Hence

```text
(B_n)_(-q)
 =x^(-q) u^ceil(q/n)(u+1)^ceil(q/(n+1)) k[u].          (7)
```

In particular every negative piece vanishes at both colors `u=0,-1`.

For homogeneous elements, direct differentiation gives

```text
J_(x,t)(x^r f(u),x^s g(u))
 =x^(r+s+n-1) W_(r,s)(f,g),
W_(r,s)(f,g)=r f g'-s f'g.                             (8)
```

This freezes both the support convolution and every valuation used below.

## 2. Homogeneous and one-by-arbitrary nonentry

A homogeneous scalar bracket would require

```text
r+s=1-n.                                               (9)
```

If both weights are negative, `(7)` makes both coefficients divisible by
`u(u+1)`, so their Wronskian is still divisible by `u(u+1)`. If, after
swapping, `r>=0` and `s=-q<0`, then `q=n+r-1` and the negative coefficient
`g` has degree at least two. For `r=0`, the Wronskian is divisible by `g`.
For `r>0`, its leading coefficient is

```text
[r deg(g)+q deg(f)] lc(f)lc(g),                         (10)
```

which is nonzero, and its degree is positive. Thus no homogeneous bracket is
a nonzero scalar.

If one output has one retained weight, its brackets with distinct weights of
the other output occupy distinct output weights. All but the unique possible
scalar complement vanish separately; deleting them leaves the forbidden
homogeneous scalar pair. This closes `1 x m` and `m x 1` for arbitrary
finite `m`.

## 3. Uniform two-by-two nonentry

Suppose both outputs have two retained weights. A single scalar contributor
would itself be a forbidden homogeneous pair, so both contributors occur on
the anti-diagonal of the support rectangle. The two outer rows must vanish.
The sign and zero-weight cases in `(8)` show that the lower outer pair is
negative. The upper pair cannot also be negative, because then both scalar
summands still vanish at `u=0`; and a mixed-sign outer pair cannot have zero
Wronskian. Hence the upper pair is nonnegative. Thus, for `R,T>=n-1`, write

```text
P=x^(-R)f+x^(T+1-n)F,
Q=x^(-T)g+x^(R+1-n)G.                                  (11)
```

The lower and upper outer equations are

```text
R f g'-T f'g=0,
(T+1-n)F G'-(R+1-n)F'G=0.                              (12)
```

Evaluate the scalar row at the distinguished color `u=0`. In the
`(-R,R+1-n)` summand, `(7)` gives

```text
ord_0(f)>=ceil(R/n).                                   (13)
```

For a nonzero scalar value, `f` must be simple and `G` a unit. The endpoint
`R=n-1` does not survive: its complement has weight zero and the Wronskian
`-R fG'` remains divisible by `f`. Since `R>=n-1`, the only surviving value
is `R=n`. Symmetrically the other scalar summand can survive only for `T=n`.
Hence

```text
R=n                         or                         T=n. (14)
```

Assume `R=n`; the other case is symmetric. Dividing the first equation in
`(12)` by `fg` and using unique factorization gives, with `d=gcd(n,T)`,

```text
f=A h^(n/d),              g=B h^(T/d),                 (15)
```

for `A,B in k*`. If `T!=n`, the `f` summand is the sole possible scalar
address and is simple at `u=0`. If `T=n`, simplicity of either summand makes
the same common base simple. Thus `(n/d)ord_0(h)=1`, so

```text
d=n,                 T=nk,
f=A h,               g=B h^k.                          (16)
```

Membership of the weight-`-n` coefficient at both colors forces

```text
u(u+1) divides h,                deg(h)>=2.             (17)
```

Put `p_0=n(k-1)+1`. The upper equation in `(12)` has weights `p_0,1`, so

```text
F=L K^p_0,                    G=M K                     (18)
```

with `L,M in k*`. Up to an overall sign, the scalar row is

```text
S=f'G+n fG'-nkF'g-p_0Fg'.                              (19)
```

Substitution of `(16)` and `(18)` factors it exactly:

```text
S=(K h'+n hK')
  [AM-k p_0 LB K^(p_0-1)h^(k-1)].                      (20)
```

The first factor has leading coefficient

```text
[deg(h)+n deg(K)]lc(h)lc(K) !=0                        (21)
```

and positive degree by `(17)`. It cannot divide the unit `1`. This proves
the uniform `2 x 2` nonentry.

## 4. Height two inherits the complete exponent-two few-support gates

Now specialize to `n=2`. Four exact properties are the complete interface
used in the symbolic proofs of THM-3569, THM-3579, THM-3583, and the repaired
THM-3592:

```text
(I)   the coefficient ring k[u] is a characteristic-zero UFD;
(II)  the bracket is x^(r+s+1)W_(r,s);
(III) at u=0, a weight -q coefficient has order >=ceil(q/2);
(IV)  every retained negative coefficient also vanishes at u=-1.         (22)
```

Properties `(I)--(III)` reproduce verbatim the singleton rows, support-sum
collisions, scalar-arm gate `(-2,1)/(1,-2)`, common-power equations,
first-order bridges, hidden gcd branches, and arm-order mismatches in those
three theorems. Property `(IV)` supplies the one global input that their
terminal Euler factors require: whenever a negative endpoint has the UFD
form `A h^m`, both irreducibles `u` and `u+1` divide `h`, so

```text
deg(h)>=2.                                               (23)
```

No multi-arm alternation is being imported silently. A ladder is selected
at the single distinguished exponent-two color `u=0`; after its absolute
weights are fixed, the displayed polynomial identities in the cited theorem
contradict the scalar row globally. The second color is used only in `(23)`.
In the repaired THM-3592 `d=3` hooked overlap, the low singleton already ties
the two candidate negative coefficients to one common base at `u=0`, so its
no-alternation seam remains valid.

For auditability, the transfer ledger is:

| dependency | imported exact mechanism | conclusion on `B_2` |
|---|---|---|
| THM-3569 | two support complements; lower/upper ladders; UFD bridges; Euler factors | no `2 x 3` or `3 x 2` |
| THM-3579 | equal-step three-by-three ladder | no equal-AP `3 x 3` |
| THM-3583 | four-vertex component path; six `R=2/T=2` ladders; square/gcd-five/gcd-four modes | no `2 x 4` or `4 x 2` |
| repaired THM-3592 | remaining three-point catalogue; component deletion; hooked/reflected/Euclidean ladders | no remaining `3 x 3` |

Every use of `Sigma^ceil(q/2)` in those proofs is a use of the local order in
`(III)`. Every use of `Sigma|h` at the terminal step is replaced by
`u(u+1)|h` from `(IV)`. No coefficient identity, support symmetry, or
orientation is changed.

It follows that a height-two Darboux pair has at least seven retained
nonconstant pieces. At total support seven, the only cells not excluded are

```text
(2,5),              (3,4),              (4,3),              (5,2). (24)
```

## 5. Positive control and exact residual

The support theorem is not a bracket-span obstruction. Put

```text
w=2u+1=2z-1.                                           (25)
```

For every `n>=2`, direct use of `(8)` and `w^2=1+4u(u+1)` gives the exact
two-bracket identity

```text
1=J(-(2z-1)p,x)+J((6/(n-1))xp,z).                      (26)
```

The first bracket is `1+6u(u+1)` and the second is `-6u(u+1)`. Thus the
constant already has bracket length at most two. The theorem proves that
compressing it to one bracket requires substantially mixed weight support:

```text
all n>=2: first unclosed cell is 2 x 3 (or 3 x 2);
n=2:     first unclosed cells are exactly those in (24).                  (27)
```

For `n>2`, the exponent-two transfer is not asserted: the distinguished
order becomes `ceil(q/n)`, changing the scalar and ladder arithmetic. For
`n=2`, cells of total support at least seven remain open. No unrestricted
Darboux pair, finite cubic map, or counterexample to `JC(2)` is proved or
refuted here. **QED.**

## Reproduction

```bash
python3 04-computation/jc2_height_tower_weight_support_thm3974.py
python3 -O 04-computation/jc2_height_tower_weight_support_thm3974.py
python3 agents/check_docs.py
```
