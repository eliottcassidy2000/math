---
id: THM-3781
title: "Scalar-centred common-step-three Danielewski Darboux nonentry"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  On every
  exponent-two squarefree Danielewski surface with at least two arms, no
  Darboux pair can have three nonconstant homogeneous weight pieces in each
  output, common step three, and scalar-centred support.  Up to exchanging
  the outputs only two support placements survive the endpoint signs.
  Their adjacent-weight equations force, respectively, a polynomial Mobius
  contradiction and a cubic-over-linear divisibility contradiction before
  the scalar equation is reached.  This closes two complete three-by-three
  cells, not planar JC(2).
source: root / planar-Jacobian Danielewski Darboux session, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The support convolution, endpoint power
  laws, both adjacent-weight integrations, arm-valuation thresholds, and
  terminal degree contradictions were derived by hand and checked by an
  exact companion.  Normal, optimized, and frozen transcript comparison
  remains required together with an independent hostile audit.
depends_on:
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
related:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
script: 04-computation/jc2_common_step_three_by_three_danielewski_thm3781.py
output: 05-knowledge/results/jc2_common_step_three_by_three_danielewski_thm3781.out
script_sha256: 893d3f5f701c2008c50270839dc0949451ee60804e0a212e968c14eb6fd1d2fc
output_sha256: e8ac756513895f4d36676245d4e1e4c6f2fc30371f6c1965cf51f511960dffc0
semantic_sha256: 34f8df9c6b069bd5f36939f38109908fde5ab2d4dd7374b019445add3e400982
hash_basis: raw LF bytes
---

# THM-3781 -- two complete three-by-three Darboux cells are empty

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  The
exponent-two Danielewski route first survives the known support census at
three weight pieces in each output.  This theorem closes the next complete
common-step family, rather than one symmetric representative.

Work over `C`.  Let `Sigma in C[b]` be squarefree with at least two distinct
roots, and put

```text
A_Sigma=C[b,c,e]/(c^2e-Sigma(b)).                       (1)
```

Use the Poisson bracket and weights

```text
{b,c}=c^2,          {c,e}=-Sigma'(b),
{b,e}=-2ce,
wt(b)=0,            wt(c)=1,             wt(e)=-2.     (2)
```

The weight-`r` part is

```text
(A_Sigma)_r={c^r f(b):
  Sigma^ceil(-r/2) divides f when r<0}.                 (3)
```

For homogeneous pieces the exact coefficient rule is

```text
{c^r f,c^s g}=c^(r+s+1)(s f'g-r f g').                (4)
```

Subtract constants from both proposed outputs.  Suppose each then has
exactly three nonzero weight pieces, its weights form an arithmetic
progression of step three, and the middle weights `p,q` obey

```text
p+q+1=0.                                               (5)
```

Thus output weight zero is the central coefficient of the support
convolution.  The theorem says no such `P,Q` satisfy

```text
{P,Q}=1.                                               (6)
```

## 1. Only two support placements survive

The extreme output weights occur once, so their brackets must vanish.  A
negative and a positive homogeneous piece cannot have zero bracket.  Indeed,
at a root of `Sigma`, write their coefficient orders as `A>0,B>=0`.  If
`r<0<s`, the first surviving coefficient in `(4)` is

```text
sA-rB>0.                                               (7)
```

A zero-weight endpoint would force its coefficient to be constant and hence
removable.  Consequently the two lower endpoints must both be negative and
the two upper endpoints both positive.  With `(5)`, this leaves

```text
p in {-2,-1,0,1}.                                      (8)
```

Exchanging `P,Q` identifies `p` with `q=-p-1`.  Hence only the following two
cells remain:

```text
I:   supp(P)=(-4,-1,2),       supp(Q)=(-3,0,3),
II:  supp(P)=(-5,-2,1),       supp(Q)=(-2,1,4).         (9)
```

Both have output-weight multiplicities

```text
1,2,3,2,1  at weights  -6,-3,0,3,6.                  (10)
```

## 2. Cell I forces its zero-weight piece to be constant

Write

```text
P=c^-4 f+c^-1 a+c^2 F,
Q=c^-3 g+h+c^3 H.                                    (11)
```

The unique weight `-6` and `6` equations are

```text
-3f'g+4fg'=0,                 3F'H-2FH'=0.            (12)
```

Unique factorization and coprime exponent pairs give

```text
f=A k^4,       g=B k^3,       F=L ell^2,       H=M ell^3,              (13)
```

for nonzero constants `A,B,L,M` and nonzero polynomials `k,ell`.  The
negative-weight membership conditions in `(3)` force `Sigma|k`, so `k` and

```text
z=k ell                                                   (14)
```

are nonconstant.

The output weights `-3` and `3` now reduce exactly to

```text
(a/k)'=(4A/(3B))h',              (a ell)'=(2L/(3M))h'. (15)
```

After integration, for constants `C,D` and nonzero constants `lambda,mu`,

```text
a=k(lambda h+C),                 a ell=mu h+D.         (16)
```

Eliminating `a` yields

```text
(lambda z-mu)h=D-Cz.                                  (17)
```

If `h` were nonconstant, the left side would have degree
`deg(z)+deg(h)`, while the right side has degree at most `deg(z)`.  Thus `h`
is constant.  But constants were subtracted before declaring the support,
so the weight-zero piece disappears and `(11)` is not a three-by-three
cell.  Cell I is empty.  Notice that the scalar weight equation was never
used.

## 3. Cell II creates an impossible cubic denominator

Now write

```text
P=c^-5 f+c^-2 a+cF,
Q=c^-2 g+ch+c^4 H.                                    (18)
```

The endpoint equations give

```text
f=A k^5,       g=B k^2,       F=L ell,       H=M ell^4.                (19)
```

Again `(3)` forces `Sigma|k`, so `z=k ell` is nonconstant.  The two adjacent
output equations integrate to

```text
a/k^2=lambda k h+C,               h/ell=mu a ell^2+D,                 (20)
```

where

```text
lambda=5A/(2B),                    mu=4M/L.             (21)
```

Substitution of the second equation into the first gives the entire
obstruction in one line:

```text
a(1-lambda mu z^3)=k^2(lambda D z+C).                 (22)
```

The factor `1-lambda mu z^3` is coprime to `k`: modulo every factor of `k`
it is one.  It must therefore divide `lambda D z+C`.  Its degree is
`3deg(z)`, while the latter polynomial has degree at most `deg(z)`.  Hence
the numerator is zero.  Nonconstancy of `z` forces `C=D=0`, after which
`(22)` forces `a=0`, contradicting the declared middle weight.  Cell II is
empty, again before the scalar equation.

## 4. What this adds and what remains

THM-3572 closed the parity-compatible step-two cell

```text
(-3,-1,1) x (-2,0,2).                                 (23)
```

The present theorem closes both inequivalent scalar-centred step-three cells
uniformly for every squarefree multiarm polynomial `Sigma`.  The mechanism
is new: step two ended in a common-factor contradiction, while step three
turns the adjacent equations into a Mobius relation in Cell I and a
cubic-over-linear divisibility in Cell II.

This is not a classification of all three-by-three supports.  Unequal step
sizes, non-centred scalar placements, gapped supports, and common steps at
least four remain open.  It constructs no Darboux pair and no planar
Jacobian counterexample.  **QED, conditional on independent hostile audit.**
