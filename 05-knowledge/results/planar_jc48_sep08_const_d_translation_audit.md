# Independent audit of the finite-point constant-D translation

**Status: independent analytic/source audit PASS; normal, optimized,
and frozen 51-gate output agree.** This accepts the extension to
every finite position of the declared sixfold boundary zero. It
does not claim an automorphism of the fixed DG surface or remove
the prescribed multiplicity and nonzero-M hypotheses.

Audited [primary proof](planar_jc48_sep08_const_d_translation.md),
[source](../../04-computation/planar_jc48_sep08_const_d_translation.py),
and [output](planar_jc48_sep08_const_d_translation.out). The
[fixed-zero theorem](planar_jc48_sep08_const_d_quartic.md) is now
PROVED and independently audited on the current worktree. The
producer's earlier message describing that dependency as still
reserved was stale; its final proof has repaired the status text.

## 1. Exact statement and complete family recovery

The original surface is `W=(P1_x x P1_z)\{z=x^2}`, with
`s=z-x^2`, `t=1/s`. The hypothesis is a global `H=N/s^2` in
`L_2`, a global `L=M/s` in `L_1`, and `F=H^2+L`, with

    N|S=a(x-p)^6,   a!=0,   M(p,p^2)!=0,
    F|D constant,                     p in C finite.

The octic has exactly the finite sixfold point and the infinity
twofold point. No claim about a different partition is hidden
in the notation.

Set `w=x-p` but keep the original graph `z=(w+p)^2`.
Translation of this projective base coordinate preserves the
section degree boxes. Write

    N=a w^6+s B(w)+s^2 C(w),   s=z-(w+p)^2.

I independently recovered the two upper-coefficient constraints.
In the original `z` expansion, its linear coefficient is
`B-2(w+p)^2 C`, while its constant coefficient is
`a w^6-(w+p)^2 B+(w+p)^4 C`. Their simultaneous degree-four
bounds first force `C4=B6=0` and then `C3=B5=0`.
The next two coefficients give

    B4=a+C2,
    B3=C1+2p(C2-a).

Thus the primary's entire preconstant seven-parameter `N`
family is exhaustive, including all lower coefficients of `B`.
For `M`, the same direct expansion gives the stated six-parameter
preconstant family, with normal coefficient
`n0+n1 w+n2 w^2` and tangential coefficients
`m0,m1,d,n1+2p n2,n2`.

The actual infinity substitution remains

    w=1/r-p,       t=-r^2-r^4 b_D.

The coefficient of `b_D` in `H|D` is `a-C2`; that in `L|D`
is `-n2`. Constancy of `H|D^2+L|D` forces each restriction
constant, since a nonconstant square of degree at most four
cannot be cancelled by a polynomial of degree at most one.
Consequently `C2=a`, `n2=0`. The complete family becomes

    N=a w^6+s(beta0+beta1 w+b w^2+c1 w^3+2a w^4)
                             +s^2(c0+c1 w+a w^2),
    M=m0+m1 w+d w^2+n1 w^3+s(n0+n1 w),

exactly the full fixed-zero family in `(w,s)`. In particular,
the lower normal jets have not been dropped. The disappearance
of `p` follows from the constant-boundary constraints, not from
an assumed translation symmetry of all global functions.

The source's all-parameter rank checks are valid. There are nine
independent constraints on the fifteen-dimensional `N` box.
The selected constraint minor is a constant unit for every `p`,
and the displayed six-vector basis has determinant one in its
chosen six coefficient positions. The `M` basis likewise has
a constant determinant one and satisfies its sole constant-D
constraint. These checks rule out special complex values of `p`
where an argument using only generic rank would be insufficient.

## 2. The map and the exact no-mate consumer

Substituting `s=1/t` gives literal polynomials `H0(w,t),L0(w,t)`
with

    H(w+p,t)=H0(w,t),       L(w+p,t)=L0(w,t).

The source-plane transformation `(w,t)->(w+p,t)` has determinant
one. It therefore transports an unrestricted rational mate by
`G0(w,t)=G(w+p,t)` and preserves the nonzero constant Jacobian.
There is no restriction on the denominator or degree of `G`.

The same full polynomial pair is global on the auxiliary surface
`W0=(P1_w x P1_zeta)\{zeta=w^2}`, using
`zeta=w^2+1/t`. Its boundary octic is `a w^6`, and its finite
M-value is the original nonzero value `m0`. Its auxiliary
boundary constants are `c0-b` and `n0-d`. Thus every hypothesis
of the already proved fixed-zero theorem is paid, with no
further local or genus calculation needed. That theorem excludes
the translated rational mate and hence the original one.

The original boundary constants are instead

    H|D=c0-b+2p c1+4a p^2,
    L|D=n0-d+2p n1.

I checked that these agree with the actual old chart. Their change
under the auxiliary interpretation is harmless: the consumer
needs constancy, not equality of the two surfaces' boundary values.
The underlying source polynomial and its original fibre parameter
are unchanged by substituting `x=w+p`.

This is a complete intersection-of-function-spaces argument inside
the common rational field. It is stronger than merely moving a
formal local root and then applying a theorem whose global
hypotheses have not been checked. The fixed-zero proof handles
all its geometric generic components where needed; no new
irreducibility assumption enters this translation.

## 3. Nonextendability and precise failure boundaries

For the actual global function `b_D=-x^2-x^4 t`, the opposite
plane translation produces

    -(x-p)^2-(x-p)^4 t
       =-2p/r+5p^2-4p^3 r+p^4 r^2+(1-pr)^4 b_D.

The simple pole coefficient `-2p` is nonzero for `p!=0`.
This verifies the primary's explicit warning: the plane map
does not extend to an automorphism of `W`. Its inverse also
fails, so using the opposite translation direction in this
hostile is legitimate. The proof transports only the special
family after proving its new globality.

The literal identity for the inherited rational mate

    h=x^2+x^4t, F=h^4+h,
    G=1/[3x^3(4h^3+1)], J(F,G)=1

was independently replayed. Its octic has multiplicity eight,
outside the present pattern. The theorem also retains the
nonzero-M requirement and does not include an infinity-six
point, an infinity-eight point, or arbitrary other boundary
multiplicities. The named real/nonreal translations test the
symbolic formulas; they are not the proof of completeness.

## 4. Source audit and frozen replays

I read the entire frozen source, including its final lower-jet
and outside-family controls. It imports no inherited mathematical
implementation. Its rational chain-rule test is only a control;
the general determinant-one argument is supplied analytically
above. The symbolic minors prove all-point completeness, so
the named locations `0,1,-3/2,i` do not function as a census.

Independent commands:

```sh
python3 -B 04-computation/planar_jc48_sep08_const_d_translation.py
python3 -B -O 04-computation/planar_jc48_sep08_const_d_translation.py
```

Both pass **51 always-active gates**, and both outputs are
byte-identical to the frozen **373-byte** output. The source is
**7,273 bytes**, SHA256
`69e58278cb6b84ae791c004b22d59f2b9aebd0e191bb6c08f39d422bcdb53dce`.
All three output copies have SHA256
`62a3a4180658aa0fcdc9f0d0b7723fff387c5a87625be0f6f4bea33f50b0a0f5`.
The semantic digest is
`396216fbb707e9ba5633c9bbacc67cf1f4e58b2c14a5d4876ddf5780909d1950`.

The updated final primary text was read after it removed the stale
dependency-status conditional. No mathematical or source correction
remains. This audit accepts the theorem for promotion by its owner;
it makes no changes to the frozen source/output or prior theorem.
