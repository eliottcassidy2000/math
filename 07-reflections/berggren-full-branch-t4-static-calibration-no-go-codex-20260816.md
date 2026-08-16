# The full ternary Berggren packet and the T4 atlas have different connections

**Status: PROVED in THM-3497 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**
The finite claims are reproduced by
`04-computation/berggren_full_branch_t4_static_calibration_no_go_20260816.py`.
This note constructs no physical LRC current, proves no Jacobian statement,
and does not identify Berggren ancestry with either problem.

## 1. Inheritance pass

The closest proved mechanism is THM-3487's distinction between a finite atlas
and a connection on that atlas.  THM-3339 supplies the three parameter
matrices

```text
A=[0 1;-1 2],       B=[0 1;1 2],       C=[1 0;2 1],       (1)
```

and their true action on the nonzero `V4=F_2^2` directions:

```text
Abar=Bbar=J=[0 1;1 0],              Cbar=I.               (2)
```

The canonical hostile example is `AGL(2,2)=S4`.  It says that every
permutation of four points is affine after some affine structure is chosen,
so bare four-point affineness cannot distinguish the two packets.  The
corrected near miss is therefore the tempting cardinality argument

```text
6*|P1(F_3)|=6*|V4|=24  =>  same branch-equivariant T4.     (false) (3)
```

The least-used load-bearing sidecar is the prescribed linear part (2).  Once
it is retained, translations cannot freely repair the action.

## 2. The typed calibration problem

Use

```text
P1(F_3)=(x0,x1,x2,x3)
       =([1:0],[0:1],[1:1],[1:2]).                        (4)
```

Reduction of (1) modulo three gives

```text
A_3=(0 1 3), fixing 2;
B_3=(0 1 3 2);
C_3=(0 3 2), fixing 1.                                    (5)
```

Let `Ord(V4*)` be the six orders of the three nonzero directions.  The two
twenty-four-state carriers are

```text
X=Ord(V4*) x P1(F_3),
Y=Ord(V4*) x V4.                                          (6)
```

Both sets index all labelled transitive tournaments on four vertices: the
first factor orders the three edges out of an owner and the second factor
selects the owner.  This is a set-level atlas statement.

For a branch letter `h`, source transport on (6) is

```text
(pi,x) |-> (hbar pi,h_3 x),                               (7)
```

whereas an affine owner lift with arbitrary translation `t_h` is

```text
(pi,u) |-> (hbar pi,hbar u+t_h).                          (8)
```

A static point calibration is a bijection
`f:P1(F_3)->V4` intertwining (5) and (8).  The more permissive
frame-preserving calibration allows six bijections `f_pi` and asks

```text
f_(hbar pi) h_3=(hbar(-)+t_h) f_pi                        (9)
```

for every frame.  Equation (9) is the finite descent demanded after the
free-tree recursion is quotiented to the six mod-two frames.

## 3. Exact generator census

The exhaustive census of all `24` point gauges gives:

| letter | type on `P1(F_3)` | true linear part | possible affine-owner types | static gauges, allowing `t` |
|---|---:|---:|---:|---:|
| `A` | `1*3` | `J` | `1^2*2` or `4` | `0` |
| `B` | `4` | `J` | `1^2*2` or `4` | `8` |
| `C` | `1*3` | `I` | `1^4` or `2^2` | `0` |

For THM-3339's pinned lifts

```text
Atilde(u)=Ju+p,       Btilde(u)=Ju+p,       Ctilde(u)=u+p, (10)
```

the counts are `(0,4,0)`.  In the displayed point orders

```text
x0,x1,x2,x3  <->  0,p,q,r,                               (11)
```

the `B` actions agree without any further gauge:

```text
B_3(u)=Ju+p.                                              (12)
```

Thus `B` is a genuine positive branch transplant, not merely a nonvanishing
count.  The four pinned gauges are the centralizer of its four-cycle; allowing
either of the two four-cycle translations gives eight.

The stronger twenty-four-state cycle types are

```text
source A:  2^3 6^3;       target A: 2^12 or 4^6;
source B:  4^6;           target B: 2^12 or 4^6;
source C:  1^6 3^6;       target C: 1^24 or 2^12.         (13)
```

So allowing the calibration to depend on the six-state frame still cannot
repair `A` or `C`.  For each lawful four-cycle translation, however, `B`
admits exactly

```text
8^3=512                                                     (14)
```

frame-preserving calibration families: one of eight compatible return gauges
on each of the three two-cycles of the frame swap.

## 4. The obstruction is the missing 3-primary linear mode

The failure is not lack of an affine model.  Under the identity point gauge,
each permutation in (5) has a unique affine decomposition because
`AGL(2,2)=S4`.  The linear parts required by the two three-cycles `A_3` and
`C_3` have order three.  They disagree with their actual direction actions
`J` and `I` in (2).

Translations in `V4` can change fixed points and turn a `J` lift into either
a transposition or a four-cycle.  They cannot manufacture an order-three
mode over a fixed `J` or `I` quotient.  Equivalently,

```text
<A_3,B_3,C_3>=S4                 has order 24,
<Atilde,Btilde,Ctilde>=V4 semidirect <J> has order 8.      (15)
```

This is a prime-support tariff: the projective `3`-torsion discarded by the
mod-two direction quotient must be restored as a sidecar if one wants a
single full-branch finite connection.

## 5. Tree, recurrence, and tournament boundary

There is no contradiction with recursive transport on the full Berggren
tree.  A node has one parent, so after choosing a root gauge one may define
each child gauge from (9).  Ordinary `H^1` of the tree is zero.  The no-go
appears only when gauges must depend on the six-state frame rather than the
entire ancestry address, or when a single static gauge is demanded.  That
finite quotient introduces return equations, already violated by `A` and
`C`.

The surviving `B` spine is especially concrete.  Its parameter recurrence is

```text
(m_(n+1),n_(n+1))=(n_n,m_n+2n_n),                        (16)
```

the Pell/silver-ratio recurrence.  Hence the exact positive transplant lives
on a distinguished binary-free ray of the ternary tree, while the other two
letters require a moving ancestry calibration or a new `C3` sidecar.

The tournament statement is consequently precise:

```text
same 24 transitive T4s:       yes;
same B-branch connection:     yes, after one of the gauges above;
same full ternary connection: no.                         (17)
```

THM-3487's repaired Fibonacci loop is time-dependent and composite, so its
positive transplant does not imply generatorwise compatibility with (5).

## 6. Connection ledger and next decisive object

- source: the projective branch action on `P1(F_3)`;
- target: affine owner transport on `V4`;
- proposed map: a static point gauge, or the six-gauge family (9);
- preserved predicate: branch action and the true mod-two channel order;
- destroyed information if (2) is dropped: the actual matching/direction
  quotient;
- missing sidecar: an order-three mode or the full ancestry-dependent gauge;
- cheapest decisive test: cycle type at one generator, strengthened by the
  twenty-four-state return census (13).

THM-3497 completes the next finite object proposed here: the regular language
of branch words `w` for which
the *composite* projective action `w_3` has a cycle type realizable by an
affine map with the true composite linear part `wbar`.  It identifies exactly
which subset of the ternary tree admits wordwise T4 calibration and turns its
level counts into a recurrence and a harmonic-subseries statement, while
retaining the quantifier boundary above.
