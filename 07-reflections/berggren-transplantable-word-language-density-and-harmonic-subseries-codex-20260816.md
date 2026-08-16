# One third of the ternary Berggren tree is wordwise T4-calibratable

**Status: PROVED in THM-3497 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**
The exact companion is
`04-computation/berggren_transplantable_word_language_harmonic_probe_20260816.py`.
The result classifies individual composite addresses.  It does not construct
one branch-monoid-equivariant calibration, a physical current, or a Jacobian
map.

## 1. From the failed generator transplant to a word language

The preceding full-branch probe proves that the projective actions of `A` and
`C` cannot individually be affine owner actions with their true mod-two
linear parts, while `B` can.  Composition creates a different question.

For a branch word `w`, write

```text
rho_3(w) in S(P1(F_3))=S4,
L_2(w) in <J>={I,J}.                                     (1)
```

Call `w` **wordwise T4-calibratable** if there exist a point bijection
`f_w:P1(F_3)->V4` and a translation `t_w in V4` such that

```text
f_w rho_3(w) f_w^(-1)(u)=L_2(w)u+t_w.                    (2)
```

Both `f_w` and `t_w` may depend on the entire word.  This quantifier is
essential: (2) is weaker than a common calibration for all letters and
stronger than mere equality of four-point cardinalities.

## 2. The twelve-state perfect-matching automaton

Let

```text
q:S4 -> S3                                                   (3)
```

be the action on the three perfect matchings of `K4`.  Its kernel is the
normal translation subgroup `V4`.  The three Berggren letters have matching
actions

```text
q(A)=q(C)=r=(0 1 2),              q(B)=s=(0 1),           (4)
```

and true linear bits

```text
epsilon(A)=epsilon(B)=1,          epsilon(C)=0,            (5)
```

where bit one means `J` and bit zero means `I`.  Thus every word has the
twelve-state address

```text
Phi(w)=(mu(w),epsilon(w)) in S3 x C2.                     (6)
```

All twelve states occur.  Indeed `(r,0)` is supplied by `C`, `(r,1)` by `A`,
so their quotient supplies the pure bit `(1,1)`; after stripping the bit from
`B`, `r` and `s` generate `S3`.

The exact acceptance rule is

```text
w passes (2)
  iff epsilon(w)=0 and mu(w)=1,
   or epsilon(w)=1 and mu(w) is one of the three reflections. (7)
```

To prove (7), first suppose `L_2(w)=I`.  The four affine maps `u->u+t` are
the identity and three double transpositions, exactly `ker(q)=V4`.  Suppose
instead `L_2(w)=J`.  Its four translates consist of two transpositions and
two four-cycles.  Those are precisely the two odd conjugacy types in `S4`,
and every odd permutation maps to a reflection under (3).  Conjugacy in
`S4` is determined by cycle type, proving both directions.  The companion
checks this equivalence on all `48` pairs in `S4 x <J>`.

There are four accepting states out of twelve:

```text
(1,0),       (s_1,1), (s_2,1), (s_3,1).                  (8)
```

This is the exact tournament content.  The only quotient needed to decide
wordwise calibratability is the action on the three `K4` matchings, together
with the one true-linear-action bit.

Partition refinement leaves all twelve states inequivalent.  Hence this DFA
is minimal, and its transition/syntactic monoid is the regular left action of
the same twelve-element group `S3 x C2`.  Acceptance is not multiplicative:
`B` and `BC` pass while `BBC` fails.  Thus “language,” rather than
“submonoid,” is load-bearing here.

## 3. Exact level count and recurrence

Let `a_n` be the number of accepted words in `{A,B,C}^n`.  The first values
are

```text
n:    0  1  2   3   4   5    6    7     8     9
a_n:  1  1  3  11  25  81  251  715  2193  6593.         (9)
```

The twelve-state transfer matrix has characteristic polynomial

```text
(lambda-3)(lambda-1)^4(lambda+1)^3
              (lambda^2+2lambda+3)^2.                   (10)
```

An exact Krylov-space calculation gives the smaller scalar recurrence

```text
a_n=2a_(n-1)+2a_(n-2)+6a_(n-3)-9a_(n-4),                (11)
```

for `n>=4`, with initials `(1,1,3,11)`.  Equivalently,

```text
sum_(n>=0) a_n x^n
 = (1-x-x^2-3x^3)/((1-x)(1-3x)(1+2x+3x^2)).             (12)
```

There is also an integral closed form.  Define

```text
c_0=1,       c_1=-1,       c_n=-2c_(n-1)-3c_(n-2).       (13)
```

Then

```text
a_n=(3^n+1+c_n)/3,
c_n=Re((-1+i sqrt(2))^n).                                (14)
```

Since `|c_n|<=3^(n/2)`, (14) proves

```text
a_n/3^n = 1/3+O(3^(-n/2)).                               (15)
```

Thus the individually calibratable addresses have exact asymptotic density
one third in the ternary ancestry levels.  This density is structural, not a
numerological ratio of two finite sets.

## 4. The Fibonacci three-ray restriction is half of the mod-six clock

THM-3339 places consecutive Fibonacci parameters on the three address rays

```text
u_(3r+2):       (BA)^r,
u_(3r+3):       A(BA)^r,
T u_(3r+4):     C(BC)^r.                                 (16)
```

Applying (7) gives the exact period-two gate

```text
(BA)^r passes       iff r is even;
A(BA)^r passes      iff r is odd;
C(BC)^r passes      iff r is odd.                         (17)
```

Equivalently, the Fibonacci index passes exactly when

```text
k mod 6 in {0,1,2}.                                      (18)
```

So the full ternary tree has density `1/3`, while the distinguished
Fibonacci three-ray family meets the language on exactly half of its index
clock.  These are different universes and must not be blended.

On Fibonacci **indices**, (18) gives

```text
sum_(k<=N, k mod 6 in {0,1,2}) 1/k
       =(1/2) log N+C_F+o(1).                             (19)
```

By contrast, weighting by the Fibonacci values or by the exponentially
growing triple coordinates yields a convergent series.  The harmonic claim
is about the index subset, not geometric mass on the triples.

## 5. Shortlex realization as a subset of the natural numbers

Enumerate the ternary tree in shortlex order: the empty word is `1`, and the
`3^n` words of level `n` occupy the next consecutive block.  Let `S` be the
set of indices of words accepted by (7).

The transfer matrix is annihilated by the squarefree polynomial

```text
(lambda-3)(lambda-1)(lambda+1)(lambda^2+2lambda+3).       (20)
```

Its Perron eigenvalue is the simple eigenvalue `3`; every other eigenvalue
has modulus at most `sqrt(3)`.  The same estimate holds uniformly after any
fixed prefix.  Decomposing a shortlex initial segment into at most two full
ternary cylinders of each remaining depth therefore gives

```text
|S intersect [1,N]|=N/3+O(sqrt(N)).                       (21)
```

Partial summation now yields

```text
sum_(m<=N, m in S) 1/m
       =(1/3) log N+C_S+O(N^(-1/2)).                      (22)
```

This is the useful refinement of “every subset of the naturals is a subset
of the harmonic series.”  The Boolean encoding alone is tautological; the
automaton supplies the nontrivial logarithmic coefficient, discrepancy, and
finite recurrence.  Ordinary harmonic mass forgets which ancestry words
produce that coefficient.  Recovering the language itself requires the
twelve-state automaton or an equivalent character bank.

## 6. Boundaries and next tests

- (7) permits a new gauge for each composite word; it does not contradict the
  generatorwise no-go.
- Acceptance need not be closed under concatenation, so this is a regular
  language, not automatically a submonoid.
- The matching quotient is intrinsic to `K4`; orienting a tournament adds an
  order sidecar that (7) does not use.
- The `1/3` tree coefficient and `1/2` Fibonacci-index coefficient are not
  LRC densities and carry no Jacobian consequence.

The next high-value refinement is to repeat the count for THM-3339's *fixed*
affine cocycle rather than allowing `t_w` to vary.  That retains more current
data and is expected to expose a parity-split rather than the clean density
(15).
