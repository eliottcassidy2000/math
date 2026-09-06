# Independent audit of the ternary double-pair Smith law

**Verdict: PASS.** I reviewed the all-depth argument in
[the double-pair report](overnight2_20260906_smith_double_pair.md), inspected
its complete producer, and independently reconstructed all 923 polynomial
minors with SymPy's polynomial-domain determinant engine. All 14,501
nonzero coefficient terms satisfy the stated universal lower bounds, and
all thirteen frontier forms have actual minors with a unique least term.
The thirteen displayed factorizations were also checked individually,
including the permitted determinant signs. No producer code was imported.

The reviewed source SHA-256 is
`bd2b5d2a9e080dd7a5422892ef5a22a4798282c7ad593376712d03e1aa186a6d`;
its frozen output SHA-256 is
`78fbfe52e3c8e7c51eda03b8b449beaffbe23469848ecc34d2b6c9a478a16fdd`.
Both matched the actual files. This review does not depend on the other
three-node candidate: it was provisional at the producer's inheritance
pass and was subsequently promoted as THM-4429 on the incoming remote.

## Module, scaling, and arbitrary lifts

The source is the full eight-dimensional module with monomial basis
`1,x,...,x^7`. First Hasse derivatives equal ordinary first derivatives,
so the rows really contain `q*x^(q-1)`; the factors divisible by three
must be retained. At the node zero, the value and derivative rows are the
first two coordinate rows. Subtracting their integral multiples from all
other rows gives an actual integral decomposition `I_2` plus the residual
six-by-six observer. No degree is removed beyond those two cleared pivots.

For normalized nodes `1,A,1+B`, with `A=3^d*a`, `B=3^r*b`, the minor
weight is exactly `3^(e*(sum selected degrees - number of derivative rows))`.
This is a factorization of each minor, not an assertion that dividing
columns by powers of three is unimodular. Every residual entry is divisible
by three when `e>=1`, and the degree-two derivative at normalized node one
has valuation exactly `e`. Hence there are precisely two initial unit
Smith factors and the next exponent is `e`.

All four cross-pair differences are units after removing the outer scale;
the two within-pair valuations are `d,r`. Translation and a unit affine
coordinate change preserve the full degree-below-eight source lattice
and its value/derivative observer up to integral invertible changes.
They therefore normalize any configuration with this distance pattern
to the stated form. Swapping the pair depths is equally lawful: on the
normalized variable, `u -> 1-u`, followed by reordering nodes, sends
`(0,1,A,1+B)` to `(0,1,-B,1-A)`. This interchanges `d,r` and changes
only unit factors. Signed or rational three-adic unit lifts cause no problem.

## Universal lower bounds and attained minima

For each polynomial term, the valuation is
`w*e+i*d+j*r+v_3(c)`. Coordinatewise dominance after the translation

```text
(w,i,j,c) -> (w,i,j,w+i+j+c)
```

proves a bound throughout `e,d,r>=1`, because the variables become
`e-1,d-1,r-1>=0`. The proof uses all terms of all minors, so it supplies
the lower bound on every determinantal ideal. Cancellation between terms
can only raise a valuation and cannot invalidate that direction.

For the upper bounds, each listed witness factors as its displayed powers
of `A,B`, its explicit integer constant, and factors that are units modulo
three. In the independent reconstruction, its distinguished monomial has
strictly smaller valuation than every other monomial for all positive
depths. Thus it cannot cancel, for any unit choices. The thirteen witnesses
attain every competing affine form in the six lower envelopes.

At a switching equality, two different witness minors may attain the same
minimum. Their terms are not being added: a determinantal ideal contains
each minor individually. Equality therefore introduces no cancellation
or unhandled unit-digit condition. The argument proves the universal
quantifiers directly; the producer's sampled specializations are controls.

The last witness gives determinant valuation `24e+4d+4r`, also equal to
four times the sum of the six pairwise node-difference valuations.
This is a useful mass check, and does not substitute for the intermediate
minor analysis.

## Sorted spectrum, equality boundaries, and precision

Put `s=min(d,r)`, `t=max(d,r)`,
`lambda=min(e,s+1)`, `mu=min(e,t+1)`, `nu=min(t,e+3s)`.
The six divisor valuations simplify exactly to

```text
D_1=e,
D_2=3e+lambda,
D_3=6e+s+mu,
D_4=11e+s+nu,
D_5=17e+4s+t,
D_6=24e+4s+4t.
```

Taking successive differences gives the reported eight-entry spectrum.
For the nontrivial sorting checks, `lambda<=mu<=e` and `nu>=s` imply

```text
(3e+s+mu-lambda)-(2e+lambda) >= s,
(5e+nu-mu)-(3e+s+mu-lambda) >= nu-s+lambda >0.
```

The next difference is
`e+3s+t-2nu+mu>=mu>0`, in both branches of `nu`; the final difference
is `e+2t-3s+nu>=e>0`. The earlier difference is also positive. These
inequalities include all switching equalities, including `d=r`,
`s+1=e`, `t+1=e`, and `t=e+3s`.

At the separate boundary `e=0`, unit-separated CRT gives the union of
the two pair lists `(0,0,d,3d)` and `(0,0,r,3r)`. Formula (1) in the
producer reduces to their sorted union because `nu=min(t,3s)`.
This boundary has its own proof; the positive-scale lower-bound
certificate is not improperly evaluated outside its domain.

The largest exponent is `7e+3t`, so the stated uniform precision loss
follows from invertible Smith coordinate changes over the three-adic
integers. It is sharp: the last Smith coordinate times `3^(N-1)` is
nonzero modulo `3^N` while its image is zero modulo `3^(N+7e+3t-1)`.
This concerns all eight source coefficients in the declared degree box.

## Independent finite symbolic replay

The following self-contained code was run in the workspace. It imports
neither the producer nor its polynomial arithmetic, constructs the
observer directly, and checks all 923 minors over the polynomial ring.
The coefficient/frontier test carries the unbounded-parameter proof;
it does not evaluate a finite list of depths.

```python
import sympy as s
from itertools import combinations
A, B = s.symbols('A B')
M = s.Matrix([
    [u**q if order == 0 else q*u**(q-1) for q in range(2, 8)]
    for u in (s.Integer(1), A, 1+B) for order in (0, 1)
])
frontiers = [
    [(1,0,0,0)],
    [(4,0,0,0),(3,1,0,1),(3,0,1,1)],
    [(7,1,0,0),(7,0,1,0),(6,1,1,1)],
    [(11,1,1,0),(12,4,0,0),(12,0,4,0)],
    [(17,4,1,0),(17,1,4,0)],
    [(24,4,4,0)],
]
def v3(c):
    c, v = abs(int(c)), 0
    if not c:
        raise RuntimeError('zero coefficient')
    while c % 3 == 0:
        c //= 3
        v += 1
    return v
def T(f):
    return f[:3] + (sum(f),)
terms = minors = 0
witnesses = [set() for _ in range(6)]
for k in range(1, 7):
    for rows in combinations(range(6), k):
        for cols in combinations(range(6), k):
            p = s.Poly(M.extract(rows, cols).det(method='domain-ge'), A, B)
            w = sum(c+2 for c in cols) - sum(r % 2 for r in rows)
            forms = [(w,i,j,v3(c)) for (i,j),c in p.terms() if c]
            for f in forms:
                if not any(all(x <= y for x,y in zip(T(g),T(f)))
                           for g in frontiers[k-1]):
                    raise RuntimeError(('dominance', rows, cols, f))
                if f in frontiers[k-1] and all(
                    h == f or (h[1] >= f[1] and h[2] >= f[2]
                               and sum(h[1:])-sum(f[1:]) > 0)
                    for h in forms
                ):
                    witnesses[k-1].add(f)
            terms += len(forms)
            minors += 1
    if witnesses[k-1] != set(frontiers[k-1]):
        raise RuntimeError(('witness', k))
    print('rank', k, 'complete; frontier witnesses', len(witnesses[k-1]))
if (minors, terms) != (923, 14501):
    raise RuntimeError(('universe', minors, terms))
print('PASS independent SymPy domain-ge audit:',
      minors, 'minors,', terms, 'terms, 13 attained forms')
```

The observed final output was

```text
PASS independent SymPy domain-ge audit: 923 minors, 14501 terms, 13 attained forms
PASS all13 displayed witness factorizations, up to permitted sign
```

The second line records a separate direct symbolic comparison of each
of the thirteen factored expressions in the report with its specified
row/column determinant. No mathematical correction was required.
