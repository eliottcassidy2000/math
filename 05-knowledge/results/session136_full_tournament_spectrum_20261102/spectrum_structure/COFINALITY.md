# Multiplicative gap amplification: relative coverage is not completeness

**PROVED elementary consequence of inherited SCC multiplication.** This is a
secondary framework, not a new construction of any missing target and not an
external-priority claim. Full spectrum completeness remains **OPEN**.

## A two-value amplifier

Suppose `1<a<b` are attained Hamiltonian-path counts. Let m0 be the least
positive integer with `b^m0 >= a^(m0+1)`. Then, for every real `x>=a^m0`,

```
S intersects [x, (b/a)x).
```

**Proof.** Pick m with `a^m<=x<a^(m+1)`; then m>=m0 and
`a^(m+1)<=b^m`. The actual order-join counts

```
a^m, a^(m-1)b, ..., ab^(m-1), b^m
```

have constant adjacent ratio b/a and bracket x. The first one at least x is
strictly below (b/a)x. Every value has an actual witness as a join of the
chosen a/b witnesses. QED.

The recipe uses unrelated joins, not a nested sequence of induced tournament
extensions. Thus it does not contradict Session134's zero-density theorem
for a fixed nested ray.

## Relative gaps vanish already in a very sparse subfamily

Let M be the submonoid `{3^u 5^v:u,v>=0}`, using C3 and the strong four-vertex
H=5 tournament. Then consecutive elements of M, and therefore consecutive
attained values, have ratio tending to one.

For any epsilon>0, irrationality of `ln3/ln5` and the elementary pigeonhole
approximation of its multiples give positive integers u,v for which
`0<|u ln3-v ln5|<ln(1+epsilon)`. The two distinct attained numbers `3^u`
and `5^v`, in increasing order a<b, have b/a<1+epsilon. The amplifier gives
an eventual (1+epsilon)-mesh. As epsilon is arbitrary, the ratio tends to one.

For completeness, the approximation follows by placing N+1 residues
`0,theta,...,N theta (mod1)` into N equal subintervals. A pair yields
`|q theta-p|<1/N`; choose N arbitrarily large. For sufficiently small error
both q,p are positive. Irrationality excludes equality, since equality would
say a positive power of three equals a positive power of five.

Yet `#(M intersect[1,X]) <= (1+ln X/ln3)(1+ln X/ln5)`, so M has natural density
zero and contains no prime other than 3 or 5. Therefore **relative gap
vanishing is much weaker than finite complement or full odd completeness**.
One can have vanishing relative gaps and infinitely many, indeed almost all,
odd numbers missing. This is a structural hostile to “cofinal and
multiplicatively dense, therefore complete.”

## An explicit arithmetic alternative

Session134's actual families realize `a=2^k+1` and `b=2^k+3` for every k>=6:
`05-knowledge/results/session134_constructive_spectrum_frontier_20261101/descent/THEORY.md`.
Thus the amplifier yields eventual relative mesh `1+2/(2^k+1)` with completely
explicit graph recipes. For example a=65,b=67 has m0=138, hence every
`x>=65^138` lies below an attained value smaller than `(67/65)x`. This huge
threshold is an exact illustration, not a useful numerical coverage claim.

Connection contract: two actual counts with graph witnesses -> a finite
geometric mesh of products -> overlapping scale intervals. H and witness
recipes are preserved. Exact target identity, induced containment and graph
strongness are not. The cheapest decisive hostile is M=<3,5>: eventual
relative approximation holds but any new odd prime is absent from M.
