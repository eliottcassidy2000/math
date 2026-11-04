# Independent audit: the realized uniform-context rank

**Verdict: PROVED, with the scope qualifications below.** The independent
finite checker passes. It was written before inspecting the producer's
checker, using the proposed statement and inherited proof sources. A later
inspection of `../profiles/PROOF.md` and `../profiles/check.py` found their
proof, stated scopes, and implementation consistent with this audit.
Historical novelty has not been established.

## Inheritance and what changes

The closest proved mechanism is [THM-3181, tournament half-grid reciprocity
and repeated-join recurrence](../../../../01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md),
especially its parity theorem and half-grid recovery. That older theorem
expressly declines a sharp minimality statement on realized profiles.
[Session176's information proof](../../session176_tournament_information_20261104/information/PROOF.md)
adds actual uniform probes, their triangular transfer, contextual equivalence,
and the realized lower bound `floor(n/3)+1`. The present construction closes
the gap between that lower bound and the older parity upper bound.

The hostile is the inherited pair of order-five profiles with `H=15` whose
common uniform wrapper produces 123 and 135. The boundary to retain is the
two-vertex reversed edge: its Hamilton count is 1, so the auxiliary
`1+2^(m-2)` formula applies only for `m>=3`. The older parity sidecar is the
key information missed by using the formal `n`-coordinate space alone.

Targeted searches in the canon and result notes for path-colour/path-cover
rank, realized rank, dimensions involving `ceil`, and the relevant correction
terms found no earlier exact-dimension statement. This is an inheritance
check, not an exhaustive literature search or a priority claim.

## 1. Exact linear and affine spaces

For fixed `n>=1`, write `d=ceil(n/2)`, and let `q_T(k)` count ordered spanning
directed path covers. Define

    Q_T(t) = sum_(k=1)^n q_T(k) binom(t,k).

The proposed conclusion is correct:

    span_Q{Q_T : |T|=n}
      = span_Q{t^n,t^(n-2),...,t^(n-2d+2)},

and the affine hull consists of the polynomials in this space whose
coefficient of `t^n` is 1. Its linear dimension is `d`, and its affine
dimension is `d-1`. The constant term is absent even when n is even.

For the upper bound, THM-3181 gives pure order parity. There is no constant
term, since every binomial `binom(t,k)` in the definition vanishes at zero.
The leading coefficient is 1, because `q_T(n)=n!`.

For the lower bound, let `B_m` be the transitive tournament on
`0,...,m-1`, except that its extreme edge is reversed to `m-1 -> 0`, with
`m>=3`. Every Hamilton path omitting that reversed edge is the unique
increasing word. Every path using it is obtained by choosing which of the
`m-2` interior vertices precede `m-1`; both the chosen vertices and the
remaining vertices have forced increasing orders. Thus

    H(B_m) = 1 + 2^(m-2),        m>=3.

To compute Q, colour the vertices by one of t ordered colours and count a
Hamilton path independently on each nonempty colour class. Compared with the
transitive tournament, nothing changes when the two extremes have different
colours. If they have the same colour and exactly j interior vertices share
it, the extra count is `2^j` for `j>=1`, and zero for `j=0`. Consequently

    Q_(B_m)(t)
      = t^m + t sum_(j=1)^(m-2) binom(m-2,j) 2^j (t-1)^(m-2-j)
      = t^m + t[(t+1)^(m-2) - (t-1)^(m-2)].

The calculation holds for all positive integers t, hence as a polynomial
identity. It also respects `t=0`. No odd-cycle formula or unproved
realizability premise is needed.

For each `1<=j<=d-1`, take the actual order-n tournament obtained by appending
`n-(2j+1)` singleton vertices in an ordered sum after `B_(2j+1)`. The
path-colour polynomial is multiplied by that many factors t. Its lowest
nonzero power is `t^(n-2j)`, with coefficient 2. Add the transitive order-n
tournament, with polynomial `t^n`. These d distinct lowest powers prove
linear independence. Their common leading coefficient and linear
independence also give affine independence. This establishes both equalities
of spaces and both exact dimensions.

For n=1 and n=2 the family has only the transitive polynomial; dimensions
are 1 and 0. Order zero is not in the statement. The witnesses need not be
strong after singleton padding, so this proof does not establish the same
rank in the strong-tournament subclass.

## 2. Recovering the profile using actual contexts

Session176's alternating transitive exterior `P_r` satisfies

    h_r = H(P_r[T]) = sum_(k=1)^(r+1) c_r(k) q_T(k),
    c_r(r+1)=r!.

Its proof includes endpoint legality and counts actual tournament contexts.
Triangular inversion of `h_0,...,h_(s-1)` recovers `q_T(1),...,q_T(s)`.
Those determine `Q_T(1),...,Q_T(s)` by the binomial formula.

Taking `s=d` recovers every polynomial in the d-dimensional linear space:
divide by `t^epsilon`, where epsilon is 1 for odd n and 2 for even n, and
interpolate at the distinct squares `1^2,...,d^2`.

For actual fixed-order profiles the leading coefficient is already 1.
Subtract `t^n` and take only `s=d-1`. Interpolation now recovers a polynomial
of degree at most `d-2` in `t^2`, after division by `t^epsilon`. Thus the
first `d-1` actual observations `h_0,...,h_(d-2)` suffice with affine
decoding. For d=1 there are no observations; the known order gives `Q=t^n`.

The rank bound is sharp for linear reconstruction on the full linear span,
and the affine bound is sharp for affine reconstruction on the affine hull.
If affine encoder and decoder compose to the identity on all actual
profiles, that identity extends to their affine hull, so the intermediate
coordinate space must have dimension at least `d-1`.

These are **not** unrestricted query, bit, or nonlinear encoding lower
bounds. In particular, a single suitable linear scalar can distinguish a
finite set of actual profiles if an arbitrary nonlinear decoder is allowed.
Nor does affine-hull equality say every monic parity polynomial is realizable
by a tournament. Nonuniform attachments can distinguish equal q profiles,
as the inherited order-two hostile shows.

## 3. Independent finite verification

`check.py` uses literal vertex permutations and mandatory cuts at every
backward adjacency to compute q; it does not use endpoint-cover recurrences.
Hamilton counts of the wrapped tournaments are checked separately by an
endpoint subset dynamic program. Skeleton rows are counted by literal
exterior permutations with inserted hole slots. Rational rank and
interpolation use exact `Fraction` arithmetic.

The explicit universe is:

- all 1,099 labelled tournaments of orders 1 through 5, with no filtering;
- every proposed lower-bound witness at orders 1 through 8;
- actual shortened uniform probes for every tournament in those two sets;
- skeleton rows for r=0,1,2,3, including their factorial diagonals;
- symbolic triangular supports at orders 1 through 128, distinguished from
  graph enumeration;
- the equal-H order-five hostile giving 123 versus 135, the nonuniform
  order-two hostile giving 3 versus 1, and the `m=2` auxiliary-formula
  boundary.

The census ranks are `1,1,2,2,3`; affine ranks are `0,0,1,1,2`.
Every full profile is recovered from the stated actual observations.
Normal and optimized runs use explicit validation rather than `assert` and
produce byte-identical receipts. The source is hashed at startup and checked
again before output; no expected receipt is regenerated inside the checker.

From the repository root:

```text
python3 05-knowledge/results/proof_principles_20261103/rank_audit/check.py --output /temporary/path/audit.json
python3 -O 05-knowledge/results/proof_principles_20261103/rank_audit/check.py --output /temporary/path/audit.optimized.json
```

Compare those files with the committed `receipt.json` and
`receipt.optimized.json`. The universal dimension and recovery statements
rest on the proofs above, not on extrapolating these finite checks.

## 4. Fresh gate-flow consequence

[The gate-flow packet, CF2 and CF6--CF7](../../circuit_flow_tournaments_20261103/PROOFS.md)
proves that the readouts of all L/U words of depth at most r span exactly
the evaluations `Q(1),...,Q(r+1)`. Its triangular basis is already given by
the prefix contexts `1,L,...,L^r`; the newest evaluation has nonzero
coefficient `product_(j=1)^r j(j+1)(2j+1)`.

Restricting these functionals to the now-known realized span V_n gives
rank `min(r+1,d)`. Indeed, after the same nonzero row scaling, their matrix
is a Vandermonde evaluation matrix on the distinct squares
`1^2,...,(r+1)^2`, with d polynomial coefficients. On differences of monic
profiles the leading coefficient vanishes, leaving d-1 coefficients, so
the affine variation rank is `min(r+1,d-1)`.

It follows that depth `d-1` is sharp for linear recovery on V_n. Depth
`d-2` is sharp for affine recovery on the fixed-order affine hull when
`d>=2`; when `d=1` the profile is already known without observations.
These depth statements refer to the collection of readouts up to that
depth, equivalently the displayed prefix probes, rather than a single
terminal Hamilton count. They are confined to the L/U language and to
linear or affine recovery on the specified spaces.

This corollary passes the algebraic audit independently of its finite
operator tests: it is a restriction of the inherited triangular observer
span, followed by a full-rank Vandermonde evaluation argument. The uniform
probes in section2 and these four-exterior-vertex gate probes are different
context families; their shared information rank does not identify them as
the same graph construction.
