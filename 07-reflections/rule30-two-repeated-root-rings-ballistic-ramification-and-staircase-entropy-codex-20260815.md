> **STATUS:** session synthesis and idea provenance.  Current truth is proved
> [THM-3476](../01-canon/theorems/THM-3476-rule30-depth-four-transverse-jet-barrier-and-slack-pascal-atlas.md),
> [THM-3480](../01-canon/theorems/THM-3480-rule30-staircase-transducer-entropy-and-nonrectangular-macroblock-compiler.md),
> and [THM-3481](../01-canon/theorems/THM-3481-rule30-cyclic-arc-norm-rank-and-marked-innovation-spectrum.md).
> None of the three Rule 30 prizes is claimed here.

# Rule 30 next targets: two repeated-root rings, ballistic ramification, and staircase entropy

This session began with three debts left by THM-3471:

1. how much of the Green-slack marker must be retained;
2. what a moving terminal arc actually forgets; and
3. whether the rectangular macroblock halo is a real state requirement.

Each debt now has an exact answer.  The answers point in different directions,
but they share one useful warning: a full profile can be algebraically hard or
even lossless while one calibrated coefficient remains uncontrolled.

## 1. The slack tower is a genuine local coordinate, but an unbounded one

Package the physical radial source before transport as

```text
S(z,X)=sum_(u,d) alpha_u(d)z^uX^d.
```

The unmarked Motzkin compiler evaluates `X` at the small root `W` of

```text
P(z,X)=X^2+(1+z)X+z^2.
```

The complete evaluation kernel is exactly `(P)`.  After deforming the slack
marker to `q=1+epsilon`, one `P`-factor becomes exactly one `epsilon`-factor.
Therefore jets through order `m` have kernel `(P^(m+1))`, and the infinite
jet tower is faithful.

Faithful is not bounded.  At physical source depth four, for every power of
two `M>=32`,

```text
[z^(5M/4+10)]R_4(z,q)=q^(M/4-7)(1+q)^M.
```

The first live transverse layer is `M`, linear in the target time.  This is a
real Rule 30 strip, not a formal `P^M` packet.

The closest older analogue is THM-2043: a finite Hasse chart can be complete
for local residues while remaining blind to magnitude and owner height.  Here
the first `2^m` jets recover exactly slack residues modulo `2^m`; higher blocks
recursively recover quotient carries.  The full Pascal tower is a binary
carry tree.

## 2. Ballistic transport ramifies that filtration by two

The slack tower is itself redundant.  At fixed source depth `u` and target
`t`, let `B_(u,t)(X)` record the source distances which survive the Green
kernel.  If `R=t-u-1`, the target slack polynomial is

```text
C_(u,t)(q)=q^R B_(u,t)(q^(-2)).
```

Because `q^(-2)+1=q^(-2)(q+1)^2`,

```text
ord_(q+1) C = 2 ord_(X+1) B.
```

More strongly, after reversing `B` into `H`,

```text
C(q)=q^delta H(q^2),
J_(2j)=h_j,
J_(2j+1)=delta h_j.
```

So half the slack jets are forced duplicates or forced zeros.  The smaller
prize-facing carrier is the Green-selected **distance packet**, with
source-depth parity retained.  The depth-four hostile is simply

```text
B=X^6(1+X)^(M/2)  ->  C=q^(M/4-7)(1+q)^M.
```

This is the same operation-level lesson as factorial boundary Stokes: a
labelled transverse layer survives a scalar cancellation.  It is also the
same circuit lesson as the three Berggren siblings: forgetting a label can
turn reconstructive data into zero.  No factorial or Berggren theorem is
being imported; the shared map is labelled boundary transport.

## 3. A terminal arc is a cyclic norm, not a smoothing heuristic

On a phase cycle of size `p=2^d`, summing a current over `m` consecutive
phases is multiplication by

```text
N_m(S)=1+S+...+S^(m-1)
```

in the repeated-root ring `F_2[S]/((S+1)^p)`.  If `m=2^a u`, `u` odd, then

```text
N_m=(S+1)^(2^a-1) times a unit,
rank N_m=max(p-2^a+1,0).
```

Odd arcs are invertible: they lose no phase-current information.  At odd
innovative Rule 30 depths, the entire physical terminal profile therefore has
maximal ANF degree; beyond the one-variable boundary it has no zero Walsh
coefficient.

Yet the center is still only

```text
c_k=A_(k,k)(-k).
```

Rotating the mark preserves every unpointed cyclic invariant and realizes
both bit values.  The obstruction is not bulk smoothing; it is calibrated
point evaluation.  This sharpens the density target from “show the current is
random-looking” to “control one moving open-line product at the physical
basepoint.”

The two local rings should not be conflated:

| ring | uniformizer | exact loss | missing datum |
|---|---|---|---|
| Motzkin source | `P(z,X)` | finite slack jets | unbounded distance carry |
| phase cycle | `S+1` | `2^nu2(m)-1` Hasse moments | physical phase mark |

Ballistic transport links their filtrations with ramification index two, but
it does not identify source distance with phase.

## 4. The rectangular halo was not the minimal macro state

One Rule 30 row has a three-state right-to-left transducer:

```text
A=00,  B=01,  C=1*.
```

Cascading `h` copies turns a two-sided rectangular halo into one ordered
staircase boundary.  From the physical zero right exterior, only a reachable
language `S_h` occurs.  A 5,873-vertex exact graph certificate proves

```text
|S_h| <= 41(117/64)^h,
entropy(S_h)<7/8 bit per row.
```

The resulting charged chunk table preserves exact Rule 30 evolution and gives

```text
Q(n)=(7/2+o(1))n^2/log^2(n)
```

main lookups, versus tariff `8` for the raw rectangular one-word cone model.
The total word-RAM upper bound remains `O(n^2/log^2 n)`: this improves the
carrier and leading lookup geometry, not the exponent and not Prize 3.

The surprising analogy is with transfer matrices and tree-width.  The right
quantity is not area of the eliminated region, but entropy of its reachable
boundary response.  The exact integer super-eigenvector is the analogue of a
certified transfer-matrix pressure bound.

## 5. The next three decisive probes

### Prize 1: repeated-root depth of the inward completion

The boundary depths `u=0,1,2` form a null circuit.  The physical depth-four
family has unbounded distance-packet order.  The next useful object is

```text
minimum (X+1)-order surviving after summing u>=3
within each source-parity channel.
```

A bounded value would expose a finite Cartier carrier; an unbounded value
with a proved lowest surviving grade would be a plausible nonperiodicity
route.  Scalar nonperiodic summands are no longer enough.

### Prize 2: pointed diagonal product correlation

Odd terminal integration is lossless, and innovation profiles are spectrally
full.  Therefore another phase-average or marginal density is unlikely to be
the missing theorem.  The target is discrepancy of the marked arc
`[-k,0)`, jointly across changing depths and calibrated to the seed boundary.
The cheapest hostile must preserve the full cyclic profile while moving its
mark.

### Prize 3: lower staircase entropy or a smaller carrier

The exact upper bound is now below `7/8`; finite Myhill--Nerode audits find no
state merging through height 13.  The fork is clean:

1. prove positive entropy for `S_h`, yielding a restricted diagonal-lookup
   obstruction; or
2. find a subexponential sufficient boundary quotient, which would point to a
   stronger time compression.

Universal block rank does not decide this fork.  The staircase compiler is
already a hostile: full fresh-bit permutivity coexists with a compressed
reachable boundary.

## 6. One reusable research rule

Before calling an average, window, or eliminated halo “simple,” factor its
operator in the relevant local ring and retain the boundary response state.
Then ask separately:

1. what is the exact kernel or repeated-root order;
2. whether the observer commutes with that filtration; and
3. whether a marked coefficient can bypass the full-profile obstruction.

That sequence of questions produced all three advances here.  It also explains
why none, by itself, settles a Rule 30 prize.
