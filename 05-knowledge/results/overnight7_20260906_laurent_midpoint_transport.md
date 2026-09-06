# Three grouped responses for the entire A=2, B=3 carry progression

**Status: PROVED for the all-parameter identities and degree reduction;
FINITE-EXACT for named grouped signs; universal signed transport OPEN.**
[Independent analytic and exact audit: PASS](overnight7_20260906_laurent_independent_audit.md). This is a structural specialization of
the incoming missed-midpoint mechanism, not a new claim of its general
path inclusion or virtual SOS theorem. No theorem ID is claimed.

The companion [endpoint-27 certificate](overnight7_20260906_laurent_quartic_carry.md)
closes one five-first-channel family. The purpose of this note is to identify
an object that can be studied uniformly in channel count, rather than add
another list of fixed-endpoint certificates.

## 1. Inheritance, the full source and the retained board

Before deriving I read incoming `6f450d8cb`'s
`05-knowledge/results/nc2_hadamard_transport_overnight_hexagon_sep05.md`
directly from `origin/main`. It proves: the full path-factor virtual row is
negative at every first root; every actual coefficient dominates its virtual
coefficient; and actual-minus-virtual is negative at the first roots of the
endpoint-15 family. Its universal signed transport is explicitly OPEN.

The first sign uses
[THM-4440, signed duplication SOS](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md)
on a single ordinary real-rooted carrier. Its coefficient condition must be
preserved. Individual real-rootedness of another carrier is not enough.
The corrected near misses remain deletion of negative-index auxiliary
support, confusing parameter doubling with auxiliary squaring, and inferring
negative-phase inequalities from coefficientwise domination.

The board is: parity midpoint sections; skipped levels; full Laurent path
sequences; contiguous Euler operators; a common zero-coefficient constraint;
quotient characteristic coefficients. The exact source here is

```text
h>=1, x>=1 integers, y=3h, g=x+3h+1,
a=6h+3,
support (-a,2g-a,3g-a).
```

The first-return interpretation additionally requires `gcd(g,a)=1`.
The path identities below do not use that gcd condition.

Define complete binomial sequences and their Laurent generating functions:

```text
O_j=binom(g,2j+1),         j>=0,
E_j=binom(g,2j),           j>=0,
B_j=binom(x+y-2j,x+j),     -x<=j<=h,
C_j=binom(x+y-1-2j,x+j),   -x<=j<=h-1,
D_j=binom(x+y-2-2j,x+j),   -x<=j<=h-1.
```

All sequences are zero outside their actual count domains. Products such
as `B^2` mean ordinary Laurent multiplication. `star` means coefficientwise
Hadamard multiplication, retaining the exponent label. The first row is

```text
P_raw=O star B.
```

Every negative index of `B,C,D` is retained before any multiplication.

## 2. Exact two-class and three-class midpoint partitions

The alpha path factor has the elementary parity identity

```text
A_double = O^2+t^-1 E^2.                           (1)
```

Indeed `[t^j]A_double=binom(2g,2+2j)`. Splitting a length-`2g` path at
its unique level-`g` vertex gives odd/odd midpoint heights, counted by
`[t^j]O^2`, or even/even heights, counted by `[t^(j+1)]E^2`.
This is also the even-part identity obtained by squaring
`(1+s)^g=E(s^2)+s O(s^2)`.

The beta factor has the useful exact identity

```text
B_double = B^2+2t C D,                             (2)
[t^j]B_double=binom(2x+2y-2j,2x+j).
```

Here is the complete path proof. The doubled path ends at
`(2y-3j,2x+j)`. The functional `X+3Y` runs strictly upward from zero
to twice `H=y+3x`. A path either hits level `H` at a vertex, or crosses
it by a vertical edge from level `H-1` or `H-2`. It cannot have two such
events. The hit class is counted by `B^2`, since all integral vertices on
the selected level are `(y-3l,x+l)`.

For a jump from `H-a`, where `a=1,2`, the prefix ends at
`(y-a-3l,x+l)`. After the forced vertical step the suffix count is that
of a path to `(y+a-3k,x-1+k)`, with `k=j-l`. Thus the two jump classes
are the convolutions

```text
beta(x,y-1) beta(x-1,y+1),
beta(x,y-2) beta(x-1,y+2).
```

Their full Laurent generating functions satisfy

```text
beta(x-1,y+1)=t D,
beta(x-1,y+2)=t C.
```

Both jump classes therefore equal `t C D`, giving `(2)`. Equality of
these grouped classes is not a pairing of individual negative frequency
products; it is an exact support-preserving reindexing of path counts.

## 3. Actual-minus-virtual has three grouped responses

The actual anchored doubled row and the virtual row are exactly

```text
Q_raw=(O^2+t^-1 E^2) star (B^2+2t C D),
V=O^2 star B^2.
```

Consequently

```text
Q_raw-V = G1+G2+G3,
G1=(t^-1 E^2) star B^2,
G2=2 O^2 star(t C D),
G3=2(t^-1 E^2) star(t C D).                        (3)
```

Every coefficient in all three groups is nonnegative. The lower carry
survives in exponent `-1`. Literal multinomial coefficients are

```text
[t^j]Q_raw=(2g)!/[(2x+j)!(2y-3j)!(2+2j)!],
            -1<=j<=2h.
```

For the full incoming factors, `V(rho)<0` at every first root `rho<0`.
Thus `G1(rho),G2(rho),G3(rho)<0` would be a sufficient uniform transport
theorem. It is stronger than necessary: only the combined signed response
is needed, and coefficientwise nonnegativity supplies none of these value
inequalities on its own. Universal individual-group or combined-group
negativity remains OPEN in this note.

The exact named bank uses

```text
(h,g)=(1,5),(1,7),(1,11),
      (2,8),(2,11),(2,16),
      (3,11),(3,13),(3,17),
      (4,14),(4,16),(4,19).
```

All 30 first roots and all three grouped responses at each root are covered.
Exact polynomial inversion retains `t^-1`, and rational interval Horner
evaluation of the remainders proves every tested group strictly negative.
The intervals are refined until the whole value interval is negative;
neither midpoint sampling nor floating roots certify these signs.

This is a small mechanism-directed bank, not evidence promoted to an
all-`h` claim. The family proof in the companion uses its own exact
characteristic certificate and does not depend on these signs.

## 4. Same-source contiguous operators and the first failed shortcut

Let `Theta=t d/dt`, and put `N=x+y`. Direct coefficient identities give

```text
g E=[1+(g-1)t] O+2t(1-t)O',
(N-2Theta) C=(y-3Theta)B,
(N-1-2Theta)D=(y-1-3Theta)C.                        (4)
```

The first follows either from differentiating the even/odd binomial split
or from its two coefficient ratios. For the other two, each common index
has `N_j=x+y-2j`, and

```text
C_j/B_j=(y-3j)/N_j,
D_j/C_j=(y-1-3j)/(N_j-1)
```

where the equations are read after cross-multiplying at vanishing boundary
coefficients. On the full support of `B`, `N_j>=x+h>=2`; the denominators
`N_j` and `N_j-1` never vanish. Thus `(4)` gives well-defined diagonal
Euler inverses on this support. It introduces no deleted endpoint.

These are operators on the same source factors, but they do not separately
preserve the root condition required by THM-4440. The earliest named control
already proves the failure. At `h=1,g=5,x=1,y=3`, the primitive source is
`(-9,1,6)` and `rho=-2` is the first-row root. Exactly

```text
(O star B)(rho)=0,
(O star C)(rho)=15,
(O star D)(rho)=10,
(E star B)(rho)=-16.                               (5)
```

Thus applying separate virtual SOS statements after replacing a factor by
`E,C,D` loses the hypothesis. First failed implication: a contiguous
operation was treated as preserving the particular vanishing coefficient.
The strongest surviving operation is `(4)` together with the common
constraint in `(5)`. A useful next target is a coupled quadratic response
of the same real-rooted carrier under that constraint, or a sign certificate
for the whole three-group sum. Independent signs or real-rootedness of
replacement carriers do not suffice. This is the recorded stopping point;
no universal transport is claimed.

## 5. General characteristic normalization without another endpoint scan

For the entire `a=6h+3` progression, first and doubled canonical fibers are

```text
(g-3h-1+j,3h-3j,1+2j),           j=0,...,h,
(2g-6h-3+j,6h+3-3j,2j),         j=0,...,2h+1.
```

With `X=alpha^(g-3h-1) beta^(3h) gamma`, define

```text
L=(g)_(2h+1)/(3h)!,              K=(2g)_(4h+2),
P_hg(t)=sum_(j=0)^h [(3h)!/((3h-3j)!(1+2j)!)]
                       *(g-2h-1)_(h-j)*t^j,
Qbar_hg(t)=sum_(j=0)^(2h+1) (2g-4h-2)_(2h+1-j)*t^j
                       /[(6h+3-3j)!(2j)!].
```

Then `CT(f^g)=X L P_hg`, and
`CT(f^(2g))=X^2 K t^-1 Qbar_hg`. Write `p_hg` for the monic first row.
The inverse-carry denominator cancels uniformly:

```text
Qbar_hg(0)/p_hg(0)
 = [2^(h+1)(3h)! / ((2h+1)!(6h+3)!)]
   *(g-3h-1)*product_(r=0)^(h-1) (2g-4h-3-2r).
```

This follows by separating the even and odd factors in the `2h+1`-term
falling factorial. With weight one assigned to `g,t`, monic reduction
preserves weighted degree. The inverse term has weight at most `2h`, as
do all positive-index terms after division by `t`. Therefore

```text
R=t^-1 Qbar_hg mod p_hg=sum_(j=0)^(h-1) R_j(g)t^j,
deg_g R_j<=2h-j,
deg_g [z^(h-k)]det(zI-R(C))<=2hk.                   (6)
```

The last inequality follows from Newton sums of the monic companion
matrix. Formula `(6)` is an all-`h` algebraic degree bound and audit
interface. It does not imply that the characteristic coefficients are
positive for all `h`. Shifted coefficient positivity at `g=s+3h+2`
remains a precise possible stronger conjecture. No `h>=5` certificate
or census was computed in this task.

## 6. Verification, connection contract and files

The companion imports no repository producer. In the 12 named controls it
checks `(1)--(4)` over the complete Laurent supports, then reconstructs
the beta partition by an independent dynamic program on 164 path endpoints.
Each path is classified at its unique hit or skipped-level crossing.
The actual row is compared with literal multinomial coefficients. All 90
group/root signs are certified by exact interval arithmetic. The source
returns **767 explicit gates**; normal and optimized outputs match.

```text
python -B 04-computation/overnight7_20260906_laurent_midpoint_transport.py
python -B -O 04-computation/overnight7_20260906_laurent_midpoint_transport.py
```

```text
source 890abb5a5bf598aebceacff6b228a0c70b76832d8d4d9812d7975133438f0649
output 1a41c93f466a3e21d65c95482b55e489fdf9737a44c5aee3ff5e282fe60ce996
semantic 4f855b7007e6f793fb65d33f076990a7b5db226838eaf52700befc6f878b59dd
```

Connection contract: complete first path factors map to a three-group
partition of the actual doubled row. The map preserves all coefficients,
negative indices and the actual first-anchor normalization. Passing to the
three groups forgets the finer missed-midpoint pairing but loses no signed
coefficient data. The necessary sidecar is the original common condition
`(O star B)(rho)=0`; the contiguous-carrier hostile proves that it cannot be
silently transported. Formula `(5)` is the cheapest failed shortcut,
and the exact 90-sign bank is the positive test for the surviving object.
