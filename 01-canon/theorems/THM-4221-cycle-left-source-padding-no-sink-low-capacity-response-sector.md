---
id: THM-4221
title: "Cycle-left source-padding no-sink low-capacity response sector"
status: >
  PROVED exact cycle-left source-padding response identity, marked-Hamilton-
  arc avoidance injection, joint endpoint-chirality sidecar, and an all-order
  no-sink low-capacity sector with strict clean sub-sector; PROVED implication
  from the open three-eighths endpoint gap and CONDITIONAL HYP-9081 iterated
  source-padding consequence; VERIFIED-EXACT on all 1,099 labelled tournaments
  through order five. The universal no-sink three-eighths gap, universal
  response monotonicity, and general (OS+) remain OPEN.
source: codex-endpoint-gap-session-20260826
depends_on:
  - THM-4208-cycle-prefix-arbitrary-context-recurrence-endpoint-energy-and-eventual-positivity
  - THM-4219-no-sink-endpoint-energy-floor-and-near-ordinal-sharpness
related:
  - HYP-9081-strong-tournament-five-copy-endpoint-energy-inequality
script: 04-computation/tournament_cycle_left_source_padding_response_sector_thm4221.py
output: 05-knowledge/results/tournament_cycle_left_source_padding_response_sector_thm4221.out
script_sha256: 206f2e0c58e1e4a3de50f0d71f80477e0cf4b9949c0dad5137f98d1d9a1da0ac
output_sha256: 488ec22c1eb87b2a37ebdc2bab18bb9ec9769bcdb99cba16c6e2dde23e51f2ae
hash_basis: raw LF bytes
primary_audit: >
  PASS. The exact exposed-path engine checks source-padding transfer and both
  forms of the response increment on every 1,099 labelled tournament through
  order five. It separately rebuilds both R_+ values from literal ordinal
  children, enumerates every marked Hamilton-path arc, and checks the improved
  endpoint sidecar and sector bounds on all 738 no-sink rows. Normal and
  optimized streams byte-match the frozen output.
referee_audit: >
  ACCEPT. An independent derivation reproduced the source-padding coordinate
  table, the response identity, marked-arc injection, low-capacity root, and
  conditional constants. The final 45H^2 endpoint improvement strengthens
  the referee's original 36H^2 bound by retaining the no-sink endpoint-
  chirality sidecar.
---

# THM-4221 -- cycle-left source-padding no-sink low-capacity response sector

**PROVED + VERIFIED-EXACT; CONDITIONAL HYP-9081 CONSEQUENCE.**

THM-4208 leaves open whether `R_+(C3,Y)>0` for every no-sink right factor.
This theorem studies a different but adjacent operation: add one universal
source to the right factor and ask how much the cycle-left response changes.
The answer is an exact endpoint-energy identity. A marked-Hamilton-arc
injection and a joint start/end chirality bound then give a new all-order
positive sector. The full no-sink statement remains conditional on the open
three-eighths endpoint gap.

## 1. Notation and statements

Let `C` be a nonempty tournament of order `n`. Retain THM-4208's Hamilton
count, capacity mass, rooted state, capacity degree, and incoming capacity
mass, and abbreviate

```text
h=H(C),                         z=W(C),
v_i=V_i^0+V_i^1,               t_i=Start_i(C),
e_i=End_i(C),                  d_i=capacity degree of i,
r_i=incoming capacity mass at i,
m=sum_i v_i^2,                 tau=sum_i t_i^2,
Delta_V=(z+h)(z+3h)-3m.                              (1)
```

Thus THM-4208's rooted collapse gives `v_i=r_i+t_i`, while
`sum_i v_i=z+h` and `sum_i t_i=sum_i e_i=h`. Define the source-padding
increment

```text
I(C)=R_+(C3,P_1 triangleright C)-R_+(C3,C)             (2)
```

and the auxiliary marked-avoidance sum

```text
A(C)=sum_i (z-d_i)(27v_i-9t_i).                        (3)
```

> **Theorem 1 (exact source-padding response).** Every nonempty tournament
> satisfies
>
> ```text
> I(C)=324z^2+1188zh+744h^2-945m-15tau
>      -27sum_i d_i v_i+9sum_i d_i t_i                (4)
>
>     =-18z^2-90zh-201h^2+315Delta_V-15tau+A(C).      (5)
> ```

For a no-sink tournament of order `n>=4`, put `x=z/h` and

```text
rho_n=
 (9n-57+sqrt(2601n^2-3402n-10391))/12.                (6)
```

> **Theorem 2 (all-order low-capacity sector).** If `C` has no sink, then
>
> ```text
> n>=4 and W(C)/H(C)<=rho_n
>       implies I(C)>=1480H(C)^2.                      (7)
> ```
>
> The simpler sufficient condition
>
> ```text
> n>=4 and W(C)<=(5n-11)H(C)                          (8)
> ```
>
> gives the strict inequality `I(C)>1480H(C)^2`. At order three, no-sink
> forces `C=C3`, and
>
> ```text
> I(C3)=13320=1480H(C3)^2.                            (9)
> ```

The radical in `(6)` is the exact positive root supplied by the scalar
bounds in this proof. It is not asserted to be the actual boundary of
source-padding monotonicity.

> **Theorem 3 (proved gap implication).** If a no-sink tournament `C`
> satisfies
>
> ```text
> Delta_V(C)>=3(W(C)+2H(C))^2/8,                      (10)
> ```
>
> then `I(C)>=1480H(C)^2`, with equality exactly when `C=C3`.

HYP-9081 would imply `(10)` for every no-sink tournament. Consequently,
conditional on HYP-9081, for every no-sink `C` and integer `a>=1`,

```text
R_+(C3,P_a triangleright C)-R_+(C3,C)
 >=1480aH(C)^2,                                       (11)
```

with equality exactly at `a=1,C=C3`.

## 2. Proof of the exact increment

Write `C^+=P_1 triangleright C` and call its new universal source `q`.
The exact ordinal transfer in THM-4184, together with THM-4208's rooted fan
collapse, gives

```text
H(C^+)=h,                         W(C^+)=2(z+h),
V_q(C^+)=(0,h),                   V_i(C^+)=(v_i,v_i),
r_q(C^+)=0,                       r_i(C^+)=2v_i,
d_q(C^+)=z+2h,                    d_i(C^+)=d_i+r_i+2t_i.
                                                               (12)
```

For example, the cross capacity at old vertex `i` is `r_i+2t_i`; adding
its inherited incoming mass `r_i` gives `r_i(C^+)=2(r_i+t_i)=2v_i`.

THM-4208 equation (43), with `y_i=V_i^0+2V_i^1`, is

```text
R_+(C3,Y)
 =108H(Y)W(Y)+120H(Y)^2
  +18sum_i y_i(W(Y)-d_i-4r_i)
  +9(3W(Y)+4H(Y))^2-84sum_i y_i^2.                    (13)
```

In the parent `C`,

```text
y_i=(3v_i+t_i)/2,                 r_i=v_i-t_i.         (14)
```

In `C^+`, the new-source coordinate is `y_q=2h` and every old coordinate
is `y_i^+=3v_i`. Substitute `(12)--(14)` into `(13)` and subtract. The
global terms and the new-source term combine with the old-vertex terms to
give `(4)`. Finally,

```text
A=27z(z+h)-9zh-27sum_i d_i v_i+9sum_i d_i t_i,        (15)
```

and substituting the definition of `Delta_V` in `(1)` turns `(4)` into
`(5)`. This proves Theorem 1.

## 3. Two retained sidecars

### 3.1 Marked Hamilton arcs

Fix a vertex `i` of a no-sink tournament. In each Hamilton path mark one
consecutive arc whose two endpoints avoid `i`. Cutting the path at the
marked arc gives an ordered two-path cover exposed at that arc. This map is
injective: concatenation across the marked arc recovers both the Hamilton
path and its mark. The exposed objects on edges avoiding `i` are counted by
`z-d_i`.

A Hamilton path has `n-3` arcs avoiding `i` when `i` is internal and one
extra when `i` is either endpoint. Therefore

```text
z-d_i >=(n-3)h+t_i+e_i.                               (16)
```

The start and end terms in `(16)` are genuine retained coordinates; dropping
them weakens the final sector.

### 3.2 Joint endpoint chirality

If `i` is not a universal source, delete `i` from a Hamilton path starting
at `i` and reinsert it immediately after the last in-neighbor of `i` in the
remaining path. Every later vertex is an outneighbor of `i`, so the result
is a Hamilton path not starting at `i`. Deletion recovers the input, hence

```text
i not a universal source  implies  t_i<=h/2.           (17)
```

Reversing all arcs gives the dual statement

```text
i not a universal sink    implies  e_i<=h/2.           (18)
```

Now suppose `C` has no sink. If it has no universal source, `(17)--(18)`
give

```text
tau<=h^2/2,                   sum_i t_i e_i<=h^2/2.
```

If it has a universal source `s`, every Hamilton path starts at `s`, while
none ends there; thus `tau=h^2` and `sum_i t_i e_i=0`. In both cases,

```text
tau+sum_i t_i e_i<=h^2.                                (19)
```

This joint estimate is stronger here than applying only the elementary
separate bounds `tau<=h^2` and `sum_i t_i e_i<=h^2`.

## 4. Proof of the low-capacity sector

THM-4219's Hamilton-path split gives `v_i>=h`. Since
`27v_i-9t_i>=0`, equations `(16)` and `(19)` imply

```text
A
 >=(n-3)h sum_i(27v_i-9t_i)
   +sum_i(t_i+e_i)(27v_i-9t_i)

 >=9(n-3)h(3z+2h)
   +54h^2-9(tau+sum_i t_i e_i)

 >=9(n-3)h(3z+2h)+45h^2.                              (20)
```

Use `(20)`, `tau<=h^2`, and THM-4219's no-sink floor
`Delta_V>=n(n-1)h^2` in `(5)`. With `x=z/h`, this gives

```text
I(C)/h^2 >=L_n(x),                                    (21)

L_n(x)=-18x^2+(27n-171)x+315n^2-297n-225.             (22)
```

The one-defect identity in THM-4219 gives `z>=(n-1)h`, so `x>=n-1`.
On this feasible ray,

```text
L_n'(x)<=-36(n-1)+27n-171=-9n-135<0.                  (23)
```

Moreover,

```text
L_n(n-1)=324n^2-459n-72>1480,              n>=4.      (24)
```

Solving `L_n(x)=1480` for its positive root gives exactly `(6)`.
Monotonicity `(23)` now proves `(7)`. The clean threshold follows from

```text
L_n(5n-11)-1480=531n-2002>0,               n>=4.      (25)
```

For `n=3`, no-sink forces `C3`. Its packet is

```text
(h,z,Delta_V,tau,A)=(3,6,54,3,432),                    (26)
```

and `(5)` gives `(9)`. This proves Theorem 2.

## 5. Proof of the three-eighths implication

Assume `(10)`. Equations `(5)` and `(20)`, together with `tau<=h^2`, give

```text
I(C)/h^2
 >=(801/8)x^2+(765/2)x+513/2
   +9(n-3)(3x+2)+45.                                  (27)
```

For `n>=4`, the right side is increasing in the feasible variable
`x>=n-1`; at its left endpoint it equals

```text
(1017n^2+738n+369)/8 >=19593/8>1480.                  (28)
```

The order-three case is `(26)`, proving Theorem 3 and its equality claim.

Under HYP-9081, its strong-terminal ordinal reduction supplies `(10)` for
every no-sink tournament. Apply Theorem 3 successively to
`C,P_1 triangleright C,...,P_(a-1) triangleright C`; each has Hamilton count
`h`. Summing the increments proves `(11)`. An increment can be equal only
when its input is `C3`, so equality in the sum occurs exactly for the single
step `a=1,C=C3`.

## 6. Exact audit and controls

The reproducible audit covers all labelled tournaments through order five:

```text
all rows by orders 1--5:                 1,099,
no-sink rows by orders 1--5:             0,0,2,32,704,
minimum I/H^2 at no-sink orders 3--5:    1480,4716,8004.
                                                               (29)
```

For every row it reconstructs the source-padded tournament literally,
checks `(12)`, evaluates both sides of `(4)--(5)`, and independently rebuilds
both `R_+` terms from the literal ordinal children. On every no-sink row it
enumerates Hamilton paths and marked arcs, checks `(16)`, `(19)--(22)`, and
the sector implications. Equality occurs only on the two labelled copies of
`C3`.

The transitive order-three tournament is a filter-hostile control:

```text
I(P_3)=-1044.                                          (30)
```

Thus no-sink is load-bearing. Conversely, THM-4219's strong tower `T(9,1)`
lies outside the clean sector while still having positive increment; hence
`(7)--(8)` are sufficient sectors, not equivalences.

Replay:

```bash
python3 -B \
  04-computation/tournament_cycle_left_source_padding_response_sector_thm4221.py
python3 -O -B \
  04-computation/tournament_cycle_left_source_padding_response_sector_thm4221.py
```

Both streams byte-match the frozen output.

## 7. Scope firewall and connection contract

The source object is a no-sink tournament with its oriented capacity tensor
and rooted Hamilton start/end data. Source padding preserves the Hamilton
count but transforms the response coordinates by `(12)`. Projecting to
`(n,H,W,Delta_V)` destroys the marked-avoidance sum `A` and the joint endpoint
term in `(19)`; `(16)` and endpoint chirality are the required sidecars.

This theorem proves neither the universal gap `(10)` nor universal
source-padding monotonicity. HYP-9081 remains open, and the abstract endpoint-
matrix hostile recorded there shows why first-order degree floors alone
cannot supply the missing gap. General `(OS+)` and the base sign
`R_+(C3,C)>0` for every no-sink `C` also remain open.
