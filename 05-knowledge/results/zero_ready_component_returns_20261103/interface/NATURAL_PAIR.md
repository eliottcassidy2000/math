# The natural pair has two different exact observation sizes

**PROVED:** an O/D reader has a closed ten-coordinate carrier for the full
saved natural-pair endpoint/P2/W2 profile. Its future two-port H responses
factor further through three numbers at each TT-prefix wrapper boundary.
The independent spare6 does not share this reduction.

This is a refinement of [NOTE.md](NOTE.md), recovering the earlier universal
root-child factorization rather than treating every pair as a generic17-state
source. Exact matrix identities and the two actual 21-vertex controls are
recorded in `natural_pair.json`.

## 1. Ten coordinates close the natural pair's full saved profile

Use the component polynomials and D/O operations from NOTE section2. Extend
the closed eight-state reader vector by E3=E(3) and Q3=Q(3). Thus

```
v=(h,e,S2,F,Q2,q,P2,S3,E3,Q3).
```

The first eight updates remain those in NOTE. For the final two coordinates:

```
D: (E3,Q3) -> (2E3,3Q3),
O: (E3,Q3) -> (2E3+Q3,3Q3).
```

Proof: evaluate `E_D(z)=(z-1)E(z)`,
`E_O(z)=(z-1)E(z)+Q(z)`, and `Q_D(z)=Q_O(z)=zQ(z)` at z=3.
This is an exact closed ten-state carrier for arbitrary D/O words, on
nontrivial rooted components. Neither minimality nor closure under R/N is
claimed.

For the wrapper `T(TT_k,B,w)` with new root r, k>=1, its root fields are in
NOTE. The old mark w has

```
(e_w,s_w,q_w,F_w,S_w(2))
  = (e,q,q,2^(k+1)F-q,S2+2^k Q2).
```

Apply the R polynomial law to the old tournament `P=TT_k -> B`, rooted at w:

```
E_w,new(z)=z^(k+1)E(z)-(z-1)^k Q(z-1),
S_w,new(z)=z^k Q(z)+(z-1)^(k+1)S(z),
Q_w,new(z)=z^(k+1)Q(z).
```

Evaluation at1,2 proves the five fields. With row/column port order(r,w),
the four ordered W2 values are

```
(0,0,2^k F,0).
```

Indeed, a first-bin-ending-r and second-bin-starting-w arrangement would
force r to be a singleton: w cannot simultaneously precede r in its first
bin. The nonempty A prefix then cannot fit in the w-starting second bin,
so `W_(r,w)=0`. For `W_(w,r)`, delete the initial r from the second bin.
This gives exactly the two-bin arrangements on `TT_k -> B` whose first bin
ends at w, counted by `2^k F`; reinsertion is always legal. The diagonal
entries vanish because a vertex cannot occur in two nonempty bins.

Together with

```
P_T(2)=2[3^k E3+2^k P2-S2],
```

the ten-state reader therefore carries every entry of the saved natural
pair `(r,w)` endpoint/P2/W2 profile, and both rooted golden triples.
It does not carry an independent old spare u without additional mixed data.

## 2. For future H alone, the natural pair factors through three numbers

The existing theorem is **R1, ROOTED_RESETS.md**, in
[`two_port_return_induction_20261103`](../../two_port_return_induction_20261103/ROOTED_RESETS.md).
After a one-inport child x is born from an arbitrary rooted P at w, retaining
ordered pair(w,x), every future two-port H response depends exactly on the
old four-profile `(H(P),s_w(P),H(P-w),E_w,P(2))`.
The original theorem used the six-command query basis. Its coefficient
factorization extends to the universal static17 basis because both span the
same full future-observation space on the universal25 carrier.

Here no change-of-basis assertion is left implicit: `natural_pair.py`
recomputes the sealing factor on the universal25 carrier with the seventeen
static query words used by `phase.json`. The five one-port cover signatures
have coefficients

```
h-s-e, e, s, q, F-e-q.
```

It applies the birth that keeps old w in slot0 and puts newborn x in slot1,
and verifies that the e coefficient cancels in all17 query columns. It then
applies the exact endpoint-label swap and repeats the check. Thus both
port orders are retained explicitly, and the resulting two integer4-by17
matrices in `(h,s,q,F)` coordinates are saved in the receipt. Their ranks
are4. These are coefficient identities, not fits to the new target values.

For `P=TT_k -> B` and k>=1, its profile at w is

```
(H(P),s_w(P),H(P-w),E_w,P(2)) = (h,0,q,2^k F).
```

The zero entry is structural: a Hamilton path cannot start in B and then
visit the nonempty prefix. Hence for the actual natural pair(r,w), every
permitted future H response is determined by only

```
(h,q,2^k F).
```

The restricted formal factor has rank3. No assertion of a minimal
three-dimensional realization on an all-size actual strong-component class
is needed here. The exact factor agrees with all17 independently gap-counted
queries of both targets:

| H | `(h,q,4F)` |
|---:|---|
| 512423 | `(49955,3,462468)` |
| 512425 | `(72701,1,439724)` |

This distinguishes two tasks that need different information. Seventeen
future-H queries at the natural root/mark boundary have a three-number
factor, but the saved full P2/endpoints and golden phase use additional
component evaluations. The independent spare6 still introduces the mixed
three-bin quantities and the `3^k` mode in NOTE section3. None of these
reductions identifies a source-changing component reset with an allowed
future birth word.

## 3. Tiny exact checks

Run from the repository root:

```
python3 05-knowledge/results/zero_ready_component_returns_20261103/interface/natural_pair.py
python3 -O 05-knowledge/results/zero_ready_component_returns_20261103/interface/natural_pair.py
```

The script checks exact integer coefficient cancellations, both ordered
profile factors and their rational ranks, all34 actual target query values,
the full natural-pair profiles, and eight D/O ten-state updates on four
order-four controls. It adds no graph census. All gates use explicit
exceptions, and all sources are pinned in `natural_pair.json`.
