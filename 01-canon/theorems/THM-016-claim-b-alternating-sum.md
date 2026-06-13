# THM-016: Claim (B) — Alternating Sum Identity for Tournament Hamiltonian Paths

**Type:** Theorem
**Certainty:** 5 — PROVED
**Status:** PROVED (inductive proof)
**Last reviewed:** kind-pasteur-2026-03-05-S10
**Disputes:** none
**Tags:** #claim-b-alternating #hamiltonian-paths #alternating-sum #induction

---

## Statement

For any tournament on m vertices W with distinguished vertex v:

```
sum_{S subset W\{v}} (-1)^|S| H(S) h_start(W\S, v) = (-1)^{m+1} h_end(W, v)
```

Equivalently (reindex R = W\S, so R always contains v):

```
sum_{R: v in R subset W} (-1)^{|R|+1} H(W\R) h_start(R, v) = h_end(W, v)
```

where:
- H(S) = total Hamiltonian path weight on S (H(empty) = 1)
- h_start(R, v) = Hamiltonian path weight on R starting at v
- h_end(W, v) = Hamiltonian path weight on W ending at v
- Tournament condition: T(a,b) + T(b,a) = 1 for all a != b

**NOTE:** This is NOT the paper's "Claim B" (THM-003, about I(Omega,2) recurrence). This is a
standalone identity about tournament Hamiltonian paths, discovered by opus-S4c and proved
by kind-pasteur-S10.

---

## Proof

By induction on m = |W|.

**Base case (m=1):** W = {v}. Only R = {v}.
LHS = (-1)^2 * H(empty) * h_start({v}, v) = 1 * 1 * 1 = 1.
RHS = h_end({v}, v) = 1. Done.

**Inductive step (m >= 2):** Assume the identity for all W' with |W'| < m.

Step 1: Expand h_start(R, v) for |R| >= 2:

```
h_start(R, v) = sum_{w in R\{v}} T(v,w) * h_start(R\{v}, w)
```

Step 2: Split the LHS into |R|=1 and |R|>=2 terms:

```
LHS = H(W\{v}) * 1                                           [R={v} term]
    + sum_{R, |R|>=2} (-1)^{|R|+1} H(W\R) * sum_w T(v,w) h_start(R\{v}, w)
```

Step 3: Exchange order of summation (collect by w in W\{v}):

```
= H(W\{v}) + sum_{w in W\{v}} T(v,w) * [inner sum over R containing {v,w}]
```

Step 4: Substitute R' = R\{v} in the inner sum. Since |R| = |R'|+1:

```
inner sum = sum_{R': w in R' subset W\{v}} (-1)^{|R'|} H((W\{v})\R') h_start(R', w)
          = -[sum_{R'} (-1)^{|R'|+1} H((W\{v})\R') h_start(R', w)]
          = -h_end(W\{v}, w)    [by inductive hypothesis on W\{v} with vertex w]
```

Step 5: Substitute back:

```
LHS = H(W\{v}) - sum_{w in W\{v}} T(v,w) h_end(W\{v}, w)
    = sum_w h_end(W\{v}, w) - sum_w T(v,w) h_end(W\{v}, w)     [since H = sum_w h_end]
    = sum_w [1 - T(v,w)] h_end(W\{v}, w)
    = sum_w T(w,v) h_end(W\{v}, w)                               [tournament: T(v,w)+T(w,v)=1]
    = h_end(W, v)                                                  [path extension]
```

The last step: a Hamiltonian path on W ending at v has some second-to-last vertex w,
giving T(w,v) * (Ham path on W\{v} ending at w). Summing over w gives h_end(W, v). QED.

---

## Dual Identity

By path reversal (T -> T^op, where T^op(a,b) = T(b,a)), h_start <-> h_end and H is preserved.
Applying to the identity gives the dual:

```
sum_{S subset W\{v}} (-1)^|S| H(S) h_end(W\S, v) = (-1)^{m+1} h_start(W, v)
```

Verified computationally (hypothesis H1/H9 in q009_creative_hypotheses.py).

---

## Verification Record

| m | Trials | Max error | Status |
|---|--------|-----------|--------|
| 1 | 100 | 0 | VERIFIED |
| 2 | 100 | 0 | VERIFIED |
| 3 | 100 | 4.4e-16 | VERIFIED |
| 4 | 100 | 1.8e-15 | VERIFIED |
| 5 | 100 | 8.0e-15 | VERIFIED |
| 6 | 100 | 2.6e-14 | VERIFIED |
| 7 | 30 | 2.1e-13 | VERIFIED |
| 8 | 30 | 1.9e-12 | VERIFIED |

All proof steps (direct computation, inductive formula, path extension, inductive hypothesis)
verified independently at each m.

Code: 04-computation/q009_claim_b_proof.py

---

## Tournament Condition is Essential

The identity FAILS for general digraphs where T(a,b) + T(b,a) != 1.
Tested at m=3: 282/500 random digraphs fail (hypothesis H2 in q009_creative_hypotheses.py).
The tournament condition is used in Step 5 where T(w,v) = 1 - T(v,w).

---

## Consequences

1. **B(Li,Rj) = B(Lj,Ri):** The alternating subset convolution symmetry holds for all n.
   (Follows by the algebraic reduction in signed-adjacency-identity.md: the s-linear
   coefficients of B(Li,Rj)-B(Lj,Ri) vanish using THM-016 for the d-independent part
   and the inductive hypothesis for the d-dependent part.)

2. **Even-odd split for all n:** sum_{|S| even} Delta(S,R) = sum_{|S| odd} Delta(S,R).

3. **Signed position identity for all n:** sum_{P: i->j} (-1)^{pos(i)} = sum_{P': j->i} (-1)^{pos(j)}.

**CAVEAT (MISTAKE-008):** These consequences are NECESSARY conditions for OCF but NOT
sufficient. The even-odd split is strictly weaker than OCF. Proving OCF requires an
additional identity connecting the odd-S sum of Delta to the cycle formula. See OPEN-Q-009.

---

## History

- **Discovery:** opus-2026-03-05-S4c (isolated as standalone identity, verified m=1,...,8, proved m=1,2,3)
- **General proof:** kind-pasteur-2026-03-05-S10 (inductive proof using h_start expansion + path extension)
- **Gap identification:** kind-pasteur-2026-03-05-S10 (the step B(Li,Rj)=B(Lj,Ri) -> OCF has a gap per MISTAKE-008)
