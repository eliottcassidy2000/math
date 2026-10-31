# Pinned Session126 macro proof excerpt

Source: ../primegame/fractran/FRACTRAN_PRIME_EXPONENT_REPORT.md, sections3–4.
The original full-source hash and this excerpt hash are in SOURCE_PINS.json.

## 3. Four counters and seven control states

On the orbit from seed 2, and more generally whenever initially all six control exponents are zero except possibly one equal to one, one has

`e_11+e_13+e_17+e_19+e_23+e_29 ∈ {0,1}`.

Proof: each applicable control-bearing instruction replaces one such prime by another or removes it. The final instruction adds 11 only when none is present, since each control prime has an earlier applicable instruction with that prime alone as denominator. The invariant follows by induction.

Thus `(a,b,c,d)=(ν_2,ν_3,ν_5,ν_7)` are unbounded data counters and the remaining primes select one of seven control states `{1,11,13,17,19,23,29}`, with 1 meaning no control factor. The reduced finite-control table is:

| Control | Prioritized condition | Update and next control |
|---|---|---|
|13|d>0|d−−; go to17 (A)|
|13|otherwise|go to11 (J)|
|17|c>0|c−−,a++,b++; go to13 (B)|
|17|else b>0|b−−; go to19 (C)|
|17|otherwise|go to1 (I)|
|19|a>0|a−−; go to23 (D)|
|19|otherwise|d++; go to11 (H)|
|11|b>0|b−−; go to29 (E)|
|11|otherwise|go to13 (K)|
|29|always|d++; go to11 (F)|
|23|always|c++; go to19 (G)|
|1|a>0|a−−,b++,c++; remain1 (L)|
|1|else d>0|d−−; remain1 (M)|
|1|otherwise|c++; go to11 (N)|

The roles are deliberately asymmetric. In particular, 7 is a data-counter prime, 17 a control prime, and 55 a product introducing a data unit and a control flag.

A separate finite structural audit agrees on 1,792 one-hot-control states with each counter in 0..3, and integer/exponent semantics agree on all 1,024 binary exponent vectors (including states outside the one-hot invariant domain).

The 10×14 integer update matrix has full integer lattice image. Columns `(I,M,J,A,N,E,F,L,C,D)` form a 10×10 minor of determinant −1. Thus no nonzero global linear conservation law, even modulo an integer m>1, is forced by these update vectors. This is not a claim that the priority graph is strongly connected: actual nonnegative reachable states have strong constraints, as the cycles below demonstrate.

## 4. Reconstructed divisor-loop proof and exact word lengths

Use A,...,N for the fourteen instructions. At a canonical trial state

`C(n,d) = 5^n 7^d 13`, with `n≥2` and `1≤d<n`, write `n=qd+r`, `0≤r<d`.

A complete quotient-removal block has instruction word

`W_d = (AB)^d J (EF)^d K`.

Starting at C(n,d), applying `W_d^q` gives `2^(qd)5^r7^d13`. Applying `(AB)^r A` then gives

`2^n 3^r 7^(d−r−1)17`.

Every indicated instruction is the first applicable one; the control table verifies this, including exhaustion tests.

* **If r>0**, the suffix `C (DG)^n H (EF)^(r−1) K` returns to `C(n,d−1)`. No pure power of 2 occurs. The complete block uses **6n+2q+2** FRACTRAN updates.
* **If r=0**, the suffix `I L^n M^(d−1) N (EF)^n K` gives `C(n+1,n)`. Immediately after I the state is `2^n7^(d−1)`. This is a pure power of 2 exactly when d=1. The complete block uses **7n+2q+d+3** updates.

Why the prime theorem follows: seed 2 follows `L N E F K` to `C(2,1)`. For each candidate n, the program tests d=n−1,n−2,... until the first divisor. It then increments n and restarts at d=n. The first divisor is 1 precisely for prime n. The only pure powers of 2 inside these blocks are the stated `2^n` events, so no composite exponent slips through.

The full instruction-word formulas, endpoints, intermediate pure-power events, and step counts were exhaustively checked for all **1,225 pairs** `2≤n≤50,1≤d<n`, involving 267,171 individual updates. These finite checks supplement the symbolic argument; they do not replace it.

More generally, starting at `2^a` with a≥1, the bootstrap word `L^a N (EF)^a K` (length 3a+2) gives `C(a+1,a)`. Subsequent pure-power exponents are exactly the primes greater than a. The special case a=0 fails: seed 1 is trapped as shown below.

