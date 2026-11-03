# Exact depth obstruction for the fixed D return

Status: proved as an integer-matrix statement; no tournament-spectrum completeness claim. This independently checks and strengthens the one-step mod-100 obstruction in `05-knowledge/results/return_channels_20261103/INCOMING.md` at commit `aa7514b988295fe0caac7fc06359a28b330dccf3`. The auxiliary script `audit.py` uses only Python's standard library, includes the six relevant atom matrices and source adjacency rows, and checks seed graphs independently by subset dynamic programming.

## Definitions and starting data

Row states are `v=(h,x,y)`. Words act by right multiplication in processing order. Define

    t = (2,-1,0),
    beta(v) = (y, h+2x-y) = v B,
    B = [[0,1],[0,2],[1,-1]].

The fixed nine-letter word `D = r24 a a a b r24 a a a` has

    P = P_D = [[336359,155087,268050],
               [672726,310170,536100],
               [560922,258618,447000]],
    P B = 3 B R,
    R = [[238350,336547],[89350,126161]],
    det R = -100,       t P = -4t.

The actual source seed states are `U=(31,14,25)`, `U2=(95,41,75)`, and `U4=(97,40,75)`, with bulks `(25,34)`, `(75,102)`, `(75,102)`. Their processing words are respectively `r5 a a`, `a a q1 a a`, `q0 a a a a a`.

The seed graph rows are included in `audit.py`; direct endpoint-subset path enumeration gives Hamilton counts `(31,95,97)`. Adding two explicit singletons, each with the newest vertex beating exactly its predecessor, gives `(H_a,H_aa)=(59,109),(177,327),(177,327)`, consistent with the state formulas `H_a=h+2x`, `H_aa=h+2(x+y)`. Atom matrix multiplication reproduces `P` and all three seed states.

## The exact return lattice at every depth

Let `L_k = Z^2 R^k` for integer `k>=1`. Write

    R^k = [[A_k,B_k],[C_k,D_k]].

Then the quotient `Z²/L_k` is cyclic of order `100^k`, and

    (a,b) ∈ L_k  iff  a D_k - b C_k ≡ 0 (mod 100^k).       (1)

Proof: `R mod 10 = [[0,7],[0,1]]` is idempotent, so `D_k ≡ 1 (mod 10)` for every positive `k`. It is therefore invertible modulo `100^k`. The determinant is `(-100)^k`. If `(a,b)=(s,t)R^k`, then `aD_k-bC_k=s det(R^k)`, proving necessity. Conversely the stated divisibility makes the first coordinate `(aD_k-bC_k)/det(R^k)` of `(a,b)(R^k)^(-1)` integral. The numerator of its second coordinate is `-aB_k+bA_k`; multiplying it by `D_k` gives `-B_k(aD_k-bC_k)+b det(R^k)`. Because `D_k` is a unit modulo `100^k`, the second coordinate is also integral. The homomorphism `(a,b) -> aD_k-bC_k mod 100^k` is surjective because `D_k` is a unit, establishing both the index and cyclic quotient.

At `k=1`, (1) simplifies to

    (a,b) ∈ L_1 iff a + 50b ≡ 0 (mod 100).               (2)

Indeed `D_1=126161` is invertible mod 100 with inverse 41, and `C_1*41 ≡ 50 (mod 100)`. In particular, every `(a,b)` in `L_1` has `50|a`. All three seed bulks violate this necessary condition since their first coordinate is odd. This elementary parity obstruction already excludes them, and (2) is the full lattice criterion.

## Exact D-depth theorem

For **any integer state** `v`, and any integers `j>=0`, `k>=1`, the equation

    w P^k = 3^k (v P^j) + m t,    w ∈ Z^3, m ∈ Z            (3)

has a solution **if and only if `beta(v) ∈ L_max(k-j,0)`**, with `L_0=Z²`. Here `w` is an unrestricted formal integer state, not required to be nonnegative, reachable or equipped with a valid port. For each of the three seeds this reduces to **if and only if `k<=j`**, since the seed bulk is outside `L_1` and the lattices are nested.

Proof. The map `beta: Z³ -> Z²` is surjective: `(a+b,0,a)` maps to `(a,b)`. Its integer kernel is precisely `Z t`. Also `beta(zP^n)=3^n beta(z)R^n`. Applying beta to (3), and cancelling `3^k`, makes the equation equivalent to

    beta(w) R^k = 3^j beta(v) R^j.                        (4)

For `k<=j`, choose `beta(w)=3^j beta(v) R^(j-k)` and lift it to an integer state. Then (4) holds, so the difference in (3) belongs to `ker beta = Z t`, yielding some integer `m`.

For `k>j`, invert `R^j` over the rationals in (4):

    beta(w) R^(k-j) = 3^j beta(v).

Since `Z²/L_(k-j)` has order a power of 100, multiplication by `3^j` is invertible on this finite quotient. Hence (4) holds only if `beta(v) ∈ L_(k-j)`. Conversely, if `beta(v)=q R^(k-j)` for some integer row `q`, choose any integer lift `w` with `beta(w)=3^j q`. Then (4) holds and its difference lifts to the required integer carry. This proves the general criterion and the seed corollary.

**Scope.** The no-solution half excludes all integer predecessors and all carries `m` for this fixed D-only equation. The solution half establishes only formal integer witnesses, not physical graph realizability, positivity, or feasible ports. It does not exclude mixed return words, re-marking, different carries outside `Zt`, or alternative seeds. It does not prove the conjectured tournament Hamiltonian-path spectrum.

## Reproduction

Run `python audit.py` and `python -O audit.py`; both modes perform explicit checks and produce identical output. The script multiplies atom matrices, reconstructs seed graphs, counts Hamilton paths and their two singleton continuations by subset DP, verifies `PB=3BR`, `tP=-4t`, seed bulks, the determinant and mod-100 gate, and checks exact-depth constructive integer witnesses for all seeds, `0<=j<=4`, `1<=k<=j`, as well as examples of the general converse with `k>j`. The power lattice proof above is symbolic for all `k`; the script additionally checks powers through eight.
