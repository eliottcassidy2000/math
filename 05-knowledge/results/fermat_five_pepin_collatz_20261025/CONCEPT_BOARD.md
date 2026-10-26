# Concept board (fourth directive)

Anchor: the Fermat line and the Pépin dictionary. Niche: the five-ladder with its tournament floor.
Wildcard: reformulations of completeness and the hostile that defeats the congruence route.

| Concept | Object / operation | Preserved | Lost / sidecar | Decisive test |
|---|---|---|---|---|
| Pépin | `3^{2^k} mod F_n`, squaring tower | primality ⟺ single unit cycle | nothing (exact) | `n = 5`: order of 3 is `11167360`, not `2^{32}` |
| affine Collatz branch mod `F_n` | `x -> 3x+1` on `Z/F_n` | conjugate to `x -> 3x`; fixed point `-1/2` | the halving branch (not a residue map) | cycle type `(1, 2^{2^n})` iff prime |
| Fermat line | `2^{2^n}+1` under Syracuse | run of `2^{n-1}-1` twos, one three, landing `(3^{2^{n-1}}+1)/2` | primality (line is prime-blind) | `F_4` and `F_5` both follow Theorem F |
| merge of consecutive lines | `2^m+1`, `2^{m+1}+1`, `m = 2L+1` | exact merge depth `L+2` (even `L`), `L+3` (`L ≡ 3 mod 4`) | `L ≡ 1 (mod 4)` merges later | `L = 2..19` receipts |
| five-ladder | index → `2^i` → `2^{2^i}+1` → `H(NT_{2^i+2})` | the tower and its square-root height | nothing (classical) | product `2^{32}-1`; 32 Sierpiński rows |
| other fives | Collatz cycles, order-27 groups, Niemeier primes, Platonic solids | each its own classification | no common carrier | the repository's prior "not the same five" ruling |

Hostiles used: the composite `F_5` (prime-blindness of the line), the prime-2-free heuristic, `x = 1`
(infinite valuation-two run; guarded), and the odd-`L` merge failure. Stopping reasons: completeness is
open and every reformulation is faithful but mechanism-free; the Paley family and the Bernoulli
denominators restate rather than attack.
