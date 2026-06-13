# Class Size Decomposition: 2^m = Odd + Even
## opus-2026-04-03-S27

### The Master Equation
2^m = sum(SC backbone tiling counts) + sum(NS paired tiling counts)
    = sum(odd terms) + sum(even terms)

| n | 2^m | SC sum | NS sum | SC sizes | NS pair sizes (×2) |
|---|-----|--------|--------|----------|---------------------|
| 3 | 2 | 2 | 0 | 1,1 | — |
| 4 | 8 | 6 | 2 | 1,5 | 2 |
| 5 | 64 | 52 | 12 | 1,1,3,5,9,9,11,13 | 2,10 |
| 6 | 1024 | 240 | 784 | 1,1,5,15,15,17,17,29,29,33,37,41 | 2,2,6,10,10,18,18,18,22,22,26,30,38,46,46,50,58,62,66,74,74,86 |

### Key Properties
1. SC sum is always even (sum of even count of odd terms; SC count is even)
2. NS sum is always even (sum of even terms)
3. SC fraction of total: 100%, 75%, 81%, 23% (vanishing at large n)
4. All SC tiling counts are ODD (THM-281)
5. All NS paired tiling counts are EVEN (2× equal partner sizes)
6. Parity alone determines type: odd = backbone, even = paired

### NOTE on |Aut|
The computation |Aut| = n! / (tiling count) is WRONG because tiling count
is the number of base-path-fixed labelings, not the total labeled tournament
count. The orbit-stabilizer theorem gives n!/|Aut| = total labeled tournaments
in the class, which is NOT the same as the tiling count.

The tiling count = number of labeled tournaments in the class that have
the specific Hamiltonian path n→(n-1)→...→1. This depends on H(T).
