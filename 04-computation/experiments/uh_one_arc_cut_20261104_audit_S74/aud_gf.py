# Bivariate generating function of the cut-arc distribution from the audit's lemma k = 3 + #U - #D on Shelor words,
# compared with the directly computed tables (aud_uh_all_9_16_16_20.out, orders 4-20) and with S74's table (orders 5-16).
import re, sympy as sp
x, z = sp.symbols('x z')
G = z**2 * (x + (1 - x) * z) / (1 - 2 * x - x * (1 - x) * z)
ser = sp.series(G, x, 0, 17).removeO()
direct = {}
for line in open('aud_uh_all_9_16_16_20.out'):
    m = re.match(r'SHELOR n=(\d+):', line)
    if m: cur = int(m.group(1))
    if 'T(n,k) direct:' in line:
        direct[cur] = {int(a): int(b) for a, b in re.findall(r'(\d+):(\d+)', line.split('|')[0].split('direct:')[1])}
bad = 0
for n in range(4, 21):
    coeff = sp.Poly(sp.expand(ser.coeff(x, n - 4)), z) if n - 4 <= 16 else None
    if coeff is None: continue
    gf = {k: int(coeff.coeff_monomial(z**k)) for k in range(0, n + 1) if coeff.coeff_monomial(z**k) != 0}
    if gf != direct[n]: bad += 1; print('MISMATCH n =', n, gf, direct[n])
print('GF vs direct tables, orders 4-20 (GF expanded to x^16 = order 20): mismatches', bad)
s74 = {5: [1,1,1], 6: [2,3,2,1], 7: [4,7,6,3,1], 8: [8,16,16,10,4,1], 9: [16,36,41,30,15,5,1], 10: [32,80,102,85,50,21,6,1],
       11: [64,176,248,231,155,77,28,7,1], 12: [128,384,592,608,456,259,112,36,8,1], 13: [256,832,1392,1560,1289,819,406,156,45,9,1],
       16: [2048,8192,16896,23552,24600,20256,13504,7380,3312,1210,352,78,12,1]}
for n, row in s74.items():
    mine = [direct[n].get(k, 0) for k in range(2, n)]
    print('S74 row n =', n, 'agrees' if mine == row else ('DIFFERS ' + str(mine)))
# column / diagonal formulas
ok = all(direct[n][2] == 2**(n-5) for n in range(5, 21)) and all(direct[n][3] * 2**7 == n * 2**n for n in range(6, 21)) \
     and all(direct[n][4] * 2**10 == (n*n + n - 8) * 2**n for n in range(7, 21)) \
     and all(direct[n][n-1] == 1 and direct[n][n-2] == n-4 and direct[n][n-3] == sp.binomial(n-3, 2) for n in range(5, 21)) \
     and all(direct[n][n-4] == n*(n-4)*(n-5)//6 for n in range(6, 21))
print('columns 2,3,4 and diagonals 1, n-4, C(n-3,2), n(n-4)(n-5)/6 hold for all computed orders <= 20:', ok)
print('row sums = F(2n-6):', all(sum(direct[n].values()) == sp.fibonacci(2*n-6) for n in range(4, 21)))
