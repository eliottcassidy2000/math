# Tournament A -> char poly of adjacency = a BINARY FORM phi_T in Sym^n P^1.
# Transitive <=> phi_T = x^n (rational normal curve vertex, nilpotent). Intransitivity first seen at
# trace(A^3)=3*(#3-cycles). Newton's identities give char poly EXACTLY from power sums trace(A^k).
from itertools import product, permutations
from fractions import Fraction as Fr
def matpow_traces(A,n):
    # trace(A^k), k=1..n, exact integer
    M=[row[:] for row in A]; tr=[sum(M[i][i] for i in range(n))]
    P=[row[:] for row in A]
    for k in range(2,n+1):
        P=[[sum(P[i][j]*A[j][l] for j in range(n)) for l in range(n)] for i in range(n)]
        tr.append(sum(P[i][i] for i in range(n)))
    return tr  # [tr A^1,...,tr A^n]
def charpoly_from_traces(p,n):
    # Newton: e_k from power sums p_k. char poly = x^n - e1 x^{n-1}+ e2 x^{n-2} - ...
    e=[Fr(1)]
    for k in range(1,n+1):
        s=Fr(0)
        for i in range(1,k+1): s+=(-1)**(i-1)*e[k-i]*p[i-1]
        e.append(s/k)
    # coefficients of char poly: c_j = (-1)^j e_j (coeff of x^{n-j})
    return [ (-1)**j*e[j] for j in range(n+1) ]  # [1, -e1, e2, ...] = coeffs high->low
def tour(n,bits):
    A=[[0]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
            idx+=1
    return A
print("=== tournament -> char poly (binary form); transitive = x^n; c3 = trace(A^3)/3 ===")
for n in (3,4,5):
    polys=set(); trans_poly=None; c3s=set()
    for bits in range(1<<(n*(n-1)//2)):
        A=tour(n,bits); tr=matpow_traces(A,n); cp=charpoly_from_traces(tr,n)
        polys.add(tuple(cp)); c3s.add(tr[2]//3 if n>=3 else 0)
        # transitive = the one with all traces 0
        if all(t==0 for t in tr): trans_poly=cp
    print(f" n={n}: {2**(n*(n-1)//2)} tournaments, {len(polys)} distinct char polys (cospectral classes);"
          f" transitive char poly = {[str(x) for x in trans_poly]} (=x^{n}); #3-cycles range {sorted(c3s)}")
# Paley n=7 (QR circulant): eigenvalues on Re=-1/2
print("\n=== Paley tournament (n=7, QR): the OTHER pole -- eigenvalues on Re=-1/2 ===")
QR={1,2,4}  # quadratic residues mod 7
A7=[[1 if ((j-i)%7) in QR else 0 for j in range(7)] for i in range(7)]
import numpy as np
ev=np.linalg.eigvals(np.array(A7,dtype=float))
print("  eigenvalues:", sorted([complex(round(z.real,3),round(z.imag,3)) for z in ev], key=lambda z:(z.real,z.imag)))
print("  => Perron 3, and (-1 +- i*sqrt7)/2 (Re=-1/2, the 'RH line'); char poly = (x-3)(x^2+x+2)^3")
print("\nDICTIONARY: transitive = X^n (rat'l normal curve vertex, nilpotent, GIT-nullcone deepest point);")
print("Paley/regular = spectrum on Re=-1/2 (cyclotomic/Gauss-sum). Symmetry => spectral degeneracy.")
