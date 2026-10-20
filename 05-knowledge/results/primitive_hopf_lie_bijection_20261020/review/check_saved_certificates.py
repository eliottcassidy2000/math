#!/usr/bin/env python3
"""Independent bounded checks of saved species certificates; no producer imports."""
from fractions import Fraction
from pathlib import Path
import json,csv,math,sys
if len(sys.argv)!=2:raise SystemExit('Usage: check_saved_certificates.py PATH_TO_SPECIES_PACKET')
BASE=Path(sys.argv[1])
checks=0
def ck(x,message):
 global checks
 checks+=1
 if not x:raise ValueError(message)
def det_laplace(a):
 if not a:return Fraction(1)
 return sum(((-1)**j)*x*det_laplace([row[:j]+row[j+1:] for row in a[1:]]) for j,x in enumerate(a[0]) if x)
rows=json.loads((BASE/'primitive_walsh.json').read_text())
for r in rows:
 a=[[Fraction(x) for x in row] for row in r['matrix']]
 d=det_laplace(a)
 ck(d==Fraction(r['determinant']),f'determinant n={r["n"]}')
 ck(all(a[i][j] for i,j in enumerate(r['matching'])),'nonzero matching')
 ck(len(set(r['matching']))==len(a),'matching distinct')
 print('independent Laplace determinant',r['n'],d)
ck(any(Fraction(x).denominator%3==0 for row in rows[-1]['matrix'] for x in row),'odd denominator hostile')
records=list(csv.DictReader((BASE/'primitive_enumeration.csv').open()))
N=10
r=records[:N+1]
def product_coefficients(signed):
 a=[1]+[0]*N
 for n in range(1,N+1):
  factors=[(int(r[n]['even_edges']),1),(int(r[n]['odd_edges']),-1 if signed else 1)]
  for count,sign in factors:
   if count:
    coeff=[sign**j*math.comb(count+j-1,j) for j in range(N//n+1)]
    a=[sum(a[k-j*n]*coeff[j] for j in range(k//n+1)) for k in range(N+1)]
 return a
ck(product_coefficients(False)==[int(x['tournaments']) for x in r],'ordinary product through 10')
ck(product_coefficients(True)==[int(x['self_converse']) for x in r],'signed product through 10')
for n in range(N+1):
 ck(int(r[n]['even_edges'])+int(r[n]['odd_edges'])==int(r[n]['connected_sign_even']),'eigenspace dimension')
 ck(int(r[n]['even_edges'])-int(r[n]['odd_edges'])==int(r[n]['primitive_trace']),'eigenspace trace')
print('Independent ordinary/signed product coefficients through n=10 PASS')
print('CHECKS',checks,'PASS')
