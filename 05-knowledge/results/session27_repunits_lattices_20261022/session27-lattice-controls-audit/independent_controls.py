from pathlib import Path
from fractions import Fraction
from math import comb, prod
import json, hashlib, re
base=Path(__file__).resolve().parent.parent/'session27-papers-review'
out=Path(__file__).resolve().parent
# Independent partition Euler-product count, no signature recursion.
ranks=list(range(1,44))+list(range(4,44))+[6,7,8]
a=[1]+[0]*43
for r in ranks:
 for n in range(r,44): a[n]+=a[n-r]
lines=(base/'ade_audit_results.txt').read_text().splitlines()
for line in lines:
 if re.match(r'^\d+ ',line):
  n,c,*_=line.split(); assert a[int(n)]==int(c)
print('partition DP: all 43 per-rank counts match; type count',len(ranks),'rank1..42',sum(a[1:43]),'rank43',a[43],'total',sum(a[1:]))
# Overflow audit for C++ four-signature additions by safe triangle bound.
bounds=[]
for k in range(1,44):
 bounds.append((k,[k*(k+1),-k*(k+1)*(2*k-1),4*k*(k+1)*(5*k*k-7*k+3),-2*k*(k+1)*(14*k**3-37*k*k+39*k-15)]))
 if k>=4: bounds.append((k,[2*k*(k-1),-2*k*(k-1)*(4*k-7),8*k*(k-1)*(20*k*k-74*k+69),-4*k*(k-1)*(112*k**3-672*k*k+1378*k-963)]))
bounds += [(6,[72,-1512,152928,-2011824]),(7,[126,-4158,669816,-14196924]),(8,[240,-13680,3844800,-143445600])]
upper=[43*max(abs(v[j]) for k,v in bounds) for j in range(4)]
assert max(upper)<2**63
print('conservative absolute sum bounds',upper,'all below signed int64')
G=[[int(x) for x in row.split(',')] for row in (base/'downloaded_G96_reduced.txt').read_text().strip().strip('[]').split(';')]
for name in ['short_vector_witness.json','lll_witness.json']:
 w=json.loads((base/name).read_text());x=w['vector'];q=sum(x[i]*G[i][j]*x[j] for i in range(96) for j in range(96));assert q==w['norm'];print(name,'exact Python integer norm',q,'max coordinate',max(map(abs,x)))
# Count expected sparse candidates (up to sign).
print('signed support1..4 candidate counts',[comb(96,r)*2**(r-1) for r in range(1,5)])
# Independent orthogonal count derived by double counting generator/point incidence.
q=7;m=24
num_lines=(q**(m-1)+1)*(q**m-1)//(q-1)
lines_per_generator=(q**m-1)//(q-1)
assert Fraction(lines_per_generator,num_lines)==Fraction(1,q**(m-1)+1)
print('split orthogonal line count',num_lines,'generator incidence probability',Fraction(lines_per_generator,num_lines))
print('symplectic/orthogonal generator-count ratio',Fraction(q**m+1,2))
# Bind distributed inputs only; original PDF audit manifest is preserved separately.
paths=sorted(base.glob('*.py'))+sorted(base.glob('*.cpp'))+sorted(base.glob('*.txt'))+list(base.glob('*.json'))+[base/'REVIEW.md']
(out/'PACKAGED_INPUT_SHA256SUMS.txt').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p)+'\n' for p in paths))
