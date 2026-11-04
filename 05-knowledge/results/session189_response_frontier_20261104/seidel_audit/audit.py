from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
import json,math,subprocess,hashlib
D=Path(__file__).resolve().parent

def need(c,m):
 if not c:raise RuntimeError(m)
def det(A):
 n=len(A);z=0
 for p in permutations(range(n)):
  inv=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n));z+=(-1)**inv*math.prod(A[i][p[i]] for i in range(n))
 return z

def cramer(A):
 d=det(A);need(d%2==1,'odd determinant');n=len(A);out=[]
 for j in range(n):out.append(F(det([[1 if k==j else A[i][k] for k in range(n)] for i in range(n)]),d))
 return out

def data(Q):
 n=len(Q);out=[]
 for mask in range(1<<n):
  ids=[i for i in range(n) if mask>>i&1];odd=len(ids)%2
  if not ids:out.append((0,ids,[]));continue
  A=[[0 if i==j else 1 if Q[i]>>j&1 else -1 for j in ids] for i in ids]
  if odd:A=[[x+1 for x in row] for row in A]
  r=cramer(A);need(sum(r)==odd,'even/odd endpoint normalization');out.append((odd,ids,r))
 return out

def response(Q,a,c):
 n=len(Q);even=F(0);odd=F(0)
 for parity,ids,r in data(Q):
  A=sum(v for i,v in zip(ids,r) if a>>i&1);C=sum(v for i,v in zip(ids,r) if c>>i&1);sign=-1 if parity else 1;v=(1+2*sign*A)**2*(1+2*sign*C)**2
  if parity:odd+=v
  else:even+=v
 return [((even if h==1 else 2*even if h==2 else 0)+(-2)**(h-1)*math.comb(h+1,2)*odd)/2**n for h in range(1,6)]

def mod(x):need(x.denominator%2==1,'integral2-adic response after cancellation');return (x.numerator*pow(x.denominator,-1,256))%256

def main():
 rows=json.loads((D/'quotients.json').read_text());need(len(rows)==7,'seven small canonical types')
 def canonical(T):
  n=len(T);return min(sum(((T[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate((i,j) for i in range(n) for j in range(i+1,n))) for p in permutations(range(n)))
 for n in (2,3,4):
  observed=[canonical(T) for T in rows if len(T)==n];universe=set()
  for code in range(1<<(n*(n-1)//2)):
   T=[0]*n;k=0
   for i in range(n):
    for j in range(i+1,n):
     if code>>k&1:T[i]|=1<<j
     else:T[j]|=1<<i
     k+=1
   universe.add(canonical(T))
  need(len(observed)==len(set(observed)) and set(observed)==universe,'independent full labeled small-quotient coverage')
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'geometric.cpp'),'-o',str(D/'geometric')],check=True);inp=str(len(rows))+'\n'+''.join(str(len(Q))+'\n'+' '.join(map(str,Q))+'\n' for Q in rows);raw=subprocess.check_output([str(D/'geometric')],input=inp.encode()).decode();(D/'run_bank.txt').write_text(raw);lines=[list(map(int,l.split())) for l in raw.splitlines()];at=0;values=[]
 for Q in rows:
  n=len(Q);sides=[]
  for side in range(2):
   T=Q if side==0 else [((1<<n)-1)^(1<<i)^r for i,r in enumerate(Q)];row=lines[at];at+=1;need(row[:3]==[side,-1,-1] and row[3:]==[(-2)**h%256 for h in range(5)],'all-cover geometric scalar fixed point')
   bank={}
   for a in range(1<<n):
    for c in range(1<<n):
     row=lines[at];at+=1;need(row[:3]==[side,a,c],'literal record ordering');f=response(T,a,c);need([mod(x) for x in f]==row[3:],'exact rational Gram versus independent cut-aware run functional');bank[a,c]=f
   sides.append(bank)
  energy=F(0);shared=[]
  for i in range(n):
   for j in range(n):
    if i==j:continue
    K=F(0)
    for mask in range(8):
     a=((mask&1)<<i)|(((mask>>1)&1)<<j);c=((mask>>2)&1)<<i;K+=(-1)**(3-mask.bit_count())*(sides[0][a,c][0]-sides[1][a,c][0])
    energy+=(1 if Q[i]>>j&1 else -1)*K;shared.append([i,j,str(K)])
  squares=sum(sum(v*v for v in r) for parity,ids,r in data(Q) if not parity);squares*=F(2**max(6-n,0),2**max(n-6,0));need(energy==squares and energy>0,'independent cubic extraction and positive energy')
  for q in sides[0]:
   delta=[x-y for x,y in zip(sides[0][q],sides[1][q])];need(delta[1]==2*delta[0] and not any(delta[2:]),'full contrast is linear in cut variable')
  values.append(dict(rows=Q,energy=str(energy),shared_cubics=shared,full_source=[[a,c,*map(str,sides[0][a,c])] for a,c in sides[0]],full_target=[[a,c,*map(str,sides[1][a,c])] for a,c in sides[1]]))
 need(at==len(lines),'no omitted direct run lines');(D/'rational_banks.json').write_text(json.dumps(values,separators=(',',':'))+'\n');print(json.dumps(dict(status='PASS',quotient_types=7,source_and_target_query_coordinates=sum(2*4**len(Q) for Q in rows),ordered_cover_coordinates=5,modulus=256,run_bank_sha256=hashlib.sha256(raw.encode()).hexdigest(),rational_bank_sha256=hashlib.sha256((D/'rational_banks.json').read_bytes()).hexdigest()),sort_keys=True,indent=2))
if __name__=='__main__':main()
