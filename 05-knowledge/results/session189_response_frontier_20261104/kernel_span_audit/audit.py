#!/usr/bin/env python3
"""Independent dense squarefree jets, symbolic span, and literal named controls."""
from pathlib import Path
from itertools import combinations,permutations
from math import comb,factorial
from fractions import Fraction as F
import hashlib,json,sys,subprocess,random
P=Path(__file__).resolve().parent

def need(b,s):
 if not b:raise RuntimeError(s)
def pmul(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out
def kernels(m):
 L=m-2;k=(m-1)//2;out=[]
 for r in list(range(2,k+1))+[m-1]:
  u=r-1;p=[0]*(L+1)
  for j in range(u+1):p[L-u+j]+=comb(u,j)
  for j in range(L-u+1):p[u+j]-=comb(L-u,j)
  vals=[sum(a*s**j for j,a in enumerate(p))for s in range(L)]
  x=[2*sum((-1)**(h-j)*comb(h,j)*vals[j]for j in range(h+1))for h in range(L)]
  out.append((r,p,x))
 return out

def inverse(A):
 n=len(A);B=[[F(x)for x in row]+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
 for j in range(n):
  q=next(i for i in range(j,n)if B[i][j]);B[j],B[q]=B[q],B[j];v=B[j][j];B[j]=[x/v for x in B[j]]
  for i in range(n):
   if i!=j:
    v=B[i][j];B[i]=[x-v*y for x,y in zip(B[i],B[j])]
 return [r[n:]for r in B]

def polynomials():
 result=[]
 for m in range(3,52,2):
  k=(m-1)//2;ks=kernels(m);wbase=[]
  for h in range(k):wbase.append([comb(h,j-h)if h<=j<=2*h else 0 for j in range(2*k-1)])
  converted=[]
  for r,p,x in ks:
   refl=[sum(p[j]*(-1)**j*comb(j,h)for j in range(h,len(p)))for h in range(len(p))];need(refl==p,'s reflection')
   refl=[sum(x[j]*(-1)**j*comb(j,h)for j in range(h,len(x)))for h in range(len(x))];need(refl==x,'t reflection')
   coeff=[];rest=x[:]
   for h in range(k):
    a=rest[h];coeff.append(a);rest=[v-a*z for v,z in zip(rest,wbase[h])]
   need(not any(rest),'full w representation');converted.append(coeff)
  inv=inverse(converted)
  # Matrix rows are kernels in the w basis: transpose inverse extracts every w^h.
  certificates=[]
  for h in range(k):
   coeff=inv[h];recovered=[sum(coeff[j]*ks[j][2][z]for j in range(k))for z in range(2*k-1)]
   need(recovered==wbase[h],'exact rational physical-kernel certificate');certificates.append([str(x)for x in coeff])
  result.append({'m':m,'selected_r':[r for r,_,_ in ks],'kernels':[x for _,_,x in ks],'w_coordinates':converted,'w_certificates':certificates})
 jets=[]
 for d in range(1,31):
  for j in range(1,d+1):
   out=[F(0)]*(d+1);coeff=[]
   for h in range(j,d+1):
    a=(-1)**(h-j)*F(j,h)*comb(2*h-j-1,h-j);coeff.append(str(a))
    for z in range(h,min(2*h,d)+1):out[z]+=a*comb(h,z-h)
   need(out==[F(z==j)for z in range(d+1)],'inverse jet formula')
   jets.append({'d':d,'j':j,'coefficients':coeff})
 return {'span':result,'inverse_jets':jets}

def mul(a,b):
 N=len(a);out=[0]*N
 for s in range(N):
  t=s
  while True:
   out[s]+=a[t]*b[s^t]
   if not t:break
   t=(t-1)&s
 return out
def add(a,b,c=1):return [x+c*y for x,y in zip(a,b)]
def paths(rows):
 n=len(rows);N=1<<n;G=[[[0]*N for _ in rows]for _ in rows]
 # Endpoint subset recurrence, independent from producer's permutation enumerator.
 for start in range(n):
  dp=[[0]*n for _ in range(N)];dp[1<<start][start]=1
  for mask in range(N):
   for end,v in enumerate(dp[mask]):
    if not v:continue
    G[start][end][mask]+=v;avail=rows[end]&~mask
    while avail:
     b=avail&-avail;avail-=b;dp[mask|b][b.bit_length()-1]+=v
 return G

def maximum(rows,O,I):
 # Exhaustive matching lists, not producer's recursive cardinality optimizer.
 pairs=[(u,v)for u in O for v in I if rows[u]>>v&1];valid=[()]
 for e in pairs:
  valid += [M+(e,)for M in valid if all(e[0]!=a and e[1]!=b for a,b in M)]
 q=max(map(len,valid));return q,[M for M in valid if len(M)==q]

def control(rows,into):
 n=len(rows);N=1<<n;full=N-1;d=(n+1)//2;O=[i for i in range(n)if not(into>>i&1)];I=[i for i in range(n)if into>>i&1];q,matchings=maximum(rows,O,I);G=paths(rows);zero=[0]*N;one=[1]+[0]*full
 S=[[sum(G[i][j][s]for j in range(n))for s in range(N)]for i in range(n)]
 E=[[sum(G[j][i][s]for j in range(n))for s in range(N)]for i in range(n)]
 alpha=[one[s]+sum(E[i][s]for i in I)for s in range(N)];beta=[one[s]+sum(S[i][s]for i in O)for s in range(N)]
 tau=[sum(G[u][v][s]for u in O for v in I)for s in range(N)]
 ell=[[sum(G[i][v][s]for v in I)for s in range(N)]for i in range(n)];rho=[[sum(G[u][i][s]for u in O)for s in range(N)]for i in range(n)];ab=mul(alpha,beta)
 powers=[one]
 for j in range(d+2):powers.append(mul(powers[-1],tau))
 need(not any(powers[q+1]),'matching nilpotence')
 moments=[]
 for i in range(n):
  A=add(mul(alpha,S[i]),mul(beta,E[i]),-1);B=mul(ab,add(ell[i],rho[i],-1))
  need(all(not A[s]for s in range(N)if s.bit_count()<2),'A minimum support')
  c=[mul(A,powers[j])[full]+(j*mul(B,powers[j-1])[full]if j else 0)for j in range(d+2)]
  need(c[-1]==0,'finite jet cutoff');moments.append(c[:-1])
 perfect=2*q==n
 if not perfect:
  witness=q+1
  unmatched={i for M in matchings for i in range(n)if all(i not in e for e in M)}
  for i in range(n):
   sign=1 if i in I else -1;val=moments[i][witness]*sign;need(val>=0,'nonperfect sign')
   if i in unmatched:need(val>=factorial(q+1),'nonperfect constructive witness')
   expected=(q+1)*mul(mul(ab,ell[i]if i in I else rho[i]),powers[q])[full]
   need(val==expected,'nonperfect exact formula')
 elif q>=2:
  witness=q;J=zero[:]
  for u in O:
   for v in I:
    if rows[u]>>v&1:J[(1<<u)|(1<<v)]=1
  JP=one
  for _ in range(q-1):JP=mul(JP,J)
  for i in range(n):
   back=[v for v in I if rows[v]>>i&1]if i in O else[u for u in O if rows[i]>>u&1]
   value=q*(sum(JP[full^(1<<i)^(1<<v)]for v in back)+powers[q-1][full^(1<<i)])
   need(moments[i][q]==(-value if i in O else value),'perfect exact formula')
  for side in [O,I]:
   need(any(moments[i][q]for i in side),'strict each side')
   for i,j in combinations(side,2):need(abs(moments[i][q])+abs(moments[j][q])>=factorial(q),'constructive any pair strict')
 else:
  witness=0;need([moments[i][0]for i in O]==[2]and[moments[i][0]for i in I]==[-2],'perfect one-edge exception')
 # Reconstruct all moments from the w battery; choose enough A names to saturate.
 m=2*d+3;ks=kernels(m);quartics=[[sum(a*b for a,b in zip(x,c))for _,_,x in ks]for c in moments]
 need(any(v for row in quartics for v in row),'selected quartic nonzero')
 for side in [O,I]:
  sidequartic=[sum(quartics[i][r]for i in side)for r in range(len(ks))]
  sideraw=[sum(moments[i][j]for i in side)for j in range(d+1)]
  if any(sideraw):need(any(sidequartic),'homogeneous maximal grouping nonzero')
 return {'rows':rows,'into':into,'matching':q,'maximum_matchings':[[list(e)for e in M]for M in matchings],'perfect':perfect,'raw_moments':moments,'witness_index':witness,'saturated_A_order':m,'selected_quartics':quartics}

def graph(n,mask):
 rows=[0]*n
 for b,(i,j)in enumerate(combinations(range(n),2)):
  if mask>>b&1:rows[i]|=1<<j
  else:rows[j]|=1<<i
 return rows

def physical(m,rows,into,reverse=False):
 out=[sum(1<<j for j in range(m)if (j<i if reverse else j>i))for i in range(m)]+[r<<m for r in rows]
 for i in range(m):
  for j in range(len(rows)):
   if into>>j&1:out[m+j]|=1<<i
   else:out[i]|=1<<(m+j)
 return out

def main():
 src=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else P.parent/'session189-kernel-span';manifest=json.loads((src/'RELEASE_FILES.json').read_text());pins={f:hashlib.sha256((src/f).read_bytes()).hexdigest()for f in manifest};poly=polynomials();(P/'polynomial_certificates.json').write_text(json.dumps(poly,separators=(',',':'))+'\n');records=[]
 for n in range(1,5):
  for z in range(1<<(n*(n-1)//2)):
   rows=graph(n,z)
   for into in range(1<<n):records.append(control(rows,into))
 rng=random.Random(189)
 for n in range(5,9):
  for z in [0,(1<<(n*(n-1)//2))-1]+[rng.getrandbits(n*(n-1)//2)for _ in range(12)]:
   rows=graph(n,z)
   for into in sorted({0,(1<<n)-1,(1<<(n//2))-1,sum(1<<i for i in range(n)if i%2),rng.randrange(1<<n)}):records.append(control(rows,into))
 (P/'moment_bank.json').write_text(json.dumps(records,separators=(',',':'))+'\n')
 # Reproduce the producer's physical examples, retaining complete literal polynomial banks.
 producer=json.loads((src/'receipt.json').read_text());examples=producer['physical_checks'];graphs=[];contrasts=[]
 for ex in examples:
  m=ex['m'];erows=ex['exterior']['rows'];into=ex['exterior']['into'];c=control(erows,into);need(c['raw_moments']==ex['exterior']['raw_moments'],'producer moment controls independent');graphs.extend([physical(m,erows,into),physical(m,erows,into,True)])
 text=str(len(graphs))+'\n'+''.join(str(len(r))+'\n'+' '.join(map(str,r))+'\n'for r in graphs);(P/'physical_graphs.txt').write_text(text);subprocess.run(['g++','-O3','-std=c++17',str(P/'literal_polynomial.cpp'),'-o',str(P/'literal_polynomial')],check=True);data=subprocess.check_output([str(P/'literal_polynomial')],input=text,text=True);(P/'physical_banks.json').write_text(data);banks=json.loads(data);total=0
 for j,ex in enumerate(examples):
  left,right=banks[2*j:2*j+2];need(left['rows']==graphs[2*j]and right['rows']==graphs[2*j+1],'physical names fixed');total+=len(left['coefficients'])+len(right['coefficients']);ix={a:i for i,a in enumerate(left['masks'])};L=len(ix);m=ex['m'];checks=[]
  for check in ex['checks']:
   r=check['r'];i=check['i'];pos=ix[3]*L+ix[(1<<r)|(1<<i)];value=left['coefficients'][pos]-right['coefficients'][pos];need(value==check['quartic_difference'],'literal physical complete bank quartic');checks.append({'r':r,'i':i,'source':left['coefficients'][pos],'target':right['coefficients'][pos],'difference':value})
  # Maximal same-incidence grouping is coefficient summation, with saturated c degree.
  into=ex['exterior']['into'];grouped=[]
  for side in [0,1]:
   vals={r:sum(v['difference']for v in checks if v['r']==r and bool(into>>(v['i']-m)&1)==bool(side))for r in sorted({c['r']for c in checks})};grouped.append(vals)
  need(any(v for g in grouped for v in g.values()),'literal grouped selected quartics');contrasts.append({'m':m,'physical_source_rows':left['rows'],'physical_target_rows':right['rows'],'checks':checks,'grouped_differences':grouped})
 (P/'physical_contrasts.json').write_text(json.dumps(contrasts,sort_keys=True,indent=2)+'\n');need(pins=={f:hashlib.sha256((src/f).read_bytes()).hexdigest()for f in manifest},'producer source unchanged')
 result={'status':'PASS','independent_exterior_cases':len(records),'kernel_orders':list(range(3,52,2)),'inverse_jet_maximum_order':30,'literal_physical_graphs':len(banks),'literal_full_named_coefficients':total,'literal_physical_quartic_contrasts':sum(len(x['checks'])for x in contrasts),'producer_pins':pins,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['audit.py','literal_polynomial.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['polynomial_certificates.json','moment_bank.json','physical_graphs.txt','physical_banks.json','physical_contrasts.json']}}
 (P/'receipt.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
