#!/usr/bin/env python3
"""Exact mixed-kernel and full physical-bank controls. Standard library only."""
import itertools,json,math,hashlib,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)

def derivative(p):return [(j+1)*p[j+1]for j in range(len(p)-1)]
def formula(m,r):
 L=m-2;u=r-1
 vals=[(s+1)**u*s**(L-u)-s**u*(s+1)**(L-u)for s in range(m-2)]
 out=[]
 while vals:
  out.append(2*vals[0]);vals=[b-a for a,b in zip(vals,vals[1:])]
 return out

def weights(m):
 w={m-1:1}
 if m%2:
  for r in range(2,(m-1)//2+1):w[r]=(-1)**r*math.comb(m-2,r-1)
 else:
  for r in range(2,(m-2)//2+1):w[r]=(-1)**r*(math.comb(m-3,r-1)-math.comb(m-3,r-2))
 return w

def literal_inner(n):
 keys=[(p,q,r)for p in range(n)for q in range(p+1,n)for r in range(n)if r not in(p,q)]
 X={k:[0]*(n-1)for k in keys};Y={k:[0]*(n-2)for k in keys}
 def poly(p,cuts):
  out=[1]
  for j,(a,b)in enumerate(zip(p,p[1:]),1):
   if j in cuts:continue
   nxt=[0]*(len(out)+1)
   for d,v in enumerate(out):
    nxt[d+1]+=v
    if a<b:nxt[d]+=v
   out=nxt
  return out
 def add(D,p,q,r,z,sign):
  if len({p,q,r})<3:return
  if p>q:p,q=q,p
  row=D[p,q,r]
  for j,x in enumerate(z):row[j]+=sign*x
 for p in itertools.permutations(range(n)):
  for c in range(1,n):
   z=poly(p,{c});add(X,p[c-1],p[c],p[-1],z,1);add(X,p[c-1],p[c],p[0],z,-1)
  for c,d in itertools.combinations(range(1,n),2):
   z=poly(p,{c,d})
   add(Y,p[d-1],p[d],p[-1],z,1);add(Y,p[d-1],p[d],p[c],z,-1)
   add(Y,p[c-1],p[c],p[d],z,-1);add(Y,p[c-1],p[c],p[-1],z,1)
 return keys,X,Y

def plus(a,b,scale=1):
 out=dict(a)
 for s,v in b.items():out[s]=out.get(s,0)+scale*v
 return {s:v for s,v in out.items()if v}
def mul(a,b):
 out={}
 for s,v in a.items():
  for t,w in b.items():
   if not(s&t):out[s|t]=out.get(s|t,0)+v*w
 return {s:v for s,v in out.items()if v}
def exterior(rows,into):
 n=len(rows);G=[[{}for _ in rows]for _ in rows]
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(rows[a]>>b&1 for a,b in zip(p,p[1:])):
    mask=sum(1<<i for i in p);g=G[p[0]][p[-1]];g[mask]=g.get(mask,0)+1
 S=[{}for _ in rows];E=[{}for _ in rows]
 for u in range(n):
  for v in range(n):S[u]=plus(S[u],G[u][v]);E[v]=plus(E[v],G[u][v])
 alpha={0:1};beta={0:1};tau={};ells=[];rhos=[]
 for v in range(n):
  if into>>v&1:alpha=plus(alpha,E[v])
  else:beta=plus(beta,S[v])
 for u in range(n):
  for v in range(n):
   if not(into>>u&1)and into>>v&1:tau=plus(tau,G[u][v])
 for i in range(n):
  ell={};rho={}
  for v in range(n):
   if into>>v&1:ell=plus(ell,G[i][v])
   else:rho=plus(rho,G[v][i])
  ells.append(ell);rhos.append(rho)
 return alpha,beta,tau,S,E,ells,rhos

def endpoints(rows,start=None):
 n=len(rows);N=1<<n;dp=[[0]*n for _ in range(N)]
 if start is None:
  for i in range(n):dp[1<<i][i]=1
 else:dp[1<<start][start]=1
 for mask in range(1,N):
  if start is not None and not(mask>>start&1):continue
  for i,x in enumerate(dp[mask]):
   if not x:continue
   opts=rows[i]&~mask
   while opts:
    b=opts&-opts;opts-=b;dp[mask|b][b.bit_length()-1]+=x
 return dp

def capacities(rows):
 n=len(rows);full=(1<<n)-1;E=endpoints(rows);rev=[sum(1<<j for j in range(n)if rows[j]>>i&1)for i in range(n)];S=endpoints(rev)
 return [[sum(E[a][i]*S[full^a][j]for a in range(full+1))for j in range(n)]for i in range(n)]

class ThreeCover:
 def __init__(self,rows):
  self.rows=rows;self.n=len(rows);self.full=(1<<self.n)-1;self.E=endpoints(rows);rev=[sum(1<<j for j in range(self.n)if rows[j]>>i&1)for i in range(self.n)];self.S=endpoints(rev);self.J={};self.cache={}
 def D(self,u,v,w,z):
  key=(u,v,w,z)
  if key in self.cache:return self.cache[key]
  if v not in self.J:self.J[v]=endpoints(self.rows,v)
  free=self.full^sum(1<<i for i in key);b=free;ans=0
  while True:
   middle=b|(1<<v)|(1<<w);j=self.J[v][middle][w]
   if j:
    rest=free^b;a=rest
    while True:
     ans+=self.E[a|(1<<u)][u]*j*self.S[(rest^a)|(1<<z)][z]
     if not a:break
     a=(a-1)&rest
   if not b:break
   b=(b-1)&free
  self.cache[key]=ans;return ans
 def quartic(self,p,q,r,i):
  return sum(self.D(u,v,w,z)+self.D(w,z,u,v)for u,v in[(p,q),(q,p)]for w,z in[(r,i),(i,r)])

def load_bank():
 paths=sorted(HERE.glob('bank.part*.json'))
 if paths:return sum([json.loads(p.read_text())for p in paths],[])
 return json.loads((HERE/'bank.json').read_text())

def main():
 if len(sys.argv)>1 and sys.argv[1]=='split':
  raw=json.loads((HERE/(sys.argv[2]if len(sys.argv)>2 else'bank.json')).read_text())
  for j in range(0,len(raw),20):(HERE/f'bank.part{j//20}.json').write_text(json.dumps(raw[j:j+20],separators=(',',':'))+'\n')
  print(json.dumps({'parts':(len(raw)+19)//20,'graphs':len(raw)}));return
 inner=[]
 for n in range(3,9):
  keys,X,Y=literal_inner(n)
  for key in keys:need(derivative(X[key])==Y[key],f'derivative kernel n={n},key={key}')
  for r in range(2,n):need(X[0,1,r]==formula(n,r)+[0],f'closed Eulerian kernel n={n},r={r}')
  inner.append({'m':n,'keys':keys,'X':[X[k]for k in keys],'Y':[Y[k]for k in keys]})
 (HERE/'inner_coefficients.json').write_text(json.dumps(inner,separators=(',',':'))+'\n')
 algebra=[]
 for m in range(3,32):
  out=[0]*(m-2)
  for r,w in weights(m).items():
   for j,z in enumerate(formula(m,r)):out[j]+=w*z
  want=[2]+([0]*(m-3)if m%2 else[4]+[0]*(m-4))
  need(out==want,'odd/even binomial kernel')
  algebra.append({'m':m,'weights':weights(m),'combined_X':out,'combined_derivative':derivative(out)})
 (HERE/'binomial_controls.json').write_text(json.dumps(algebra,indent=2)+'\n')
 controls=json.loads((HERE/'controls.json').read_text());bank=load_bank();need(len(bank)==len(controls['graphs']),'missing graph')
 for B,G in zip(bank,controls['graphs']):
  need(B['id']==G['id']and B['rows']==G['rows'],'graph identity')
  n=len(B['rows']);expected=[0]+[1<<i for i in range(n)]+[s for s in range(1<<n)if s.bit_count()==2];need(B['masks']==expected,'canonical full named basis')
  L=len(expected);need(len(B['coefficients'])==L*L and len(B['values'])==L*L,'incomplete bank')
  ix={a:i for i,a in enumerate(expected)}
  for ia,a in enumerate(expected):
   for ic,c in enumerate(expected):
    value=0;s=a
    while True:
     t=c
     while True:
      value+=B['coefficients'][ix[s]*L+ix[t]]
      if not t:break
      t=(t-1)&c
     if not s:break
     s=(s-1)&a
    need(value==B['values'][ia*L+ic],'Boolean inversion')
 physical=[];literal=[]
 for case in controls['cases']:
  m=case['m'];erows=case['exterior_rows'];e=len(erows);into=case['into_mask'];alpha,beta,tau,S,E,ell,rho=exterior(erows,into);full=(1<<e)-1
  powers=[{0:1}]
  for j in range(m-2):powers.append(mul(powers[-1],tau))
  theta=[mul(mul(alpha,beta),p).get(full,0)for p in powers]
  if case['label']=='strong_C3_TT_m_TT3_TT3':need(theta[3]>0,'missing genuine third excursion')
  Q=[0]+[r<<1 for r in erows]
  for j in range(e):
   if into>>j&1:Q[j+1]|=1
   else:Q[0]|=1<<(j+1)
  C=capacities(Q);B0=bank[case['graph_ids'][0]];B1=bank[case['graph_ids'][1]];masks=B0['masks'];ix={a:i for i,a in enumerate(masks)};L=len(masks)
  records=[]
  for i in range(e):
   Kpre=plus(mul(alpha,S[i]),mul(beta,E[i]),-1);gamma=mul(mul(alpha,beta),plus(ell[i],rho[i],-1));K=C[0][i+1]-C[i+1][0];need(Kpre.get(full,0)==K,'contracted capacity sign')
   contrasts={}
   for r in range(2,m):
    col=ix[3]*L+ix[(1<<r)|(1<<(m+i))];delta=B0['coefficients'][col]-B1['coefficients'][col];phi=formula(m,r)
    pred=sum(x*mul(Kpre,powers[j]).get(full,0)for j,x in enumerate(phi))+sum(x*mul(gamma,powers[j]).get(full,0)for j,x in enumerate(derivative(phi)))
    need(delta==pred,'physical two-channel quartic identity');contrasts[r]=delta
   combo=sum(w*contrasts[r]for r,w in weights(m).items());expected=2*K if m%2 else 2*K+4*mul(tau,Kpre).get(full,0)+4*gamma.get(full,0)
   need(combo==expected,'physical odd/even binomial identity')
   records.append({'physical_i':m+i,'individual_r_contrasts':contrasts,'K_i':K,'binomial_contrast':combo,'expected':expected})
  record={'case':case['id'],'m':m,'theta':theta,'records':records}
  if 'witness_groups'in case:
   X,Y=case['witness_groups'];KX=sum(records[i-m]['K_i']for i in X);KY=sum(records[i-m]['K_i']for i in Y);need(KX>KY,'ordered exterior skew')
   BQ=bank[case['graph_ids'][2]];grouped=[]
   for group in [X,Y]:
    local=target=0
    for r,w in weights(m).items():
     for i in group:
      col=ix[3]*L+ix[(1<<r)|(1<<i)];local+=w*(B0['coefficients'][col]-B1['coefficients'][col]);target+=w*(B0['coefficients'][col]-BQ['coefficients'][col])
    need(local==target,'whole quotient/all-cover restoration control');grouped.append({'group':group,'local_contrast':local,'quotient_contrast':target})
   need(grouped[0]['quotient_contrast']-grouped[1]['quotient_contrast']==2*(KX-KY)>0,'grouped quotient-only separation')
   record['witness']={'K_X':KX,'K_Y':KY,'grouped':grouped,'positive_difference':2*(KX-KY)}
   # Independent normal form with complete joint support, not attachment DP.
   for gid in case['graph_ids'][:2]:
    B=bank[gid];engine=ThreeCover(B['rows']);checks=[]
    for r in weights(m):
     for i in [X[0],Y[0]]:
      exact=engine.quartic(0,1,r,i);col=ix[3]*L+ix[(1<<r)|(1<<i)];need(exact==B['coefficients'][col],'independent joint three-cover check');checks.append({'r':r,'i':i,'coefficient':exact})
    literal.append({'graph_id':gid,'checks':checks})
  physical.append(record)
 (HERE/'physical_receipt.json').write_text(json.dumps({'cases':physical,'independent_joint_three_cover_controls':literal},indent=2)+'\n')
 files=['controls.json','graphs.txt','inner_coefficients.json','binomial_controls.json','physical_receipt.json']+[p.name for p in sorted(HERE.glob('bank.part*.json'))]
 receipt={'status':'PASS','literal_inner_orders':list(range(3,9)),'symbolic_binomial_orders':list(range(3,32)),'physical_cases':len(physical),'full_named_graphs':len(bank),'complete_named_coordinates':sum(len(b['coefficients'])for b in bank),'old_orders':sorted({len(b['rows'])for b in bank}),'all_named_cubic_quartic_zeros_retained':True,'joint_three_cover_graphs':len(literal),'joint_three_cover_coordinates':sum(len(z['checks'])for z in literal),'genuine_higher_excursion_cases':[z['case']for z in physical if len(z['theta'])>3 and z['theta'][3]>0],'quotient_only_controls':[z['case']for z in physical if'witness'in z],'data_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest()for f in files},'scope':'Exact identities with an ordered-witness root-family corollary; no arbitrary-prime closure.'}
 print(json.dumps(receipt,indent=2,sort_keys=True))
if __name__=='__main__':main()
