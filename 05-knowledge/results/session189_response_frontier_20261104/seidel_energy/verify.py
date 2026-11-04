#!/usr/bin/env python3
"""Exact rational producer checks; independent physical/run audits are separate."""
from fractions import Fraction as F
from pathlib import Path
import itertools,math,json,hashlib
ROOT=Path(__file__).resolve().parent

def need(ok,*msg):
 if not ok:raise RuntimeError(msg)
def solve(A,b):
 n=len(A);M=[[F(x) for x in row]+[F(y)] for row,y in zip(A,b)]
 for j in range(n):
  choices=[i for i in range(j,n) if M[i][j]];need(choices,'singular matrix')
  p=choices[0];M[j],M[p]=M[p],M[j];v=M[j][j];M[j]=[x/v for x in M[j]]
  for i in range(n):
   if i!=j and M[i][j]:
    v=M[i][j];M[i]=[x-v*y for x,y in zip(M[i],M[j])]
 return [row[-1] for row in M]
def inverse(A):
 n=len(A);columns=[solve(A,[int(i==j) for i in range(n)]) for j in range(n)]
 return [[columns[j][i] for j in range(n)] for i in range(n)]
def bilinear(x,M,y):return sum(x[i]*M[i][j]*y[j] for i in range(len(x)) for j in range(len(y)))
def encode(x):x=F(x);return [x.numerator,x.denominator]
def v2(x):
 x=F(x);need(x,'valuation of zero')
 return (abs(x.numerator)&-abs(x.numerator)).bit_length()-1-((x.denominator&-x.denominator).bit_length()-1)
def mod(x,M):
 x=F(x);need(x.denominator%2,'even reduced denominator');return x.numerator*pow(x.denominator,-1,M)%M

def converse(Q):return tuple(((1<<len(Q))-1)^(1<<i)^r for i,r in enumerate(Q))
def valid(Q):return all(not(r>>i&1) and all(((r>>j)&1)+((Q[j]>>i)&1)==1 for j in range(len(Q)) if i!=j) for i,r in enumerate(Q))
def tournament(n,mask):
 Q=[0]*n;b=0
 for i in range(n):
  for j in range(i+1,n):
   if mask>>b&1:Q[i]|=1<<j
   else:Q[j]|=1<<i
   b+=1
 return tuple(Q)
def carousel(n):return tuple(sum(1<<((i+d)%n) for d in range(1,(n+1)//2)) for i in range(n))
def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def subsets(Q):
 m=len(Q);out=[]
 for mask in range(1<<m):
  ids=[i for i in range(m) if mask>>i&1];n=len(ids)
  S=[[0 if i==j else (1 if Q[i]>>j&1 else -1) for j in ids] for i in ids]
  if not n:r=[]
  elif n%2:r=solve([[S[i][j]+1 for j in range(n)] for i in range(n)],[1]*n)
  else:r=solve(S,[1]*n)
  need(all(x.denominator%2 for x in r),'odd inverse denominators')
  need(sum(r)==n%2,'inverse vector sum')
  if n:
   target=0 if n%2 else 1
   need(all(sum(S[i][j]*r[j] for j in range(n))==target for i in range(n)),'Seidel equation')
  out.append((mask,ids,S,r))
 return out

def square_coeff(ids,r,sign):
 out={0:F(1)}
 for i,x in zip(ids,r):out[1<<i]=sign*4*x+4*x*x
 for a in range(len(ids)):
  for b in range(a+1,len(ids)):out[(1<<ids[a])|(1<<ids[b])]=8*r[a]*r[b]
 return out

def gram_banks(Q,retain=True):
 m=len(Q);half=[x for x in range(1<<m) if x.bit_count()<=2];keys=list(itertools.product(half,repeat=2));position={p:i for i,p in enumerate(keys)}
 source=[F(0)]*len(keys);target=[F(0)]*len(keys);odd=[F(0)]*len(keys);E=F(0);sub=subsets(Q)
 for mask,ids,S,r in sub:
  f=square_coeff(ids,r,-1 if len(ids)%2 else 1)
  opp=f if len(ids)%2 else square_coeff(ids,r,-1)
  for a,x in f.items():
   for c,y in f.items():source[position[a,c]]+=x*y
  for a,x in opp.items():
   for c,y in opp.items():target[position[a,c]]+=x*y
  if len(ids)%2:
   for a,x in f.items():
    for c,y in f.items():odd[position[a,c]]+=x*y
  elif ids:E+=sum(x*x for x in r)
 scale=F(1,1<<m);source=[x*scale for x in source];target=[x*scale for x in target];odd=[x*scale for x in odd];delta=[x-y for x,y in zip(source,target)]
 E*=F(64,1<<m)
 need(source[0]==target[0]==1,'Hamilton scalar fixed point')
 need(all(x.denominator%2 for x in source+target),'2-adic integrality')
 cubic=[[F(0)]*m for _ in range(m)]
 for i in range(m):
  for j in range(m):
   if i!=j:cubic[i][j]=delta[position[(1<<i)|(1<<j),1<<i]]
 contraction=sum((1 if Q[i]>>j&1 else -1)*cubic[i][j] for i in range(m) for j in range(m) if i!=j)
 need(E==contraction,'positive energy identity')
 if m>=2:need(E>0 and v2(E)>=0,'positive integral energy')
 need(all(not x for (a,c),x in zip(keys,delta) if a.bit_count()+c.bit_count()==4),'complete quartic cancellation')
 # All ordered covers are encoded by the same full coefficient banks.
 for h in range(1,8):
  evenfactor=1 if h==1 else 2 if h==2 else 0
  oddfactor=(-2)**(h-1)*math.comb(h+1,2)
  for s,t,o,d in zip(source,target,odd,delta):
   sh=evenfactor*(s-o)+oddfactor*o;th=evenfactor*(t-o)+oddfactor*o
   need(sh-th==evenfactor*d,'complete cover contrast')
   need(sh.denominator%2 and th.denominator%2,'complete cover integrality')
 if not retain:return E,cubic,sub
 def value(coeff,a,c):return sum(x for (u,v),x in zip(keys,coeff) if u&a==u and v&c==v)
 witness=None
 if m>=2:
  i,j=min(((i,j) for i in range(m) for j in range(m) if i!=j and cubic[i][j]),key=lambda p:(v2(cubic[p[0]][p[1]]),p))
  L=v2(cubic[i][j])+1;cube=[];fd=F(0)
  for x,y,z in itertools.product((0,1),repeat=3):
   a=(x<<i)|(y<<j);c=z<<i;d=value(delta,a,c);fd+=(-1)**(3-x-y-z)*d
   cube.append(dict(a=a,c=c,source=encode(value(source,a,c)),target=encode(value(target,a,c)),contrast=encode(d),residue=mod(d,1<<L)))
  need(fd==cubic[i][j],'literal eight-corner cubic')
  need(any(v['residue'] for v in cube),'nonzero actual corner')
  witness=dict(i=i,j=j,cubic=encode(cubic[i][j]),sufficient_modulus_exponent=L,corner_values=cube)
 record=dict(rows=list(Q),half_masks=half,coefficient_order='lexicographic (a,c), both masks of weight at most two',source_coefficients=[encode(x) for x in source],converse_coefficients=[encode(x) for x in target],common_odd_subset_coefficients=[encode(x) for x in odd],cover_reconstruction='h=1: E+O; h=2: 2E-6O; h>=3: (-2)^(h-1)*binom(h+1,2)*O, where E=Hamilton-O',complete_shared_cubics=[[encode(x) for x in row] for row in cubic],energy=encode(E),energy_modulus_exponent=None if m<2 else v2(E)+1,witness=witness,subset_vectors=[dict(mask=mask,names=ids,vector=[encode(x) for x in r]) for mask,ids,S,r in sub])
 return record,E,cubic,sub

def check_resolvents(Q,sub):
 total=0
 for mask,ids,S,r in sub:
  n=len(ids)
  if not n:continue
  for t in (F(0),F(1,2),F(1)):
   v=1+2*t;K=inverse([[S[i][j]+v for j in range(n)] for i in range(n)])
   ones=[F(1)]*n
   for a0,c0 in ((0,0),(1,(1<<n)-1),((1<<n)-1,0),(sum(1<<i for i in range(n) if i%2),sum(1<<i for i in range(n) if i%3))):
    a=[F((a0>>i)&1) for i in range(n)];c=[F((c0>>i)&1) for i in range(n)]
    La=1-2*bilinear(ones,K,[x+t for x in a]);Lc=1-2*bilinear(ones,K,[x+t for x in c])
    Ra=1-2*bilinear([1+t-x for x in a],K,ones);Rc=1-2*bilinear([1+t-x for x in c],K,ones)
    xy=1+t-2*bilinear([1+t-x for x in a],K,[x+t for x in c]);yx=t-2*bilinear([1+t-x for x in c],K,[x+t for x in a])
    A=sum(x*y for x,y in zip(a,r));C=sum(x*y for x,y in zip(c,r));expected=(v*(1+2*A)**2*(1+2*C)**2 if n%2==0 else -(1-2*A)**2*(1-2*C)**2/v**3)
    need(La*xy*Rc+Lc*yx*Ra==expected,'direct two-port resolvent')
    total+=1
 return total

def carousel_formula_controls():
 rows=[]
 for n in range(3,20,2):
  k=(n-1)//2;Q=carousel(n);order=[0]+[v for i in range(1,k+1) for v in (k+i,i)];d=[(-1)**(i+1) for i in range(n)]
  need(all((0 if i==j else (1 if Q[order[i]]>>order[j]&1 else -1))==d[i]*d[j]*(0 if i==j else (1 if i<j else -1)) for i in range(n) for j in range(n)),'carousel switching order')
  count=1<<(n-3);cubic_sum=0;second=0
  for mask in range(count):
   z=(n-3)-2*mask.bit_count();cubic_sum+=(z-1)**2*(z+1);second+=z*z
  need(second==(n-3)*count,'independent-sign second moment')
  coefficient=F(64,1<<n)*cubic_sum
  need(coefficient==-8*(n-4) and mod(coefficient,16)==8,'all-size carousel formula control')
  rows.append(dict(order=n,switch_order=order,subsets=count,cubic_sum=cubic_sum,coefficient=encode(coefficient)))
 return rows

def transitive_congruences():
 checks=0
 for L in range(1,11):
  M=1<<L
  for N in (M,2*M,3*M):
   for k in range(1,2*L+7):
    actual=sum((-1)**(k-j)*math.comb(k,j)*pow(j,N,M) for j in range(k+1))%M
    need(actual==pow(-2,k-1,M),'dyadic TT congruence');checks+=1
 for N in range(4,101,4):
  for k in range(1,16):
   actual=sum((-1)**(k-j)*math.comb(k,j)*pow(j,N,16) for j in range(k+1))%16
   need(actual==pow(-2,k-1,16),'TT multiple-four profile');checks+=1
 return checks

def main():
 small=[];resolvent=0
 for n in range(1,5):
  for mask in range(1<<(n*(n-1)//2)):
   Q=tournament(n,mask);E,K,sub=gram_banks(Q,False);resolvent+=check_resolvents(Q,sub)
   small.append(dict(rows=list(Q),energy=encode(E),shared_cubics=[[encode(x) for x in row] for row in K]))
 named=[('edge',transitive(2)),('C3',carousel(3)),('TT4',transitive(4)),('carousel5',carousel(5)),('Paley7',tuple(sum(1<<((i+d)%7) for d in (1,2,4)) for i in range(7))),('carousel7',carousel(7)),('endpoint_flat8',(238,4,240,118,195,146,34,74)),('carousel9',carousel(9)),('carousel11',carousel(11))]
 banks=[]
 for label,Q in named:
  need(valid(Q),'valid named tournament');record,E,K,sub=gram_banks(Q);record['name']=label;banks.append(record)
  if label.startswith('carousel') or label=='C3':need(K[0][(len(Q)-1)//2]==-8*(len(Q)-4),'exact carousel named coefficient')
  if len(Q)<=9:resolvent+=check_resolvents(Q,sub)
 (ROOT/'rational_banks.json').write_text(json.dumps(dict(scope='Complete quotient-named block-uniform geometric response coefficient banks, not full physical banks for large inflations.',records=banks),sort_keys=True,separators=(',',':'))+'\n')
 (ROOT/'small_energy_bank.json').write_text(json.dumps(small,sort_keys=True,separators=(',',':'))+'\n')
 carousel_checks=carousel_formula_controls();(ROOT/'carousel_controls.json').write_text(json.dumps(carousel_checks,sort_keys=True,indent=2)+'\n')
 receipt=dict(status='PASS',scope='Structural all-order proof with exact rational producer controls; independent physical/run audits are separate.',small_labeled_quotients=len(small),named_complete_banks=len(banks),named_coefficient_positions_per_side=sum(len(r['source_coefficients']) for r in banks),direct_resolvent_checks=resolvent,transitive_congruence_checks=transitive_congruences(),carousel_orders=[r['order'] for r in carousel_checks],named_summaries=[dict(name=r['name'],rows=r['rows'],energy=r['energy'],energy_modulus_exponent=r['energy_modulus_exponent'],witness=r['witness']) for r in banks])
 receipt['data_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['rational_banks.json','small_energy_bank.json','carousel_controls.json']}
 print(json.dumps(receipt,sort_keys=True,indent=2))
if __name__=='__main__':main()
