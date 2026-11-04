#!/usr/bin/env python3
"""Independent exact controls; Python standard library only."""
import json, math, hashlib, itertools
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def need(b,msg):
 if not b: raise RuntimeError(msg)

def fcoef(n,d,r):
 return (r+1)*sum((-1)**(r-s)*math.comb(r,s)*(s+1)**(n-d-1)*(s**(d-1)+(s+2)**(d-1)) for s in range(r+1))

def literal_distance(n):
 out={(i,j):[0]*(n-1) for i in range(n) for j in range(i+1,n)}
 for p in itertools.permutations(range(n)):
  for cut in range(1,n):
   poly=[1]
   for k in range(n-1):
    if k==cut-1:continue
    nxt=[0]*(len(poly)+1)
    for r,x in enumerate(poly):
     nxt[r+1]+=x
     if p[k]<p[k+1]:nxt[r]+=x
    poly=nxt
   z=out[tuple(sorted((p[cut-1],p[cut])))]
   for r,x in enumerate(poly):z[r]+=x
 return out

def endpoint_tables(rows):
 n=len(rows);N=1<<n
 end=[[0]*n for _ in range(N)]
 for i in range(n):end[1<<i][i]=1
 for S in range(1,N):
  for i,x in enumerate(end[S]):
   if not x:continue
   w=rows[i]&~S
   while w:
    b=w&-w;w-=b;end[S|b][b.bit_length()-1]+=x
 rev=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
 start=[[0]*n for _ in range(N)]
 for i in range(n):start[1<<i][i]=1
 for S in range(1,N):
  for i,x in enumerate(start[S]):
   if not x:continue
   w=rev[i]&~S
   while w:
    b=w&-w;w-=b;start[S|b][b.bit_length()-1]+=x
 return end,start

def capacity(rows):
 n=len(rows);E,S=endpoint_tables(rows);N=(1<<n)-1
 C=[[0]*n for _ in range(n)]
 for mask in range(1,N):
  for j in range(n):
   if E[mask][j]:
    for k in range(n):C[j][k]+=E[mask][j]*S[N^mask][k]
 return C

def delete(rows,i):
 names=[v for v in range(len(rows)) if v!=i]
 return [sum(1<<k for k,v in enumerate(names) if rows[u]>>v&1) for u in names],{v:k for k,v in enumerate(names)}

def shared(rows,m):
 ans=[]
 for i in range(m):
  r,ix=delete(rows,i);C=capacity(r)
  for j in range(m):
   for k in range(j+1,m):
    if i not in (j,k):ans.append([i,j,k,C[ix[j]][ix[k]]+C[ix[k]][ix[j]]])
 return ans

def mul(a,b):
 out={}
 for S,x in a.items():
  for T,y in b.items():
   if not S&T:out[S|T]=out.get(S|T,0)+x*y
 return out

def exterior_theta(rows,into,rmax):
 e=len(rows);alpha={0:1};beta={0:1};tau={}
 for k in range(1,e+1):
  for p in itertools.permutations(range(e),k):
   if all(rows[p[j]]>>p[j+1]&1 for j in range(k-1)):
    s=sum(1<<v for v in p)
    if into>>p[-1]&1:alpha[s]=alpha.get(s,0)+1
    if not(into>>p[0]&1):beta[s]=beta.get(s,0)+1
    if into>>p[-1]&1 and not(into>>p[0]&1):tau[s]=tau.get(s,0)+1
 base=mul(alpha,beta);power={0:1};theta=[]
 for r in range(rmax+1):
  theta.append(mul(base,power).get((1<<e)-1,0));power=mul(power,tau)
 return theta

def make_rows(m,erows,into,order=None):
 order=list(range(m)) if order is None else list(order)
 out=[0]*(m+len(erows))
 for j,u in enumerate(order):
  for v in order[j+1:]:out[u]|=1<<v
 for j,r in enumerate(erows):
  out[m+j]=r<<m
  for i in range(m):
   if into>>j&1:out[m+j]|=1<<i
   else:out[i]|=1<<(m+j)
 return out

def recover(records,m):
 sums={(j,k):0 for j in range(m) for k in range(j+1,m)}
 for i,j,k,x in records:sums[j,k]+=x
 ma=max(sums.values());ends=[p for p,x in sums.items() if x==ma]
 need(len(ends)==1,'endpoints not unique')
 left,right=ends[0]
 def score(u,v):return sums[tuple(sorted((u,v)))]
 interior=[v for v in range(m) if v not in (left,right)]
 interior.sort(key=lambda v:(score(left,v),-score(right,v)))
 return [left]+interior+[right]

def literal_polynomial(rows):
 n=len(rows);x=n;y=n+1;out={}
 for p in itertools.permutations(range(n+2)):
  terms={(0,0):1};valid=True
  for u,v in zip(p,p[1:]):
   if u<n and v<n:
    if not(rows[u]>>v&1):valid=False;break
   elif u>=n and v>=n:
    if not(u==x and v==y):valid=False;break
   else:
    old=u if u<n else v;port=v if u<n else u
    var=(1<<old,0) if port==x else (0,1<<old)
    new={}
    for (a,c),z in terms.items():
     key=(a|var[0],c|var[1])
     new[key]=new.get(key,0)+(z if u<n else -z)
     if u>=n:new[a,c]=new.get((a,c),0)+z
    terms=new
  if valid:
   for key,z in terms.items():out[key]=out.get(key,0)+z
 return out

def main():
 receipt={'scope':'exact finite controls; no collision or unrestricted reconstruction claimed'}
 polys=[]
 for n in range(2,8):
  data=literal_distance(n)
  for (i,j),v in data.items():need(v==[fcoef(n,j-i,r) for r in range(n-1)],f'literal polynomial n={n},i={i},j={j}')
  polys.append({'n':n,'distance_polynomials':[[fcoef(n,d,r) for r in range(n-1)] for d in range(1,n)]})
 for n in range(3,26):
  for d in range(1,n-1):
   diff=[fcoef(n,d+1,r)-fcoef(n,d,r) for r in range(n-1)]
   need(all(x>=0 for x in diff),'coefficient monotonicity')
   need((all(x==0 for x in diff))==(d==1),'strict distance difference')
 (ROOT/'distance_polynomials.json').write_text(json.dumps(polys,separators=(',',':'))+'\n')
 receipt['literal_distance_orders']=list(range(2,8));receipt['positivity_formula_orders']=list(range(3,26))
 census=[]
 for m in range(5,8):
  classes={};records=[]
  for p in itertools.permutations(range(m)):
   pos={v:i for i,v in enumerate(p)};W=lambda d:2 if d<=2 else 2**(d-1)
   data=[]
   for i in range(m):
    for j in range(m):
     for k in range(j+1,m):
      if i in (j,k):continue
      lo,hi=sorted((pos[j],pos[k]));d=hi-lo-int(lo<pos[i]<hi)
      data.append([i,j,k,W(d)])
   recovered=recover(data,m)
   need(tuple(recovered) in (p,p[::-1]),'order reconstruction')
   sig=tuple(z[3] for z in data);classes.setdefault(sig,[]).append(list(p))
   records.append({'order':list(p),'quartics':[z[3] for z in data]})
  need(len(classes)==math.factorial(m)//2,'wrong class count')
  need(all(len(c)==2 and c[0]==c[1][::-1] for c in classes.values()),'nonconverse class')
  chunks=[records] if m<7 else [records[:2520],records[2520:]]
  for chunk_no,chunk in enumerate(chunks):
   fname=f'order_census_{m}.json' if m<7 else f'order_census_{m}.part{chunk_no}.json'
   (ROOT/fname).write_text(json.dumps({'canonical_index':'i ascending, j<k ascending, exclude i=j,k','record_start':0 if m<7 else 2520*chunk_no,'records':chunk},separators=(',',':'))+'\n')
  census.append({'m':m,'orders':len(records),'classes':len(classes),'size_of_every_class':2})
 receipt['complete_physical_order_census']=census
 ambient=[]
 # Exhaust every exterior tournament on at most three names and every root incidence.
 for e in range(4):
  edges=list(itertools.combinations(range(e),2))
  for mask in range(1<<len(edges)):
   erows=[0]*e
   for b,(i,j) in enumerate(edges):
    if mask>>b&1:erows[i]|=1<<j
    else:erows[j]|=1<<i
   for into in range(1<<e):
    theta=exterior_theta(erows,into,2);need(theta[0]>0,'zero constant transfer')
    W={d:sum(theta[r]*fcoef(4,d,r) for r in range(3)) for d in range(1,4)}
    cases=[]
    for order in (list(range(5)),[0,2,1,3,4]):
     rows=make_rows(5,erows,into,order);data=shared(rows,5);pos={v:i for i,v in enumerate(order)}
     for i,j,k,z in data:
      lo,hi=sorted((pos[j],pos[k]));d=hi-lo-int(lo<pos[i]<hi)
      need(z==W[d],'ambient distance transfer')
     need(recover(data,5) in (order,order[::-1]),'ambient recovery')
     cases.append({'order':order,'rows':rows,'shared_quartics':data})
    ambient.append({'exterior_rows':erows,'into_mask':into,'theta':theta,'cases':cases})
 (ROOT/'ambient_controls.json').write_text(json.dumps(ambient,separators=(',',':'))+'\n')
 receipt['exhaustive_ambient_rooted_patterns']=len(ambient);receipt['actual_capacity_cases']=2*len(ambient)
 bank=json.loads((ROOT/'bank.json').read_text());literal_ids=[]
 for z in bank['cases']:
  rows=z['rows'];m=z['m'];masks=z['masks'];L=len(masks);ix={a:i for i,a in enumerate(masks)};co=z['coefficients'];val=z['values']
  need(len(co)==L*L and len(val)==L*L,'incomplete bank')
  for i,a in enumerate(masks):
   for j,c in enumerate(masks):
    reconstructed=sum(co[ii*L+jj] for ii,aa in enumerate(masks) if aa&a==aa for jj,cc in enumerate(masks) if cc&c==cc)
    need(reconstructed==val[i*L+j],'finite-difference inversion')
  data=shared(rows,m)
  for i,j,k,x in data:need(co[ix[(1<<i)|(1<<j)]*L+ix[(1<<i)|(1<<k)]]==x,'full bank shared quartic')
  need(recover(data,m) in (z['order'],z['order'][::-1]),'full bank reconstruction')
  for seed in range(4):
   n=len(rows);full=(1<<n)-1;a=((0x9e3779b1*(seed+1))^(z['id']*37))&full;c=((0x85ebca6b*(seed+3))^(z['id']*53))&full
   aug=list(rows)+[1<<(n+1),0]
   for v in range(n):
    if a>>v&1:aug[v]|=1<<n
    else:aug[n]|=1<<v
    if c>>v&1:aug[v]|=1<<(n+1)
    else:aug[n+1]|=1<<v
   E,_=endpoint_tables(aug);actual=sum(E[-1])
   pred=sum(co[ii*L+jj] for ii,aa in enumerate(masks) if aa&a==aa for jj,cc in enumerate(masks) if cc&c==cc)
   need(pred==actual,'high-weight actual physical query')
  if z['id'] in (0,12):
   literal=literal_polynomial(rows)
   for i,a in enumerate(masks):
    for j,c in enumerate(masks):need(literal.get((a,c),0)==co[i*L+j],f'literal full polynomial case{z["id"]}')
   need(all(a.bit_count()<=2 and c.bit_count()<=2 for (a,c),v in literal.items() if v),'bidegree violation')
   literal_ids.append(z['id'])
 receipt['complete_named_response_banks']={'cases':len(bank['cases']),'old_orders':sorted(set(len(z['rows']) for z in bank['cases'])),'coordinates':sum(len(z['coefficients']) for z in bank['cases']),'literal_word_independent_ids':literal_ids,'all_named_cubic_quartic_zeros_retained':True,'additional_high_weight_queries':4*len(bank['cases'])}
 files=['bank.json','distance_polynomials.json','ambient_controls.json']+[f'order_census_{m}.json' for m in range(5,7)]+[f'order_census_7.part{i}.json' for i in range(2)]
 receipt['data_sha256']={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files}
 print(json.dumps(receipt,indent=2,sort_keys=True))
if __name__=='__main__':main()
