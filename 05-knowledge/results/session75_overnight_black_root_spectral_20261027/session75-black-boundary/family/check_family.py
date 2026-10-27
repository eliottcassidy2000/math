#!/usr/bin/env python3
"""Portable exact verification of rooted-converse black-distance obstruction.
No external imports, no assertions, no input generators from producer packets.
"""
import itertools,json,hashlib

def need(x,msg):
 if not x:raise RuntimeError(msg)
def arc(t,u,v):return bool(t[u] & (1<<v))
def build(a):
 m=a+4;r=m-1;u=[0]*m
 for i in range(a):
  for j in range(i+1,a):u[i]|=1<<j
  for j in range(a,a+3):u[i]|=1<<j
 for i in range(3):u[a+i]|=1<<(a+(i+1)%3);u[a+i]|=1<<r
 for i in range(a):u[r]|=1<<i
 full=(1<<m)-1
 v=[full^(1<<i)^u[i] for i in range(m)]
 t=[x|(full<<m) for x in v]+[x<<m for x in u]
 t[r]^=1<<(m+r);t[m+r]^=1<<r
 return tuple(t),tuple(u),r

def paths(t):
 n=len(t);full=(1<<n)-1
 def go(p,used):
  if used==full:yield tuple(p);return
  z=t[p[-1]]&~used&full
  while z:
   b=z&-z;z-=b;v=b.bit_length()-1
   yield from go(p+[v],used|b)
 for i in range(n):yield from go([i],1<<i)
def anti_failure(t,p):
 n=len(t)
 for i in range(n):
  for j in range(i+1,n):
   if arc(t,p[i],p[j])!=arc(t,p[n-1-j],p[n-1-i]):return (p[i],p[j],p[n-1-j],p[n-1-i])
 return None
def apply(t,p):
 need(sorted(p)==list(range(len(t))),'not permutation')
 need(all(arc(t,u,v) for u,v in zip(p,p[1:])),'not directed')
 z=list(t)
 for u,v in zip(p,p[1:]):z[u]^=1<<v;z[v]^=1<<u
 return tuple(z)
def converse(t):
 f=(1<<len(t))-1
 return tuple(f^(1<<i)^r for i,r in enumerate(t))
def strong(t,vs=None):
 vs=list(range(len(t))) if vs is None else list(vs)
 if not vs:return False
 s=set(vs);root=vs[0]
 for backwards in (False,True):
  seen={root};todo=[root]
  while todo:
   x=todo.pop()
   for y in s-seen:
    if arc(t,y,x) if backwards else arc(t,x,y):seen.add(y);todo.append(y)
  if seen!=s:return False
 return True
def cuts(t):
 n=len(t);f=(1<<n)-1
 return {s:sum((t[v]&(f^s)).bit_count() for v in range(n) if s>>v&1) for s in range(1,f)}
def mincut_score_certificate(t):
 ds=[r.bit_count() for r in t];order=sorted(range(len(t)),key=lambda i:(ds[i],i));prefix=0;vv=[]
 for k,v in enumerate(order[:-1],1):
  prefix+=ds[v];vv.append(prefix-k*(k-1)//2)
 l=min(vv);sizes=[i+1 for i,v in enumerate(vv) if v==l]
 unique=len(sizes)==1 and ds[order[sizes[0]-1]]<ds[order[sizes[0]]]
 shore=sum(1<<v for v in order[:sizes[0]]) if unique else None
 return {'lambda':l,'size_minima':sizes,'unique':unique,'shore':shore,'prefix_cut_minima':vv}
def literal_k(t,p):
 pe={frozenset((u,v)) for u,v in zip(p,p[1:])};n=len(t)
 return tuple(sum(1<<v for v in range(n) if u!=v and (arc(t,u,v)==(frozenset((u,v)) in pe))) for u in range(n))
def route_for(a):
 m=a+4;r=m-1
 p=tuple(range(a+1,-1,-1))+tuple(range(m,2*m))+(r,a+2)
 q=tuple(range(2*m-1,m-1,-1))+tuple(range(a))+(a,a+2,r,a+1)
 return p,q

def check(a,exhaustive=False):
 t,u,r=build(a);m=a+4;n=2*m;shore=((1<<m)-1)<<m
 need(strong(t),'initial not strong');need(strong(u),'U not strong')
 ds=[x.bit_count() for x in t];delta=min(min(d,n-1-d) for d in ds)
 need(delta==2,'delta is not2');mc=mincut_score_certificate(t)
 need(mc['lambda']==1 and mc['unique'] and mc['shore']==shore,'unique min shore failure')
 # The exactly three root-ending U paths are constructed separately.
 endings=[tuple(range(a))+tuple(a+(i+j)%3 for j in range(3))+(r,) for i in range(3)]
 need(all(apply(u,p) for p in endings),'root-end directedness')
 dec=[tuple(x+m for x in p)+q[::-1] for p in endings for q in endings]
 for p in dec:
  apply(t,p);need(anti_failure(t,p) is None,'geodesic path is black')
 p,q=route_for(a);curr=t;moves=[];states=[t]
 for j,w in enumerate((p,q)):
  need(anti_failure(curr,w) is not None,'route step blue')
  nxt=apply(curr,w);need(anti_failure(nxt,w[::-1]) is not None,'reverse route step blue')
  ks=curr if j%2==0 else converse(curr);kw=w if j%2==0 else w[::-1]
  kn=literal_k(ks,kw);expected=converse(nxt) if j%2==0 else nxt
  need(kn==expected,'literal K conversion');need(anti_failure(ks,kw) is not None,'K source blue');need(anti_failure(kn,kw) is not None,'K target blue')
  moves.append({'R_word':w,'R_source':curr,'R_target':nxt,'failed_anti':anti_failure(curr,w),'K_word':kw,'K_source':ks,'K_target':kn})
  curr=nxt;states.append(curr)
 need(strong(states[1]),'intermediate not strong');need(not strong(states[2]),'terminal strong')
 need(sum((curr[v]&((1<<m)-1)).bit_count() for v in range(m,2*m))==0,'U is not terminal zero shore')
 # After step1, induced U and V cannot be converse-isomorphic.
 mid=states[1];su=sum(((mid[v]>>m)&((1<<m)-1)).bit_count()**2 for v in range(m,2*m));sv=sum((m-1-(mid[v]&((1<<m)-1)).bit_count())**2 for v in range(m))
 need(su!=sv and abs(su-sv)==2,'lost block degree-square obstruction')
 report={'a':a,'n':n,'initial_rows':t,'scores':ds,'delta':delta,'minimum_cut':mc,'decreasing_words':dec,'all_decreasing_words_blue':True,'route':moves,'intermediate_strong':True,'terminal_strong':False,'block_score_square_pair':(su,sv),'black_nonstrong_distance':2}
 if exhaustive:
  allp=list(paths(t));good=[p for p in allp if shore>>p[0]&1 and not(shore>>p[-1]&1)]
  need(set(good)==set(dec),'9-word cone incomplete')
  black=[p for p in allp if anti_failure(t,p) is not None]
  need(all(strong(apply(t,p)) for p in black),'one-step black nonstrong exit')
  c=cuts(t);need([s for s,v in c.items() if v==1]==[shore],'explicit mincut mismatch')
  need(sorted(paths(u)),'U has no paths');need(set(p for p in paths(u) if p[-1]==r)==set(endings),'rooted path completeness')
  report['exhaustive']={'all_global_hamilton_paths':len(allp),'black_paths':len(black),'all_black_neighbors_strong':True,'decreasing_paths':len(good),'direct_cuts_checked':len(c)}
 return report

def main():
 cases=[check(a,exhaustive=(a==1)) for a in range(1,21)]
 print(json.dumps({'status':'PASS','scope':'all-order proof separately; exact controls a=1..20, n10 complete Hamilton-word check','cases':cases},sort_keys=True,indent=2))
if __name__=='__main__':main()
