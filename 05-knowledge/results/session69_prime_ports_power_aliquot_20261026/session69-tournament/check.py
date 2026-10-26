from itertools import permutations,product
from collections import Counter
from math import comb
from functools import reduce
from operator import mul
import json,time

def child(bit):
 a=[[0]*4 for _ in range(4)]
 for i,j in [(1,2),(2,3),(3,1)]:a[i][j]=1
 for i in range(1,4):a[0][i]=1-bit;a[i][0]=bit
 return a

def profiles(a):
 n=len(a);scores=list(map(sum,a));F=[0]*(n+1);S=[Counter() for _ in F];E=[Counter() for _ in F];B=[Counter() for _ in F]
 for p in permutations(range(n)):
  bad=sum((not a[p[i]][p[i+1]])<<i for i in range(n-1))
  for cuts in range(1<<(n-1)):
   if cuts&bad!=bad:continue
   c=1+cuts.bit_count();F[c]+=1;S[c][scores[p[0]]]+=1;E[c][scores[p[-1]]]+=1;B[c][scores[p[0]],scores[p[-1]]]+=1
 return F,S,E,B

def compute():
 prof=[profiles(child(b)) for b in range(2)];F=prof[0][0]
 if F!=prof[1][0]:raise RuntimeError('F mismatch')
 m=7;base=5;pw=[base**i for i in range(m)];N=base**m
 digits=[tuple((z//pw[i])%base for i in range(m)) for z in range(N)]
 Q=[[int(i!=j and (j-i)%7 in (1,2,4)) for j in range(m)]for i in range(m)]
 dp=[[0]*m for _ in range(N)];dp[pw[0]][0]=1
 for z,c in enumerate(digits):
  for last,val in enumerate(dp[z]):
   if not val:continue
   for j in range(m):
    if c[j]<4 and Q[last][j]:dp[z+pw[j]][j]+=val
 diag=[Counter(),Counter()]
 cross={(j,a,b):Counter() for j in (1,3) for a,b in product(range(2),repeat=2)}
 for z,c in enumerate(digits):
  if 0 in c:continue
  allF=reduce(mul,(F[ci] for ci in c),1)
  v=dp[z][0]
  if v:
   scale=v*(allF//F[c[0]])
   for a in range(2):
    for (s,t),count in prof[a][3][c[0]].items():diag[a][s-t-1]+=scale*count
  for j in (1,3):
   v=dp[z][j]
   if not v:continue
   scale=v*(allF//F[c[0]]//F[c[j]])
   for a,b in product(range(2),repeat=2):
    for s,sc in prof[a][1][c[0]].items():
     for t,tc in prof[b][2][c[j]].items():cross[j,a,b][s-t-1]+=scale*sc*tc
 out={}
 for k in range(8):
  hist=Counter()
  for a,count in enumerate([7-k,k]):
   for delta,num in diag[a].items():hist[delta]+=count*num
  pair={(0,0):comb(7-k,2),(1,1):comb(k,2),(0,1):k*(7-k)//2,(1,0):k*(7-k)//2}
  for (a,b),count in pair.items():
   for j in (1,3):
    for delta,num in cross[j,a,b].items():hist[delta]+=count*num
  out[k]={d:hist[d] for d in sorted(hist) if hist[d]}
 return {'profile':F,'histograms':out,'H_by_weight':{k:sum(h.values()) for k,h in out.items()},'diagonal':[{str(d):v for d,v in sorted(h.items())}for h in diag],'cross':{','.join(map(str,k)):{str(d):v for d,v in sorted(h.items())}for k,h in cross.items()}}


def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def fmul(a,b,poly=11):
 r=0;top=1<<(poly.bit_length()-1)
 while b:
  if b&1:r^=a
  a<<=1
  if a&top:a^=poly
  b>>=1
 return r

def fpow(a,n,poly=11):
 r=1
 while n:
  if n&1:r=fmul(r,a,poly)
  a=fmul(a,a,poly);n>>=1
 return r

def trace(a,degree=3,poly=11):
 r=0
 for _ in range(degree):r^=a;a=fmul(a,a,poly)
 return r

def matrix(n,code):
 a=[[0]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):a[i][j]=(code>>k)&1;a[j][i]=1-a[i][j];k+=1
 return a

def code(a,p=None):
 n=len(a);p=range(n) if p is None else p;r=0;k=0
 for i in range(n):
  for j in range(i+1,n):r|=a[p[i]][p[j]]<<k;k+=1
 return r

def autos(a):
 n=len(a)
 return [p for p in permutations(range(n)) if all(a[i][j]==a[p[i]][p[j]] for i in range(n) for j in range(n))]

from functools import lru_cache
@lru_cache(None)
def canonical(n,bits):
 a=matrix(n,bits)
 return min(code(a,p) for p in permutations(range(n)))

def rooted(a,marks):
 best=None
 for p in permutations(range(len(a))):
  inv={v:i for i,v in enumerate(p)};key=(tuple(inv[x] for x in marks),code(a,p))
  if best is None or key<best:best=key
 return best

def ham(a):
 n=len(a)
 def rec(path,left):
  if not left:yield tuple(path);return
  for v in left:
   if not path or a[path[-1]][v]:yield from rec(path+[v],left-{v})
 return list(rec([],set(range(n))))

def hdp(a):
 n=len(a);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for s in range(1<<n):
  for i,v in enumerate(D[s]):
   if v:
    for j in range(n):
     if not(s>>j)&1 and a[i][j]:D[s|(1<<j)][j]+=v
 return sum(D[-1])

def c3(a):
 n=len(a)
 return sum(a[i][j]*a[j][k]*a[k][i]+a[i][k]*a[k][j]*a[j][i] for i in range(n) for j in range(i+1,n)for k in range(j+1,n))

def partner(a,path):
 u=[[a[j][i] for j in range(len(a))]for i in range(len(a))]
 for x,y in zip(path,path[1:]):u[x][y]=a[x][y];u[y][x]=a[y][x]
 return u

def field_check():
 alpha=[fpow(2,i) for i in range(7)];need(len(set(alpha))==7,'F8 primitive')
 Q=[[int(i!=j and (j-i)%7 in (1,2,4))for j in range(7)]for i in range(7)]
 need(all(Q[i][j]==int(trace(fmul(alpha[j],fpow(alpha[i],6)))==0) for i in range(7)for j in range(7)if i!=j),'trace tournament')
 G=autos(Q);aff=[tuple((t+e*i)%7 for i in range(7))for t in range(7)for e in (1,2,4)]
 need(set(G)==set(aff),'full group')
 nontrivial=[]
 for mask in range(1,1<<7):
  S=[i for i in range(7)if mask>>i&1]
  if 1<len(S)<7 and all(len({Q[v][s] for s in S})==1 for v in range(7)if v not in S):nontrivial.append(S)
 need(not nontrivial,'prime')
 fam=[]
 for base in ((1,2,4),(3,5,6)):
  fam.append({tuple(sorted((i+t)%7 for i in base))for t in range(7)})
  counts=Counter(pair for L in fam[-1]for pair in __import__('itertools').combinations(L,2))
  need(set(counts.values())=={1} and len(counts)==21,'Fano')
 cycles={v for v in __import__('itertools').combinations(range(7),3) if sum(Q[i][j] for i in v for j in v)==3 and all(sum(Q[i][j]for j in v)==1 for i in v)}
 need(len(cycles)==14 and fam[0].isdisjoint(fam[1])and fam[0]|fam[1]==cycles,'two Fanos')
 def word(a,b,c):return tuple(a^trace(fmul(b,x)^fmul(c,fpow(x,3)))for x in alpha)
 params={word(a,b,c):(a,b,c)for a,b,c in product(range(2),range(8),range(8))}
 need(len(params)==128,'trace bijection')
 pointset=set();ordinary={};coef_action_checks=0;flip_checks=0
 for w,(a,b,c) in params.items():
  need(a==sum(w)%2,'parity coordinate')
  orbit={tuple(w[p[i]] for i in range(7))for p in G};rep=min(orbit)
  s=21//len(orbit);ordinary.setdefault(rep,{'weight':sum(w),'root_stabilizer':s,'parameters':params[rep]})
  for u in alpha:
   expect=tuple(a^trace(fmul(b,fmul(u,x))^fmul(c,fpow(fmul(u,x),3)))for x in alpha)
   need(word(a,fmul(b,u),fmul(c,fpow(u,3)))==expect,'scale action');coef_action_checks+=1
  need(word(a,fpow(b,4),fpow(c,4))==tuple(w[(2*i)%7]for i in range(7)),'Frobenius action')
  need(word(a^1,fpow(c,2),fpow(b,4))==tuple(1-w[(-i)%7]for i in range(7)),'converse action')
  for iz,z in enumerate(alpha):
   B=fmul(b,z);C=fmul(c,fpow(z,3));point=min((a,fpow(B,2**j),fpow(C,2**j))for j in range(3));pointset.add(point)
   delta=word(1,fpow(z,6),fpow(z,4));need(delta==tuple(int(i==iz)for i in range(7)),'delta polynomial')
   new=w[:iz]+(1-w[iz],)+w[iz+1:];aa,bb,cc=params[new]
   need((aa,fmul(bb,z),fmul(cc,fpow(z,3)))==(a^1,B^1,C^1),'pointed flip');flip_checks+=1
 need(len(ordinary)==12 and Counter(x['root_stabilizer']for x in ordinary.values())=={21:2,3:6,1:4},'ordinary strata')
 need(len(pointset)==48,'pointed count')
 pointedges=set()
 for a,B,C in pointset:
  target=min((a^1,fpow(B^1,2**j),fpow(C^1,2**j))for j in range(3))
  need(target in pointset,'point flip closure');pointedges.add(tuple(sorted(((a,B,C),target))))
 need(len(pointedges)==24,'unoriented flip ports')
 # A hostile to extending the inverse-trace rule from F8 to F32.
 need(len({fpow(2,k,37)for k in range(31)})==31,'F32 primitive control')
 hostile=next((z,fpow(z,30,37),trace(z,5,37))for z in range(2,32)if trace(z,5,37)==trace(fpow(z,30,37),5,37))
 P=ham(Q);need(len(P)==189==hdp(Q),'Paley H')
 prep={min(tuple(p[x]for x in path)for p in G)for path in P};need(len(prep)==9,'Paley presentation count')
 ports=[]
 for p in sorted(prep):
  U=partner(Q,p);need(c3(U)==13,'Paley line flux');need(partner(U,p)==Q,'partner involution')
  target_code=canonical(7,code(U));V=matrix(7,target_code);modules=[]
  for mask in range(1,1<<7):
   S=[i for i in range(7)if mask>>i&1]
   if 1<len(S)<7 and all(len({V[v][t]for t in S})==1 for v in range(7)if v not in S):modules.append(S)
  ports.append({'path':p,'target_canonical_code':target_code,'target_H':hdp(U),'target_Aut':len(autos(U)),'target_proper_nontrivial_modules':modules})
 multiplicity_witnesses=[]
 for i,x in enumerate(ports):
  for y in ports[i+1:]:
   if x['target_canonical_code']!=y['target_canonical_code']:continue
   U=partner(Q,x['path']);V=partner(Q,y['path'])
   isom=next(p for p in permutations(range(7))if all(U[a][b]==V[p[a]][p[b]]for a in range(7)for b in range(7)))
   need(isom not in G,'target merger must not be source automorphism')
   multiplicity_witnesses.append({'source_paths':[x['path'],y['path']],'target_isomorphism':isom,'target_code':x['target_canonical_code']})
 need(len(multiplicity_witnesses)==2,'two literal doubled target classes')
 return {'field_elements_alpha_powers':alpha,'trace_zero_exponents':[i for i,x in enumerate(alpha)if trace(x)==0],'root_Aut':len(G),'C3_count':len(cycles),'ordinary_classes':[dict(word=''.join(map(str,k)),**v)for k,v in sorted(ordinary.items())],'pointed_classes':48,'modular_flip_unoriented_ports':24,'coefficient_scale_checks':coef_action_checks,'pointed_flip_checks':flip_checks,'F32_inverse_trace_hostile':hostile,'Paley_H':189,'Paley_tc':9,'Paley_literal_ports':ports,'Paley_target_merger_witnesses':multiplicity_witnesses,'Paley_distinct_target_classes':len({p['target_canonical_code']for p in ports})}

def covers(a):
 n=len(a);out={c:[]for c in range(1,n+1)}
 for p in permutations(range(n)):
  for cuts in range(1<<(n-1)):
   if any(not a[p[i]][p[i+1]]and not(cuts>>i)&1 for i in range(n-1)):continue
   paths=[];start=0
   for i in range(n-1):
    if cuts>>i&1:paths.append(p[start:i+1]);start=i+1
   paths.append(p[start:]);out[len(paths)].append(tuple(paths))
 return out

def sub(Q,children):
 offsets=[];n=0
 for A in children:offsets.append(n);n+=len(A)
 T=[[0]*n for _ in range(n)]
 for i,A in enumerate(children):
  for j,B in enumerate(children):
   for u in range(len(A)):
    for v in range(len(B)):T[offsets[i]+u][offsets[j]+v]=A[u][v]if i==j else Q[i][j]
 return T,offsets

def local_cover_partner(S,C):
 U=[[S[j][i]for j in range(len(S))]for i in range(len(S))]
 for path in C:
  for x,y in zip(path,path[1:]):U[x][y]=S[x][y];U[y][x]=S[y][x]
 return U

def walks(Q,content):
 out=[]
 def rec(w,c):
  if not any(c):out.append(tuple(w));return
  for i,num in enumerate(c):
   if num and (not w or Q[w[-1]][i]):
    d=list(c);d[i]-=1;rec(w+[i],d)
 rec([],content);return out

def glue(Q,local,marks,w):
 opp=[[Q[j][i]for j in range(len(Q))]for i in range(len(Q))];U,off=sub(opp,local)
 seen=[0]*len(Q);parts=[]
 for i in w:
  r=seen[i];seen[i]+=1;parts.append((off[i]+marks[i][2*r],off[i]+marks[i][2*r+1]))
 for (_,end),(start,_) in zip(parts,parts[1:]):U[end][start]^=1;U[start][end]^=1
 return U

def compiler_check():
 K1=[[0]];T2=[[0,1],[0,0]];C3=[[0,1,0],[0,0,1],[1,0,0]];choices=[K1,T2,C3]
 cv=[covers(S)for S in choices];auts=[len(autos(S))for S in choices];response=[]
 for S,cs,aut in zip(choices,cv,auts):
  r={}
  for c,Cs in cs.items():
   raw=Counter(rooted(local_cover_partner(S,C),tuple(v for p in C for v in (p[0],p[-1])))for C in Cs)
   need(all(v%aut==0 for v in raw.values()),'response integrality')
   r[c]={rho:v//aut for rho,v in raw.items()}
  response.append(r)
 rows=[];totalpaths=0;compressed_cases=0
 for Q in (T2,C3):
  for ids in product(range(3),repeat=len(Q)):
   children=[choices[i]for i in ids];T,off=sub(Q,children);n=len(T);direct=ham(T);need(len(direct)==hdp(T),'source independent H')
   directset=set(direct);produced=set();rawtargets=Counter()
   K=sum(all(ids[p[i]]==ids[i]for i in range(len(Q)))for p in autos(Q));aut=K*reduce(mul,(auts[i]for i in ids),1)
   if n<=6:need(len(autos(T))==aut,'full source group')
   combined=Counter()
   for content in product(*(range(1,len(S)+1)for S in children)):
    Ws=walks(Q,content)
    for w in Ws:
     for Cs in product(*(cv[i][c]for i,c in zip(ids,content))):
      local=[local_cover_partner(S,C)for S,C in zip(children,Cs)];marks=[tuple(v for p in C for v in (p[0],p[-1]))for C in Cs]
      seen=[0]*len(Q);path=[]
      for i in w:
       path.extend(off[i]+v for v in Cs[i][seen[i]]);seen[i]+=1
      path=tuple(path);need(path in directset,'compiler legal path');need(path not in produced,'run bijection uniqueness');produced.add(path)
      U=glue(Q,local,marks,w);need(U==partner(T,path),'literal glue exact')
      if n<=6:rawtargets[canonical(n,code(U))]+=1
     if n<=6:
      for rr in product(*(list(response[i][c].items())for i,c in zip(ids,content))):
       local=[matrix(len(S),rho[1])for S,(rho,count)in zip(children,rr)];marks=[rho[0]for rho,count in rr]
       target=canonical(n,code(glue(Q,local,marks,w)));combined[target]+=reduce(mul,(count for rho,count in rr),1)
   need(produced==directset,'complete compiler universe')
   if n<=6:
    need(all(v%aut==0 for v in rawtargets.values()),'ordinary target raw normalization')
    need(all(v%K==0 for v in combined.values()),'top aggregate normalization')
    need({k:v//aut for k,v in rawtargets.items()}=={k:v//K for k,v in combined.items()},'P2 target incidence convolution');compressed_cases+=1
   rows.append({'quotient_order':len(Q),'child_types':ids,'n':n,'H':len(direct),'Aut':aut,'target_classes':len(rawtargets)if n<=6 else None});totalpaths+=len(direct)
 return {'substitutions':len(rows),'labelled_paths_compiled':totalpaths,'full_compressed_target_convolutions':compressed_cases,'cases':rows}


def port_hostile_check():
 Q=[[0,1,0],[0,0,1],[1,0,0]];S=matrix(4,7);one=[[0]]
 C=((0,1),(3,2));D=((0,2,1),(3,));w=(0,1,2,0)
 us=[local_cover_partner(S,X)for X in (C,D)]
 outer=[rooted(U,(X[0][0],X[-1][-1]))for U,X in zip(us,(C,D))]
 need(outer[0]==outer[1], 'outer endpoint hostile requires same rooted response')
 targets=[]
 for U,X in zip(us,(C,D)):
  marks=[tuple(v for path in X for v in (path[0],path[-1])),(0,0),(0,0)]
  V=glue(Q,[U,one,one],marks,w);targets.append(canonical(6,code(V)))
 need(targets==[81,83],'intermediate component marks hostile')
 S5=matrix(5,10);aut=len(autos(S5));cv=covers(S5);responses={}
 for c,Cs in cv.items():
  raw=Counter(rooted(local_cover_partner(S5,C),tuple(v for path in C for v in (path[0],path[-1])))for C in Cs)
  need(all(v%aut==0 for v in raw.values()),'hostile child orbit integrality')
  responses[c]={rho:v//aut for rho,v in raw.items()}
 weighted=0;unweighted=0
 for c,rr in responses.items():
  num=len(walks(Q,(c,1,1)));weighted+=num*sum(rr.values());unweighted+=num*len(rr)
 T,off=sub(Q,[S5,one,one]);need(hdp(T)==weighted==93 and unweighted==88 and aut==1,'response multiplicities are necessary')
 return {'missing_intermediate_marks':{'child_code':7,'covers':[C,D],'common_outer_rooted_response':outer[0],'quotient_walk':w,'distinct_target_codes':targets},'lost_response_multiplicities':{'child_code':10,'child_Aut':aut,'parent_H':weighted,'wrong_unweighted_total':unweighted,'nonunit_responses':[{'c':c,'marks':rho[0],'local_target_code':rho[1],'R':v}for c,rr in responses.items()for rho,v in sorted(rr.items())if v>1]}}

def walk_tuples(Q,content,d):
 out=[];length=sum(content)
 def rec(word,left):
  if not any(left):
   bad=sum((word[i]==word[i+1] or not Q[word[i]][word[i+1]])<<i for i in range(length-1))
   for cuts in range(1<<(length-1)):
    if cuts.bit_count()!=d-1 or cuts&bad!=bad:continue
    parts=[];start=0
    for i in range(length-1):
     if cuts>>i&1:parts.append(tuple(word[start:i+1]));start=i+1
    parts.append(tuple(word[start:]));out.append(tuple(parts))
   return
  for i,v in enumerate(left):
   if v:
    rest=list(left);rest[i]-=1;rec(word+[i],rest)
 rec([],content);return out

def glue_tuple(Q,local,marks,ww):
 opp=[[Q[j][i]for j in range(len(Q))]for i in range(len(Q))];U,off=sub(opp,local);seen=[0]*len(Q);outer=[]
 for w in ww:
  parts=[]
  for i in w:
   r=seen[i];seen[i]+=1;parts.append((off[i]+marks[i][2*r],off[i]+marks[i][2*r+1]))
  for (_,end),(start,_)in zip(parts,parts[1:]):U[end][start]^=1;U[start][end]^=1
  outer.extend((parts[0][0],parts[-1][1]))
 return U,tuple(outer)

def recursive_check():
 K1=[[0]];T2=[[0,1],[0,0]];C3=[[0,1,0],[0,0,1],[1,0,0]];choices=[K1,T2,C3]
 cases=[(T2,(2,0)),(C3,(0,0,0)),(C3,(1,0,0)),(C3,(2,0,0)),(C3,(1,1,0))]
 response=[]
 for S in choices:
  aut=len(autos(S));rr={}
  for c,Cs in covers(S).items():
   raw=Counter(rooted(local_cover_partner(S,C),tuple(v for path in C for v in(path[0],path[-1])))for C in Cs)
   need(all(v%aut==0 for v in raw.values()),'recursive child integrality');rr[c]={rho:v//aut for rho,v in raw.items()}
  response.append(rr)
 out=[];channels=0;total_cover_orbits=0
 for Q,ids in cases:
  children=[choices[i]for i in ids];T,off=sub(Q,children);n=len(T);aut=len(autos(T));K=sum(all(ids[p[i]]==ids[i]for i in range(len(Q)))for p in autos(Q))
  direct={}
  for d,Cs in covers(T).items():
   raw=Counter(rooted(local_cover_partner(T,C),tuple(v for path in C for v in(path[0],path[-1])))for C in Cs)
   need(all(v%aut==0 for v in raw.values()),'parent cover action free');direct[d]={rho:v//aut for rho,v in raw.items()}
  combined={d:Counter()for d in range(1,n+1)}
  for content in product(*(range(1,len(S)+1)for S in children)):
   for d in range(1,sum(content)+1):
    for ww in walk_tuples(Q,content,d):
     for rr in product(*(list(response[i][c].items())for i,c in zip(ids,content))):
      local=[matrix(len(S),rho[1])for S,(rho,count)in zip(children,rr)];marks=[rho[0]for rho,count in rr]
      U,outer=glue_tuple(Q,local,marks,ww);rho=rooted(U,outer);combined[d][rho]+=reduce(mul,(count for rho,count in rr),1)
  for d,dist in combined.items():
   need(all(v%K==0 for v in dist.values()),'recursive quotient group integrality')
   need({rho:v//K for rho,v in dist.items()}==direct[d],'all-cover response recursion P3')
   channels+=len(dist);total_cover_orbits+=sum(direct[d].values())
  out.append({'quotient_order':len(Q),'child_types':ids,'n':n,'Aut':aut,'response_types_by_d':{d:len(h)for d,h in direct.items()},'cover_orbits_by_d':{d:sum(h.values())for d,h in direct.items()}})
 return {'cases':out,'response_channels':channels,'ordered_cover_orbits':total_cover_orbits}

if __name__=='__main__':
 field=field_check();compiler=compiler_check();hostiles=port_hostile_check();recursive=recursive_check();endpoint=compute()
 need(set(endpoint['H_by_weight'].values())=={5307704039306577787515},'Paley28 common H')
 H=next(iter(endpoint['H_by_weight'].values()))
 need(all(H%(3**7*s)==0 for s in (1,3,21)),'Paley28 port normalization')
 for k,h in endpoint['histograms'].items():
  need(h==endpoint['histograms'][7-k],'converse flux histogram')
 # Endpoint transfer reduction to one-dimensional cut balance: second difference is constant.
 for delta in range(-4,3):
  a=endpoint['histograms'][0].get(delta,0);b=endpoint['histograms'][1].get(delta,0)
  need(all(6*(h.get(delta,0)-a)==k*(7-k)*(b-a)for k,h in endpoint['histograms'].items()),'quadratic weight flux law')
 result={'status':'FINITE-EXACT producer checks; proof audit separate','field_and_ports':field,'literal_compiler':compiler,'port_state_hostiles':hostiles,'recursive_all_cover_response':recursive,'paley28_endpoint_transfer':endpoint,'paley28_tc':{s:H//(3**7*s)for s in (1,3,21)}}
 print(json.dumps(result,sort_keys=True,indent=2))
