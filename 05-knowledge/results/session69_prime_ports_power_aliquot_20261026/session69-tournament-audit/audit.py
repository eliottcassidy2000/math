#!/usr/bin/env python3
"""Independent exact Session69 field/action, response compiler and endpoint audit.
No imports from the producer. Every check remains active under python -O.
"""
import itertools as it, collections as co, math, json, hashlib, pathlib, time
OUT=pathlib.Path(__file__).parent
checks=0

def ck(ok, label):
 global checks
 checks+=1
 if not ok: raise RuntimeError(label)

def fieldmul(x,y):
 bits=0
 for j in range(3):
  if y>>j&1: bits^=x<<j
 for j in range(4,2,-1):
  if bits>>j&1: bits^=11<<(j-3)
 return bits

def fp(x,k):
 r=1
 for _ in range(k):r=fieldmul(r,x)
 return r

def trace(x):return x^fp(x,2)^fp(x,4)

def paley():return tuple(tuple(int(i!=j and (j-i)%7 in (1,2,4)) for j in range(7)) for i in range(7))

def autos(A):
 n=len(A)
 return [p for p in it.permutations(range(n)) if all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(n))]

def ctri(A):
 return {p for p in it.combinations(range(len(A)),3) if all(sum(A[i][j] for j in p)==1 for i in p)}

def ispath(A,p):return all(A[i][j] for i,j in zip(p,p[1:]))

def field_audit():
 A=paley();alpha=[fp(2,j) for j in range(7)];G=autos(A)
 ck(len(set(alpha))==7 and sorted(alpha)==list(range(1,8)), 'field generator')
 ck([i for i,x in enumerate(alpha) if trace(x)==0]==[1,2,4], 'trace orbit')
 ck(all(A[i][j]==int(i!=j and trace(fieldmul(alpha[j],fp(alpha[i],6)))==0) for i in range(7) for j in range(7)), 'trace tournament')
 fieldG={tuple(alpha.index(fieldmul(u,fp(x,2**j))) for x in alpha) for u in range(1,8) for j in range(3)}
 ck(set(G)==fieldG and len(G)==21, 'full field group')
 modules=[]
 for k in range(2,7):
  for M in it.combinations(range(7),k):
   if all(len({A[x][y] for y in M})==1 for x in set(range(7))-set(M)):modules.append(M)
 ck(not modules,'Paley prime')
 cycles=ctri(A)
 planes=[{tuple(sorted((x+t)%7 for x in B)) for t in range(7)} for B in [(1,2,4),(3,5,6)]]
 ck(len(cycles)==14 and planes[0]|planes[1]==cycles and not planes[0]&planes[1], 'two planes exact')
 for P in planes:
  paircount=co.Counter(pair for line in P for pair in it.combinations(line,2))
  ck(len(P)==7 and set(paircount.values())=={1} and len(paircount)==21,'Fano pairs')
 hyperG=[p for p in it.permutations(range(7)) if {tuple(sorted(p[i] for i in C)) for C in cycles}==cycles]
 ck(len(hyperG)==42,'hypergraph full aut')
 ck(sum(ispath(A,p) for p in it.permutations(range(7)))==189,'Paley H')
 coeff={}
 for a,b,c in it.product(range(2),range(8),range(8)):
  w=tuple(a^trace(fieldmul(b,x))^trace(fieldmul(c,fp(x,3))) for x in alpha)
  ck(w not in coeff,'polynomial injective')
  ck(sum(w)%2==a,'parity')
  coeff[w]=(a,b,c)
 ck(len(coeff)==128,'polynomial surjective')
 def cls(a,b,c):
  if not b:return (a,0,int(c!=0))
  r=fieldmul(c,fp(b,4));return (a,1,min(r,fp(r,2),fp(r,4)))
 buckets=co.defaultdict(set)
 for w,abc in coeff.items():buckets[cls(*abc)].add(w)
 orbitdata=[]
 for param, ws in sorted(buckets.items()):
  w=min(ws);orb={tuple(w[p[i]] for i in range(7)) for p in G}
  ck(ws==orb,'coefficient orbit exact')
  st=21//len(ws);ck(st in [1,3,21],'stabilizer')
  orbitdata.append({'parameter':param,'weight':sum(w),'stabilizer':st,'canonical_word':''.join(map(str,w))})
 ck(len(buckets)==12,'12 classes');ck(co.Counter(x['stabilizer'] for x in orbitdata)=={21:2,3:6,1:4},'stabilizer strata')
 for w,(a,b,c) in coeff.items():
  dual=tuple(1-w[(-i)%7] for i in range(7))
  ck(coeff[dual]==(a^1,fp(c,2),fp(b,4)),'converse polynomial')
  for zi,z in enumerate(alpha):
   B=fieldmul(b,z);C=fieldmul(c,fp(z,3));flipped=list(w);flipped[zi]^=1
   aa,bb,cc=coeff[tuple(flipped)]
   ck((aa,fieldmul(bb,z),fieldmul(cc,fp(z,3)))==(a^1,B^1,C^1),'pointed flip')
 pointkeys={(a,min((B,C),(fp(B,2),fp(C,2)),(fp(B,4),fp(C,4)))) for a,B,C in it.product(range(2),range(8),range(8))}
 ck(len(pointkeys)==48,'pointed classes')
 pointst=co.Counter(3 if B in [0,1] and C in [0,1] else 1 for a,(B,C) in pointkeys)
 ck(pointst=={3:8,1:40},'point stabilizers')
 # Exhibit distinct same-stabilizer weight-three classes from the two planes.
 reps=[min(tuple(int(i in line) for i in range(7)) for line in P) for P in planes]
 ck(cls(*coeff[reps[0]])!=cls(*coeff[reps[1]]),'Fano decorations distinct')
 return {'group':21,'hypergraph_group':42,'cyclic_triangles':14,'H':189,'classes':orbitdata,'pointed_classes':48,'point_stabilizers':dict(pointst),'same_stabilizer_weight3_words':[''.join(map(str,w)) for w in reps]}

def child(bit):
 A=[[0]*4 for _ in range(4)]
 for i,j in [(1,2),(2,3),(3,1)]:A[i][j]=1
 for i in range(1,4):A[0][i]=1-bit;A[i][0]=bit
 return tuple(map(tuple,A))

def covers(A):
 n=len(A)
 for perm in it.permutations(range(n)):
  for cut in range(1<<(n-1)):
   parts=[];start=0
   for j in range(n-1):
    if cut>>j&1:parts.append(perm[start:j+1]);start=j+1
   parts.append(perm[start:])
   if all(ispath(A,p) for p in parts):yield tuple(parts)

def local_profiles(A):
 # Independent subset-of-vertices composition rather than permutation/cuts.
 n=len(A);scores=list(map(sum,A));P=co.defaultdict(list)
 for k in range(1,n+1):
  for subset in it.combinations(range(n),k):
   for p in it.permutations(subset):
    if ispath(A,p):P[sum(1<<x for x in subset)].append(p)
 rows=co.Counter()
 def visit(mask, cnt, start, end):
  if mask==(1<<n)-1:rows[cnt,start,end]+=1;return
  rem=((1<<n)-1)^mask;s=rem
  while s:
   for p in P[s]:visit(mask|s,cnt+1,scores[p[0]] if cnt==0 else start,scores[p[-1]])
   s=(s-1)&rem
 visit(0,0,None,None)
 F=co.Counter();L=co.Counter();E=co.Counter()
 for (c,s,t),v in rows.items():F[c]+=v;L[c,s]+=v;E[c,t]+=v
 return F,L,E,rows

def signperm(p):return -1 if sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2 else 1

def minor_poly(A,rm=None,cm=None):
 n=len(A);rows=[i for i in range(n) if i!=rm];cols=[i for i in range(n) if i!=cm];poly=co.Counter()
 for p in it.permutations(range(len(cols))):
  coeff=signperm(p);mask=0
  for k,ci in enumerate(p):
   i,j=rows[k],cols[ci]
   if i!=j:
    if not A[i][j]:coeff=0;break
    coeff=-coeff;mask|=1<<j
  poly[mask]+=coeff
 return {m:c for m,c in poly.items() if c}

def rational_walks(A,cap,first,last):
 # Cramer's rule: W_first,last = x_first Cofactor(last,first)/det(I-A X).
 m=len(A);base=cap+1;N=base**m;pw=[base**i for i in range(m)]
 D=minor_poly(A);ck(D[0]==1,'denominator unit')
 terms=[(tuple(i for i in range(m) if mask>>i&1),sum(pw[i] for i in range(m) if mask>>i&1),v) for mask,v in D.items() if mask]
 inv=[0]*N;inv[0]=1;digits=[None]*N
 for z in range(N):
  c=tuple(z//pw[i]%base for i in range(m));digits[z]=c
  if z:inv[z]=-sum(v*inv[z-offset] for inds,offset,v in terms if all(c[i] for i in inds))
 M=minor_poly(A,last,first);ns=[]
 for mask,v in M.items():
  inds=[i for i in range(m) if mask>>i&1];inds.append(first)
  exp=tuple(inds.count(i) for i in range(m));ns.append((exp,sum(pw[i]*exp[i] for i in range(m)),v*((-1)**(first+last))))
 W=[0]*N
 for z,c in enumerate(digits):
  W[z]=sum(v*inv[z-off] for exp,off,v in ns if all(c[i]>=exp[i] for i in range(m)))
 return W,digits,D

def endpoint_audit():
 p=[local_profiles(child(b)) for b in range(2)]
 ck(p[0][0]==p[1][0]=={1:3,2:18,3:36,4:24},'independent local profile')
 Q=paley();tables={};den=None
 for end in (0,1,3):
  W,digits,den=rational_walks(Q,4,0,end)
  for a,b in it.product(range(2),repeat=2):
   if end==0 and b!=a:continue
   hist=co.Counter()
   for z,c in enumerate(digits):
    if not W[z] or 0 in c:continue
    if end==0:
     scale=W[z]*math.prod(p[0][0][c[i]] for i in range(1,7))
     for (cc,s,t),v in p[a][3].items():
      if cc==c[0]:hist[s-t-1]+=scale*v
    else:
     scale=W[z]*math.prod(p[0][0][c[i]] for i in range(7) if i not in (0,end))
     for (cc,s),v in p[a][1].items():
      if cc!=c[0]:continue
      for (dd,t),u in p[b][2].items():
       if dd==c[end]:hist[s-t-1]+=scale*v*u
   tables[end,a,b]=hist
 ck(tables[0,0,0]==tables[0,1,1],'diagonal equality')
 for end in (1,3):ck(tables[end,0,0]==tables[end,1,1],'same colour equality')
 histograms={}
 for k in range(8):
  h=co.Counter()
  for a,num in enumerate((7-k,k)):
   for d,v in tables[0,a,a].items():h[d]+=num*v
  counts={(0,0):math.comb(7-k,2),(1,1):math.comb(k,2),(0,1):k*(7-k)//2,(1,0):k*(7-k)//2}
  for end in (1,3):
   for (a,b),num in counts.items():
    for d,v in tables[end,a,b].items():h[d]+=num*v
  histograms[k]={d:v for d,v in sorted(h.items()) if v}
  ck(sum(h.values())==5307704039306577787515,'independent rational H')
  ck(histograms[k]==histograms.get(7-k,histograms[k]),'complement flux equality')
 # Producer is consulted only after all independent arithmetic has finished.
 produced=json.loads((OUT/'inputs'/'producer_endpoint_reference.json').read_text())
 ck({str(k):{str(d):v for d,v in h.items()} for k,h in histograms.items()}==produced['histograms'],'full independent rational histogram match')
 for (end,a,b),h in tables.items():
  expected=produced['diagonal'][a] if end==0 else produced['cross'][f'{end},{a},{b}']
  ck({str(d):v for d,v in h.items()}==expected,'independent endpoint channel match')
 return {'method':'Cramer determinant/cofactor reciprocal polynomial, independent of last-slot DP','denominator_monomials':len(den),'denominator':den,'histograms':histograms,'H':5307704039306577787515,'local_profiles':[dict(x[0]) for x in p]}


from functools import lru_cache

def enc(A,order=None):
 n=len(A);p=range(n) if order is None else order
 return sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(it.combinations(range(n),2)))

def dec(n,z):
 A=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(it.combinations(range(n),2)):A[i][j]=z>>k&1;A[j][i]=1-A[i][j]
 return tuple(map(tuple,A))

@lru_cache(None)
def canon(n,z):
 A=dec(n,z)
 return min(enc(A,p) for p in it.permutations(range(n)))

def rootedcanon(A,marks):
 # Lexicographically minimal marked labels are forced by first occurrence.
 flat=tuple(x for pair in marks for x in pair);fixed=tuple(dict.fromkeys(flat));rank={v:i for i,v in enumerate(fixed)}
 rest=[v for v in range(len(A)) if v not in rank]
 return tuple(rank[x] for x in flat),min(enc(A,fixed+tail) for tail in it.permutations(rest))

def target(A,C):
 U=[list(col) for col in zip(*A)]
 for p in C:
  for x,y in zip(p,p[1:]):U[x][y]^=1;U[y][x]^=1
 return tuple(map(tuple,U))

def substitution(Q,children):
 slots=[];loc=[];offset=[]
 for i,A in enumerate(children):
  offset.append(len(slots));slots.extend([i]*len(A));loc.extend(range(len(A)))
 n=len(slots)
 return tuple(tuple(children[slots[x]][loc[x]][loc[y]] if slots[x]==slots[y] else Q[slots[x]][slots[y]] for y in range(n)) for x in range(n)),offset

@lru_cache(None)
def responses(A):
 raw=co.Counter();n=len(A)
 for C in covers(A):
  key=(len(C),rootedcanon(target(A,C),[(p[0],p[-1]) for p in C]))
  raw[key]+=1
 aut=len(autos(A));rows=co.defaultdict(list)
 for (c,(marks,z)),v in raw.items():
  ck(v%aut==0,'local orbit integrality')
  rows[c].append((dec(n,z),tuple(zip(marks[::2],marks[1::2])),v//aut))
 return dict(rows)

def qwalks(Q,content):
 def rec(w,left):
  if not any(left):yield w;return
  for i in range(len(Q)):
   if left[i] and (not w or Q[w[-1]][i]):
    nxt=list(left);nxt[i]-=1;yield from rec(w+(i,),nxt)
 yield from rec((),content)

def gluer(Q,choice,w):
 local=[r[0] for r in choice];U,offs=substitution(tuple(zip(*Q)),local);U=[list(row) for row in U]
 occurrence=[0]*len(Q);ports=[]
 for i in w:
  p=choice[i][1][occurrence[i]];occurrence[i]+=1;ports.append((offs[i]+p[0],offs[i]+p[1]))
 for left,right in zip(ports,ports[1:]):
  a,b=left[1],right[0];U[a][b]^=1;U[b][a]^=1
 return tuple(map(tuple,U))

def compiler_audit():
 one=((0,),);edge=((0,1),(0,0));triangle=((0,1,0),(0,0,1),(1,0,0));prime5=tuple(tuple(int(i!=j and (j-i)%5 in (1,2)) for j in range(5)) for i in range(5))
 palette=[one,edge,triangle,child(0),child(1)]
 controls=[]
 # Every C3 child tuple over the named palette of total order <=6.
 for ids in it.product(range(5),repeat=3):
  if sum(len(palette[i]) for i in ids)<=6:controls.append((triangle,tuple(palette[i] for i in ids),ids))
 # Non-C3 prime actions and all-singleton Paley endpoint/ordinary incidence.
 controls.extend((prime5,tuple(edge if i==j else one for i in range(5)),('prime5',j)) for j in range(5))
 controls.extend([(prime5,(one,)*5,('prime5','singleton')),(paley(),(one,)*7,('Paley7','singleton')),(triangle,(dec(5,10),one,one),('multiplicity5',10))])
 out=[];flux_tests=0;lost_endpoint=None
 for Q,children,ids in controls:
  T,off=substitution(Q,children);n=len(T);G=autos(T)
  childkeys=[(len(A),canon(len(A),enc(A))) for A in children]
  K=sum(all(childkeys[i]==childkeys[g[i]] for i in range(len(Q))) for g in autos(Q))
  ck(len(G)==K*math.prod(len(autos(A)) for A in children),'source automorphism formula')
  raw=co.Counter();endpoint=co.Counter();port_reps=set();scores=list(map(sum,T))
  for path in it.permutations(range(n)):
   if not ispath(T,path):continue
   U=target(T,(path,));Z=canon(n,enc(U));raw[Z]+=1
   endpoint[path[0],path[-1]]+=1
   ck(target(U,(path,))==T,'target involution')
   ck(len(ctri(U))-len(ctri(T))==scores[path[0]]-scores[path[-1]]-1,'literal triangle flux');flux_tests+=1
   port_reps.add(min(tuple(g[v] for v in path) for g in G))
  ck(len(port_reps)*len(G)==sum(raw.values()),'presentation freeness')
  combined=co.Counter();unweighted=co.Counter()
  for content in it.product(*(range(1,len(A)+1) for A in children)):
   for w in qwalks(Q,content):
    for choice in it.product(*(responses(A)[c] for A,c in zip(children,content))):
     U=gluer(Q,choice,w);Z=canon(n,enc(U));combined[Z]+=math.prod(r[2] for r in choice);unweighted[Z]+=1
  ck(all(v%K==0 for v in combined.values()),'top integer normalization')
  expected={z:v//len(G) for z,v in raw.items()}
  ck(all(v%len(G)==0 for v in raw.values()),'raw target normalization')
  ck({z:v//K for z,v in combined.items()}==expected,'independent compressed target-class equality')
  if ids==('multiplicity5',10):ck(dict(unweighted)!=dict(combined),'response multiplicity cannot be discarded')
  # Independent endpoint rational coefficients, all source/end slots, for tiny substitutions.
  if len(Q)==3:
   prof=[local_profiles(A) for A in children]
   predicted=co.Counter()
   for i,j in it.product(range(3),repeat=2):
    W,digits,D=rational_walks(Q,max(map(len,children)),i,j)
    for z,c in enumerate(digits):
     if not W[z] or any(c[h]<1 or c[h]>len(children[h]) for h in range(3)):continue
     if i==j:
      scale=W[z]*math.prod(prof[h][0][c[h]] for h in range(3) if h!=i)
      for (cc,s,t),v in prof[i][3].items():
       if cc==c[i]:predicted[s+sum(len(children[h])*Q[i][h] for h in range(3))-t-sum(len(children[h])*Q[j][h] for h in range(3))-1]+=scale*v
     else:
      scale=W[z]*math.prod(prof[h][0][c[h]] for h in range(3) if h not in (i,j))
      for (cc,s),v in prof[i][1].items():
       if cc!=c[i]:continue
       for (dd,t),u in prof[j][2].items():
        if dd==c[j]:predicted[s+sum(len(children[h])*Q[i][h] for h in range(3))-t-sum(len(children[h])*Q[j][h] for h in range(3))-1]+=scale*v*u
   actual=co.Counter()
   for (i,j),v in endpoint.items():actual[scores[i]-scores[j]-1]+=v
   ck(predicted==actual,'tiny rational endpoint full histogram')
  out.append({'id':ids,'order':n,'H':sum(raw.values()),'Aut':len(G),'K':K,'ordinary_ports':len(port_reps),'ordinary_target_incidence':expected,'unweighted_response_count':sum(unweighted.values()),'weighted_response_count':sum(combined.values())})
 # Scalar path-cover data are identical for A0 and A1 but ordinary target incidence differs in the same prime context.
 scalar_pairs=[]
 for r in out:
  if r['id']==(3,0,0):left=r
  if r['id']==(4,0,0):right=r
 ck(local_profiles(child(0))[0]==local_profiles(child(1))[0],'scalar hostile equal full cover data')
 ck(left['ordinary_target_incidence']!=right['ordinary_target_incidence'],'scalar profiles fail literal target incidence')
 return {'cases':out,'case_count':len(out),'literal_flux_path_checks':flux_tests,'scalar_hostile':{'left':left,'right':right},'multiplicity_hostile':next(r for r in out if r['id']==('multiplicity5',10))}

def loss_hostiles():
 Q=((0,1,0),(0,0,1),(1,0,0));one=((0,),);S=dec(4,7)
 C=((0,1),(3,2));D=((0,2,1),(3,));w=(0,1,2,0)
 U=target(S,C);V=target(S,D)
 ck(ispath(S,C[0]) and ispath(S,C[1]) and ispath(S,D[0]) and ispath(S,D[1]),'hostile covers legal')
 ru=rootedcanon(U,((C[0][0],C[-1][-1]),));rv=rootedcanon(V,((D[0][0],D[-1][-1]),))
 ck(ru==rv==((0,1),4),'outer endpoint target response collision')
 g=[]
 for A,cover in [(U,C),(V,D)]:
  choice=((A,tuple((p[0],p[-1]) for p in cover),1),(one,((0,0),),1),(one,((0,0),),1))
  glued=gluer(Q,choice,w);g.append(canon(6,enc(glued)))
 ck(g==[81,83],'ordered inner-port loss hostile')
 # Response marginal has one equal tuple, while the quotient walk has three choices.
 ck(len(list(qwalks(Q,(1,1,1))))==3 and len(autos(Q))==3,'bare marginal division hostile')
 return {'port_loss':{'child_order':4,'child_code':7,'covers':[C,D],'common_outer_response':ru,'walk':w,'glued_target_codes':g},'bare_marginal':{'child_response_tuples':1,'K':3,'quotient_walks':3,'ordinary_ports':1}}

def walkfamilies(Q,content,d):
 def rec(family,left):
  if not any(left):
   if len(family)==d:yield family
   return
  for i,num in enumerate(left):
   if not num:continue
   nxt=list(left);nxt[i]-=1
   if not family:yield from rec(((i,),),nxt)
   else:
    if Q[family[-1][-1]][i]:yield from rec(family[:-1]+(family[-1]+(i,),),nxt)
    if len(family)<d:yield from rec(family+((i,),),nxt)
 yield from rec((),content)

def gluefamily(Q,choice,family):
 U,offs=substitution(tuple(zip(*Q)),[r[0] for r in choice]);U=[list(r) for r in U]
 used=[0]*len(Q);globalmarks=[]
 for walk in family:
  ports=[]
  for i in walk:
   p=choice[i][1][used[i]];used[i]+=1;ports.append((offs[i]+p[0],offs[i]+p[1]))
  globalmarks.append((ports[0][0],ports[-1][1]))
  for (_,a),(b,_) in zip(ports,ports[1:]):U[a][b]^=1;U[b][a]^=1
 return tuple(map(tuple,U)),globalmarks

def recursive_response_audit():
 Q=((0,1,0),(0,0,1),(1,0,0));one=((0,),);edge=((0,1),(0,0))
 cases=[(Q,(one,one,one)),(Q,(edge,one,one)),(Q,(Q,one,one)),(Q,(edge,edge,one)),(Q,(edge,edge,edge)),(edge,(Q,one)),(edge,(edge,edge))]
 rows=[]
 for Q,children in cases:
  T,off=substitution(Q,children);aut=len(autos(T));n=len(T);keys=[(len(A),canon(len(A),enc(A))) for A in children]
  K=sum(all(keys[i]==keys[p[i]] for i in range(len(Q))) for p in autos(Q))
  expected=co.Counter()
  for C in covers(T):expected[len(C),rootedcanon(target(T,C),tuple((p[0],p[-1]) for p in C))]+=1
  compiled=co.Counter()
  for content in it.product(*(range(1,len(A)+1) for A in children)):
   for d in range(1,n+1):
    for family in walkfamilies(Q,content,d):
     for choice in it.product(*(responses(A)[c] for A,c in zip(children,content))):
      U,marks=gluefamily(Q,choice,family)
      compiled[d,rootedcanon(U,marks)]+=math.prod(r[2] for r in choice)
  ck(all(v%aut==0 for v in expected.values()),'recursive raw orbit divisibility')
  ck(all(v%K==0 for v in compiled.values()),'recursive K divisibility')
  ck({k:v//aut for k,v in expected.items()}=={k:v//K for k,v in compiled.items()},'all-d recursive response identity')
  rows.append({'quotient_order':len(Q),'child_sizes':list(map(len,children)),'order':n,'Aut':aut,'K':K,'all_ordered_covers':sum(expected.values()),'response_classes':len(expected)})
  rows[-1]['by_d']={d:sum(v for (dd,r),v in expected.items() if dd==d)//aut for d in range(1,n+1)}
 return {'cases':rows,'all_ordered_covers':sum(r['all_ordered_covers'] for r in rows),'all_response_classes':sum(r['response_classes'] for r in rows)}

if __name__=='__main__':
 result={'input_pin':'61d7cca5d8df80169f939268269239c2f67ea724','field':field_audit(),'endpoint':endpoint_audit(),'compiler':compiler_audit(),'loss_hostiles':loss_hostiles(),'recursive_response':recursive_response_audit(),'checks':checks}
 (OUT/'audit_results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(json.dumps(result,indent=2,sort_keys=True))
