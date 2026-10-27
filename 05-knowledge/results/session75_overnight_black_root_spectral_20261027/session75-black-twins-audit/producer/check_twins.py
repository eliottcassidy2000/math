#!/usr/bin/env python3
"""Two exact labelled tournaments only; stdlib; explicit failures survive -O."""
import itertools,functools,json,math,hashlib
OLD=(966,1004,1016,1009,995,193,384,768,544,96)
NEW=(1008,1001,995,997,494,448,640,768,576,48)
OLD_CORE=((2,3,4,6,7,8,9,5,0,1),(3,0,2,5,9,8,7,6,4,1))
NEW_CORE=((0,5,6,7,8,9,4,3,2,1),(2,3,7,9,8,6,5,0,4,1))
def need(x,m):
 if not x:raise RuntimeError(m)
def edge(a,u,v):return bool(a[u]>>v&1)
def opposite(a):return tuple(((1<<len(a))-1)^(1<<i)^r for i,r in enumerate(a))
def hps(a):
 n=len(a)
 def go(p,unused):
  if not unused:yield tuple(p);return
  for v in sorted(unused):
   if not p or edge(a,p[-1],v):yield from go(p+[v],unused-{v})
 yield from go([],set(range(n)))
def reversal_blue(a,p):
 image={p[i]:p[-1-i]for i in range(len(a))}
 return all(edge(a,u,v)==edge(a,image[v],image[u])for u in range(len(a))for v in range(u+1,len(a)))
def rmove(a,p):
 b=list(a)
 for u,v in zip(p,p[1:]):b[u]^=1<<v;b[v]^=1<<u
 return tuple(b)
def strong(a,vs=None):
 vs=set(range(len(a)))if vs is None else set(vs)
 if not vs:return False
 for b in(a,opposite(a)):
  seen={min(vs)}
  while True:
   new=seen|{v for u in seen for v in vs if edge(b,u,v)}
   if new==seen:break
   seen=new
  if seen!=vs:return False
 return True
def core(a):return any(a[v].bit_count()in(0,len(a)-1)and strong(a,set(range(len(a)))-{v})for v in range(len(a)))
def cuts(a):
 n=len(a)
 return {s:sum(edge(a,u,v)for u in range(n)if s>>u&1 for v in range(n)if not(s>>v&1))for s in range(1,(1<<n)-1)}
def maps(a,b,anti=False):
 n=len(a);ds=[r.bit_count()for r in a];es=[r.bit_count()for r in b];groups={d:[i for i in range(n)if ds[i]==d]for d in sorted(set(ds))};targets={d:[i for i in range(n)if es[i]==(n-1-d if anti else d)]for d in groups}
 if any(len(groups[d])!=len(targets[d])for d in groups):return []
 out=[]
 for choices in itertools.product(*(itertools.permutations(targets[d])for d in groups)):
  p=[None]*n
  for d,seq in zip(groups,choices):
   for u,v in zip(groups[d],seq):p[u]=v
  if all(edge(a,u,v)==(edge(b,p[v],p[u])if anti else edge(b,p[u],p[v]))for u in range(n)for v in range(u+1,n)):out.append(tuple(p))
 return out

def odd_cycle_profile(a):
 """F(x3,x5,...): coefficient counts vertex-disjoint directed odd-cycle collections.
 Cycles are directed cyclic words modulo rotation (not reversal); unused vertices allowed.
 Root each simple cycle at its smallest label, count words by subset/end DP.
 """
 n=len(a);full=(1<<n)-1;cyc={}
 for start in range(n):
  states={(1<<start,start):1}
  for length in range(1,n):
   nxt={}
   for (mask,last),count in states.items():
    if length>=3 and length%2 and edge(a,last,start):cyc[mask]=cyc.get(mask,0)+count
    for v in range(start+1,n):
     if not(mask>>v&1)and edge(a,last,v):nxt[mask|1<<v,v]=nxt.get((mask|1<<v,v),0)+count
   states=nxt
  # n is even here, so every relevant odd length was inspected in1..n-1.
 lengths=tuple(range(3,n+1,2));zero=(0,)*len(lengths)
 @functools.lru_cache(None)
 def F(mask):
  if not mask:return {zero:1}
  bit=mask&-mask;out=dict(F(mask^bit));sub=mask
  while sub:
   if sub&bit and sub in cyc:
    j=lengths.index(sub.bit_count())
    for profile,count in F(mask^sub).items():
     q=list(profile);q[j]+=1;q=tuple(q);out[q]=out.get(q,0)+cyc[sub]*count
   sub=(sub-1)&mask
  return out
 profile=F(full)
 return {'cycle_lengths':lengths,'collection_coefficients':[{'multiplicities':p,'count':c}for p,c in sorted(profile.items())],'single_cycle_counts':{str(k):sum(c for m,c in cyc.items()if m.bit_count()==k)for k in lengths},'H_evaluation_at_all_variables_2':sum(c*2**sum(p)for p,c in profile.items())}
def independent_cycle_profile(a):
 # Separate literal-word enumeration and unordered disjoint-packing recursion.
 n=len(a);lengths=tuple(range(3,n+1,2));cycles=[]
 def go(start,p,left):
  if len(p)>=3 and len(p)%2 and edge(a,p[-1],start):cycles.append(tuple(p))
  if len(p)==n-1:return
  for v in sorted(left):
   if edge(a,p[-1],v):go(start,p+[v],left-{v})
 for start in range(n):go(start,[start],set(range(start+1,n)))
 need(len(cycles)==len(set(cycles)),"duplicated directed cycle")
 masks=[sum(1<<v for v in c)for c in cycles];profile={}
 def pack(first,used,p):
  profile[p]=profile.get(p,0)+1
  for i in range(first,len(cycles)):
   if not(masks[i]&used):
    q=list(p);q[lengths.index(len(cycles[i]))]+=1;pack(i+1,used|masks[i],tuple(q))
 pack(0,0,(0,)*len(lengths))
 return profile,cycles

def route(a,pp):
 n=len(a);out=[]
 for j,p in enumerate(pp):
  need(sorted(p)==list(range(n))and all(edge(a,u,v)for u,v in zip(p,p[1:])),'route directedness');need(not reversal_blue(a,p),'blue witness');b=rmove(a,p);need(not reversal_blue(b,p[::-1]),'R blackness at target')
  s=opposite(a)if j%2 else a;q=p[::-1]if j%2 else p;pathpairs={frozenset((u,v))for u,v in zip(q,q[1:])}
  k=tuple(sum(1<<v for v in range(n)if v!=u and edge(s,u,v)==(frozenset((u,v))in pathpairs))for u in range(n));expected=b if j%2 else opposite(b)
  need(k==expected and not reversal_blue(s,q)and not reversal_blue(k,q),'literal K edge')
  bad=next([p[i],p[t],p[-1-t],p[-1-i]]for i in range(n)for t in range(i+1,n)if edge(a,p[i],p[t])!=edge(a,p[-1-t],p[-1-i]))
  out.append({'R_path':p,'R_source_rows':a,'R_target_rows':b,'failed_anti_named_pairs':bad,'K_path':q,'K_source_rows':s,'K_target_rows':k,'R_target_strong':strong(b)});a=b
 return out,a

def analyse(name,a,pp):
 n=len(a);need(strong(a),'source not strong');values=cuts(a);lam=min(values.values());shores=[s for s,c in values.items()if c==lam];delta=min(min(r.bit_count(),n-1-r.bit_count())for r in a);words=list(hps(a));bb=[p for p in words if reversal_blue(a,p)];kk=[p for p in words if not reversal_blue(a,p)];exits=[p for p in kk if not strong(rmove(a,p))];desc=[p for p in words if any(s>>p[0]&1 and not(s>>p[-1]&1)for s in shores)];aut=maps(a,a);anti=maps(a,a,True);moves,end=route(a,pp);need(core(end),'core route failed');prof=odd_cycle_profile(a);need(prof['H_evaluation_at_all_variables_2']==len(words),'odd-cycle identity control');independent,cyclewords=independent_cycle_profile(a);need(independent=={tuple(x['multiplicities']):x['count']for x in prof['collection_coefficients']},'independent cycle/packing disagreement')
 report={'name':name,'rows':a,'scores':[r.bit_count()for r in a],'sorted_scores':sorted(r.bit_count()for r in a),'lambda':lam,'delta':delta,'minimum_shores':shores,'H':len(words),'blue_words':len(bb),'black_words':len(kk),'automorphism_order':len(aut),'tc_H_div_Aut':len(words)//len(aut),'automorphisms':aut,'anti_automorphisms':anti,'anti_automorphism_count':len(anti),'involutive_anti_automorphism_count':sum(all(p[p[i]]==i for i in range(n))for p in anti),'descent_words':len(desc),'blue_descent_words':sum(reversal_blue(a,p)for p in desc),'black_descent_words':sum(not reversal_blue(a,p)for p in desc),'black_nonstrong_neighbors_counted_by_words':len(exits),'black_nonstrong_distance':1 if exits else 2,'black_core_distance':2,'core_route':moves,'odd_cycle_collection_polynomial':prof,'independently_enumerated_directed_odd_cycle_words':cyclewords}
 if exits:report['one_black_nonstrong_route']=route(a,[exits[0]])[0]
 return report

def incoming():
 # Independently transcribe THM4811 E(C3;K1): D=C3[1,C3,1], rootx0=0.
 D=(14,20,24,18,1);m=5;a=[0]*10
 for u in range(m):
  a[u]=D[u]|sum(1<<v for v in range(m,2*m));a[m+u]=sum(1<<(m+v)for v in range(m)if edge(D,v,u))
 a[0]^=1<<m;a[m]^=1;return tuple(a)
def main():
 old=analyse('Phase1: two cyclic R5 blocks, one reversed bridge',OLD,OLD_CORE);new=analyse('Phase2: rooted-converse E10',NEW,NEW_CORE);inc=incoming();iso=maps(NEW,inc);need(iso,'incoming E10 not isomorphic');need(not maps(OLD,NEW),'old and new unexpectedly isomorphic')
 expected=(old['H'],old['blue_words'],old['black_words'],old['automorphism_order'],old['black_nonstrong_distance'],new['H'],new['blue_words'],new['black_words'],new['automorphism_order'],new['black_nonstrong_distance']);need(expected==(657,33,624,1,1,657,99,558,9,2),'suspected counts differ')
 need(cuts(OLD)==cuts(NEW),'named cut functions differ');need(old['scores']==new['scores'],'named degrees differ');need(old['sorted_scores']==new['sorted_scores']and old['lambda']==new['lambda']==1 and old['delta']==new['delta']==2,'compression twin mismatch')
 print(json.dumps({'status':'PASS','scope':'exactly two labelled sources; no expanded search range','old':old,'new':new,'incoming_THM4811':{'commit':'87d61f67673','rows_from_source_recipe':inc,'isomorphism_new_label_to_incoming_label':iso[0],'credit':'Concurrent independent construction/classification; Session75 does not claim first discovery of E10 or the rooted-converse obstruction.'},'same_named_cut_function_all_1022_shores':True,'same_named_cut_function_sha256':hashlib.sha256(json.dumps(cuts(OLD),sort_keys=True).encode()).hexdigest(),'same_full_odd_cycle_collection_polynomial':old['odd_cycle_collection_polynomial']==new['odd_cycle_collection_polynomial']},sort_keys=True,indent=2))
if __name__=='__main__':main()
