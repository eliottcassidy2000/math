"""Independent separator referee: deletion-split blocks and literal subset costs.

No mathematical producer is imported or executed. Complete two-regular boards
are generated as all row-pair products followed by exact column counts.
"""
from pathlib import Path
from itertools import combinations,product
from functools import lru_cache
from collections import defaultdict,Counter
from hashlib import sha256
import argparse,json,sys
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
DEST=HERE.parent.parent/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
gates=0
def need(ok,msg):
 global gates
 gates+=1
 if not ok:raise ArithmeticError(msg)

def original_lines(points,slopes):
 groups=defaultdict(list)
 for i,(r,c) in enumerate(points):
  for a in slopes:groups[(a,c-a*r)].append(i)
 return [(k,tuple(v)) for k,v in sorted(groups.items()) if len(v)>2]

def components(vertices,edges,removed=None):
 unseen=set(vertices)-{removed};parts=[]
 while unseen:
  part={min(unseen)}
  while True:
   nxt=part|{v for a,b in edges for v in (a,b) if a!=removed and b!=removed and (a in part or b in part)}
   if nxt==part:break
   part=nxt
  parts.append(part);unseen-=part
 return parts

def split_blocks(ls):
 # Recursive articulation deletion, not DFS lowlinks or the producer Tarjan route.
 n=len(ls);active=sorted({p for line in ls for p in line});encode={p:n+i for i,p in enumerate(active)}
 decode={v:k for k,v in encode.items()}
 edges={tuple(sorted((i,encode[p]))) for i,line in enumerate(ls) for p in line}
 vertices={v for e in edges for v in e};blocks=[]
 def recurse(vs,es):
  for v in sorted(vs):
   parts=components(vs,es,v)
   if len(parts)>1:
    for part in parts:
     sub=part|{v};recurse(sub,{e for e in es if set(e)<=sub})
    return
  blocks.append((vs,es))
 for part in components(vertices,edges):recurse(part,{e for e in edges if set(e)<=part})
 need(sum(len(es) for vs,es in blocks)==len(edges),'deletion blocks partition every original incidence edge')
 groups=[set(vs) for vs,es in blocks]
 while True:
  pair=next(((i,j) for i in range(len(groups)) for j in range(i) if any(v<n for v in groups[i]&groups[j])),None)
  if pair is None:break
  i,j=pair;groups[j]|=groups[i];groups.pop(i)
 bags=sorted([(tuple(sorted(v for v in G if v<n)),tuple(sorted(decode[v] for v in G if v>=n))) for G in groups])
 memberships=defaultdict(list)
 for b,(owned,cells) in enumerate(bags):
  need(set(cells)==set().union(*(set(ls[a]) for a in owned)),'whole owned lines determine every bag cell')
  contained=[(vs,es) for vs,es in blocks if set(v for v in vs if v<n)<=set(owned)]
  excess=sum(2*(len(es)-len(vs)+1)-2 for vs,es in contained if len(es)>1)
  local_high=0
  for p in cells:
   touched=[(vs,es) for vs,es in contained if encode[p] in vs]
   need(len(touched)==1,'a cell touches only one ordinary block within its merged bag')
   local_degree=sum(p in ls[a] for a in owned)
   if local_degree>=3:
    local_high+=1;vs,es=touched[0]
    need(len(es)>1 and sum(encode[p] in e for e in es)==local_degree,'each local exceptional cell lies wholly in one nonbridge block')
  need(local_high<=excess,'complete cycle-excess upper bound on the local branching count')
  if excess==0:need(local_high==0,'incidence cactus blocks eliminate all local branching')
  for p in cells:memberships[p].append(b)
 seps={p:bs for p,bs in memberships.items() if len(bs)>1}
 need(sorted(a for owned,cells in bags for a in owned)==list(range(n)),'each complete original line has exactly one owner')
 initial=len(components(vertices,edges))
 articulations={p for p in active if len(components(vertices,edges,encode[p]))>initial}
 need(set(seps)==articulations,'remaining separators are exactly original articulation cells')
 # Bipartite quotient: bag IDs followed by fresh separator IDs.
 node={p:len(bags)+i for i,p in enumerate(seps)}
 qedges={(b,node[p]) for p,bs in seps.items() for b in bs}
 qvs=set(range(len(bags)+len(seps)))
 need(len(qedges)==len(qvs)-len(components(qvs,qedges)),'the actual bag-cell quotient is a forest')
 return bags,seps

def whole_deletion(ls,forced_deleted=(),forced_retained=()):
 active=sorted({p for a in ls for p in a}|set(forced_deleted)|set(forced_retained))
 delete=set(forced_deleted);retain=set(forced_retained)
 need(not delete&retain,'whole original forced states are consistent')
 free=[p for p in active if p not in delete|retain]
 for size in range(len(free)+1):
  for extra in combinations(free,size):
   removed=delete|set(extra)
   if all(len(set(a)-removed)<=2 for a in ls):return len(removed)
 return 10**9

def localized(ls,root_hint=None,inspect=True):
 bags,seps=split_blocks(ls)
 owned=[b[0] for b in bags];cells=[b[1] for b in bags]
 inc=[{p:tuple(j for j,a in enumerate(aa) if p in ls[a]) for p in cc} for aa,cc in bags]
 hs=[sum(len(ii)>=3 for ii in D.values()) for D in inc]
 calls=0;negative=[];parent_high=0;states_checked=0
 @lru_cache(None)
 def cm(p,parent,state):return 1-state+sum(bm(b,p,state) for b in seps[p] if b!=parent)
 @lru_cache(None)
 def bm(b,parent,state):
  nonlocal calls,parent_high,states_checked
  aa,cc=bags[b];free=[p for p in cc if p!=parent]
  unary={p:(cm(p,b,0),cm(p,b,1)) if p in seps else (1,0) for p in free}
  capacity=[2-int(state and parent in ls[a]) for a in aa]
  if parent is not None and len(inc[b][parent])>=3:parent_high+=1
  weight={p:unary[p][0]-unary[p][1] for p in free}
  negative.extend((b,p,w) for p,w in weight.items() if w<0)
  # Literal full bag subsets: no graph compiler and no shared-cell projection.
  best=10**9
  for bits in range(1<<len(free)):
   selected=[p for j,p in enumerate(free) if bits>>j&1]
   if any(sum(j in inc[b][p] for p in selected)>capacity[j] for j in range(len(aa))):continue
   cost=sum(unary[p][int(p in selected)] for p in free)
   best=min(best,cost)
  # Separate high-cell enumeration plus integral residual edge-set oracle.
  high=[p for p in free if len(inc[b][p])>=3];ordinary=[p for p in free if p not in high]
  compiled=10**9;base=sum(unary[p][0] for p in free)
  for hb in range(1<<len(high)):
   chosen=[p for j,p in enumerate(high) if hb>>j&1]
   residual=[capacity[j]-sum(j in inc[b][p] for p in chosen) for j in range(len(aa))]
   if min(residual,default=0)<0:continue
   calls+=1;reward=-10**9
   for eb in range(1<<len(ordinary)):
    kept=[p for j,p in enumerate(ordinary) if eb>>j&1]
    if any(sum(j in inc[b][p] for p in kept)>residual[j] for j in range(len(aa))):continue
    reward=max(reward,sum(weight[p] for p in kept))
   compiled=min(compiled,base-sum(weight[p] for p in chosen)-reward)
  states_checked+=1
  need(compiled==best,'every conditional bag state: weighted integral compiler equals literal owned-line subsets')
  need(best<10**9,'both parent states remain feasible')
  return best
 seen=set();roots=[]
 order=list(range(len(bags)))
 if root_hint is not None:order=[root_hint]+[b for b in order if b!=root_hint]
 for b in order:
  if b in seen:continue
  roots.append(b);front={b}
  while front:
   seen|=front;front={child for current in front for p in cells[current] if p in seps for child in seps[p]}-seen
 total=sum(bm(b,None,0) for b in roots)
 need(calls<=2*sum(2**h for h in hs),'exact theorem graph-call upper bound for the chosen rooting')
 return dict(tau=total,bags=len(bags),local_h=hs,separators={str(p):bs for p,bs in seps.items()},calls=calls,negative_rewards=negative,fixed_high_parent=parent_high,states=states_checked)

def audit_geometry(record,slopes,reroot=True,brute=True):
 points=record['points'];pairs=original_lines(points,slopes);ls=[line for key,line in pairs]
 need([[list(k),list(v)] for k,v in pairs]==record['lines'],'literal coordinate geometry reconstructs every and only overfull line')
 ans=localized(ls)
 for key in ['tau','bags','local_h','separators','calls','fixed_high_parent']:
  need(ans[key]==record[key],'independent default-root result matches '+key)
 need(sorted(map(tuple,record['negative_rewards']))==sorted(ans['negative_rewards']),'all negative continuation rewards retained')
 if brute:need(ans['tau']==whole_deletion(ls),'localized conditional recurrence equals direct whole-board deletion')
 high_total=ans['fixed_high_parent']
 if reroot:
  high_total=0
  for root in range(ans['bags']):
   rr=localized(ls,root)
   need(rr['tau']==ans['tau'],'every bag rooting retains the exact global optimum')
   high_total+=rr['fixed_high_parent']
 return ls,ans,high_total

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--producer',type=Path,default=DEST if HERE.parent.name=='04-computation' else Path('C:/w/continuing13_20260908_no3'))
 args=ap.parse_args()
 path=args.producer/'continuing13_20260908_no3_cell_separators_certificate.json'
 raw=path.read_bytes();J=json.loads(raw)
 need(J['total_gates']==45383,'declared frozen producer gate universe')
 need(sha256(raw).hexdigest()=='418f449018069ae13fbc7c5c121a3a80a518a952667104b012511b89cb23c271','exact frozen producer certificate')
 totals={};totalboards=0
 for n,expected in [(2,1),(3,6),(4,90),(5,2040)]:
  count=0;neg=0;globalH=0;calls=0;hist=Counter()
  pairs=list(combinations(range(n),2))
  for rows in product(pairs,repeat=n):
   if any(sum(c in row for row in rows)!=2 for c in range(n)):continue
   count+=1;pts=[(r,c) for r,row in enumerate(rows) for c in row]
   ls=[v for k,v in original_lines(pts,(1,-1,2))]
   ans=localized(ls);direct=whole_deletion(ls)
   need(ans['tau']==direct,'entire two-regular board: separator messages equal original deletion')
   need(max(ans['local_h'],default=0)==0,'complete small-board local core is empty')
   incidence=Counter(p for a in ls for p in a)
   exceptional={p for p,degree in incidence.items() if degree>=3}
   need(exceptional<=set(map(int,ans['separators'])),'three-direction exceptional articulation corollary in complete small bank')
   neg+=bool(ans['negative_rewards']);globalH+=bool(exceptional);calls+=ans['calls'];hist[str(max(ans['local_h'],default=0))]+=1
  need(count==expected,'complete row-product and column-filter universe count')
  row=dict(boards=count,negative_reward_boards=neg,with_global_exception=globalH,total_oracle_calls=calls,max_local_h_histogram=dict(hist))
  need(row==J['complete_boards'][str(n)],'complete board aggregate certificate')
  totals[str(n)]=row;totalboards+=count
 need(totalboards==2137,'entire declared geometric board universe')
 abstract={}
 for n,expected in [(0,1),(1,1),(2,2),(3,9),(4,96)]:
  shared=list(combinations(range(n),2))+list(combinations(range(n),3))
  count=0;hist=Counter()
  for mask in range(1<<len(shared)):
   chosen=[v for j,v in enumerate(shared) if mask>>j&1]
   covered=Counter(pair for v in chosen for pair in combinations(v,2))
   if any(v>1 for v in covered.values()) or any(sum(a in v for v in chosen)>3 for a in range(n)):continue
   ls=[[] for _ in range(n)]
   for p,v in enumerate(chosen):
    for a in v:ls[a].append(p)
   p=len(chosen)
   for line in ls:
    while len(line)<3:line.append(p);p+=1
   ls=list(map(tuple,ls));ans=localized(ls);count+=1;hist[str(max(ans['local_h'],default=0))]+=1
   need(ans['tau']==whole_deletion(ls),'complete abstract linear triple system matches direct original-cell deletion')
   for root in range(ans['bags']):need(localized(ls,root)['tau']==ans['tau'],'every abstract-system bag root preserves exact conditional optimum')
  need(count==expected,'complete independently enumerated abstract linear-incidence universe')
  row=dict(systems=count,max_local_h_histogram=dict(hist))
  need(row==J['abstract_systems'][str(n)],'complete abstract local-core histogram')
  abstract[str(n)]=row
 for record in J['connected_chain_family']:
  ls,ans,high=audit_geometry(record,(1,-1,2),brute=False)
  k=record['k'];pts=record['points'];need(len(pts)==4*k+3 and len(ls)==2*k+1,'connected family exact cardinalities')
  centres=set(range(k))
  need(all(len(set(a)-centres)<=2 for a in ls),'deleting all actual chain centres attains the upper bound')
  disjoint=[set(line) for key,line in original_lines(pts,(-1,))]
  need(len(disjoint)==k and sum(map(len,disjoint))==len(set().union(*disjoint)),'k pairwise disjoint overfull lines give matching lower bound')
  need(ans['tau']==k and ans['calls']==4*k+1 and all(h==0 for h in ans['local_h']),'complete connected-family formula and zero local branching')
  need(len({p[0] for p in pts})==len(pts)==len({p[1] for p in pts}),'all occupied chain rows and columns contain one point')
 joint=J['geometric_joint_choice'];ls,ans,high=audit_geometry(joint,(1,-1,2))
 for state,wanted in [('00',3),('01',4),('10',4),('11',2)]:
  deleted=[j for j,s in enumerate(state) if s=='1'];retained=[j for j,s in enumerate(state) if s=='0']
  need(whole_deletion(ls,deleted,retained)==wanted==joint['forced_deletion_cost'][state],'joint exceptional deletion table 3,4,4,2')
 naive=J['naive_component_hostile'];ls,ans,high=audit_geometry(naive,(1,-1,2))
 need(ans['bags']==4 and ans['local_h']==[0]*4,'whole cyclic core retained as one bag')
 # Reconstruct the incorrect remove-all-articulations quotient and its cycle rank.
 n=len(ls);active=sorted({p for a in ls for p in a});code={p:n+i for i,p in enumerate(active)}
 edges={(a,code[p]) for a,line in enumerate(ls) for p in line};vertices={v for e in edges for v in e}
 articulation=set(map(int,ans['separators']));removed={code[p] for p in articulation}
 parts=components(vertices-removed,{e for e in edges if not set(e)&removed})
 owner={v:i for i,part in enumerate(parts) for v in part};qe={(owner[a],len(parts)+p) for a,b in edges for p in articulation if b==code[p]}
 qv=set(range(len(parts)))|{len(parts)+p for p in articulation}
 need(len(qe)-len(qv)+len(components(qv,qe))==1,'naive articulation deletion quotient really contains a cycle')
 four=J['four_direction_boundary'];ls,ans,high=audit_geometry(four,(1,-1,2,-2))
 need(high==four['all_root_fixed_high_states']==4,'fixed high-incidence parent is exercised under all rootings')
 need(sorted(ans['local_h'])==[0,0,2] and set(ans['separators'])=={'0','1'},'articulation cells can remain exceptional inside a four-direction bag')
 for name,record in J['inherited_controls'].items():
  pts=[(r,c) for r,row in enumerate(record['rows']) for c in row]
  for mode,target in record['modes'].items():
   slopes=(1,-1,2,-2) if '-2' in mode else (1,-1,2)
   rec=dict(target,points=pts);audit_geometry(rec,slopes)
 core=J['inherited_controls']['first_saturated_core']
 need(core['rows']==[[0,1],[0,1],[2,3],[4,5],[3,5],[2,4]],'actual smallest-dimension saturated core board')
 need(core['modes']['(1, -1, 2)']['local_h']==[1] and core['modes']['(1, -1, 2)']['tau']==4,'dimension-six local core and exact original cost')
 # Non-iff: the subdivided K4 incidence system has cycle excess four but no high cell.
 K4=[[] for _ in range(4)]
 for p,(a,b) in enumerate(combinations(range(4),2)):K4[a].append(p);K4[b].append(p)
 kk=localized(list(map(tuple,K4)))
 need(kk['bags']==1 and kk['local_h']==[0],'noncactus high-cycle block can still have zero exceptional-cell count')
 cert=dict(status='Independent complete finite audit PASS supporting separate analytic proof',
 producer_certificate_sha256=sha256(raw).hexdigest(),complete_boards=totals,complete_abstract_systems=abstract,chain_lengths=list(range(1,9)),smallest_three_direction_two_regular_local_core_dimension=6,
 retained_controls=['joint exceptional choices 3/4/4/2','naive articulation quotient cycle','four-direction local-high articulation and fixed parent','all inherited integral matching and concurrency controls','cycle-excess bound and non-iff cactus boundary'],
 method='Deletion-based recursive block split; literal local subset costs; independent residual integral edge subsets; direct whole original deletion',
 gates=gates,scope='All2137two-regular boards n2..5 plus declared geometric controls; no probabilistic or extremal promotion')
 out=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()+b'\n'
 (DEST/(HERE.stem+'_certificate.json')).write_bytes(out)
 print('INDEPENDENT: deletion-split blocks, complete owned lines, conditional subset minima, weighted residual edge subsets.')
 print('COMPLETE_BANKS 2137 two-regular boards and109 linear triple systems; 8 chains; all geometric hostiles including the first saturated core.')
 print('PASS: shared-cell costs paid once, negative rewards retained, high parent fixed before branching, exact graph-call bound.')
 print('CERTIFICATE_SHA256',sha256(out).hexdigest())
 print('Always-active exact gates:',gates)
if __name__=='__main__':main()
