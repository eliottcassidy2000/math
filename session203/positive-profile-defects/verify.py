#!/usr/bin/env python3
"""Literal finite controls for the independently proved all-order statements."""
import importlib.util,itertools,json,pathlib,sys
from counting import J,gate,labelled,nine,relabel
from compiler import move,base,inventory,compile_target,COST
ROOT=pathlib.Path(__file__).resolve().parent

def require(ok,detail):
 if not ok:raise RuntimeError(detail)

def child(a,op):
 n=len(a)
 if op=='O':return tuple([x|(1<<n)for x in a]+[0])
 if op=='Z':return tuple([a[0]]+[x|(1<<n)for x in a[1:]]+[1])
 if op=='J':return J(a)
 return gate(a,int(op[-1]))

def coefficients(p):
 h,s,e,L,R,d,C,b,T=p
 return dict(tau=2*T+2*L-C-2*R-b+2*d,g=C-3*L-h,
  f0=R-h+e-d,q0=b-L-d,a0=C-2*L-b,w0=T-b-s,p0=d-e)

def forests(a):
 A=tuple(range(1,len(a)))
 if not A:return {()}
 result=set()
 for perm in itertools.permutations(A):
  for mask in range(1<<(len(A)-1)):
   parts=[];start=0
   for i in range(len(A)):
    if i==len(A)-1 or mask>>i&1:
     part=perm[start:i+1];start=i+1
     if not all(a[x]>>y&1 for x,y in zip(part,part[1:])):break
     parts.append(part)
   else:
    if len(parts)<=3:result.add(tuple(sorted(parts)))
 return result

def path_type(a,path):return (int(bool(a[0]>>path[0]&1)),int(bool(a[path[-1]]&1)))
def forest_coefficient(a,F):
 z=[path_type(a,p)for p in F];ell=len(z)
 if ell==0:return 1
 if ell==1:return 2*(1-z[0][0])
 if ell==2:return 2*(1-z[0][1]*z[1][0]-z[1][1]*z[0][0])
 if ell==3:return 2*sum(z[j][0]*(2-sum(z[i][1]for i in range(3)if i!=j))for j in range(3))
 return 0

def verify():
 totals=dict(rooted_inputs=0,literal_children=0,forests=0,negative_forests=0,cancellation_images=0)
 for n in range(1,6):
  for a in labelled(n):
   totals['rooted_inputs']+=1;p=nine(a);h,s,e,L,R,d,C,b,T=p;z=coefficients(p)
   require(L>=h and R>=h and R>=d and d>=e and h>=s and z['f0']>=0 and C>=2*L and C>=2*R and T>=R and T>=d,('basic',a,p))
   require(z['tau']>=0,('tau',a,p))
   ff=forests(a);totals['forests']+=len(ff)
   require(sum(forest_coefficient(a,F)for F in ff)==z['tau'],('forest expansion',a,p))
   images=set()
   for F in ff:
    c=forest_coefficient(a,F)
    if c>=0:continue
    totals['negative_forests']+=1
    require(c==-2 and len(F)==2 and all(path_type(a,P)==(1,1)for P in F),('negative shape',a,F,c))
    P=F[0];cuts=[i for i in range(len(P)-1)if a[0]>>P[i]&1 and a[P[i+1]]&1]
    require(bool(cuts),('missing cut',a,F));i=cuts[0]+1
    G=tuple(sorted((P[:i],P[i:],F[1])))
    require(G in ff and forest_coefficient(a,G)==2,('image',a,F,G))
    require(G not in images,('injection collision',a,F,G));images.add(G)
    pre=next(Q for Q in G if path_type(a,Q)==(1,0));suf=next(Q for Q in G if path_type(a,Q)==(0,1));other=next(Q for Q in G if path_type(a,Q)==(1,1))
    require(tuple(sorted((pre+suf,other)))==F,('left inverse',a,F,G))
   totals['cancellation_images']+=len(images)
   if n==1:continue
   require(z['w0']>=0,('w0',a,p))
   if a[0]:
    require(z['q0']>=0 and z['a0']>=0,('outneighbor cone',a,p))
    A=tuple(x>>1 for x in a[1:]);terms=0
    for i in range(1,n):
     if a[0]>>i&1:terms+=nine(relabel(A,[i-1]+[j for j in range(n-1)if j!=i-1]))[3]
    require(b-L==terms,('q endpoint identity',a,p,terms))
   for op in ('O','Z','J','K1','K2'):
    t=nine(child(a,op));totals['literal_children']+=1
    require(t==move(p,op),('physical map',a,op,p,t,move(p,op)))
    g1=coefficients(t)['g'];B=base(p);D=base(t)-2*B
    if op=='O':
     require(g1==2*z['g']+3*L+h and base(t)==B,('O identities',a))
     continue
    if op=='Z':gexp=z['tau']+z['f0']+z['q0']+3*z['a0'];dexp=s+z['f0']+2*z['p0']+z['q0']
    if op=='J':gexp=z['tau']+z['f0']+z['q0'];dexp=z['g']+R-e+s
    if op=='K1':gexp=z['tau']+2*(C-L)+2*T-R-h+e;dexp=2*z['q0']+z['a0']+2*s+2*z['w0']+3*R-h-e
    if op=='K2':gexp=2*z['tau']+5*(C-2*L)+12*T-3*h-2*R-2*d+2*e;dexp=h+2*(R-e)+4*(R-d)+2*d+C+2*(C-2*L)+2*T+2*(T-b)
    require(g1==gexp and D==dexp,('defect identities',a,op,p,t,gexp,dexp))
    if op in ('J','Z')and a[0]:require(g1>=0,('JZ cone',a,op))
    if op in ('K1','K2'):require(g1>=0,('gate cone',a,op))
    if s>0:
     if op in ('Z','K1','K2')or z['g']>=0:require(D>0,('positive defect',a,op,p))
     if op=='J'and z['g']>=0:require(D>=2*s,('J 2s defect',a,p))
     if op=='Z':
      zz=coefficients(t);require((t[1],zz['f0'],zz['q0'],zz['p0'])==(s,z['f0'],z['q0'],0),('Z invariants',a,p))
      delta=s+z['f0']+z['q0'];cur=p
      for k in range(1,6):
       cur=move(cur,'Z');require(base(cur)==2**k*(B+delta+z['p0'])-delta,('Z ray',a,k))
 require(totals['negative_forests']>0,('hostile control absent',totals))
 spec=importlib.util.spec_from_file_location('old_compiler',ROOT/'dependencies/COMPILER_BEFORE_CAP.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
 comparisons=[]
 for bound in (139,997,9999):
  ai,aw=inventory(bound);bi,bw=old.inventory(bound)
  require(ai['complete'] and bi['complete'],('incomplete comparison',bound))
  require(set(aw)==set(bw),('cap support difference',bound,sorted(set(aw)^set(bw))))
  for B,(word,p)in aw.items():
   require(2+sum(COST[x]for x in word)<=max(2,B.bit_length()+1),('target cap',B,word))
   require(coefficients(p)['g']>=0,('grammar cone',B,word))
  comparisons.append(dict(bound=bound,new_states=ai['states'],old_states=bi['states'],support=len(aw),new_order_cap=ai['order_cap'],old_order_cap=bi['order_cap']))
 positive=compile_target(137);negative=compile_target(139);incomplete=compile_target(139,1)
 require(positive['status']=='PHYSICAL_START_TWO_WITNESS'and negative['status']=='ABSENT_START_TWO_PROFILE'and incomplete['status']=='INCOMPLETE',('compiler controls',positive,negative,incomplete))
 require(nine(tuple(positive['rows']))==tuple(positive['profile']),('literal target witness',positive))
 return dict(status='PASS',scope='Finite independent controls of all-order proofs; not staggered occupancy.',totals=totals,cap_comparisons=comparisons,compiler_controls=dict(positive=positive['status'],negative=negative['status'],resource=incomplete['status']))

if __name__=='__main__':print(json.dumps(verify(),indent=2,sort_keys=True))
