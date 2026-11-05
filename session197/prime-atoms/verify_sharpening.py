#!/usr/bin/env python3
from verify import *
from verify_factor_bridge import star, converse


def compress(rows,r):
 vs=[x for x in range(len(rows)) if x!=r]
 S=induced(rows,vs);cs=components(S)
 O={i for i,x in enumerate(vs) if rows[r]>>x&1}
 require(O,'positive U required')
 first=next(i for i,C in enumerate(cs) if O.intersection(C))
 keep=[x for i,C in enumerate(cs) if i>=first or len(C)>1 for x in C]
 retained=[vs[x] for x in keep]+[r]
 return induced(rows,retained),len(retained)-1,sum(len(C) for C in cs[first+1:])


def check_compression(rows,r,z=None):
 if z is None:z=profiles(rows)[0][r]
 h,e,R,d=z;U=R-d;f=R+e-h-d
 require(U>0,'positive U')
 require(h-e>=d and U>=d+f,'full deletion insertion bound')
 Q,q,suffix=compress(rows,r)
 require(profiles(Q)[0][q]==z,'prefix compression profile')
 require(U>=(1<<suffix),'suffix exponential bound')
 require(len(Q)<=d+U.bit_length()-1+2,'sharp order bound')
 return Q,q


def main():
 inequalities=0;compressions=0
 for n in range(1,7):
  for rows in labelled(n):
   for r,z in enumerate(profiles(rows)[0]):
    h,e,R,d=z;U=R-d
    if U==0:continue
    require(h-e>=d and U>=d+(R+e-h-d),'complete insertion universe')
    inequalities+=1
    if n<=5:check_compression(rows,r,z);compressions+=1
 rng=random.Random(197029)
 for n in [6,7,8,9,10]:
  for _ in range(30):
   rows=[0]*n
   for i,j in itertools.combinations(range(n),2):
    a,b=(i,j) if rng.randrange(2) else (j,i);rows[a]|=1<<b
   for r,z in enumerate(profiles(rows)[0]):
    if z[2]>z[3]:check_compression(rows,r,z);compressions+=1
 # Nontrivial components on both sides of neutral singleton chains.
 def join(blocks):
  n=sum(map(len,blocks));out=[];off=0
  for block in blocks:
   later=((1<<n)-1)^((1<<(off+len(block)))-1)
   out.extend((x<<off)|later for x in block);off+=len(block)
  return out
 S=join([[4,1,2],[0],[0],[0],[4,1,2],[0],[0],[0]])
 for outs in [[6],[6,7],[0,7],[3,10]]:
  T,r=add_apex(S,sum(1<<x for x in outs))
  for _ in range(3):
   perm=list(range(len(T)));rng.shuffle(perm)
   Q=induced(T,perm);q=perm.index(r)
   check_compression(Q,q);compressions+=1
 sources={};partners={}
 for n in range(1,5):
  for rows in labelled(n):
   ps,_=profiles(rows);qs,_=profiles(converse(rows))
   for r,z in enumerate(ps):
    if z[2]>z[3]:sources.setdefault(z,(rows,r))
    g,s,L,t=qs[r];partners.setdefault((g,t,L),(rows,r))
 pairs=actual=0
 for z,(X,r) in sources.items():
  h,e,R,d=z;U=R-d;f=R+e-h-d
  for (g,t,L),(Y,s) in partners.items():
   beta=d*g-f*t
   if beta not in [-1,1]:continue
   require(2*d<=U+beta,'sharp general factor budget')
   NX,nr=check_compression(X,r,z)
   require(len(NX)<=(U+beta)//2+U.bit_length()-1+2,'sharp source certificate')
   if L==1:NY,ns=[0],0
   else:
    require(L-t>0,'nonconstant converse partner')
    NW,ns=check_compression(converse(Y),s)
    NY=converse(NW)
   ng,_,nL,nt=profiles(converse(NY))[0][ns]
   require((ng,nL,nt)==(g,L,t),'partner triple')
   require(len(NY)<=(1 if L==1 else L+1),'sharp partner certificate')
   for m in [1,2]:
    T,_=star(NX,nr,NY,ns,m)
    require(profiles(T)[0][0][0]==U*L*(1<<m)+beta,'sharp graph compiler')
    actual+=1
   pairs+=1
 bases={(1,1):[1],(1,-1):[3,5],(3,1):[6],(3,-1):[9,15],
        (5,1):[11,12,16],(5,-1):[15,39],(7,1):[],(7,-1):[],
        (9,1):[20,24,32,44],(9,-1):[26,27,34,35,41,45,99],
        (11,1):[29,37,47],(11,-1):[],(13,1):[],(13,-1):[]}
 allowed=[];budgets=[]
 for p in [3,5,7,11,13,17,19,23,29]:
  for beta in [-1,1]:
   D=(p+beta)//2
   known=[d for d in range(1,D+1,2) if (d,beta) in bases and p in bases[d,beta]]
   unknown=[d for d in range(1,D+1,2) if (d,beta) not in bases]
   require(known or not unknown,'no absent case uses unknown deletion data')
   if known:allowed.append((p,beta))
   budgets.append({'prime':p,'offset':beta,'deletion_bound':D,
                   'known_witness_deletions':known,'unneeded_unclassified_deletions':unknown})
 require(allowed==[(3,-1),(5,-1),(11,1),(29,1)],'prime amplitude classification through29')
 X=[2,44,49,53,3,17];r=5
 require(profiles(X)[0][r]==(23,4,40,11),'amplitude29 source')
 vals=[]
 for m in range(1,6):
  val=profiles(k_graph(X,r,m))[0][r][0]
  require(val==29*(1<<m)+1,'amplitude29 actual family');vals.append(val)
 print(json.dumps({'status':'PASS','insertion_inequality_root_checks_through_order6':inequalities,
                   'exact_profile_compression_checks':compressions,
                   'sharp_unit_pair_certificates':pairs,'actual_star_graph_checks':actual,
                   'prime_amplitude_budgets':budgets,'allowed_odd_prime_amplitudes_through29':allowed,
                   'amplitude29_positive_actual_counts':vals,
                   'scope':'Sharper exact inverse certificate; full tournament spectrum remains open.'},
                  indent=2,sort_keys=True))


if __name__=='__main__':main()
