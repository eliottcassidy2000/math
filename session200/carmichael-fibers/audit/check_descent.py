from collections import defaultdict
import json
N=5001;t=(1,0,2,1)
def need(q,m):
 if not q:raise RuntimeError(m)
def move(v,a):
 D,E,U,d=v
 return ((D,0,2*U,d),(D-E+U,d,2*U-d,d),(D-E+U,D-E,D+U,D))[a]
S={t};todo=[t]
for s in todo:
 for a in range(3):
  q=move(s,a)
  if q[2]<N and q not in S:S.add(q);todo.append(q)
F=lambda s:s[0]-s[1]+s[2]
G=lambda s:s[0]+2*s[2]
P={(G(s),2*s[2]-s[3]) for s in S}
bal={3}|{G(s) for s in S if G(s)<=N};fib=defaultdict(list)
for s in S:
 if F(s)<=N:fib[F(s)].append(s)
family=set(bal)
for m,h in P:
 k=1
 while m+2*(2**k-1)*h<=N:family.add(m+2*(2**k-1)*h);k+=1
need(family==set(fib),'all-k equality')
nonbalanced=[];miss=[]
for n in sorted(set(fib)-bal):
 ks=[]
 for D,E,U,d in fib[n]:
  v=U-d;k=(v&-v).bit_length()-1;h=v//2**k;m=n-2*(2**k-1)*h
  need(k>=1 and (m,h) in P and 3<=h<m<2*h and 2*m<n,'joint parent')
  need(n<2**(k+1)*h and n>(2**(k+1)-1)*h,'h interval')
  ks.append(k)
 nonbalanced.append(n)
 if min(ks)>1:miss.append(n)
need(miss==[443,617],'one-step misses')
for D,E,U,d in S:
 m=G((D,E,U,d));h=2*U-d
 need(max((h+1)//2,(m+2)//3)<=U<=min((m+h)//4,2*h//3),'parent pruning')
need(set(fib[443])=={(219,1,225,1),(191,1,253,1)},'443')
print(json.dumps({'status':'PASS','target_cap':N,'physical_profiles':len(S),'nonbalanced_targets':nonbalanced,'one_step_exceptions':miss,'F443':sorted(fib[443]),'all_k_compiler_equals_F':True,'parent_pruning_profiles_checked':len(S)},sort_keys=True,indent=2))
