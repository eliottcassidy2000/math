from interface import *
import json

def require(x,s):
 if not x:raise RuntimeError(s)
def potential(s,k):
 return 2*k-sum(a<k for a,b in s)-sum(b<k for a,b in s)+int(any(a==k for a,b in s))
def main():
 formal=[];edgecount=0
 for k in range(5):
  ss=shapes(k);entries=0
  for mask in range(1<<k):
   for i,tr in enumerate(transitions(k,mask,k)):
    diag=0
    for j,c in tr:
     entries+=1
     require(potential(ss[j],k)>=potential(ss[i],k),('potential',k,mask,i,j))
     require(i==j or potential(ss[j],k)>potential(ss[i],k),('equal level',k,mask,i,j))
     if i==j:diag+=c
    d=sum(b==k for a,b in ss[i]);require(diag==d,('diagonal',k,mask,i,diag,d))
  formal.append({'ports':k,'shapes':len(ss),'transition_entries':entries});edgecount+=entries
 physical=0;profiles=0
 for n in range(1,5):
  edges=list(combinations(range(n),2))
  for bits in range(1<<len(edges)):
   rows=[0]*n
   for z,(u,v)in enumerate(edges):
    if bits>>z&1:u,v=v,u
    rows[u]|=1<<v
   for k in range(min(n,3)+1):
    for ports in permutations(range(n),k):
     old=brute(rows,list(ports));profiles+=1
     require(H(old,k)>0,'positive actual H')
     for mask in range(1<<k):
      rr,pp=physical_birth(rows,list(ports),mask,k)
      require(pp==list(ports),'same physical ports')
      require(step(old,k,mask,k)==brute(rr,pp),('physical transfer',rows,ports,mask))
      physical+=1
 return {'status':'PASS','formal':formal,'formal_transition_entries':edgecount,
         'actual_profiles':profiles,'physical_birth_checks':physical}
physical_birth=physical
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
