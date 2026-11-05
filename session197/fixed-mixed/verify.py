from interface import *
import graphs as g
from collections import Counter
import itertools,json
need=g.require

def potential(s,k,side):
 P=sum(v<k for c in s for v in c)
 return 2*k-P+sum(c[0 if side=='R'else 1]==k for c in s)

def compose(A,B):
 out=[]
 for row in B:
  d=Counter()
  for i,x in row:
   for j,y in A[i]:d[j]+=x*y
  out.append(tuple(sorted(d.items())))
 return tuple(out)

def main():
 formal=[]
 for k in range(5):
  ss=shapes(k);entries=0
  for side in ('R','L'):
   for mask in range(1<<k):
    for i,row in enumerate(trans(k,side,mask)):
     need(0<=potential(ss[i],k,side)<=2*k+1,'potential range');diag=0
     for j,c in row:
      need(c>0,'positive coefficient');entries+=1
      need(i==j or potential(ss[j],k,side)>potential(ss[i],k,side),'strict off-diagonal potential')
      if i==j:diag+=c
     expected=sum(c[1 if side=='R'else 0]==k for c in ss[i])
     need(diag==expected and expected<=k+1,'mask-independent diagonal')
  formal.append({'ports':k,'coordinates':len(ss),'transition_entries':entries})
 commutations=0
 for k in range(4):
  for a,b in itertools.product(range(1<<k),repeat=2):
   need(compose(trans(k,'R',a),trans(k,'L',b))==compose(trans(k,'L',b),trans(k,'R',a)),'opposite-side formal commutation');commutations+=1
 profiles=births=0
 for n in range(1,5):
  pairs=list(itertools.combinations(range(n),2))
  for bits in range(1<<len(pairs)):
   rows=[0]*n
   for j,(u,v)in enumerate(pairs):
    if bits>>j&1:u,v=v,u
    rows[u]|=1<<v
   for k in range(min(n,3)+1):
    for ports in itertools.permutations(range(n),k):
     old=brute(rows,ports);profiles+=1
     need(H(old,k)==g.hp(rows),'Hamilton readout')
     for side,mask in itertools.product(('R','L'),range(1<<k)):
      out=physical(rows,list(ports),side,mask)
      need(step(old,k,side,mask)==brute(out,ports),'literal mixed cover transfer');births+=1
 rows=g.SEED[:];ports=g.PORTS[:];mixed=[]
 for n in range(5):
  core=g.clone(g.SEED,7,n+1);t=n;expected=[0]*(len(core)+t)
  for i in range(t):expected[i]=sum(1<<j for j in range(i+1,len(expected)))
  for i,r in enumerate(core):expected[t+i]=r<<t
  mapping=[8+2*(n-1-i)for i in range(n)]+list(range(7))+[7]+[9+2*i for i in range(n)]
  actual=[0]*len(expected)
  for i,r in enumerate(expected):
   for j in range(len(expected)):
    if r>>j&1:actual[mapping[i]]|=1<<mapping[j]
  need(actual==rows,'mixed clock all-edge identity control')
  h=g.hp(rows);need(h==60*3**(n+1)-1,'direct mixed clock count')
  mixed.append({'periods':n,'order':len(rows),'H':h,'ports':ports[:]})
  if n<4:
   for side,mask in [('L',0),('R',7)]:rows,ports=g.birth(rows,ports,side,mask,3)
   need(ports==g.PORTS,'fixed actual ports')
 return {'status':'PASS','formal':formal,'commuting_mask_pairs':commutations,
         'actual_profiles':profiles,'literal_birth_checks':births,'mixed_plus2':mixed}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
