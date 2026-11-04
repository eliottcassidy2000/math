"""Independent exact convention maps for the primary three-family classification."""
from pathlib import Path
import hashlib,json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def matrix(names,arc):return tuple(sum(1<<j for j,v in enumerate(names) if arc(u,v)) for u in names)
def main():
 records=[]
 for r in range(2,16):
  n=2*r+1
  # Primary definitions on printed page 198, kept distinct from packet code.
  source1=matrix(list(range(n)),lambda u,v: 0<(v-u)%n<=r)
  own1=tuple(sum(1<<((u+d)%n) for d in range(1,r+1)) for u in range(n))
  need(source1==own1,'family one convention')
  names=[('a',i) for i in range(r+1)]+[('b',i) for i in range(1,r+1)]
  def arc(u,v):
   a,i=u;b,j=v
   return (a==b=='b' and i<j) or (a=='b' and b=='a' and i<=j) or (a=='a' and b=='b' and i<j) or (a==b=='a' and i>j)
  source2=matrix(names,arc);f=[2*r-i if a=='a' else r-i for a,i in names]
  def stair(u,v):
   if u<r and v<r:return u>v
   if u>=r and v>=r:return u<v
   if u<r<=v:return v-r<=u
   return u-r>v
  own2=matrix(list(range(n)),stair)
  need(sorted(f)==list(range(n)) and all(bool(source2[i]>>j&1)==bool(own2[f[i]]>>f[j]&1) for i in range(n) for j in range(n)),'family two exact label map')
  names3=[('b',0)]+[('a',i) for i in range(1,2*r+1)]
  source3=matrix(names3,lambda u,v: (u[0]==v[0]=='a' and u[1]<v[1]) or(u[0]=='b' and v[0]=='a' and v[1]%2==1) or(u[0]=='a' and v[0]=='b' and u[1]%2==0))
  g=[2*r]+list(range(2*r));own3=matrix(list(range(n)),lambda u,v: (u<2*r and v<2*r and u<v) or(u==2*r and v<2*r and v%2==0) or(u<2*r and v==2*r and u%2==1))
  need(all(bool(source3[i]>>j&1)==bool(own3[g[i]]>>g[j]&1) for i in range(n) for j in range(n)),'family three exact label map')
  records.append({'r':r,'order':n,'family2_map':f,'family3_map':g})
 root=Path(__file__).resolve().parent;paths=[root/'PROOF.md',root/'verify.py']+[root.parent/d/'PROOF.md' for d in ('ordered_module_skew','alternating_prime_quotient','staircase_prime_quotient','marked_pair_quotient')]
 print(json.dumps({'status':'PASS','family_maps':records,'primary_source_url':'https://trotter.math.gatech.edu/papers/83.pdf','primary_definition_page':198,'primary_classification_page':203,'source_sha256':{str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}},indent=2,sort_keys=True))
if __name__=='__main__':main()
