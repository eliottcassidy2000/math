"""Independent U-family audit. No imports from producer or response checkers."""
from pathlib import Path
from itertools import product
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'staircase_prime_quotient'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def induce(T,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs)
def endpoint(T):
 n=len(T);N=1<<n;P=op(T);D=[[0]*n for _ in range(N)]
 for u in range(n):D[1<<u][u]=1
 for M in range(1,N):
  rem=M
  while rem:
   bit=rem&-rem;rem-=bit;v=bit.bit_length()-1;S=M^bit
   if S:
    pred=S&P[v];z=0
    while pred:
     b=pred&-pred;pred-=b;z+=D[S][b.bit_length()-1]
    D[M][v]=z
 return D
def hp(T):return sum(endpoint(T)[-1]) if T else 1
def extend(T,a,c):
 n=len(T);out=[]
 for u,row in enumerate(T):out.append(row|((bool(a>>u&1))<<n)|((bool(c>>u&1))<<(n+1)))
 out.append(((1<<n)-1)^a|1<<(n+1));out.append(((1<<n)-1)^c);return tuple(out)
def uquot(r):
 n=2*r+1
 return tuple(sum(1<<v for v in range(n) if u!=v and (((v-u)%n<=r) ^ (u<r and v<r))) for u in range(n))
def mod(T,M):return all((T[v]&M) in (0,M) for v in range(len(T)) if not M>>v&1)
def sides(r,z):
 if z<r:return [r+z],[r+z+1]
 j=z-r
 if j==0:return [v for v in range(2*r+1) if v not in (z,0)],[0]
 if j==r:return [r-1],[v for v in range(2*r+1) if v not in (z,r-1)]
 return [j],[j-1]
def inflate(Q,children):
 sizes=list(map(len,children));off=[sum(sizes[:i]) for i in range(len(sizes)+1)]
 which=[i for i,m in enumerate(sizes) for _ in range(m)];out=[]
 for u,i in enumerate(which):
  out.append(0)
  for v,j in enumerate(which):
   arc=(children[i][u-off[i]]>>(v-off[i])&1) if i==j else (Q[i]>>j&1)
   if arc:out[u]|=1<<v
 return tuple(out),off

def main():
 paths=[SRC/'PROOF.md',SRC/'verify.py',BASE/'marked_pair_quotient'/'PROOF.md',BASE/'marked_pair_quotient'/'verify.py',BASE/'ordered_module_skew'/'PROOF.md']
 pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};data=json.loads((SRC/'normal.json').read_text());structures=[]
 for r in range(1,9):
  T=uquot(r);n=len(T);cycle=list(range(r-1,-1,-1))+list(range(r,2*r+1))
  need(all(T[u]>>v&1 for u,v in zip(cycle,cycle[1:]+cycle[:1])),'spanning directed cycle')
  need(all(not mod(T,M) for M in range(1,(1<<n)-1) if M.bit_count()>1),'independent exhaustive primeness')
  for z in range(n):
   A,B=sides(r,z);M=sum(1<<v for v in A+B)
   need(A and B and z not in A+B and len(set(A+B))==len(A+B),'disjoint named exterior sides')
   need(all(T[u]>>v&1 for u in A for v in B),'physical ordered sides')
   need(all(T[v]&M in (0,M) for v in range(n) if v!=z and not M>>v&1),'all deletion modules')
   structures.append({'r':r,'deleted':z,'ordered_sides':[A,B]})
 results=[];queries=0
 for rec in data['positive_literal_batteries']+data['zero_quartic_boundary_batteries']:
  Q=tuple(rec['quotient']);r=(len(Q)-1)//2;need(Q==uquot(r),'independent cyclic-reversal definition');z=rec['pair_index'];children=rec['children'];T,off=inflate(Q,children);U,_=inflate(op(Q),children)
  need(T==tuple(rec['rows_T']) and U==tuple(rec['rows_U']),'independent physical blow-ups');I,J=sides(r,z);need([I,J]==rec['marked_child_sets'],'independent sides');p,q=range(off[z],off[z+1]);A=[u for i in I for u in range(off[i],off[i+1])];B=[u for i in J for u in range(off[i],off[i+1])];am=sum(1<<u for u in A);cm=sum(1<<u for u in B)
  keep=[u for u in range(len(T)) if u not in (p,q)];idx={v:i for i,v in enumerate(keep)};R=induce(T,keep);E=endpoint(R);S=endpoint(op(R));N=(1<<len(R))-1;Ai=[idx[v] for v in A];Bi=[idx[v] for v in B]
  skew=sum(sum(E[M][u] for u in Bi)*sum(S[N^M][v] for v in Ai)-sum(E[M][u] for u in Ai)*sum(S[N^M][v] for v in Bi) for M in range(1,N));values=[];coefficient=0
  for bits,old in zip(product((0,1),repeat=4),rec['literal_battery']):
   a=(bits[0]<<q)|(am if bits[1] else 0);c=(bits[2]<<p)|(cm if bits[3] else 0);f,g=hp(extend(T,a,c)),hp(extend(U,a,c));sign=(-1)**sum(bits)
   need((a,c,f,g,sign)==(old['a'],old['c'],old['T'],old['U'],old['sign']),'independent physical query equality');coefficient+=sign*(f-g);queries+=1;values.append([a,c,f,g])
  need(coefficient==skew==rec['quartic_coefficient']==rec['cover_skew'],'independent ordered-cover / full-query coefficient')
  need((coefficient>0)==(rec in data['positive_literal_batteries']),'positive versus zero boundary scope');results.append({'r':r,'root':z,'children':children,'quartic':coefficient,'literal_queries':values})
 normal=subprocess.check_output([sys.executable,str(SRC/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SRC/'verify.py')]);need(normal==optimized==(SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer receipt replay')
 need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'frozen producer sources')
 print(json.dumps({'status':'PASS','structural_deletions_through_order17':structures,'actual_queries_per_graph':queries,'independent_batteries':results,'producer_modes_equal':True,'source_sha256':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()
