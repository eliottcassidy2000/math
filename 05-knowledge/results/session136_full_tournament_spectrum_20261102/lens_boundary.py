"""Small independent GF(2)-similarity / H-information-loss certificate."""
import itertools,json

def graph(n,code):
 a=[0]*n
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  u,v=(i,j) if code>>k&1 else (j,i)
  a[u]|=1<<v
 return a

def product(a,b):
 ans=[]
 for row in a:
  v=0
  for j in range(len(a)):
   if row>>j&1:v^=b[j]
  ans.append(v)
 return ans

def ranks(a):
 n=len(a);p=[1<<i for i in range(n)];rs=[]
 for _ in range(n+1):
  basis={}
  for x in p:
   while x:
    i=x.bit_length()-1
    if i not in basis:basis[i]=x;break
    x^=basis[i]
  rs.append(len(basis));p=product(p,a)
 return rs

def h(a):
 n=len(a);d={1<<i:{i:1} for i in range(n)}
 for mask in range(1,1<<n):
  for v,c in d.get(mask,{}).copy().items():
   for w in range(n):
    if a[v]>>w&1 and not mask>>w&1:
     row=d.setdefault(mask|1<<w,{})
     row[w]=row.get(w,0)+c
 return sum(d[(1<<n)-1].values())

def main():
 a=graph(6,32767); wanted=[6,5,4,3,2,1,0]
 if ranks(a)!=wanted or h(a)!=1:raise RuntimeError('transitive control')
 for code in range(32768):
  b=graph(6,code)
  if ranks(b)==wanted and h(b)>1:
   literal=sum(all(b[u]>>v&1 for u,v in zip(p,p[1:])) for p in itertools.permutations(range(6)))
   if literal!=h(b):raise RuntimeError('literal mismatch')
   return {'status':'verified','interpretation':'Both adjacency matrices over GF(2) have one nilpotent Jordan block; H differs. This verifies inherited lens blindness, not a new classification theorem.', 'order':6,'rank_powers':wanted,'transitive':{'rows':a,'H':h(a)},'hostile':{'code':code,'rows':b,'H':literal},'encoding':'LSB lexicographic pairs, bit1 means i->j'}
 raise RuntimeError('no hostile')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
