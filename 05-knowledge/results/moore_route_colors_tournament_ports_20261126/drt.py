from itertools import combinations
from collections import Counter
# McKay's two order-15 DR tournaments (upper triangle, row-major).
ROWS=(
"111111100000001110001111000100111010011101101011101010011101000100111010010111110100010110001011001010101",
"111111100000001110001111000100111010011110100111100110101101001000111001010111110100100110101011001010101",
)
def decode(word,n):
 a=[0]*n
 for (i,j),b in zip(combinations(range(n),2),word,strict=True):
  if b=='0':i,j=j,i
  a[i]|=1<<j
 return a
def paths(a):
 n=len(a); U=(1<<n)-1
 D=[[0]*n for _ in range(U+1)]
 for v in range(n):D[1<<v][v]=1
 for s in range(1,U):
  c=s
  while c:
   t=c&-c;c^=t;j=t.bit_length()-1;w=D[s][j]
   x=a[j]&~s
   while x:
    t=x&-x;x^=t;i=t.bit_length()-1
    D[s|t][i]+=w
 return sum(D[-1]),Counter(D[-1]),D
def main():
 for w in ROWS:
  a=decode(w,15)
  assert all(x.bit_count()==7 for x in a)
  assert all((a[i]&a[j]).bit_count()==3 for i,j in combinations(range(15),2))
  for i in range(15):
   for j in range(15):
    assert sum((a[i]>>h&1)*(a[h]>>j&1) for h in range(15))+(a[i]>>j&1)==4*(i!=j)
  h,ends,D=paths(a)
  print("H",h,"per-end histogram",sorted(ends.items()),
    "vertex deletion H histogram", sorted(Counter(sum(D[-1^(1<<j)]) for j in range(15)).items()))

def edit_catalog(a,D):
 n=len(a); U=(1<<n)-1
 ar=[sum((a[j]>>i&1)<<j for j in range(n)) for i in range(n)]
 _,_,R=paths(ar)
 E=[[0]*n for _ in range(n)]
 for s in range(1,U):
  t=U^s;j=s
  while j:
   b=j&-j;j^=b;u=b.bit_length()-1;x=D[s][u]
   if not x:continue
   k=t
   while k:
    b=k&-k;k^=b;v=b.bit_length()-1;E[u][v]+=x*R[t][v]
 res=Counter(E[v][u]-E[u][v] for u in range(n) for v in range(n) if a[u]>>v&1)
 return res

if __name__=='__main__':
 main()
 arc_spectra=[]
 for idx,w in enumerate(ROWS):
  a=decode(w,15);h,_,D=paths(a);r=edit_catalog(a,D)
  arc_spectra.append(r)
  print("arc effects",idx,
    len(r),"distinct", "min/max",min(r),max(r),"near zero",sorted(r.items(),key=lambda p:abs(p[0]))[:3])
  for u,v in [(0,1),(3,8)]:
   if not(a[u]>>v&1):u,v=v,u
   a[u]^=1<<v;a[v]^=1<<u
   h2,_,_=paths(a)
   assert h2-h in r
   a[u]^=1<<v;a[v]^=1<<u
 assert arc_spectra[0] == arc_spectra[1]
