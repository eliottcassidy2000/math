#!/usr/bin/env python3
from verify import *

def induced(r,A):
 keep=list(bits(A));return [sum(1<<b for b,j in enumerate(keep)if r[i]>>j&1)for i in keep]
def main():
 parents=paths=cycles=strongparents=0
 for n in range(2,6):
  full=(1<<n)-1
  for code in range(1<<(n*(n-1)//2)):
   r=rows_from_code(n,code);H,C,h=boundary(r);hs=[hp(induced(r,A))for A in range(1<<n)]
   P=[[0]*n for _ in r];q=[[0]*n for _ in r]
   for A in range(1,1<<n):
    if A.bit_count()%2:continue
    vs=list(bits(A))
    for p in permutations(vs):
     if all(r[a]>>b&1 for a,b in zip(p,p[1:])):
      P[p[0]][p[-1]]+=hs[full^A];paths+=1
      if p[0]==min(vs) and r[p[-1]]>>p[0]&1:
       for a,b in zip(p,p[1:]+p[:1]):q[a][b]+=hs[full^A]
       cycles+=1
   for i,j in combinations(range(n),2):need(C[i][j]==2*(P[i][j]+P[j][i]),'capacity path expansion')
   for i in range(n):
    for j in bits(r[i]):need(q[i][j]==P[j][i],'even cycle circulation expansion')
    need(sum(q[i])==sum(q[j][i]for j in range(n)),'zero divergence')
   if strong(r):
    vals=[flow(r,C,0,v)[0]for v in range(1,n)]+[flow(r,C,v,0)[0]for v in range(1,n)]
    need(min(vals)==min(cut(r,C,S)for S in range(1,full)),'2(n-1) mincut cover');strongparents+=1
   parents+=1
 out={'status':'PASS','parents':parents,'even_paths':paths,'even_cycles':cycles,'strong_parents_checked':strongparents}
 raw=json.dumps(out,sort_keys=True)+'\n';Path(__file__).with_name('paths_receipt.json').write_text(raw);print(raw,end='')
if __name__=='__main__':main()
