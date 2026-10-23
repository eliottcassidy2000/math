"""Independent triangle enumeration verifies aligned ternary decoding."""
from pathlib import Path
from itertools import combinations
import argparse,csv,json,hashlib
HERE=Path(__file__).resolve().parent
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise AssertionError(msg)
def matrix(code,n):
 A=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):A[i][j]=code>>(n*(n-1)//2-1-k)&1;A[j][i]=1-A[i][j]
 return A
def check(A):
 n=len(A);c=[0]*n;e=[[0]*n for _ in range(n)];c3=0
 for i,j,k in combinations(range(n),3):
  if A[i][j]==A[j][k]==A[k][i]:
   c3+=1
   for v in (i,j,k):c[v]+=1
   for u,v in ((i,j),(i,k),(j,k)):e[u][v]+=1;e[v][u]+=1
 need(sum(c)==3*c3,'triangle total');r=[0]*n;rs=[0]*n;d=[0]*n;ds=[0]*n
 for u in range(n):
  for v in range(n):
   if A[u][v]:
    r[u]+=(c[v]+e[u][v])%3;rs[v]+=(c[u]+e[u][v])%3
    d[u]+=(c[v]-e[u][v])%3;ds[v]+=(c[u]-e[u][v])%3
    need(((c[v]+e[u][v])-(c[u]+e[u][v]))%3==(c[v]-c[u])%3,'gradient')
 rec=[(r[u]+rs[u])%3 for u in range(n)]
 for u in range(n):
  need(rec[u]==c[u]%3,'aligned c recovery');need(r[u]%3==(c[u]+sum(A[u][v]*c[v] for v in range(n)))%3,'forward sum operator');need(rs[u]%3==(c[u]+sum(A[v][u]*c[v] for v in range(n)))%3,'converse sum operator');need((d[u]+ds[u])%3==0,'difference decoder collapse')
  for v in range(n):
   if A[u][v]:need(((c[v]+e[u][v])%3-rec[v])%3==e[u][v]%3,'aligned edge recovery')
 return c,r,rs,d,ds

def main():
 p=argparse.ArgumentParser();p.add_argument('--classes',type=Path,default=HERE.parent/'wiring/classes_input_n8.tsv');p.add_argument('--write',action='store_true');args=p.parse_args()
 if not args.classes.exists() and args.classes==HERE.parent/'wiring/classes_input_n8.tsv':args.classes=HERE.parent/'session38-wiring/classes_input_n8.tsv'
 levels=[]
 for n in range(2,7):
  N=1<<(n*(n-1)//2)
  for code in range(N):check(matrix(code,n))
  levels.append({'n':n,'all_labeled_tournaments':N})
 R=list(csv.DictReader(args.classes.open(),delimiter='\t'));need(len(R)==6880,'class count')
 for r in R:check(matrix(int(r['canon']),8))
 A=[[0,1,0],[0,0,1],[1,0,0]];c,r,rs,d,ds=check(A);need(c==[1,1,1] and all(x%3==0 for x in d+ds),'minus collapse hostile')
 O={'scope':'PROVED all-size identities; exhaustive all labeled tournaments n2..6 and all6880 order8 class reps','classes_sha256':hashlib.sha256(args.classes.read_bytes()).hexdigest(),'levels':levels,'order8_classes':6880,'cyclic_triangle_hostile':{'outgoing_rows':[2,4,1],'vertex_triangle_counts':c,'plus_weighted_row_sums_mod3':[x%3 for x in r],'plus_converse_weighted_row_sums_mod3':[x%3 for x in rs],'minus_weighted_row_sums_mod3':[x%3 for x in d],'minus_converse_weighted_row_sums_mod3':[x%3 for x in ds]},'checks':checks};dest=HERE/'three_color_receipt.json'
 if args.write:dest.write_text(json.dumps(O,indent=2)+'\n')
 else:need(json.loads(dest.read_text())==O,'frozen three-color receipt differs')
 print(json.dumps(O,indent=2))
if __name__=='__main__':main()
