import itertools as it,json
from pathlib import Path

def fixed(T):
 m=len(T)//2; d=[[0]*(2*m) for _ in range(1<<m)]
 for v in range(2*m):d[1<<(v//2)][v]=1
 for mask in range(1,1<<m):
  for v in range(2*m):
   for w in range(2*m):
    if not (mask>>(w//2)&1) and T[v]>>w&1:d[mask|1<<(w//2)][w]+=d[mask][v]
 return sum(d[-1][::2])
def total(T):
 n=len(T);d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if d[mask][v]:
    for w in range(n):
     if not(mask>>w&1) and T[v]>>w&1:d[mask|1<<w][w]+=d[mask][v]
 return sum(d[-1])
def cent(T):
 n=len(T);m=n//2;c=0
 for p in it.permutations(range(m)):
  pp=[2*p[v//2]+v%2 for v in range(2*m)]+([2*m] if n%2 else [])
  c+=all(((T[v]>>w)&1)==((T[pp[v]]>>pp[w])&1) for v in range(n) for w in range(n))
 return c
def main():
 import argparse
 parser=argparse.ArgumentParser()
 parser.add_argument('--witnesses',type=Path,default=Path(__file__).resolve().parent.parent/'session29-coherent-halfpaths'/'blue_checks.json')
 args=parser.parse_args()
 records=json.loads(args.witnesses.read_text())
 out=[]
 for r in records:
  T=r['outmasks']; f=fixed(T);c=cent(T)
  assert f==r['fixed_hp'] and c==r['centralizer']
  row={'n':len(T),'f':f,'centralizer':c,'blue':f//c}
  if len(T)<=11:row['total_H']=total(T)
  out.append(row)
 Path(__file__).with_name('literal_results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(out)

if __name__=='__main__':main()
