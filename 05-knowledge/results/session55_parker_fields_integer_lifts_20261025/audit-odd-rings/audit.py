import json,itertools,math,hashlib
from pathlib import Path
D=Path(__file__).resolve().parent
LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
PAIRS=list(itertools.combinations(range(9),2))
BASE=[5,7,11,13,17,19,23,25,31,43,47,49,67]
ODD=sorted([1,3,9]+BASE+[3*x for x in BASE])
def chk(b,msg):
 if not b: raise RuntimeError(msg)
def arrays(n):
 Q={x*x%n for x in range(n)}
 # Independent construction from top row a,b,c and center e, using b=3e-a-c.
 for e in Q:
  for a in Q:
   i=(2*e-a)%n
   if i not in Q:continue
   for c in Q:
    g=(2*e-c)%n
    if g not in Q:continue
    b=(3*e-a-c)%n
    if b not in Q:continue
    h=(2*e-b)%n
    if h not in Q:continue
    d=(3*e-a-g)%n
    if d not in Q:continue
    f=(2*e-d)%n
    if f in Q:
     yield(a,b,c,d,e,f,g,h,i)
def masks(n):
 out={};count=0
 for m in arrays(n):
  chk(len({sum(m[i] for i in line)%n for line in LINES})==1,'magic')
  count+=1
  mask=sum(1<<k for k,(i,j) in enumerate(PAIRS) if m[i]==m[j])
  out.setdefault(mask,m)
 small={}
 for x in sorted(out,key=int.bit_count):
  if not any(x&y==y for y in small):small[x]=out[x]
 return count,out,small
def crt(a,b,m,n):return (a+m*((b-a)*pow(m,-1,n)%n))%(m*n)
def verify(m,n):
 Q={x*x%n for x in range(n)}
 chk(len(set(m))==9 and all(x in Q for x in m),'witness square/distinct')
 chk(len({sum(m[i] for i in line)%n for line in LINES})==1,'witness lines')
def main():
 tabs={};results={};P=[];global_masks={}
 for n in sorted(set(ODD+[2,4,8,16,32,64])):
  ct,full,mini=masks(n);tabs[n]=mini
  chk(0 not in full,'claimed Parker factor '+str(n))
  results[n]={'arrays':ct,'all_masks':len(full),'minimal_masks':{str(k):v for k,v in mini.items()}}
 for odd in ODD:
  for k in range(7):
   even=2**k;n=even*odd
   witness=None
   intersections={a&b for a in tabs[even] for b in tabs[odd]}
   minimum=[]
   for x in sorted(intersections,key=int.bit_count):
    if not any(x&y==y for y in minimum):minimum.append(x)
   global_masks[n]=sorted(minimum)
   for a,ma in tabs[even].items():
    for b,mb in tabs[odd].items():
     if a&b==0:
      witness=[crt(x,y,even,odd) for x,y in zip(ma,mb)]
      verify(witness,n);break
    if witness:break
   if witness is None:P.append(n)
 out={'odd_Parker_including_1':ODD,'all_Parker_candidates_that_pass':sorted(P),'largest':max(P),'local_tables':results}
 producer=D.parent/'crt-wiring/classification_candidates.json'
 if producer.exists():
  reference=json.loads(producer.read_text())
  chk(sorted(P)==reference['parker_moduli'],'producer Parker list differs')
  for n,t in reference['local_antichains'].items():
   chk(set(tabs[int(n)])=={int(row['mask'],16) for row in t},'local antichain differs '+n)
  for row in reference['rows']:
   chk(global_masks[row['n']]==sorted(int(x,16) for x in row['minimal_masks']),'global antichain differs '+str(row['n']))
  out['producer_comparison']='All local and global minimal-mask sets and all203 statuses agree.'
 (D/'audit_results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'odd':ODD,'all_Parker':sorted(P),'count':len(P),'largest':max(P)},indent=2))
if __name__=='__main__':main()
