"""Independent complete cut-polynomial suffix DP, with no producer imports."""
from pathlib import Path
from functools import lru_cache
import hashlib,json,math
D=Path(__file__).resolve().parent
B=64; MASK=(1<<B)-1

def need(c,s):
 if not c:raise RuntimeError(s)

def tournament(n,code):
 rows=[0]*n;p=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>p&1:rows[i]|=1<<j
   else:rows[j]|=1<<i
   p+=1
 return rows

def kernels(R,i,k):
 n=len(R); other=((1<<n)-1)^(1<<i); N=n+k-1
 need(math.factorial(N)*(1<<(N-1)) < 1<<B,'packed polynomial no-carry bound')
 @lru_cache(None)
 def suffix(S,r,last):
  if S==0 and r==0:return 1
  out=0
  if r:
   child=suffix(S,r-1,i);out+=child<<B
   if last!=i and R[last]>>i&1:out+=child
  left=S
  while left:
   b=left&-left;left-=b;v=b.bit_length()-1;child=suffix(S^b,r,v);out+=child<<B
   if R[last]>>v&1:out+=child
  return out
 poly=suffix(other,k-1,i)
 left=other
 while left:
  b=left&-left;left-=b;poly+=suffix(other^b,k,b.bit_length()-1)
 return [(poly>>(B*d))&MASK for d in range(N)]

def main():
 data=[];coords=0
 for n in range(2,6):
  for code in range(1<<(n*(n-1)//2)):
   R=tournament(n,code)
   for i in range(n):
    for k in [2,3,5]:
     K=kernels(R,i,k);want={2:3,3:3,5:5}[k]
     for h,value in enumerate(K,1):
      if h>=2:need(value%2==int(h==want),'universal slot parity '+str((R,i,k,h,value)))
      elif k in [3,5]:need(value%2==0,'return-segment even Hamilton parity')
     need(K[-1]==math.factorial(len(K))//math.factorial(k),'all-singleton ordered slot count')
     data.append({'rows':R,'distinguished':i,'multiplicity':k,'K_1_through_length':K});coords+=len(K)
 out=json.dumps(data,separators=(',',':'))+'\n';p=D/'kernel_recount.json'
 if p.exists():need(p.read_text()==out,'retained complete kernel bank stable')
 else:p.write_text(out)
 return {'status':'PASS','method':'Independent suffix recurrence with one packed exact coefficient per cut count; no producer import.','complete_kernel_records':len(data),'exact_cut_coordinates':coords,'maximum_word_length':9,'labeled_tournaments_with_nonempty_complement':sum(1<<(n*(n-1)//2) for n in range(2,6)),'h_beyond_word_length':'zero by nonempty path slots','source_sha256':hashlib.sha256((D/'kernel_audit.py').read_bytes()).hexdigest(),'data_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
