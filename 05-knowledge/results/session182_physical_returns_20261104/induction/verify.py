import itertools,json
from pathlib import Path
OUT=Path(__file__).parent

def require(x,msg):
 if not x: raise RuntimeError(msg)

def tournaments(n):
 for mask in range(1<<(n*(n-1)//2)):
  rows=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if mask>>k&1: rows[i]|=1<<j
    else: rows[j]|=1<<i
    k+=1
  yield rows

def sig(rows,W):
 n=len(rows);h=p=f=g=j=e=s=0
 for perm in itertools.permutations(range(n)):
  bad=[k for k in range(n-1) if not(rows[perm[k]]>>perm[k+1]&1)]
  if not bad:
   h+=1;e+=bool(W>>perm[-1]&1);s+=bool(W>>perm[0]&1)
  for cut in range(n+1):
   if any(k!=cut-1 for k in bad):continue
   p+=1
   if cut and W>>perm[cut-1]&1:f+=1
   if cut<n and W>>perm[cut]&1:g+=1
   if 0<cut<n and W>>perm[cut-1]&1 and W>>perm[cut]&1:j+=1
 return (h,p,f,g,j,e,s)

def combine(x,y):
 h,p,f,g,j,e,s=x;H,P,F,G,J,E,S=y
 return(h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H)

def ordinal(a,b):
 n=len(a);m=len(b)
 return [r+(((1<<m)-1)<<n) for r in a]+[r<<n for r in b]

def strong(rows):
 n=len(rows)
 for v in range(n):
  seen=1<<v
  while True:
   nxt=seen
   for u in range(n):
    if seen>>u&1:nxt|=rows[u]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True

def terminal_met(rows,W):
 n=len(rows)
 for v in range(n):
  seen=1<<v
  while True:
   nxt=seen
   for u in range(n):
    if seen>>u&1:nxt|=rows[u]
   if nxt==seen:break
   seen=nxt
  if not seen&W:return False
 return True

def attach(rows,W,a=1):
 n=len(rows);r=a+n
 out=[]
 for i in range(a):out.append(sum(1<<k for k in range(i+1,a+n)))
 for i,z in enumerate(rows):out.append((z<<a)+((1<<r) if W>>i&1 else 0))
 out.append(((1<<a)-1)+sum(1<<(a+i) for i in range(n) if not(W>>i&1)))
 return out

def main():
 checks=0;sc=0;ev=0
 banks={n:[(r,w,sig(r,w)) for r in tournaments(n) for w in range(1<<n)] for n in range(1,5)}
 for n in range(1,4):
  for m in range(1,4):
   if n+m>5:continue
   for a,u,x in banks[n]:
    for b,v,y in banks[m]:
     require(combine(x,y)==sig(ordinal(a,b),u+(v<<n)), 'ordinal transfer');checks+=1
 for n,bank in banks.items():
  for b,w,x in bank:
   t=attach(b,w)
   require(strong(t)==terminal_met(b,w),'strongness');sc+=1
   require(sig(t,1<<(n+1))[0]==x[0]-x[4]+2*x[2],'multiport H');ev+=1
   z=sig(t,1<<(n+1));require(z[2]-z[0]==x[4],'defect')
 binary_checks=0
 for m in range(1,6):
  rows=[sum(1<<k for k in range(i+1,m)) for i in range(m)]
  for W in range(1<<m):
   f=g=j=e=0
   for i in range(m):
    b=W>>i&1
    f,g,j,e=f+b*(1<<i),2*g+b,j+b*(g+e),b
   z=sig(rows,W)
   require((f,g,j,e)==(z[2],z[3],z[4],z[5]), "binary transfer")
   binary_checks+=1
 require(sig([6,4,0],7)==(1,8,7,7,6,1,1),"source signature")
 require(sig([30,28,24,16,0],21)==(1,32,21,21,12,1,1),"target signature")
 examples=[]
 for sclock in range(1,5):
  a=[sum(1<<j for j in range(i+1,3)) for i in range(3)]
  b=[sum(1<<j for j in range(i+1,5)) for i in range(5)]
  # DP instead of factorial enumeration for full large examples.
  def hp(rows):
   n=len(rows);dp=[[0]*n for _ in range(1<<n)]
   for i in range(n):dp[1<<i][i]=1
   for mask in range(1,1<<n):
    for i,z in enumerate(dp[mask]):
     if z:
      for k in range(n):
       if not mask>>k&1 and rows[i]>>k&1:dp[mask|1<<k][k]+=z
   return sum(dp[-1])
  aa=attach(a,7,sclock);bb=attach(b,21,sclock)
  ha=hp(aa);hb=hp(bb)
  require(hb==3*ha+4 and strong(aa) and strong(bb),'uniform +4')
  examples.append({'clock':sclock,'source_H':ha,'target_H':hb,'source_rows':aa,'target_rows':bb})
 result={'binary_signature_checks':binary_checks,'ordinal_signature_checks':checks,'strongness_checks':sc,'attachment_and_defect_checks':ev,'uniform_plus4_examples':examples}
 (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
