import itertools,json
from verify import SEEDS,inflate,reverse

def hp(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m in range(1,1<<n):
  for i in range(n):
   if not d[m][i]:continue
   opts=a[i]&~m
   while opts:
    b=opts&-opts;opts-=b;d[m|b][b.bit_length()-1]+=d[m][i]
 return sum(d[-1])
def extend(t,a,c):
 n=len(t);out=list(t)+[1<<(n+1),0]
 for i in range(n):
  for p,m in [(n,a),(n+1,c)]:
   if m>>i&1:out[i]|=1<<p
   else:out[p]|=1<<i
 return out
if __name__=='__main__':
 records=[]
 for r in range(1,5):
  for j,s in enumerate(SEEDS):
   t=inflate(s,r);rev=reverse(t);n=len(t)
   for a,c in [(0,0)]+[(1<<i,0) for i in range(n)]+[(0,1<<i) for i in range(n)]:
    f=hp(extend(t,a,c));g=hp(extend(rev,a,c))
    if f!=g:
     records.append({'r':r,'state':j,'a':a,'c':c,'F':f,'F_op':g});break
   else:raise RuntimeError('needs higher query')
 print(json.dumps({'literal_two_port_witnesses':records},indent=2))
