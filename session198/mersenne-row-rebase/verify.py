"""Standalone exact source replay and literal tournament subset checks."""
import json,math,hashlib

def require(x,m):
 if not x:raise RuntimeError(m)
def move(v,c):
 D,E,U,d=v
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 if c=='J':return D-E+U,D-E,D+U,D
 raise RuntimeError('letter')
def run(word):
 v=(1,0,2,1)
 for c in word:v=move(v,c)
 return v
def B(v):return v[0]+2*v[2]
def words(b):
 require(b>=5,'exponent domain')
 if b==5:return 'OZJZO'
 if b%2==0:
  i=(b-4)//2;return 'O'*i+'JOO'+'Z'*i+'J'
 i=(b-5)//2;return 'O'*i+'JOZO'+'Z'*i+'J'
def graph(word):
 rows=[2,0];root=0
 for c in word:
  n=len(rows);rows.append(0)
  for v in range(n):
   if c!='O' and v==root:rows[n]|=1<<v
   else:rows[v]|=1<<n
  if c=='J':root=n
 return rows,root

def literal(rows,root):
 n=len(rows);size=1<<n
 require(all(not rows[i]>>i&1 for i in range(n)),'loops')
 require(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n)for j in range(i+1,n)),'tournament')
 end=[[0]*n for _ in range(size)];start=[[0]*n for _ in range(size)]
 for v in range(n):end[1<<v][v]=1
 start[1<<root][root]=1
 H=[0]*size;H[0]=1
 for mask in range(1,size):
  for v in range(n):
   if not mask>>v&1:continue
   old=mask^(1<<v)
   if old:
    end[mask][v]=sum(end[old][u]for u in range(n)if old>>u&1 and rows[u]>>v&1)
    start[mask][v]=sum(start[old][u]for u in range(n)if old>>u&1 and rows[u]>>v&1)
  H[mask]=sum(end[mask])
 full=size-1
 profile=(H[full],end[full][root],sum(H[m]*sum(start[full^m])for m in range(size)if not m>>root&1),H[full^(1<<root)])
 cycles=sum(start[full][v]for v in range(n)if rows[v]>>root&1)
 return profile,sum(start[full]),cycles

def main():
 closed=0
 for i in range(1,65):
  x=1<<i
  for kind in('even','odd'):
   b=2*i+(4 if kind=='even'else 5);q=run(words(b))
   want=(16*x*x-1,8*x*x-3*x-2,16*x*x,8*x*x-3*x-1)if kind=='even'else(32*x*x-1,16*x*x-5*x-2,32*x*x,16*x*x-5*x-1)
   require(q==want,'source profile')
   require(B(q)==3*(1<<b)-1 and math.gcd(q[0],q[2])==1,'target and primitive')
   out=move(move(q,'J'),'Z');require(out==(B(q),q[0],B(q),q[0]),'balanced map')
   require(len(words(b))+4==b+4,'order');closed+=1
 polynomial=0
 for i in range(1,11):
  for j in range(1,11):
   x,y=1<<i,1<<j
   require(B(run('O'*i+'JOO'+'Z'*j+'J'))==48*x*y-18*x+18*y-1,'even polynomial')
   require(B(run('O'*i+'JOZO'+'Z'*j+'J'))==96*x*y-30*x+30*y-1,'odd polynomial')
   polynomial+=2
 require(run(words(5))==(19,0,38,5),'seed profile')
 records=[]
 for b in range(5,11):
  word=words(b);q=run(word);rows,root=graph(word);prof,s,c=literal(rows,root)
  require(prof==q and s==1,'literal source')
  rows,root=graph(word+'JZ');prof,s,c=literal(rows,root)
  target=3*(1<<b)-1
  require(prof==(target,q[0],target,q[0])and s==1 and c==1,'literal balanced graph')
  require(len(rows)==b+4,'literal order')
  records.append(dict(b=b,H=target,C=c,order=len(rows),source_word=word,source_profile=q,balanced_root=root,balanced_rows=rows))
 print(json.dumps(dict(status='PASS',closed_profile_checks=closed,two_run_polynomial_checks=polynomial,literal_balanced_graphs=records),indent=2,sort_keys=True))
if __name__=='__main__':main()
