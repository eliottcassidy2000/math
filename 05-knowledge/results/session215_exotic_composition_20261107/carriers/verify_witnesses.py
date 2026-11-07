from inherited import PAIRS,profile,remark,move,graph_move,base,need
import json

def hp(rows):
 # Forward subset DP, independent of inherited ends/starts profile engine.
 n=len(rows);full=(1<<n)-1;dp=[[0]*n for _ in range(full+1)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,full+1):
  for u,a in enumerate(dp[mask]):
   if not a:continue
   avail=rows[u]&~mask
   while avail:
    bit=avail&-avail;avail-=bit;v=bit.bit_length()-1
    dp[mask|bit][v]+=a
 return sum(dp[full])

def compiler(rows):
 n=len(rows);out=list(rows)+[0]*3;core=(2,8,3,5);labels=[0,n,n+1,n+2]
 for i in range(4):
  for j in range(4):
   if core[i]>>j&1:out[labels[i]]|=1<<labels[j]
 for v in labels[1:]:
  for u in range(1,n):out[v]|=1<<u
 return tuple(out)

if __name__ == '__main__':
 witnesses=[(1547,'A','JJ'),(759,'B','ZOOZ'),(677,'A','ZZJJJ'),(677,'B','ZOOJJ')]
 result=[]
 for num,side,w in witnesses:
  pair=next(z for z in PAIRS if z[0]==num);_,js,root,a,b=pair
  j=js[side=='B'];rows=remark((a,b)[side=='B'],root);p=profile(rows)
  for c in w:rows=graph_move(rows,c);p=move(p,c)
  need(profile(rows)==p,'literal complete nine-state mismatch')
  x=base(p,j);r=1+(side=='B');need(p[1]==r+3*j,'start label')
  child=compiler(rows);count=hp(child);need(count==3*x+2*r,'compiler count')
  result.append(dict(seed=num,side=side,word=w,j=j,r=r,rows=rows,profile=p,X=x,compiled_rows=child,compiled_count=count))
 print(json.dumps(result,indent=2))
