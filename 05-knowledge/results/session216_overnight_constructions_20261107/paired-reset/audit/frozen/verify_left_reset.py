from inherited import profile,need,graph_move
from compiler_control import compiler,hp
import json,itertools

def Q(rows):
 n=len(rows);return tuple([rows[0]|(1<<n)]+list(rows[1:])+[(1<<n)-2])
def P(rows):
 n=len(rows);return tuple(list(rows)+[(1<<n)-1])
def seven(p):return tuple(p[i] for i in [0,1,2,3,4,5,7])
def mv(p,c):
 h,s,e,L,R,d,b=p
 if c=='O':return(h,s,0,L,2*R,d,2*b)
 if c=='P':return(h,0,e,2*L,R,d,2*b)
 if c=='Z':return(h+R-e,s,d,L+b,2*R-d,d,2*b)
 return(h+L-s,d,e,2*L-d,R+b,d,2*b)
def reduced(p):
 h,s,e,L,R,d,b=p;need(s==d,'reset condition');r=d%3;j=(d-r)//3
 return(h+L-d+2*j,L-d,b,R+b-e,d-e)
def redmove(v,c):
 X,a,b,w,t=v
 return (X+2*a,2*a,2*b,w+2*b,t) if c=='Q' else (X+w,a+b,2*b,2*w-2*t,0)
controls=0
for n in range(2,5):
 for bits in range(1<<(n*(n-1)//2)):
  g=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if bits>>k&1:g[i]|=1<<j
    else:g[j]|=1<<i
    k+=1
  p=seven(profile(g))
  for c in 'OZPQ':
   child=Q(g) if c=='Q' else P(g) if c=='P' else graph_move(g,c)
   need(seven(profile(child))==mv(p,c),'literal transfer');controls+=1
  q=mv(p,'Q');v=reduced(q)
  for c in 'QZ':need(reduced(mv(q,c))==redmove(v,c),'reduced transfer');controls+=1
x=list(map(int,open('../session216-carriers/BANK.txt').read().split()));i=1;bank=[]
for _ in range(x[0]):n=x[i];bank.append(x[i+1:i+1+n]);i+=n+1
pair=json.load(open('RESET_PAIRED_RAYS.json'))[0];records=[]
for side in ['A','B']:
 z=pair[side];g=bank[z['ix']];mask=z['mask'];rows=tuple([mask<<1]+[(g[j]<<1)|(0 if mask>>j&1 else 1) for j in range(len(g))]);need(tuple(z['profile'])==profile(rows),'initial source')
 rows=Q(rows);state=seven(profile(rows));need(reduced(state)[:2]==(69,36),'paired reset')
 rec={'side':side,'initial_rows':list(rows),'profile':profile(rows),'j':z['j'],'counts':[]}
 for k in range(5):
  p=seven(profile(rows));v=reduced(p);need(v[0]==72*(1<<k)-3 and v[1]==36*(1<<k),'ray')
  child=compiler(rows);h=hp(child);r=p[1]%3;need(h==3*v[0]+2*r,'literal compiler');rec['counts'].append({'k':k,'X':v[0],'H':h,'order':len(rows)})
  rows=Q(rows)
 records.append(rec)
print(json.dumps({'status':'PASS','transfer_controls':controls,'paired_ray':records},indent=2))
