"""Independent exact DP controls for the three-state twin contrast bridge."""
import itertools,json
SEEDS=((16,5,33,39,46,3),(16,5,49,39,42,3),(18,4,49,39,42,3))
def need(ok,msg):
 if not ok: raise RuntimeError(msg)
def inflate(a,r):
 out=[0]*(5+r)
 for u in range(5):
  out[u]=a[u]&31
  for v in range(5,5+r):
   if a[u]>>5&1:out[u]|=1<<v
   else:out[v]|=1<<u
 for u,v in itertools.combinations(range(5,5+r),2):out[u]|=1<<v
 return out
def reverse(a,subset=None):
 subset=range(len(a)) if subset is None else subset
 b=list(a)
 for u,v in itertools.combinations(subset,2):b[u]^=1<<v;b[v]^=1<<u
 return b
def packet(a):
 n=len(a);full=(1<<n)-1
 ends=[[0]*n for _ in range(1<<n)];starts=[[0]*n for _ in range(1<<n)]
 for i in range(n):ends[1<<i][i]=starts[1<<i][i]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if not(mask>>v&1):continue
   rest=mask^(1<<v)
   if not rest:continue
   ends[mask][v]=sum(ends[rest][u] for u in range(n) if rest>>u&1 and a[u]>>v&1)
   starts[mask][v]=sum(starts[rest][u] for u in range(n) if rest>>u&1 and a[v]>>u&1)
 c=[[0]*n for _ in range(n)]
 for mask in range(1,full):
  other=full^mask
  for u in range(n):
   for v in range(n):c[u][v]+=ends[mask][u]*starts[other][v]
 return sum(ends[full]),c

def contrast(c,a=5,b=6,w=7):return c[b][w]-c[a][w]+c[w][a]-c[w][b]
if __name__=='__main__':
 records=[]
 for r in range(3,8):
  for j,seed in enumerate(SEEDS):
   t=inflate(seed,r);h,c=packet(t);expected=2*(4*3**(r-2)+9*2**(r-2)-3+2*j)
   need(contrast(c)==expected,('contrast',r,j))
   for name,u in [('whole_converse',reverse(t)),('module_converse',reverse(t,range(5,5+r)))]:
    hh,cc=packet(u);need(hh==h,('H',name,r,j));need(contrast(cc)==-expected,('reverse contrast',name,r,j))
   records.append({'r':r,'state':j,'H':h,'L':expected,'converse_L':-expected})
 print(json.dumps({'status':'PASS','controls':records},indent=2,sort_keys=True))
