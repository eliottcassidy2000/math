import itertools as it,collections as co,json
from pathlib import Path

def require(x,m):
 if not x:raise RuntimeError(m)
def edges(n):return list(it.combinations(range(n),2))
def arc(m,n,i,j):return (m>>edges(n).index((i,j)))&1 if i<j else 1-arc(m,n,j,i)
def induce(m,n,p):return sum(arc(m,n,p[i],p[j])<<b for b,(i,j) in enumerate(edges(len(p))))
def complement(m,n):
 keep=sum(1<<b for b,(i,j) in enumerate(edges(n)) if j==i+1)
 return m^(((1<<len(edges(n)))-1)^keep)
def hdata(m,n):
 out=[sum(arc(m,n,i,j)<<j for j in range(n) if i!=j) for i in range(n)]
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for v in range(n):
   x=dp[s][v]
   if not x:continue
   todo=out[v]&~s
   while todo:
    bit=todo&-todo;todo-=bit;w=bit.bit_length()-1;dp[s|bit][w]+=x
 return [sum(row) for row in dp],dp

def tilelist(m):
 return sorted({induce(m,5,p) for p in it.permutations(range(5)) if all(arc(m,5,p[i],p[i+1]) for i in range(4))})
def ports(m):
 h,e=hdata(m,5);_,s=hdata(m^1023,5)
 qe=[sum(e[sub][v]*h[31^sub] for sub in range(1,31) if(sub>>v)&1) for v in range(5)]
 qs=[sum(s[sub][v]*h[31^sub] for sub in range(1,31) if(sub>>v)&1) for v in range(5)]
 return {'H':h[31],'E':e[31],'S':s[31],'QE':qe,'QS':qs}
def join(ts):
 n=5*len(ts)
 return sum((arc(ts[i//5],5,i%5,j%5) if i//5==j//5 else 1)<<b for b,(i,j) in enumerate(edges(n)))
def mirror(m,n):return induce(m^((1<<len(edges(n)))-1),n,list(range(n-1,-1,-1)))
T={m:tilelist(m) for m in(40,41)}
require([len(T[m]) for m in(40,41)]==[5,13],'seed tile counts')
P={t:ports(complement(t,5)) for m in T for t in T[m]}
rows=[];seen=set();hist=co.Counter();bluehist=co.Counter()
for letter1,letter2 in it.product((40,41),repeat=2):
 for t,u in it.product(T[letter1],T[letter2]):
  orig=join((t,u));target=complement(orig,10)
  require(target not in seen,'distinct target tile');seen.add(target)
  require(arc(target,10,9,0)==1,'Hamiltonian cycle closure')
  hx=P[t];hy=P[u];a=4;b=0
  pred=hx['H']*hy['H']-hy['E'][b]*hx['S'][a]+(hx['E'][a]+hx['QE'][a])*(hy['S'][b]+hy['QS'][b])
  actual=hdata(target,10)[0][-1];require(pred==actual,'two block seam formula')
  blue=mirror(orig,10)==orig
  require(blue==(mirror(target,10)==target),'complement preserves mirror')
  hist[actual]+=1
  if blue:bluehist[actual]+=1
  rows.append({'word':[letter1,letter2],'tiles':[t,u],'target':target,'H':actual,'blue':blue})
require(len(rows)==324 and sum(bluehist.values())==18,'outward line budgets')
result={'tile_alphabet':T,'port_data':P,'edges':324,'blue_edges':18,'black_edges':306,'H_histogram':sorted(hist.items()),'blue_H_histogram':sorted(bluehist.items()),'sum_H':sum(k*v for k,v in hist.items()),'sum_blue_H':sum(k*v for k,v in bluehist.items()),'rows':rows}
Path(__file__).with_name('seamed-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in('rows','port_data')},indent=2))
