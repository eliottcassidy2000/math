from pathlib import Path
import json,itertools
P=Path(__file__).parent
# Extract independent definitions without running that file's enumeration driver.
s=(P/'independent.py').read_text();ns={'__file__':str(P/'independent.py')};exec(s[:s.index('hist={};num=0')],ns);exec(s[s.index('def paths(a):'):s.index('small=[a')],ns);build=ns['build']
def table(a):
 n=len(a);N=1<<n;E=[[0]*n for S in range(N)];H=[0]*N;H[0]=1
 for j in range(n):E[1<<j][j]=1
 for S in range(1,N):
  H[S]=sum(E[S])
  for j in range(n):
   if not E[S][j]:continue
   z=a[j]&~S
   while z:
    b=z&-z;z-=b;E[S|b][b.bit_length()-1]+=E[S][j]
 return H,E

def state(a,cut):
 H,E=table(a);n=len(a);rev=[sum((a[j]>>i&1)<<j for j in range(n)) for i in range(n)];_,S=table(rev);F=len(H)-1
 starts=[sum(row[j] for j in range(n) if cut>>j&1) for row in S];ends=[sum(row[j] for j in range(n) if cut>>j&1) for row in E]
 return H[F],sum(H[F^z]*starts[z] for z in range(F+1)),sum(ends[z]*starts[F^z] for z in range(F+1)),ends[F]
seeds=json.loads((P/'frozen/MENU.json').read_text())['seeds'];graphs=[];states=0
for i,rec in enumerate(seeds):
 a=rec['rows'];cut=rec['cut'];expected=(9,245,150+2*i,9)
 if state(a,cut)!=expected:raise RuntimeError('initial state')
 core=build(a,[],cut,0,[0]);l=len(a);mapping={j:j+2 for j in range(l)};mapping[l]=1;mapping[l+1]=0;cy=[mapping[j] for j in rec['core_cycle']]
 if sorted(cy)!=list(range(len(core))) or not all(core[x]>>y&1 for x,y in zip(cy,cy[1:]+cy[:1])):raise RuntimeError('cycle')
 for m in range(4):
  for word in itertools.product((0,1),repeat=m):
   b=list(a);cc=cut;st=expected
   for eps in word:
    n=len(b);b=[row|(1<<n) for row in b]+[0];cc|=eps<<n
    h,al,ta,e=st;st=(h,2*al+eps*h,ta+eps*(al+e),eps*h)
   if state(b,cc)!=st:raise RuntimeError(('append update',word))
   states+=1
   for U,hU,qU in [([0],1,2),([2,0],1,4),([2,4,1],3,12)]:
    graphs.append((st[1]*qU+(9-st[2])*hU,build(b,[],cc,0,U),True))
# Empty L and all small cut-labeled tournaments exercise the update independently.
emptychecks=0
for n in range(4):
 for a in ns['tournaments'](n):
  for cut in range(1<<n):
   h,al,ta,e=state(a,cut)
   for eps in (0,1):
    b=[row|(1<<n) for row in a]+[0];cc=cut|(eps<<n)
    if state(b,cc)!=(h,2*al+eps*h,ta+eps*(al+e),eps*h):raise RuntimeError('small update')
    emptychecks+=1
(P/'MENU_GRAPHS.txt').write_text(''.join(f'{h} {int(strong)} {len(a)} '+ ' '.join(map(str,a))+'\n' for h,a,strong in graphs))
print(json.dumps({'seed_states':4,'seed_cycles':4,'append_states_through_length3':states,'arbitrary_small_update_checks':emptychecks,'literal_menu_composites':len(graphs),'max_order':max(len(a) for h,a,strong in graphs)},indent=2))
