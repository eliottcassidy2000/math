"""Sharp graceful rainbow-forest completion by a consecutive star buffer."""
from pathlib import Path
from itertools import product,combinations
from collections import Counter,deque
import json,gzip,sys
HERE=Path(__file__).resolve().parent
def need(p,m):
    if not p:raise RuntimeError(m)
class DSU:
    def __init__(self,n):self.p=list(range(n))
    def root(self,a):
        while a!=self.p[a]:a=self.p[a]
        return a
    def join(self,a,b):
        a,b=self.root(a),self.root(b)
        if a==b:return False
        self.p[a]=b;return True

def initial(m,fixed,n):
    d=DSU(n+1)
    for a,b in list(fixed.values())+[(0,x) for x in range(m+1,n+1)]:need(d.join(a,b),'frozen graph remains a forest')
    return d

def complete(m,fixed,C,n=None):
    n=m+max(C)-1 if n is None else n
    d=initial(m,fixed,n);added={}
    for c in C:
        e=next(((a,a+c) for a in range(n-c+1) if d.root(a)!=d.root(a+c)),None)
        if e is None:return None
        need(d.join(*e),'greedy crossing edge');added[c]=e
    need(len({d.root(v) for v in range(n+1)})==1,'complete spanning tree')
    return added

def fiber(m,fixed,C,n):
    out=[]
    for starts in product(*(range(n-c+1) for c in C)):
        d=initial(m,fixed,n)
        if all(d.join(a,a+c) for a,c in zip(starts,C)):out.append(starts)
    # Width two: key retains every coordinate except a chosen pair.
    d=DSU(len(out));bins={}
    for i,t in enumerate(out):
        for changed in combinations(range(len(C)),min(2,len(C))):
            key=(changed,tuple(x for j,x in enumerate(t) if j not in changed))
            if key in bins:d.join(i,bins[key])
            else:bins[key]=i
    sizes=sorted(Counter(d.root(i) for i in range(len(out))).values())
    return out,sizes

def run():
    checks=0;failures_before=[]
    # Every partially chosen color tuple through m=5, including empty fibers.
    for m in range(1,6):
        for choices in product(*(range(-1,m-c+1) for c in range(1,m+1))):
            fixed={c:(a,a+c) for c,a in enumerate(choices,1) if a>=0};C=[c for c,a in enumerate(choices,1) if a<0]
            d=DSU(m+1)
            if not all(d.join(*e) for e in fixed.values()) or not C:continue
            need(complete(m,fixed,C) is not None,'universal consecutive-buffer construction');checks+=1
    sharp=[];stabilization=[]
    for m in range(3,31):
        F={d:(0,d) for d in range(1,m-1)};F[m-1]=(1,m)
        need(all(complete(m,F,[m],m+b) is None for b in range(m-1)),'sharp family has no earlier entry')
        e=complete(m,F,[m]);need(e=={m:(m-1,2*m-1)},'first entry is forced')
        sharp.append({'m':m,'first_buffer':m-1,'forced_edge':e[m]})
        # The last possible free edge really can enter at n=m+D=2m.
        G={m-i:(i,m) for i in range(1,m)}
        counts=[len(fiber(m,G,[m],m+b)[0]) for b in range(m+2)]
        need(counts==list(range(1,m+2))+[m+1],'sharp stabilization family')
        stabilization.append({'m':m,'counts':counts})
    old=json.loads((HERE.parent/'graceful_fixed_color_counterexample_20261019/witness.json').read_text(encoding='utf8'))
    F={int(c):tuple(e) for c,e in old['fixed_forest'].items()};C=old['free'];sequence=[]
    for n in range(12,25):
        states,sizes=fiber(12,F,C,n)
        sequence.append({'order_minus_one':n,'states':len(states),'component_sizes':sizes,'first_state':states[0] if states else None})
    need(sequence[0]['component_sizes']==[2,3],'inherited disconnected five-state fiber')
    need(all(r['component_sizes']==[r['states']] for r in sequence if r['order_minus_one']>=22),'all contact graphs connected by bound')
    # All valid new edges use an old endpoint, so state sets stabilize by m+D.
    a,_=fiber(12,F,C,23);b,_=fiber(12,F,C,24);need(a==b,'literal stabilization after m+max(C)')
    original,_=fiber(12,F,C,12);expanded,_=fiber(12,F,C,13)
    # Find a literal route connecting old states from different old components.
    def route(states,start,finish):
        prev={start:None};queue=deque([start])
        while queue:
            x=queue.popleft()
            if x==finish:
                path=[]
                while x is not None:path.append(x);x=prev[x]
                return path[::-1]
            for y in states:
                if y not in prev and sum(a!=b for a,b in zip(x,y))<=2:prev[y]=x;queue.append(y)
        return None
    pair=next((x,y) for x,y in combinations(original,2) if route(original,x,y) is None)
    bridge=route(expanded,*pair);need(bridge is not None,'one new frozen star leaf joins old components')
    out={'scope':'All-order proof is in GRACEFUL.md; finite controls enumerate every partial-color forest through m=5, sharp families m3..30, and the inherited m12 obstruction with star padding through n24. Abstract tree shape may change.','partial_forest_controls':checks,'sharp_family_controls':sharp,'sharp_stabilization_controls':stabilization,'frozen_m12_padding':sequence,'one_leaf_bridge':{'free_colors':C,'edge_start_tuples':bridge}}
    tag='' if __debug__ else '.optimized';(HERE/('graceful_padding'+tag+'.json.gz')).write_bytes(gzip.compress((json.dumps(out,separators=(',',':'))+'\n').encode(),mtime=0))
    print('controls',checks,'padding',sequence,flush=True)
if __name__=='__main__':run()
