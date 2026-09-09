"""Independent literal-matching and block-state referee; no producer imports."""
from collections import Counter,defaultdict
from functools import lru_cache
from itertools import combinations,product
from pathlib import Path
import hashlib,json,sys
import networkx as nx

sys.stdout.reconfigure(newline='\n')
G=Counter()
def need(ok,label):
    G[label]+=1
    if not ok:raise ArithmeticError(label)
def brute_matching(edges,n):
    adj=defaultdict(list)
    for (u,v),w in edges.items():adj[u].append((v,w));adj[v].append((u,w))
    @lru_cache(None)
    def visit(mask):
        if not mask:return 0,()
        bit=mask&-mask;u=bit.bit_length()-1;rest=mask^bit;best=visit(rest)
        for v,w in adj[u]:
            if rest>>v&1:
                score,chosen=visit(rest^(1<<v));candidate=(score+w,chosen+((min(u,v),max(u,v)),))
                if candidate[0]>best[0]:best=candidate
        return best
    return visit((1<<n)-1)
def compiled(cells,profiles,exact=None):
    if exact:
        a,r=exact
        if r<0 or r>=len(profiles[a]):return None
        profiles=dict(profiles);profiles[a]=profiles[a][:r+1]
    ids={};edges={};pc={};degree={}
    def node(label):
        if label not in ids:ids[label]=len(ids)
        return ids[label]
    def edge(u,v,w,ports=0,line=None):
        key=tuple(sorted((node(u),node(v))));edges[key]=w;pc[key]=ports;degree[key]=int(exact is not None and line==exact[0])
    slots={};inccharge={};constant=0
    for a,R in profiles.items():
        cap=len(R)-1;slots[a]=[('slot',a,j) for j in range(cap)]
        for v in slots[a]:node(v)
        if cap==2:
            L=min(R[1]-R[2],R[0]-R[1]);leaf=R[1]-R[2]-L;internal=R[0]-R[2]-2*L
            need(internal>=leaf>=0,'arbitrary profile gadget signs')
            edge(slots[a][0],slots[a][1],internal)
            for v in slots[a]:edge(v,('leaf',a),leaf)
            constant+=R[2]+2*L;inccharge[a]=-L
        elif cap==1:constant+=R[0];inccharge[a]=R[1]-R[0]
        else:constant+=R[0];inccharge[a]=0
    ports=[]
    for p,inc,reward in cells:
        u=('cell',p,0);v=('cell',p,1);ports.extend((node(u),node(v)));edge(u,v,0,2)
        for j,a in enumerate(inc):
            for slot in slots[a]:edge(u if j==0 else v,slot,inccharge[a]+(reward if j==0 else 0),1,a)
        if len(inc)==1:edge(v,('free',p),0,1)
    N=sum(map(abs,edges.values()))+1
    adjusted={e:w+N*degree[e] for e,w in edges.items()}
    M=sum(map(abs,adjusted.values()))+1
    score,matching=brute_matching({e:w+M*pc[e] for e,w in adjusted.items()},len(ids))
    covered={v for e in matching for v in e}
    need(set(ports)<=covered,'literal maximum matching covers all original cell ports')
    chosen=[];deg={a:0 for a in profiles}
    for p,inc,reward in cells:
        pair=tuple(sorted((ids[('cell',p,0)],ids[('cell',p,1)])))
        if pair not in matching:
            chosen.append(p)
            for a in inc:deg[a]+=1
    true=sum(w for p,inc,w in cells if p in chosen)+sum(profiles[a][deg[a]] for a in profiles)
    recovered=score-M*len(ports)-N*sum(degree[e] for e in matching)+constant
    need(true==recovered,'ordinary matching objective decodes to exact original-cell objective')
    if exact and deg[exact[0]]!=exact[1]:return None
    return true
def literal(cells,profiles,exact=None):
    best=None
    for bits in product((0,1),repeat=len(cells)):
        d=Counter();score=0
        for bit,(p,inc,w) in zip(bits,cells):
            if bit:score+=w;d.update(inc)
        if any(d[a]>=len(R) for a,R in profiles.items()):continue
        if exact and d[exact[0]]!=exact[1]:continue
        score+=sum(R[d[a]] for a,R in profiles.items());best=score if best is None else max(best,score)
    return best
def overfull(points):
    buckets=defaultdict(list)
    for i,(r,c) in enumerate(points):
        for a in (1,-1,2):buckets[a,c-a*r].append(i)
    return {a:tuple(p) for a,p in sorted(buckets.items()) if len(p)>=3}
def direct_repair(lines):
    cells=sorted(set().union(*map(set,lines))) if lines else []
    best=len(cells)
    for bits in range(1<<len(cells)):
        loss=len(cells)-bits.bit_count()
        if loss>=best:continue
        kept={p for j,p in enumerate(cells) if bits>>j&1}
        if all(len(kept.intersection(a))<=2 for a in lines):best=loss
    return best
def block_solver(lines,root_choice=0):
    X=nx.Graph()
    for a,line in enumerate(lines):
        for p in line:X.add_edge(('L',a),('P',p))
    if not X:return 0,[],0
    edgeblocks=list(nx.biconnected_component_edges(X));blocks=[];appear=defaultdict(list)
    for es in edgeblocks:
        inc=defaultdict(list)
        for u,v in es:
            line,cell=(u,v) if u[0]=='L' else (v,u);inc[cell[1]].append(line[1])
        block={p:tuple(sorted(aa)) for p,aa in inc.items()};blocks.append(block)
        for p in block:appear['P',p].append(len(blocks)-1)
        for a in set(a for aa in block.values() for a in aa):appear['L',a].append(len(blocks)-1)
    separators={v for v,bb in appear.items() if len(bb)>1}
    need(separators==set(nx.articulation_points(X)),'independent biconnected components retain both separator types')
    forest=nx.Graph();forest.add_nodes_from(('B',i) for i in range(len(blocks)))
    for s in separators:
        for b in appear[s]:forest.add_edge(s,('B',b))
    need(nx.is_forest(forest),'unmerged block-cut incidence is a forest')
    @lru_cache(None)
    def sep_cost(s,parent,state):
        children=[b for b in appear[s] if b!=parent]
        if s[0]=='P':return sum(bag(b,s,state) for b in children)
        values={0:0}
        for b in children:
            new={}
            for used,cost in values.items():
                for r in range(3-used):
                    new[used+r]=min(new.get(used+r,float('inf')),cost+bag(b,s,r))
            values=new
        return min(c for r,c in values.items() if r<=2-state)
    @lru_cache(None)
    def bag(b,parent,state):
        inc=blocks[b];cells=sorted(inc);best=float('inf')
        for bits in product((0,1),repeat=len(cells)):
            chosen=dict(zip(cells,bits));degree=Counter(a for p in cells if chosen[p] for a in inc[p])
            if any(v>2 for v in degree.values()):continue
            if parent and parent[0]=='P' and chosen[parent[1]]!=state:continue
            if parent and parent[0]=='L' and degree[parent[1]]!=state:continue
            # Charge shared cells in their nearest-to-root owning block, never in a separator node.
            cost=sum(1-chosen[p] for p in cells if parent!=('P',p))
            nodes=[('P',p) for p in cells]+[('L',a) for a in set(a for aa in inc.values() for a in aa)]
            for s in nodes:
                if s in separators and s!=parent:cost+=sep_cost(s,b,chosen[s[1]] if s[0]=='P' else degree[s[1]])
            best=min(best,cost)
        return best
    total=0
    for component in nx.connected_components(forest):
        bs=sorted(v[1] for v in component if v[0]=='B');root=root_choice if root_choice in bs else bs[0]
        total+=bag(root,None,0)
    hs=[sum(len(aa)>=3 for aa in block.values()) for block in blocks]
    return total,hs,len(blocks)
def boards(n):
    for rows in product(tuple(combinations(range(n),2)),repeat=n):
        if Counter(c for row in rows for c in row)==Counter({c:2 for c in range(n)}):yield [(r,c) for r,row in enumerate(rows) for c in row]
def abstract(v):
    candidates=list(combinations(range(v),2))+list(combinations(range(v),3))
    for bits in range(1<<len(candidates)):
        shared=[s for j,s in enumerate(candidates) if bits>>j&1]
        if max(Counter(a for s in shared for a in s).values(),default=0)>3:continue
        if max(Counter(p for s in shared for p in combinations(s,2)).values(),default=0)>1:continue
        ls=[[] for _ in range(v)]
        for p,s in enumerate(shared):
            for a in s:ls[a].append(p)
        p=len(shared)
        for line in ls:
            while len(line)<3:line.append(p);p+=1
        yield ls
def main():
    comparisons=0
    patterns=[[(0,(0,),-4),(1,(0,),0),(2,(0,),9)],[(0,(0,1),-4),(1,(0,1),9),(2,(1,),1)]]
    for R in product((-7,2,11),repeat=3):
        for cells in patterns:
            prof={0:R}
            if any(1 in inc for p,inc,w in cells):prof[1]=(3,-5,7)
            for r in (None,0,1,2):
                exact=None if r is None else (0,r)
                need(compiled(cells,prof,exact)==literal(cells,prof,exact),'literal weighted matching equals exhaustive retained-subset optimum');comparisons+=1
    need(compiled([],{0:(0,4,-3)},(0,2)) is None,'impossible exact parent degree rejected')
    need(compiled([(0,(0,1),10)],{0:(0,0),1:(0,)})==0,'coverage dominance excludes half-retained original cell')
    # With reward profile (0,1,0), two separate leaves give4 instead of the required free-two reward2.
    wrong,_=brute_matching({(0,2):2,(1,3):2,(0,1):2},4)
    right,_=brute_matching({(0,2):2,(1,2):2,(0,1):2},3)
    need((wrong,right)==(4,2),'shared-leaf identity is load-bearing for nonconvex profiles')
    census={}
    for n in range(2,5):
        count=0
        for pts in boards(n):
            ls=list(overfull(pts).values());value,hs,b=block_solver(ls)
            need(value==direct_repair(ls),'independent exact separator DP equals original retained subsets');count+=1
        need(count=={2:1,3:6,4:90}[n],'independent full Cartesian two-regular board universe');census[n]=count
    abstract_counts={}
    for v in range(5):
        count=0
        for ls in abstract(v):
            expected=direct_repair(ls);value,hs,b=block_solver(ls)
            need(value==expected,'all abstract original-cell systems agree')
            for root in range(b):need(block_solver(ls,root)[0]==expected,'every abstract root preserves exact original cost')
            count+=1
        need(count=={0:1,1:1,2:2,3:9,4:96}[v],'complete independent abstract incidence universe');abstract_counts[v]=count
    theta=[(6,12),(12,12),(9,15),(10,8),(0,0),(8,10),(7,17),(11,11),(13,19),(15,18),(14,28)]
    first=overfull(theta);common=first[1,0];phi={r:99 for r in range(3)}
    for bits in range(1<<11):
        kept={i for i in range(11) if bits>>i&1};r=len(kept.intersection(common))
        if r<3 and all(len(kept.intersection(line))<=2 for label,line in first.items() if label!=(1,0)):phi[r]=min(phi[r],11-len(kept))
    need(phi=={0:4,1:3,2:2},'complete independent theta contribution table')
    family=[]
    for k in range(1,7):
        pts=[(1000**j*(r+100),1000**j*(c+100)) for j in range(k) for r,c in theta]
        named=overfull(pts);ls=list(named.values());value,hs,b=block_solver(ls)
        need(len(pts)==11*k and len(ls)==5*k+1 and len(named[1,0])==3*k,'all actual geometric family lines and cells')
        need(value==4*k-2 and b==7*k and Counter(hs)==Counter({0:6*k,2:k}),'exact geometric optimum and unmerged parameter')
        high=[p for p,n in Counter(p for line in ls for p in line).items() if n>=3]
        valid=sum(all(len(set(p for p,keep in zip(high,choice) if keep).intersection(line))<=2 for line in ls) for choice in product((0,1),repeat=len(high)))
        need(valid==2**k*(1+k+k*(k-1)//2),'exact old retained exceptional branch count')
        family.append(dict(k=k,tau=value,blocks=b,old_valid_branches=valid,new_bound=3*sum(2**h for h in hs)))
    one=[(100,100),(101,102),(102,104),(103,103),(104,102)]
    double=one+[(1000*r,1000*c) for r,c in one];ls=list(overfull(double).values())
    need(direct_repair(ls)==block_solver(ls)[0]==3,'actual ten-point nonconvex original repair')
    other=[(0,1,2),(2,3,4)];phi2={r:99 for r in range(3)}
    for bits in range(32):
        kept={i for i in range(5) if bits>>i&1};r=len(kept.intersection({0,3}))
        if all(len(kept.intersection(line))<=2 for line in other):phi2[r]=min(phi2[r],5-len(kept))
    message=[min(phi2[r] for r in range(3-t)) for t in range(3)]
    need(phi2=={0:2,1:2,2:1} and message==[1,2,2],'literal nonconvex branch and continuation profiles')
    cert=dict(status='PASS independent analytic referee controls',matching_comparisons=comparisons,
      matching_oracle='literal exhaustive ordinary matching by vertex-mask recursion; no blossom call',
      board_counts=census,abstract_counts=abstract_counts,theta_costs=phi,nonconvex_costs=phi2,nonconvex_message=message,
      family=family,gates=dict(sorted(G.items())),total_gates=sum(G.values()),scope='Finite controls support the separately proved all-board theorem; no producer imports and no general polynomial-time hyperedge claim.')
    dest=Path(__file__).resolve().parent
    if dest.name=='04-computation':dest=dest.parent/'05-knowledge/results'
    path=dest/(Path(__file__).stem+'_certificate.json');path.write_bytes((json.dumps(cert,sort_keys=True,indent=2)+'\n').encode())
    print('Independent unmerged no3 line/cell states: PASS')
    print('216 literal ordinary-matching comparisons;97 complete geometric boards;109 complete abstract systems, every root')
    print('Actual theta family k1..6; exact nonconvex continuation; mandatory coverage and shared-leaf hostiles: PASS')
    print('Always-active exact gates:',sum(G.values()))
    print('Certificate:',path.name)
if __name__=='__main__':main()
