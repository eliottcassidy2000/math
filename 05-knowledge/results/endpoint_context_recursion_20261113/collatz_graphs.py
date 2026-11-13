"""Minimax reachability in precisely three explicitly bounded graph models."""
from heapq import heappop,heappush
from fractions import Fraction
from context_head import need,save


def forward(n):return n//2 if n%2==0 else 3*n+1


def predecessors(n,mode):
    p={2*n}
    if n%6==4:p.add((n-1)//3)
    if mode=='bidirectional':p.add(forward(n))
    if mode=='literal':
        if n%4==0:p.add(n//2)
        if n%2==0:p.add(3*n+1)
    return p


def allowed(a,b,mode):
    if forward(a)==b:return True
    if mode=='bidirectional':return forward(b)==a
    if mode=='literal':return (a%2==0 and b==2*a) or (a>1 and a%2 and a%3==1 and b==(a-1)//3)
    return False


def minimax(cap,mode):
    inf=cap+1;cost=[inf]*(cap+1);nxt=[0]*(cap+1);cost[1]=1;heap=[(1,1)]
    while heap:
        c,n=heappop(heap)
        if cost[n]!=c:continue
        for p in predecessors(n,mode):
            if 1<=p<=cap:
                cc=max(c,p)
                if cc<cost[p]:cost[p]=cc;nxt[p]=n;heappush(heap,(cc,p))
    return cost,nxt


def run():
    cap=1_000_000;target=100_000;rows=[];costs={}
    # Exhaustively audit the reverse relation against the literal forward definition.
    for mode in ('ordinary','bidirectional','literal'):
        for n in range(1,301):
            for p in range(1,301):need((p in predecessors(n,mode))==allowed(p,n,mode),'reverse edge guard')
        d,nxt=minimax(cap,mode);costs[mode]=d
        missing=[n for n in range(1,target+1) if d[n]>cap]
        witnessed=[n for n in range(1,target+1) if d[n]<=cap]
        worst=max(witnessed,key=lambda n:Fraction(d[n],n))
        paths={}
        for n in sorted(set([3,7,27,31,703,9663,worst])):
            if d[n]>cap:continue
            path=[n];seen={n}
            while n!=1:
                m=nxt[n];need(allowed(n,m,mode) and m not in seen,'witness is simple valid path');seen.add(m);path.append(m);n=m
            paths[str(path[0])]=path
            need(max(path)==d[path[0]],'bottleneck witness')
        rec=dict(mode=mode,cap=cap,starts=target,reached=len(witnessed),unresolved_count=len(missing),first_unresolved=missing[:20],largest_peak_ratio=dict(start=worst,peak=d[worst],ratio=str(Fraction(d[worst],worst))),paths=paths)
        rows.append(rec);print('COLLATZ',mode,'REACHED',len(witnessed),'FIRST_MISSING',missing[:5],'WORST',rec['largest_peak_ratio'],flush=True)
    # A functional graph's undirected components do not shorten the minimax threshold.
    need(costs['ordinary']==costs['bidirectional'],'all finite-cap basin connectivity equivalence')
    need(all(b<=a for a,b in zip(costs['ordinary'][1:],costs['literal'][1:])),'literal contains ordinary arrows')
    improved=[n for n in range(1,target+1) if costs['literal'][n]<costs['ordinary'][n]]
    # A separate actual path strategy; the induced-cap misses above remain misses.
    def greedy(n):
        if n%2==0:return n//2
        if n>1 and n%3==1:return (n-1)//3
        return 3*n+1
    known={1:(0,1)};failures=[]
    for start in range(2,target+1):
        n=start;path=[];seen=set()
        while n not in known and n not in seen and len(path)<10000:
            seen.add(n);path.append(n);n=greedy(n)
        if n not in known:failures.append(start);continue
        for v in reversed(path):
            w=greedy(v);need(allowed(v,w,'literal'),'literal greedy edge')
            steps,peak=known[w];known[v]=steps+1,max(v,peak)
    greedy_known=[n for n in range(1,target+1) if n in known]
    peak_start=max(greedy_known,key=lambda n:known[n][1]);step_start=max(greedy_known,key=lambda n:known[n][0])
    greedy_data=dict(starts=target,certified=len(greedy_known),unresolved=failures,visited_rank_states=len(known),largest_peak=dict(start=peak_start,peak=known[peak_start][1]),longest_path=dict(start=step_start,steps=known[step_start][0]))
    forced=[]
    for K in range(2,61,2):
        chain=[3**j*2**(K-j)-1 for j in range(K)]
        need(all(x>1 and x%2 and x%3!=1 for x in chain),'inverse guard unavailable')
        need(all(3*x+1==2*y for x,y in zip(chain,chain[1:])),'forced odd chain')
        need(3*chain[-1]+1==2*3**K-2,'forced native peak')
        forced.append(dict(K=K,start=chain[0],forced_peak=2*3**K-2))
    save('collatz_graphs',dict(status='FINITE-EXACT',scope='Positive vertices1..1,000,000, induced graphs, exact minimum maximum vertex on a path to1. Claims about1..100,000; unresolved means no path within this cap, not divergence. Separate deterministic literal strategy with10000-step cap and30 algebraic forced-excursion controls.',reverse_guard_checks=270000,ordinary_bidirectional_costs_identical=True,literal_strict_improvements=len(improved),rows=rows,greedy=greedy_data,forced_excursions=forced))
    print('LITERAL GREEDY',greedy_data,flush=True)


if __name__=='__main__':run()
