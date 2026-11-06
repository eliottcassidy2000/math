"""Independent full-path enumeration and actual rooted insertion maps."""
from itertools import combinations,product
from pathlib import Path
from collections import deque
import json,random
CORE=[102,0,10,83,135,30,182,47]
OLD=[0,2,3,4,6,7]
VISITS=0

def need(b,s):
    if not b:raise RuntimeError(s)

def valid(a):
    n=len(a)
    need(all(0<=x<1<<n and not(x>>i&1) for i,x in enumerate(a)),'range/loop')
    need(all((a[i]>>j&1)+(a[j]>>i&1)==1 for i in range(n) for j in range(i)),'orientation')

def induced(a,vs):return [sum(1<<j for j,w in enumerate(vs) if a[v]>>w&1) for v in vs]
def tournament(n,code):
    a=[0]*n
    for bit,(i,j) in enumerate(combinations(range(n),2)):
        if code>>bit&1:a[i]|=1<<j
        else:a[j]|=1<<i
    return a

def kernel(body,color):
    # Independent convention: body first, followed by old labels 0,2,3,4,6,7.
    n=len(body);a=body[:]+[0]*6
    for i,u in enumerate(OLD):
        for j,v in enumerate(OLD):
            if CORE[u]>>v&1:a[n+i]|=1<<(n+j)
        for j in range(n):
            win=(color>>j&1) if u==3 else (CORE[u]>>5&1)
            if win:a[n+i]|=1<<j
            else:a[j]|=1<<(n+i)
    valid(a);return a

def reach(a,r):
    seen={r};todo=[r]
    while todo:
        v=todo.pop()
        for w in range(len(a)):
            if a[v]>>w&1 and w not in seen:seen.add(w);todo.append(w)
    return seen

def paths(a,root=None):
    global VISITS
    full=(1<<len(a))-1
    def dfs(word,mask):
        global VISITS
        VISITS+=1
        if mask==full:yield tuple(word);return
        todo=a[word[-1]]&~mask
        while todo:
            bit=todo&-todo;todo-=bit
            yield from dfs(word+[bit.bit_length()-1],mask|bit)
    for r in (range(len(a)) if root is None else [root]):yield from dfs([r],1<<r)

def starts(a):
    ans=[0]*len(a)
    for p in paths(a):ans[p[0]]+=1
    return ans

def injection(a,subset,r):
    n=len(a);S=set(subset)
    queue=deque([(r,)]);seen={r};entry=None
    while queue:
        p=queue.popleft();v=p[-1]
        if v in S:entry=p;break
        for w in range(n):
            if a[v]>>w&1 and w not in seen:seen.add(w);queue.append(p+(w,))
    need(entry is not None,'entry path');v=entry[-1]
    need(set(entry[:-1]).isdisjoint(S),'first entry')
    # Freeze a vertex schedule independently of the initial Hamilton word.
    K=S|set(entry);schedule=[]
    while len(K)<n:
        w=next((w for w in range(n) if w not in K and any(a[x]>>w&1 for x in K)),None)
        need(w is not None,'root reachability supplies next vertex');schedule.append(w);K.add(w)
    sub=induced(a,subset);all_outputs=set();count=0
    for initial in paths(sub,subset.index(v)):
        original=tuple(subset[i] for i in initial);word=list(entry[:-1]+original)
        for w in schedule:
            j=max(i for i,x in enumerate(word) if a[x]>>w&1)
            need(j==len(word)-1 or a[w]>>word[j+1]&1,'root-preserving legal gap')
            word.insert(j+1,w)
        need(word[0]==r and len(word)==n and len(set(word))==n,'root/output vertex set')
        need(all(a[x]>>y&1 for x,y in zip(word,word[1:])),'directed full output')
        need(tuple(x for x in word if x in S)==original,'inverse deletion recovers seed path')
        all_outputs.add(tuple(word));count+=1
    need(len(all_outputs)==count,'injectivity')
    return count,len(entry),len(schedule)

def main():
    here=Path(__file__).resolve().parent
    frozen=json.loads((here/'PRODUCER_COLORED_KERNELS.json').read_text());records={};literal=0
    for case in frozen['kernels']:
        need(case['body_rows']==tournament(3,case['body_code']),'every labeled body code')
        a=kernel(case['body_rows'],case['color']);ss=starts(a);need(all(len(reach(a,r))==9 for r in range(9)),'strong nine-kernel')
        order=[(label[1] if label[0]=='b' else 3+OLD.index(label[1])) for label in case['labels']]
        need(induced(a,order)==case['rows'],'archived nine-kernel rows')
        need([ss[i] for i in order]==case['starts'],'archived every-root counts')
        need(min(ss)==case['minimum']>=27 and sum(ss)==case['H'],'nine-kernel floor')
        records[(case['body_code'],case['color'])]=(a,ss);literal+=1
    need(set(records)==set(product(range(8),repeat=2)),'all labeled 3-body cases')
    reps=[]
    for body,colors in [([6,4,0],range(8)),([2,4,1],(0,1,3,7))]:
        reps.append([min(starts(kernel(body,c))) for c in colors])
    need(reps==[[27,35,31,39,29,37,33,41],[45,51,57,63]],'displayed representative minima')
    four=[];embeddings=0
    for case in json.loads((here/'PRODUCER_KERNELS.json').read_text())['kernels']:
        seed=case['seed'];a=kernel([30,28,24,16,0],seed);ss=starts(a)
        # The archived labeling keeps old0,2,3,4, then the body, then old6,7.
        order=[5,6,7,8,0,1,2,3,4,9,10]
        need(induced(a,order)==case['rows'] and [ss[i] for i in order]==case['starts'],'four-word kernel rows/counts')
        need(min(ss)==case['minimum_start'] and sum(ss)==case['H'],'four-word minima')
        need(all(len(reach(a,r))==11 for r in range(11)),'strong eleven-kernel')
        for k in range(26):
            m=5+2*k;body=[sum(1<<j for j in range(i+1,m)) for i in range(m)];big=kernel(body,(seed<<k)+(1<<k)-1)
            need(induced(big,list(range(k,k+5))+list(range(m,m+6)))==a,'all-scale retained middle kernel');embeddings+=1
        four.append({'seed':seed,'H':sum(ss),'minimum_start':min(ss)});literal+=1
    need([c['seed'] for c in four]==[10,14,17,21] and [c['minimum_start'] for c in four]==[119,127,133,141] and [c['H'] for c in four]==[7213,7221,7227,7235],'four displayed finite kernels')
    # Hostile body-size boundary: enumerate all colored bodies of order 2.
    small=[]
    for code in range(2):
        for color in range(4):
            a=kernel(tournament(2,code),color);ss=starts(a)
            small.append({'code':code,'color':color,'minimum_start':min(ss),'starts':ss})
    need(min(x['minimum_start'] for x in small)<27,'order-three hypothesis is substantive')
    # Actual maps, not just inequalities, on every rooted small input/kernel.
    maps=images=outside=nonempty_schedules=0;max_entry=0;host_graphs=0
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            a=tournament(n,code);valid(a);host_graphs+=1
            roots=[r for r in range(n) if len(reach(a,r))==n]
            for mask in range(1,1<<n):
                subset=[i for i in range(n) if mask>>i&1]
                for r in roots:
                    count,entry,schedule=injection(a,subset,r);images+=count;maps+=1;outside+=r not in subset;nonempty_schedules+=schedule>0;max_entry=max(max_entry,entry)
    print(json.dumps({'status':'PASS','literal_primary_kernels':literal,'nine_kernel_cases':64,'representative_minima':reps,'four_word_kernels':four,'retained_middle_embeddings':embeddings,'scales_checked':[0,25],'small_body_boundary':small,'injection_host_graphs_through5':host_graphs,'actual_injection_maps':maps,'injected_paths':images,'outside_root_maps':outside,'maps_with_new_vertex_insertion':nonempty_schedules,'longest_entry_path':max_entry,'directed_prefix_visits':VISITS,'scope':'Only tournaments containing a specified strong kernel/interface; unrelated constructions and F remain open.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
