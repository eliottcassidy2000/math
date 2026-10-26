"""(C3, score function) determines T iff no prime node of the modular tree is size-balanced: compares the exception
lists of c3score (files exc{n}.txt in DATA, one class index per line) with a direct modular decomposition."""
import os, subprocess
DATA=os.environ.get("DATA",".")
def balanced_nodes(n,out):
    nodes=[]
    def scc(Vm):
        vs=sorted([v for v in range(n) if Vm>>v&1], key=lambda v: bin(out[v]&Vm).count('1'))
        comps=[]; cum=0; cur=0
        for i,v in enumerate(vs):
            cum+=bin(out[v]&Vm).count('1'); cur|=1<<v
            if cum==(i+1)*i//2: comps.append(cur); cur=0
        return comps
    def closure(Vm,S):
        ch=True
        while ch:
            ch=False
            for z in range(n):
                if Vm>>z&1 and not S>>z&1:
                    a=out[z]&S
                    if a and a!=S: S|=1<<z; ch=True
        return S
    def rec(Vm):
        if bin(Vm).count('1')==1: return
        comps=scc(Vm)
        if len(comps)>=2:
            for c in comps: rec(c)
            return
        vs=[v for v in range(n) if Vm>>v&1]; M={}
        for x in vs:
            m=1<<x
            for y in vs:
                if y!=x:
                    C=closure(Vm,(1<<x)|(1<<y))
                    if C!=Vm: m|=C
            M[x]=m
        blocks=sorted(set(M.values())); reps=[(b&-b).bit_length()-1 for b in blocks]; sizes=[bin(b).count('1') for b in blocks]
        bal=all(sum(sizes[j] for j,rj in enumerate(reps) if j!=i and out[r]>>rj&1)==sum(sizes[j] for j,rj in enumerate(reps) if j!=i and not out[r]>>rj&1) for i,r in enumerate(reps))
        nodes.append(bal)
        for b in blocks: rec(b)
    rec((1<<n)-1); return nodes
for n in [3,4,5,6,7]:
    lines=subprocess.run(['gentourng',str(n)],capture_output=True,text=True).stdout.split()
    exc=set(int(x) for x in open(os.path.join(DATA,f'exc{n}.txt')).read().split())
    pred=set()
    for idx,line in enumerate(lines,1):
        out=[0]*n; p=0
        for i in range(n):
            for j in range(i+1,n):
                if line[p]=='1': out[i]|=1<<j
                else: out[j]|=1<<i
                p+=1
        if any(balanced_nodes(n,out)): pred.add(idx)
    print(f'n={n}: (C3,score function) exceptions={len(exc)}; classes with a size-balanced prime node={len(pred)}; equal sets: {exc==pred}')
