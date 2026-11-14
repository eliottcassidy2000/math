# check that the Monsky nullity is additive over strong components, for every non-strong class of orders 2-8
# needs monsky.py (this folder) and classgen lists cls_n.txt in the working directory
src=open('monsky.py').read().split('strong0=set()')[0]
ns={}; exec(src,ns); decode=ns['decode']; nullity=ns['nullity']
def sccs(out):
    n=len(out); reach=[0]*n
    for v in range(n):
        seen=fr=1<<v
        while fr:
            nf=0
            for w in range(n):
                if fr>>w&1: nf|=out[w]
            nf&=~seen; seen|=nf; fr=nf
        reach[v]=seen
    cs=[]; done=0
    for v in range(n):
        if done>>v&1: continue
        c=sum(1<<w for w in range(n) if reach[v]>>w&1 and reach[w]>>v&1); cs.append(c); done|=c
    return cs
def sub(out,mask):
    idx=[v for v in range(len(out)) if mask>>v&1]; pos={v:i for i,v in enumerate(idx)}
    return [sum(1<<pos[w] for w in idx if out[v]>>w&1) for v in idx]
bad=tot=0
for n in range(2,9):
    for line in open('cls_%d.txt'%n):
        out=decode(int(line.split()[0]),n); cs=sccs(out)
        if len(cs)==1: continue
        tot+=1
        if nullity(out)!=sum(nullity(sub(out,c)) for c in cs): bad+=1
print('non-strong classes checked (orders 2-8):',tot,' additivity violations:',bad)
