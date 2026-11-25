from pathlib import Path
from runpy import run_path
from sys import argv
from struct import pack
base=run_path(str(Path(__file__).with_name('gewirtz-transport.py')))
globals().update({k:v for k,v in base.items() if not k.startswith('__')})
def descriptor(csets):
    bits=[sum(1<<i for i in L) for L in csets]
    index={m:i for i,m in enumerate(bits)}
    def act(p):
        return bytes(index[sum(1<<p[i] for i in a)] for a in csets)
    return act
point=descriptor([[j for j,h in enumerate(orbit) if h>>i&1] for i in range(21)])
line=descriptor([[j for j,h in enumerate(orbit) if not h&v] for v in line_bits])
binary_path=argv[1] if len(argv)>1 else 'gewirtz-seven.bin'
ders=[p for p in G1 if all(i!=j for i,j in enumerate(p))]
assert len(ders)==19080
left=set(ders); rep=[]; sizes=[]
gen=G+[pl_perm(frob)]
inv=lambda p:bytes(p.index(i) for i in range(56))
while left:
    g=min(left);rep.append(g);left.remove(g);V=[g]; seen={g}
    for p in V:
        for h in gen:
            k=inv(h)
            q=bytes(h[p[k[i]]] for i in range(56))
            if q not in seen:seen.add(q);V.append(q);left.remove(q)
    sizes.append(len(V))
assert sorted(sizes)==[2520,2880,2880,2880,2880,5040]

from re import findall
from pathlib import Path
def check(ids, swaps_only):
    T=[bytes(range(56))]+[ders[i] for i in ids]
    assert all(len({p[i] for p in T})==7 for i in range(56))
    which=[(point(p),line(p)) for p in T]
    assert all(all(sum(bool(line_bits[l]&(1<<p[i])) for p in A)%2 for l in range(21))
               for A in zip(*which) for i in range(21))
    print('seed',ids,'subgroup closed',all(bytes(p[q[i]] for i in range(56)) in T for p in T for q in T),'Fano-source counts', [sum(len({p[i] for p in A})==7 for i in range(21)) for A in zip(*which)])
    count=[{p[i] for p in T} for i in range(56)]
    g1=set(G1)
    options=[p for p in G2 if (not swaps_only or p not in g1) and all(p[i] not in count[i] for i in range(56))]
    cells=[i*56+j for i in range(56) for j in range(56) if j not in count[i]]
    avail={c:0 for c in cells}; cover=[]
    for i,p in enumerate(options):
        D=[56*j+k for j,k in enumerate(p)]
        cover.append(D)
        for c in D:avail[c]|=1<<i
    assert len(options)==(819 if swaps_only else 1407)
    assert min(a.bit_count() for a in avail.values())==(8 if swaps_only else 24)
    print('exact-cover remaining rows',len(options),'cells',len(cells),'minimum cell options',min(a.bit_count() for a in avail.values()),flush=True)
    conf=[0]*len(options)
    for i,D in enumerate(cover):
        for c in D:conf[i]|=avail[c]
    visits=0
    def exact(C,A):
        nonlocal visits
        visits+=1
        if not C:return []
        b=0; few=len(options)+1
        for c in C:
            x=A & avail[c]; n=x.bit_count()
            if n<few:
                few=n;b=x
                if not b:return None
                if n==1:break
        while b:
            x=b&-b;b^=x;i=x.bit_length()-1
            next=exact(C.difference(cover[i]),A&~conf[i])
            if next is not None:return [i]+next
        return None
    ans=exact(set(cells),(1<<len(options))-1)
    print('completion',ans,'visits',visits)
    assert ans is None
ids=[]
for f in argv[1:]:
    ids.extend(tuple(map(int,x.split())) for x in findall(r'^compatible seven-set:((?: \d+){6})$',Path(f).read_text(),flags=8))
assert len(ids)==2
for a in ids:
    for swaps_only in (True,False):check(a,swaps_only)
