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
index={p:i for i,p in enumerate(ders)}
with open(binary_path,'wb') as fp:
    fp.write(pack('<II',len(ders),len(rep)))
    fp.write(pack('<6I',*(index[p] for p in rep)))
    fp.write(pack('<21I',*line_bits))
    for p in ders:fp.write(p+point(p)+line(p))
print('exact group input',len(ders),'class sizes',sizes)
