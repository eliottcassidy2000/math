"""Check all physical edges and the unique cycle. --exhaustive settles the missing edges."""
import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, product
from math import isqrt, prod
from pathlib import Path
import subprocess
import tempfile

from square36 import all_factors, components

HERE = Path(__file__).resolve().parent
VERTICES = {100*v+s for v in range(1,37) for s in (-1,1)}
ARITHMETIC = {(a,b) for a,b in combinations(sorted(VERTICES),2) if isqrt(a+b)**2 == a+b}
DEEP_EDGES = {(199,701), (299,3301), (399,3201), (501,3099), (599,3001), (699,1801)}


def path_count(rows):
    n=len(rows)
    assert all(r>=0 and r<1<<n and not (r>>i&1) for i,r in enumerate(rows))
    assert all(((rows[i]>>j)^(rows[j]>>i))&1 for i,j in combinations(range(n),2))
    full=(1<<n)-1
    states=[[0]*n for _ in range(full+1)]
    for i in range(n):
        states[1<<i][i]=1
    for s,endpoints in enumerate(states):
        for i,count in enumerate(endpoints):
            b=rows[i]&~s
            while count and b:
                v=b&-b
                states[s|v][v.bit_length()-1]+=count
                b-=v
    return sum(states[full])


def physical_edges():
    edges={}
    for line in (HERE/'physical_edges.rows').read_text().splitlines():
        if line.startswith('#'):
            continue
        a,b,u,v,*rows=map(int,line.split())
        assert (a,b) in ARITHMETIC and u!=v and (a,b) not in edges
        assert path_count(rows)==a
        rows[u]^=1<<v
        rows[v]^=1<<u
        assert path_count(rows)==b
        edges[a,b]=len(rows)
    return edges


@dataclass
class Core:
    h: int
    q: int
    start: list
    end: list
    starts: list
    ends: list
    flips: set

    @property
    def n(self):
        return len(self.start)


SINGLE=Core(1,2,[1],[1],[1],[1],set())


def strong_catalogue():
    bank=defaultdict(list)
    catalogue_sizes={}
    with tempfile.TemporaryDirectory() as directory:
        binary=Path(directory)/'strong_ports'
        subprocess.run(['c++','-O3','-std=c++17',str(HERE/'strong_ports.cpp'),'-o',str(binary)],check=True)
        graphs=b'&BP_\n'  # directed triangle in digraph6
        for n in range(3,14):
            catalogue_sizes[n]=len(graphs.splitlines())
            records=subprocess.run([binary,'read'],input=graphs,capture_output=True,check=True).stdout
            for line in records.splitlines():
                code,h,q,*out=line.split()
                assert len(code)>2 and code[1]-63==n
                out=list(map(int,out))
                arrays=[out[i*n:(i+1)*n] for i in range(4)]
                bank[int(h)].append(Core(int(h),int(q),*arrays,set(out[4*n:])))
            if n==13:
                break
            raw=subprocess.run([binary,'extend'],input=graphs,capture_output=True,check=True).stdout
            canonical=subprocess.run(['labelg','-q'],input=raw,capture_output=True,check=True).stdout
            graphs=b'\n'.join(sorted(set(canonical.splitlines())))+b'\n'
    print('strong classes with H <= 699 by order',catalogue_sizes)
    assert list(catalogue_sizes.values())==[1,1,6,35,353,6008,17748,12194,5369,1150,23]
    return bank


def factor_words(h, bank):
    if h==1:
        yield ()
        return
    for d in bank:
        if h%d==0:
            for tail in factor_words(h//d,bank):
                yield (d,)+tail


def between(left, middle, right, delta):
    """Yield exact possible numbers of intervening singleton SCCs."""
    inner=prod(k.q for k in middle)
    for u,a in enumerate(left.starts):
        for v,b in enumerate(right.ends):
            base=a*b*inner
            if not middle and base-left.end[u]*right.start[v]==delta:
                yield 0
            if delta%base:
                continue
            ratio=delta//base
            if ratio&(ratio-1):
                continue
            k=ratio.bit_length()-1
            if k or middle:
                yield k


def search_pair(lo,hi,bank):
    least=None
    core_orders=set()
    def record(order,marked_cores):
        nonlocal least
        least=order if least is None else min(least,order)
        if lo==199:
            only=[block.n for block in marked_cores if block.h>1]
            assert len(only)==1
            core_orders.add(only[0])

    minima={h:min(c.n for c in classes) for h,classes in bank.items()}
    for word in factor_words(lo,bank):
        for pos,h in enumerate(word):
            g=lo//h
            if hi%g:
                continue
            outside=sum(minima[d] for d in word[:pos]+word[pos+1:])
            for c in bank[h]:
                if hi//g in c.flips:
                    record(outside+c.n,[c])
        for i in range(len(word)+1):
            for j in range(i,len(word)+1):
                segment, exterior=word[i:j], word[:i]+word[j:]
                g=prod(exterior)
                if (hi-lo)%g:
                    continue
                outside=sum(minima[h] for h in exterior)
                for add_left,add_right in product((False,True),repeat=2):
                    if len(segment)+add_left+add_right<2:
                        continue
                    slots=([(SINGLE,)] if add_left else [])+[bank[h] for h in segment]
                    if add_right:
                        slots.append((SINGLE,))
                    for blocks in product(*slots):
                        middle=blocks[1:-1]
                        for k in between(blocks[0],middle,blocks[-1],(hi-lo)//g):
                            record(outside+sum(c.n for c in blocks)+k,blocks)
    return least,core_orders


def cycle_order(factor):
    adjacent=defaultdict(set)
    for a,b in factor:
        adjacent[a].add(b)
        adjacent[b].add(a)
    path=[min(adjacent)]
    path.append(min(adjacent[path[0]]))
    while len(path)<len(adjacent):
        path.append(next(iter(adjacent[path[-1]]-{path[-2]})))
    return path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exhaustive',action='store_true')
    args=parser.parse_args()
    positives=physical_edges()
    assert len(ARITHMETIC)==118 and len(positives)==101
    if args.exhaustive:
        bank=strong_catalogue()
        goals=(ARITHMETIC-positives.keys())|DEEP_EDGES|{(199,201)}
        for lo,hi in sorted(goals):
            order,core_orders=search_pair(lo,hi,bank)
            if (lo,hi) in positives:
                assert order==positives[lo,hi], (lo,hi,order)
            else:
                assert order is None, (lo,hi,order)
            if (lo,hi) in {(199,201),(199,701)}:
                assert core_orders==({8} if hi==201 else {10})
            print('pair',lo,hi,'least order',order,'199-core orders',sorted(core_orders))
    histogram=Counter()
    for factor in all_factors(VERTICES,set(positives),2):
        count=len(components(VERTICES,factor))
        histogram[count]+=1
        if count==1:
            assert (399,3201) in factor
            print('unique Hamilton cycle',cycle_order(factor))
    assert histogram=={1:1,2:6,3:8,4:14,5:4,6:4}
    print('arithmetic edges',len(ARITHMETIC),'physical edges',len(positives),
          '2-factors by component count',sorted(histogram.items()))


if __name__=='__main__':
    main()
