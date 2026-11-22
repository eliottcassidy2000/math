"""All 29-step transitive inflations of the primitive 27-to-29 arcs at exponent 61."""
from collections import Counter
from functools import cache
from itertools import combinations, permutations
from math import comb, prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'chanpak_arc_fiber_grammar_20261122'))
from fibers import automorphisms, bounded_classes, paths, strong, strong_cap

CENTER = 2**60*(2**61-1)


def edit(rows, a, b):
    other=list(rows)
    other[a]^=1<<b
    other[b]^=1<<a
    return tuple(other)


def seeds():
    for n, rows, h in bounded_classes(strong_cap(27), 27):
        if h != 27:
            continue
        arcs={arc for arc in combinations(range(n),2) if paths(edit(rows,*arc))==29}
        if not arcs:
            continue
        autos=[perm for perm in permutations(range(n))
               if all((rows[i]>>j&1)==(rows[perm[i]]>>perm[j]&1) for i,j in combinations(range(n),2))]
        assert len(autos)==automorphisms(rows)
        while arcs:
            a,b=min(arcs)
            orbit={tuple(sorted((perm[a],perm[b]))) for perm in autos}
            assert orbit<=arcs and len(orbit)==len(autos)
            arcs-=orbit
            other=edit(rows,a,b)
            assert strong(rows) and strong(other)
            yield rows,other,(a,b)


@cache
def surjections(n,r):
    return sum((-1)**(r-j)*comb(r,j)*j**n for j in range(r+1))


def has_cycle(rows, ports):
    return any(
        all(sum(rows[i]>>j&1 for j in triple)==1 for i in triple)
        for triple in combinations(ports,3))


@cache
def signature(rows, ports):
    """Count histories of maximal runs. The vertices of ports induce a transitive tournament."""
    assert not has_cycle(rows,ports)
    full=(1<<len(rows))-1
    ordinary=full ^ sum(1<<p for p in ports)
    layer=Counter()
    result=Counter()
    for i in range(len(rows)):
        layer[ordinary & (1<<i), tuple(int(i==p) for p in ports),i]+=1
    while layer:
        nxt=Counter()
        for (mask,runs,last),mult in layer.items():
            if mask==ordinary and all(runs):
                result[runs]+=mult
            rest=rows[last]&~mask
            while rest:
                bit=rest&-rest
                rest^=bit
                j=bit.bit_length()-1
                visits=runs
                if j in ports:
                    visits=tuple(t+(j==p) for t,p in zip(runs,ports))
                nxt[mask|(bit&ordinary),visits,j]+=mult
        layer=nxt
    return result


def value(coeff,widths):
    return sum(w*prod(surjections(n,r) for n,r in zip(widths,runs)) for runs,w in coeff.items())


def displacements(parts,budget):
    if not parts:
        yield ()
        return
    for a in range(1,budget-parts+2):
        for rest in displacements(parts-1,budget-a):
            yield (a,)+rest


def check():
    output=Counter()
    samples=[]
    assert strong_cap(CENTER-1)==157 and surjections(30,3)**3>CENTER+1
    for rows,other,arc in seeds():
        n=len(rows)
        allowed=tuple(i for i in range(n) if i not in arc)
        cap=(157-n)//29
        tally=Counter()
        closest=None
        for size in range(len(allowed)+1):
            for ports in combinations(allowed,size):
                assignments=list(displacements(size,cap))
                if has_cycle(rows,ports):
                    tally['cyclic']+=len(assignments)
                    continue
                a,b=signature(rows,ports),signature(other,ports)
                assert value(a,(1,)*size)==27 and value(b,(1,)*size)==29
                for increments in assignments:
                    widths=tuple(1+29*t for t in increments)
                    lower,upper=value(a,widths),value(b,widths)
                    tally['evaluated']+=1
                    tally['unit in marked direction']+=upper-lower==2
                    tally['unit in reverse direction']+=lower-upper==2
                    tally['exact lower']+=CENTER-1 in (lower,upper)
                    tally['exact upper']+=CENTER+1 in (lower,upper)
                    tally['exact pair']+={lower,upper}=={CENTER-1,CENTER+1}
                    if abs(lower-upper)==2:
                        sample=(abs((lower+upper)//2-CENTER),(lower+upper)//2-CENTER,ports,widths)
                        if closest is None or sample<closest:
                            closest=sample
        assert tally['evaluated']+tally['cyclic']==comb(len(allowed)+cap,cap)
        output.update(tally)
        samples.append((rows,arc,dict(tally),closest))
    assert len(samples)==4
    assert (output['evaluated'], output['cyclic'], output['unit in marked direction']) == (714,42,28)
    assert not any(output[k] for k in ('unit in reverse direction','exact lower','exact upper','exact pair'))
    for row in samples: print(row)
    print('all types',dict(output))


if __name__=='__main__':
    check()
