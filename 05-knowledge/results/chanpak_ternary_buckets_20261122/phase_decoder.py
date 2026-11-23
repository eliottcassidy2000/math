"""Decode one of five physical source families and its centered digit from the odd target."""
from argparse import ArgumentParser

RAYS = ((21,-13),(33,7),(55,7),(280,7),(1019,219))
MODULUS = 27*7*19


def residues(modulus):
    phases=[]
    v=2%modulus
    seen=set()
    while v not in seen:
        phases.append(v)
        seen.add(v)
        v=2*v%modulus
    sources={(a*v+b)%modulus for a,b in RAYS for v in phases}
    children={(3*s+d)%modulus for s in sources for d in (-2,0,2)}
    return sources,children


def decoder():
    table={}
    for a,b in RAYS:
        for s in range(1,19):
            for digit in (-2,0,2):
                key=(3*(a*2**s+b)+digit)%MODULUS
                assert key not in table
                table[key]=(a,b,s,digit)
    return table


def read_target(target,table):
    """Return the exact arithmetic address, if it exists (physical domains in the note)."""
    candidate=table.get(target%MODULUS)
    if candidate is None:
        return None
    a,b,phase,digit=candidate
    q,rem=divmod(target-digit-3*b,3*a)
    if rem or q<2 or q&(q-1):
        return None
    s=q.bit_length()-1
    assert (s-phase)%18==0
    return a,b,s,digit


def main(target):
    table=decoder()
    sources,children=residues(MODULUS)
    assert (len(sources),len(children),len(table))==(90,270,270)
    assert set(table)==children and sources.isdisjoint(children)
    assert all(read_target(3*(a*2**s+b)+d,table)==(a,b,s,d)
               for a,b in RAYS for s in (1,2,13,18,19,61,1000) for d in (-2,0,2))
    for p in range(1,18,2):
        u=pow(2,p-1,MODULUS)
        mark=table.get((u*(2*u-1)-1)%MODULUS)
        assert mark==((280,7,13,0) if p==7 else None)
    print("parent residues",len(sources),"child residues",len(children),"period",18)
    if target is not None:
        print("arithmetic address",read_target(target,table))


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("target",nargs="?",type=int)
    main(parser.parse_args().target)
