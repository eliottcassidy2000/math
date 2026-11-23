"""Exact count changes along shortest physical edits between two growing constructions."""
from argparse import ArgumentParser
from itertools import combinations
from pathlib import Path
from subprocess import run

from physical_navigation import EXAMPLES, align, flipped, pad, sampled


def relabel(rows,verts):
    return tuple(sum(((rows[i]>>v)&1)<<j for j,v in enumerate(verts)) for i in verts)


def aligned(parent,child,s,perm):
    p,pl,pr=sampled(parent,s)
    c,cl,cr=sampled(child,s)
    L,R=max(len(pl),len(cl)),max(len(pr),len(cr))
    p=pad(p,pl,pr,L,R)
    c=pad(c,cl,cr,L,R)
    P=pl+pr
    C=cl+cr
    pt=[i for i in range(len(p)) if i not in P]
    ct=[i for i in range(len(c)) if i not in C]
    return relabel(p,P+pt),relabel(c,[C[j] for j in perm]+ct),R


def cube_heights(binary,parent,child,perm):
    values=[]
    edges=None
    for s in range(1,5):
        p,c,R=aligned(parent,child,s,perm)
        assert R<=3
        changes=[(i,j) for i,j in combinations(range(len(p)),2) if (p[i]^c[i])>>j&1]
        if edges is None: edges=changes
        assert edges==changes and all(j<len(p)-s for i,j in edges)
        corners=[p]
        for edge in edges:
            corners+= [flipped(rows,*edge) for rows in corners]
        assert corners[-1]==c
        data="".join(f"{len(rows)} "+' '.join(map(str,rows))+"\n" for rows in corners)
        output=run([binary,"--heights"],input=data,text=True,capture_output=True,check=True)
        values.append(list(map(int,output.stdout.split())))
    return edges,list(zip(*values))


def affine(vals):
    difference = vals[1] - vals[0]
    if difference % 2:
        return None
    a = difference // 2
    b = 2 * vals[0] - vals[1]
    return (a,b) if all((a<<s)+b==v for s,v in enumerate(vals,1)) else None


def less(p,q):
    return p[0]<=q[0] and 2*p[0]+p[1]<=2*q[0]+q[1]


def allowed_path(charts,increasing):
    best={0:[]}
    for m,z in enumerate(charts):
        if m==0 or z is None or not (less(charts[0],z) and less(z,charts[-1])):continue
        for i in range(m.bit_length()):
            j=m^(1<<i)
            if j<m and j in best and (not increasing or less(charts[j],z)):
                best[m]=best[j]+[i]
                break
    return best.get(len(charts)-1)


def main(binary):
    for parent,child,_,a,b in EXAMPLES:
        _,dist,perm=align(parent,child,1)
        edges,heights=cube_heights(binary,parent,child,perm)
        charts=[affine(v) for v in heights]
        assert charts[0]==(a,b) and charts[-1]==(3*a,3*b) and len(edges)==dist
        if a in (280, 1019):
            assert allowed_path(charts, False) is not None
        for increasing in (True,False):
            path=allowed_path(charts,increasing)
            print("source",a,b,"increasing",increasing,"path",path)
            if path is None: continue
            m=0
            for i in path:
                m^=1<<i
                print(edges[i],charts[m])
            break

if __name__=="__main__":
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("binary",type=Path)
    main(ap.parse_args().binary.resolve())
