"""Two transitive modules, exact run polynomials and marked golden observations."""
from boundary_search import *
from carrier import bin_values

BASE=(2,4,9,19,39,15)
VARIABLE=(0,5)
FIXED=(1,2,3,4)


def bivariate(kind='H',root=1):
    """H; E1 paths ending at root; S2 two bins with first starting at root."""
    position={v:i for i,v in enumerate(FIXED)};full=(1<<len(FIXED))-1
    separators=int(kind=='S2')
    @lru_cache(None)
    def walk(mask,last,used):
        out=Counter()
        if mask==full and used==separators and (kind!='E1' or last==root):out[0,0]=1
        if used<separators and last!=-1:
            out.update(walk(mask,-1,used+1))
        for v in range(6):
            if v in position and mask>>position[v]&1:continue
            if last!=-1 and not BASE[last]>>v&1:continue
            nextmask=mask|(1<<position[v]) if v in position else mask
            shift=(int(v==VARIABLE[0]),int(v==VARIABLE[1]))
            for (i,j),weight in walk(nextmask,v,used).items():out[i+shift[0],j+shift[1]]+=weight
        return out
    return walk(1<<position[root],root,0) if kind=='S2' else walk(0,-1,0)


def surjections(k,limit):
    v=[1]+[0]*limit
    for _ in range(k):v=[0]+[j*(v[j]+v[j-1]) for j in range(1,limit+1)]
    return v


def evaluate2(c,m,n):
    a=surjections(m,max(i for i,j in c));b=surjections(n,max(j for i,j in c))
    return sum(w*a[i]*b[j] for (i,j),w in c.items())


def source2(m,n):
    addresses=[(v,j) for v,c in enumerate((m,1,1,1,1,n)) for j in range(c)]
    T=tuple(sum(1<<j for j,(w,l) in enumerate(addresses) if
                (v==w and k<l) or (v!=w and BASE[v]>>w&1)) for v,k in addresses)
    ports=tuple(addresses.index((v,0)) for v in FIXED)
    return T,ports


def phase(m,n,root=1):
    h=evaluate2(bivariate('H'),m,n)
    e=evaluate2(bivariate('E1',root),m,n)
    s=evaluate2(bivariate('S2',root),m,n)
    require((s-e)%2==0,'observer parity')
    v=(h,(s-e)//2,(s+e)//2)
    return v,field(v)


def main():
    polys={k:bivariate(k) for k in ('H','E1','S2')}
    print('polynomials',{k:{str(ij):w for ij,w in c.items() if min(ij)>0} for k,c in polys.items()},flush=True)
    for m in range(1,5):
        for n in range(1,5):
            T,p=source2(m,n);h=tables(T)[0][-1]
            require(h==2*(2**m+1)*(2**n+1)-1==evaluate2(polys['H'],m,n),'bivariate Hamilton')
            P,E,S=bin_values(T,p[0],2)
            require(E[1]==evaluate2(polys['E1'],m,n) and S[2]==evaluate2(polys['S2'],m,n),'independent marked cuts')
            tables.cache_clear()
    rows=[dict(m=m,n=n,observer=phase(m,n)[0],phase=phase(m,n)[1]) for m,n in ((1,1),(3,1),(3,3),(2,2),(1,6),(4,3))]
    print(rows,flush=True)
    (HERE/'clone_clock.json').write_text(json.dumps(dict(polynomials={k:[[i,j,w] for (i,j),w in sorted(c.items())] for k,c in polys.items()},samples=rows),indent=2)+'\n')


if __name__=='__main__':main()
