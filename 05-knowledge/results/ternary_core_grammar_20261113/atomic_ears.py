"""Root-preserving near-transitive ear certificates on a complete small bank."""
from itertools import combinations
from functools import lru_cache
from triangle_probe import bank,need,save,HERE
from peeling import induced,strong,near

@lru_cache(None)
def parse(T,r):
    atom=near(T)
    if atom is not None:return dict(kind='atom',rows=T,root=r,reversed_edge=atom)
    n=len(T);others=[v for v in range(n) if v!=r]
    for count in range(2,n-2):
        for removed in combinations(others,count):
            keep=[v for v in range(n) if v not in removed];U=induced(T,keep)
            if not strong(U):continue
            for v in keep:
                av=(v,)+removed;A=induced(T,av)
                if not strong(A) or near(A) is None:continue
                parent=parse(U,keep.index(r))
                if parent is not None:return dict(kind='ear',rows=T,root=r,removed=removed,shared=v,atom_rows=A,keep=keep,parent=parent)
    return None

def pure_ends(T):
    n=len(T);D=[[0]*n for _ in range(1<<n)]
    for v in range(n):D[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if D[mask][v]:
                for w in range(n):
                    if not(mask>>w&1) and T[v]>>w&1:D[mask|1<<w][w]+=D[mask][v]
    return D

def profile(T,r):
    n=len(T);end=pure_ends(T);H=[sum(row) for row in end];full=(1<<n)-1
    rev=tuple(full^(1<<v)^T[v] for v in range(n));s=sum(pure_ends(rev)[full][r:r+1])
    F=sum(end[S][r]*H[full^S] for S in range(1,1<<n) if S>>r&1)
    # Empty complement is one path choice.
    F+=end[full][r]
    q=H[full^(1<<r)] if n>1 else 1
    return (H[full],s,q,F)

def check_tree(tree):
    T=tuple(tree['rows']);r=tree['root'];need(strong(T),'strong entry at every tree node')
    h=sum(pure_ends(T)[-1])
    if tree['kind']=='atom':
        need(h==(1<<(len(T)-2))+1,'near-transitive atom count');return (1,0)
    U=tuple(tree['parent']['rows']);A=tuple(tree['atom_rows']);ha=sum(pure_ends(A)[-1]);hu=sum(pure_ends(U)[-1])
    need(h>=ha*hu and ha>=3 and hu<h,'target-decreasing actual entry')
    need(tree['keep'][tree['parent']['root']]==r,'same physical root retained')
    leaves,depth=check_tree(tree['parent']);return leaves+1,depth+1

def run():
    records=[];hist={};twin25=[]
    for rec in bank(7):
        T=tuple(rec['rows']);r=len(T)-1;tree=parse(T,r)
        need(tree is not None,'rooted atomic ear exists')
        leaves,depth=check_tree(tree)
        hist[len(T)]=hist.get(len(T),0)+1
        records.append(dict(n=len(T),ordinary_code=rec['ordinary_code'],rows=T,root=r,leaves=leaves,depth=depth,tree=tree))
        z=profile(T,r)
        if z[0]==25:twin25.append(dict(n=len(T),rows=T,code=rec['ordinary_code'],profile=z,readiness=z[3]-z[0]))
    save('atomic_ears',dict(scope='Every rooted strong class through7; literal roots, cross edges and atom order retained; independent positive integer path count verifies each descent',levels=hist,records=records,H25=twin25))
    print('ROOTED LEVELS',hist,'TOTAL',len(records))
    print('H25 PROFILES',twin25)

if __name__=='__main__':run()
