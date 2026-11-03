"""All-level affine clocks on a recursive source/triangle tournament family."""
from itertools import combinations
from support import need,save
def first_unit(n):
    while n%3==0:n//=3
    return n%3
def build(k,flipped=False):
    N=3**k;B=N//3;n=N+B;rows=[0]*n
    block=lambda i:i%B if i<N else i-N
    for i,j in combinations(range(n),2):
        b,c=block(i),block(j)
        if b!=c:arc=first_unit(c-b)==1
        else:
            arc=(i>=N) if i>=N or j>=N else ((j//B-i//B)%3==1)
            if flipped and b==0:arc=not arc
        rows[i if arc else j]|=1<<(j if arc else i)
    return tuple(rows)
def affine(k,a,b):
    N=3**k;B=N//3
    return tuple((a*i+b)%N for i in range(N))+tuple(N+(a*i+b)%B for i in range(B))
def triangles(T):return {t for t in combinations(range(len(T)),3) if all(sum(T[u]>>v&1 for v in t)==1 for u in t)}
def auto(T,p):return all(((T[i]>>j)&1)==((T[p[i]]>>p[j])&1) for i,j in combinations(range(len(T)),2))
def main():
    checks=[]
    for k in range(1,5):
        N=3**k;B=N//3;A=build(k);T=build(k,True);C=triangles(A)
        need(C==triangles(T),'all-level identical labelled cyclic triples')
        for p in (affine(k,1,1),affine(k,2,0)):
            need({tuple(sorted(p[i] for i in t)) for t in C}==C,'translation/doubling preserve actual triples')
        counts=None
        if k<=3:
            GA=GB=0
            for a in range(N):
                if a%3==0:continue
                for b in range(N):
                    p=affine(k,a,b);aa=auto(A,p);bb=auto(T,p)
                    need(aa==(a%3==1),'exact oriented uniform affine subgroup')
                    need(bb==(a%3==1 and b%B==0),'exact one-flip affine subgroup')
                    GA+=aa;GB+=bb
            need((GA,GB)==(3**(2*k-1),3**k),'all affine intersection counts')
            counts=[GA,GB]
        checks.append({'k':k,'vertices':len(A),'inner_clock':N,'cyclic_triples':len(C),
            'affine_group':2*3**(2*k-1),'checked_affine_intersections':counts,
            'proved_full_aut_orders':[3**((N-1)//2),3**((N-1)//2-(k-1))],
            'proved_full_aut_ratio':B})
    save('ternary_family.json',{'status':'FINITE_EXACT_GENERATOR_CONTROLS_FOR_PROVED_ALL_LEVEL_FAMILY','levels':checks,
        'scope':'full automorphism orders use the stated unique maximal-module recursion; enumerated affine subgroup checks do not count all automorphisms at high levels'})
    print('PASS same-triple families through108 vertices, affine intersections through36; all-size Aut recurrence in proof',flush=True)
if __name__=='__main__':main()
