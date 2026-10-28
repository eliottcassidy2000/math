"""Frobenius and gain transport on actual Paley-seven strong-atom attachments."""
from itertools import product
from pathlib import Path
import json
import kernel125 as k125

CHECKS=0
def need(b,msg):
    global CHECKS
    CHECKS+=1
    if not b: raise RuntimeError(msg)

def mul8(a,b):
    z=0
    while b:
        if b&1:z^=a
        a<<=1;b>>=1
        if a&8:a^=11
    return z
def pow8(a,k):
    z=1
    for _ in range(k):z=mul8(z,a)
    return z
def trace8(a):return a^pow8(a,2)^pow8(a,4)
def profile(x,y):return tuple(trace8(mul8(pow8(x,1<<j),y)) for j in range(3))
G=list(product(range(7),range(3)))
def compose(g,h):
    b,a=g;c,d=h
    return ((b+pow(2,a,7)*c)%7,(a+d)%3)
def inverse(g):
    b,a=g
    return ((-pow(2,(-a)%3,7)*b)%7,(-a)%3)
def perm(g):
    b,a=g
    return tuple((pow(2,a,7)*u+b)%7 for u in range(7))
def conv(g):return (-g[0]%7,g[1])
def difference(g,h):return compose(h,inverse(g))
P=tuple(sum(1<<((i+d)%7) for d in (1,2,4)) for i in range(7))

def attached(g):
    h=perm(g)
    rows=[P[i]|(1<<(7+h[i])) for i in range(7)]
    rows+=[(P[j]<<7)|sum(1<<i for i in range(7) if h[i]!=j) for j in range(7)]
    return tuple(rows)

def main():
    need(k125.strong(P),'Paley source strong')
    images=set()
    for g in G:
        p=perm(g);images.add(p)
        need(all((P[i]>>j&1)==(P[p[i]]>>p[j]&1) for i in range(7) for j in range(7)),'actual Paley automorphism')
        need(compose(g,inverse(g))==(0,0),'inverse')
        need(k125.strong(attached(g)),'actual matching attachment strong')
        for u in range(7):
            lhs=pow8(2,p[u]); rhs=mul8(pow8(2,g[0]),pow8(pow8(2,u),1<<g[1]))
            need(lhs==rhs,'semilinear field intertwining')
        for h in G:
            need(perm(compose(g,h))==tuple(perm(g)[v] for v in perm(h)),'group action')
            need(conv(compose(g,h))==compose(conv(g),conv(h)),'converse normalizes gains')
    # Exhaust all permutations to establish that these are *all* chart actions.
    from itertools import permutations
    exact={p for p in permutations(range(7)) if all((P[i]>>j&1)==(P[p[i]]>>p[j]&1) for i in range(7) for j in range(7))}
    need(exact==images,'full automorphism group21')
    need([g for g in G if perm(g)[0]==0]==[(0,0),(0,1),(0,2)],'owner pin leaves Frobenius3')
    orbit=[pow8(2,1<<j) for j in range(3)]
    need(orbit[0]^orbit[1]^orbit[2]==0,'actual Frobenius orbit is nonnormal')
    words=[profile(2,y) for y in orbit]
    need(len(set(words))==3,'restricted orbit profile separates phases')
    need(profile(2,2)==profile(2,3),'unrestricted field collision')
    for a,b,c in product(range(3),repeat=3):
        delta=(b-a)%3; epsilon=(c-b)%3
        need(profile(pow8(2,1<<a),pow8(2,1<<c))==profile(2,pow8(2,1<<((delta+epsilon)%3))),'actual pinned Frobenius LCA')
    # Rooted tree 0-1,0-2,1-3, fixed root chart identity. All 21^3 states.
    pair_checks=0
    for g1,g2,g3 in product(G,repeat=3):
        charts=[(0,0),g1,g2,g3]
        tree=[g1,g2,difference(g1,g3)]
        recovered=[(0,0),tree[0],tree[1],compose(tree[2],tree[0])]
        need(recovered==charts,'semilinear tree reconstruction')
        for i in range(4):
            for j in range(i+1,4):
                need(perm(difference(charts[i],charts[j]))==perm(difference(recovered[i],recovered[j])),'actual off-tree matching')
                pair_checks+=1
    complement=k125.native_complement(P,tuple(range(7)))
    need(sorted(row.bit_count() for row in complement)!=[3]*7,'native complement exits Paley family')
    out={'status':'FINITE-EXACT','actual_strong_atom_order':7,'full_chart_group':21,'pinned_owner_chart_group':3,'matching_attachments':21,'all_attachments_strong_order':14,'F8_polynomial_mask':11,'actual_frobenius_orbit':orbit,'restricted_profiles':words,'unrestricted_profile_collision':[2,3],'four_atom_relative_gain_states':21**3,'off_tree_matching_checks':pair_checks,'native_complement_scores':sorted(row.bit_count() for row in complement),'checks':CHECKS}
    Path(__file__).with_name('strong_ports.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out))

if __name__=='__main__':main()
