"""Tree structure of the certified census of hostile points x = -p/3^j (output of ..._dim_z13scan.py).

For each census point x != -1 its PARENT is the census point h of smaller 3-adic exponent j_h < j_x that
agrees with x to the largest 2-adic depth i = v_2(x-h).  Reported:
  * agreement-depth deficit 1.585*j_x - i                    (hypothesis (iii) of the dimension note);
  * relative cost rho = (|x|-|h|) * R_min^h(i),  R_min^h(i) = 3^alpha_h(i)/2^i  (alpha_h exact, dim_alpha);
    rho = +1 means x = h - 2^i/3^alpha_h(i) (the canonical Theorem-P perturbation);
  * whether every child's scale i exceeds its parent's own scale;
  * the generation (depth) distribution of the resulting tree, overall and per j.
  * canonical steps x = h - 2^i/3^alpha_h(i) failing condition (A) of Theorem P, and (given the
    Theorem-P closure file) whether census points outside the closure have an all-canonical,
    (A)-satisfying ancestry.
usage: python3 ..._dim_census_tree.py CENSUS_FILE ALPHA_BIN [P0] [CLOSURE_FILE]
"""
import sys, subprocess
from fractions import Fraction as F
from collections import Counter

def j3(x):
    q=x.denominator; e=0
    while q%3==0: q//=3; e+=1
    return e
def v2(n):
    n=abs(n); return (n & -n).bit_length()-1

def main():
    cen=sorted([F(l.strip()) for l in open(sys.argv[1]) if l.strip()], key=lambda x: x.denominator)
    alpha_bin=sys.argv[2]; P0=sys.argv[3] if len(sys.argv)>3 else '26'
    rows=[]
    for x in cen:
        if x==-1: continue
        jx=j3(x); best=None
        for h in cen:
            if j3(h)>=jx: continue
            i=v2((x-h).numerator)
            if best is None or i>best[0]: best=(i,h)
        rows.append((x,jx,best[1],best[0]))
    need={}
    for x,jx,h,i in rows: need[h]=max(need.get(h,0),i)
    hs=list(need)
    out=subprocess.run([alpha_bin,P0,str(max(need.values()))]+[t for h in hs for t in (str(h.numerator),str(h.denominator))],
                       capture_output=True,text=True,check=True).stdout.strip().splitlines()
    al={h:[int(t) for t in line.split('alpha:')[1].split()] for h,line in zip(hs,out)}
    rho=Counter(); gap=Counter()
    for x,jx,h,i in rows:
        R=F(3**al[h][i-1],2**i)
        rho[(( -x)-(-h))*R]+=1
        gap[round(1.585*jx-i)]+=1
    par={x:h for x,jx,h,i in rows}; own={x:i for x,jx,h,i in rows}; own[F(-1)]=0
    nest=sum(1 for x,jx,h,i in rows if i<=own[h])
    def depth(x):
        d=0
        while x!=F(-1): x=par[x]; d+=1
        return d
    D=Counter(depth(x) for x in cen)
    byj={}
    for x in cen: byj.setdefault(j3(x),[]).append(depth(x))
    kids=Counter(par.values())
    print(f"census points: {len(cen)}; with parent: {len(rows)}")
    print("agreement-depth deficit 1.585 j - i:", dict(sorted(gap.items())))
    print("relative cost rho = (|x|-|h|) R_min^h(i):", {str(k):v for k,v in sorted(rho.items())})
    print("children whose scale does not exceed the parent's own scale:", nest)
    print("generation distribution:", dict(sorted(D.items())))
    print("max generation per j:", {j:max(v) for j,v in sorted(byj.items())})
    print(f"nodes with children: {len(kids)}; children of -1: {kids[F(-1)]}; max children: {max(kids.values())}")
    canon=[(x,jx,h,i) for x,jx,h,i in rows if x==h-F(2**i,3**al[h][i-1])]
    failA=[x for x,jx,h,i in canon if not (F(1,2)-F(1,2*3**jx) > (-h-1)+F(2**i,3**jx))]
    print(f"canonical steps x = h - 2^i/3^alpha_h(i): {len(canon)}; failing (A): {len(failA)} {[str(x) for x in failA]}")
    cs=set(t[0] for t in canon)
    nonc=Counter(jx for x,jx,h,i in rows if x not in cs)
    print("non-canonical steps by j:", dict(sorted(nonc.items())),
          "; largest j with a non-canonical step:", max(nonc) if nonc else None,
          "; census points above it (all canonical):", sum(1 for x,jx,h,i in rows if jx>max(nonc or [0])))
    if len(sys.argv)>4:
        clos={F(l.split()[0]) for l in open(sys.argv[4]) if l.strip()}
        good={x for x,jx,h,i in canon if x not in failA}
        def chain_ok(x):
            while x!=F(-1):
                if x not in good: return False
                x=par[x]
            return True
        out=[x for x in cen if x!=F(-1) and x not in clos]
        print(f"census points outside the Theorem-P closure: {len(out)}; with an all-canonical (A)-satisfying ancestry: {sum(chain_ok(x) for x in out)}")

if __name__=="__main__":
    main()
