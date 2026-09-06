#!/usr/bin/env python3
"""Independent universal LRC local-bound audit, no repository math imports.

Recomputes all 308 coefficient slopes by continuous knapsack on error slices,
then optionally compares a native interval TSV against a fresh raw integer box.
The H601 head is not rerun; its source and manifest are audited in the report.
"""
import argparse
import csv
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path
import sys

R = Q(3,14)
CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def cross(u,v):
    return (u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])


def speed_vertices(v):
    vertices = set()
    for free in range(3):
        if v[free] == 0:
            continue
        other = [j for j in range(3) if j != free]
        for values in product((Q(0),Q(1)),repeat=2):
            w = [Q(0)]*3
            for j,x in zip(other,values):
                w[j] = x
            w[free] = -sum(v[j]*w[j] for j in other)/v[free]
            if 0 <= w[free] <= 1:
                vertices.add(tuple(w))
    need(bool(vertices), ("nonempty normalized speed polytope",v))
    return sorted(vertices)


def fill_value(items, amount, reverse):
    total = Q(0)
    for ratio, capacity in sorted(items, reverse=reverse):
        take = min(capacity,amount)
        total += take*ratio
        amount -= take
    need(amount == 0, "continuous resource is filled exactly")
    return total


def slice_width(v,w,delta):
    """LP width via fractional knapsack; no error polygon vertices or clipping."""
    pivot = next(j for j in range(3) if v[j])
    j,k = (pivot+1)%3,(pivot+2)%3
    coefficients = [Q(0)]*3
    coefficients[j],coefficients[k] = -Q(w[k],v[pivot]),Q(w[j],v[pivot])
    items,free_width = [],Q(0)
    for p,ell in zip(v,coefficients):
        if p:
            # Change e to sign(p)e, then shift [-R,R] to [0,2R].
            # The constant objective shift cancels between max and min.
            slope = ell*(1 if p>0 else -1)/abs(p)
            items.append((slope,2*R*abs(p)))
        else:
            free_width += 2*R*abs(ell)
    amount = Q(delta)+R*sum(map(abs,v))
    need(0 <= amount <= 2*R*sum(map(abs,v)), "error slice resource range")
    return free_width+fill_value(items,amount,True)-fill_value(items,amount,False)


def pattern_slope(p):
    S = sum(p)
    defects = tuple(d for d in range(-3*S,3*S+1)
                    if 14*abs(d)<3*S and ((d%3 == 0) == all(x%3 for x in p)))
    unit = all(x%3 for x in p)
    rho = Q(2 if unit else 1,3)
    intercept = Q(4,3)*len(defects) if unit else Q(len(defects))
    maximum = Q(0)
    for negative in range(3):
        if p[negative] == 0:
            continue
        v = tuple(-x if j == negative else x for j,x in enumerate(p))
        for w in speed_vertices(v):
            value = rho*sum((slice_width(v,w,d) for d in defects),Q(0))
            maximum = max(maximum,value)
    return defects,maximum,intercept


def coefficient_box():
    # Full Cartesian cube, followed by symmetry quotient; no norm-20 list.
    patterns = sorted({tuple(sorted(p)) for p in product(range(19),repeat=3)
                       if sum(x!=0 for x in p)>=2 and sum(p)%2==0
                       and gcd(gcd(p[0],p[1]),p[2])==1
                       and sum(x%3==0 for x in p)<=1
                       and tuple(sorted(p)) not in ((0,1,1),(1,1,2))})
    need(len(patterns)==308, "complete coefficient box")
    need(Counter(sum(x!=0 for x in p) for p in patterns)=={2:15,3:293}, "support sizes")
    digest = sha256()
    equalities,max_defect = [],0
    for p in patterns:
        defects,slope,intercept = pattern_slope(p)
        need(slope<=Q(15,98), ("universal finite-box slope",p,slope))
        if slope==Q(15,98):
            equalities.append(p)
        max_defect=max(max_defect,max(map(abs,defects)))
        digest.update(repr((p,defects,slope,intercept)).encode())
    need(equalities==[(1,7,8)], "complete slope equality pattern")
    need(digest.hexdigest()=="3be81c2522a1df6a146e50634754620103f9d2d8840d17f34c9e9a4954e849f7",
         "all exact pattern/defect/slope/intercept rows match frozen semantic digest")
    need(max(sum(p) for p in patterns)==52 and max_defect==11, "actual norm and defect extremes")
    norm4=pattern_slope((1,1,2))[1]
    norm3=pattern_slope((1,1,1))[1]
    need(norm4==Q(2,7) and norm3==Q(2,7), "norm-four and missing-parity norm-three slope hostiles")
    return len(patterns),equalities,digest.hexdigest(),max_defect


def residue_fibers():
    cases=0
    for w in product((1,2),repeat=3):
        kernel=[v for v in product(range(3),repeat=3) if sum(a*b for a,b in zip(v,w))%3==0]
        for v in kernel:
            if v==(0,0,0):
                continue
            unit=all(v)
            need(sum(x==0 for x in v)==(0 if unit else 1), "primitive residue type dichotomy")
            for delta in range(3):
                fiber=[C for C in kernel if all((x-delta*y)%3==0 for x,y in zip(cross(v,C),w))]
                need(len(fiber)==3, "complete affine residue line")
                live=sum(all(C) for C in fiber)
                expected=(2 if delta==0 else 0) if unit else (0 if delta==0 else 1)
                need(live==expected, "exact owner word density and allowed defect")
                cases+=1
    return cases


def analytic_controls():
    need(Q(6,49)+Q(4,7*19)==Q(142,931)<Q(15,98), "large-M strict slope")
    need(Q(2,11)-Q(15,98)==Q(31,1078), "count-slope target gap")
    need(Q(308,31)*56+Q(4312,93)==Q(56056,93)<603, "S<=56 cutoff")
    g58=Q(3*58*58,16)-Q(308*58,31)-Q(4312,93)
    need(g58==Q(3023,372)>0 and Q(3*58,8)-Q(308,31)>0, "S>=58 quadratic cutoff")
    for S in range(2,201,2):
        limit=Q(3*S,14)
        integers=[d for d in range(-S,S+1) if abs(d)<limit]
        Du=[d for d in integers if d%3==0]
        Dz=[d for d in integers if d%3]
        need(Q(4*len(Du),3)<Q(2*S,7)+Q(4,3), "unit intercept bound")
        need(len(Dz)<Q(2*S,7)+Q(4,3), "one-zero intercept bound")
    # Polynomial difference has degree at most two in each of a,b. Nine
    # exact evaluations therefore certify the supplied identity identically.
    for a,b in product(map(Q,range(3)),repeat=2):
        left=8*(a*b+a+b)-3*(a+b)*(a+1)*(b+1)
        right=3*a*(1-a)*(b+1)+3*b*(1-b)*(a+1)+2*a*(1-b)+2*b*(1-a)
        need(left==right, "positive hexagon-area remainder identity")
    for w in ((1,5,11),(19,23,29),(1,137,277),(541,547,601)):
        a,b,c=map(Q,w)
        vertices=((c/(a+c),Q(0)),(Q(0),c/(b+c)),(-b/(a+b),a/(a+b)),
                  (-c/(a+c),Q(0)),(Q(0),-c/(b+c)),(b/(a+b),-a/(a+b)))
        area=abs(sum(p[0]*q[1]-p[1]*q[0] for p,q in zip(vertices,vertices[1:]+vertices[:1])))/2
        expected=2*c*(a*b+a*c+b*c)/((a+b)*(a+c)*(b+c))
        need(area==expected and area>Q(3,4), "direct six-vertex area and strict gain")
    # Actual-zero endpoint section is nonzero in the closed cube. The
    # strict-defect sampling convention sets it to zero without changing area.
    v,w=(13,-1,0),(1,13,19)
    need(slice_width(v,w,3)==Q(3,7), "support-two endpoint jump hostile")


def raw(w):
    a,b,c=w
    roofs=[3*(sum(w)-wi) for wi in w]
    found=[]
    for y in range(-roofs[1]//14-1,roofs[1]//14+2):
        if not y%3 or 14*abs(y)>=roofs[1]:
            continue
        for z in range(-roofs[2]//14-1,roofs[2]//14+2):
            if not z%3 or 14*abs(z)>=roofs[2] or (-b*y-c*z)%a:
                continue
            x=(-b*y-c*z)//a
            if x%3 and 14*abs(x)<roofs[0]:
                found.append((x,y,z))
    D=42*a*b*c
    q=18*a*b
    E=[0,0,0]
    mass=0
    for C in found:
        terms=[min(q,3*w[i]*(roofs[i]-14*abs(C[i]))) for i in range(3)]
        need(min(terms)>0, "strict raw projection summands")
        for i in range(3):
            E[i]+=terms[i]
        mass+=min(terms)
    return D,tuple(E),mass,len(found)


def native_comparison(path):
    expected={w for w in combinations([n for n in range(1,98) if n%2 and n%3],3)
              if gcd(gcd(w[0],w[1]),w[2])==1}
    observed=set()
    equalities=[]
    norm4=0
    with Path(path).open(newline="",encoding="utf-8") as f:
        for line in csv.DictReader(f,delimiter="\t"):
            row={k:int(v) for k,v in line.items()}
            w=tuple(row[k] for k in ("a","b","c"))
            need(w in expected and w not in observed, "native fresh-head row coverage")
            observed.add(w)
            D,E,mass,N=raw(w)
            need(D==row["denominator"], "common 42abc denominator")
            need(E==tuple(row["E"+str(i)+"_numerator"] for i in range(3)), "all three native/raw projection values")
            need(mass==row["mass_numerator"] and 3*N==row["contacts"], "physical mass and three inverse branches")
            need(77*min(E)<=6*D, "fresh-head selected bound")
            if 77*min(E)==6*D:
                equalities.append(w)
            a,b,c=w
            has_norm4=c in (2*a+b,a+2*b) or 2*b==a+c
            norm4+=has_norm4
            if not has_norm4:
                need(77*max(E)<6*D, "fresh-head strict bound outside norm four")
    need(observed==expected and len(observed)==5409, "complete independently generated H97 universe")
    need(equalities==[(1,5,11)] and norm4==241, "fresh-head equality and exact norm-four partition")
    return len(observed),norm4,equalities


def full_head_universe():
    speeds=[n for n in range(1,602) if n%2 and n%3]
    rows=norm4=0
    for a,b,c in combinations(speeds,3):
        if gcd(gcd(a,b),c)!=1:
            continue
        rows+=1
        norm4+=c in (2*a+b,a+2*b) or 2*b==a+c
    need((rows,norm4)==(1317935,9201), "independent H601 universe and norm-four identity counts")
    return rows,norm4


def additive_parity_controls():
    rows=[]
    for w in ((2,5,7),(1,7,8),(1,10,11),(1,61,62),(1,121,122)):
        a,b,c=w
        D,E,mass,N=raw(w)
        exact=Q(mass,D)
        integral=2*(Q(3,7*c)*Q(3*a,14)
                    +(Q(3,7*c)+Q(3,14*b))*Q(3*(b-a),28)
                    +Q(3,14*b)*Q(3*a,28))
        need(integral==Q(27,196), "norm-three physical tent area")
        need(abs(exact-Q(9,98))<=Q(6,7*c), "norm-three lattice quadrature control")
        need(exact>Q(6,77), "mixed-parity additive hostile")
        rows.append((w,str(exact),N))
    need(Q(9,98)-Q(6,77)==Q(15,1078) and Q(6,7*62)<Q(15,1078),
         "cofinite additive physical obstruction cutoff")
    D,E,mass,N=raw((2,11,20))
    need(Q(mass,D)==Q(11,140)>Q(6,77) and N==4,
         "mixed-parity norm-four nonadditive hostile")
    return rows


def main():
    sys.stdout.reconfigure(newline="\n")
    parser=argparse.ArgumentParser()
    parser.add_argument("--native-tsv")
    args=parser.parse_args()
    box=coefficient_box()
    fibers=residue_fibers()
    analytic_controls()
    full_head=full_head_universe()
    native=native_comparison(args.native_tsv) if args.native_tsv else None
    additive=additive_parity_controls()
    print("status=INDEPENDENT_EXACT_AUDIT; analytic proof and H601 evidence scope in companion report")
    print("coefficient_box=(patterns,equality,semantic_sha256,max_abs_defect):"+str(box))
    print("width_method=continuous_knapsack; no_error_vertex_or_polygon_clipping_imports")
    print("coefficient_supports=293_full+15_support_two; actual_max_norm=52")
    print("finite_field_owner_fibers="+str(fibers))
    print("analytic_controls=zonotope_scale,intercepts,hexagon_area,strict_56_58_cutoff,support_two_endpoint")
    print("fresh_H97_native_vs_raw=(rows,norm4,equalities):"+str(native))
    print("H601_universe_only=(rows,norm4):"+str(full_head))
    print("mixed_parity_additive_controls="+str(additive))
    print("mixed_parity_additive_family=mu>=9/98-6/(7c)>6/77_for_c>=62")
    print("additive_only_extension=FALSE; nonadditive_(2,11,20)_physical_mass=11/140")
    print("H601_full_census=SOURCE_AND_TRANSCRIPT_AUDITED_NOT_RERUN")
    print("checks="+str(CHECKS)+"; RESULT=PASS; LRC14=OPEN")


if __name__=="__main__":
    main()
