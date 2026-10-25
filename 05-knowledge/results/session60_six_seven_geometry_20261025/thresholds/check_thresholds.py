#!/usr/bin/env python3
"""Exact bounded controls for the root-cover threshold and its singularity guard.

Only Python's standard library is used. These computations check formulas and
explicit matrix instances; they do not replace the proofs in ROOT_THRESHOLDS.md.
"""
from fractions import Fraction
from itertools import combinations
from math import comb, gcd
from functools import reduce
import json
from pathlib import Path


def check(condition):
    """A validation gate that stays active under python -O."""
    if not condition:
        raise AssertionError("exact validation gate failed")


def canonical(d, n, r=2):
    return (r - 1) * n - r * (d + 1)


def hilbert_koszul(d, n, r, m):
    if m < 0:
        return 0
    c = n - d - 1
    return sum((-1)**j * comb(c,j) * comb(n-1+m-r*j,n-1)
               for j in range(min(c,m//r)+1))


def hilbert_free_basis(d, n, r, m):
    if m < 0:
        return 0
    p = [1]
    for _ in range(n-d-1):
        q = [0] * (len(p)+r-1)
        for i, coefficient in enumerate(p):
            for e in range(r):
                q[i+e] += coefficient
        p = q
    return sum(coefficient * comb(d+m-e,d)
               for e,coefficient in enumerate(p) if e <= m)


def top_chern(d,n,r):
    c = n-d-1
    return sum(comb(n,d-j) * (-r)**j * (comb(c+j-1,j) if c else int(j==0))
               for j in range(d+1)) * r**c


def rank(matrix):
    if not matrix or not matrix[0]:
        return 0
    a = [list(map(Fraction,row)) for row in matrix]
    row = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(row,len(a)) if a[i][col]),None)
        if pivot is None:
            continue
        a[row],a[pivot] = a[pivot],a[row]
        p = a[row][col]
        a[row] = [x/p for x in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                p = a[i][col]
                a[i] = [x-p*y for x,y in zip(a[i],a[row])]
        row += 1
        if row == len(a):
            break
    return row


def invert(matrix):
    n = len(matrix)
    a = [[Fraction(x) for x in matrix[i]]+[Fraction(i==j) for j in range(n)]
         for i in range(n)]
    for col in range(n):
        pivot = next(i for i in range(col,n) if a[i][col])
        a[col],a[pivot] = a[pivot],a[col]
        p = a[col][col]
        a[col] = [x/p for x in a[col]]
        for i in range(n):
            if i != col:
                p = a[i][col]
                a[i] = [x-p*y for x,y in zip(a[i],a[col])]
    return [row[n:] for row in a]


def relation_rows(forms):
    d1 = len(forms[0])
    inv = invert(forms[:d1])
    result = []
    for j in range(d1,len(forms)):
        # L_j = coefficient row times matrix of first d+1 forms.
        coefficients = [sum(forms[j][k]*inv[k][i] for k in range(d1))
                        for i in range(d1)]
        row = [-x for x in coefficients]+[Fraction(0)]*(len(forms)-d1)
        row[j] = 1
        result.append(row)
    return result


def smooth_support_check(d,n):
    # L_t = y_0+t*y_1+...+t^d*y_d, t=0,...,n-1.
    forms = [[t**i for i in range(d+1)] for t in range(n)]
    minors = 0
    for indices in combinations(range(n),d+1):
        check(rank([forms[i] for i in indices]) == d+1)
        minors += 1
    relations = relation_rows(forms)
    c = n-d-1
    supports = 0
    for number_zero in range(d+1):
        for zero in combinations(range(n),number_zero):
            nonzero = [j for j in range(n) if j not in zero]
            check(rank([[row[j] for j in nonzero] for row in relations]) == c)
            supports += 1
    # Restrict to L_0=0: remove zero first column and divide row L_t by t.
    restricted = [[Fraction(forms[t][i],t) for i in range(1,d+1)]
                  for t in range(1,n)]
    for indices in combinations(range(n-1),d):
        check(rank([restricted[i] for i in indices]) == d)
    return {"d":d,"n":n,"full_size_nonzero_minors":minors,
            "smooth_jacobian_supports":supports,
            "coordinate_restriction_general_position":True}


def row(d,n,r=2):
    a = canonical(d,n,r)
    degree = r**(n-d-1)
    answer = {"d":d,"n":n,"r":r,"codimension":n-d-1,
              "canonical_coefficient":a,"degree":degree,
              "cover_degree":r**(n-1),
              "geometric_genus":hilbert_koszul(d,n,r,a),
              "canonical_top_intersection":a**d*degree,
              "euler_characteristic":top_chern(d,n,r)}
    if d == 1:
        g2 = 1 + Fraction(a*degree,2)
        check(g2.denominator == 1)
        check(answer["geometric_genus"] == g2)
        answer["curve_genus"] = int(g2)
    return answer



def cross(a,b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def primitive(vector):
    divisor = reduce(gcd,(abs(x) for x in vector))
    check(divisor)
    v = tuple(x//divisor for x in vector)
    return v if next(x for x in v if x) > 0 else tuple(-x for x in v)


def incidences(lines):
    points = {primitive(cross(a,b)) for a,b in combinations(lines,2)}
    counts = {}
    for point in sorted(points):
        m = sum(sum(x*y for x,y in zip(line,point)) == 0 for line in lines)
        counts[m] = counts.get(m,0)+1
    check(sum(comb(m,2)*t for m,t in counts.items()) == comb(len(lines),2))
    return counts


def arrangement_control(name,lines):
    n = len(lines)
    counts = incidences(lines)
    complement = 3-2*n+sum((m-1)*t for m,t in counts.items())
    open_lines = 2*n-sum(m*t for m,t in counts.items())
    e_strata = 2**(n-1)*complement + 2**(n-2)*open_lines
    for m,t in counts.items():
        if m == 2:
            e_strata += 2**(n-3)*t
        else:
            g = 1+(m-4)*2**(m-3)
            e_strata += t*2**(n-m-1)*(2-2*g)
    e_closed = 2**(n-3)*(12-4*n+counts.get(2,0)+sum(m*t for m,t in counts.items() if m>=3))
    check(e_strata == e_closed)
    k2 = 2**(n-3)*((n-6)**2-sum(t*(m-3)**2 for m,t in counts.items() if m>=4))
    smooth = row(2,n)
    check(e_closed == smooth['euler_characteristic']-2**(n-4)*sum(t*m*(m-3) for m,t in counts.items() if m>=4))
    chi = Fraction(e_closed+k2,12)
    check(chi.denominator == 1)
    chi_corrected = Fraction(1+smooth['geometric_genus'])-Fraction(2)**(n-6)*sum(t*(m-3)*(m-2) for m,t in counts.items() if m>=4)
    check(chi == chi_corrected)
    positive_genus_quotients = [1+(m-4)*2**(m-3) for m in counts if m>=4]
    return {'name':name,'n':n,'intersection_multiplicities':counts,'resolved_K_squared':k2,
            'resolved_euler':e_closed,'resolved_chi_O':int(chi),
            'positive_genus_concurrent_quotients':positive_genus_quotients,
            'obstruction_mechanism':'positive-genus curve quotient' if positive_genus_quotients else ('positive geometric genus' if n>=6 else 'not asserted')}


def main():
    checks = 0
    for d in range(1,7):
        for r in range(2,7):
            for n in range(d+1,2*d+6):
                a = canonical(d,n,r)
                for m in sorted(set([0,1,r-1,r,r+1,max(0,a),max(0,2*a)])):
                    check(hilbert_koszul(d,n,r,m) == hilbert_free_basis(d,n,r,m))
                    checks += 1
                if a == 0:
                    check(hilbert_koszul(d,n,r,a) == 1)
                if 0 <= a < r:
                    check(hilbert_koszul(d,n,r,a) == comb(n-1+a,n-1))
    threshold_rows = []
    for d in range(1,7):
        for n in (2*d+1,2*d+2,2*d+3):
            threshold_rows.append(row(d,n))
        cy,gt = threshold_rows[-2:]
        check(cy["canonical_coefficient"] == 0)
        check(cy["degree"] == 2**(d+1))
        check(cy["geometric_genus"] == 1)
        check(gt["canonical_coefficient"] == 1)
        check(gt["geometric_genus"] == 2*d+3)
        check(gt["canonical_top_intersection"] == 2**(d+2))
    surface_controls = [row(2,n) for n in range(4,10)]
    for r0 in surface_controls:
        n = r0["n"]
        check(r0["euler_characteristic"] == 2**(n-4)*(n*n-9*n+24))
        check(r0["euler_characteristic"] == 12*(1+r0["geometric_genus"])-r0["canonical_top_intersection"])
    cy_equal_degree = [row(d,n,r) for d in range(1,7) for r in range(2,d+3)
                       for n in range(d+2,2*d+3) if canonical(d,n,r)==0]
    check([(x["r"],x["n"]) for x in cy_equal_degree if x["d"]==2] == [(2,6),(4,4)])
    check([(x["r"],x["n"]) for x in cy_equal_degree if x["d"]==3] == [(2,8),(3,6),(5,5)])
    local_cones = []
    for r in range(2,6):
        for m in range(3,9):
            degree = r**(m-2)
            genus = 1+Fraction((r*(m-2)-m)*degree,2)
            check(genus.denominator == 1)
            check((genus == 0) == ((r,m) == (2,3)))
            check(genus == hilbert_koszul(1,m,r,canonical(1,m,r)))
            discrepancy = m-1-r*(m-2)
            check((discrepancy+1)*degree == -(2*genus-2))
            local_cones.append({"r":r,"concurrent_lines":m,"exceptional_degree":degree,
                                "exceptional_genus":int(genus),"discrepancy":discrepancy})
    hostile_cones = []
    for n in range(5,11):
        m = n-1
        degree = 2**(n-3)
        genus = 1+(n-5)*2**(n-4)
        discrepancy = 4-n
        kres2 = ((n-6)**2-discrepancy**2)*degree
        check(kres2 == 8*(1-genus))
        hostile_cones.append({"n":n,"concurrent_lines":m,"dualizing_coefficient":n-6,
                              "base_curve_genus":genus,"resolution_geometric_genus":0,
                              "resolution_irregularity":genus,"resolution_Kodaira_dimension":"-infinity",
                              "discrepancy":discrepancy,"resolution_K_squared":kres2})
    restriction_ladders = {
      "three_quadrics_genus_five":[row(d,d+4) for d in (3,2,1)],
      "four_quadrics_degree_sixteen":[row(d,d+5) for d in (3,2,1)]}
    for ladder in restriction_ladders.values():
        for before,after in zip(ladder,ladder[1:]):
            check(after["canonical_coefficient"] == before["canonical_coefficient"]+1)
            check(after["degree"] == before["degree"])
    arrangements = [
        arrangement_control('seven GP lines', [[1,t,t*t] for t in range(7)]),
        arrangement_control('magic nine', [[1,0,0],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[1,1,1],[1,-1,-1],[1,1,-1],[1,-1,1]]),
        arrangement_control('hourglass seven', [[1,0,0],[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[1,1,1],[1,-1,-1]]),
        arrangement_control('four concurrent plus two transverse', [[0,1,t] for t in range(4)]+[[1,0,0],[1,1,7]]),
        arrangement_control('four concurrent plus one transverse', [[0,1,t] for t in range(4)]+[[1,0,0]]),
        arrangement_control('six concurrent plus one transverse', [[0,1,t] for t in range(6)]+[[1,0,0]])]
    canonical_surfaces = [row(2,n,r) for r in range(2,7) for n in range(4,10) if canonical(2,n,r)==1]
    check([(x['r'],x['n'],x['geometric_genus'],x['canonical_top_intersection']) for x in canonical_surfaces] == [(2,7,7,16),(3,5,5,9),(5,4,4,5)])
    result = {"status":"FINITE-EXACT controls, not a replacement for the proof",
              "hilbert_formula_comparisons":checks,
              "smoothness_support_controls":[smooth_support_check(d,2*d+3) for d in range(1,5)],
              "square_root_threshold_rows":threshold_rows,
              "surface_controls":surface_controls,
              "canonical_surface_classification":canonical_surfaces,
              "all_incidence_arrangement_controls":arrangements,
              "equal_degree_Calabi_Yau_rows_d1_through_d6":cy_equal_degree,
              "coordinate_section_ladders":restriction_ladders,
              "concurrent_local_cones":local_cones,
              "hostile_projective_cones":hostile_cones}
    out = Path(__file__).with_name("threshold_results.json")
    out.write_text(json.dumps(result,indent=2)+"\n")
    print(f"PASS: {checks} independent Hilbert comparisons; general-position matrix controls d=1..4; thresholds d=1..6; cone guards. Wrote {out.name}")

if __name__ == "__main__":
    main()
