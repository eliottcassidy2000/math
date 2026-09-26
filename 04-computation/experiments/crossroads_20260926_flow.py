"""Exact bad-ancestor incidence energy; stdlib, no Collatz oracle.

Reproduce: python3 04-computation/experiments/crossroads_20260926_flow.py
Each cell retains the exact endpoint congruence and affine slope.
"""
from fractions import Fraction
from math import gcd
import json
from collections import Counter


def shortcut(n, q=3, sign=1):
    return (q*n+sign)//2 if n & 1 else n//2


def cells(depth, q=3, sign=1, cap=None):
    """All prefixes k<depth of slope-undecided parity classes mod 2^depth."""
    rows = []
    bad = []
    for source in range(1 << depth):
        v, e, history, admissible = source, 0, [], True
        for k in range(depth):
            history.append((k, e, v))
            e += v & 1
            v = shortcut(v, q, sign)
            if q**e <= 2**(k+1) or (cap is not None and q**e > cap*2**(k+1)):
                admissible = False
                break
        if admissible:
            bad.append(source)
            for k, e, v in history:
                modulus = 2**(depth-k) * q**e
                rows.append((k, e, v % modulus, modulus))
    return bad, rows


def energy(depth, q=3, sign=1, cap=None):
    bad, rows = cells(depth, q, sign, cap)
    numerator = 0
    # All pair-intersection contributions have common denominator 2^(2L).
    # Keep duplicate occurrences: capacities count each possible hit time.
    for i, (k, e, residue, modulus) in enumerate(rows):
        for j in range(i + 1):
            ell, f, other, other_modulus = rows[j]
            if (residue-other) % gcd(modulus, other_modulus) == 0:
                term = q**min(e, f) * 2**(depth-max(k, ell))
                numerator += term if i == j else 2*term
    rho = Fraction(len(bad), 2**depth)
    second = Fraction(numerator, 2**(2*depth))
    mean = depth*rho
    return {"L": depth, "bad_classes": len(bad), "cells": len(rows),
            "rho": str(rho), "mean": str(mean), "second": str(second),
            "bound": str(rho*rho/second) if second else "0",
            "energy_over_L2_rho": float(second/(depth*depth*rho)) if rho else 0}, rows


def direct_audit(depth, q=3, sign=1):
    data, rows = energy(depth, q, sign)
    period = 2**depth * q**(depth-1)
    total = total2 = Fraction(0)
    for v in range(period):
        value = sum((Fraction(q**e, 2**k) for k,e,r,m in rows if v % m == r), Fraction(0))
        total += value
        total2 += value*value
    assert total/period == Fraction(data["mean"])
    assert total2/period == Fraction(data["second"])
    return period


def cohomology_probes():
    """Periodic residue rank obstruction, real cycles, and affine guard sidecar."""
    for power in range(2, 13):
        n = 2**power-1
        x = n
        for j in range(power-1):
            assert (3*x+1) % 4 == 2
            x = (3*x+1)//2
        assert x == 2*3**(power-1)-1
    # Word (1) closes every odd residue quotient along -1, but not in N.
    # Exact affine fixed point is -1 for plus, +1 for minus, -1/3 for q=5.
    centers = {"3n+1": Fraction(-1), "3n-1": Fraction(1), "5n+1": Fraction(-1,3)}
    for q, sign, cycle, exponents in [(3,-1,[5,7],[1,2]), (5,1,[13,33,83],[1,1,5])]:
        A,C,B = 1,0,1
        for k in exponents:
            A,C,B = q*A, q*C+sign*B, B*2**k
        assert Fraction(-C,A-B) == cycle[0]
        for x,y in zip(cycle, cycle[1:]+cycle[:1]):
            z=q*x+sign
            while z%2 == 0: z//=2
            assert z == y
    return {k:str(v) for k,v in centers.items()}


def incidence_harmonic(depth, q=3, sign=1):
    """Exact joint source/hub harmonic dual on its full CRT period."""
    bad, rows = cells(depth, q, sign)
    period = 2**depth*q**(depth-1)
    source_period = period*2**(depth-1)
    denominator = 2**(depth-1)
    capacities = [0]*period
    for k,e,residue,modulus in rows:
        weight = q**e*2**(depth-1-k)
        for v in range(residue,period,modulus):
            capacities[v] += weight
    histogram = Counter()
    sum_along_sources=0
    for residue in bad:
        for n in range(residue,source_period,2**depth):
            v, maximum = n, 0
            for _ in range(depth):
                maximum = max(maximum,capacities[v % period])
                sum_along_sources += capacities[v % period]
                v = shortcut(v,q,sign)
            assert maximum > 0
            histogram[maximum] += 1
    harmonic = sum((Fraction(count*denominator,maximum*source_period) for maximum,count in histogram.items()),Fraction(0))
    assert sum(histogram.values()) == len(bad)*q**(depth-1)*2**(depth-1)
    data,_=energy(depth,q,sign)
    assert Fraction(sum_along_sources,denominator*source_period)==Fraction(data["second"])
    assert harmonic >= Fraction(data["bound"])
    return {"L":depth,"endpoint_period":period,"source_period":source_period,"distinct_maxima":len(histogram),"harmonic":str(harmonic),"decimal":float(harmonic),"over_CS":float(harmonic/Fraction(data["bound"]))}


def tilted_diagonal(depth):
    """Exact p=3/4 survival and the last-time diagonal lower bound."""
    masses={0:1}
    for k in range(1,depth):
        new=Counter()
        for e,weight in masses.items():
            if 3**e > 2**k: new[e] += weight
            if 3**(e+1) > 2**k: new[e+1] += 3*weight
        masses=new
    survival=Fraction(sum(masses.values()),4**(depth-1))
    assert survival >= Fraction(1,4)
    return {"L":depth,"tilted_survival":float(survival),"last_diagonal_lower":float(survival/2)}


def trimmed_probe(depth,cap):
    data,rows=energy(depth,cap=cap)
    slope_count=1
    while 3**slope_count <= cap:
        slope_count += 1
    uniform_bound=cap*slope_count*sum(k//3+1 for k in range(depth))
    data["cap"]=cap
    data["uniform_capacity_bound"]=uniform_bound
    data["capacity_price_bound"]=str(Fraction(data["rho"])/uniform_bound)
    assert Fraction(data["second"]) <= uniform_bound*Fraction(data["mean"])
    # Affine offset h=c/w in every retained source prefix. This is not
    # a density approximation; it bounds each actual integer ancestor.
    for source in cells(depth,cap=cap)[0]:
        v,e=source,0
        for k in range(depth):
            w=Fraction(3**e,2**k)
            h=Fraction(v,w)-source
            assert 1 <= w <= cap and 0 <= h <= Fraction(k,3)
            e += v & 1
            v=shortcut(v)
    return data


def inverse_multiplicity_audit():
    """Independent exact inverse-tree enumeration; no source-cylinder code."""
    groups_checked=ancestors_checked=0
    largest_ratio=Fraction(0)
    for q,sign in [(3,1),(3,-1),(5,1)]:
        for hub in range(1,129):
            ancestors={hub}
            for k in range(1,19):
                next_ancestors={2*n for n in ancestors}
                next_ancestors.update((2*n-sign)//q for n in ancestors
                    if (2*n-sign)%q==0 and (2*n-sign)//q>0 and ((2*n-sign)//q)&1)
                ancestors=next_ancestors
                by_exponent=Counter()
                for source in ancestors:
                    v,e,valid=source,0,True
                    for j in range(1,k+1):
                        e += v&1
                        v=shortcut(v,q,sign)
                        if q**e <= 2**j:
                            valid=False
                            break
                    if not valid: continue
                    assert v==hub
                    h=Fraction(hub*2**k,q**e)-source
                    assert 0 <= sign*h <= Fraction(k,q)
                    by_exponent[e] += 1
                    ancestors_checked += 1
                for e,count in by_exponent.items():
                    groups_checked += 1
                    bound=k//q+1
                    assert count <= bound,(q,sign,hub,k,e,count,bound)
                    largest_ratio=max(largest_ratio,Fraction(count,bound))
    return {"hubs_per_sheet":128,"max_depth":18,"groups_checked":groups_checked,
            "admissible_ancestors_checked":ancestors_checked,"largest_count_over_bound":str(largest_ratio)}


def parameter_tube_audit():
    """Forward-source implementation for A<1 and signed intercepts."""
    groups_checked=0
    old_bound_hostile=None
    for q,sign,A in [(3,1,Fraction(1,16)),(3,-1,Fraction(1,4)),(5,1,Fraction(1,2))]:
        groups={}
        for source in range(1,8193):
            v,e=source,0
            for k in range(1,19):
                e += v&1
                v=shortcut(v,q,sign)
                w=Fraction(q**e,2**k)
                if w < A: break
                h=Fraction(v,w)-source
                assert 0 <= sign*h <= Fraction(k,q)/A
                groups.setdefault((v,k,e),[]).append(source)
        for (hub,k,e),sources in groups.items():
            bound=int(Fraction(k,q)/A)+1
            assert len(sources)<=bound
            groups_checked+=1
            if len(sources)>k//q+1 and old_bound_hostile is None:
                old_bound_hostile={"q":q,"sign":sign,"A":str(A),"hub":hub,"k":k,"e":e,"sources":sources,"old_bound":k//q+1,"repaired_bound":bound}
    assert old_bound_hostile is not None
    return {"sources_per_sheet":8192,"max_depth":18,"groups_checked":groups_checked,"A1_bound_hostile":old_bound_hostile}


def main():
    print("Exact incidence, not a random-orbit model. Badness is the slope-prefix predicate.")
    print("Direct pointwise energy audits:", sum(direct_audit(L,q,s) for q,s in [(3,1),(3,-1),(5,1)] for L in range(1,5)))
    print("Affine/cycle controls:", json.dumps(cohomology_probes(), sort_keys=True))
    print("Independent inverse multiplicity audit:",json.dumps(inverse_multiplicity_audit(),sort_keys=True))
    print("Parameterized affine tubes:",json.dumps(parameter_tube_audit(),sort_keys=True))
    plus = []
    for q,s,maxdepth in [(3,1,14),(3,-1,10),(5,1,10)]:
        print("SHEET", q, s)
        for L in range(1,maxdepth+1):
            result,_=energy(L,q,s)
            print(json.dumps(result, sort_keys=True))
            if q == 3 and s == 1: plus.append(result)
            if q == 3 and s == -1:
                assert result == plus[L-1], "sign conjugacy failed"
    print("Joint-incidence harmonic dual:")
    for L in range(1,7):
        print(json.dumps(incidence_harmonic(L),sort_keys=True))
    print("Tilted diagonal obstruction:")
    for L in [2,8,16,32,64,128,256]:
        print(json.dumps(tilted_diagonal(L),sort_keys=True))
    print("Trimmed slope tube and affine-offset capacity:")
    for cap in [2,4,8,16,32]:
        for L in [6,8,10,12,14]:
            print(json.dumps(trimmed_probe(L,cap),sort_keys=True))


if __name__ == "__main__":
    main()
