#!/usr/bin/env python3
"""Independent exact overlap audit. Does not enumerate external frames."""
from fractions import Fraction
from itertools import permutations
from math import factorial, comb
from pathlib import Path
import argparse, hashlib, json

def coefficient_formula(k):
    result = Fraction(0)
    for j in range(k // 2 + 1):
        r = k - 2*j + 1
        coefficient = sum(comb(r,l) * comb(r+j-l-1,j-l)
                          for l in range(min(r,j)+1))
        result += Fraction(factorial(k-2*j)*coefficient,4**j)
    return factorial(k)*result

def direct_overlap(k):
    result = Fraction(0)
    for q in permutations(range(k)):
        pos = {v:i for i,v in enumerate(q)}
        shared = [(i,1 if pos[i+1] == pos[i]+1 else -1)
                  for i in range(k-1) if abs(pos[i+1]-pos[i]) == 1]
        same_end = q[-1] == k-1
        numerator = 0
        for d in range(1 << k):
            if same_end and (d >> (k-1)) & 1:
                continue
            product = 1
            for i,epsilon in shared:
                a,b = (d >> i)&1,(d >> (i+1))&1
                if a == b:
                    if (a == 0) == (epsilon == 1):
                        product *= 2
                    else:
                        product = 0
                        break
            numerator += product
        result += Fraction(numerator,1 << (k-1 if same_end else k))
    return factorial(k)*result

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--census-dir',type=Path)
    p.add_argument('--max-k',type=int,default=7)
    args = p.parse_args()
    if args.max_k < 1:
        p.error("--max-k must be at least 1")
    results = []
    for k in range(1,args.max_k+1):
        formula = coefficient_formula(k)
        overlap = direct_overlap(k)
        if formula != overlap:
            raise ValueError(f"k={k}: coefficient formula {formula} differs from overlap sum {overlap}")
        census_checked = False
        census_sha256 = None
        if args.census_dir and k <= 5:
            source = args.census_dir/f'exact_k{k}.json'
            raw = source.read_bytes()
            data = json.loads(raw)
            if data['k'] != k or data['external_frames'] != 1 << (k*(k-1)):
                raise ValueError(f'{source}: incorrect k or complete-frame universe size')
            census_moment = Fraction(data['sum_h_squared'],data['external_frames'])
            if formula != census_moment:
                raise ValueError(f'k={k}: coefficient formula {formula} differs from census moment {census_moment}')
            census_sha256 = hashlib.sha256(raw).hexdigest()
            census_checked = True
        results.append({'k':k,'mean_h_squared':str(formula),
                        'variance_h':str(formula-factorial(k)**2),
                        'overlap_pair_check':True,'external_census_check':census_checked,
                        'external_census_sha256':census_sha256})
    print(json.dumps(results,indent=2))

if __name__ == '__main__':
    main()
