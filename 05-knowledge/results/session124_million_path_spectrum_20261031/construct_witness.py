#!/usr/bin/env python3
"""Produce an actual adjacency matrix for any allowed odd target <=999999.

This uses an exact finite word tree (no depth cutoff), the frozen four block
repairs, and the fourteen ordinal-product recipes. It imports no DP verifier.
"""
import argparse
from itertools import combinations
import json
from math import prod
from pathlib import Path

BASE = Path(__file__).resolve().parent
FACTORS = {34869:[3,59,197],85885:[5,89,193],116587:[23,37,137],153451:[173,887],245923:[31,7933],535405:[5,13,8237],595543:[13,61,751],609873:[3,11,18481],618723:[49,9,23,61],767135:[5,153427],782103:[3,260701],824119:[149,5531],826373:[353,2341],868367:[157,5531]}

def marked(seed, blocks):
    n, code, cut = seed
    A = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(combinations(range(n),2)):
        A[i][j] = (code>>k)&1
        A[j][i] = 1-A[i][j]
    bits = [(cut>>i)&1 for i in range(n)]
    for block in blocks:
        if block not in ('0','1','C0','C1','C2','C3'):
            raise ValueError('Invalid block')
        old = len(A); size = 3 if block.startswith('C') else 1
        A = [row+[1]*size for row in A]+[[0]*(old+size) for _ in range(size)]
        if size == 3:
            for i in range(3): A[old+i][old+(i+1)%3] = 1
            bits.extend(int(i<int(block[1])) for i in range(3))
        else:
            bits.append(int(block))
    return [row+[1-b] for row,b in zip(A,bits)]+[bits+[0]]

def singleton_recipe(target):
    roots = json.loads((BASE/'constructors/seed_roots.json').read_text())
    for root in roots:
        todo = [(tuple(root['state']), ())]
        while todo:
            state, word = todo.pop(); h,F,q,t = state
            if F > target: continue
            if F == target: return {'seed':root['witness'],'blocks':list(word)}
            if q > target: continue
            # Remove only exact self-loops; positive q doubles, so this is finite.
            for child, bit in (((h,F,2*q+h,0),'1'),((h,F+q+h-t,2*q,h),'0')):
                if child != state: todo.append((child,word+(bit,)))
    raise ValueError('Target not in singleton support')

def construct(target):
    if target in FACTORS:
        factors = FACTORS[target]
        if prod(factors) != target: raise RuntimeError('Invalid product')
        pieces = [construct(f) for f in factors]; A=[]
        for B,_ in pieces:
            old=len(A); size=len(B)
            A=[row+[1]*size for row in A]+[[0]*old+row for row in B]
        return A, {'target':target,'ordinal_product':[r for _,r in pieces]}
    repairs = {r['target']:r for r in map(json.loads,(BASE/'constructors/block_target_witnesses.jsonl').read_text().splitlines())}
    recipe = repairs.get(target) or singleton_recipe(target)
    A = marked(recipe['seed'],recipe['blocks'])
    return A, {'target':target,'seed':recipe['seed'],'blocks':recipe['blocks']}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('target',type=int)
    p.add_argument('--adjacency',type=Path,required=True)
    p.add_argument('--recipe',type=Path)
    args=p.parse_args(); t=args.target
    if not (1<=t<=999999 and t%2 and t not in (7,21)):
        p.error('Require a positive odd target <=999999 other than 7,21')
    A, recipe=construct(t); n=len(A)
    if any(A[i][i] for i in range(n)) or any(A[i][j]+A[j][i]!=1 for i in range(n) for j in range(i+1,n)):
        raise RuntimeError('Invalid constructed tournament')
    args.adjacency.write_text(f'{n} {t}\n'+'\n'.join(''.join(map(str,row)) for row in A)+'\n')
    if args.recipe: args.recipe.write_text(json.dumps(recipe,sort_keys=True)+'\n')
    print(f'constructed order {n} target {t}')
if __name__=='__main__':main()
