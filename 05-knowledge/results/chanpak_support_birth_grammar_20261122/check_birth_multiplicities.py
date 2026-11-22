"""Class census from exact ordered output of birth_classes.cpp."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path


def strong(A):
    full = (1<<len(A))-1
    forward = backward = 1
    for _ in A:
        for i, arc in enumerate(A):
            if forward>>i&1:
                forward |= arc
            if arc&backward:
                backward |= 1<<i
    return forward==backward==full


def self_converse(A):
    n = len(A)
    scores = [arc.bit_count() for arc in A]
    if sorted(scores)!=sorted(n-1-d for d in scores):
        return False
    pi = [-1]*n

    def visit(unassigned, fixed):
        if not unassigned:
            return True
        u = min(unassigned)
        for v in sorted(unassigned):
            if scores[u]+scores[v]!=n-1 or (u==v and fixed):
                continue
            candidate = {u:v,v:u}
            if any((A[x]>>y&1)==(A[z]>>pi[y]&1)
                   for x,z in candidate.items() for y in range(n) if pi[y]>=0):
                continue
            # A new swapped pair reverses its internal arc automatically.
            pi[u], pi[v] = v, u
            if visit(unassigned-set(candidate), fixed or u==v):
                return True
            pi[u] = pi[v] = -1
        return False
    return visit(set(range(n)), False)


parser = argparse.ArgumentParser()
parser.add_argument("classes", type=Path)
parser.add_argument("--fibers", type=Path)
args = parser.parse_args()
by = defaultdict(Counter)
fixed = defaultdict(Counter)
decomposable = defaultdict(Counter)
birth = {}
orders = Counter()
for line in args.classes.open():
    n,h,text_rows = line.split()
    n,h = int(n),int(h)
    orders[n] += 1
    birth.setdefault(h,n)
    if birth[h]!=n:
        continue
    A = tuple(map(int,text_rows.split(",")))
    by[n][h] += 1
    fixed[n][h] += self_converse(A)
    if not strong(A):
        decomposable[n][h] += 1

known = [1,1,2,4,12,56,456,6880,191536]
if sorted(orders)!=list(range(1,max(orders)+1)) or any(orders[n]!=known[n-1] for n in orders):
    raise ValueError(orders)

print("n values classes one-class-values self-converse decomposable maxima(H,classes,self-converse)")
for n in sorted(orders):
    c = by[n]
    if any((v-fixed[n][h])%2 for h,v in c.items()):
        raise ValueError(n)
    maxima = [(h,v,fixed[n][h]) for h,v in sorted(c.items()) if v==max(c.values())]
    print(n,len(c),sum(c.values()),sum(v==1 for v in c.values()),sum(fixed[n].values()),
          sum(decomposable[n].values()),maxima)

if args.fibers:
    with args.fibers.open("w") as out:
        print("H\tfirst_order\tclasses\tself_converse\tconverse_pairs\tdecomposable",file=out)
        for h,n in sorted(birth.items()):
            a,f = by[n][h],fixed[n][h]
            print(h,n,a,f,(a-f)//2,decomposable[n][h],sep="\t",file=out)
