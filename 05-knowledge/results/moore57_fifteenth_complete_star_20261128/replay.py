"""Replay the 15-star Moore partial graph. No solver is needed."""
from collections import Counter
from pathlib import Path
import argparse
import contextlib
import sys

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "moore57_quadratic_tensors_20261125"
if not OLD.is_dir():
    OLD = HERE
sys.path.insert(0, str(OLD))
from nested_core_check import orbit_chart, graph
from planted_stars import full_graph
from three_stars import IDENTITY, SLOTS, inverse, routes, symmetric_group_certificate
from translated_slots import contacts

def build(data_path):
    C, S = orbit_chart("period_fourteen_site.txt", 4, ((10, 45, 2), (27, 49, 0)))
    data = [list(map(int, s.split())) for s in data_path.read_text().splitlines() if s.strip()]
    num, = data[0]
    trades = data[1:1+num]
    raw = [tuple(x) for x in data[1+num:]]
    assert len(raw) == 52 and len(set(raw)) == 52
    old = {(h,y): contacts(C,h,y) for h in SLOTS for y in IDENTITY}
    forward0 = {(h,j): C[h,j] for h in SLOTS for j in (0,28)}
    assert len(set(old.values())) == 3024
    source, marker = (1, 1), (2, 7)
    old_star_triples = {u: [old[v] for v in vs] for u, vs in S.items()}
    flags = {v:u for u, vs in S.items() for v in vs}
    assert len(flags) == 14 * 53
    assert marker not in flags and not (set(raw) & set(flags))
    constraints = [r[source[1]] for j in (56,0,28) for r in routes(C,source[0],j)]
    raw_signature = Counter(h for h,_ in raw)
    counts = Counter()
    for h,k,y in trades:
        assert h!=k and h in SLOTS and k in SLOTS and not ({h,k} & {1,2})
        old_p,old_q = inverse(C[k,0]),inverse(C[k,28])
        comp = {y}
        check = [y]
        for t in check:
            for w in (old_p[C[h,0][t]],old_q[C[h,28][t]]):
                if w not in comp: comp.add(w);check.append(w)
        assert {flags[h,t] for t in comp if (h,t) in flags} == {
            flags[k,t] for t in comp if (k,t) in flags}
        counts[len(comp)] += 1
        for t in comp:
            a,b = flags.pop((h,t),None),flags.pop((k,t),None)
            if a is not None: flags[k,t]=a
            if b is not None: flags[h,t]=b
        for j in (0,28):
            p,q = list(C[h,j]),list(C[k,j])
            for t in comp:p[t],q[t]=q[t],p[t]
            C[h,j],C[k,j]=tuple(p),tuple(q)
            C[j,h],C[j,k]=inverse(p),inverse(q)
    fresh = {(h,y): contacts(C,h,y) for h in SLOTS for y in IDENTITY}
    assert set(fresh.values())==set(old.values())
    assert fresh[source]==old[source] and fresh[marker]==old[marker]
    assert constraints==[r[source[1]] for j in (56,0,28) for r in routes(C,source[0],j)]
    renamed = {v:u for u,v in fresh.items()}
    oldstars = {u:[renamed[v] for v in triples] for u,triples in old_star_triples.items()}
    assert flags == {v:u for u, vs in oldstars.items() for v in vs}
    newstar = [renamed[old[v]] for v in [*raw, marker]]
    allstars = {**oldstars,source:newstar}
    assert set(allstars)==set(S)|{source}
    assert len({v for vs in allstars.values() for v in vs})==15*53
    old_moves = sum(len(set(S[u])-set(oldstars[u])) for u in S)
    centered_moves = sum(a!=b for key in forward0 for a,b in zip(forward0[key], C[key]))
    for (h,y),targets in allstars.items():
        assert {k for k,_ in targets}==set(SLOTS)-{h}
        for j in (56,0,28):
            reached = [r[y] for r in routes(C,h,j)]
            reached += [C[k,j][z] for k,z in targets]
            assert sorted(reached)==list(IDENTITY)
    print("Trades",len(trades),"component sizes",sorted(counts.items()))
    print("Raw branch occupancy profile", sorted(Counter(raw_signature.get(b,0) for b in SLOTS if b not in (1,2)).items()))
    print("Old exterior contacts moved", old_moves, "center contacts moved", centered_moves)
    print("New contacts",sorted(newstar))
    assert full_graph(C,allstars)==(13284,241)
    A=graph(C, allstars)
    labels=[0]*len(A)
    for _ in range(10):
        ids={}
        new=[ids.setdefault((labels[v],tuple(sorted(labels[w] for w in row))),len(ids))
             for v,row in enumerate(A)]
        if len(set(new))==len(A):break
        labels=new
    assert len(set(new))==3250
    print("Unmarked color refinement is discrete")
    words={(0,28):(6,10,39),(0,56):(10,20,3),(28,56):(4,20,11)}
    for pair,word in words.items():
        print("Symmetric group", pair, symmetric_group_certificate(routes(C,*pair),word))
    return C, allstars

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("witness",nargs="?",type=Path,default=HERE/"witness.txt")
    parser.add_argument("--adjacency",action="store_true",help="format for frozen_hall.cpp")
    args=parser.parse_args()
    with contextlib.redirect_stdout(sys.stderr if args.adjacency else sys.stdout):
        C, S = build(args.witness)
    if args.adjacency:
        A=graph(C,S)
        print(len(A))
        for row in A: print(len(row),*sorted(row))

if __name__=="__main__": main()
