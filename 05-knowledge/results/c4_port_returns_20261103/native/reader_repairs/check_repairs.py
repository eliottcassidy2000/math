"""Actual retained-name reader repairs, with bounded guarding and exact replay."""
from pathlib import Path
from functools import lru_cache
from collections import Counter
from itertools import combinations
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
PACKET = HERE.parents[1]
RESULTS = PACKET.parent
sys.path.insert(0, str(HERE.parent))
from check_native import native, legal, flip, add_named
sys.path.insert(0, str(RESULTS/"recursive_seed_switching_20261103/transport"))
import fast_gap
from fast_gap import Engine
spec = importlib.util.spec_from_file_location("inherited_traced_native", RESULTS/"adaptive_metagraph_reformulation_20261103/metagraph_probe.py")
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


class Budget(Exception):
    pass


def bounded_path(rows, u, v, budget=1000):
    full = (1 << len(rows))-1
    count = 0
    @lru_cache(None)
    def go(mask, last):
        nonlocal count
        count += 1
        if count > budget:
            raise Budget
        if mask == full:
            return ()
        allowed = ((1 << v) & ~mask) if last == u else ((full if last < 0 else rows[last]) & ~mask & ~(1 << v))
        while allowed:
            bit = allowed & -allowed
            allowed -= bit
            w = bit.bit_length()-1
            tail = go(mask | bit, w)
            if tail is not None:
                return (w,)+tail
        return None
    try:
        return go(0, -1)
    except Budget:
        return None


def insertion_path(rows, vertices):
    path = []
    for v in vertices:
        at = next((i for i, w in enumerate(path) if rows[v] >> w & 1), len(path))
        path.insert(at, v)
    need(all(rows[a] >> b & 1 for a,b in zip(path,path[1:])), "insertion-built path")
    return path


def compile_arc(rows, names, a, b, helper):
    u, v = names.index(a), names.index(b)
    if not rows[u] >> v & 1:
        u, v = v, u
        a, b = b, a
    source = tuple(rows)
    for inverse in (False, True):
        graph = old.flip(source,u,v) if inverse else source
        x,y = (v,u) if inverse else (u,v)
        word = bounded_path(graph,x,y)
        if word is None:
            continue
        U,V,W,p,q,case = old.traced_gadget(graph,x,y,word)
        if inverse:
            U,V,W,p,q = W,V,U,q,p
        need(tuple(r & ((1 << len(names))-1) for r in U[:-1]) == source, "actual helper birth source")
        local = names+[helper]
        return {"old_arc": [a,b], "case": case, "inverse": inverse,
                "births": [{"name": helper, "incoming": [names[i] for i in range(len(names)) if not U[-1] >> i & 1]}],
                "guards": [[local[i] for i in word] for word in (p,q)], "deletions": [helper]}
    return bridge_program(rows,names,a,b,helper)


def bridge_program(rows,names,a,b,helper):
    """Guaranteed fallback: two helpers, no Hamilton-language search."""
    u,v=names.index(a),names.index(b)
    need(rows[u] >> v & 1,"bridge input is the actual oriented arc")
    z,h = helper,helper+1
    rest = insertion_path(rows,[i for i in range(len(names)) if i not in (u,v)])
    word = [a,b,z]+[names[i] for i in rest]
    return {"old_arc": [a,b], "case": "bridge_then_append", "inverse": False,
            "births": [{"name":z,"incoming":[b]},
                       {"name":h,"incoming":[w for w in names+[z] if w != a]}],
            "guards": [word+[h], [a,h]+word[1:]], "deletions":[h,z]}


def replay(rows, names, program):
    original, original_names = list(rows), list(names)
    for item in program["births"]:
        need(item["name"] not in names, "fresh helper")
        rows,names = add_named(rows,names,item["name"],item["incoming"])
    for word in program["guards"]:
        rows = native(rows,names,word)
    need(program["deletions"] == list(reversed(names[len(original_names):])), "only helper names deleted")
    full = (1 << len(original_names))-1
    rows = [r & full for r in rows[:len(original_names)]]
    need(rows == flip(original,original_names,*program["old_arc"]), "exact old-arc repair")
    return rows


def feedback(rows):
    n = len(rows)
    @lru_cache(None)
    def search(removed,budget):
        active = [i for i in range(n) if not removed >> i & 1]
        triangle = None
        for i,a in enumerate(active):
            for j in range(i+1,len(active)):
                b = active[j]
                u,v = (a,b) if rows[a] >> b & 1 else (b,a)
                for c in active[j+1:]:
                    if rows[v] >> c & 1 and rows[c] >> u & 1:
                        triangle = (a,b,c)
                        break
                if triangle: break
            if triangle: break
        if triangle is None:
            return tuple(i for i in range(n) if removed >> i & 1)
        if budget == 0:
            return None
        for v in triangle:
            answer = search(removed | (1 << v),budget-1)
            if answer is not None:
                return answer
        return None
    for size in range(n+1):
        result = search(0,size)
        if result is not None:
            return result
    raise RuntimeError("feedback search exhausted")


def sidecar(engine,rows,names,ports):
    rows = tuple(rows)
    fb = feedback(rows)
    h = engine.count(rows,fb)
    observations = {}
    for port in ports:
        r = names.index(port)
        e = engine.count(rows,fb,1,(),((0,r),))
        s2 = engine.count(rows,fb,2,((0,r),),())
        f = engine.count(rows,fb,2,(),((0,r),))
        need((s2-e) % 2 == 0, "integral golden decoder")
        x,y = (s2-e)//2,(s2+e)//2
        a,b,c = (x+4*(h-y)) % 7,y % 7,(h-y) % 7
        observations[str(port)] = {"E2":f,"H_X_Y":[h,x,y],"F49_plus_charge":[a,b,c],
                                   "F49_frobenius_plus_charge":[(a+4*b)%7,-b%7,c],
                                   "F343_vector":[h%7,x%7,y%7],
                                   "F343_frobenius":[(h+4*x+6*y)%7,2*x%7,(3*x+4*y)%7]}
    pairs = [(99,0),(0,99),(99,2),(2,99),(0,4),(4,0)]
    j = {f"{a},{b}":engine.count(rows,fb,2,((1,names.index(b)),),((0,names.index(a)),)) for a,b in pairs}
    return {"H":h,"feedback_names":[names[i] for i in fb],"ports":observations,"ordered_J":j}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=HERE/"receipt.json")
    parser.add_argument("--clocks",type=int,nargs="*",default=list(range(1,11)))
    args = parser.parse_args()
    witness_path = PACKET/"reader/witnesses.json"
    witnesses = json.loads(witness_path.read_text())
    records=[]
    bridge_controls=0
    for n in range(2,5):
        ns=list(range(n))
        for mask in range(1 << (n*(n-1)//2)):
            rr=old.decode(n,mask)
            for a,b in combinations(ns,2):
                if not rr[a] >> b & 1:
                    a,b=b,a
                replay(rr,ns,bridge_program(rr,ns,a,b,99000))
                bridge_controls+=1
    need(bridge_controls==410,"all labeled order2..4 arc controls")
    temporary=tempfile.TemporaryDirectory(prefix="reader_gap_")
    binary=Path(temporary.name)/"gap.exe"
    subprocess.run(["g++","-O3","-std=c++17",str(fast_gap.SOURCE),"-o",str(binary)],check=True)
    fast_gap.executable=lambda:binary
    engine=Engine()
    cases=Counter()
    try:
        for number,record in enumerate(witnesses["records"]):
            if record["t"] not in args.clocks:
                continue
            source,target=record["source"],record["target"]
            rows,names=list(source["rows"]),list(source["names"])
            original_names=list(names)
            target_index={name:i for i,name in enumerate(target["names"])}
            births=[]
            for new in target["fresh_names"]:
                incoming=[old for old in names if target["rows"][target_index[old]] >> target_index[new] & 1]
                rows,names=add_named(rows,names,new,incoming)
                births.append({"name":new,"incoming":incoming})
            expected=[sum(1 << j for j,b in enumerate(names) if target["rows"][target_index[a]] >> target_index[b] & 1) for a in names]
            programs=[]
            for i,a in enumerate(names):
                for j in range(i+1,len(names)):
                    b=names[j]
                    if ((rows[i]^expected[i]) >> j)&1:
                        program=compile_arc(rows,names,a,b,30000+1000*number+2*len(programs))
                        rows=replay(rows,names,program)
                        programs.append(program)
                        cases[program["case"]+("_inverse" if program["inverse"] else "_forward")]+=1
            need(rows==expected,"exact named target endpoint")
            need(set(original_names)<=set(names) and 99 in names,"all old names/root retained")
            ports=[99,0,1,2,3,4]
            src_sidecar=sidecar(engine,source["rows"],source["names"],ports)
            dst_sidecar=sidecar(engine,rows,names,ports)
            need(src_sidecar["H"]==source["H"],"independent actual source count")
            need(dst_sidecar["H"]==record["target_H"]==3*src_sidecar["H"]+4,"independent actual target count")
            result={"family":record["family"],"t":record["t"],"source":source,
                    "permanent_births":births,"programs":programs,"target_rows":rows,"target_names":names,
                    "declared_target":target,"source_sidecar":src_sidecar,"target_sidecar":dst_sidecar,
                    "native_complements":2*len(programs),"deleted_old_names":[],
                    "peak_order":len(names)+max((len(p["births"]) for p in programs),default=0)}
            records.append(result)
            print(f'{record["family"]} t{record["t"]}: {src_sidecar["H"]}->{dst_sidecar["H"]}; {len(programs)} arc gadgets; peak{result["peak_order"]}',flush=True)
    finally:
        engine.close()
        temporary.cleanup()
    pins=[Path(__file__),witness_path,HERE.parent/"check_native.py",
          RESULTS/"adaptive_metagraph_reformulation_20261103/metagraph_probe.py",
          RESULTS/"recursive_seed_switching_20261103/transport/fast_gap.cpp"]
    report={"status":"FINITE-EXACT actual reader repairs with independent native replay and gap counts",
            "scope":"Finite clocks only; no all-clock common word, branch closure, fixed-order black connectivity, or physical Frobenius action.",
            "cases":dict(cases),"records":records,
            "exhaustive_bridge_arc_controls_orders2_to4":bridge_controls,
            "sha256":{str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}}
    args.output.write_text(json.dumps(report,indent=2)+"\n")


if __name__=="__main__":
    main()
