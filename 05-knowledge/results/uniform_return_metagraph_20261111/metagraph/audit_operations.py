"""Independent forward-operation and sharp-return controls for cap9.

Uses a start-end recurrence different from the producer's forward extension
DP. Imports only the literal family/dictionary data and canonical SCC word.
"""
import importlib.util
import itertools as it
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("cap9", HERE / "classify_cap9.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

@lru_cache(None)
def count(rows):
    n = len(rows)
    @lru_cache(None)
    def paths(mask, first):
        if mask == 1 << first:
            return 1
        rest = mask ^ (1 << first)
        return sum(paths(rest, v) for v in range(n)
                   if (rest >> v) & 1 and (rows[first] >> v) & 1)
    if not n:
        return 1
    return sum(paths((1 << n)-1, first) for first in range(n))

def grow(rows, mask):
    n = len(rows)
    return tuple(row | (0 if mask >> v & 1 else 1 << n)
                 for v, row in enumerate(rows)) + (mask,)

def main():
    br = m.birth_rules()
    ar = m.arc_rules()
    support_birth = {(r["source"],r["target"]) for r in br}
    support_arc = {(r["source"],r["target"]) for r in ar}
    birth_checks = valid_births = arc_checks = valid_arcs = 0
    for n in range(0, 8):
        for word in m.family_words(n):
            rows = m.from_word(word)
            for mask in range(1 << n):
                child = grow(rows, mask)
                h = count(child)
                birth_checks += 1
                if h > 9:
                    continue
                valid_births += 1
                target = m.word_of(child)
                need(target in set(m.family_words(n+1)), "forward birth missed family")
                parts = m.sccs(child)
                own = next(p for p in parts if n in p)
                if len(own) == 1:
                    need(any(target == word[:j]+("x",)+word[j:]
                             for j in range(len(word)+1)), "singleton insertion")
                else:
                    keep = tuple(v for v in own if v != n)
                    source_word = m.word_of(m.induced(rows, keep))
                    target_symbol = m.word_of(m.induced(child, own))
                    need(len(target_symbol)==1, "strong child symbol")
                    need(("".join(source_word),target_symbol[0]) in support_birth,
                         "forward birth absent from rooted dictionary")
    for n in range(2, 11):
        for word in m.family_words(n):
            rows = m.from_word(word)
            components = m.sccs(rows)
            owner = {v:i for i,part in enumerate(components) for v in part}
            for u,v in it.combinations(range(n),2):
                target = list(rows)
                target[u] ^= 1 << v
                target[v] ^= 1 << u
                target = tuple(target)
                h = count(target)
                arc_checks += 1
                if h > 9:
                    continue
                valid_arcs += 1
                i,j = owner[u],owner[v]
                keep = tuple(w for part in components[i:j+1] for w in part)
                need(len(keep)<=5, "unbounded affected interval under cap")
                left = m.word_of(m.induced(rows,keep))
                right = m.word_of(m.induced(target,keep))
                need(("".join(left),"".join(right)) in support_arc,
                     "missing local arc rule")
                expected = word[:i]+right+word[j+1:]
                need(m.word_of(target)==expected,"local arc changed exterior")
    # The exact upward excursion and an immutable-four-vertex alternative.
    aa = m.from_word(("A","A"))
    lifted = list(aa)
    lifted[0] ^= 1 << 3
    lifted[3] ^= 1
    lifted = tuple(lifted)
    returned = m.induced(lifted,[0,1,3,4,5])
    need((count(aa),count(lifted),count(returned))==(9,17,9),"9-17-9 return")
    need(m.word_of(returned)==("C",),"return is S5a")
    native = list(m.reverse(aa))
    for v in range(5):
        native[v] ^= 1 << (v+1)
        native[v+1] ^= 1 << v
    need(count(tuple(native))==17,"AA native escape")
    # First remove two vertices of the first triangle: xA remains literally.
    step1 = m.induced(aa,[0,1,3,4,5])
    step2 = m.induced(aa,[0,3,4,5])
    lower = []
    for mask in range(1 << 4):
        child = grow(step2,mask)
        if count(child)==9 and m.word_of(child)==("C",):
            lower.append(dict(mask=mask,rows=child))
    need(lower,"immutable four-vertex birth repair")
    need(count(step1)==count(step2)==3,"downward return")
    # Neutral kernel quotient directly from complete rule dictionary.
    edges=set()
    for rules in (br,ar):
        for rule in rules:
            src=tuple(rule["source"]);dst=tuple(rule["target"])
            if count(m.from_word(src))==count(m.from_word(dst)):
                edges.add(("".join(c for c in src if c!="x"),
                           "".join(c for c in dst if c!="x")))
    for rule in m.native_lines():
        if rule["source_h"]==rule["target_h"]:
            edges.add((rule["source"].replace("x",""),
                       rule["target"].replace("x","")))
    need(all(a==b or {a,b}=={"C","D"} for a,b in edges),
         "unexpected neutral kernel edge")
    need(("C","D") in edges or ("D","C") in edges,"missing C-D connection")
    result=dict(status="PASS",birth_queries=birth_checks,births_under_cap=valid_births,
                arc_queries=arc_checks,arcs_under_cap=valid_arcs,
                native_lines=Counter((str(x["n"])+"_"+x["color"]) for x in m.native_lines()),
                upward_return=dict(rows=[aa,lifted,returned],counts=[9,17,9],deleted_name=2),
                downward_return=dict(old_rows=aa,deleted_names=[1,2],
                                     retained_names=[0,3,4,5],
                                     step1=step1,step2=step2,births=lower,counts=[9,3,3,9]),
                neutral_kernel_edges=sorted(edges),
                scope="Finite independent operation audit; all-order completeness and height17 lower bound are written proofs.")
    (HERE/"operation_audit.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:result[k] for k in ("status","birth_queries","births_under_cap","arc_queries","arcs_under_cap")}))

if __name__=="__main__":
    main()
