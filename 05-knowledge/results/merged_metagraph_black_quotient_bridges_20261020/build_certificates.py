#!/usr/bin/env python3
"""Build compact exact certificates after running census.py through n=8."""
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def main():
    rows=[];hashes=[]
    for n in range(3,9):
        rel=Path("n8")/f"census_n{n}.json" if n==8 else Path(f"census_n{n}.json")
        raw=(ROOT/rel).read_bytes()
        hashes.append(hashlib.sha256(raw).hexdigest()+"  "+str(rel))
        d=json.loads(raw)
        nodes=d["nodes"]
        graphs={}
        for name,g in d["graphs"].items():
            graphs[name]={k:g[k] for k in ("vertices","edges","loops","support_pairs","component_count","cycle_rank")}
            graphs[name]["active_components"]=sum(c["edges"]>0 for c in g["components"])
            graphs[name]["component_vertex_sizes"]=dict(Counter(len(c["vertices"])for c in g["components"]))
            graphs[name]["bridges"]=[{"endpoints":[nodes[u],nodes[v]],
                                     "has_quotient_leaf":any(nodes[z]["reflection_pairs"]==1 for z in (u,v))}
                                    for u,v in g["bridges"]]
        blue_spectrum=dict(sorted(Counter(r["blue"]for r in nodes if r["SC"]).items()))
        max_b=max(blue_spectrum)
        rows.append({"n":n,"m":d["m"],"tilings":len(d["pi"]),"classes":len(d["classes"]),
                     "merged":len(nodes),"SC":sum(r["SC"]for r in nodes),
                     "node_types":dict(Counter(r["type"]for r in nodes)),
                     "graphs":graphs,"SC_blue_spectrum":blue_spectrum,
                     "SC_blue_spectrum_by_triangle_parity":{
                         p:dict(sorted(Counter(r["blue"]for r in nodes if r["SC"]and r["c3"]%2==p).items()))
                         for p in (0,1)},
                     "maximum_blue_witnesses":[{"node":r,"classes":[d["classes"][k]for k in r["classes"]]}
                                               for r in nodes if r["blue"]==max_b],
                     "all_Q_bridges_have_a_leaf":all(any(nodes[z]["reflection_pairs"]==1 for z in (u,v))
                                                   for u,v in d["graphs"]["paired_black"]["bridges"])})
    (ROOT/"compact_certificates.json").write_text(json.dumps(rows,indent=2)+"\n")
    (ROOT/"FULL_ATLAS_SHA256SUMS").write_text("\n".join(hashes)+"\n")


if __name__=="__main__":
    main()
