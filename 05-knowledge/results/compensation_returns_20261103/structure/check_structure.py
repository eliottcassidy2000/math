"""Bounded structural controls: deletion caps, persistent endpoints, migration."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]
OLD=RESULTS/"c4_port_returns_20261103/native/reader_repairs/check_entry_scope.py"
spec=importlib.util.spec_from_file_location("inherited_small_endpoint_dp",OLD)
dp=importlib.util.module_from_spec(spec);spec.loader.exec_module(dp)


def need(ok,why):
    if not ok: raise RuntimeError(why)


def birth(rows,out):
    n=len(rows)
    return [r | (0 if out >> i & 1 else 1 << n) for i,r in enumerate(rows)]+[out]


def c2(rows,aux):
    n=len(rows)
    return birth(birth(rows,1 << aux),1 << n)


def induced(rows,vertices):
    return [sum(1 << j for j,b in enumerate(vertices) if rows[a] >> b & 1) for a in vertices]


def flip(rows,a,b):
    out=list(rows);out[a]^=1 << b;out[b]^=1 << a;return out


def join(rows,a,b):
    full=(1 << len(rows))-1
    h,e=dp.tables(rows)
    _,s=dp.tables([full ^ r ^ (1 << i) for i,r in enumerate(rows)])
    return sum(e[mask][a]*s[full ^ mask][b] for mask in range(full+1))


def inspect(rows,names,mark=0,root=99):
    p=dp.profile(rows)
    _,ends=dp.tables(rows)
    w=names.index(mark)
    reduced=induced(rows,[i for i in range(len(rows)) if i!=w])
    return {"H":p["H"],"E_root":p["E2"][names.index(root)],
            "q":dp.tables(reduced)[0][-1],"e_mark":ends[-1][w],
            "start_defects":dict(zip(map(str,names),p["start_defects"]))}


def cap_graph(m,tail,bits):
    # Deleted carrier is TT_m -> C3 -> optional singleton tail.
    n=m+3+int(tail)
    rows=[0]*n
    for i in range(m): rows[i]=sum(1 << j for j in range(i+1,n))
    for i in range(3):
        rows[m+i]=1 << (m+(i+1)%3)
        if tail: rows[m+i]|=1 << (n-1)
    return birth(rows,((1 << n)-1)^bits)


def initial_cap(tail,out):
    n=3+tail
    rows=[(1 << ((i+1)%3)) | sum(1 << j for j in range(3,n)) for i in range(3)]
    rows += [sum(1 << j for j in range(i+1,n)) for i in range(3,n)]
    return birth(rows,out)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=HERE/"receipt.json")
    args=parser.parse_args()
    cap_checks=[]
    for tail in (False,True):
        for m in range(2,5):
            n=m+3+int(tail); tested=0;ready=0
            for bits in range(1 << n):
                if bits&1 or not bits&2: continue
                k=sum(bits >> (m+i)&1 for i in range(3))
                if tail:
                    if not bits >> (n-1)&1: continue
                elif not k: continue
                rows=cap_graph(m,tail,bits);p=dp.profile(rows);tested+=1
                need(p["E2"][0]==p["H"],"root first and second deletion vertex dominates mark")
                if tail:
                    need(p["start_defects"][n-1]==k,"one terminal singleton exposes cap endpoint mass")
                elif m>=3:
                    need(not any(x in (1,2) for x in p["start_defects"]),"terminal triangle star obstruction")
                ready+=any(x in (1,2) for x in p["start_defects"])
            cap_checks.append({"tail":tail,"TT_length":m,"stars":tested,"with_start1_or2":ready})
    need(cap_checks[0]["with_start1_or2"]>0,"m2 boundary is real")
    initial_checks=[]
    for b in range(2,5):
        tested=ready=0
        for out in range(1 << (b+3)):
            if not out&7 or out >> (b+2)&1:continue
            p=dp.profile(initial_cap(b,out));tested+=1
            actual=[i for i,x in enumerate(p["start_defects"]) if x in (1,2)]
            exception=(out&7).bit_count()==2 and out>>3==0
            need(actual==([b+3] if exception else []),"initial-cap exact apex exception classification")
            ready+=bool(actual)
        initial_checks.append({"tail":b,"strong_stars":tested,"apex_ready2":ready})
    apex_hostiles=[]
    for b in (1,2,3,4,9):
        p=dp.profile(initial_cap(b,3))
        need(p["H"]==6*2**b+1 and p["start_defects"][-1]==2,"arbitrary-length initial-cap apex hostile")
        apex_hostiles.append({"tail":b,"H":p["H"],"apex_start_defect":p["start_defects"][-1]})
    modules=[];rows=[0];aux=0
    for j in range(6):
        h,e=dp.tables(rows)
        need(h[-1]==3**j and e[-1][0]==1,"unbounded cap count and fixed anchor endpoint")
        modules.append({"j":j,"rows":rows,"anchor":0,"auxiliary":aux,"H":h[-1],"e_anchor":e[-1][0]})
        n=len(rows);rows=c2(rows,aux);aux=n
    wrong=c2(modules[2]["rows"],0)
    wrong_e=dp.tables(wrong)[1][-1][0]
    need(wrong_e==3,"at depth2 the endpoint anchor is not the transported auxiliary")
    path=RESULTS/"c4_port_returns_20261103/entry_selector/receipt.json"
    data=json.loads(path.read_text())
    records=[]
    for rec in data["records"][:2]:
        child=next(s["graph"] for s in rec["stages"] if s["stage"]=="C2_child")
        ns=child["names"];rr=child["rows"]
        source=inspect(rr,ns)
        need(source["H"] in (1409,1427) and source["q"]==3 and source["e_mark"]==1,"literal child source")
        stages=[]
        for cap in (201,401,400):
            a,b=ns.index(200),ns.index(cap)
            need(rr[a] >> b&1,"old reader200 precedes cap")
            before=inspect(rr,ns)
            jold,jnew=join(rr,a,b),join(rr,b,a)
            rr=flip(rr,a,b)
            after=inspect(rr,ns)
            need(after["H"]-before["H"]==jnew-jold,"fresh joint join-price of each flip")
            need(after["E_root"]==after["H"],"protected outer root stays zero-ready")
            stages.append({"arc":[200,cap],"old_J":jold,"reverse_J":jnew,"state":after})
        moved=inspect(rr,ns)
        need(moved["q"]==3 and moved["start_defects"]["200"]==1,"relocation prepares terminal reader200")
        records.append({"family":rec["family"],"source":source,"source_graph":child,
                        "three_flip_stages":stages,"relocated_graph":{"rows":rr,"names":ns},
                        "relocated":moved,"delta_H":moved["H"]-source["H"]})
    report={"status":"PROVED mechanisms with FINITE-EXACT bounded controls",
            "cap_star_controls":cap_checks,"zero_channel_modules":modules,
            "initial_cap_controls":initial_checks,"initial_cap_apex_hostiles":apex_hostiles,
            "wrong_auxiliary_anchor_endpoint_count":wrong_e,"literal_relocations":records,
            "sha256":{str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (Path(__file__),OLD,path)}}
    args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"cap_controls":cap_checks,"wrong_anchor_e":wrong_e,
                      "actual_relocations":[{"family":r["family"],"source_H":r["source"]["H"],
                          "new_H":r["relocated"]["H"],"delta":r["delta_H"],
                          "sigma200":r["relocated"]["start_defects"]["200"]} for r in records]},indent=2))


if __name__=="__main__":main()
