#!/usr/bin/env python3
"""Independent subset DP, inherited-census checks, and reproducible summaries."""
from itertools import permutations
from pathlib import Path
import argparse, array, json, math

HERE=Path(__file__).resolve().parent

def adjacency(k, code):
    """No diagonal matching arcs needed; vertices 2i=x_i, 2i+1=y_i."""
    pairs=[(i,j) for i in range(k) for j in range(i+1,k)]
    m=len(pairs)
    a=[[0]*(2*k) for _ in range(2*k)]
    for z,(i,j) in enumerate(pairs):
        v=(code>>z)&1;b=(code>>(m+z))&1
        a[2*i][2*j]=v;a[2*j][2*i]=1-v
        a[2*i+1][2*j+1]=1-v;a[2*j+1][2*i+1]=v
        a[2*i][2*j+1]=a[2*j][2*i+1]=b
        a[2*j+1][2*i]=a[2*i+1][2*j]=1-b
    return a

def dp_h(k, code):
    a=adjacency(k,code)
    dp=[[0]*(2*k) for _ in range(1<<k)]
    for i in range(2*k):dp[1<<(i//2)][i]=1
    for mask in range(1,1<<k):
        for i,x in enumerate(dp[mask]):
            if not x:continue
            for j in range(2*k):
                if not mask>>(j//2)&1 and a[i][j]:
                    dp[mask|1<<(j//2)][j]+=x
    return sum(dp[-1][::2])

def frame_aut(k,code):
    a=adjacency(k,code)
    return [p for p in permutations(range(k))
            if all(a[2*i][2*j]==a[2*p[i]][2*p[j]]
                   and a[2*i][2*j+1]==a[2*p[i]][2*p[j]+1]
                   for i in range(k) for j in range(k) if i!=j)]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inherited",type=Path,
                        help="Optional directory containing pairing_n3.json through pairing_n9.json")
    args=parser.parse_args()
    inherited_dir=args.inherited
    if inherited_dir is None:
        candidates=[HERE.parent/"odd_blue_spectrum_hollow_cubes_20261020",
                    HERE.parent/"session14-odd-blue-pairing"]
        inherited_dir=next((p for p in candidates if (p/"pairing_n3.json").exists()),None)
    elif not inherited_dir.is_dir():
        parser.error("--inherited must name an existing census directory")
    reports=[]
    for k in range(1,6):
        data=json.loads((HERE/f"exact_k{k}.json").read_text())
        raw=array.array('H');raw.frombytes((HERE/f"h_k{k}.bin").read_bytes())
        assert len(raw)==1<<(k*(k-1))
        assert sum(raw)==math.factorial(k)*len(raw)
        checks=0
        if k<=4 and inherited_dir is not None:
            inherited=inherited_dir/f"pairing_n{2*k+1}.json"
            for r in json.loads(inherited.read_text())["rows"]:
                code=r["mask"]>>k
                assert raw[code]==r["symmetric_HP"]
                assert dp_h(k,code)==r["symmetric_HP"]
                checks+=1
        witness=data["max_masks"][0]
        assert dp_h(k,witness)==data["max_h"]
        aut=frame_aut(k,witness)
        regular=next(a for a in range(1<<k)
                     if sum(all((a>>i&1)==(a>>p[i]&1) for i in range(k)) for p in aut)==1)
        report={key:data[key] for key in ("k","external_frames","sum_h","sum_h_squared","max_h","elapsed_seconds")}
        report.update(witness=witness,witness_H_order=len(aut),regular_matching=regular,
                      inherited_class_checks=checks,maximizing_labelled_frames=len(data["max_masks"]),
                      inherited_census_available=inherited_dir is not None,
                      mean_h=math.factorial(k),variance_h=data["sum_h_squared"]/len(raw)-math.factorial(k)**2,
                      fixed_A_cross_mean_min=min(data["cross_sums_by_A"])/(1<<(k*(k-1)//2)),
                      fixed_A_cross_mean_max=max(data["cross_sums_by_A"])/(1<<(k*(k-1)//2)),
                      fixed_transitive_A_cross_mean=data["cross_sums_by_A"][0]/(1<<(k*(k-1)//2)))
        reports.append(report)
    (HERE/"verification.json").write_text(json.dumps(reports,indent=2)+"\n")
    print(json.dumps(reports,indent=2))

if __name__=="__main__":main()
