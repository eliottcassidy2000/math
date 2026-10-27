#!/usr/bin/env python3
"""Ensure the independent optimizer rejects four corrupt certificates."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("certificate")
    args=parser.parse_args()
    with open(args.certificate) as f:
        original=json.load(f)
    variants=[]
    d=copy.deepcopy(original); d["cases"][-1]["minimum"]+=1
    variants.append(("final_minimum",d))
    d=copy.deepcopy(original); d["cases"][-1]["subset_lower_bounds"][31]+=1
    variants.append(("interior_subset_optimum",d))
    d=copy.deepcopy(original); d["cases"][1]["witness"]["orientation"]=0
    variants.append(("illegal_quotient_orientation",d))
    d=copy.deepcopy(original); d["cases"][1]["witness"]["children"][1]["vertex"]=2
    variants.append(("duplicated_vertex",d))
    script=Path(__file__).with_name("independent_distances.py")
    rejected=[]
    with tempfile.TemporaryDirectory(prefix="spectral-corruption-") as tmp:
        for name,data in variants:
            path=Path(tmp)/(name+".json")
            path.write_text(json.dumps(data))
            for opt in (False,True):
                command=[sys.executable]+(["-O"] if opt else [])+["-B",str(script),"--certificate",str(path)]
                result=subprocess.run(command,capture_output=True,text=True)
                require(result.returncode!=0 and "RuntimeError" in result.stderr,"corruption accepted or wrong failure: "+name)
            rejected.append(name)
    print(json.dumps({"status":"PASS","corruptions_rejected_normal_and_optimized":rejected},indent=2,sort_keys=True))


if __name__=="__main__":
    main()
