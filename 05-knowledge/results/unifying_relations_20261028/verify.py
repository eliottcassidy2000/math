"""Validate frozen receipts and optionally replay both Python execution modes."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def need(ok,message):
    if not ok: raise RuntimeError(message)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--replay",action="store_true")
    args=parser.parse_args()
    count=0
    for line in (HERE/"SHA256SUMS").read_text(encoding="utf-8").splitlines():
        digest,name=line.split("  ",1)
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,"hash: "+name)
        count+=1
    for doc in HERE.glob("*.md"):
        for link in re.findall(r"\]\(([^)]+)\)",doc.read_text(encoding="utf-8")):
            if re.match(r"[a-zA-Z]+:",link):continue
            path=link.split("#")[0]
            if path:need((doc.parent/path).exists(),"local link: "+link)
    pins=json.loads((HERE/"SOURCES.json").read_text(encoding="utf-8"))
    for row in pins["repository_sources"]:
        blob=subprocess.check_output(["git","show",row["commit"]+":"+row["path"]],cwd=ROOT)
        need(hashlib.sha256(blob).hexdigest()==row["sha256"],"source pin: "+row["path"])
    normal=(HERE/"normal.json").read_bytes()
    need(normal==(HERE/"optimized.json").read_bytes(),"normal/optimized frozen equality")
    if args.replay:
        for flags in ([],["-O"]):
            output=subprocess.check_output([sys.executable,"-B",*flags,str(HERE/"check.py")],cwd=ROOT)
            need(output==normal,"live replay "+str(flags))
    print(json.dumps({"status":"PASS","packet_hashes":count,
                      "source_pins":len(pins["repository_sources"]),
                      "normal_optimized_equal":True,"replayed":args.replay}))

if __name__=="__main__": main()
