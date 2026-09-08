#!/usr/bin/env python3
"""Read-only historical pin verification; this does not replay mathematical proofs."""
from collections import Counter, defaultdict
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

sys.stdout.reconfigure(newline="\n")
BASE = "47bf7e6225ce8f1aa5a4eb23b9f74974a6f598f5"
SNAP = "840f1e47cc5c59a18bfd358737c8937c0e282920"
MANIFEST = "05-knowledge/results/planar_jc48_sep06_manifest.json"
HERE = Path(__file__).resolve()
STEM = HERE.stem
GATES = Counter()


def check(name, condition):
    GATES[name] += 1
    if not condition:
        raise RuntimeError("Historical verification failed: " + name)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


anchor = HERE.parent.parent if HERE.parent.name == "04-computation" else Path.cwd()
ROOT = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], cwd=anchor).decode().strip())


def git(*args, data=None):
    return subprocess.run(["git", *args], cwd=ROOT, input=data,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          check=True).stdout


for commit in (BASE, SNAP):
    check("pinned_commit_identity", git("rev-parse", commit+"^{commit}").decode().strip() == commit)

manifest_raw = {commit:git("show", commit+":"+MANIFEST) for commit in (BASE,SNAP)}
manifests = {commit:json.loads(raw) for commit,raw in manifest_raw.items()}
tables = {}
for commit, count in ((BASE,375),(SNAP,395)):
    entries = manifests[commit]["artifacts"]
    table = {x["path"]:x for x in entries}
    check("unique_manifest_paths", len(table)==len(entries))
    check("declared_manifest_count", len(table)==count==manifests[commit]["artifact_count"])
    tables[commit] = table

# Check filename namespaces AND YAML IDs in the pinned tree, not current HEAD.
all_names = git("ls-tree","-r","--name-only",SNAP,"01-canon/theorems").decode().splitlines()
pattern = r'''^id:[[:space:]]*['"]?THM-44(59|60|61|62)['"]?[[:space:]]*$'''
yaml_hits = git("grep","-l","-E",pattern,SNAP,"--","01-canon/theorems").decode().splitlines()
yaml_paths = [x.removeprefix(SNAP+":") for x in yaml_hits]
check("yaml_search_snapshot_prefix", all(x.startswith(SNAP+":") for x in yaml_hits))
for num in range(4459,4463):
    check("unique_filename_namespace", sum(Path(x).name.startswith(f"THM-{num}-") for x in all_names)==1)


def batch_blobs(requests):
    requests = list(dict.fromkeys(requests))
    check("safe_batch_names", all("\n" not in x and "\r" not in x for x in requests))
    raw = git("cat-file","--batch",data=("\n".join(requests)+"\n").encode())
    position = 0
    result = {}
    for requested in requests:
        end = raw.index(b"\n",position)
        header = raw[position:end].decode().split()
        check("snapshot_blob_type", len(header)==3 and header[1]=="blob")
        size = int(header[2]); start = end+1
        content = raw[start:start+size]
        check("snapshot_blob_frame", len(content)==size and raw[start+size:start+size+1]==b"\n")
        result[requested] = content
        position = start+size+1
    check("snapshot_batch_consumed", position==len(raw))
    return result


requests = [commit+":"+p for commit in (BASE,SNAP) for p in sorted(tables[commit])]
requests += [SNAP+":"+p for p in yaml_paths]
blobs = batch_blobs(requests)
verified_artifacts = {}
for commit in (BASE,SNAP):
    checked = []
    for path, expected in sorted(tables[commit].items()):
        raw = blobs[commit+":"+path]
        check("artifact_sha256", sha(raw)==expected["sha256"])
        check("artifact_byte_count", len(raw)==expected["bytes"])
        checked.append(dict(path=path,bytes=len(raw),sha256=sha(raw)))
    verified_artifacts[commit] = checked

old,new = tables[BASE],tables[SNAP]
changed = [dict(path=p,old_sha256=old[p]["sha256"],new_sha256=new[p]["sha256"])
           for p in sorted(old.keys() & new.keys()) if old[p]["sha256"]!=new[p]["sha256"]]
added,removed = sorted(new.keys()-old.keys()),sorted(old.keys()-new.keys())
check("snapshot_transition_counts", len(changed)==4 and len(added)==20 and not removed)

# The four old primary migrations are exact status/footer edits, not theorem edits.
migrations = {
    "planar_jc48_sep07_supplier_genus.md":(
        "The source and output are frozen. Independent analytic audit is pending.",
        "The source and output are frozen; the independent analytic audit below is accepted."),
    "planar_jc48_sep08_const_d_translation.md":(
        "The source and output are frozen. Independent analytic/source\nreview is pending; the primary remains RESERVED until that review\nis complete. The fixed-zero theorem is PROVED and independently\naudited.",
        "The source and output are frozen. The independent analytic/source\nreview below is accepted; this all-point theorem and its fixed-zero\nsupplier are PROVED and independently audited."),
    "planar_jc48_sep08_polynomial_weight_gate.md":(
        "The source/output are frozen. The proof and applications remain\nRESERVED until independent audit and parent-owned status promotion.",
        "The source/output are frozen. The proof and applications are PROVED\nfollowing the independent audit and parent-owned promotion below."),
    "planar_jc48_sep08_seven_one_two_locations.md":(
        "The complete independent analytic/source audit remains pending.\nThis candidate is outside the proved dependency graph until accepted.",
        "The complete independent analytic/source audit below is accepted.\nThis theorem is in the proved dependency graph."),
}
check("status_migration_path_set", {Path(x["path"]).name for x in changed}==set(migrations))
for row in changed:
    path=row["path"]; before,after=migrations[Path(path).name]
    a,b=blobs[BASE+":"+path],blobs[SNAP+":"+path]
    check("status_migration_unique_text", a.count(before.encode())==1)
    check("status_only_migration", a.replace(before.encode(),after.encode())==b)

canonical = []
by_id = defaultdict(list)
for path in yaml_paths:
    raw=blobs[SNAP+":"+path]; text=raw.decode("utf-8")
    match=re.search(r'''^id:\s*['"]?(THM-\d+)['"]?\s*$''',text,re.M)
    check("canonical_yaml_id_present", match is not None)
    by_id[match.group(1)].append(path)
for num in range(4459,4463):
    label=f"THM-{num}"
    check("unique_yaml_id", len(by_id[label])==1)
    path=by_id[label][0]; raw=blobs[SNAP+":"+path]; text=raw.decode("utf-8")
    check("yaml_matches_filename", Path(path).name.startswith(label+"-"))
    row=dict(id=label,path=path,bytes=len(raw),sha256=sha(raw),pins=[])
    for prefix in ("","spatial_"):
        for typ in ("script","output"):
            a=re.search(r"^"+prefix+typ+r": (.+)$",text,re.M)
            b=re.search(r"^"+prefix+typ+r"_sha256: ([0-9a-f]{64})$",text,re.M)
            check("canonical_pin_pair_complete", (a is None)==(b is None))
            if a is not None:
                row["pins"].append(dict(path=a.group(1).strip(),expected_sha256=b.group(1)))
    check("canonical_pin_count", len(row["pins"])==(4 if num==4459 else 2))
    canonical.append(row)

pins = batch_blobs([SNAP+":"+p["path"] for row in canonical for p in row["pins"]])
for row in canonical:
    for pin in row["pins"]:
        raw=pins[SNAP+":"+pin["path"]]
        check("canonical_artifact_sha256", sha(raw)==pin["expected_sha256"])
        pin.update(bytes=len(raw),actual_sha256=sha(raw))

record = dict(status="VERIFIED historical pins and exact status-only migrations; no mathematical proof replay",
    baseline_commit=BASE,reviewed_commit=SNAP,current_worktree_verification=False,
    manifest_sha256={commit:sha(raw) for commit,raw in manifest_raw.items()},
    baseline_planar_count=375,reviewed_planar_count=395,unchanged_old_count=371,
    changed_old=changed,added=added,removed=removed,verified_artifacts=verified_artifacts,
    canonical=canonical,gates=dict(sorted(GATES.items())),gate_total=sum(GATES.values()),
    source_sha256=sha(HERE.read_bytes()))
destination=HERE.parent
if destination.name=="04-computation":
    destination=destination.parent/"05-knowledge"/"results"
certificate=destination/(STEM+"_certificate.json")
certificate.write_bytes((json.dumps(record,indent=2,sort_keys=True)+"\n").encode("utf-8"))
print("PASS continuing11 historical incoming pins")
print("Baseline:",BASE)
print("Reviewed snapshot:",SNAP)
print("Planar artifacts: baseline 375; reviewed 395; unchanged old 371; status-only migrations 4; added 20; removed 0")
print("Canonical namespaces: 4 unique filenames and YAML IDs; producer/output pins: 10")
print("Current worktree contents: not verified; mathematical proofs: not replayed")
print("Always-active gates:",sum(GATES.values()))
print("Certificate:",certificate.name)
