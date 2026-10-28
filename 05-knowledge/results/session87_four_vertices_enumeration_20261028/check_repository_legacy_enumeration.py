#!/usr/bin/env python3
"""Repo-relative regression for Session87's two narrowly scoped legacy repairs.

Default inputs are the live 04-computation sources and frozen Session22 baseline.
Explicit source overrides allow testing proposed copies before publication.
No obsolete whole-source snapshots are needed; C probes build only in /tmp.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import tempfile


def need(test, message):
    if not test:
        raise RuntimeError(message)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parts(n, upper=None):
    if not n:
        yield ()
        return
    top = min(n, n if upper is None else upper)
    top -= (top + 1) % 2
    for size in range(top, 0, -2):
        for m in range(1, n // size + 1):
            for tail in parts(n - size * m, size - 2):
                yield (size,) * m + tail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path)
    ap.add_argument("--python-source", type=Path)
    ap.add_argument("--c-source", type=Path)
    ap.add_argument("--baseline", type=Path)
    args = ap.parse_args()
    repo = args.repo or next((p for p in Path(__file__).resolve().parents
                             if (p / "04-computation").is_dir()), None)
    need(repo is not None, "Pass --repo for a non-repository invocation")
    py = (args.python_source or repo / "04-computation/a000568_speedup.py").resolve()
    c = (args.c_source or repo / "04-computation/k_uniform_fast_enum.c").resolve()
    base = (args.baseline or repo / "05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py").resolve()
    fixed = load(py, "repaired_legacy")
    refs, _ = load(base, "degree_local_baseline").local_dp(30)
    for n in range(31):
        emitted = list(fixed.partitions_odd_parts(n))
        got = {tuple(p) for p in emitted}
        need(got == set(parts(n)) and len(got) == len(emitted), f"Partition bijection n={n}")
        need(all(sum(p) == n and all(s > 0 and s % 2 for s in p)
                 and list(p) == sorted(p, reverse=True) for p in emitted), f"Domain n={n}")
        need(fixed.a000568_davis(n) == refs[n], f"Davis baseline n={n}")
        if n <= 20:
            need(fixed.a000568_dirichlet(n) == refs[n], f"Dirichlet n={n}")
    need(fixed.a000568_davis(5) == 12, "Parity repair hostile")
    with tempfile.TemporaryDirectory(prefix="session87-regression-") as tmp:
        tmp = Path(tmp)
        binary = tmp / "guard"
        subprocess.run(["gcc", "-O2", "-std=gnu11", str(c), "-lgmp", "-lpthread", "-o", str(binary)], check=True)
        # Recreate the discarded legacy generator only to authenticate its known failure.
        driver = tmp / "probe.c"
        driver.write_text('#define main legacy_main\n#include ' + json.dumps(str(c)) + '''
#undef main
static void old_capped(int r,int *prefix,int plen) {
  if(!r) { if(num_comps>=10000)return; memcpy(comp_data[num_comps],prefix,plen*sizeof(int));
    comp_len[num_comps++]=plen; return; }
  for(int v=r;v>=1;v--) {prefix[plen]=v;old_capped(r-v,prefix,plen+1);}
}
int main(void) {
  int prefix[MAX_K];
  for(int k=1;k<=14;k++) {
    global_k=k;num_comps=0;gen_compositions(k,prefix,0);
    if(num_comps!=(1<<(k-1))) return 2;
  }
  global_k=15;num_comps=0;old_capped(15,prefix,0);
  int ps[3]={1,2,3};parts_list=ps;
  int pk[3]={2,1,0},pm[3]={1,1,12};
  long long got=compute_ck(pk,pm,3);
  int pk_asc[3]={0,1,2},pm_asc[3]={12,1,1};
  long long reordered=compute_ck(pk_asc,pm_asc,3);
  printf("k1..14 composition counts complete; legacy k15 descending=%lld ascending=%lld true93\\n",got,reordered);
  return !(num_comps==10000 && got==92 && reordered==93);
}
''')
        probe = tmp / "probe"
        subprocess.run(["gcc", "-O2", "-std=gnu11", str(driver), "-lgmp", "-lpthread", "-o", str(probe)], check=True)
        witness = subprocess.run([str(probe)], capture_output=True, text=True, check=True).stdout.strip()
        cases = [(1, n, n + 1) for n in range(1, 7)]
        cases += [(2, n, x) for n, x in zip(range(2, 8), [2, 4, 11, 34, 156, 1044])]
        cases += [(3, 3, 2), (3, 4, 5), (3, 5, 34), (14, 14, 2), (15, 14, 1)]
        for k, n, expected in cases:
            run = subprocess.run([str(binary), str(k), str(n), "1"], capture_output=True, text=True, check=True, timeout=30)
            match = re.search(r"a\(\d+,\d+\) = (\d+)", run.stdout)
            need(match is not None and int(match.group(1)) == expected, f"C exact control {k,n}")
        rejects = [("15", "15", "1"), ("15", "17", "1"), ("20", "20", "1"),
                   ("0", "0", "1"), ("-1", "3", "1"), ("3", "-1", "1"),
                   ("3", "5", "0"), ("3", "5", "-2"), ("3", "oops", "1"),
                   ("3", "5", "9999999999999999999999999999999999")]
        for argv in rejects:
            run = subprocess.run([str(binary), *argv], capture_output=True, text=True, timeout=10)
            need(run.returncode != 0 and "ERROR" in run.stderr and "a(" not in run.stdout,
                 f"Unsafe input not rejected: {argv}")
    print(json.dumps({"status": "PASS", "python_generator_and_davis_orders": [0, 30],
                      "python_dirichlet_orders": [0, 20], "c_exact_controls": len(cases),
                      "c_rejected_controls": len(rejects), "compiled_probe": witness,
                      "scope": "Parity generator; C composition capacity and input validation only",
                      "not_validated": ["legacy CRT", "C arithmetic overflow", "C bucket caps"],
                      "sha256": {"python_source": hashlib.sha256(py.read_bytes()).hexdigest(),
                                 "c_source": hashlib.sha256(c.read_bytes()).hexdigest(),
                                 "baseline": hashlib.sha256(base.read_bytes()).hexdigest()}}, indent=2))


if __name__ == "__main__":
    main()
