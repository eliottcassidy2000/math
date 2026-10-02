#!/usr/bin/env python3
"""procgen_extcert_20261001_lrc15_archive_check.py -- the orchestrator's integrity and consistency audit of the 71 gate
certificates of J. Allikvere, "Fifteen lonely runners" (Zenodo record 22667683; repo LRC(15) = 14 speeds).

Own code; the package's audit script and pipeline are not read or run. Only data is used: the archives
scratch/lrc_certs/lrc15/evidence_<p>_s0.tgz (downloaded with the owner's permission, MD5-checked against Zenodo) and
SHA256SUMS.txt. Each archive is STREAMED (never unpacked); memory stays small.

Per archive:
  1. SHA-256 of the archive equals its SHA256SUMS.txt line.
  2. MANIFEST_SHA256.json against the archive contents. Every member file of the subdirectories (rows_ir/, rows_km1/,
     filt/, kills/) is hashed. The multiset {(path, sha256)} is compared with the manifest's by count and by an
     order-independent digest: the sum mod 2^256 of sha256(path NUL hash). Top-level files are compared one by one.
     The manifest's tool entries (../bgk14, ../cascade_k14, ../tight_lift_15) are collected and must agree across all
     gates.
  3. SUMMARY.json: k = 14, p = the file's prime, status GATE_CLOSED, smax = 12, variant 'decomposition' iff no cover
     on <= 12 classes was found (covers_upto_smax = 0), claim type, persistent-orbit count = number of kill records.
  4. Every kill record kills/kill_orbit*.txt: header P = p, K = 14, L = 15. 'J(K,p) empty' gates need
     improper_after_gcd = 0 in every record; few-exception-shift gates need improper_after_neargcd = 0 in every record.
     Strict gates: the only persistent orbit is (1, 2, .., 14), with 3 / 5 witness-free lifts at levels 3 / 5,
     15 CRT candidates, 7 witness-free level-15 lifts, all settled by the gcd clause (paper, Section 4.3).
usage: procgen_extcert_20261001_lrc15_archive_check.py [workers]      (default 2)
"""
import hashlib
import json
import multiprocessing as mp
import os
import re
import sys
import tarfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, '..', '..', 'scratch', 'lrc_certs', 'lrc15')
MOD = 1 << 256
PAIR = re.compile(rb'"([^"\\]+)"\s*:\s*"([0-9a-f]{64})"')
SUBDIRS = ('rows_ir/', 'rows_km1/', 'filt/', 'kills/')
KILLRE = re.compile(r'P=(\d+) K=(\d+) L=(\d+) .*?witness_free_level_3=(\d+) witness_free_level_5=(\d+) crt_candidates=(\d+)\s+'
                    r'nodes=\d+ witness_free_level_15=(\d+) improper_after_gcd=(\d+) improper_after_neargcd=(\d+)', re.S)


def pair_digest(path, h):
    return int.from_bytes(hashlib.sha256(path.encode() + b'\0' + h.encode()).digest(), 'big')


def sha256_file(fn):
    h = hashlib.sha256()
    with open(fn, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def parse_manifest_stream(f):
    """flat JSON dict path -> sha256, parsed in chunks"""
    top, tools = {}, {}
    dig, n = 0, 0
    buf = b''
    while True:
        chunk = f.read(1 << 20)
        if not chunk:
            break
        buf += chunk
        last = 0
        for m in PAIR.finditer(buf):
            path, h = m.group(1).decode(), m.group(2).decode()
            last = m.end()
            if path.startswith('../'):
                tools[path] = h
            elif path.startswith(SUBDIRS):
                dig = (dig + pair_digest(path, h)) % MOD
                n += 1
            else:
                top[path] = h
        buf = buf[last:]
    return top, tools, dig, n


def summary_head(txt):
    """the scalar fields of a huge SUMMARY.json, from its first megabytes (they precede the per-orbit lists)"""
    def grab(pat, conv=str):
        mm = re.search(pat, txt, re.S)
        return conv(mm.group(1)) if mm else None
    g = re.search(r'"gate"\s*:\s*\{\s*"k"\s*:\s*(\d+)\s*,\s*"p"\s*:\s*(\d+)', txt)
    pre = re.search(r'"precondition"\s*:\s*\{(.*?)\}', txt, re.S).group(1)
    reps = re.search(r'"representatives"\s*:\s*\[\s*"([^"]*)"', txt)
    return {'gate': {'k': int(g.group(1)), 'p': int(g.group(2))},
            'claim_type': grab(r'"claim_type"\s*:\s*"([^"]*)"'), 'status': grab(r'"status"\s*:\s*"([^"]*)"'),
            'precondition': {'smax': int(re.search(r'"smax"\s*:\s*(\d+)', pre).group(1)),
                             'covers_upto_smax': int(re.search(r'"covers_upto_smax"\s*:\s*(\d+)', pre).group(1)),
                             'variant': re.search(r'"variant"\s*:\s*"([^"]*)"', pre).group(1)},
            'base_family': {'irredundant_rows': grab(r'"irredundant_rows"\s*:\s*(\d+)', int)},
            'persistent_orbits': {'count': grab(r'"persistent_orbits"\s*:\s*\{\s*"count"\s*:\s*(\d+)', int),
                                  'representatives': [reps.group(1)] if reps else []}}


def check_archive(args):
    name, expected_sha = args
    p = int(name.split('_')[1])
    t0 = time.time()
    res = {'p': p, 'errors': []}
    err = res['errors'].append
    fn = os.path.join(PKG, name)
    if sha256_file(fn) != expected_sha:
        err('archive sha256 differs from SHA256SUMS.txt')
    prefix = f'out_{p}/'
    fdig, fn_sub = 0, 0
    top_files = {}
    man = None
    summary_raw = None
    kills = {'n': 0, 'bad_header': 0, 'max_gcd': 0, 'max_near': 0, 'records': []}
    with tarfile.open(fn, 'r|gz') as tf:
        for m in tf:
            if not m.isfile():
                continue
            if not m.name.startswith(prefix):
                err(f'member outside {prefix}: {m.name}')
                continue
            rel = m.name[len(prefix):]
            f = tf.extractfile(m)
            if rel == 'MANIFEST_SHA256.json':
                man = parse_manifest_stream(f)
                continue
            if rel == 'SUMMARY.json':
                # not in the manifest (covered by the archive hash); up to 3.6 GB (per-orbit lists), so only the head
                # is read when large: every scalar field precedes the long lists
                summary_raw = f.read(m.size if m.size <= (64 << 20) else (8 << 20))
                summary_full = m.size <= (64 << 20)
                top_files[rel] = None
                continue
            data = f.read()
            h = hashlib.sha256(data).hexdigest()
            if rel.startswith(SUBDIRS):
                fdig = (fdig + pair_digest(rel, h)) % MOD
                fn_sub += 1
                if rel.startswith('kills/') and rel.endswith('.txt'):
                    km = KILLRE.search(data.decode(errors='replace'))
                    kills['n'] += 1
                    if not km or int(km.group(1)) != p or int(km.group(2)) != 14 or int(km.group(3)) != 15:
                        kills['bad_header'] += 1
                        continue
                    g, nr = int(km.group(8)), int(km.group(9))
                    kills['max_gcd'] = max(kills['max_gcd'], g)
                    kills['max_near'] = max(kills['max_near'], nr)
                    if len(kills['records']) < 3:
                        base = re.search(r'base: ([\d ]+)', data.decode()).group(1).strip()
                        kills['records'].append((base,) + tuple(int(km.group(i)) for i in range(4, 10)))
            else:
                top_files[rel] = h
    if man is None:
        err('no MANIFEST_SHA256.json')
        return res
    mtop, tools, mdig, mn = man
    res['tools'] = tools
    if mn != fn_sub:
        err(f'manifest lists {mn} subdirectory files, archive has {fn_sub}')
    if mdig != fdig:
        err('subdirectory digest differs from the manifest')
    for path, h in mtop.items():
        if path == 'SUMMARY.json' or top_files.get(path) != h:
            err(f'top-level manifest entry {path} missing or different')
    res['top_not_in_manifest'] = sorted(set(top_files) - set(mtop))
    res['manifest_entries'] = mn + len(mtop) + len(tools)
    if summary_raw is None:
        err('no SUMMARY.json')
        return res
    if summary_full:
        S = json.loads(summary_raw)
    else:
        S = summary_head(summary_raw.decode(errors='replace'))
        res['summary_head_only'] = True
    del summary_raw
    g, pre = S.get('gate', {}), S.get('precondition', {})
    if g.get('k') != 14 or g.get('p') != p:
        err(f'gate fields {g}')
    if S.get('status') != 'GATE_CLOSED':
        err(f"status {S.get('status')}")
    if pre.get('smax') != 12:
        err(f"smax {pre.get('smax')}")
    variant = pre.get('variant', '')
    decomp = variant.startswith('decomposition')
    if decomp != (pre.get('covers_upto_smax') == 0):
        err(f"variant {variant} vs covers_upto_smax {pre.get('covers_upto_smax')}")
    ct = S.get('claim_type', '')
    strict = ct.startswith('J(K,p) empty')
    near = ct.startswith('divisibility')
    if strict == near:
        err(f'claim type {ct}')
    if strict != decomp:
        err('claim type and variant disagree (strict gates use the decomposition variant)')
    po = S.get('persistent_orbits', {})
    if po.get('count') != kills['n']:
        err(f"persistent_orbits.count {po.get('count')} vs {kills['n']} kill records")
    if kills['bad_header']:
        err(f"{kills['bad_header']} kill records with a bad header")
    if strict:
        if kills['max_gcd'] != 0:
            err('a strict gate has an improper lift after the gcd clause')
        reps = po.get('representatives', [])
        if po.get('count') != 1 or reps != ['1 2 3 4 5 6 7 8 9 10 11 12 13 14']:
            err(f'strict gate persistent orbits {po.get("count")} {reps[:2]}')
        if kills['records'] and kills['records'][0][1:] != (3, 5, 15, 7, 0, 0):
            err(f"strict gate kill record {kills['records'][0]}")
    else:
        if kills['max_near'] != 0:
            err('an improper lift remains after the few-exception shift lemma')
    res.update(strict=strict, variant=variant, covers=pre.get('covers_upto_smax'), orbits=kills['n'],
               max_gcd=kills['max_gcd'], max_near=kills['max_near'], sample=kills['records'][:1],
               irr_rows=S.get('base_family', {}).get('irredundant_rows'), secs=round(time.time() - t0, 1))
    return res


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    sums = [line.split() for line in open(os.path.join(PKG, 'SHA256SUMS.txt')) if line.strip()]
    jobs = sorted(((nm, h) for h, nm in sums), key=lambda x: -os.path.getsize(os.path.join(PKG, x[0])))
    t0 = time.time()
    results = []
    with mp.Pool(workers) as pool:
        for r in pool.imap_unordered(check_archive, jobs):
            results.append(r)
            print(f"p = {r['p']}: {'OK' if not r['errors'] else 'ERRORS ' + '; '.join(r['errors'])} | "
                  f"{'strict' if r.get('strict') else 'shift'} | variant {r.get('variant')} | covers<=12: {r.get('covers')} | "
                  f"irredundant rows {r.get('irr_rows')} | persistent orbits {r.get('orbits')} | max improper after gcd "
                  f"{r.get('max_gcd')}, after shift {r.get('max_near')} | manifest entries {r.get('manifest_entries')} | "
                  f"top-level files not in manifest {r.get('top_not_in_manifest')} | {r.get('secs')} s", flush=True)
    results.sort(key=lambda r: r['p'])
    tools = {tuple(sorted(r.get('tools', {}).items())) for r in results}
    nstrict = sum(1 for r in results if r.get('strict'))
    bad = [r['p'] for r in results if r['errors']]
    print(f'==== {len(results)} archives; errors at {bad}; strict gates {nstrict}, few-exception-shift gates {len(results) - nstrict}')
    print(f'==== shift gates: {[r["p"] for r in results if not r.get("strict")]}')
    print(f'==== persistent orbits in all: {sum(r.get("orbits") or 0 for r in results)}; manifest entries in all: '
          f'{sum(r.get("manifest_entries") or 0 for r in results)}')
    print(f'==== tool hashes identical across gates: {len(tools) == 1}: {sorted(tools)[0] if tools else None}')
    print(f'elapsed {time.time() - t0:.0f} s')
    ok = not bad and len(results) == 71 and nstrict == 52 and len(tools) == 1
    print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')


if __name__ == '__main__':
    main()
