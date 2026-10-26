#!/usr/bin/env python3
"""Independent finite audit. Never imports or executes audited source code.
Usage: python3 -B check.py [path/to/session65-twentyfive]
"""
import hashlib
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
SOURCE = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE.parent.parent / 'session65-twentyfive'
CODE = '1RB1RD_1LC1RB_0RA1LC_0LQ1RE_0LF1RG_0LC1LF_0LF0LH_1LQ1LI_0RJ0LI_1RK0LJ_0RL0RS_1RL0RM_1RN1RM_0LO0LU_0LP1LO_1RH1LX_1LR1LQ_0RK0LT_1LR1RS_---1RC_1LV1LU_0LW0LJ_0RK0LW_1RY1LX_1RE1RY'
RENAMING = 'A D E G H I V J W K L M N U O P Q R X Y F S B C T'.split()

def need(ok, why):
    if not ok:
        raise ValueError(why)

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def sealed_integrity():
    allowed = (SOURCE / 'ALLOWLIST.txt').read_text().splitlines()
    actual = sorted(str(p.relative_to(SOURCE)) for p in SOURCE.rglob('*') if p.is_file())
    need(sorted(allowed) == actual and len(allowed) == len(set(allowed)) == 48, 'exact 48-file source allowlist')
    for name in allowed:
        p = SOURCE / name
        need(not p.is_symlink() and p.resolve().is_relative_to(SOURCE.resolve()), 'safe regular source path')
    sums = (SOURCE / 'SHA256SUMS').read_text().splitlines()
    for line in sums:
        expected, name = line.split('  ', 1)
        need(name in allowed and digest(SOURCE / name) == expected, 'source SHA256: ' + name)
    manifest = json.loads((SOURCE / 'MANIFEST.json').read_text())
    for item in manifest['files']:
        p = SOURCE / item['path']
        need(len(p.read_bytes()) == item['bytes'] and digest(p) == item['sha256'], 'source manifest: ' + item['path'])
    receipt = json.loads((SOURCE / 'RECEIPT.json').read_text())
    tree = json.loads((SOURCE / 'SOURCE_TREE.json').read_text())
    blobs = {x['path']: x['sha'] for x in tree['tree'] if x['type'] == 'blob'}
    for item in receipt['upstream_byte_checks']:
        p = SOURCE / item['path']
        raw = p.read_bytes()
        blob = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
        need(blob == item['expected_git_blob'] == blobs[item['path'].removeprefix('inputs/goldbach_tm/')], 'upstream Git blob consistency')
    need((SOURCE/'checks.json').read_bytes() == (SOURCE/'checks-O.json').read_bytes(), 'source normal/-O equality')
    return {'source_files': len(allowed), 'sha256_entries': len(sums), 'manifest_entries': len(manifest['files']), 'git_blob_checks': len(receipt['upstream_byte_checks'])}

def lean_table():
    text = (SOURCE / 'inputs/goldbach_tm/GoldbachTm/Tm25/TuringMachine25.lean').read_text()
    body = text.split('def machine : Machine\n', 1)[1].split('\ndef nth_cfg', 1)[0]
    table = {}
    for line in body.splitlines():
        if not line.startswith('| ⟨') or '_+25' in line:
            continue
        left, right = line.split('=>', 1)
        q = int(re.search(r'\d+', left).group())
        read = {'zero':0, 'one':1}[re.search(r'Γ\.(\w+)', left).group(1)]
        need((q, read) not in table, 'duplicate Lean branch')
        if right.strip() == 'none':
            table[q, read] = None
        else:
            target = int(re.search(r'\d+', right).group())
            write = {'zero':0, 'one':1}[re.search(r'Γ\.(\w+)', right).group(1)]
            direction = {'left':-1, 'right':1}[re.search(r'Turing.Dir\.(\w+)', right).group(1)]
            table[q, read] = (target, write, direction)
    need(set(table) == {(q,b) for q in range(25) for b in (0,1)}, 'all Lean slots exactly once')
    return table

def wiki_table():
    rows = CODE.split('_')
    need(len(rows) == 25 and all(len(x) == 6 for x in rows), 'wiki dimensions')
    out = {}
    for q, row in enumerate(rows):
        for bit in (0,1):
            v = row[3*bit:3*bit+3]
            out[q, bit] = None if v == '---' else (ord(v[2])-65, int(v[0]), {'L':-1, 'R':1}[v[1]])
    return out

def same_table(a, b):
    rho = [ord(x)-65 for x in RENAMING]
    if sorted(rho) != list(range(25)):
        return False
    for (q, r), inst in a.items():
        want = None if inst is None else (rho[inst[0]], inst[1], inst[2])
        if b[rho[q], r] != want:
            return False
    return True

def one(table, cfg):
    if cfg is None:
        return None
    q, h, cells = cfg
    instruction = table[q, int(h in cells)]
    if instruction is None:
        return None
    target, write, direction = instruction
    cells = set(cells)
    if write:
        cells.add(h)
    else:
        cells.discard(h)
    return target, h + direction, frozenset(cells)

def steps(table, cfg, k):
    for _ in range(k):
        cfg = one(table, cfg)
    return cfg

def padded(a, branch):
    need(a[branch] is not None, 'only ordinary instructions are padded')
    target, write, direction = a[branch]
    b = dict(a)
    b[branch] = (25, write, direction)
    for bit in (0,1):
        b[25, bit] = (26, bit, 1)
        b[26, bit] = (target, bit, -1)
    return b

def all_padding_guards(a):
    need(all(v is None or 0 <= v[0] < 25 for v in a.values()), 'fresh state labels do not collide with old targets')
    try:
        padded(a,(24,0))
    except ValueError:
        rejected_halt=True
    else:
        rejected_halt=False
    need(rejected_halt, 'reject padding an undefined transition')
    cases = unchanged = 0
    for branch, inst in a.items():
        if inst is None:
            continue
        b = padded(a, branch)
        for h in (-11, 0, 17):
            for mask in range(32):
                cells = {h+j-2 for j in range(5) if mask >> j & 1}
                cells.update((h-1000, h+1000))
                bit = int(h in cells)
                cfg = (branch[0], h, frozenset(cells))
                k = 3 if bit == branch[1] else 1
                need(steps(b, cfg, k) == one(a, cfg), 'padding restores full configuration')
                cases += int(k == 3)
                unchanged += int(k == 1)
        need(one(b, (24,0,frozenset())) is None, 'halt branch untouched')
    chosen = padded(a, (0,0))
    start = (0, 0, frozenset())
    need(one(chosen,start)[0] == 25 and steps(chosen,start,2)[0] == 26, 'fresh states reached on blank input')
    # Coupled original/padded run: time changes, every old boundary agrees.
    original = expanded = start
    replacements = elapsed = 0
    for _ in range(10000):
        q,h,c = original
        k = 3 if (q,int(h in c)) == (0,0) else 1
        replacements += int(k == 3)
        original = one(a,original)
        expanded = steps(chosen,expanded,k)
        elapsed += k
        need(original == expanded and original is not None, 'finite padded conjugacy at old boundaries')
    # Finite halting toy: one selected right-going step followed by a halt.
    toy = {(0,0):(1,1,1), (0,1):(1,0,-1), (1,0):None, (1,1):None}
    b = padded(toy,(0,0))
    need(steps(toy,start,2) is None and steps(b,start,4) is None, 'positive halting toy')
    need(steps(toy,start,1) is not None and steps(b,start,3) is not None, 'exact halting attempt clocks')
    return {'ordinary_branches_tested':49, 'guarded_local_configurations':cases, 'unselected_read_controls':unchanged, 'translations':[-11,0,17], 'far_tape_sentinels_preserved':True, 'undefined_selection_rejected':rejected_halt, 'blank_new_state_path':[25,26], 'coupled_old_steps':10000, 'padded_steps':elapsed, 'selected_instruction_visits':replacements, 'halting_toy_old_new_attempts':[2,4]}

def trace(a, w):
    c = (0,0,frozenset())
    d = (0,0,frozenset())
    rho = [ord(x)-65 for x in RENAMING]
    macro=[]
    seen=set()
    for t in range(30001):
        q,h,cells = c
        seen.add(q)
        need(d == (rho[q],h,cells), 'independent literal configurations conjugate')
        if q == 21 and cells and h not in cells and cells == frozenset(range(h-len(cells),h)):
            if not macro or macro[-1]['n'] != len(cells):
                macro.append({'n':len(cells),'step':t})
        if t == 30000:
            break
        c,d=one(a,c),one(w,d)
        need(c is not None and d is not None, 'finite trace does not halt')
    need(macro[:4] == [{'n':4,'step':45},{'n':6,'step':132},{'n':8,'step':301},{'n':10,'step':596}], 'macro clocks')
    original = json.loads((SOURCE/'checks.json').read_text())['trace']['macro_blocks']
    need(macro == original[:len(macro)], 'small independent macro prefix versus reported output')
    return {'microsteps':30000, 'visited_states':len(seen), 'macro_blocks':macro, 'infinite_nonhalting_conclusion':False}

def arithmetic():
    # Eratosthenes sieve, independent of source trial division.
    sieve=[True]*501
    sieve[0]=sieve[1]=False
    for p in range(2,23):
        if sieve[p]:
            for n in range(p*p,501,p):
                sieve[n]=False
    primes=[i for i in range(501) if sieve[i]]
    counts={n:0 for n in range(4,501,2)}
    for p in primes:
        for q in primes:
            if p+q in counts:
                counts[p+q]+=1
    for n,c in counts.items():
        need(c == sum(sieve[p] and sieve[n-p] for p in range(2,n-1)) and c > 0, 'finite pair coefficient')
    pairs=[(p,q) for p in primes if 3<=p<=200 for q in primes if 3<=q<=200]
    for p,q in pairs:
        E,O=p+q,p+2*q
        need((2*E-O,O-E)==(p,q) and q+2*p==3*E-O and ((p==q)==(2*O==3*E)), 'faithful pair interface')
    need({x for x,y in {(0,0),(1,1)}}=={x for x,y in {(0,1),(1,0)}} and not ({(0,0),(1,1)} & {(0,1),(1,0)}), 'marginal-loss hostile')
    return {'even_targets':len(counts),'largest_even':500,'odd_prime_pairs':len(pairs),'first_pair_counts':{n:c for n,c in counts.items() if n<=20}}

def source_static():
    root=SOURCE/'inputs/goldbach_tm'
    pending=['GoldbachTm.Tm25.Content']; visited=set(); external=set()
    while pending:
        name=pending.pop()
        if name in visited:
            continue
        visited.add(name)
        p=root/(name.replace('.','/')+'.lean')
        need(p.is_file(), 'complete saved local import cone: '+name)
        text=p.read_text()
        need(not re.search(r'\b(sorry|admit|axiom|sorryAx|native_decide)\b',text), 'static flagged token')
        for module in re.findall(r'^import\s+(\S+)',text,re.M):
            if module.startswith('GoldbachTm.'):
                pending.append(module)
            else:
                external.add(module)
    need('import GoldbachTm.Tm25.Content' in (root/'GoldbachTm.lean').read_text(), 'root imports endpoint')
    content=(root/'GoldbachTm/Tm25/Content.lean').read_text()
    need("theorem halt_lemma':" in content and 'theorem halt_lemma_rev :' in content, 'both reduction endpoints')
    need('never_none : ∀ i, nth_cfg i ≠ none' in content, 'never_halt has nonhalting hypothesis')
    return {'local_import_cone':sorted(visited), 'external_imports_not_built':sorted(external), 'root_endpoint_import':True, 'lean_build':'NOT RUN', 'endpoint_axiom_audit':'NOT RUN'}

def main():
    integrity=sealed_integrity()
    a,b=lean_table(),wiki_table()
    need(same_table(a,b),'all 50 slots relabel exactly')
    need(sum(v is not None for v in a.values())==49 and a[24,0] is None,'49 instructions, unique halt')
    need(one(a,(24,0,frozenset())) is None and one(a,(24,0,frozenset({0}))) is not None,'state alone loses halt semantics')
    mutation_failures=0
    for change in [(1,1,1),(22,0,1),(22,1,-1),None]:
        bad=dict(a);bad[0,0]=change
        need(not same_table(bad,b),'reject altered target/write/move/halt')
        mutation_failures+=1
    out={'status':'PASS FINITE-EXACT and static source audit only','integrity':integrity,'table':{'slots':50,'ordinary_instructions':49,'renaming':RENAMING,'hostile_mutations_rejected':mutation_failures},'trace':trace(a,b),'padding':all_padding_guards(a),'arithmetic':arithmetic(),'source_static':source_static()}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
