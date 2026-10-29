import re, glob, collections
rows = []
for f in sorted(glob.glob('out/census_n*_k*.out')):
    n, k = map(int, re.findall(r'n(\d+)_k(\d+)', f)[0])
    lines = [l for l in open(f) if re.match(r'^[01]{3,}\s', l)]
    c = collections.Counter()
    for l in lines:
        strong = ' strong ' in l
        c['classes'] += 1
        c['strong'] += strong
        m = re.search(r' M=(.*?)( unglued=| P=|$)', l)
        if m:
            snf = m.group(1).strip()
            uni = all(t.startswith('1^') for t in snf.split())
            sing = '0^' in snf
            c['unimodular'] += uni
            c['singular'] += sing
            c['strong_unimodular'] += uni and strong
            c['strong_singular'] += sing and strong
        else:
            m2 = re.search(r'rankM=(\d+)/(\d+)', l)
            if m2:
                c['singular'] += int(m2.group(1)) < int(m2.group(2))
                c['strong_singular'] += strong and int(m2.group(1)) < int(m2.group(2))
        sp = re.search(r'spur2=(\d+) spur3=(\d+)', l)
        if sp:
            c['spurious'] += int(sp.group(1)) + int(sp.group(2)) > 0
            c['spur2>0'] += int(sp.group(1)) > 0
            c['spur3>0'] += int(sp.group(2)) > 0
        c['killsBdry'] += 'killsBdry=True' in l
    rows.append((n, k, c))
for n, k, c in rows:
    print('n=%d k=%d: %s' % (n, k, ', '.join('%s=%d' % kv for kv in c.items())))
