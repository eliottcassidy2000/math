import re, sys
for fn in sys.argv[1:]:
    L = open(fn).read().splitlines()
    print('==', fn)
    print(' oct | E_V/N mean max | G0/N mean | G1/G0 mean [lo,hi] | (G0-E0)/N mean max | Gam1/N mean max | HL/N mean max | E+/N mean | E_V/E+ ')
    for i, l in enumerate(L):
        m = re.match(r'octave 2\^\s*(\d+):', l)
        if not m: continue
        k = int(m.group(1))
        ev = re.search(r'E_V max ([\d.]+) mean ([\d.]+)', l); g0 = re.search(r'E_clk max ([\d.]+) mean ([\d.]+)', l)
        r1 = re.search(r'G1/G0 mean ([+-][\d.]+) range \[([+-][\d.]+), ([+-][\d.]+)\]', l)
        hl = re.search(r'higher lags \(E_V-G0-2G1\)/N mean ([+-][\d.]+) range \[([+-][\d.]+), ([+-][\d.]+)\]', l)
        ep = re.search(r'unsigned E\+/N mean ([\d.]+).*dip E_V/E\+ mean ([\d.]+)', L[i+1])
        lg = re.search(r'same-clock \(G0-E0\)/N mean ([+-][\d.]+) max ([+-][\d.]+) ; other lag-1 Gam1/N = \(G1\+SIB\)/N mean ([+-][\d.]+) max ([+-][\d.]+)', L[i+2])
        print(' %2d | %.3f %.3f | %.3f | %s [%s, %s] | %s %s | %s %s | %s %s | %s | %s' % (k, float(ev.group(2)), float(ev.group(1)), float(g0.group(2)),
              r1.group(1), r1.group(2), r1.group(3), lg.group(1), lg.group(2), lg.group(3), lg.group(4), hl.group(1), hl.group(3), ep.group(1), ep.group(2)))
