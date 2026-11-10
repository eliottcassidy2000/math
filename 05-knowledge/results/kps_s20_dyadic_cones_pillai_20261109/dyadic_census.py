# Dyadic-port census from the S19 port statistics (all classes of order <= 8, all ports).
# Input line: n code u h2 E2 O2 hu2 h3 E3 O3 hu3 ; base-2 triple (h, e, R) = (h2, O2 - E2, E2 + O2).
import sys, collections, math
odd=lambda r: r>>((r&-r).bit_length()-1)
ispow=lambda r: r>0 and r&(r-1)==0
oddparts=collections.defaultdict(set); best={}; vals=set(); triples=set()
for line in open(sys.argv[1]):
    p=list(map(int,line.split())); n,code,u,h,E,O,hu=p[:7]; R=E+O; e=O-E
    vals.add(h); oddparts[h].add(odd(R))
    if ispow(R):
        triples.add((h,e,R)); j=R.bit_length()-1
        if h not in best or j<best[h][0]: best[h]=(j,n,code,u)
small=[v for v in sorted(vals) if v<=300]
print('attained odd values <= 300 (order <= 8):',len(small),' with a dyadic port:',sum(1 for v in small if v in best))
print('distinct dyadic triples (h,e,R):',len(triples))
for h in (5,13):
    print('response odd parts at H = %d:'%h, sorted(oddparts[h]))
print('values <= 300 with a dyadic port and least exponent j:')
print(' '.join('%d:%d'%(h,best[h][0]) for h in sorted(best) if h<=300))
