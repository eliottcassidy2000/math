# emit one-centre tournaments TT_m + hub (vertex m) in gentourng upper-triangle format
import sys, itertools
M=int(sys.argv[1])
for m in range(1,M+1):
    n=m+1
    for x in itertools.product('01',repeat=m):
        if m>1 and x[0]=='0': continue   # leading zeros redundant
        bits=[]
        for i in range(n):
            for j in range(i+1,n):
                if j<m: bits.append('1')            # TT arc i->j
                else: bits.append(x[i])             # j==m is hub: bit 1 means i->hub ; x_i=1 means hub->i
        # hub->i iff x_i=1, so i->hub iff x_i=0
        bits=[ (('1' if x[i]=='0' else '0') if j==m else '1') for i in range(n) for j in range(i+1,n)]
        print(''.join(bits))
