import sys
pairs=set()
for n in range(3,10):
    for line in open(f"hl{n}.txt"):
        parts=line.split(); H=int(parts[1])
        for item in parts[2:]:
            v,rest=item.split(':'); h,L,E,S=map(int,rest.split(','))
            pairs.add((H,L))
for t in range(0,40): pairs.add((1,2**t))
LIM=int(sys.argv[1])
cov=bytearray(LIM//2+2)
for (H,L) in pairs:
    j=1
    while H+(L<<j)<=LIM:
        cov[(H+(L<<j))>>1]=1; j+=1
miss=[N for N in range(1,LIM+1,2) if not cov[N>>1]]
print("odd N<=%d that are not cuts h+2^j L of rooted tournaments of order<=9: %d"%(LIM,len(miss)))
print("first misses:",miss[:25])
big=[N for N in miss if N>67]
print("first miss above 67:", big[:1], " count of misses above 67:",len(big))
