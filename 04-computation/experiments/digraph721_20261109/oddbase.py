# least order of an ODD realiser for each odd k < 512, from the exhaustive census of connected graphs with d <= 511
# (all orders).  A graph is odd iff some component is odd; d is multiplicative; any value m has a realiser on L(m)
# vertices (threshold), and 1 needs no vertices (empty union).
INF=99
minconn_odd={}; 
for line in open('dom_K511.out'):
    p=line.split(); n=int(p[0]); d=int(p[2])
    if p[3]=='O': minconn_odd[d]=min(minconn_odd.get(d,INF),n)
def L(k): return k.bit_length()
def Lany(m): return 0 if m==1 else L(m)
minodd={}
for k in range(1,512,2):
    best=INF
    for a,o in minconn_odd.items():
        if k%a==0: best=min(best,o+Lany(k//a))
    minodd[k]=best
miss=[k for k in minodd if minodd[k]>=INF]
need1=[k for k in minodd if minodd[k]==L(k)+1]
over=[k for k in minodd if L(k)+1<minodd[k]<INF]
print('no odd realiser:',miss)
print('least order L(k)+1 (not L(k)):',need1)
print('least order > L(k)+1:',over)
