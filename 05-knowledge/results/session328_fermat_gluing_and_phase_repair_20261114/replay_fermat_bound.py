from profile_tools import canon, extend, profile
T={(0,)}
for n in range(2,9):
    T={canon(extend(A,c)) for A in T for c in range(1<<(n-1))}
M=max(profile(A)[1] for A in T)
if len(T)!=6880 or M!=8848:
    raise RuntimeError('the full order-eight Q2 census failed')
print('6880 order-eight classes; max Q2=8848; j<=2 cannot exceed 44240')
