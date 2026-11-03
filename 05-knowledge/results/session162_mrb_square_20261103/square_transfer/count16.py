def H(n,back):
 out=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if (j,i) in back:out[j]|=1<<i
   else:out[i]|=1<<j
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for m in range(1,1<<n):
  for j in range(n):
   if dp[m][j]:
    x=out[j]&~m
    while x:
     bit=x&-x;dp[m|bit][bit.bit_length()-1]+=dp[m][j];x-=bit
 return sum(dp[-1])
for j,want in enumerate((6435,9945,13455)):
 back={(5,i) for i in range(5)}|{(6,5),(8,6),(15,8)}
 if j>=1:back.add((2,0))
 if j>=2:back.add((3,0))
 actual=H(16,back)
 if actual!=want:raise AssertionError((j,want,actual))
 print(j,actual)
