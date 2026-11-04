#include <array>
#include <iostream>
#include <vector>
using namespace std;
// Over F2, determinants count directed cycle covers. Isolate the unique cycle
// containing vertex0 by subset convolution; no path/endpoint recurrence is used.
int determinant(const vector<int>& rows,int mask){
 int a[11],n=0;for(unsigned bits=mask;bits;bits&=bits-1)a[n++]=rows[__builtin_ctz(bits)]&mask;
 int pivot=0;for(unsigned bits=mask;bits;bits&=bits-1){int b=bits&-bits,k=pivot;while(k<n && !(a[k]&b))k++;if(k==n)return 0;swap(a[k],a[pivot]);for(int i=pivot+1;i<n;i++)if(a[i]&b)a[i]^=a[pivot];pivot++;}return 1;
}
int cycles(const vector<int>&r){
 int N=1<<r.size();vector<unsigned char>cover(N),one(N);cover[0]=1;
 for(int s=1;s<N;s++)cover[s]=determinant(r,s);
 for(int s=1;s<N;s+=2){int value=cover[s];for(int a=(s-1)&s;a;a=(a-1)&s)if(a&1)value^=one[a]&cover[s^a];one[s]=value;}
 return one[N-1];
}
int main(){int cases;cin>>cases;for(int z=0;z<cases;z++){int node,a,c;cin>>node>>a>>c;vector<int>q(9),r(11);for(int &v:q)cin>>v;cout<<node<<'\t'<<a<<'\t'<<c;for(int op=0;op<2;op++){for(int i=0;i<9;i++)r[i]=(op?(511^q[i]^(1<<i)):q[i])|(((a>>i)&1)<<9)|(((c>>i)&1)<<10);r[9]=(511^a)|1024;r[10]=511^c;cout<<'\t'<<cycles(r);}cout<<'\n';if(z%1000==0)cerr<<z+1<<'/'<<cases<<'\n';}}
