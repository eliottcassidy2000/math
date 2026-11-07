#include <bits/stdc++.h>
using namespace std;using I=long long;using Key=array<I,9>;struct Rec{int ix,mask;I s;};
int main(int argc,char**argv){ifstream in(argv[1]);int total;in>>total;map<Key,array<Rec,2>> bank;I orientations=0,retained=0,hits=0;
for(int ix=0;ix<total;ix++){int n;in>>n;vector<int>a(n);for(auto&v:a)in>>v;int N=1<<n,F=N-1;vector<I>en(N*n),st(N*n),H(N),C(N);H[0]=1;
for(int i=0;i<n;i++)en[(1<<i)*n+i]=st[(1<<i)*n+i]=1;
for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){int A=S^(1<<v);for(int w=0;w<n;w++)if(A>>w&1){if(a[w]>>v&1)en[S*n+v]+=en[A*n+w];else st[S*n+v]+=st[A*n+w];}H[S]+=en[S*n+v];}
for(int S=0;S<N;S++){int A=S;while(true){C[S]+=H[A]*H[S^A];if(!A)break;A=(A-1)&S;}}
vector<I>EL(n),SR(n),ST(n),QH(n*n),QC(n*n);
for(int A=1;A<N;A++){int rem=F^A;for(int u=0;u<n;u++){EL[u]+=en[A*n+u]*H[rem];SR[u]+=st[A*n+u]*H[rem];ST[u]+=st[A*n+u]*C[rem];}
for(int B=rem;B;B=(B-1)&rem){I h=H[rem^B];for(int u=0;u<n;u++)if(en[A*n+u])for(int v=0;v<n;v++)if(st[B*n+v]){I x=en[A*n+u]*st[B*n+v];QC[u*n+v]+=x*h;if(B==rem)QH[u*n+v]+=x;}}}
for(int mask=0;mask<N;mask++){orientations++;I s=0,e=0,L=H[F],R=H[F],T=C[F],h=0,c=H[F];for(int u=0;u<n;u++)if(mask>>u&1){s+=st[F*n+u];R+=SR[u];T+=ST[u];c+=SR[u];}else{e+=en[F*n+u];L+=EL[u];c+=EL[u];}
if(s<=0||s>26||s%3==0)continue;retained++;h=s+e;for(int u=0;u<n;u++)if(!(mask>>u&1))for(int v=0;v<n;v++)if(mask>>v&1){h+=QH[u*n+v];c+=QC[u*n+v];}
cout<<ix<<" "<<mask<<" "<<n+1<<" "<<h<<" "<<s<<" "<<e<<" "<<L<<" "<<R<<" "<<H[F]<<" "<<2*c<<" "<<C[F]<<" "<<T<<"\n";

}
}
cerr<<"orientations "<<orientations<<" retained "<<retained<<" keys "<<bank.size()<<" hits "<<hits<<"\n";
}
