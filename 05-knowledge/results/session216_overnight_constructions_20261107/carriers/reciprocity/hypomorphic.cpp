#include <bits/stdc++.h>
using namespace std;using I=long long;
int main(int argc,char**argv){int maxk=argc>1?atoi(argv[1]):8;I count=0,first=0,second=0,hits=0;
for(int k=1;k<=maxk;k++){int n=k+4,N=1<<n,F=N-1;for(int dist=0;dist<(1<<k);dist++)for(int bits=1;bits<7;bits++){
vector<int>a(n);for(int i=0;i<=k;i++)for(int j=i+1;j<=k;j++)if(dist>>(j-i-1)&1)a[i]|=1<<j;else a[j]|=1<<i;
for(int i=0;i<3;i++)a[k+1+i]|=1<<(k+1+(i+1)%3);
for(int i=0;i<=k;i++)for(int j=0;j<3;j++)if(bits>>((j+i)%3)&1)a[i]|=1<<(k+1+j);else a[k+1+j]|=1<<i;
vector<I>en(N*n),st(N*n),H(N);H[0]=1;for(int i=0;i<n;i++)en[(1<<i)*n+i]=st[(1<<i)*n+i]=1;
for(int S=1;S<N;S++)for(int v=0;v<n;v++)if(S>>v&1){int A=S^(1<<v);int old=A;while(old){int w=__builtin_ctz((unsigned)old);old&=old-1;if(a[w]>>v&1)en[S*n+v]+=en[A*n+w];else st[S*n+v]+=st[A*n+w];}H[S]+=en[S*n+v];}
count++;I s0=st[F*n],sk=st[F*n+k];if(s0<1||sk<1||s0>26||sk>26||abs(s0-sk)!=2||en[F*n]!=en[F*n+k])continue;
I hi=max(s0,sk),lo=min(s0,sk);if(hi%3!=1||lo%3!=2)continue;first++;
I l0=0,lk=0,r0=0,rk=0;for(int S=1;S<N;S++){l0+=en[S*n]*H[F^S];lk+=en[S*n+k]*H[F^S];r0+=st[S*n]*H[F^S];rk+=st[S*n+k]*H[F^S];}
cout<<"CANDIDATE "<<k<<" "<<dist<<" "<<bits<<" starts "<<s0<<" "<<sk<<" end "<<en[F*n]<<" L "<<l0<<" "<<lk<<" R "<<r0<<" "<<rk<<"\n"<<flush;
if(l0!=lk||r0!=rk)continue;second++;
vector<I>C(N);for(int S=0;S<N;S++){int A=S;while(true){C[S]+=H[A]*H[S^A];if(!A)break;A=(A-1)&S;}}
I t0=0,tk=0;for(int S=1;S<N;S++){t0+=st[S*n]*C[F^S];tk+=st[S*n+k]*C[F^S];}if(t0!=tk)continue;hits++;cout<<"UNIVERSAL";for(auto row:a)cout<<" "<<row;cout<<"\n"<<flush;
}cerr<<"k "<<k<<" graphs "<<count<<" first "<<first<<" second "<<second<<" hits "<<hits<<"\n";}
}
