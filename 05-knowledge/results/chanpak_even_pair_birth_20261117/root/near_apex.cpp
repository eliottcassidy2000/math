#include<bits/stdc++.h>
using namespace std;
using I=long long;
using V=vector<int>;
constexpr int N=13,U=1<<N;
I d[U][N],h[U],W[N][N],cut[U];
int select_mask(const V&a,I target){
    int n=a.size(),M=1<<n;
    auto dp=[&](){
        for(int s=1;s<M;s++)
            for(int b=s;b;b&=b-1){
                int i=__builtin_ctz(b); I z=d[s][i];if(!z)continue;
                for(int c=a[i]&~s;c;c&=c-1){int x=c&-c;d[s|x][__builtin_ctz(x)]+=z;}
            }
    };
    memset(d,0,M*sizeof(d[0]));
    for(int i=0;i<n;i++)d[1<<i][i]=1;dp();h[0]=1;
    for(int s=1;s<M;s++)h[s]=accumulate(d[s],d[s]+n,0LL);
    I r=target-h[M-1];
    if(r<0 || r%2)return -1;r/=2;
    memset(W,0,sizeof(W));
    for(int i=0;i<n;i++){
        memset(d,0,M*sizeof(d[0]));d[1<<i][i]=1;dp();
        for(int s=1;s<M;s++)if((s>>i&1)&&!(__builtin_popcount(s)&1)){
            I rem=h[M-1-s];
            for(int b=s;b;b&=b-1){int j=__builtin_ctz(b);W[i][j]+=d[s][j]*rem;}
        }
    }
    I out[N]={};
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)out[i]+=W[i][j];
    cut[0]=0;
    for(int s=1;s<M;s++){
        int i=__builtin_ctz(s),p=s&(s-1);I v=cut[p]+out[i];
        for(int b=p;b;b&=b-1){int j=__builtin_ctz(b);v-=W[i][j]+W[j][i];}
        cut[s]=v;
        if(v==r)return s;
    }
    return r==0?0:-1;
}
int main(int argc,char**argv){
    if(argc!=4)return 1;
    int n=atoi(argv[1]),attempts=atoi(argv[2]);I t=atoll(argv[3]);
    if(n<3||n>N)return 2;
    mt19937_64 rng(97213);
    for(int it=0;it<attempts;it++){
        V a(n);
        for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
            int v=j-i;bool forward=(it%2?(v==1||v==3||v==4||v==5||v==9):(2*v<n));
            if(forward)a[i]|=1<<j;else a[j]|=1<<i;
        }
        int flips=(it/24)%16;
        for(int f=0;f<flips;f++){
            int i=rng()%n,j=rng()%n;if(i==j)continue;
            a[i]^=1<<j;a[j]^=1<<i;
        }
        int q=select_mask(a,t);
        if(q>=0){
            cout<<n+1<<" "<<t<<" "<<it<<" "<<q<<" deletion ";
            for(int z:a)cout<<z<<",";
            cout<<" literal ";
            for(int i=0;i<n;i++)cout<<(a[i]|((1-(q>>i&1))<<n))<<",";
            cout<<q<<"\n";return 0;
        }
    }
    cout<<"SAMPLED MISS "<<n+1<<" "<<attempts<<" target "<<t<<"\n";
}
