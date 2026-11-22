#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

int n,d;
long long code(int n,int mask){
    int t=1,ans=0;
    for (int i=0;i<n;i++,t*=3) ans+=(mask>>i&1)*t;
    return ans;
}
int height(int n,int mask){
    std::array<int,6> in{};
    int j=0;
    for(int a=0;a<n;a++)for(int b=a+1;b<n;b++,j++){
        if(mask>>j&1)in[a]|=1<<b; else in[b]|=1<<a;
    }
    int dp[64][6]{};
    for(int i=0;i<n;i++)dp[1<<i][i]=1;
    for(int s=1;s<(1<<n);s++){
        for(int b=0;b<n;b++)if((s>>b&1)&&s!=(1<<b)){
            for(int a=0;a<n;a++)if((in[b]>>a&1))dp[s][b]+=dp[s^(1<<b)][a];
        }
    }
    int h=0; for(int b=0;b<n;b++)h+=dp[(1<<n)-1][b];
    return h;
}
int main(){
    for(n=3;n<=6;n++){
        d=n*(n-1)/2;
        std::vector<int> pow3={1};for(int i=0;i<d;i++)pow3.push_back(3*pow3.back());
        std::vector<uint16_t> weight(pow3.back()); std::vector<uint8_t> dim(pow3.back());
        int maxheight=0;
        for(int a=0;a<(1<<d);a++){
            int h=height(n,a);
            weight[code(d,a)]=h;
            maxheight=std::max(maxheight,h);
        }
        if((n==5 && maxheight!=15)||(n==6 && maxheight!=45))return 2;
        std::array<std::array<int,16>,721> g{};
        for(int code=0;code<pow3.back();code++){
            int q=code,j=0;
            for(;q&&q%3!=2;q/=3,j++);
            if(q){int i=pow3[j];weight[code]=std::max(weight[code-2*i],weight[code-i]);dim[code]=dim[code-i]+1;}
            g[weight[code]][dim[code]]++;
        }
        for(int L:{1,3,5,7,9,19,21,23,45}){
            std::cout<<n<<" cap "<<L<<":";
            int chi=0;
            for(int i=0;i<=d;i++){
                int sum=0;
                for(int j=1;j<=L;j++)sum+=g[j][i];
                std::cout<<" "<<sum;
                chi+=i%2?-sum:sum;
            }
            std::cout<<" chi "<<chi<<"\n";
        }
    }
    return 0;
}
