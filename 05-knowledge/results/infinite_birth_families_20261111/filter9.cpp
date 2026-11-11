#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
using namespace std;
int main(int argc,char**argv){
    if(argc!=3)return 2;
    ifstream in(argv[1]); ofstream out(argv[2]);
    uint64_t code; int aut; long long tested=0,kept=0;
    while(in>>code>>aut){
        array<int,9> rows{}; int b=35;
        for(int i=0;i<9;i++)for(int j=i+1;j<9;j++,b--)
            if(code>>b&1)rows[i]|=1<<j;else rows[j]|=1<<i;
        array<int,9> reach=rows;
        for(int i=0;i<9;i++)reach[i]|=1<<i;
        for(int k=0;k<9;k++)for(int i=0;i<9;i++)if(reach[i]>>k&1)reach[i]|=reach[k];
        tested++;bool strong=true;for(int x:reach)strong &= x==511;
        if(!strong)continue;
        array<array<int,9>,512> dp{};
        for(int i=0;i<9;i++)dp[1<<i][i]=1;
        for(int s=1;s<512;s++)for(int i=0;i<9;i++)if(s>>i&1){
            int next=rows[i] & (511^s);
            while(next){int bit=next&-next,j=__builtin_ctz(next);next^=bit;dp[s|bit][j]+=dp[s][i];}
        }
        int h=0;for(int x:dp[511])h+=x;
        if(h<125){out<<code<<" "<<aut<<" "<<h;for(int x:rows)out<<" "<<x;out<<"\n";kept++;}
    }
    cout<<tested<<" "<<kept<<"\n";
    return tested==191536?0:3;
}
