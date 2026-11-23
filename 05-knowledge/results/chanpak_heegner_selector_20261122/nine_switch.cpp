// Read the seven marked-apex path counts on the full switching orbit and its arc neighbours.
#include <array>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using Rows = std::array<int,9>;
using Chart = std::array<int,7>; // h,q,c,t,e,r,k

Chart chart(const Rows& rows,int apex) {
    int idx[8],m=0,tag=0,out[8]{};
    for(int i=0;i<9;i++) if(i!=apex) idx[m++]=i;
    for(int i=0;i<8;i++) {
        if(rows[apex]>>idx[i]&1) tag|=1<<i;
        for(int j=0;j<8;j++) if(rows[idx[i]]>>idx[j]&1) out[i]|=1<<j;
    }
    int start[256][8]{},all[256]{},end[256][8]{};
    all[0]=1;
    for(int i=0;i<8;i++) start[1<<i][i]=end[1<<i][i]=1;
    for(int s=1;s<256;s++) for(int b=s;b;b&=b-1) {
        int i=__builtin_ctz((unsigned)b),rest=s^(1<<i);
        for(int t=rest&out[i];t;t&=t-1) start[s][i]+=start[rest][__builtin_ctz((unsigned)t)];
        for(int t=rest&~out[i];t;t&=t-1) end[s][i]+=end[rest][__builtin_ctz((unsigned)t)];
        all[s]+=end[s][i];
    }
    Chart z{all[255],0,0,0,0,0,0};
    for(int s=0;s<256;s++) {
        int u=255^s,starts=0,ends=0;
        for(int t=s&tag;t;t&=t-1) starts+=start[s][__builtin_ctz((unsigned)t)];
        for(int t=s&~tag;t;t&=t-1) ends+=end[s][__builtin_ctz((unsigned)t)];
        int parity=(__builtin_popcount((unsigned)s)%2)?1:-1;
        z[1]+=all[s]*all[u];
        z[2]+=starts*all[u]; z[3]+=ends*all[u];
        z[4]+=parity*starts*all[u]; z[5]+=parity*ends*all[u];
        if(!s || s==255) z[6]+=starts+ends;
        else for(int t=u&~tag;t;t&=t-1) z[6]+=end[u][__builtin_ctz((unsigned)t)]*starts;
    }
    return z;
}
int main(){
    constexpr std::uint64_t code=0x7103d59dcULL;
    for(int cut=0;cut<256;cut++) {
        Rows rows{};
        int bit=0;
        for(int i=0;i<9;i++) for(int j=i+1;j<9;j++,bit++)
            if(((code>>bit)^(cut>>i)^(cut>>j))&1) rows[i]|=1<<j;
            else rows[j]|=1<<i;
        for(int a=0;a<9;a++) {
            Chart z=chart(rows,a);
            for(int i=0;i<9;i++) for(int j=i+1;j<9;j++) {
                Rows changed=rows;
                changed[i]^=1<<j; changed[j]^=1<<i;
                Chart zz=chart(changed,a);
                if(z[0]!=zz[0]) continue;
                std::cout<<cut<<" "<<a<<" "<<i<<" "<<j;
                for(int v:z) std::cout<<" "<<v;
                for(int k=0;k<7;k++) std::cout<<" "<<zz[k]-z[k];
                std::cout<<"\n";
            }
        }
    }
}
