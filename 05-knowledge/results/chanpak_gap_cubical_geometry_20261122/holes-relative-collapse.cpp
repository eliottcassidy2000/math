#define main cubical_census_main
#include "holes-cubical.cpp"
#undef main

bool union_strong(int code) {
    int out[6]={},in[6]={};
    for(int i=0;i<6;i++)for(int j=i+1;j<6;j++){
        int t=code%3;code/=3;
        if(t!=1){out[i]|=1<<j;in[j]|=1<<i;}
        if(t!=0){out[j]|=1<<i;in[i]|=1<<j;}
    }
    for(auto arr:{out,in}){
        int s=1,old=0;
        while(s!=old){old=s;for(int i=0;i<6;i++)if((s>>i)&1)s|=arr[i];}
        if(s!=63)return false;
    }
    return true;
}
int main(){
    constexpr int D=15;
    std::vector<int> p={1};for(int j=0;j<D;j++)p.push_back(3*p.back());
    std::vector<uint16_t> height_max(p.back()); std::vector<uint8_t> dim(p.back()),active(p.back()),deg(p.back());
    for(int v=0;v<(1<<D);v++)height_max[code(D,v)]=height(6,v);
    for(int c=0;c<p.back();c++){
        int z=c,j=0;
        while(z&&z%3!=2){z/=3;j++;}
        if(z){int w=p[j];height_max[c]=std::max(height_max[c-2*w],height_max[c-w]);dim[c]=dim[c-w]+1;}
    }
    for(int L:{15,17,19,21,23}){
        std::fill(active.begin(),active.end(),0);
        std::fill(deg.begin(),deg.end(),0);
        int tot=0,k=0;
        std::array<int,D+1> hist0{},hist1{};
        std::vector<int> pending;
        for(int c=0;c<p.back();c++)if(height_max[c]<=L&&union_strong(c)){
            active[c]=1;tot++;hist0[dim[c]]++;
        }
        auto cofaces=[&](int c,auto fn){
            for(int j=0,z=c;j<D;j++,z/=3)if(z%3<2&&active[c+(2-z%3)*p[j]])fn(c+(2-z%3)*p[j]);
        };
        for(int c=0;c<p.back();c++)if(active[c]){
            cofaces(c,[&](int){deg[c]++;});
            if(deg[c]==1)pending.push_back(c);
        }
        for(size_t i=0;i<pending.size();i++){
            int c=pending[i];
            if(!active[c]||deg[c]!=1)continue;
            int c0=-1;cofaces(c,[&](int g){c0=g;});if(c0<0)return 3;
            if(deg[c0])return 4;
            active[c]=active[c0]=0;k++;
            for(int g:{c,c0})for(int j=0,z=g;j<D;j++,z/=3)if(z%3==2)
                for(int t:{g-p[j],g-2*p[j]})if(active[t]&&--deg[t]==1)pending.push_back(t);
        }
        for(int c=0;c<p.back();c++)if(active[c])hist1[dim[c]]++;
        if(tot!=2*k)return 5;
        std::cout<<"cap "<<L<<" relative "<<tot<<" collapses "<<k<<" left "<<tot-2*k<<"\npre ";
        for(int a:hist0)std::cout<<a<<" ";
        std::cout<<"\nleft ";for(int a:hist1)std::cout<<a<<" ";
        std::cout<<"\n";
    }
}
