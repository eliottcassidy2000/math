#define main declared_family_main
#include "twelve.cpp"
#undef main
int main(int argc,char**argv) {try {
    if(argc!=2)throw std::runtime_error("control input");
    std::ifstream f(argv[1]);int count;f>>count;
    std::cout<<"[";
    for(int i=0;i<count;++i) {
        int n;f>>n;if(n<2||n>12||(n&1))throw std::runtime_error("even control order");
        Rows a{};for(int j=0;j<n;++j)f>>a[j];if(!f)throw std::runtime_error("truncated input");
        DP e=endpoints(a,n),s=endpoints(reverse(a,n),n);Rows c=cubic_columns(s,n);
        if(i)std::cout<<",";
        std::cout<<"["<<e[(1<<n)-1]<<","<<s[(1<<n)-1]<<","<<cover_gate(a,n);
        for(int j=0;j<n;++j)std::cout<<","<<c[j];std::cout<<"]";
    }
    std::cout<<"]\n";
} catch(const std::exception&e) {std::cerr<<e.what()<<"\n";return 1;}}
