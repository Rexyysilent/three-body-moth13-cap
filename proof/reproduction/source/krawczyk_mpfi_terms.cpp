#include <boost/multiprecision/mpfi.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
namespace bm=boost::multiprecision;
using R=bm::number<bm::backends::mpfr_float_backend<80>>;
using H=bm::number<bm::backends::mpfi_float_backend<80>>;
constexpr int Q=8;
using VH=std::array<H,Q>;
using MH=std::array<std::array<H,Q>,Q>;
using DBH=std::array<std::array<H,Q>,3>;
static H rangeH(const R&a,const R&b){H x;mpfi_interv_fr(x.backend().data(),a.backend().data(),b.backend().data());return x;}
static R lo(const H&x){return bm::lower(x);}static R hi(const H&x){return bm::upper(x);}
static R mag(const H&x){return std::max(R(abs(lo(x))),R(abs(hi(x))));}
static H sym(const R&r){return rangeH(-r,r);}
static R point(const std::string&s,bool ld){return ld?R(std::strtold(s.c_str(),nullptr)):R(s);}
static H readH(std::istream&f,bool ld){std::string a,b;f>>a>>b;return rangeH(point(a,ld),point(b,ld));}
struct Mesh{int n;R tau,rtau,rstate;std::vector<R>den;std::vector<VH> node;std::vector<std::array<std::string,Q>> rawNode;std::string rawTau;};
static Mesh readMesh(const char*p){
 std::ifstream f(p);if(!f)throw std::runtime_error("mesh open");Mesh m;std::string s;
 f>>m.n>>m.rawTau;m.tau=R(m.rawTau);f>>s;m.rtau=R(s);f>>s;m.rstate=R(s);
 m.den.resize(m.n);for(auto&x:m.den){f>>s;x=R(s);}
 m.node.resize(m.n+1);m.rawNode.resize(m.n+1);
 for(int k=0;k<=m.n;k++)for(int j=0;j<Q;j++){f>>s;m.rawNode[k][j]=s;m.node[k][j]=H(R(s));}
 return m;
}
struct Rec{VH xc{},xb{};MH M{};R mind2=0;};
struct Log{std::vector<Rec>r;R globalD2=0;};
static Log readLog(const char*p,int n,bool ld){
 std::ifstream f(p);if(!f)throw std::runtime_error("log open");Log L;L.r.resize(n);std::string t,s;
 for(int k=0;k<n;k++){
  int idx,ok;f>>t>>idx>>t>>ok;if(!f||idx!=k||ok!=1)throw std::runtime_error("bad REC");
  f>>t;f>>s>>s;f>>s;L.r[k].mind2=point(s,ld);f>>s>>s;
  f>>t;for(auto&x:L.r[k].xc)x=readH(f,ld);
  f>>t;for(auto&x:L.r[k].xb)x=readH(f,ld);
  f>>t;for(auto&row:L.r[k].M)for(auto&x:row)x=readH(f,ld);
 }
 int ok,count;f>>t>>ok>>t>>count>>t>>s;if(!f||ok!=1||count!=n)throw std::runtime_error("bad footer");
 L.globalD2=point(s,ld);
 for(const auto& rec:L.r){
  if(!(rec.mind2>0))throw std::runtime_error("nonpositive segment collision bound");
  if(!(L.globalD2<=rec.mind2))throw std::runtime_error("footer collision bound exceeds segment bound");
 }
 return L;
}
static VH field(const VH&x){
 H mu=H(R(1001)/1000),q3x=-(x[0]+x[2])/mu,q3y=-(x[1]+x[3])/mu;
 auto g=[](H a,H b){H r2=a*a+b*b;if(!(lo(r2)>0))throw std::runtime_error("field collision");H r3=r2*sqrt(r2);return std::array<H,2>{a/r3,b/r3};};
 auto a=g(x[2]-x[0],x[3]-x[1]),b=g(q3x-x[0],q3y-x[1]);
 auto c=g(x[0]-x[2],x[1]-x[3]),d=g(q3x-x[2],q3y-x[3]);
 return {x[4],x[5],x[6],x[7],a[0]+mu*b[0],a[1]+mu*b[1],c[0]+mu*d[0],c[1]+mu*d[1]};
}
static std::array<H,3> boundary(const VH&x){return {x[0]+x[2],x[1]+x[3],x[0]*(x[4]-x[6])+x[1]*(x[5]-x[7])};}
static DBH boundaryD(const VH&x){
 DBH d{};for(auto&r:d)for(auto&v:r)v=H(0);d[0][0]=1;d[0][2]=1;d[1][1]=1;d[1][3]=1;
 d[2][0]=x[4]-x[6];d[2][1]=x[5]-x[7];d[2][4]=x[0];d[2][5]=x[1];d[2][6]=-x[0];d[2][7]=-x[1];return d;
}
static std::vector<H> denseApply(const std::vector<long double>&C,int D,const std::vector<H>&v){
 std::vector<H>o(D,H(0));for(int i=0;i<D;i++){H s(0);auto row=&C[(std::size_t)i*D];for(int j=0;j<D;j++)if(!(lo(v[j])==0&&hi(v[j])==0))s+=H(R(row[j]))*v[j];o[i]=s;}return o;
}
int main(int ac,char**av){
 try{
  if(ac!=8){std::cerr<<"mesh mpfi-center.log ld-box.log C.bin point.bin R.bin cert\n";return 2;}
  Mesh mesh=readMesh(av[1]);int N=mesh.n,D=Q*N+3,bb=Q*N;
  int n200=0,n3200=0,n25600=0;bool denOK=true;
  for(const auto& d:mesh.den){if(d==R(200))++n200;else if(d==R(3200))++n3200;else if(d==R(25600))++n25600;else denOK=false;}
  bool meshTime=(N==244&&denOK&&n200==198&&n3200==30&&n25600==16);
  bool initialSlice=(mesh.node[0][0]==H(-1)&&mesh.node[0][1]==H(0)&&mesh.node[0][2]==H(1)&&mesh.node[0][3]==H(0)&&mesh.node[0][4]==mesh.node[0][6]&&mesh.node[0][5]==mesh.node[0][7]);
  bool structural=meshTime&&initialSlice&&(mesh.tau>mesh.rtau);
  Log cen=readLog(av[2],N,false),box=readLog(av[3],N,true);
  std::ifstream fc(av[4],std::ios::binary);int dc;fc.read((char*)&dc,sizeof(dc));if(dc!=D)throw std::runtime_error("C dim");
  std::vector<long double>C((std::size_t)D*D);fc.read((char*)C.data(),sizeof(long double)*C.size());if(!fc)throw std::runtime_error("C read");
  std::ifstream fp(av[5],std::ios::binary);int np;fp.read((char*)&np,sizeof(np));if(np!=N)throw std::runtime_error("point dim");
  std::vector<std::array<std::array<long double,Q>,Q>>M0(N);std::vector<std::array<long double,Q>>c0(N);std::array<std::array<long double,Q>,3>db0{};
  for(int k=0;k<N;k++){fp.read((char*)M0[k].data(),sizeof(long double)*Q*Q);fp.read((char*)c0[k].data(),sizeof(long double)*Q);}fp.read((char*)db0.data(),sizeof(long double)*3*Q);if(!fp)throw std::runtime_error("point read");
  std::ifstream fr(av[6],std::ios::binary);int dr;long double qRld;fr.read((char*)&dr,sizeof(dr));fr.read((char*)&qRld,sizeof(qRld));if(dr!=D)throw std::runtime_error("R dim");
  std::vector<long double>Rradld(D);fr.read((char*)Rradld.data(),sizeof(long double)*D);if(!fr)throw std::runtime_error("R read");
  std::ofstream out(av[7]);if(!out)throw std::runtime_error("cert open");out<<std::scientific<<std::setprecision(70);
  H tenth=H(1)/H(10),root2=sqrt(H(2));bool mpfi=(lo(tenth)<hi(tenth)&&lo(root2)<hi(root2));
  out<<"MPFI_OUTWARD_OK "<<mpfi<<" TENTH "<<lo(tenth)<<" "<<hi(tenth)<<" SQRT2 "<<lo(root2)<<" "<<hi(root2)<<"\n";
  R outerState=R(std::strtold("1.001e-14",nullptr)),outerTau=R(std::strtold("1.001e-13",nullptr));
  R ldTau=R(std::strtold(mesh.rawTau.c_str(),nullptr));bool domain=abs(mesh.tau-ldTau)+mesh.rtau<=outerTau;R maxOffset=abs(mesh.tau-ldTau);
  for(int k=0;k<N;k++)for(int j=0;j<Q;j++){R ld=R(std::strtold(mesh.rawNode[k][j].c_str(),nullptr)),off=abs(R(mesh.rawNode[k][j])-ld);maxOffset=std::max(maxOffset,off);domain=domain&&(off+mesh.rstate<=outerState);}
  out<<"MESH_TIME_SUM_EXACT "<<meshTime<<" INITIAL_SLICE_EXACT "<<initialSlice<<" STRUCTURAL_CHECKS "<<structural<<"\n";
  out<<"DOMAIN_CONTAINMENT "<<domain<<" MAX_CENTER_OFFSET "<<maxOffset<<"\n";
  out<<"COLLISION_D2_LOWER "<<box.globalD2<<"\n";
  std::vector<H>F(D,H(0));R maxF=0;int maxFIndex=-1;
  for(int k=0;k<N;k++)for(int i=0;i<Q;i++){F[Q*k+i]=cen.r[k].xc[i]-mesh.node[k+1][i];if(mag(F[Q*k+i])>maxF){maxF=mag(F[Q*k+i]);maxFIndex=Q*k+i;}}
  auto bf=boundary(mesh.node[N]);for(int i=0;i<3;i++){F[bb+i]=bf[i];if(mag(bf[i])>maxF){maxF=mag(bf[i]);maxFIndex=bb+i;}}
  auto CF=denseApply(C,D,F);
  std::vector<R>rad(D,mesh.rstate);rad[2]=mesh.rtau;
  std::vector<R>Dmag(D,R(0));
  for(int k=0;k<N;k++){
   auto ci=field(box.r[k].xb);for(int i=0;i<Q;i++)ci[i]/=H(mesh.den[k]);
   for(int i=0;i<Q;i++){H s(0);
    if(k==0){
     H ja=box.r[k].M[i][4]+box.r[k].M[i][6],jb=box.r[k].M[i][5]+box.r[k].M[i][7];
     long double ja0=M0[k][i][4]+M0[k][i][6],jb0=M0[k][i][5]+M0[k][i][7];
     s+=(ja-H(R(ja0)))*sym(rad[0]);s+=(jb-H(R(jb0)))*sym(rad[1]);
    }else for(int j=0;j<Q;j++)s+=(box.r[k].M[i][j]-H(R(M0[k][i][j])))*sym(mesh.rstate);
    s+=(ci[i]-H(R(c0[k][i])))*sym(mesh.rtau);Dmag[Q*k+i]=mag(s);
   }
  }
  VH lastBox;for(int i=0;i<Q;i++)lastBox[i]=mesh.node[N][i]+sym(mesh.rstate);auto dbi=boundaryD(lastBox);
  for(int i=0;i<3;i++){H s(0);for(int j=0;j<Q;j++)s+=(dbi[i][j]-H(R(db0[i][j])))*sym(mesh.rstate);Dmag[bb+i]=mag(s);}
  std::vector<H>Ebox(D);for(int i=0;i<D;i++)Ebox[i]=sym(Dmag[i]);auto CE=denseApply(C,D,Ebox);
  R inflate=R("1.000000000001"),qR=R(qRld)*inflate,maxCorr=0,maxDeriv=0,maxK=0,minMargin=R("1e100"),qTot=0;int worst=-1;bool inclusion=true;
  for(int i=0;i<D;i++){
   R rr=R(Rradld[i])*inflate;H K=-CF[i]+sym(rr)-CE[i];R ratio=mag(K)/rad[i],margin=std::min(lo(K)+rad[i],rad[i]-hi(K));
   maxCorr=std::max(maxCorr,mag(CF[i])/rad[i]);maxDeriv=std::max(maxDeriv,mag(CE[i])/rad[i]);qTot=std::max(qTot,(rr+mag(CE[i]))/rad[i]);
   if(ratio>maxK){maxK=ratio;worst=i;}minMargin=std::min(minMargin,margin);if(!(lo(K)>-rad[i]&&hi(K)<rad[i]))inclusion=false;
  }
  bool contraction=qTot<1;
  bool nonsingular=qR<1,proof=mpfi&&structural&&domain&&(box.globalD2>0)&&nonsingular&&contraction&&inclusion;
  out<<"MAX_F_INTERVAL_MAG "<<maxF<<" MAX_F_INDEX "<<maxFIndex<<"\nQR_NORM_BOUND "<<qR<<"\nPRECONDITIONER_NONSINGULAR "<<nonsingular<<"\n";
  out<<"MAX_CORRECTION_RATIO "<<maxCorr<<"\nMAX_DERIVATIVE_RATIO "<<maxDeriv<<"\nTOTAL_CONTRACTION_BOUND "<<qTot<<" CONTRACTION_LT_ONE "<<contraction<<"\n";
  out<<"KRAWCZYK_MAX_RATIO "<<maxK<<" WORST_INDEX "<<worst<<"\nKRAWCZYK_MIN_MARGIN "<<minMargin<<"\nKRAWCZYK_STRICT_INCLUSION "<<inclusion<<"\nPROOF_OK "<<proof<<"\n";
  return proof?0:10;
 }catch(std::exception&e){std::cerr<<"ERROR "<<e.what()<<"\n";return 99;}
}