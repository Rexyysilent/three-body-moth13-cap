#pragma STDC FENV_ACCESS ON
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/transc.hpp>
#include <boost/numeric/interval/rounded_transc.hpp>
#include <array>
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace bn=boost::numeric;
namespace il=boost::numeric::interval_lib;
using RR=il::save_state<il::rounded_transc_std<long double>>;
using PP=il::policies<RR,il::checking_base<long double>>;
using I=bn::interval<long double,PP>;
static constexpr int Q=8;
using VI=std::array<I,Q>;
using VP=std::array<long double,Q>;
using MI=std::array<std::array<I,Q>,Q>;
using MP=std::array<std::array<long double,Q>,Q>;
using P3=std::array<long double,3>;
using I3=std::array<I,3>;
using DBI=std::array<std::array<I,Q>,3>;
using DBP=std::array<std::array<long double,Q>,3>;
using P83I=std::array<std::array<I,3>,Q>;

static long double lo(const I&x){return bn::lower(x);}
static long double hi(const I&x){return bn::upper(x);}
static long double mag(const I&x){return std::max(std::fabs(lo(x)),std::fabs(hi(x)));}
static long double mid(const I&x){return lo(x)+(hi(x)-lo(x))/2;}
static long double parse(const std::string&s){return std::strtold(s.c_str(),nullptr);}
static I sym(long double r){return I(-r,r);}
static bool zeroI(const I&x){return lo(x)==0 && hi(x)==0;}

struct Mesh{
  int n=0;
  long double tau=0,rtau=0,rstate=0;
  std::vector<long double> den;
  std::vector<VP> node;
};
static Mesh readMesh(const std::string&path){
  std::ifstream f(path); if(!f)throw std::runtime_error("cannot open mesh");
  Mesh m; std::string z;
  f>>m.n>>z; m.tau=parse(z);
  f>>z; m.rtau=parse(z); f>>z; m.rstate=parse(z);
  m.den.resize(m.n);
  for(auto&d:m.den){f>>z;d=parse(z);}
  m.node.resize(m.n+1);
  for(int k=0;k<=m.n;k++)for(int j=0;j<Q;j++){f>>z;m.node[k][j]=parse(z);}
  if(!f)throw std::runtime_error("truncated mesh");
  return m;
}
struct Rec{
  VI xc{},xb{}; MI M{};
  long double steps=0,minpow=0,mind2=0,margin=0,maxmw=0;
};
struct Log{
  std::vector<Rec> r;
  long double globalD2=0;
};
static I readInterval(std::istream&f){
  std::string a,b;f>>a>>b;return I(parse(a),parse(b));
}
static Log readLog(const std::string&path,int n){
  std::ifstream f(path);if(!f)throw std::runtime_error("cannot open log "+path);
  Log L;L.r.resize(n);std::string t,z;
  for(int k=0;k<n;k++){
    int idx,ok;
    f>>t;if(t!="REC")throw std::runtime_error("expected REC "+std::to_string(k)+" got "+t);
    f>>idx>>t>>ok;if(t!="OK"||idx!=k||ok!=1)throw std::runtime_error("bad REC");
    f>>t;if(t!="META")throw std::runtime_error("expected META");
    f>>z;L.r[k].steps=parse(z);f>>z;L.r[k].minpow=parse(z);
    f>>z;L.r[k].mind2=parse(z);f>>z;L.r[k].margin=parse(z);f>>z;L.r[k].maxmw=parse(z);
    f>>t;if(t!="XC")throw std::runtime_error("expected XC");
    for(int j=0;j<Q;j++)L.r[k].xc[j]=readInterval(f);
    f>>t;if(t!="XB")throw std::runtime_error("expected XB");
    for(int j=0;j<Q;j++)L.r[k].xb[j]=readInterval(f);
    f>>t;if(t!="MB")throw std::runtime_error("expected MB");
    for(int i=0;i<Q;i++)for(int j=0;j<Q;j++)L.r[k].M[i][j]=readInterval(f);
  }
  int ok,count;
  f>>t>>ok>>t>>count>>t>>z;
  if(t!="GLOBAL_D2"||ok!=1||count!=n)throw std::runtime_error("bad BATCH_OK footer");
  L.globalD2=parse(z);
  return L;
}
static VI field(const VI&x){
  I mu=I(1001)/I(1000);
  I q3x=-(x[0]+x[2])/mu,q3y=-(x[1]+x[3])/mu;
  auto gv=[](I rx,I ry){
    I r2=rx*rx+ry*ry;
    if(!(lo(r2)>0))throw std::runtime_error("field collision");
    I r3=r2*sqrt(r2);
    return std::array<I,2>{rx/r3,ry/r3};
  };
  auto g12=gv(x[2]-x[0],x[3]-x[1]);
  auto g13=gv(q3x-x[0],q3y-x[1]);
  auto g21=gv(x[0]-x[2],x[1]-x[3]);
  auto g23=gv(q3x-x[2],q3y-x[3]);
  VI f;
  f[0]=x[4];f[1]=x[5];f[2]=x[6];f[3]=x[7];
  f[4]=g12[0]+mu*g13[0];f[5]=g12[1]+mu*g13[1];
  f[6]=g21[0]+mu*g23[0];f[7]=g21[1]+mu*g23[1];
  return f;
}
static I3 boundary(const VI&x){
  return {x[0]+x[2],x[1]+x[3],
          x[0]*(x[4]-x[6])+x[1]*(x[5]-x[7])};
}
static DBI boundaryD(const VI&x){
  DBI d;for(auto&r:d)for(auto&v:r)v=I(0);
  d[0][0]=I(1);d[0][2]=I(1);
  d[1][1]=I(1);d[1][3]=I(1);
  d[2][0]=x[4]-x[6];d[2][1]=x[5]-x[7];
  d[2][4]=x[0];d[2][5]=x[1];d[2][6]=-x[0];d[2][7]=-x[1];
  return d;
}
static std::array<std::array<long double,3>,3>
inverse3(std::array<std::array<long double,3>,3>a){
  std::array<std::array<long double,6>,3> w{};
  for(int i=0;i<3;i++){for(int j=0;j<3;j++)w[i][j]=a[i][j];w[i][3+i]=1;}
  for(int k=0;k<3;k++){
    int p=k;for(int i=k+1;i<3;i++)if(std::fabs(w[i][k])>std::fabs(w[p][k]))p=i;
    if(w[p][k]==0)throw std::runtime_error("singular Schur midpoint");
    if(p!=k)std::swap(w[p],w[k]);
    long double v=w[k][k];for(int j=0;j<6;j++)w[k][j]/=v;
    for(int i=0;i<3;i++)if(i!=k){v=w[i][k];for(int j=0;j<6;j++)w[i][j]-=v*w[k][j];}
  }
  std::array<std::array<long double,3>,3> b{};
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)b[i][j]=w[i][3+j];
  return b;
}
static I det3(const std::array<std::array<I,3>,3>&a){
  return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
       - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
       + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
}
static VP pmatvec(const MP&A,const VP&x){
  VP y{};for(int i=0;i<Q;i++)for(int j=0;j<Q;j++)y[i]+=A[i][j]*x[j];return y;
}
static P3 dbvec(const DBP&A,const VP&x){
  P3 y{};for(int i=0;i<3;i++)for(int j=0;j<Q;j++)y[i]+=A[i][j]*x[j];return y;
}
struct Data{
  int n=0,D=0;
  std::vector<MP>M;
  std::vector<VP>c;
  DBP db{};
  std::array<std::array<long double,3>,3> sinv{};
};
static std::vector<long double>
pointSolve(const Data&d,const std::vector<VP>&y,const P3&yb){
  VP w{};
  for(int k=0;k<d.n;k++){
    VP z=pmatvec(d.M[k],w);
    for(int i=0;i<Q;i++)z[i]-=y[k][i];
    w=z;
  }
  P3 rhs=yb,dbw=dbvec(d.db,w);
  for(int i=0;i<3;i++)rhs[i]-=dbw[i];
  P3 p{};
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)p[i]+=d.sinv[i][j]*rhs[j];
  std::vector<long double> out(d.D,0);
  for(int i=0;i<3;i++)out[i]=p[i];
  VP cur{};cur[4]=p[0];cur[5]=p[1];cur[6]=p[0];cur[7]=p[1];
  for(int k=0;k<d.n;k++){
    VP z=pmatvec(d.M[k],cur);
    for(int i=0;i<Q;i++)z[i]+=d.c[k][i]*p[2]-y[k][i];
    for(int i=0;i<Q;i++)out[3+Q*k+i]=z[i];
    cur=z;
  }
  return out;
}
static std::vector<I>
denseApply(const std::vector<long double>&C,int D,const std::vector<I>&v){
  std::vector<I>o(D,I(0));
  for(int i=0;i<D;i++){
    I s(0);const long double*row=&C[(std::size_t)i*D];
    for(int j=0;j<D;j++)if(!zeroI(v[j]))s+=I(row[j])*v[j];
    o[i]=s;
  }
  return o;
}
int main(int argc,char**argv){
  try{
    if(argc!=5){std::cerr<<"usage: krawczyk mesh center.log box.log certificate.log\n";return 2;}
    std::ofstream cert(argv[4]);if(!cert)throw std::runtime_error("cannot open certificate output");
    cert<<std::setprecision(20)<<std::scientific;
    I tenth=I(1)/I(10),root2=sqrt(I(2));
    bool rounding=(lo(tenth)<hi(tenth)&&lo(root2)<hi(root2)&&lo(root2)*lo(root2)<=2&&hi(root2)*hi(root2)>=2);
    cert<<"ROUNDING_OK "<<rounding<<" TENTH "<<lo(tenth)<<" "<<hi(tenth)
        <<" SQRT2 "<<lo(root2)<<" "<<hi(root2)<<"\n";
    if(!rounding)return 3;
    Mesh mesh=readMesh(argv[1]);
    Log cen=readLog(argv[2],mesh.n),box=readLog(argv[3],mesh.n);
    if(!(box.globalD2>0))throw std::runtime_error("nonpositive collision bound");
    const int N=mesh.n,D=Q*N+3;
    Data d;d.n=N;d.D=D;d.M.resize(N);d.c.resize(N);
    std::vector<MI> Mi(N);std::vector<VI> Ci(N);
    for(int k=0;k<N;k++){
      Mi[k]=box.r[k].M;
      for(int i=0;i<Q;i++)for(int j=0;j<Q;j++)d.M[k][i][j]=mid(Mi[k][i][j]);
      VI fi=field(box.r[k].xb);
      for(int i=0;i<Q;i++){Ci[k][i]=fi[i]/I(mesh.den[k]);d.c[k][i]=mid(Ci[k][i]);}
    }
    VI lastBox;
    for(int i=0;i<Q;i++)lastBox[i]=I(mesh.node[N][i])+sym(mesh.rstate);
    DBI dbi=boundaryD(lastBox);
    for(int i=0;i<3;i++)for(int j=0;j<Q;j++)d.db[i][j]=mid(dbi[i][j]);

    std::array<std::array<long double,3>,Q> Pp{};
    Pp[4][0]=1;Pp[5][1]=1;Pp[6][0]=1;Pp[7][1]=1;
    for(int k=0;k<N;k++){
      std::array<std::array<long double,3>,Q> z{};
      for(int i=0;i<Q;i++)for(int j=0;j<3;j++){
        for(int l=0;l<Q;l++)z[i][j]+=d.M[k][i][l]*Pp[l][j];
        if(j==2)z[i][j]+=d.c[k][i];
      }
      Pp=z;
    }
    std::array<std::array<long double,3>,3> Sm{};
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int l=0;l<Q;l++)Sm[i][j]+=d.db[i][l]*Pp[l][j];
    d.sinv=inverse3(Sm);
    std::array<std::array<I,3>,3> SI{};
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)SI[i][j]=I(Sm[i][j]);
    I sdet=det3(SI);
    cert<<"DIM "<<D<<" SEGMENTS "<<N<<"\n";
    cert<<"COLLISION_D2_LOWER "<<box.globalD2<<"\n";
    cert<<"SCHUR_POINT";
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)cert<<" "<<Sm[i][j];
    cert<<"\nSCHUR_POINT_DET_ENC "<<lo(sdet)<<" "<<hi(sdet)<<"\n";
    cert<<"SCHUR_INV";
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)cert<<" "<<d.sinv[i][j];
    cert<<"\n";

    cert<<"BUILD_C_BEGIN\n";
    std::vector<long double>C((std::size_t)D*D);
    std::vector<VP> y(N);P3 yb{};
    for(int col=0;col<D;col++){
      for(auto&v:y)v.fill(0);yb.fill(0);
      if(col<Q*N)y[col/Q][col%Q]=1;else yb[col-Q*N]=1;
      auto z=pointSolve(d,y,yb);
      for(int row=0;row<D;row++)C[(std::size_t)row*D+col]=z[row];
    }
    cert<<"BUILD_C_END\n";
    {std::ofstream q(std::string(argv[4])+".C.bin",std::ios::binary);q.write((char*)&D,sizeof(D));q.write((char*)C.data(),sizeof(long double)*C.size());}
    {std::ofstream q(std::string(argv[4])+".point.bin",std::ios::binary);q.write((char*)&N,sizeof(N));for(int k=0;k<N;k++){q.write((char*)d.M[k].data(),sizeof(long double)*Q*Q);q.write((char*)d.c[k].data(),sizeof(long double)*Q);}q.write((char*)d.db.data(),sizeof(long double)*3*Q);}


    std::vector<std::vector<std::pair<int,long double>>> Jcol(D);
    for(int k=0;k<N;k++){
      int rb=Q*k;
      if(k==0){
        for(int i=0;i<Q;i++){
          long double va=d.M[k][i][4]+d.M[k][i][6];
          long double vb=d.M[k][i][5]+d.M[k][i][7];
          Jcol[0].push_back({rb+i,va});Jcol[1].push_back({rb+i,vb});
        }
      }else{
        int xb=3+Q*(k-1);
        for(int j=0;j<Q;j++)for(int i=0;i<Q;i++)Jcol[xb+j].push_back({rb+i,d.M[k][i][j]});
      }
      for(int i=0;i<Q;i++)Jcol[2].push_back({rb+i,d.c[k][i]});
      int nb=3+Q*k;
      for(int i=0;i<Q;i++)Jcol[nb+i].push_back({rb+i,-1});
    }
    int xbN=3+Q*(N-1),bb=Q*N;
    for(int j=0;j<Q;j++)for(int i=0;i<3;i++)if(d.db[i][j]!=0)Jcol[xbN+j].push_back({bb+i,d.db[i][j]});

    std::vector<I> F(D,I(0));
    long double maxF=0;
    for(int k=0;k<N;k++)for(int i=0;i<Q;i++){
      F[Q*k+i]=cen.r[k].xc[i]-I(mesh.node[k+1][i]);
      maxF=std::max(maxF,mag(F[Q*k+i]));
    }
    VI lastPoint;for(int i=0;i<Q;i++)lastPoint[i]=I(mesh.node[N][i]);
    I3 bf=boundary(lastPoint);
    for(int i=0;i<3;i++){F[bb+i]=bf[i];maxF=std::max(maxF,mag(bf[i]));}
    cert<<"MAX_F_INTERVAL_MAG "<<maxF<<"\n";
    auto CF=denseApply(C,D,F);

    std::vector<long double> rad(D,mesh.rstate);
    rad[2]=mesh.rtau;
    std::vector<long double>Dmag(D,0);
    for(int k=0;k<N;k++){
      int rb=Q*k;
      for(int i=0;i<Q;i++){
        I s(0);
        if(k==0){
          I ja=Mi[k][i][4]+Mi[k][i][6];
          I jb=Mi[k][i][5]+Mi[k][i][7];
          long double ja0=d.M[k][i][4]+d.M[k][i][6];
          long double jb0=d.M[k][i][5]+d.M[k][i][7];
          s+=(ja-I(ja0))*sym(rad[0]);
          s+=(jb-I(jb0))*sym(rad[1]);
        }else{
          for(int j=0;j<Q;j++)s+=(Mi[k][i][j]-I(d.M[k][i][j]))*sym(mesh.rstate);
        }
        s+=(Ci[k][i]-I(d.c[k][i]))*sym(mesh.rtau);
        Dmag[rb+i]=mag(s);
      }
    }
    for(int i=0;i<3;i++){
      I s(0);
      for(int j=0;j<Q;j++)s+=(dbi[i][j]-I(d.db[i][j]))*sym(mesh.rstate);
      Dmag[bb+i]=mag(s);
    }
    std::vector<I> Ebox(D,I(0));
    for(int j=0;j<D;j++)Ebox[j]=sym(Dmag[j]);
    auto CE=denseApply(C,D,Ebox);

    std::vector<long double> Rrad(D,0);
    long double qR=0,maxCorrRatio=0,maxERatio=0,maxKRatio=0,minMargin=std::numeric_limits<long double>::infinity();
    int worst=-1;
    std::vector<I>K(D,I(0));
    for(int i=0;i<D;i++){
      I rr(0);
      const long double*cr=&C[(std::size_t)i*D];
      for(int j=0;j<D;j++){
        I q(i==j?1:0);
        for(const auto&e:Jcol[j])q-=I(cr[e.first])*I(e.second);
        rr+=q*sym(rad[j]);
      }
      Rrad[i]=mag(rr);
      qR=std::max(qR,Rrad[i]/rad[i]);
      K[i]=-CF[i]+rr-CE[i];
      long double ratio=mag(K[i])/rad[i];
      long double margin=std::min(lo(K[i])+rad[i],rad[i]-hi(K[i]));
      maxCorrRatio=std::max(maxCorrRatio,mag(CF[i])/rad[i]);
      maxERatio=std::max(maxERatio,mag(CE[i])/rad[i]);
      if(ratio>maxKRatio){maxKRatio=ratio;worst=i;}
      minMargin=std::min(minMargin,margin);
    }
    {std::ofstream q(std::string(argv[4])+".R.bin",std::ios::binary);q.write((char*)&D,sizeof(D));q.write((char*)&qR,sizeof(qR));q.write((char*)Rrad.data(),sizeof(long double)*Rrad.size());}
    bool nonsingular=(qR<1);
    bool inclusion=true;
    for(int i=0;i<D;i++)if(!(lo(K[i])>-rad[i]&&hi(K[i])<rad[i]))inclusion=false;
    cert<<"QR_NORM_BOUND "<<qR<<"\n";
    cert<<"PRECONDITIONER_NONSINGULAR "<<nonsingular<<"\n";
    cert<<"MAX_CORRECTION_RATIO "<<maxCorrRatio<<"\n";
    cert<<"MAX_DERIVATIVE_RATIO "<<maxERatio<<"\n";
    cert<<"KRAWCZYK_MAX_RATIO "<<maxKRatio<<" WORST_INDEX "<<worst<<"\n";
    cert<<"KRAWCZYK_MIN_MARGIN "<<minMargin<<"\n";
    cert<<"KRAWCZYK_STRICT_INCLUSION "<<inclusion<<"\n";
    int shown=0;
    for(int i=0;i<D&&shown<20;i++)if(mag(K[i])/rad[i]>.5L){
      cert<<"WIDE "<<i<<" "<<lo(K[i])<<" "<<hi(K[i])<<" RAD "<<rad[i]<<"\n";shown++;
    }
    bool proof=rounding&&(box.globalD2>0)&&nonsingular&&inclusion;
    cert<<"PROOF_OK "<<proof<<"\n";
    return proof?0:10;
  }catch(const std::exception&e){
    std::cerr<<"ERROR "<<e.what()<<"\n";return 99;
  }
}
