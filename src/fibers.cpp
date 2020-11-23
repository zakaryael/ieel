//
//  fibers.h
//  i-eel
//
//  Created by Jeremie Bec on 23/11/2020.
//

#include "channelflow/fibers.h"

Fiber::Fiber() {
#ifdef HAVE_MPI
    cfmpi_ = &CfMPI::getInstance();
#endif
}

Fiber::Fiber(int Ns, Real L, Real mu, Real beta, Real gamma, CfMPI* cfmpi) {
#ifdef HAVE_MPI
    if (cfmpi == nullptr)
        cfmpi = &CfMPI::getInstance();
#endif
    assert(Ns >= 0);
    Ns_ = Ns;
    assert(L>0);
    L_ = L;
    ds_ = L_/(Real)Ns_;
    assert(mu>0);
    mu_ = mu;
    assert(beta>0);
    beta_ = beta;
    assert(gamma>0);
    gamma_ = gamma;
    Fiber::alloc();
}

void Fiber::read(int N, Real L, Real mu, Real beta, Real gamma, const std::string& filebase, int ifib, FlowField& u, FlowField& du, CfMPI* cfmpi) {
#ifdef HAVE_MPI
    if (cfmpi == nullptr)
        cfmpi = &CfMPI::getInstance();
#endif
    // Random initial orientation
    Ns_ = N;
    assert(L>0);
    L_ = L;
    ds_ = L_/(Real)Ns_;
    assert(mu>0);
    mu_ = mu;
    assert(beta>0);
    beta_ = beta;
    assert(gamma>0);
    gamma_ = gamma;
    
    Fiber::alloc();
    
    char cname[512];
    strcpy(cname, filebase.c_str());
    FILE *fin = fopen(cname,"r");
    if(fseek(fin,(Ns_+1)*9*ifib*sizeof(Real),SEEK_SET)>0) {
        cerr<<"Error in reading "<<filebase<<endl;
        exit(-1);
    }
    for(int is=0; is<=Ns_; is++) {
        Real xtmp[3];
        if(fread(&xtmp, sizeof(Real), 3, fin)==0) {
            cerr<<"Error in reading "<<filebase<<endl;
            exit(-1);
        }
        X_(is,0) = xtmp[0];
        X_(is,1) = xtmp[1];
        X_(is,2) = xtmp[2];
    }
    fclose(fin);
    
    Xold_ = X_;
    diff();
    interp_U(u,du);
    calc_tension();
    D1T_ = D1_*T_;
    for(int dim=0; dim<3; dim++)
        Gold_.col(dim) = Uf_.col(dim) + (2.0/mu_)*(D1T_%D1X_.col(dim)) + (1.0/mu_)*(T_%D2X_.col(dim));
}

void Fiber::init(int N, Real L, Real mu, Real beta, Real gamma, FlowField& u, FlowField& du, std::default_random_engine& rng, CfMPI* cfmpi) {
#ifdef HAVE_MPI
    if (cfmpi == nullptr)
        cfmpi = &CfMPI::getInstance();
#endif
    // Random initial orientation
    Ns_ = N;
    assert(L>0);
    L_ = L;
    ds_ = L_/(Real)Ns_;
    assert(mu>0);
    mu_ = mu;
    assert(beta>0);
    beta_ = beta;
    assert(gamma>0);
    gamma_ = gamma;
    
    Fiber::alloc();
    
    MyCol p(3);
    std::uniform_real_distribution<Real> unif(-1,1);
    Real rr = 2;
    while(rr>1) {
        p(0) = unif(rng);
        p(1) = 0; //unif(rng);
        p(2) = unif(rng);
        rr = p(0)*p(0)+p(1)*p(1)+p(2)*p(2);
    }
    rr = sqrt(rr);
    p(0) /= rr; p(1) /= rr; p(2) /= rr;
    
    std::uniform_real_distribution<Real> unifX(0,u.Lx());
    std::uniform_real_distribution<Real> unifY(max(u.a(),u.a()-L_*p(0)),min(u.b(),u.b()-L_*p(0)));
    std::uniform_real_distribution<Real> unifZ(0,u.Lz());
    X_(0,0) = unifX(rng);
    X_(0,1) = unifY(rng);
    X_(0,2) = unifZ(rng);
    Xold_(0,0) = X_(0,0); Xold_(0,1) = X_(0,1); Xold_(0,2) = X_(0,2);
    for(int is=1; is<=Ns_; is++) {
        X_.row(is) = X_.row(is-1)+ds_*p.t();
        Xold_.row(is) = X_.row(is);
        if((X_(is,1)<u.a())||(X_(is,1)>u.b())) {
            cerr<<"Error in initializing the fiber!"<<endl;
            cerr<<"is = "<<is<<", Y(is) = "<<X_(is,1)<<endl;
            exit(-1);
        }
    }
    diff();
    interp_U(u,du);
    calc_tension();
    D1T_ = D1_*T_;
    for(int dim=0; dim<3; dim++)
        Gold_.col(dim) = Uf_.col(dim) + (2.0/mu_)*(D1T_%D1X_.col(dim)) + (1.0/mu_)*(T_%D2X_.col(dim));
}

void Fiber::alloc() {
    X_.set_size(Ns_+1,3);
    Xold_.set_size(Ns_+1,3);
    Gold_.set_size(Ns_+1,3);
    D1X_.set_size(Ns_+1,3);
    D2X_.set_size(Ns_+1,3);
    D3X_.set_size(Ns_+1,3);
    D4X_.set_size(Ns_+1,3);
    
    T_.set_size(Ns_+1);
    D1T_.set_size(Ns_+1);
    NormXi_.set_size(Ns_+1);
    
    Uf_.set_size(Ns_+1,3);
    D1Uf_.set_size(Ns_+1,3);
    
    D1_.set_size(Ns_+1,Ns_+1);
    D1_.diag(1)  += 1;
    D1_.diag(-1) -= 1;
    D1_(0,0) = -3; D1_(0,1) = 4; D1_(0,2) = -1;
    D1_(Ns_,Ns_) = 3; D1_(Ns_,Ns_-1) = -4; D1_(Ns_,Ns_-2) = 1;
    D1_ /= 2.0*ds_;
    
    D2_.set_size(Ns_+1,Ns_+1);
    D2_.diag(1)  += 1;
    D2_.diag()   -= 2;
    D2_.diag(-1) += 1;
    D2_(0,0)=2; D2_(0,1)=-5; D2_(0,2)=4; D2_(0,3)=-1;
    D2_(Ns_,Ns_)=2; D2_(Ns_,Ns_-1)=-5; D2_(Ns_,Ns_-2)=4; D2_(Ns_,Ns_-3)=-1;
    D2_ /= ds_*ds_;
    
    LapDirich_.set_size(Ns_-1,Ns_-1);
    LapDirich_.diag() -= 2;
    LapDirich_.diag(1) += 1;
    LapDirich_.diag(-1) += 1;
    LapDirich_ /= ds_*ds_;
    
    D3_.set_size(Ns_+1,Ns_+1);
    D3_.diag(2)  += 1;
    D3_.diag(1)  -= 2;
    D3_.diag(-1) += 2;
    D3_.diag(-2) -= 1;
    D3_(0,0)=-5; D3_(0,1)=18; D3_(0,2)=-24; D3_(0,3)=14; D3_(0,4)=-3;
    D3_(1,0)=-3; D3_(1,1)=10; D3_(1,2)=-12; D3_(1,3)=6; D3_(1,4)=-1;
    D3_(Ns_,Ns_)=5; D3_(Ns_,Ns_-1)=-18; D3_(Ns_,Ns_-2)=24; D3_(Ns_,Ns_-3)=-14; D3_(Ns_,Ns_-4)=3;
    D3_(Ns_-1,Ns_)=3; D3_(Ns_-1,Ns_-1)=-10; D3_(Ns_-1,Ns_-2)=12; D3_(Ns_-1,Ns_-3)=-6; D3_(Ns_-1,Ns_-4)=1;
    D3_ /= 2.0*ds_*ds_*ds_;
    
    D4_.set_size(Ns_+1,Ns_+1);
    D4_.diag(2)  += 1;
    D4_.diag(1)  -= 4;
    D4_.diag()   += 6;
    D4_.diag(-1) -= 4;
    D4_.diag(-2) += 1;
    D4_(0,0)=3; D4_(0,1)=-14; D4_(0,2)=26; D4_(0,3)=-24; D4_(0,4)=11; D4_(0,5)=-2;
    D4_(1,0)=2; D4_(1,1)=-9; D4_(1,2)=16; D4_(1,3)=-14; D4_(1,4)=6; D4_(1,5)=-1;
    D4_(Ns_,Ns_)=3; D4_(Ns_,Ns_-1)=-14; D4_(Ns_,Ns_-2)=26; D4_(Ns_,Ns_-3)=-24; D4_(Ns_,Ns_-4)=11; D4_(Ns_,Ns_-5)=-2;
    D4_(Ns_-1,Ns_)=2; D4_(Ns_-1,Ns_-1)=-9; D4_(Ns_-1,Ns_-2)=16; D4_(Ns_-1,Ns_-3)=-14; D4_(Ns_-1,Ns_-4)=6; D4_(Ns_-1,Ns_-5)=-1;
    D4_ /= ds_*ds_*ds_*ds_;
    
    Op4_.set_size(Ns_-3,Ns_-3);
    Op4_.diag(2)  += 1;
    Op4_.diag(1)  -= 4;
    Op4_.diag()   += 6;
    Op4_.diag(-1) -= 4;
    Op4_.diag(-2) += 1;
    Op4_(0,0) = 2.0/11.0; Op4_(0,1) = -4.0/11.0; Op4_(0,2) = 2.0/11.0;
    Op4_(1,0) = -16.0/11.0; Op4_(1,1) = 43.0/11.0; Op4_(1,2) = -38.0/11.0;
    Op4_(Ns_-5,Ns_-4) = -16.0/11.0; Op4_(Ns_-5,Ns_-5) = 43.0/11.0; Op4_(Ns_-5,Ns_-6) = -38.0/11.0;
    Op4_(Ns_-4,Ns_-4) = 2.0/11.0; Op4_(Ns_-4,Ns_-5) = -4.0/11.0; Op4_(Ns_-4,Ns_-6) = 2.0/11.0;
    Op4_ /= ds_*ds_*ds_*ds_;
    
}

void Fiber::save(const std::string& filebase) const {
    char cname[512];
    strcpy(cname, filebase.c_str());
    FILE *fout = fopen(cname,"a");
    
    for(int is=0; is<=Ns_; is++) {
        fwrite(&X_(is,0), sizeof(Real), 1, fout);
        fwrite(&X_(is,1), sizeof(Real), 1, fout);
        fwrite(&X_(is,2), sizeof(Real), 1, fout);
    }
    
    MyMat V = (X_-Xold_)/dt_;
    for(int is=0; is<=Ns_; is++) {
        fwrite(&V(is,0), sizeof(Real), 1, fout);
        fwrite(&V(is,1), sizeof(Real), 1, fout);
        fwrite(&V(is,2), sizeof(Real), 1, fout);
    }
    
    for(int is=0; is<=Ns_; is++) {
        fwrite(&Uf_(is,0), sizeof(Real), 1, fout);
        fwrite(&Uf_(is,1), sizeof(Real), 1, fout);
        fwrite(&Uf_(is,2), sizeof(Real), 1, fout);
    }
    
    fclose(fout);
}

void Fiber::diff() {
    D1X_    = D1_*X_;
    D2X_    = D2_*X_;
    D3X_    = D3_*X_;
    D4X_    = D4_*X_;
    NormXi_ = sum(D1X_%D1X_,1);
}

void Fiber::calc_tension() {
    MySpMat Lap = 2.0*LapDirich_;
    Lap.diag() -= sum(D2X_.rows(1,Ns_-1)%D2X_.rows(1,Ns_-1),1);
    
    MyCol A = -mu_*sum(D1X_.rows(1,Ns_-1)%D1Uf_.rows(1,Ns_-1),1);
    A += 7.0*sum(D2X_.rows(1,Ns_-1)%D4X_.rows(1,Ns_-1),1)
    + 6.0*sum(D3X_.rows(1,Ns_-1)%D3X_.rows(1,Ns_-1),1);
    A += (mu_*beta_)*(1-NormXi_.rows(1,Ns_-1));
    
    MyCol  Ttmp;
    if(spsolve(Ttmp,Lap,A,"superlu")==false || isnan(Ttmp(0))) {
        cerr<<"No solution for the tension"<<endl;
        X_.save("X.dat",raw_ascii);
        exit(-1);
    }
    
    T_(0) = 0;
    T_.rows(1,Ns_-1) = Ttmp;
    T_(Ns_) = 0;
}

void Fiber::evol(Real dt, FlowField& u, FlowField& du) {
    
    diff();
    dt_ = dt;
    interp_U(u,du);
    
    calc_tension();
    D1T_ = D1_*T_;
    
    MyMat G(Ns_+1,3);
    for(int dim=0; dim<3; dim++)
        G.col(dim) = Uf_.col(dim) + (2.0/mu_)*(D1T_%D1X_.col(dim)) + (1.0/mu_)*(T_%D2X_.col(dim));
    
    MyMat RHS = (4.0/3.0)*X_.rows(2,Ns_-2)-(1.0/3.0)*Xold_.rows(2,Ns_-2)+(2*dt_/3.0)*(2*G.rows(2,Ns_-2)-Gold_.rows(2,Ns_-2));
    MyMat XXi = D1_*(2.0*X_-Xold_);
    XXi = XXi.rows(2,Ns_-2);
    
    Xold_ = X_;
    Gold_ = G;
    
    int nn = Ns_-3;
    MySpMat Axx(nn,nn), Ayy(nn,nn), Azz(nn,nn), Axy(nn,nn), Axz(nn,nn), Ayz(nn,nn);
    Axx.diag() = XXi.col(0) % XXi.col(0) + 1;
    Ayy.diag() = XXi.col(1) % XXi.col(1) + 1;
    Azz.diag() = XXi.col(2) % XXi.col(2) + 1;
    Axy.diag() = XXi.col(0) % XXi.col(1);
    Axz.diag() = XXi.col(0) % XXi.col(2);
    Ayz.diag() = XXi.col(1) % XXi.col(2);
    Axy = Axy*Op4_;
    Axz = Axz*Op4_;
    Ayz = Ayz*Op4_;
    
    MySpMat MM = speye<MySpMat>(3*nn,3*nn) + (2.0*dt_/(3.0*mu_))*join_cols(join_rows(Axx*Op4_,join_rows(Axy,Axz)), join_cols(join_rows(Axy,join_rows(Ayy*Op4_,Ayz)), join_rows(Axz,join_rows(Ayz,Azz*Op4_))));
    
    MyMat Xnew(3*nn,1);
    if(spsolve(Xnew,MM,vectorise(RHS),"superlu")==false || isnan(Xnew(0))) {
        cerr<<"No solution for elasticity"<<endl;
        X_.save("X.dat",raw_ascii);
        exit(-1);
    }
    
    Xnew = reshape(Xnew,nn,3);
    
    X_.row(0) = (48.0/11.0)*Xnew.row(0)-(52.0/11.0)*Xnew.row(1)+(15.0/11.0)*Xnew.row(2);
    X_.row(1) = (28.0/11.0)*Xnew.row(0)-(23.0/11.0)*Xnew.row(1)+(6.0/11.0)*Xnew.row(2);
    X_.rows(2,Ns_-2) = Xnew;
    X_.row(Ns_-1) = (28.0/11.0)*Xnew.row(Ns_-4)-(23.0/11.0)*Xnew.row(Ns_-5)+(6.0/11.0)*Xnew.row(Ns_-6);
    X_.row(Ns_)   = (48.0/11.0)*Xnew.row(Ns_-4)-(52.0/11.0)*Xnew.row(Ns_-5)+(15.0/11.0)*Xnew.row(Ns_-6);
    
}

void Fiber::interp_U(FlowField& u, FlowField& du) {
    
    const Vector YY = u.ygridpts();
    const int Nx = u.Nx();
    const int Ny = u.Ny();
    const int Nz = u.Nz();
    const Real Lx = u.Lx();
    const Real Lz = u.Lz();
    const Real a = u.a();
    const Real b = u.b();
    
    const Real c = 0.5 * (b + a);
    const Real r = 0.5 * (b - a);
    const Real piN = pi / (Real)(Ny - 1);
    const Real dx = Lx/(Real)Nx;
    const Real dz = Lz/(Real)Nz;
    
    const int nd = u.Nd();
    if(nd != 3)
        cerr << "Dimension Mismatch" << endl;
    
    Real dy, x, y, z, xmod, zmod;
    
    int i0, i1, j0, j1, k0, k1;
    
    for(int is=0; is<=Ns_; is++) {
        xmod = wrap(X_(is,0),Lx);
        i0 = (int)floor(Nx*xmod/Lx);
        i1 = (i0==Nx-1)?0:(i0+1);
        if(i0<0 || i1>= Nx) {
            cerr << endl << "Error in X position:" << endl;
            cerr << "Point " << is << " located at x = " << X_(is,0) << endl;
            exit(-1);
        }
        
        j1 = (int)floor(acos((X_(is,1)-c)/r)/piN);
        if(j1==Ny-1) j1--;
        j0 = j1+1;
        
        zmod = wrap(X_(is,2),Lz);
        k0 = (int)floor(Nz*zmod/Lz);
        k1 = (k0==Nz-1)?0:(k0+1);
        if(k0<0 || k1>= Nz) {
            cerr << endl << "Error in Z position" << endl;
            cerr << "Point " << is << " located at z = " << X_(is,2) << endl;
            exit(-1);
        }
        
        dy = YY[j1]-YY[j0];
        x = (xmod - i0*dx)/dx;
        y = (X_(is,1) - YY[j0])/dy;
        z = (zmod - k0*dz)/dz;
        
        // Trilinear interpolation
        if(j0<0 || j1>= Ny) {
            Uf_(is,0) = 0;
            Uf_(is,2) = 0;
            if(X_(is,1)<a) {
                Uf_(is,1) = -gamma_*(X_(is,1)-a);
                D1Uf_.row(is) = -gamma_*D1X_.row(is);
            }
            else if(X_(is,1)>b) {
                Uf_(is,1) = -gamma_*(X_(is,1)-b);
                D1Uf_.row(is) = -gamma_*D1X_.row(is);
            }
            else {
                cerr << endl << "Error in Y position" << endl;
                cerr << "Point " << is << " located at y = " << X_(is,1) << "has wrong behavior" << endl;
            }
        }
        else {
            for(int dim=0; dim<3; dim++) {
                Uf_(is,dim) = ((u(i0,j0,k0,dim)*(1-x)+u(i1,j0,k0,dim)*x)*(1-y)+(u(i0,j1,k0,dim)*(1-x)+u(i1,j1,k0,dim)*x)*y)*(1-z)
                    + ((u(i0,j0,k1,dim)*(1-x)+u(i1,j0,k1,dim)*x)*(1-y)+(u(i1,j1,k0,dim)*(1-x)+u(i1,j1,k0,dim)*x)*y)*z;
                D1Uf_(is,dim) = 0;
                for(int dd=0; dd<3; dd++)
                    D1Uf_(is,dim) += D1X_(is,dd)*( ((du(i0,j0,k0,i3j(dim,dd))*(1-x) +du(i1,j0,k0,i3j(dim,dd))*x)*(1-y) +(du(i0,j1,k0,i3j(dim,dd))*(1-x) +du(i1,j1,k0,i3j(dim,dd))*x)*y)*(1-z) + ((du(i0,j0,k1,i3j(dim,dd))*(1-x) +du(i1,j0,k1,i3j(dim,dd))*x)*(1-y) +(du(i1,j1,k0,i3j(dim,dd))*(1-x) +u(i1,j1,k0,i3j(dim,dd))*x)*y)*z );
            }
        }
    }
}
