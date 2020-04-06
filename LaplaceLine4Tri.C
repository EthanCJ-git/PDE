// Calculate the Voltage for a 2D system with wire at fixed voltage immediately above
// a 3-sided box at 0V

#include "TGraph2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TBox.h"
#include "TApplication.h"
#include "TImage.h"

#include <getopt.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using std::vector;
using std::cout;
using std::endl;

//static const double epsilon_0 = 625000.0/22468879468420441/M_PI;

// generic code to do one iteration of finite difference method
// Jacobi Method
double iterateJ(vector<vector<double>> &V, const double &rho, const double &delta, const int &plateTop, const int &plateBot, const int &plateLef, const int &plateRig){
  auto Vtmp = V;
  double dVmax=1e-50;
  int nx=V.size();
  int ny=V[0].size();
  double drho = 2*rho/(plateTop-plateBot);
  for (int i=1; i<nx-1; i++){
    for (int j=1; j<ny-1; j++){
      double Vnew = 0.25*(Vtmp[i+1][j]+Vtmp[i-1][j]+Vtmp[i][j+1]+Vtmp[i][j-1]);
      if(plateLef <= i && i <= plateRig && plateBot < j && j < plateTop)
	Vnew += M_PI*(rho-drho*(plateTop-j))*pow(delta,2);
      double dV=fabs(Vnew-V[i][j]);
      dVmax=std::max(dVmax,dV);    // keep track of max change in this sweep
      V[i][j] = Vnew;
    }
  }
  return dVmax;
}


// Gauss-Seidel Method
double iterateGS(vector<vector<double>> &V){
  double dVmax=1e-50;
  int nx=V.size();
  int ny=V[0].size();
  for (int i=1; i<nx-1; i++){
    for (int j=1; j<ny-1; j++){
      double Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
      double dV=fabs(Vnew-V[i][j]);
      dVmax=std::max(dVmax,dV);    // keep track of max change in this sweep
      V[i][j] = Vnew;
    }
  }
  return dVmax;
}

// fill a TGraph2D object from a vector of voltages
// delta: grid spacing
// the optional range parameter defines the subregion to plot
void fillGraph(TGraph2D* tg, const vector<vector<double>> &V, double delta, TBox *range=0){
  int nx=V.size();
  int ny=V[0].size();
  tg->Clear();                 // reset the graph
  for (int i=0; i<nx; i++){
    double x = i*delta;
    for (int j=1; j<ny; j++){
      double y = j*delta;
      if (range && range->IsInside(x,y))
	tg->SetPoint(tg->GetN(),x,y,V[i][j]);
    }
  }
}

double boundary(const double &x, const double &w) {
  if(x <= w/2)
    return 200.0*x/w;
  return 100.0*(1-1.0*x/w);
}

// Define box 0<x<L, 0<y<L
// potential on top edge at y=L
// eps: convergence criteria (max size of change at any grid point in an iteration)
// maxIter: max iterations in case of non-converence
// Npts : smoothness parameter, number of grid points in x,y
// pass a tcanvas for an animated solution, with specified max rate of frames/second
TGraph2D* LaplaceLine(int maxIter=100, double eps=0.001, int Npts=100, TCanvas *tc=0, int rate=10){
  double L=100;            // length of any side
  double Vtop=100;         // Voltage at top of box
  int maxgraphlines=200;   // max lines to draw in each direction
  
  vector<vector<double>> V(Npts, vector<double> (Npts, 0));  // create N x N vector, init to 0
  double delta = L/(Npts-1);                                 // grid spacing
  int plateTop = Npts/3*2;
  int plateBot = Npts/3;
  int plateLef = Npts/6;
  int plateRig = Npts/6*5;
  int width = Npts/12;
  double rho = 0.25;
  for (int i=plateLef; i<=plateRig; i++)
    for( int j=0; j<=width; j++) {
      V[i][plateTop+j] = Vtop;            // set voltage at wire
      V[i][plateBot-j] = -Vtop;
    }

  for(int j=0; j<Npts; j++) {
    V[0][j] = boundary(0, L);
    V[Npts-1][j] = boundary(L, L);
  }
  for(int i=1; i<Npts-1; i++) {
    double value = boundary(i*delta, L);
    V[i][0] = value;
    V[i][Npts-1] = value;
  }

  int msec = 1000/rate;                                      // milliseconds sleep between frames
  TBox *plotRange = new TBox(0,0,1.1*L,1.1*L);

  TGraph2D* tgV = new TGraph2D();                            // graph to store result
  if (Npts<50) tgV->SetLineWidth(3);                         
  tgV->SetLineColor(kRed);
  tgV->SetNpx(std::min(maxgraphlines,Npts));  tgV->SetNpy(std::min(maxgraphlines,Npts)); 
  tgV->SetTitle("Voltage;x;y;V");
  
  double dV;
  int niter=0;
  do{
    dV=iterateJ(V, rho, delta, plateTop, plateBot, plateLef, plateRig);   // iterate using Jacobi method
    //dV=iterateGS(V);   // iterate using Gauss-Seidel method
    for (int i=plateLef; i<=plateRig; i++)
      for (int j=0; j<=width; j++) {
	V[i][plateTop+j] = Vtop;
	V[i][plateBot-j] = -Vtop;
      }
    ++niter;
    if (tc) {
      tc->cd();
      fillGraph(tgV,V,delta,plotRange);
      tgV->Draw("surf");
      tc->Update();
      gSystem->Sleep(msec);
    }
  } while (dV>eps && niter<maxIter);
  
  cout << "Ended calculation with " << niter << " iterations, dVmax = " << dV << endl;
 
  fillGraph(tgV,V,delta,plotRange);
  return tgV;
}

void usage(char *prog){
  std::cerr << "Usage: " << prog << " <option(s)> SOURCES"
	    << "Options:\n"
	    << "\t-h\t\tShow this help message\n"
	    << "\t-a\t\tDisplay animation of solution"
    	    << "\t-I\t\t(max) Number of iterations [100]"
	    << "\t-e\t\tconvergence criteria [0.001]"
    	    << "\t-N\t\tNumber of points in x,y [100]"
	    << "\t-R\t\tmax frames/second with animation [10]"
	    << std::endl;
  exit(0);
}


int main(int argc, char *argv[]){
  TApplication theApp("App", &argc, argv, NULL, -1);  // -1 disables ROOT arg processing

  // defaults for LaplaceLine
  int maxIter=100;
  double eps=0.001;
  int Npts=100;
  TCanvas *tc=0;
  int rate=10;
  
  int opt;
  while ((opt = getopt(argc, argv, "haI:e:N:r:")) != -1) {
    switch (opt) {
    case 'h':
      usage(argv[0]);
      break;
    case 'a':
      tc=new TCanvas();
      break;
    case 'I':
      maxIter=atoi(optarg);
      break;
     case 'e':
      eps=atof(optarg);
      break; 
    case 'N':
      Npts=atoi(optarg);
      break;
    case 'r':
      rate=atoi(optarg);
      break;
    }
  }
 
  auto tg=LaplaceLine(maxIter,eps,Npts,tc,rate);

  // display final result
  if (!tc) tc=new TCanvas();
  tg->Draw("surf");              // explore other drawing options!
  
  TImage* img = TImage::Create();
  img->FromPad(tc);
  img->WriteImage("LaplaceLine4Tri.png");

  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}

