#include <Math/Factory.h> 
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <cmath> 
#include <vector>
#include <TAxis.h>
#include <TH1F.h> 
#include <TFile.h> 
#include <functional> 
#include <TF1.h> 
#include <TCanvas.h> 
#include <TStyle.h> 
#include <iostream> 
#include <stdexcept> 
#include <TRandom3.h> 

using namespace std; 

inline bool is_NaN(const double& x) { return x != x; }

inline double fcn_gaus(double x, double A, double mean, double sigma)
{
    double arg = (x-mean)/sigma;
    return A * exp( -0.5 * arg * arg ); 
}

double fcn_two_gaus(double x, double A1, double x1, double s1, double A2, double x2, double s2)
{
    return fcn_gaus(x,A1,x1,s1) + fcn_gaus(x,A2,x2,s2);
}

double fcn_gumbel(double x, double A, double mu, double beta)
{
    if (is_NaN(A) || is_NaN(mu) || is_NaN(beta)) {
        Error("fcn_gumbel", "Nan encountered! params: %f, %f, %f", A, mu, beta); 
        throw invalid_argument("See err above from <fcn_gumbel>.");
        return -1.; 
    }
    //printf("params: %f, %f, %f", A, mu, beta); cout << endl; 

    x = (x - mu)/beta; 
    return A * exp( -exp(-x) - x );
}

//compute the NLL for a given histogram and function
double NLL_from_hist(const TH1F* hist, const std::function<double(double)>& fcn)
{
    cout << "Computing NLL..." << endl; 
    double LL=0.; 
    auto ax = hist->GetXaxis(); 
    for (int b=1; b<=ax->GetNbins(); b++) {
        double expect  = fcn(ax->GetBinCenter(b));
        double observe = hist->GetBinContent(b);
        if (expect != expect) Error("NLL_from_hist", "NAN encountered! x = %f, bin=%i", ax->GetBinCenter(b), b);
        LL += observe * log(expect) - expect - (observe > 0. ? observe * log(observe) - observe : 0.); 
    }
    return -LL; 
}

double Chi2_from_hist(const TH1F* hist, const std::function<double(double)>& fcn)
{
    double chi2=0.;
    auto ax = hist->GetXaxis(); 
    for (int b=1; b<=ax->GetNbins(); b++) {
        double expect  = fcn(ax->GetBinCenter(b));
        double observe = hist->GetBinContent(b);

        if (observe < 1.) continue; 
        double arg = (expect - observe)/observe;
        chi2 += arg * arg; 
    }
    return chi2; 
}

int main(int argc, char* argv[])
{      
    //open the file
    auto file = new TFile("distros.root", "READ");

    auto dist1 = (TH1F*)file->Get("dist1"); 

    auto NLL_dist1_twogaus = [&dist1](const double *par)
    {
        //there should be 6 parameters
        auto test_fcn = [&par](double x){ return fcn_two_gaus(x, par[0],par[1],par[2],par[3],par[4],par[5]); };
        return NLL_from_hist(dist1, test_fcn);
    };

    //create the minuit object 
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
    
    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6); 
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(2);    

    auto f_minimizer_2g = ROOT::Math::Functor(NLL_dist1_twogaus, 6);
    
    minimizer->SetFunction(f_minimizer_2g);

    //set the list of variables
    int i_var=0; 
    minimizer->SetVariable(i_var++, "A1", 800., 1e-4);
    minimizer->SetVariable(i_var++, "x1", 80.,  1e-4);
    minimizer->SetVariable(i_var++, "s1", 7.,   1e-4);
    minimizer->SetVariable(i_var++, "A2", 100., 1e-4);
    minimizer->SetVariable(i_var++, "x2", 82., 1e-4);
    minimizer->SetVariable(i_var++, "s2", 12.5, 1e-4);

    bool fit_status = minimizer->Minimize();

    const double* par = minimizer->X(); 

    auto c = new TCanvas("c", "test", 1400, 800); 

    gStyle->SetOptStat(0); 

    c->Divide(2,1); 
    c->cd(1);
    dist1->SetTitle("Dist1 (twin-gaus fit)"); 
    dist1->DrawCopy("E"); 

    auto tf1_2g = new TF1("fit_result", [](double *x, double* p){
        return fcn_two_gaus(x[0], p[0],p[1],p[2],p[3],p[4],p[5]);
    }, dist1->GetXaxis()->GetXmin(), dist1->GetXaxis()->GetXmax(), 6); 

    tf1_2g->SetParameters(par); 
    tf1_2g->DrawCopy("SAME");

    //auto NLL_dist1_gumbel = 

    delete minimizer; 
    minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
    
    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6);
    minimizer->SetTolerance(1e-2);
    minimizer->SetPrintLevel(5); 
    
    auto f_minimizer_gumbel = ROOT::Math::Functor([&dist1](const double *par)
    {
        //there should be 6 parameters
        auto test_fcn = [&par](double x){ return fcn_gumbel(x, par[0],par[1],par[2]); };
        return Chi2_from_hist(dist1, test_fcn);
    }, 3);
    
    minimizer->SetFunction(f_minimizer_gumbel);

    //set the list of variables
    i_var=0; 
    minimizer->SetVariable(i_var++, "Amplitude", 3300., 1e-3);
    minimizer->SetVariable(i_var++, "mu",        80.,   1e-2);
    minimizer->SetVariable(i_var++, "beta",      7.5,   1e-4);


    fit_status = minimizer->Minimize();

    const double* par_gumbel = minimizer->X(); 
    
    auto tf1 = new TF1("fit_result_gumbel", [](double *x, double* p){
        return fcn_gumbel(x[0], p[0],p[1],p[2]);
    }, dist1->GetXaxis()->GetXmin(), dist1->GetXaxis()->GetXmax(), 3); 
    tf1->SetParameters(par_gumbel);

    c->cd(2); 
    dist1->SetTitle("Dist1 (gumbel fit)");
    dist1->DrawCopy("E");
    tf1->DrawCopy("SAME"); 

    c->SaveAs("ex1.png");

    return 0; 
}