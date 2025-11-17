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
#include <TH2D.h> 

using namespace std; 

inline bool is_NaN(const double& x) { return x != x; }


double Chi2_from_hist(const TH2D* hist, const std::function<double(double,double)>& fcn)
{
    double chi2=0.;
    auto xax = hist->GetXaxis();
    auto yax = hist->GetYaxis();  

    for (int bx=1; bx<=xax->GetNbins(); bx++) {
        for (int by=1; by<=yax->GetNbins(); by++) {
            double expect  = fcn(xax->GetBinCenter(bx), yax->GetBinCenter(by));
            double observe = hist->GetBinContent(bx,by);

            if (observe < 1.) continue; 
            double arg = (expect - observe)/observe;
            chi2 += arg * arg;
        } 
    }
    return chi2; 
}


double gaus_2d(double x, double y, double A, double x0, double y0, double sig_x, double sig_y)
{
    x = (x - x0)/sig_x;
    y = (y - y0)/sig_y; 
    
    return A * exp( -0.5 * ( x*x + y*y ) ); 
}

int main(int argc, char* argv[])
{      
    //open the file
    auto file = new TFile("fitInputs.root", "READ");

    auto h_data = (TH2D*)file->Get("hdata"); 
    auto h_bg = (TH2D*)file->Get("hbkg");

    //create the minuit object 
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
    
    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6); 
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(4);    

    auto f_minimizer_2g = ROOT::Math::Functor([&h_data, &h_bg](const double *par)
    {
        //these are the histograms
        int ipar=0; 
        const double& A_sig     = par[ipar++];
        const double& A_bg      = par[ipar++];
        
        const double& mean_x    = par[ipar++];
        const double& mean_y    = par[ipar++];
        const double& sigma_x   = par[ipar++];
        const double& sigma_y   = par[ipar++];

        auto fcn_fit = [&](double x, double y)
        {
            int bx = h_bg->GetXaxis()->FindBin(x); 
            int by = h_bg->GetYaxis()->FindBin(y); 

            double val = A_bg * h_bg->GetBinContent(bx, by); 

            return val + A_sig * gaus_2d(x,y, A_sig, mean_x,mean_y, sigma_x,sigma_y);
        };

        return Chi2_from_hist(h_data, fcn_fit);
    }, 6);
    
    minimizer->SetFunction(f_minimizer_2g);

    //set the list of variables
    int i_var=0; 
    minimizer->SetVariable(i_var++, "A_signal",     50.,    1e-4);
    minimizer->SetVariable(i_var++, "A_background", 10.,    1e-4);
    minimizer->SetVariable(i_var++, "mean_x",       1.2,    1e-4);
    minimizer->SetVariable(i_var++, "mean_y",       3.5,    1e-4);
    minimizer->SetVariable(i_var++, "sigma_x",      1.5,    1e-4);
    minimizer->SetVariable(i_var++, "sigma_y",      1.,     1e-4);

    bool fit_status = minimizer->Minimize();
    return 0; 
#if 0 

    const double* par = minimizer->X(); 

    auto c = new TCanvas("c", "test", 1400, 800); 
    c->Divide(2,1); 

    auto tf1_1 = new TF1("tf1_1", [par](double *x, double *p)
        {
            return exp_and_hist(x[0], p[0], p[1], p[2], p[3], p[4]); 
        }, 
        h1->GetXaxis()->GetXmin(), 
        h1->GetXaxis()->GetXmax(), 
        5
    ); 
    double params_1[] = {par[4], par[0], par[1], par[2], par[3]}; 
    tf1_1->SetParameters(params_1);
    
    c->cd(1);
    h1->Draw("E");
    tf1_1->DrawCopy("SAME");

    auto tf1_2 = new TF1("tf1_2", [par](double *x, double *p)
        {
            return exp_and_hist(x[0], p[0], p[1], p[2], p[3], p[4]); 
        }, 
        h2->GetXaxis()->GetXmin(), 
        h2->GetXaxis()->GetXmax(), 
        5
    ); 
    double params_2[] = {par[7], par[0], par[1], par[5], par[6]}; 
    tf1_2->SetParameters(params_2);

    c->cd(2);
    h2->GetYaxis()->SetRangeUser(0., h2->GetMaximum()*1.1); 
    h2->Draw("E");
    tf1_2->Draw("SAME");

    c->SaveAs("ex2.png"); 

    return 0; 
#endif 
}