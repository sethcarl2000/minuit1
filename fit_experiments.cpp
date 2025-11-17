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

double exp_and_hist(double x, double A, double sigma, double mean, double exp_amplitude, double exp_factor )
{
    double arg = (x - mean)/sigma; 
    return exp_amplitude * exp( -x * exp_factor )  + A * exp( -0.5 * arg * arg ); 
}

int main(int argc, char* argv[])
{      
    //open the file
    auto file = new TFile("experiments.root", "READ");

    auto h1 = (TH1F*)file->Get("hexp1"); 
    auto h2 = (TH1F*)file->Get("hexp2"); 


    //create the minuit object 
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
    
    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6); 
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(4);    

    auto fit_two_histograms_Chi2 = [&h1, &h2](const double *par)
    {
        //these are the histograms
        int ipar=0; 
        const double& sigma_sig = par[ipar++];
        const double& mean_sig  = par[ipar++];
        
        const double& exp_amp_1 = par[ipar++];
        const double& exp_len_1 = par[ipar++];
        const double& A_1       = par[ipar++];
        
        const double& exp_amp_2 = par[ipar++];
        const double& exp_len_2 = par[ipar++];
        const double& A_2       = par[ipar++];

        double chi2_1 = Chi2_from_hist(h1, [=](double x){ return exp_and_hist(x, A_1, sigma_sig, mean_sig, exp_amp_1, exp_len_1); });
        double chi2_2 = Chi2_from_hist(h2, [=](double x){ return exp_and_hist(x, A_2, sigma_sig, mean_sig, exp_amp_2, exp_len_2); });
        
        return chi2_1 + chi2_2;
    };

    auto f_minimizer_2g = ROOT::Math::Functor(fit_two_histograms_Chi2, 8);
    
    minimizer->SetFunction(f_minimizer_2g);

    //set the list of variables
    int i_var=0; 
    minimizer->SetVariable(i_var++, "signal_sigma", 10.,    1e-4);
    minimizer->SetVariable(i_var++, "signal_mean",  75.,    1e-4);
    minimizer->SetVariable(i_var++, "exp_amp_1",    1708.,  1e-4);
    minimizer->SetVariable(i_var++, "exp_len_1",    0.075,  1e-4);
    minimizer->SetVariable(i_var++, "A1",           15.,    1e-4);
    minimizer->SetVariable(i_var++, "exp_amp_2",    279.,   1e-4);
    minimizer->SetVariable(i_var++, "exp_len_2",    0.025,  1e-4);
    minimizer->SetVariable(i_var++, "A2",           12.5,   1e-4);

    bool fit_status = minimizer->Minimize();

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
}