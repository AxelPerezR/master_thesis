#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>

void energy_resolution() {
    gStyle->SetOptFit(0);

    // Load yo file and tree
    TFile* file = TFile::Open("eicrecon_0.1GeV_NoCu.root");
    TTree* tree = (TTree*)file->Get("events");

    // Define alias
    tree->SetAlias("eta", 
        "0.5*log((sqrt(MCParticles.momentum.x[0]*MCParticles.momentum.x[0] + MCParticles.momentum.y[0]*MCParticles.momentum.y[0] + MCParticles.momentum.z[0]*MCParticles.momentum.z[0]) + MCParticles.momentum.z[0]) / (sqrt(MCParticles.momentum.x[0]*MCParticles.momentum.x[0] + MCParticles.momentum.y[0]*MCParticles.momentum.y[0] + MCParticles.momentum.z[0]*MCParticles.momentum.z[0]) - MCParticles.momentum.z[0]))");

    tree->SetAlias("ClusterEnergy", "Max$(EcalEndcapNClusters.energy)");    

    // Create histogram
    TH1D* h_ep = new TH1D("h_ep", ";E_{reconstructed}/p_{thrown};Normalized Number of Events", 105, 0., 1.05);
    // TH1D* h_ep = new TH1D("h_ep", "E/p;E/p;Events", 101, 0., 1.01);

    // Fill histogram with max cluster energy normalized by beam energy
    double beamEnergy = 0.1; // GeV
    tree->Draw(Form("ClusterEnergy/%f >> h_ep", beamEnergy), "eta > -3.5 && eta < -2 && ClusterEnergy > 0", "goff");

    // Normalize histogram
    h_ep->Scale(1.0 / h_ep->GetMaximum());

    h_ep->SetLineColor(kBlue);       // Change line color
    h_ep->SetLineWidth(2);           // Make line thicker

    // Fit the histogram with a Crystal Ball

    // The followup fit is inspired by Dimitry's work

    int startBin = 11;
    int endBin = h_ep->GetNbinsX();
    double normGuess = h_ep->Integral(startBin, endBin);
    std::cout << "Normalization guess: " << normGuess << std::endl;

    double mu_guess = h_ep->GetMean();
    double sigma_guess = h_ep->GetRMS();
    double norm_guess = h_ep->Integral(h_ep->FindBin(0.5), h_ep->FindBin(1.1));

    TF1* fitFunc = new TF1("fitFunc", "crystalball", 0.5, 1.1);
    //instead of mu guess, use 0.95
    // fitFunc->SetParameters(norm_guess, mu_guess, sigma_guess, 2.0, 3.0); // n (normalization factor), mean, sigma, alpha (tail threshold), n (tail exponent)
    fitFunc->SetParameters(1, mu_guess, 0.05, 2.0, 3.0); // n (normalization factor), mean, sigma, alpha (tail threshold), n (tail exponent)
    
    fitFunc->SetParLimits(2, 0.001, 1.0);  // Sigma must be > 0
    fitFunc->SetParLimits(3, 0.5, 5.0);    // Alpha
    fitFunc->SetParLimits(4, 1.01, 10.0);  // n > 1
    
    h_ep->Fit(fitFunc, "R");

    // Extract resolution
    double mu = fitFunc->GetParameter(1);      // Mean
    double mu_err = fitFunc->GetParError(1);
    double sigma = fitFunc->GetParameter(2);   // Std. deviation (from Crystal Ball core)
    double sigma_err = fitFunc->GetParError(2);

    double resolution = sigma / mu;
    std::cout << "Relative Energy Resolution (sigma/E): " << resolution*100 << " %" << std::endl;

    // Determine the resolution from FWHM
    // Find maximum value of the fit function (y_max)
    double y_max = fitFunc->GetMaximum(0.5, 1.1);

    // Define half-maximum value
    double half_max = y_max / 2.0;
    // double half_max = 0.5;
    // std::cout << "Half-maximum value: " << half_max << std::endl;

    // Find x values at half-maximum
    // double x_minus = fitFunc->GetX(half_max, mu - 2*sigma, mu);
    // double x_plus  = fitFunc->GetX(half_max, mu, mu + 2*sigma);
    double x_minus = fitFunc->GetX(half_max, mu - 3*sigma, mu);
    double x_plus  = fitFunc->GetX(half_max, mu, mu + 3*sigma);

    // Compute FWHM and relative resolution
    double fwhm = x_plus - x_minus;
    double sigma_rel_fwhm = (fwhm / 2.0) / (std::sqrt(2.0 * std::log(2.0)) * mu);

    // Approximate fwhm error using derivative w.r.t sigma (conservatively)
    // since x_minus/x_plus depend on sigma range; use sigma_err × 2 as proxy
    double fwhm_err = 2.0 * sigma_err;  // rough conservative estimate

    // Error propagation
    double denom = 2.0 * std::sqrt(2.0 * std::log(2.0)) * mu;
    double rel_err_fwhm = std::sqrt(
        std::pow(fwhm_err / denom, 2) +
        std::pow((fwhm * mu_err) / (denom * mu), 2)
    );

    // Output
    std::cout << "Relative resolution (sigma/E from FWHM): " 
            << sigma_rel_fwhm * 100 << " ± " << rel_err_fwhm * 100 << " %" << std::endl;

    //std::cout << "Relative resolution (sigma/E from FWHM): " << sigma_rel_fwhm * 100 << " %" << std::endl;

    // Plotting    
    TCanvas* c1 = new TCanvas("c1", "E/p Histogram", 600, 600);
    h_ep->Draw("E1 HIST");
    //fitFunc->Draw("same");
    //c1->SetGrid();

    // Draw FWHM lines
    TLine* line_low = new TLine(x_minus, 0, x_minus, y_max);
    TLine* line_high = new TLine(x_plus, 0, x_plus, y_max);
    line_low->SetLineStyle(2);   // dotted
    line_high->SetLineStyle(2);  // dotted
    line_low->SetLineWidth(3);
    line_high->SetLineWidth(3);
    line_low->SetLineColor(kRed);
    line_high->SetLineColor(kRed);
    line_low->Draw("same");
    line_high->Draw("same");

    //gPad->SetLogy();
    //h_ep->SetMinimum(1e-3);
    h_ep->SetStats(0);

    TPaveStats* stats = nullptr;
    stats = (TPaveStats*)h_ep->GetListOfFunctions()->FindObject("stats");
    if (!stats) stats = (TPaveStats*)gPad->GetPrimitive("stats");

    if (stats) {
        // Clear & repopulate existing box
        stats->Clear();
    } else {
        // Create a new stats box if none exists
        stats = new TPaveStats(0.65, 0.60, 0.90, 0.85, "NDC");
        stats->SetBorderSize(1);
        stats->SetFillColor(kWhite);
        stats->SetTextFont(42);
        stats->SetTextSize(0.04);
    }

    // Now add your custom lines 0.033
    stats->AddText("Parameters:"); 
    stats->AddText(Form(" #mu      = %.3f", mu));
    stats->AddText(Form(" #sigma   = %.3f",  sigma));
    stats->AddText(Form(" #sigma/E (%) = %.3f ", sigma_rel_fwhm*100));

    // If we just created it, draw it; otherwise Update() will suffice
    if (!h_ep->GetListOfFunctions()->FindObject(stats->GetName())) {
        stats->Draw();
    }

    // Redraw the pad
    c1->Modified();
    c1->Update();
}