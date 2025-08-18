void plot_acceptanceGap_etaFromEvent() {
    TFile *file = TFile::Open("eicrecon_175cm_500MeV_90to178deg.root");
    TTree *events = (TTree*)file->Get("events");

    events->SetAlias("eta", 
        "0.5*log((sqrt(MCParticles.momentum.x[0]*MCParticles.momentum.x[0] + MCParticles.momentum.y[0]*MCParticles.momentum.y[0] + MCParticles.momentum.z[0]*MCParticles.momentum.z[0]) + MCParticles.momentum.z[0]) / (sqrt(MCParticles.momentum.x[0]*MCParticles.momentum.x[0] + MCParticles.momentum.y[0]*MCParticles.momentum.y[0] + MCParticles.momentum.z[0]*MCParticles.momentum.z[0]) - MCParticles.momentum.z[0]))");

    events->SetAlias("Barrel", "Sum$(EcalBarrelScFiRecHits.energy)");
    events->SetAlias("EEEMCal", "Sum$(EcalEndcapNRecHits.energy)");

    // Create a TH2F histogram
    const int nEtaBins = 50;
    TH2F *hist2D = new TH2F("hist2D", 
        "Energy reconstruction of 500MeV e- (EEEMCal at z=175 cm);#eta_{thrown}; #Sigma E_{dep}/p_{thrown}", 
        nEtaBins, -4, 0, 100, 0, 1.2);

    events->Draw("(EEEMCal + Barrel)/0.5:eta >> hist2D", "", "goff");

    // Build TGraphErrors from bin-wise projection
    TGraphErrors *graph = new TGraphErrors();
    int pointIndex = 0;

    for (int i = 1; i <= nEtaBins; ++i) {
        TH1D *projY = hist2D->ProjectionY("_py", i, i);
        if (projY->GetEntries() < 2) {
            delete projY;
            continue;
        }

        double etaCenter = hist2D->GetXaxis()->GetBinCenter(i);
        double mean = projY->GetMean();
        double stddev = projY->GetStdDev();
        //Mean error from histogram
        //double meanError = projY->GetMeanError();

        graph->SetPoint(pointIndex, etaCenter, mean);
        graph->SetPointError(pointIndex, 0, stddev);
        ++pointIndex;

        delete projY;
    }

    //Normalization
    hist2D->Scale(1.0 / hist2D->GetMaximum());

    // Style the graph
    graph->SetLineColor(kBlack);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlack);

    // Draw canvas
    TCanvas *c1 = new TCanvas("c1", "Energy reconstruction of 500 MeV e- (EEEMCal at z=175 cm)", 900, 600);
    c1->SetGrid();
    gPad->SetLogz();
    gStyle->SetPalette(104);

    hist2D->Draw("colz");
    hist2D->SetMinimum(1e-3);
    hist2D->SetStats(0);

    graph->Draw("P same");
}