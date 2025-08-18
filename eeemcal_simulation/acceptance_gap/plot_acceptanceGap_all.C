void plot_acceptanceGap_all() {
    std::vector<int> distances = {174, 175, 176, 178, 180, 200};
    const int nEtaBins = 100;

    std::vector<TGraphErrors*> graphs;
    // std::vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta-4, kOrange+7, kMagenta+3};

    // std::vector<int> colors = {
    // kAzure + 1,     // Soft, modern blue
    // kOrange + 7,    // Warm mustard orange
    // kGreen + 3,     // Fresh teal-green
    // kMagenta - 4,   // Muted violet-pink
    // kGray + 2,      // Modern neutral gray
    // kSpring + 5     // Light mint-green
    // };

    std::vector<int> colors = {
    kBlack,         // Classic black for primary data or fit
    kRed,           // Standard red – attention-grabbing
    kBlue,          // Clean blue – easy to distinguish
    kGreen + 2,     // Bright green – common for overlays
    kMagenta,       // Standard magenta – contrasts well
    kCyan + 2       // Soft cyan – good for distinction
    };

    TLegend *legend = new TLegend(0.65, 0.7, 0.88, 0.88);

    for (size_t idx = 0; idx < distances.size(); ++idx) {
        int d = distances[idx];
        TString filename = Form("eicrecon_%dcm_500MeV_90to178deg.root", d);
        TFile *file = TFile::Open(filename);
        if (!file || file->IsZombie()) {
            std::cerr << "Could not open file: " << filename << std::endl;
            continue;
        }

        TTree *events = (TTree*)file->Get("events");

        events->SetAlias("eta", 
            "0.5*log((sqrt(MCParticles.momentum.x[0]*MCParticles.momentum.x[0] + MCParticles.momentum.y[0]*MCParticles.momentum.y[0] + MCParticles.momentum.z[0]*MCParticles.momentum.z[0]) + MCParticles.momentum.z[0]) / (sqrt(MCParticles.momentum.x[0]*MCParticles.momentum.x[0] + MCParticles.momentum.y[0]*MCParticles.momentum.y[0] + MCParticles.momentum.z[0]*MCParticles.momentum.z[0]) - MCParticles.momentum.z[0]))");

        //events-SetAlias("eta", "0.5*log((sqrt(ReconstructedChargedParticles.momentum.x[0]*ReconstructedChargedParticles.momentum.x[0] + ReconstructedChargedParticles.momentum.y[0]*ReconstructedChargedParticles.momentum.y[0] + ReconstructedChargedParticles.momentum.z[0]*ReconstructedChargedParticles.momentum.z[0]) + ReconstructedChargedParticles.momentum.z[0]) / (sqrt(ReconstructedChargedParticles.momentum.x[0]*ReconstructedChargedParticles.momentum.x[0] + ReconstructedChargedParticles.momentum.y[0]*ReconstructedChargedParticles.momentum.y[0] + ReconstructedChargedParticles.momentum.z[0]*ReconstructedChargedParticles.momentum.z[0]) - ReconstructedChargedParticles.momentum.z[0]))");

        events->SetAlias("eta", "0.5*log((sqrt(ReconstructedChargedParticles.momentum.x*ReconstructedChargedParticles.momentum.x + ReconstructedChargedParticles.momentum.y*ReconstructedChargedParticles.momentum.y + ReconstructedChargedParticles.momentum.z*ReconstructedChargedParticles.momentum.z) + ReconstructedChargedParticles.momentum.z) / (sqrt(ReconstructedChargedParticles.momentum.x*ReconstructedChargedParticles.momentum.x + ReconstructedChargedParticles.momentum.y*ReconstructedChargedParticles.momentum.y + ReconstructedChargedParticles.momentum.z*ReconstructedChargedParticles.momentum.z) - ReconstructedChargedParticles.momentum.z))");

        events->SetAlias("Barrel", "Sum$(EcalBarrelScFiRecHits.energy)");
        events->SetAlias("EEEMCal", "Sum$(EcalEndcapNRecHits.energy)");

        TH2F *hist2D = new TH2F("hist2D", "", nEtaBins, -4, 0, 100, 0, 1.2);
        events->Draw("(EEEMCal + Barrel)/0.5:eta >> hist2D", "", "goff");

        TGraphErrors *graph = new TGraphErrors();
        int pointIndex = 0;

        for (int i = 1; i <= nEtaBins; ++i) {
            double etaCenter = hist2D->GetXaxis()->GetBinCenter(i);
            double etaOffset = 0.003 * idx;
            if (etaCenter < -1.8 || etaCenter > -1.55) continue;

            TH1D *projY = hist2D->ProjectionY("_py", i, i);
            if (projY->GetEntries() < 2) {
                delete projY;
                continue;
            }

            double mean = projY->GetMean();
            double stddev = projY->GetStdDev();

            graph->SetPoint(pointIndex, etaCenter+etaOffset, mean);
            graph->SetPointError(pointIndex, 0, stddev);
            pointIndex++;

            delete projY;
        }

        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(colors[idx]);
        graph->SetLineColor(colors[idx]);
        graph->SetLineWidth(2);
        graph->SetMarkerSize(1.0);
        graph->SetTitle(Form("Distance = %d cm", d));

        legend->AddEntry(graph, Form("%d cm", d), "p");
        graphs.push_back(graph);

        file->Close();
    }

    // Plot
    TCanvas *c = new TCanvas("c", "Energy Response vs Eta for Various Distances", 900, 700);
    c->SetGrid();

    for (size_t i = 0; i < graphs.size(); ++i) {
        if (i == 0){
            graphs[i]->SetMinimum(0.0); 
            graphs[i]->SetMaximum(1.2);  
            graphs[i]->Draw("AP");
            graphs[i]->GetXaxis()->SetLimits(-1.8, -1.55);
        }
        else{
            graphs[i]->Draw("P SAME");
        }
    }

    graphs[0]->SetTitle("Energy reconstruction of 500 MeV e- for EEEMCal at different z;#eta_{thrown}; #Sigma E_{dep}/p_{thrown}");

    legend->Draw();

    //c->SaveAs("eta_response_range_-2_to_-1_vs_distance_allOnSamePlot.png");
}
