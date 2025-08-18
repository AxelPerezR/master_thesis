#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

typedef unsigned int uint;

void WriteRMSCSV(const std::string& filename,
                 const std::vector<double>& channel,
                 const std::vector<double>& rms_no,
                 const std::vector<double>& rms_cab,
                 const std::vector<double>& rms_cab_hv){
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    file << "index_pb,PB_INA,channel_shihai,rms_no,rms_cab,rms_cab_hv\n";
    size_t n = channel.size();
    int index_asic = 0;
    for (size_t i = 0; i < n; ++i) {
        index_asic = i % 72; // Assuming 144 channels per ASIC
        file << i << "," << index_asic << "," << channel[i] << "," << rms_no[i] << "," << rms_cab[i] << ","
             << std::fixed << std::setprecision(3) << rms_cab_hv[i] << "\n";
    }

    file.close();
    std::cout << " \n CSV file written: " << filename << std::endl;
}

void pedestals_all_sample(int run_no_cables = 32, int run_cables = 33, int run_cables_hv = 34)
{
    // Tristan channels to ignore
    std::vector<int> ignored_channels = {8, 17, 26, 35, 44, 53, 62, 71, 80, 89, 98, 107, 116, 125, 134, 143};
    std::vector<int> non_used_channels_0to151 = {9, 18, 28, 37, 47, 56, 66, 75, 85, 94, 104, 113, 123, 132, 142, 151}; 

    const int num_kcu = 3;
    const int channels_per_kcu = 144;
    const int num_channels = num_kcu * channels_per_kcu;

    // Shihai mapping
    std::vector<int> shihai_channel_map;
    for (int ch = 1; ch <= 151; ++ch)
    {
        if (ch == 19 || ch == 38 || ch == 57 || ch == 76 || ch == 95 || ch == 114 || ch == 133)
            continue;
        shihai_channel_map.push_back(ch);
    }

    std::vector<int> fpga_ids = {208, 209, 210};
    std::vector<std::string> protoboards = {"01", "008", "06"};

    const TString input_dir = "/home/axelperezr/eic/pedestals/data/";
    const TString output_dir = "/home/axelperezr/eic/pedestals/plots/";

    TString run_no_cables_id = Form("run%03d", run_no_cables);
    TString run_cables_id = Form("run%03d", run_cables);
    TString run_cables_hv_id = Form("run%03d", run_cables_hv);

    TFile *f_no_cables = TFile::Open(input_dir + run_no_cables_id + ".root");
    TFile *f_cables = TFile::Open(input_dir + run_cables_id + ".root");
    TFile *f_cables_hv = TFile::Open(input_dir + run_cables_hv_id + ".root");

    if (!f_no_cables || f_no_cables->IsZombie() || !f_cables || f_cables->IsZombie() || !f_cables_hv || f_cables_hv->IsZombie())
    {
        std::cerr << "Cannot open input files." << std::endl;
        return;
    }

    TTree *t_no_cables = (TTree *)f_no_cables->Get("events");
    TTree *t_cables = (TTree *)f_cables->Get("events");
    TTree *t_cables_hv = (TTree *)f_cables_hv->Get("events");
    if (!t_no_cables || !t_cables || !t_cables_hv)
    {
        std::cerr << "Cannot find TTree 'events' in files." << std::endl;
        return;
    }

    uint adc_no_cables[num_channels][20];
    uint adc_cables[num_channels][20];
    uint adc_cables_hv[num_channels][20];

    uint hit_ped_no_cables[num_channels];
    uint hit_ped_cables[num_channels];
    uint hit_ped_cables_hv[num_channels];


    t_no_cables->SetBranchAddress("adc", adc_no_cables);
    t_cables->SetBranchAddress("adc", adc_cables);
    t_cables_hv->SetBranchAddress("adc", adc_cables_hv);

    t_no_cables->SetBranchAddress("hit_pedestal", hit_ped_no_cables);
    t_cables->SetBranchAddress("hit_pedestal", hit_ped_cables);
    t_cables_hv->SetBranchAddress("hit_pedestal", hit_ped_cables_hv);

    const int num_samples = 20;
    Long64_t nEntries = std::min({t_no_cables->GetEntries(), t_cables->GetEntries(), t_cables_hv->GetEntries()});

    std::vector<std::vector<std::vector<uint>>> max_adc_no_cables(nEntries, std::vector<std::vector<uint>>(num_channels, std::vector<uint>(num_samples, 0)));
    std::vector<std::vector<std::vector<uint>>> max_adc_cables(nEntries, std::vector<std::vector<uint>>(num_channels, std::vector<uint>(num_samples, 0)));
    std::vector<std::vector<std::vector<uint>>> max_adc_cables_hv(nEntries, std::vector<std::vector<uint>>(num_channels, std::vector<uint>(num_samples, 0)));

    std::vector<std::vector<std::vector<int>>> adc_minus_pedestal_no_cables(nEntries, std::vector<std::vector<int>>(num_channels, std::vector<int>(num_samples, 0)));
    std::vector<std::vector<std::vector<int>>> adc_minus_pedestal_cables(nEntries, std::vector<std::vector<int>>(num_channels, std::vector<int>(num_samples, 0)));
    std::vector<std::vector<std::vector<int>>> adc_minus_pedestal_cables_hv(nEntries, std::vector<std::vector<int>>(num_channels, std::vector<int>(num_samples, 0)));

    std::vector<TH1I*> channel_histograms_no_cables(num_channels);
    std::vector<TH1I*> channel_histograms_cables(num_channels);
    std::vector<TH1I*> channel_histograms_cables_hv(num_channels);

    for (int ch = 0; ch < num_channels; ++ch)
    {
        channel_histograms_no_cables[ch] = new TH1I(Form("h_no_cables_ch_%d", ch), Form("Channel %d - No Cables", ch), 200, 0, 200);
        channel_histograms_cables[ch] = new TH1I(Form("h_cables_ch_%d", ch), Form("Channel %d - Cables", ch), 200, 0, 200);
        channel_histograms_cables_hv[ch] = new TH1I(Form("h_cables_hv_ch_%d", ch), Form("Channel %d - Cables HV", ch), 200, 0, 200);
    }

    std::vector<TH2I*> adc_samples_histograms_no_cables(num_channels);
    std::vector<TH2I*> adc_samples_histograms_cables(num_channels);
    std::vector<TH2I*> adc_samples_histograms_cables_hv(num_channels);

    for (int ch = 0; ch < num_channels; ++ch) {
        adc_samples_histograms_no_cables[ch] = new TH2I(Form("adc_no_cables_ch%d", ch), "ADC vs Sample", num_samples, 0, num_samples, 1024, 0, 1024);
        adc_samples_histograms_cables[ch] = new TH2I(Form("adc_cables_ch%d", ch), "ADC vs Sample", num_samples, 0, num_samples, 1024, 0, 1024);
        adc_samples_histograms_cables_hv[ch] = new TH2I(Form("adc_cables_hv_ch%d", ch), "ADC vs Sample", num_samples, 0, num_samples, 1024, 0, 1024);
    }


    std::vector<std::vector<std::vector<int>>> adc_event_no_cables(nEntries, std::vector<std::vector<int>>(num_channels, std::vector<int>(num_samples)));
    std::vector<std::vector<std::vector<int>>> adc_event_cables(nEntries, std::vector<std::vector<int>>(num_channels, std::vector<int>(num_samples)));
    std::vector<std::vector<std::vector<int>>> adc_event_cables_hv(nEntries, std::vector<std::vector<int>>(num_channels, std::vector<int>(num_samples)));

    // Histogram for hit_max distribution for channel 0 on Protoboard 008
    int target_fpga = 0;
    int local_ch = 4;
    int global_ch = target_fpga * channels_per_kcu + local_ch;

    TH1F *h_adc_distribution_no_cables = new TH1F("h_adc_distribution_no_cables", "Pedestals distribution for a representative high RMS Channel (Channel 5 on Protoboard 01);ADC Value;Frequency", 100, 0, 120);
    TH1F *h_adc_distribution_cables = new TH1F("h_adc_distribution_cables", "Pedestals distribution for a representative high RMS Channel (Channel 5 on Protoboard 01);ADC Value;Frequency", 100, 0, 120);
    TH1F *h_adc_distribution_cables_hv = new TH1F("h_adc_distribution_cables_hv", "Pedestals distribution for a representative high RMS Channel (Channel 5 on Protoboard 01);ADC Value;Frequency", 100, 0, 120);

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        t_no_cables->GetEntry(i);
        t_cables->GetEntry(i);
        t_cables_hv->GetEntry(i);

        for (int sample = 0; sample < num_samples; ++sample)
        {
            h_adc_distribution_no_cables->Fill(adc_no_cables[global_ch][sample] - hit_ped_no_cables[global_ch]);
            h_adc_distribution_cables->Fill(adc_cables[global_ch][sample] - hit_ped_cables[global_ch]);
            h_adc_distribution_cables_hv->Fill(adc_cables_hv[global_ch][sample] - hit_ped_cables_hv[global_ch]);
        }

        for (int ch = 0; ch < num_channels; ++ch)
        {
            for (int sample = 0; sample < num_samples; ++sample)
            {
                max_adc_no_cables[i][ch][sample] = adc_no_cables[ch][sample];
                adc_minus_pedestal_no_cables[i][ch][sample] = adc_no_cables[ch][sample] - hit_ped_no_cables[ch];

                max_adc_cables[i][ch][sample] = adc_cables[ch][sample];
                adc_minus_pedestal_cables[i][ch][sample] = adc_cables[ch][sample] - hit_ped_cables[ch];

                max_adc_cables_hv[i][ch][sample] = adc_cables_hv[ch][sample];
                adc_minus_pedestal_cables_hv[i][ch][sample] = adc_cables_hv[ch][sample] - hit_ped_cables_hv[ch];

                int content_no_cables = adc_no_cables[ch][sample];
                int content_cables = adc_cables[ch][sample];
                int content_cables_hv = adc_cables_hv[ch][sample];

                channel_histograms_no_cables[ch]->Fill(content_no_cables);
                channel_histograms_cables[ch]->Fill(content_cables);
                channel_histograms_cables_hv[ch]->Fill(content_cables_hv);

                adc_samples_histograms_no_cables[ch]->Fill(sample, content_no_cables);
                adc_samples_histograms_cables[ch]->Fill(sample, content_cables);
                adc_samples_histograms_cables_hv[ch]->Fill(sample, content_cables_hv);

                adc_event_no_cables[i][ch][sample] = adc_no_cables[ch][sample];
                adc_event_cables[i][ch][sample] = adc_cables[ch][sample];
                adc_event_cables_hv[i][ch][sample] = adc_cables_hv[ch][sample];
            }
        }
    }

    //std::cout << "size of adc_minus_pedestal_no_cables: " << adc_minus_pedestal_no_cables.size() << std::endl;
    //std::cout << "size of adc_minus_pedestal_no_cables channel 0: " << adc_minus_pedestal_no_cables[0] << std::endl;

    TCanvas *c_adc_freq = new TCanvas("c_adc_freq", "ADC Distribution", 800, 800);
    c_adc_freq->SetGrid();

    //h_adc_distribution_cables->SetLineColor(kBlue + 1);
    // h_adc_distribution_cables->SetLineWidth(2);
    h_adc_distribution_no_cables->SetStats(0); // Disable stats box

    h_adc_distribution_cables->SetFillColor(kBlue); // 0.0 = fully transparent, 1.0 = fully opaque
    h_adc_distribution_cables->SetLineColor(kBlue - 9); // 0.0 = fully transparent, 1.0 = fully opaque
    h_adc_distribution_cables->SetFillStyle(3004); // 0.0 = fully transparent, 1.0 = fully opaque
    h_adc_distribution_cables->SetLineWidth(2);

    h_adc_distribution_no_cables->SetLineColor(kRed);
    h_adc_distribution_no_cables->SetFillStyle(3004);
    h_adc_distribution_no_cables->SetFillColor(kRed - 9);
    h_adc_distribution_no_cables->SetLineWidth(2);

    h_adc_distribution_cables_hv->SetFillStyle(3004);
    h_adc_distribution_cables_hv->SetLineColor(kGreen + 2);
    h_adc_distribution_cables_hv->SetFillColor(kGreen + 1);
    h_adc_distribution_cables_hv->SetLineWidth(2);

    h_adc_distribution_no_cables->Draw("HIST");
    h_adc_distribution_cables->Draw("HIST SAME");
    h_adc_distribution_cables_hv->Draw("HIST SAME");

    // Place the legend in the top right corner
    TLegend *leg4 = new TLegend(0.45, 0.85, 0.85, 0.70);
    leg4->SetTextSize(0.035); // Adjust text size for better visibility
    leg4->AddEntry(h_adc_distribution_no_cables, Form("No Cables (run %d)", run_no_cables), "f");
    leg4->AddEntry(h_adc_distribution_cables, Form("Cables (run %d)", run_cables), "f");
    leg4->AddEntry(h_adc_distribution_cables_hv, Form("Cables HV (run %d)", run_cables_hv), "f");
    //leg4->AddEntry((TObject*)0, Form("Total Channels Affected: 169"), "");
    leg4->Draw();

    // Add a text box with the total high RMS channels
    int total_low_rms = 169;
    int total_high_rms = 130;
    int total_average_rms = 5;
    int total_slightly_higher_rms = 80;
    TLatex *text = new TLatex(0.4, 0.65, Form("Total High RMS Channels: %d", total_high_rms));
    text->SetNDC();
    text->SetTextSize(0.035);
    text->Draw();

    c_adc_freq->SaveAs(output_dir + "adc_distribution_channel34_pb008.png");

    auto style_histogram = [](TH1 *h, const char *title, const char *ytitle)
    {
        h->SetTitle(title);
        h->GetXaxis()->SetTitle("Channel"); // It was, "Shihai Channel"
        h->GetYaxis()->SetTitle(ytitle);
        h->GetXaxis()->SetTitleSize(0.045);
        h->GetYaxis()->SetTitleSize(0.045);
        h->GetXaxis()->SetLabelSize(0.04);
        h->GetYaxis()->SetLabelSize(0.04);
        h->GetXaxis()->SetTitleOffset(1.2);
        h->GetYaxis()->SetTitleOffset(1.3);
        h->LabelsOption("h");
        h->SetLineWidth(2);
        h->SetStats(0); // Disable stats box

        h->GetXaxis()->SetTitleOffset(1.0); // Move X-axis title closer
        h->GetYaxis()->SetTitleOffset(0.6); // Move Y-axis title closer
    };

    for (int fpga = 0; fpga < num_kcu; ++fpga)
    {
        TString label = Form("fpga%d_%d_%s", fpga, fpga_ids[fpga], protoboards[fpga].c_str());

        TString out1 = output_dir + run_no_cables_id + "_vs_" + run_cables_id + "_" + label + "_max_adc_vs_channel.png";
        TString out2 = output_dir + run_cables_id + "_vs_" + run_cables_hv_id + "_" + label + "_adc_minus_pedestal_vs_channel.png";

        int ch_start = fpga * channels_per_kcu;
        int ch_end = ch_start + channels_per_kcu;

        // Always use 144 bins for all FPGAs
        TCanvas *c1 = new TCanvas(Form("c1_fpga%d", fpga), "ADC", 2000, 800);
        c1->SetGrid();
        TH1I *h1_no_cables = new TH1I(Form("h1_no_cables_protoboard%d", fpga), Form("ADC (Protoboard %s - No Cables)", protoboards[fpga].c_str()), channels_per_kcu, 0, channels_per_kcu);
        TH1I *h1_cables = new TH1I(Form("h1_cables_protoboard%d", fpga), Form("ADC (Protoboard %s - Cables)", protoboards[fpga].c_str()), channels_per_kcu, 0, channels_per_kcu);
        TH1I *h1_cables_hv = new TH1I(Form("h1_cables_hv_protoboard%d", fpga), Form("ADC (Protoboard %s - Cables HV)", protoboards[fpga].c_str()), channels_per_kcu, 0, channels_per_kcu);

        h1_no_cables->Sumw2();  // Enables storage of bin errors
        h1_cables->Sumw2();     // Enables storage of bin errors
        h1_cables_hv->Sumw2();  // Enables storage of bin errors

        // TProfile* p_no_cables = new TProfile(Form("p_no_cables_fpga%d", fpga), "Mean & StdDev per Channel", num_channels, 0, num_channels);
        // TProfile* p_cables = new TProfile(Form("p_cables_fpga%d", fpga), "Mean & StdDev per Channel", num_channels, 0, num_channels);
        // TProfile* p_cables_hv = new TProfile(Form("p_cables_hv_fpga%d", fpga), "Mean & StdDev per Channel", num_channels, 0, num_channels);

        TCanvas *c_rms = new TCanvas(Form("c_rms_fpga%d", fpga), "RMS per Channel", 2000, 800);
        c_rms->SetGrid();

        std::vector<double> x_values; // X values are the same for all
        std::vector<double> rms_no_cables;
        std::vector<double> rms_cables;
        std::vector<double> rms_cables_hv;

        // TH1F *h_rms_no_cables = new TH1F(Form("h_rms_distribution_no_cables_fpga%d", fpga), "RMS Frequency Across Channels (No Cables);RMS Value;Number of Channels", 50, 0, 30);
        // TH1F *h_rms_cables = new TH1F(Form("h_rms_distribution_cables_fpga%d", fpga), "RMS Frequency Across Channels (Cables);RMS Value;Number of Channels", 50, 0, 30);
        // TH1F *h_rms_cables_hv = new TH1F(Form("h_rms_distribution_cables_hv_fpga%d", fpga), "RMS Frequency Across Channels (Cables HV);RMS Value;Number of Channels", 50, 0, 30);

        for (int local_ch = 0; local_ch < channels_per_kcu; ++local_ch)
        {
            // if (std::find(ignored_channels.begin(), ignored_channels.end(), local_ch) != ignored_channels.end())
            //     continue;

            int global_ch = ch_start + local_ch;
            int bin = local_ch + 1;

            double sum = 0.0;
            int count = 0;

            for (Long64_t i = 0; i < nEntries; ++i)
            {
                for (int sample = 0; sample < num_samples; ++sample)
                {
                    sum += max_adc_no_cables[i][global_ch][sample];
                    ++count;
                }
            }

            double mean_adc_no_cables = (count > 0) ? (sum / count) : 0.0;

            double variance_sum = 0.0;
            for (Long64_t i = 0; i < nEntries; ++i)
            {
                for (int sample = 0; sample < num_samples; ++sample)
                {
                    double val = max_adc_no_cables[i][global_ch][sample];
                    variance_sum += (val - mean_adc_no_cables) * (val - mean_adc_no_cables);
                }
            }

            double stddev_no_cables = (count > 0) ? std::sqrt(variance_sum / count) : 0.0;

            // std::cout << "Global Channel: " << global_ch << ", Local Channel: " << local_ch << ", Mean ADC No Cables: " << mean_adc_no_cables << ", StdDev No Cables: " << stddev_no_cables << std::endl;

            sum = 0.0;
            count = 0;

            for (Long64_t i = 0; i < nEntries; ++i)
            {
                for (int sample = 0; sample < num_samples; ++sample)
                {
                    sum += max_adc_cables[i][global_ch][sample];
                    ++count;
                }
            }

            double mean_adc_cables = (count > 0) ? (sum / count) : 0.0;

            variance_sum = 0.0;

            for (Long64_t i = 0; i < nEntries; ++i)
            {
                for (int sample = 0; sample < num_samples; ++sample)
                {
                    double val = max_adc_cables[i][global_ch][sample];
                    variance_sum += (val - mean_adc_cables) * (val - mean_adc_cables);
                }
            }

            double stddev_cables = (count > 0) ? std::sqrt(variance_sum / count) : 0.0;

            sum = 0.0;
            count = 0;

            for (Long64_t i = 0; i < nEntries; ++i)
            {
                for (int sample = 0; sample < num_samples; ++sample)
                {
                    sum += max_adc_cables_hv[i][global_ch][sample];
                    ++count;
                }
            }

            double mean_adc_cables_hv = (count > 0) ? (sum / count) : 0.0;

            variance_sum = 0.0;
            for (Long64_t i = 0; i < nEntries; ++i)
            {
                for (int sample = 0; sample < num_samples; ++sample)
                {
                    double val = max_adc_cables_hv[i][global_ch][sample];
                    variance_sum += (val - mean_adc_cables_hv) * (val - mean_adc_cables_hv);
                }
            }

            double stddev_cables_hv = (count > 0) ? std::sqrt(variance_sum / count) : 0.0;

            h1_no_cables->SetBinContent(bin, mean_adc_no_cables);
            h1_no_cables->SetBinError(bin, stddev_no_cables);
            h1_cables->SetBinContent(bin, mean_adc_cables);
            h1_cables->SetBinError(bin, stddev_cables);
            h1_cables_hv->SetBinContent(bin, mean_adc_cables_hv);
            h1_cables_hv->SetBinError(bin, stddev_cables_hv);

            int shihai_label = shihai_channel_map[local_ch];

            // Check if histograms are empty
            // std::cout << "Entries for Channel " << global_ch << ": " << channel_histograms_no_cables[global_ch]->GetEntries() << std::endl;

            double rms_no = channel_histograms_no_cables[global_ch]->GetRMS();
            //printf("RMS for Channel %d (No Cables): %.2f\n", shihai_label, rms_no);
            double rms_cab = channel_histograms_cables[global_ch]->GetRMS();
            //printf("RMS for Channel %d (Cables): %.2f\n", shihai_label, rms_cab);
            double rms_hv = channel_histograms_cables_hv[global_ch]->GetRMS();
            //printf("RMS for Channel %d (Cables HV): %.2f\n", shihai_label, rms_hv);

            // h_rms_no_cables->Fill(rms_no);
            // h_rms_cables->Fill(rms_cab);
            // h_rms_cables_hv->Fill(rms_hv);

            // Only label every 10 bins
            if (bin % 10 == 1) {
                h1_no_cables->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_label));
                h1_cables->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_label));
                h1_cables_hv->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_label));
            }
            else {
                h1_no_cables->GetXaxis()->SetBinLabel(bin, "");
                h1_cables->GetXaxis()->SetBinLabel(bin, "");
                h1_cables_hv->GetXaxis()->SetBinLabel(bin, "");
            }
            
            //index_values.push_back(global_ch);
            x_values.push_back(shihai_label);
            rms_no_cables.push_back(rms_no);
            rms_cables.push_back(rms_cab);
            rms_cables_hv.push_back(rms_hv);
        }

        // Write RMS values to CSV
        std::string csv_filename = "RMS_values_fpga" + std::to_string(fpga) + ".csv";
        WriteRMSCSV(csv_filename, x_values, rms_no_cables, rms_cables, rms_cables_hv);

        h1_no_cables->SetLineColor(kRed + 1);
        h1_cables->SetLineColor(kBlue + 1);
        h1_cables_hv->SetLineColor(kGreen + 1);
        style_histogram(h1_no_cables, Form("Pedestal distribution (Protoboard %s)", protoboards[fpga].c_str()), "ADC");
        h1_no_cables->SetMaximum(150);
        //h1_no_cables->SetMaximum(300); // Set maximum for better visibility
        //h_rms->SetMaximum(6); // Set maximum for RMS histogram

        c1->cd();  // Make sure you are in c1
        h1_no_cables->Draw("HIST");
        h1_no_cables->Draw("E SAME"); // Draw the error bars for no cables
        h1_cables->Draw("HIST SAME");
        h1_cables->Draw("E SAME"); // Draw the error bars for cables
        h1_cables_hv->Draw("HIST SAME");
        h1_cables_hv->Draw("E SAME"); // Draw the error bars for cables HV

        TLegend *leg1 = new TLegend(0.15, 0.75, 0.28, 0.88);
        leg1->AddEntry(h1_no_cables, Form("No Cables (run %d)", run_no_cables), "l");
        leg1->AddEntry(h1_cables, Form("Cables (run %d)", run_cables), "l");
        leg1->AddEntry(h1_cables_hv, Form("Cables HV (run %d)", run_cables_hv), "l");
        leg1->Draw();

        // Place a text in the right corner with the number of channels, samples and events
        TLatex *text1 = new TLatex(0.55, 0.75, Form("Channels: 128, Total Samples: %d, Number of Events: %lld", num_samples, nEntries));
        text1->SetNDC();
        text1->SetTextSize(0.035);
        text1->Draw();

        //h_rms->SetLineColor(kBlack);
        //style_histogram(h_rms, "RMS per Channel", "RMS (ADC)");
        //h_rms->Draw("HIST");

        c1->SaveAs(out1);
        //c_rms->SaveAs(output_dir + Form("rms_per_channel_fpga%d.png", fpga));

        TGraph *g_no_cables = new TGraph(x_values.size(), &x_values[0], &rms_no_cables[0]);
        TGraph *g_cables = new TGraph(x_values.size(), &x_values[0], &rms_cables[0]);
        TGraph *g_cables_hv = new TGraph(x_values.size(), &x_values[0], &rms_cables_hv[0]);

        g_no_cables->SetTitle(Form("RMS per Channel - Protoboard %s", protoboards[fpga].c_str()));
        g_no_cables->SetMarkerStyle(20);
        g_no_cables->SetMarkerSize(1);
        g_no_cables->SetMarkerColor(kRed + 1);
        g_no_cables->SetLineColor(kRed + 1);

        g_cables->SetMarkerStyle(21);
        g_cables->SetMarkerSize(1);
        g_cables->SetMarkerColor(kBlue + 1);
        g_cables->SetLineColor(kBlue + 1);

        g_cables_hv->SetMarkerStyle(22);
        g_cables_hv->SetMarkerSize(1);
        g_cables_hv->SetMarkerColor(kGreen + 1);
        g_cables_hv->SetLineColor(kGreen + 1);

        c_rms->cd();  // Make sure you are in c_rms
        g_no_cables->Draw("AP");    // Draw axes and points for the first graph
        g_cables->Draw("P SAME");   // Add points for the second graph
        g_cables_hv->Draw("P SAME");// Add points for the third graph

        // Draw a vertical line at channel 103
        // TLine *line = new TLine(103, 0, 103, 10); // Adjust Y limits as needed
        // line->SetLineColor(kRed);
        // line->SetLineStyle(3); // Dashed line
        // line->Draw("SAME");

        // Draw three horizontal lines at RMS No Cables: 1.79409, RMS Cables: 2.08897, RMS Cables HV: 1.84272
        // TLine *line_no_cables = new TLine(1, 1.79409, 151, 1.79409);
        // line_no_cables->SetLineColor(kRed);
        // line_no_cables->SetLineStyle(2); // Dashed line
        // line_no_cables->Draw("SAME");

        // Compute the mean RMS values
        double mean_rms_no_cables = std::accumulate(rms_no_cables.begin(), rms_no_cables.end(), 0.0) / rms_no_cables.size();
        double mean_rms_cables = std::accumulate(rms_cables.begin(), rms_cables.end(), 0.0) / rms_cables.size();
        double mean_rms_cables_hv = std::accumulate(rms_cables_hv.begin(), rms_cables_hv.end(), 0.0) / rms_cables_hv.size();

        // Draw horizontal lines for mean RMS values
        TLine *line_mean_no_cables = new TLine(1, mean_rms_no_cables, 151, mean_rms_no_cables);
        line_mean_no_cables->SetLineColor(kRed);
        line_mean_no_cables->SetLineStyle(2); // Dashed line
        line_mean_no_cables->SetLineWidth(2);
        line_mean_no_cables->Draw("SAME");

        TLine *line_mean_cables = new TLine(1, mean_rms_cables, 151, mean_rms_cables);
        line_mean_cables->SetLineColor(kBlue);
        line_mean_cables->SetLineStyle(2); // Dashed line
        line_mean_cables->SetLineWidth(2);
        line_mean_cables->Draw("SAME");

        TLine *line_mean_cables_hv = new TLine(1, mean_rms_cables_hv, 151, mean_rms_cables_hv);
        line_mean_cables_hv->SetLineColor(kGreen);
        line_mean_cables_hv->SetLineStyle(2); // Dashed line
        line_mean_cables_hv->SetLineWidth(2);
        line_mean_cables_hv->Draw("SAME");


        // Find the most stable channels for different configurations //Between configurations the change is not significant
        // double percent_tolerance = 0.20; // 20% tolerance
        // int channels_within_tolerance = 0;

        // std::cout << "Stable Channels for Protoboard " << protoboards[fpga] << " (within " << percent_tolerance * 100 << "% tolerance):" << std::endl;

        // std::cout << "Number of channels: " << channels_per_kcu << std::endl;
        // for (int i = 0; i < channels_per_kcu; ++i)
        // {
        //     double tolerance_no_cables = percent_tolerance * rms_no_cables[i];
        //     double tolerance_cables = percent_tolerance * rms_cables[i];
        //     double tolerance_cables_hv = percent_tolerance * rms_cables_hv[i];

        //     bool is_close_no_cables = fabs(rms_no_cables[i] - rms_cables[i]) <= tolerance_no_cables;
        //     bool is_close_cables = fabs(rms_cables[i] - rms_cables_hv[i]) <= tolerance_cables;
        //     bool is_close_cables_hv = fabs(rms_cables_hv[i] - rms_no_cables[i]) <= tolerance_cables_hv;

        //     if (is_close_cables && std::find(non_used_channels_0to151.begin(), non_used_channels_0to151.end(), x_values[i]) == non_used_channels_0to151.end())
        //     {
        //         std::cout << "Channel " << x_values[i] << " is within " << percent_tolerance * 100 << "% tolerance for all configurations." << std::endl;
        //         //std::cout << "Channel index: " << i << std::endl;
        //         std::cout << "RMS No Cables: " << rms_no_cables[i] << ", RMS Cables: " << rms_cables[i] << ", RMS Cables HV: " << rms_cables_hv[i] << std::endl;
        //         //std::cout << "Tolerance No Cables: " << tolerance_no_cables << ", Tolerance Cables: " << tolerance_cables << ", Tolerance Cables HV: " << tolerance_cables_hv << std::endl;
        //         channels_within_tolerance++;
        //     }
        // }

        // std::cout << "Total channels within " << percent_tolerance * 100 << "% of each other: " << channels_within_tolerance << std::endl;

        // std::cout << "Below or equal to the mean RMS values for connected and hv cables for Protoboard " << protoboards[fpga] << ":" << std::endl;

        // int channels_below_mean = 0;

        // for (int i = 0; i < channels_per_kcu; ++i)
        // {
        //     if ((rms_cables[i] <= mean_rms_cables && rms_cables_hv[i] <= mean_rms_cables_hv) && std::find(non_used_channels_0to151.begin(), non_used_channels_0to151.end(), x_values[i]) == non_used_channels_0to151.end())
        //     {
        //         if (x_values[i] < 1 || x_values[i] > 151) continue; // Skip if x_value is out of range
        //         std::cout << x_values[i] << std::endl;
        //         //std::cout << "Channel " << x_values[i] << " is below or equal to the mean RMS values." << std::endl;
        //         //std::cout << "RMS Cables: " << rms_cables[i] << ", RMS Cables HV: " << rms_cables_hv[i] << std::endl;
        //         channels_below_mean++;
        //     }
        // }

        // std::cout << "Total channels below or equal to the mean RMS values for connected and hv cables: " << channels_below_mean << std::endl;

        // std::cout << "Above the mean RMS values for connected and hv cables for Protoboard " << protoboards[fpga] << ":" << std::endl;

        // int channels_above_mean = 0;

        // for (int i = 0; i < channels_per_kcu; ++i)
        // {
        //     if ((rms_cables[i] > mean_rms_cables || rms_cables_hv[i] > mean_rms_cables_hv) && std::find(non_used_channels_0to151.begin(), non_used_channels_0to151.end(), x_values[i]) == non_used_channels_0to151.end())
        //     {
        //         if (x_values[i] < 1 || x_values[i] > 151) continue; // Skip if x_value is out of range
        //         std::cout << x_values[i] << std::endl;
        //         //std::cout << "Channel " << x_values[i] << " is above the mean RMS values." << std::endl;
        //         channels_above_mean++;
        //     }
        // }

        // std::cout << "Total channels above the mean RMS values for connected and hv cables: " << channels_above_mean << std::endl;

        // int valid_channels = 0;

        // for (int i = 0; i < channels_per_kcu; ++i)
        // {
        //     if (std::find(non_used_channels_0to151.begin(), non_used_channels_0to151.end(), x_values[i]) != non_used_channels_0to151.end())
        //     {
        //         // Skip the channels that are not used
        //         continue;
        //     }
        //     valid_channels++;
        //     std::cout << "Channel " << x_values[i] << ": RMS No Cables = " << rms_no_cables[i] << ", RMS Cables = " << rms_cables[i] << ", RMS Cables HV = " << rms_cables_hv[i] << std::endl;
        // }

        // std::cout << "Total valid channels: " << valid_channels << std::endl;

        TLegend *leg = new TLegend(0.50, 0.75, 0.85, 0.88);
        // Increase the legend size to accommodate longer text
        leg->SetTextSize(0.035);
        leg->AddEntry(g_no_cables, Form("No Cables (run %d), Mean RMS = %.2f", run_no_cables, mean_rms_no_cables), "p");
        leg->AddEntry(g_cables, Form("Cables (run %d), Mean RMS = %.2f", run_cables, mean_rms_cables), "p");
        leg->AddEntry(g_cables_hv, Form("Cables HV (run %d), Mean RMS = %.2f", run_cables_hv, mean_rms_cables_hv), "p");
        leg->Draw();

        // Set Y Limits for better visibility
        g_no_cables->GetYaxis()->SetRangeUser(0, 20);
        g_cables->GetYaxis()->SetRangeUser(0, 20);
        g_cables_hv->GetYaxis()->SetRangeUser(0, 20);

        // Set X axis limits
        g_no_cables->GetXaxis()->SetLimits(1, 151);
        g_cables->GetXaxis()->SetLimits(1, 151);
        g_cables_hv->GetXaxis()->SetLimits(1, 151);

        // Set axis titles
        g_no_cables->GetXaxis()->SetTitle("Channel");
        g_no_cables->GetYaxis()->SetTitle("RMS (ADC)");
        g_cables->GetXaxis()->SetTitle("Channel");
        g_cables->GetYaxis()->SetTitle("RMS (ADC)");
        g_cables_hv->GetXaxis()->SetTitle("Channel");
        g_cables_hv->GetYaxis()->SetTitle("RMS (ADC)");

        c_rms->SaveAs(output_dir + Form("rms_per_channel_fpga%d.png", fpga));

        // // RMS distribution histograms
        // TCanvas *c_rms_freq = new TCanvas("c_rms_freq", "RMS Frequency Distribution", 1200, 800);
        // c_rms_freq->SetGrid();
        // h_rms_no_cables->SetTitle(Form("RMS Frequency Distribution - Protoboard %s", protoboards[fpga].c_str()));

        // h_rms_no_cables->SetLineColor(kRed + 1);
        // h_rms_no_cables->SetLineWidth(3);
        // h_rms_cables->SetLineColor(kBlue + 1);
        // h_rms_cables->SetLineWidth(3);
        // h_rms_cables_hv->SetLineColor(kGreen + 1);
        // h_rms_cables_hv->SetLineWidth(3);

        // h_rms_no_cables->Draw("HIST");
        // h_rms_cables->Draw("HIST SAME");
        // h_rms_cables_hv->Draw("HIST SAME");

        // TLegend *leg3 = new TLegend(0.15, 0.75, 0.35, 0.88);
        // leg3->AddEntry(h_rms_no_cables, "No Cables", "l");
        // leg3->AddEntry(h_rms_cables, "Cables", "l");
        // leg3->AddEntry(h_rms_cables_hv, "Cables HV", "l");
        // leg3->Draw();


        //c_rms_freq->SaveAs(output_dir + Form("rms_frequency_comparison_fpga%d.png", fpga));

        // Max ADC - Pedestal
        TCanvas *c2 = new TCanvas(Form("c2_fpga%d", fpga), "ADC - Calibration Pedestal", 2000, 800);
        c2->SetGrid();
        TH1I *h2_no_cables = new TH1I(Form("h2_no_cables_protoboard%d", fpga), Form("ADC - Calibration Pedestal (Protoboard %s - No Cables)", protoboards[fpga].c_str()), channels_per_kcu, 0, channels_per_kcu);
        TH1I *h2_cables = new TH1I(Form("h2_cables_protoboard%d", fpga), Form("ADC - Calibration Pedestal (Protoboard %s - Cables)", protoboards[fpga].c_str()), channels_per_kcu, 0, channels_per_kcu);
        TH1I *h2_cables_hv = new TH1I(Form("h2_cables_hv_protoboard%d", fpga), Form("ADC - Calibration Pedestal (Protoboard %s - Cables HV)", protoboards[fpga].c_str()), channels_per_kcu, 0, channels_per_kcu);

        // c2->SetLeftMargin(-0.09);   // Allow enough space for Y-axis title
        c2->SetBottomMargin(0.12); // Increase bottom margin for X-axis labels

        for (int local_ch = 0; local_ch < channels_per_kcu; ++local_ch)
        {
            if (std::find(ignored_channels.begin(), ignored_channels.end(), local_ch) != ignored_channels.end())
                continue;

            int global_ch = ch_start + local_ch;
            int bin = local_ch + 1;

            double sum = 0.0;
            int count = 0;

            for (Long64_t i = 0; i < nEntries; ++i)
            {
                for (int sample = 0; sample < num_samples; ++sample)
                {
                    sum += adc_minus_pedestal_no_cables[i][global_ch][sample];
                    ++count;
                }
            }

            double mean_adc_minus_pedestal_no_cables = (count > 0) ? (sum / count) : 0.0;

            sum = 0.0;
            count = 0;

            for (Long64_t i = 0; i < nEntries; ++i)
            {
                for (int sample = 0; sample < num_samples; ++sample)
                {
                    sum += adc_minus_pedestal_cables[i][global_ch][sample];
                    ++count;
                }
            }

            double mean_adc_minus_pedestal_cables = (count > 0) ? (sum / count) : 0.0;

            sum = 0.0;
            count = 0;

            for (Long64_t i = 0; i < nEntries; ++i)
            {
                for (int sample = 0; sample < num_samples; ++sample)
                {
                    sum += adc_minus_pedestal_cables_hv[i][global_ch][sample];
                    ++count;
                }
            }

            double mean_adc_minus_pedestal_cables_hv = (count > 0) ? (sum / count) : 0.0;

            h2_no_cables->SetBinContent(bin, mean_adc_minus_pedestal_no_cables);
            h2_cables->SetBinContent(bin, mean_adc_minus_pedestal_cables);
            h2_cables_hv->SetBinContent(bin, mean_adc_minus_pedestal_cables_hv);

            int shihai_label = shihai_channel_map[local_ch];

            // Only label every 10 bins
            if (bin % 10 == 1) {
                h2_no_cables->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_label));
                h2_cables->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_label));
                h2_cables_hv->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_label));
            }
            else {
                h2_no_cables->GetXaxis()->SetBinLabel(bin, "");
                h2_cables->GetXaxis()->SetBinLabel(bin, "");
                h2_cables_hv->GetXaxis()->SetBinLabel(bin, "");
            }
        }

        h2_no_cables->SetLineColor(kRed + 1);
        h2_cables->SetLineColor(kBlue + 1);
        h2_cables_hv->SetLineColor(kGreen + 1);
        style_histogram(h2_no_cables, Form("ADC - Calibration Pedestal NEW (Protoboard %s)", protoboards[fpga].c_str()), "ADC - Calibration Pedestal");
        //h2_no_cables->SetMaximum(1024);
        h2_no_cables->SetMaximum(100); // Set maximum for better visibility

        h2_no_cables->Draw("HIST");
        h2_cables->Draw("HIST SAME");
        h2_cables_hv->Draw("HIST SAME");

        TLegend *leg2 = new TLegend(0.15, 0.75, 0.28, 0.88);
        leg2->AddEntry(h2_no_cables, Form("No Cables (run %d)", run_no_cables), "l");
        leg2->AddEntry(h2_cables, Form("Cables (run %d)", run_cables), "l");
        leg2->AddEntry(h2_cables_hv, Form("Cables HV (run %d)", run_cables_hv), "l");
        leg2->Draw();

        c2->SaveAs(out2);

        // Plot channel_histrograms in a separate canvas
        TCanvas *c_hist = new TCanvas(Form("c_hist_fpga%d", fpga), Form("Channel Histograms - Protoboard %s", protoboards[fpga].c_str()), 2000, 800);
        c_hist->Divide(3, 5); // 3 columns, 5 rows
        c_hist->SetGrid();

        for (int ch = 0; ch < 1; ++ch) {
            //c_hist->cd(ch + 1);
            channel_histograms_no_cables[ch]->Draw("HIST");
        }
        c_hist->SaveAs(output_dir + Form("channel_histograms_fpga%d.png", fpga));

        // Typical Channel.
        // FPGA 0. 8, 64 // FPGA 1. 42, 126 // FPGA 2. 67, 140

        // Plot the different configurations in a single canvas for channel 39 in Shihai map for fpga 2
        // if (fpga == 1) // Channel 39 in Shihai map corresponds to local channel 36 in FPGA 1
        // {
        //     TCanvas *c_channel_39 = new TCanvas("c_channel_39", "Channel 39 Histograms", 800, 600);
        //     c_channel_39->SetGrid();

        //     // Set the title and axis labels
        //     c_channel_39->SetTitle("Pedestal distribution for a stable channel (Channel 39 in Protoboard 008)");

        //     int actual_channel = 36 + 144*fpga; // Channel 39 in Shihai map corresponds to local channel 97 in FPGA 1

        //     // Clear the default histogram title
        //     channel_histograms_no_cables[actual_channel]->SetTitle("");

        //     // Set vertical limits for better visibility
        //     //channel_histograms_no_cables[actual_channel]->GetYaxis()->SetRangeUser(0, 50); // Adjust as needed

        //     channel_histograms_no_cables[actual_channel]->GetXaxis()->SetTitle("ADC Value");
        //     channel_histograms_no_cables[actual_channel]->GetYaxis()->SetTitle("Entries");
        //     channel_histograms_no_cables[actual_channel]->SetStats(0); // Disable stats box
        //     channel_histograms_no_cables[actual_channel]->SetLineColor(kRed);
        //     channel_histograms_no_cables[actual_channel]->SetFillColor(kRed - 9);
        //     channel_histograms_no_cables[actual_channel]->SetFillStyle(3004);
        //     channel_histograms_no_cables[actual_channel]->SetLineWidth(2);

        //     channel_histograms_cables[actual_channel]->SetLineColor(kBlue - 9);
        //     channel_histograms_cables[actual_channel]->SetFillColor(kBlue);
        //     channel_histograms_cables[actual_channel]->SetFillStyle(3004);
        //     channel_histograms_cables[actual_channel]->SetLineWidth(2);

        //     channel_histograms_cables_hv[actual_channel]->SetLineColor(kGreen + 2);
        //     channel_histograms_cables_hv[actual_channel]->SetFillColor(kGreen + 1);
        //     channel_histograms_cables_hv[actual_channel]->SetFillStyle(3004);
        //     channel_histograms_cables_hv[actual_channel]->SetLineWidth(2);

        //     channel_histograms_no_cables[actual_channel]->Draw("HIST");
        //     channel_histograms_cables[actual_channel]->Draw("HIST SAME");
        //     channel_histograms_cables_hv[actual_channel]->Draw("HIST SAME");

        //     // Add a text box with the total number of channels
        //     TLatex *text_channel_1 = new TLatex(0.60, 0.70, "Total Channels: 182");
        //     text_channel_1->SetNDC();
        //     text_channel_1->SetTextSize(0.035);
        //     text_channel_1->Draw();

        //     // RMS values for each configuration // Round to 2 decimal places
        //     double rms_run_no_cables = channel_histograms_no_cables[actual_channel]->GetRMS();
        //     double rms_run_cables = channel_histograms_cables[actual_channel]->GetRMS();
        //     double rms_run_cables_hv = channel_histograms_cables_hv[actual_channel]->GetRMS();

        //     TLegend *leg_channel_1 = new TLegend(0.15, 0.75, 0.45, 0.88);
        //     leg_channel_1->AddEntry(channel_histograms_no_cables[actual_channel], Form("No Cables - RMS: %.2f (run %d)", rms_run_no_cables, run_no_cables), "f");
        //     leg_channel_1->AddEntry(channel_histograms_cables[actual_channel], Form("Cables - RMS: %.2f (run %d)", rms_run_cables, run_cables), "f");
        //     leg_channel_1->AddEntry(channel_histograms_cables_hv[actual_channel], Form("Cables HV - RMS: %.2f (run %d)", rms_run_cables_hv, run_cables_hv), "f");
        //     leg_channel_1->Draw();

        //     // Create custom title using TLatex
        //     TLatex *text_title_1 = new TLatex(0.13, 0.97, "Pedestal distribution for a channel above the mean RMS values");
        //     text_title_1->SetNDC();
        //     text_title_1->SetTextSize(0.035);
        //     text_title_1->Draw();

        //     TLatex *text_title_2 = new TLatex(0.20, 0.92, "Ch 39 in Pb 01 w/o calibration pedestal subtraction");
        //     text_title_2->SetNDC();
        //     text_title_2->SetTextSize(0.035);
        //     text_title_2->Draw();

        //     c_channel_39->SaveAs(output_dir + "channel_39_histograms.png");

            // TCanvas *c_channel_39_samples = new TCanvas("c_channel_39_samples", "Channel 39 Samples", 1800, 600);
            // c_channel_39_samples->Divide(3, 1); // 3 columns, 1 row
            // c_channel_39_samples->SetGrid();

            // adc_samples_histograms_no_cables[actual_channel]->SetTitle("No Cables (INA36 - ASIC 0 - Pb 06)");
            // adc_samples_histograms_cables[actual_channel]->SetTitle("Cables (INA36 - ASIC 0 - Pb 06)");
            // adc_samples_histograms_cables_hv[actual_channel]->SetTitle("Cables + voltage (INA36 - ASIC 0 - Pb 06)");

            // c_channel_39_samples->cd(1);
            // adc_samples_histograms_no_cables[actual_channel]->GetXaxis()->SetTitle("Sample instant");
            // adc_samples_histograms_no_cables[actual_channel]->GetYaxis()->SetTitle("ADC");
            // adc_samples_histograms_no_cables[actual_channel]->SetStats(0); // Disable stats box
            // adc_samples_histograms_no_cables[actual_channel]->GetYaxis()->SetRangeUser(0, 200); // Set Y-axis range for better visibility
            // // Get the maximum value of the histogram to set Z-axis range and store in a integer
            // int max_value = adc_samples_histograms_no_cables[actual_channel]->GetMaximum();
            // adc_samples_histograms_no_cables[actual_channel]->GetZaxis()->SetRangeUser(0, max_value); // Set Z-axis range for better visibility

            // adc_samples_histograms_no_cables[actual_channel]->Draw("colz");

            // c_channel_39_samples->cd(2);

            // adc_samples_histograms_cables[actual_channel]->GetXaxis()->SetTitle("Sample instant");
            // adc_samples_histograms_cables[actual_channel]->GetYaxis()->SetTitle("ADC");
            // adc_samples_histograms_cables[actual_channel]->SetStats(0); // Disable stats box
            // adc_samples_histograms_cables[actual_channel]->GetYaxis()->SetRangeUser(0, 200); // Set Y-axis range for better visibility
            // adc_samples_histograms_cables[actual_channel]->GetZaxis()->SetRangeUser(0, max_value); // Set Z-axis range for better visibility
            // adc_samples_histograms_cables[actual_channel]->Draw("colz");

            // c_channel_39_samples->cd(3);
            // adc_samples_histograms_cables_hv[actual_channel]->GetXaxis()->SetTitle("Sample instant");
            // adc_samples_histograms_cables_hv[actual_channel]->GetYaxis()->SetTitle("ADC");
            // adc_samples_histograms_cables_hv[actual_channel]->SetStats(0); // Disable stats box
            // adc_samples_histograms_cables_hv[actual_channel]->GetYaxis()->SetRangeUser(0, 200); // Set Y-axis range for better visibility
            // adc_samples_histograms_cables_hv[actual_channel]->GetZaxis()->SetRangeUser(0, max_value); // Set Z-axis range for better visibility
            // adc_samples_histograms_cables_hv[actual_channel]->Draw("colz");

            // TPad* titlePad = new TPad("titlePad", "titlePad", 0, 0.91, 1, 1.0); // Top 9% of canvas
            // titlePad->SetFillStyle(0);
            // titlePad->SetFrameFillStyle(0);
            // titlePad->Draw();
            // titlePad->cd();

            // TLatex* title = new TLatex(0.5, 0.5, "Pedestal Distribution - Channel 39 in Protoboard 06");
            // title->SetTextAlign(22);
            // title->SetTextSize(0.04);
            // title->Draw();

            // c_channel_39_samples->SaveAs(output_dir + "channel_39_samples_histograms.png");

        // }

        std::map<int, std::vector<int>> typical_channels = {
            {0, {7, 60}}, // 8, 64
            {1, {39, 119}}, // 42, 126
            {2, {63}} // 67
        };

        std::map<int, std::vector<int>> non_stable_channels = {
            {0, {1, 72}}, // 2, 77
            //{1, {4, 96}}, // 5, 102
            {2, {4, 21, 87}} // 5, 23, 92
        };

        std::map<int, std::vector<int>> ideal_channels = {
            {0, {13, 15, 16, 40}}, // 14, 16, 17, 43
            {2, {97}} // 103
        };

        std::map<int, std::vector<int>> good_channels = {
            {0, {7, 41}}, // 8, 44
            {1, {36, 108, 114}} // 39, 115, 121
        };

        //std::string type = "Typical Channel"; // Type of channels to plot

        std::string type = "Non Stable Channel"; // Type of channels to plot

        //std::string type = "Ideal Channel"; // Type of channels to plot

        //std::string type = "Good Channel"; // Type of channels to plot

        std::map<int, std::string> fpga_to_pb = {
            {0, "01"},
            {1, "008"},
            {2, "06"}
        };

        // Create one big canvas: 5 rows (channels), 3 columns (conditions)
        TCanvas* c_all_channels = new TCanvas("c_all_channels", "All Channels Samples", 1800, 3600);
        c_all_channels->Divide(3, 5); // 3 columns (tests), 5 rows (channels)
        c_all_channels->SetGrid();

        // Create a canvas for all channels waveforms
        std::vector<int> selected_events;
        for (int i = 0; i <= 150; ++i) selected_events.push_back(i);

        TCanvas* c_all_events_waveforms = new TCanvas("c_all_channels_waveforms", "All Channels Waveforms", 2400, 2880);
        c_all_events_waveforms->Divide(15, 10); // 15 events
        c_all_events_waveforms->SetGrid();

        int pad_index = 1;
        int pad_index_waveforms = 1;

        for (const auto& [fpga, channels] : non_stable_channels) {
            for (int channel : channels) {
                int actual_channel = channel + 144 * fpga;

                // Determine ASIC
                TString asic = (channel <= 71) ? "ASIC 0" : "ASIC 1";

                // Build base title
                TString title_base = Form("INA%d - %s - Pb %s", channel % 72,
                                        asic.Data(),
                                        fpga_to_pb[fpga].c_str());

                int max_value_no = adc_samples_histograms_no_cables[actual_channel]->GetMaximum();
                int max_value_cables = adc_samples_histograms_cables[actual_channel]->GetMaximum();
                int max_value_cables_hv = adc_samples_histograms_cables_hv[actual_channel]->GetMaximum();

                int max_value = std::max({max_value_no, max_value_cables, max_value_cables_hv});

                // No Cables
                c_all_channels->cd(pad_index++);
                adc_samples_histograms_no_cables[actual_channel]->SetTitle(type + " (" + title_base + ") - No Cables");
                adc_samples_histograms_no_cables[actual_channel]->GetXaxis()->SetTitle("Sample instant");
                adc_samples_histograms_no_cables[actual_channel]->GetYaxis()->SetTitle("ADC");
                adc_samples_histograms_no_cables[actual_channel]->SetStats(0);
                adc_samples_histograms_no_cables[actual_channel]->GetYaxis()->SetRangeUser(0, 200);
                //adc_samples_histograms_no_cables[actual_channel]->GetZaxis()->SetRangeUser(0, max_value);
                adc_samples_histograms_no_cables[actual_channel]->Draw("colz");

                double rms_run_no_cables = channel_histograms_no_cables[actual_channel]->GetRMS();

                TLatex* label_no_cables = new TLatex(0.15, 0.85, Form("RMS: %.2f (Mean RMS = %.2f)", rms_run_no_cables, mean_rms_no_cables));
                label_no_cables->SetNDC();
                label_no_cables->SetTextSize(0.04);
                label_no_cables->Draw();

                // Cables
                c_all_channels->cd(pad_index++);
                adc_samples_histograms_cables[actual_channel]->SetTitle(type + " (" + title_base + ") - Cables");
                adc_samples_histograms_cables[actual_channel]->GetXaxis()->SetTitle("Sample instant");
                adc_samples_histograms_cables[actual_channel]->GetYaxis()->SetTitle("ADC");
                adc_samples_histograms_cables[actual_channel]->SetStats(0);
                adc_samples_histograms_cables[actual_channel]->GetYaxis()->SetRangeUser(0, 200);
                //adc_samples_histograms_cables[actual_channel]->GetZaxis()->SetRangeUser(0, max_value);
                adc_samples_histograms_cables[actual_channel]->Draw("colz");

                double rms_run_cables = channel_histograms_cables[actual_channel]->GetRMS();
                TLatex* label_cables = new TLatex(0.15, 0.85, Form("RMS: %.2f (Mean RMS = %.2f)", rms_run_cables, mean_rms_cables));
                label_cables->SetNDC();
                label_cables->SetTextSize(0.04);
                label_cables->Draw();

                // Cables + HV
                c_all_channels->cd(pad_index++);
                adc_samples_histograms_cables_hv[actual_channel]->SetTitle(type + " (" + title_base + ") - Cables + V");
                adc_samples_histograms_cables_hv[actual_channel]->GetXaxis()->SetTitle("Sample instant");
                adc_samples_histograms_cables_hv[actual_channel]->GetYaxis()->SetTitle("ADC");
                adc_samples_histograms_cables_hv[actual_channel]->SetStats(0);
                adc_samples_histograms_cables_hv[actual_channel]->GetYaxis()->SetRangeUser(0, 200);
                //adc_samples_histograms_cables_hv[actual_channel]->GetZaxis()->SetRangeUser(0, max_value);
                adc_samples_histograms_cables_hv[actual_channel]->Draw("colz");
                
                double rms_run_cables_hv = channel_histograms_cables_hv[actual_channel]->GetRMS();
                TLatex* label_cables_hv = new TLatex(0.15, 0.85, Form("RMS: %.2f (Mean RMS = %.2f)", rms_run_cables_hv, mean_rms_cables_hv));
                label_cables_hv->SetNDC();
                label_cables_hv->SetTextSize(0.04);
                label_cables_hv->Draw();

                for (int evt : selected_events) {
                    c_all_events_waveforms->cd(pad_index_waveforms++);
                    TGraph* g_no_cables = new TGraph(num_samples);
                    TGraph* g_cables = new TGraph(num_samples);
                    TGraph* g_cables_hv = new TGraph(num_samples);

                    for (int s = 0; s < num_samples; ++s) {
                        g_no_cables->SetPoint(s, s, adc_event_no_cables[evt][actual_channel][s]);
                        g_cables->SetPoint(s, s, adc_event_cables[evt][actual_channel][s]);
                        g_cables_hv->SetPoint(s, s, adc_event_cables_hv[evt][actual_channel][s]);
                    }

                    g_no_cables->SetTitle(Form("%s - Event %d;Sample instant;ADC", title_base.Data(), evt));
                    //g->SetLineColor(kBlue + (evt / 10) % 9); // cycle colors
                    g_no_cables->SetLineWidth(1.5);
                    g_no_cables->SetLineColor(kRed);
                    //g_no_cables->SetMarkerStyle(20); // cycle marker styles
                    g_no_cables->SetMarkerColor(kRed);
                    //g_no_cables->SetMarkerSize(1);
                    
                    g_no_cables->GetXaxis()->SetLimits(0, num_samples);
                    g_no_cables->GetYaxis()->SetRangeUser(0, 200); // Set Y-axis range for better visibility

                    g_cables->SetLineWidth(1.5);
                    g_cables->SetLineColor(kBlue);
                    //g_cables->SetMarkerStyle(21); // cycle marker styles
                    g_cables->SetMarkerColor(kBlue);
                    //g_cables->SetMarkerSize(1);
                    // g_cables->SetFillColor(kBlue - 4);
                    // g_cables->SetFillStyle(3004);

                    g_cables_hv->SetLineWidth(1.5);
                    g_cables_hv->SetLineColor(kGreen + 2);
                    //g_cables_hv->SetMarkerStyle(22); // cycle marker styles
                    g_cables_hv->SetMarkerColor(kGreen + 2);
                    //g_cables_hv->SetMarkerSize(1);
                    // g_cables_hv->SetFillColor(kGreen + 2 - 4);
                    // g_cables_hv->SetFillStyle(3004);

                    g_no_cables->Draw("ALP");
                    g_cables->Draw("LP SAME");
                    g_cables_hv->Draw("LP SAME");

                    // Add legend
                    TLegend* legend_waveform = new TLegend(0.65, 0.75, 0.88, 0.88);
                    legend_waveform->AddEntry(g_no_cables, "No Cables", "l");
                    legend_waveform->AddEntry(g_cables, "Cables", "l");
                    legend_waveform->AddEntry(g_cables_hv, "Cables + V", "l");
                    legend_waveform->SetTextSize(0.035);
                    legend_waveform->SetBorderSize(0);
                    legend_waveform->Draw();

                }
                // Save the waveform canvas
                c_all_events_waveforms->SaveAs(output_dir + Form("waveforms_%s.png", title_base.Data()));

                pad_index_waveforms = 1; // Reset pad index for next channel
                c_all_events_waveforms->Clear();
                c_all_events_waveforms->Divide(15, 10); // Reset to 3 columns, 5 rows
            }
        }

        // Save final canvas
        c_all_channels->SaveAs(output_dir + type + "_all_channels_samples_histograms.png");

        TCanvas* c_all_nonstable = new TCanvas("c_all_nonstable", "Non-Stable Channel Pedestals", 1800, 900);
        c_all_nonstable->Divide(3, 2); // Adjust rows/cols to match # of channels
        int pad_index_distr = 1;

        for (const auto& [fpga, channels] : good_channels) {
            for (int channel : channels) {
                int actual_channel = channel + 144 * fpga;

                // Determine ASIC
                TString asic = (channel <= 71) ? "ASIC 0" : "ASIC 1";

                // Draw to current pad
                c_all_nonstable->cd(pad_index_distr++);

                // Format histograms
                channel_histograms_no_cables[actual_channel]->SetLineColor(kRed);
                channel_histograms_no_cables[actual_channel]->SetFillColor(kRed - 9);
                channel_histograms_no_cables[actual_channel]->SetFillStyle(3004);
                channel_histograms_no_cables[actual_channel]->SetLineWidth(2);
                channel_histograms_no_cables[actual_channel]->SetStats(0);
                channel_histograms_no_cables[actual_channel]->GetXaxis()->SetTitle("ADC Value");
                channel_histograms_no_cables[actual_channel]->GetYaxis()->SetTitle("Entries");
                channel_histograms_no_cables[actual_channel]->SetTitle(Form("Pedestal Distribution. INA%d - %s - Pb %s", channel % 72, asic.Data(), fpga_to_pb[fpga].c_str()));
                gPad->SetLogy();

                channel_histograms_cables[actual_channel]->SetLineColor(kBlue - 9);
                channel_histograms_cables[actual_channel]->SetFillColor(kBlue);
                channel_histograms_cables[actual_channel]->SetFillStyle(3004);
                channel_histograms_cables[actual_channel]->SetLineWidth(2);
                channel_histograms_cables[actual_channel]->SetTitle("");

                channel_histograms_cables_hv[actual_channel]->SetLineColor(kGreen + 2);
                channel_histograms_cables_hv[actual_channel]->SetFillColor(kGreen + 1);
                channel_histograms_cables_hv[actual_channel]->SetFillStyle(3004);
                channel_histograms_cables_hv[actual_channel]->SetLineWidth(2);
                channel_histograms_cables_hv[actual_channel]->SetTitle("");

                // Draw histograms overlaid
                channel_histograms_no_cables[actual_channel]->Draw("HIST");
                channel_histograms_cables[actual_channel]->Draw("HIST SAME");
                channel_histograms_cables_hv[actual_channel]->Draw("HIST SAME");

                // Compute RMS
                double rms_run_no_cables = channel_histograms_no_cables[actual_channel]->GetRMS();
                double rms_run_cables = channel_histograms_cables[actual_channel]->GetRMS();
                double rms_run_cables_hv = channel_histograms_cables_hv[actual_channel]->GetRMS();

                // Draw legend
                TLegend *leg = new TLegend(0.55, 0.70, 0.90, 0.88);
                leg->SetTextSize(0.03);
                leg->AddEntry(channel_histograms_no_cables[actual_channel], Form("No Cables - RMS: %.2f", rms_run_no_cables), "f");
                leg->AddEntry(channel_histograms_cables[actual_channel], Form("Cables - RMS: %.2f", rms_run_cables), "f");
                leg->AddEntry(channel_histograms_cables_hv[actual_channel], Form("Cables HV - RMS: %.2f", rms_run_cables_hv), "f");
                leg->Draw();

                // // Dynamic title text
                // TString title_1 = Form("Pedestal Distribution. INA%d - %s - Pb %s", channel % 72, asic.Data(), fpga_to_pb[fpga].c_str());
                // TLatex *text1 = new TLatex(0.30, 0.95, title_1);
                // text1->SetNDC();
                // text1->SetTextSize(0.035);
                // text1->Draw();
            }
        }

        // Save the full canvas
        c_all_nonstable->SaveAs(output_dir + "good_channels_pedestal_overlay.png");


    }

    f_no_cables->Close();
    f_cables->Close();
    f_cables_hv->Close();
    
    std::cout << "Total Numbers of Entries: " << nEntries << std::endl;
}
