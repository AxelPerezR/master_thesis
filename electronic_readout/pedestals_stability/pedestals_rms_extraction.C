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

// Function to write RMS values to a CSV file
void WriteRMSCSV(const std::string& filename,
                 const std::vector<int>& channel,
                 const std::vector<double>& rms_no,
                 const std::vector<double>& rms_cab,
                 const std::vector<double>& rms_cab_hv){
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    file << "index_pb,channel_shihai,PB_INA,rms_no,rms_cab,rms_cab_hv\n";
    size_t n = channel.size();
    int index_asic = 0;
    for (size_t i = 0; i < n; ++i) {
        index_asic = i % 72; // Assuming 144 channels per ASIC
        file << i << "," << channel[i] << "," << index_asic << "," << rms_no[i] << "," << rms_cab[i] << ","
             << std::fixed << std::setprecision(3) << rms_cab_hv[i] << "\n";
    }

    file.close();
    std::cout << " \n CSV file written: " << filename << std::endl;
}

// Function to analyze RMS values of pedestals
// This function reads the ADC values from ROOT files, calculates the RMS for each channel, and
// writes the results to a CSV file.
// The function takes three run numbers as input parameters, which correspond to different data sets.
// The function assumes that the input files are structured in a specific way, with ADC values stored in branches named "adc" and "hit_pedestal".
// The first run corresponds to data without cables, the second run with cables, and the third run with high voltage cables.
// Run in terminal: root 'pedestals_rms_analysis.C(run_no_cables, run_cables, run_cables_hv)', where run_no_cables, run_cables, and run_cables_hv are the run numbers for the respective data sets.
void pedestals_rms_extraction(int run_no_cables = 32, int run_cables = 33, int run_cables_hv = 34)
{
    const int num_kcu = 3; // Modify this if the number of KCUs changes
    const int channels_per_kcu = 144; // Number of readout channels per protoboard // Do not modify
    const int num_samples = 20; // Number of samples per channel
    const int num_channels = num_kcu * channels_per_kcu;

    // Shihai mapping
    std::vector<int> shihai_channel_map;
    for (int ch = 1; ch <= 151; ++ch)
    {   
        // Skip 'dead channels' in the protoboard
        // These channels are marked as 'dead' during the pedestal alignment process
        if (ch == 19 || ch == 38 || ch == 57 || ch == 76 || ch == 95 || ch == 114 || ch == 133)
            continue;
        shihai_channel_map.push_back(ch);
    }

    // FPGA IDs and protoboard names
    std::vector<int> fpga_ids = {208, 209, 210}; // IP addresses of the KCUs // Modify as needed
    std::vector<std::string> protoboards = {"01", "008", "06"}; // Protoboard names corresponding to the FPGA IDs // Modify as needed

    const TString input_dir = "/home/axelperezr/eic/pedestals/data/"; // Directory where the input ROOT files are located // Modify as needed
    const TString output_dir = "/home/axelperezr/eic/pedestals/plots/"; // Directory where the output plots and csv will be saved // Modify as needed

    // Create file names based on run numbers
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

    // Get the TTree from each file
    TTree *t_no_cables = (TTree *)f_no_cables->Get("events");
    TTree *t_cables = (TTree *)f_cables->Get("events");
    TTree *t_cables_hv = (TTree *)f_cables_hv->Get("events");
    if (!t_no_cables || !t_cables || !t_cables_hv)
    {
        std::cerr << "Cannot find TTree 'events' in files." << std::endl;
        return;
    }

    // Arrays to hold ADC values for each channel
    // Each channel has 20 samples, so we use a 2D array for each run

    uint adc_no_cables[num_channels][num_samples];
    uint adc_cables[num_channels][num_samples];
    uint adc_cables_hv[num_channels][num_samples];

    t_no_cables->SetBranchAddress("adc", adc_no_cables);
    t_cables->SetBranchAddress("adc", adc_cables);
    t_cables_hv->SetBranchAddress("adc", adc_cables_hv);

    Long64_t nEntries = std::min({t_no_cables->GetEntries(), t_cables->GetEntries(), t_cables_hv->GetEntries()});

    std::vector<TH1I*> channel_histograms_no_cables(num_channels);
    std::vector<TH1I*> channel_histograms_cables(num_channels);
    std::vector<TH1I*> channel_histograms_cables_hv(num_channels);

    for (int ch = 0; ch < num_channels; ++ch)
    {
        channel_histograms_no_cables[ch] = new TH1I(Form("h_no_cables_ch_%d", ch), Form("Channel %d - No Cables", ch), 200, 0, 200);
        channel_histograms_cables[ch] = new TH1I(Form("h_cables_ch_%d", ch), Form("Channel %d - Cables", ch), 200, 0, 200);
        channel_histograms_cables_hv[ch] = new TH1I(Form("h_cables_hv_ch_%d", ch), Form("Channel %d - Cables HV", ch), 200, 0, 200);
    }

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        t_no_cables->GetEntry(i);
        t_cables->GetEntry(i);
        t_cables_hv->GetEntry(i);

        for (int ch = 0; ch < num_channels; ++ch)
        {
            for (int sample = 0; sample < num_samples; ++sample)
            {
                int content_no_cables = adc_no_cables[ch][sample];
                int content_cables = adc_cables[ch][sample];
                int content_cables_hv = adc_cables_hv[ch][sample];

                channel_histograms_no_cables[ch]->Fill(content_no_cables);
                channel_histograms_cables[ch]->Fill(content_cables);
                channel_histograms_cables_hv[ch]->Fill(content_cables_hv);
            }
        }
    }

    for (int fpga = 0; fpga < num_kcu; ++fpga)
    {
        int ch_start = fpga * channels_per_kcu; // Starting channel for the current FPGA

        std::vector<int> shihai_values; // shihai channel values
        std::vector<double> rms_no_cables;
        std::vector<double> rms_cables;
        std::vector<double> rms_cables_hv;

        for (int local_ch = 0; local_ch < channels_per_kcu; ++local_ch)
        {

            int global_ch = ch_start + local_ch;
            int bin = local_ch + 1;
            int shihai_label = shihai_channel_map[local_ch];

            double rms_no = channel_histograms_no_cables[global_ch]->GetRMS();
            double rms_cab = channel_histograms_cables[global_ch]->GetRMS();
            double rms_hv = channel_histograms_cables_hv[global_ch]->GetRMS();
            
            //index_values.push_back(global_ch);
            shihai_values.push_back(shihai_label);
            rms_no_cables.push_back(rms_no);
            rms_cables.push_back(rms_cab);
            rms_cables_hv.push_back(rms_hv);
        }

        // Write RMS values to CSV
        std::string csv_filename = "RMS_values_fpga" + std::to_string(fpga) + ".csv";
        WriteRMSCSV(csv_filename, shihai_values, rms_no_cables, rms_cables, rms_cables_hv);
    }

    f_no_cables->Close();
    f_cables->Close();
    f_cables_hv->Close();
}
