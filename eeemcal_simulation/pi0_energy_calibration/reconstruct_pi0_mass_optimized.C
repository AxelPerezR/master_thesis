#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <vector>

// ROOT modern physics vector
#include "Math/Vector4D.h"     // For ROOT::Math::PxPyPzEVector
#include "Math/Vector3D.h"     // For ROOT::Math::XYZVector

#include "edm4eic/ClusterData.h"  // Cluster struct definition

#include <iomanip> // for std::setprecision

// Use TTreeReader
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

void WritePhotonEventCSV(const std::string& filename,
                         const std::vector<int>& event,
                         const std::vector<int>& photon,
                         const std::vector<int>& crystal_id,
                         const std::vector<float>& x,
                         const std::vector<float>& y,
                         const std::vector<float>& photon_energy) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    file << "event,photon,crystal_id,x,y,photon_energy\n";
    size_t n = event.size();
    for (size_t i = 0; i < n; ++i) {
        file << event[i] << "," << photon[i] << "," << crystal_id[i] << ","
             << std::fixed << std::setprecision(3) << x[i] << "," << y[i] << ","
             << std::fixed << std::setprecision(3) << photon_energy[i] << "\n";
    }

    file.close();
    std::cout << " \n CSV file written: " << filename << std::endl;
}

std::vector<std::tuple<int, float, float>> loadCrystalPositions(const std::string& filename) {
    std::vector<std::tuple<int, float, float>> crystals;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open crystal position file." << std::endl;
        return crystals;
    }

    std::string line;
    std::getline(file, line);  // Skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string id_str, x_str, y_str;
        std::getline(ss, id_str, ',');
        std::getline(ss, x_str, ',');
        std::getline(ss, y_str, ',');
        crystals.emplace_back(std::stoi(id_str), std::stof(x_str), std::stof(y_str));
    }

    return crystals;
}

int findClosestCrystalID(float cluster_x, float cluster_y,
                         const std::vector<std::tuple<int, float, float>>& crystals) {
    float min_dist = 1e6;
    int best_id = -1;
    for (const auto& [id, x, y] : crystals) {
        float dist = std::hypot(x - cluster_x, y - cluster_y);
        if (dist < min_dist) {
            min_dist = dist;
            best_id = id;
        }
    }
    return best_id;
}

// Loads the correction factors from the file
std::unordered_map<int, float> loadCorrections(const std::string& correction_file) {
    std::unordered_map<int, float> corrections;
    std::ifstream file(correction_file);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open correction file " << correction_file << std::endl;
        return corrections; // returns empty map
    }

    std::string line;
    // Skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string id_str, correction_str;

        std::getline(ss, id_str, ',');
        std::getline(ss, correction_str, ',');

        try {
            int id = std::stoi(id_str);
            float correction = std::stof(correction_str);
            corrections[id] = correction;
        } catch (const std::exception& e) {
            std::cerr << "Warning: Skipping malformed line: " << line << std::endl;
        }
    }

    return corrections;
}

// Applies the correction using the pre-loaded map
float applyCorrection(float photon_energy, int crystal_id, const std::unordered_map<int, float>& corrections) {
    auto it = corrections.find(crystal_id);
    if (it != corrections.end()) {
        return photon_energy * it->second;
    } else {
        std::cerr << "Warning: Crystal ID " << crystal_id << " not found in correction data." << std::endl;
        return photon_energy;
    }
}

void reconstruct_pi0_mass_optimized() {
    TStopwatch timer;
    timer.Start();

    TChain chain("events");

    std::string collisionEnergy = "18x275"; // Collision energy, can be "5x41" or "18x275"
    int energyLimit = 20; // GeV
    const int nFiles = 1394;  // Maximum number of files // 564 is the maximum number of files in the dataset DIS 18x275 // 306 at 5x41 GeV for version 25.03.1, 311 for version 25.04.1, 380 for 25.05.0, 1501 for 25.02.0, 1472 for 25.01.1, 406 for 24.10.0
    // for SIDIS datasets, the maximum number of files is 1385 in 25.03.1, 1394 in 25.04.1
    int numberOfDatasets = 22; // Number of datasets to pick from // DIS datasets: 1-5, SIDIS datasets: 1-9
    std::string version = "25.04.1"; // Version of the dataet, can be "25.03.1" or "25.04.1"
    bool DIS = false; // Set to true if you want to use DIS datasets
    bool SIDIS = true; // Set to true if you want to use SIDIS datasets
    // there are 50 runs available for SIDIS datasets, but only 9 are used in the analysis

    std::cout << "Process: " << (DIS ? "DIS" : "SIDIS") << std::endl;  
    std::cout << "Number of datasets: " << numberOfDatasets << std::endl;
    std::cout << "Number of files per dataset: " << nFiles << std::endl;
    std::cout << "Using collision energy: " << collisionEnergy << std::endl;


    if (DIS == false && SIDIS == false) {
        std::cerr << "Please set either DIS or SIDIS to true." << std::endl;
        return;
    } else if (DIS == true && SIDIS == true) {
        std::cerr << "Please set either DIS or SIDIS to true, not both." << std::endl;
        return;
    }

    if (DIS == true) {
        for (int i = 1; i <= numberOfDatasets; ++i) {
            for (int j = 1; j <= nFiles; ++j) {
                if (i==3 && j == 21 && collisionEnergy == "5x41") continue;
                if (i==1 && j == 214 && collisionEnergy == "18x275") continue;
                if (i==3 && j == 9 && collisionEnergy == "18x275") continue;
                if (i==3 && j == 416 && collisionEnergy == "18x275") continue;

                TString filename = Form(
                    "root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/25.03.1/"
                    "epic_craterlake/DIS/NC/%s/minQ2=1/"
                    "pythia8NCDIS_%s_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_%d.%04d.eicrecon.edm4eic.root",
                    collisionEnergy.c_str(), collisionEnergy.c_str(), i, j);

                std::cout << "Adding file: " << filename << std::endl;
                chain.Add(filename);
            }
        }
    } else if (SIDIS == true) {
        for (int i = numberOfDatasets-1; i <= numberOfDatasets; ++i) {
            for (int j = 1; j <= nFiles; ++j) {
                if (i==5 && j == 176 && collisionEnergy == "5x41" && version == "25.05.0") continue;
                if (i==8 && j == 127 && collisionEnergy == "5x41" && version == "25.05.0") continue;
                if (i==9 && j == 308 && collisionEnergy == "5x41" && version == "25.05.0") continue;
                if (i==8 && j == 15 && collisionEnergy == "5x41" && version == "25.04.1") continue;
                if (i==1 && j == 146 && collisionEnergy == "5x41" && version == "25.03.1") continue;

                if (i==1 && j == 536 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==1 && j == 583 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==2 && j == 786 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==2 && j == 871 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==3 && j == 678 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==4 && j == 431 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==4 && j == 602 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==4 && j == 869 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==5 && j == 1117 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==6 && j == 627 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==8 && j == 837 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==10 && j == 471 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==13 && j == 1149 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==16 && j == 198 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==17 && j == 1310 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==24 && j == 1191 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==24 && j == 518 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==25 && j == 84 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==26 && j == 38 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==26 && j == 673 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==26 && j == 813 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==26 && j == 1278 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==27 && j == 626 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==32 && j == 971 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==36 && j == 559 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==36 && j == 939 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==33 && j == 934 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==38 && j == 295 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==38 && j == 714 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==38 && j == 756 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==38 && j == 814 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==39 && j == 824 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==40 && j == 970 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==41 && j == 148 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==41 && j == 253 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==41 && j == 256 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==41 && j == 303 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==41 && j == 417 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==41 && j == 1130 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==42 && j == 361 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==43 && j == 677 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==47 && j == 239 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==47 && j == 970 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==47 && j == 1288 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==48 && j == 56 && collisionEnergy == "18x275" && version == "25.03.1") continue;
                if (i==49 && j == 1193 && collisionEnergy == "18x275" && version == "25.03.1") continue;

                if (i==19 && j == 2 && collisionEnergy == "5x41" && version == "25.03.1") continue;
                if (i==26 && j == 128 && collisionEnergy == "5x41" && version == "25.03.1") continue;
                if (i==29 && j == 121 && collisionEnergy == "5x41" && version == "25.03.1") continue;
                if (i==31 && j == 188 && collisionEnergy == "5x41" && version == "25.03.1") continue;
                if (i==35 && j == 119 && collisionEnergy == "5x41" && version == "25.03.1") continue;
                if (i==38 && j == 51 && collisionEnergy == "5x41" && version == "25.03.1") continue;
                if (i==41 && j == 131 && collisionEnergy == "5x41" && version == "25.03.1") continue;
                if (i==46 && j == 69 && collisionEnergy == "5x41" && version == "25.03.1") continue;

                TString filename = Form(
                    "root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/%s/"
                    "epic_craterlake/SIDIS/pythia6-eic/1.0.0/%s/q2_0to1/"
                    "pythia_ep_noradcor_%s_q2_0.000000001_1.0_run%d.ab.%04d.eicrecon.edm4eic.root",
                    version.c_str(), collisionEnergy.c_str(), collisionEnergy.c_str(), i, j);

                std::cout << "Adding file: " << filename << std::endl;
                chain.Add(filename);
            }
        }
    } else
    {
        std::cerr << "Invalid dataset type. Please set either DIS or SIDIS to true." << std::endl;
        return;
    }

    // Branch pointer replaced by TTreeReader
    TTreeReader reader(&chain);
    TTreeReaderValue<std::vector<edm4eic::ClusterData>> clusters(reader, "EcalEndcapNClusters");

    // Histogram for pi0 mass
    double minEnergy = 80.0; // MeV
    double maxEnergy = 180.0; // MeV
    // The x axis should be the mass of the pi0 going to gamma gamma in GeV
    TH1F* hist_mass = new TH1F("pi0_mass", "Reconstructed #pi^{0} Mass;m_{#pi^{0} (#gamma#gamma)} [MeV]; Number of #gamma#gamma", 100, minEnergy, maxEnergy);
    TH1F* hist_massCorrected = new TH1F("pi0_massCorrected", "Reconstructed #pi^{0} Mass;m_{#pi^{0} (#gamma#gamma)} [MeV]; Number of #gamma#gamma", 100, minEnergy, maxEnergy);
    // TH1F* hist_mass = new TH1F("pi0_mass", "Reconstructed #pi^{0} Mass;m_{#pi^{0} (#gamma#gamma)} [MeV]; Number of #gamma#gamma", 200, 30.0, 230.0);
    // TH1F* hist_massCorrected = new TH1F("pi0_massCorrected", "Reconstructed #pi^{0} Mass;m_{#pi^{0} (#gamma#gamma)} [MeV]; Number of #gamma#gamma", 200, 30.0, 230.0);

    // Histogram for Cluster Position
    double r_cal = 627.0;
    const int nbins = 59;

    TH2F* hxy_assoc = new TH2F("hxy_assoc", "Crystals Contributing to #pi^{0} Mass Reconstruction;x [mm];y [mm]", nbins, -r_cal, r_cal, nbins, -r_cal, r_cal);

    using Vec4 = ROOT::Math::PxPyPzEVector;
    using Vec3 = ROOT::Math::XYZVector;

    Long64_t nEntries = chain.GetEntries();
    std::cout << "Total events: " << nEntries << std::endl;

    int eventNumber = 1;

    // Vectors to store photon events
    std::vector<int>    vEvent;
    std::vector<int>    vPhoton;
    std::vector<int>    vCrystal;
    std::vector<float>  vX;
    std::vector<float>  vY;
    std::vector<float>  vPhotonEnergy;

    // Vectors to store two sigma photon events
    std::vector<int>    vEventTwoSigma;
    std::vector<int>    vPhotonTwoSigma;
    std::vector<int>    vCrystalTwoSigma;
    std::vector<float>  vXTwoSigma;
    std::vector<float>  vYTwoSigma;
    std::vector<float>  vPhotonEnergyTwoSigma;

    vEvent.reserve(nEntries);
    vPhoton.reserve(nEntries);
    vCrystal.reserve(nEntries);
    vX.reserve(nEntries);
    vY.reserve(nEntries);
    vPhotonEnergy.reserve(nEntries);

    // in MeV
    double massPDG = 134.9766; // PDG value for pi0 mass in MeV
    double sigma_sim = 4.42015; // Energy resolution for the calorimeter

    // double massRaw = 127.758; // Raw mass of the pi0 in MeV // For 18x275 GeV data set
    // double sigmaRaw = 4.42913; // Raw sigma of the pi0 in MeV // For 18x275 GeV data set
    double massRaw = 127.148; // Raw mass of the pi0 in MeV // For 5x41 GeV data set
    double sigmaRaw = 4.71922; // Raw sigma of the pi0 in MeV // For 5x41 GeV data set
    double upperThreeSigma = massPDG + 3*sigma_sim; // Upper limit for the three sigma range
    double lowerThreeSigma = massPDG - 3*sigma_sim; // Lower limit
    double upperTwoSigma = massPDG + 2*sigma_sim; // Upper limit for the two sigma range
    double lowerTwoSigma = massPDG - 2*sigma_sim; // Lower limit

    double lowerTwoSigmaRaw = massRaw - 2*sigmaRaw; // Lower limit for the two sigma range for raw mass
    double upperTwoSigmaRaw = massRaw + 2*sigmaRaw; // Upper limit for the two sigma range for raw mass
   
    // Load the correction factors from the file
    auto corrections = loadCorrections("calibration_constants_pi0_final_10it.csv");

    // Load crystal positions from the CSV file
    auto crystalPositions = loadCrystalPositions("crystal_id_xy.csv");

    Long64_t entryNum = 0;

    while (reader.Next()) {
        ++entryNum;

        if (entryNum % 100000 == 0) {
            std::cout << "\nProcessing event " << entryNum << std::endl;
            std::cout << "\rProgress: " << (entryNum * 100.0 / nEntries) << "%";
            std::cout.flush();
        }

        if (clusters->size() != 2) continue;

        auto c1 = clusters->at(0);
        auto c2 = clusters->at(1);

        Vec3 dir1(c1.position.x, c1.position.y, c1.position.z);
        Vec3 dir2(c2.position.x, c2.position.y, c2.position.z);

        dir1 = dir1.Unit();
        dir2 = dir2.Unit();

        Vec4 q1(c1.energy * dir1.X(), c1.energy * dir1.Y(), c1.energy * dir1.Z(), c1.energy);
        Vec4 q2(c2.energy * dir2.X(), c2.energy * dir2.Y(), c2.energy * dir2.Z(), c2.energy);

        Vec4 pi0 = q1 + q2;
        double mass = pi0.M()*1000; // MeV
        //double calibratedMass = mass * 1.05758; // Calibration factor, adjust as needed //Other calibration factors: 1.05758, 1.0595, 1.05374

        hist_mass->Fill(mass);

        // if (mass < lowerTwoSigmaRaw || mass > upperTwoSigmaRaw) {
        //     continue; // Skip events outside the energy range
        // }

        if (mass < 80 || mass > 180) {
            continue; // Skip events outside the energy range
        }

        // Find the best matching crystal IDs for the clusters
        int matched_id_1 = findClosestCrystalID(c1.position.x, c1.position.y, crystalPositions);
        int matched_id_2 = findClosestCrystalID(c2.position.x, c2.position.y, crystalPositions);


        // Your existing matched crystal IDs
        int matched_id[2] = { matched_id_1, matched_id_2 };

        // Hit positions from clusters
        float x[2] = { c1.position.x, c2.position.x };
        float y[2] = { c1.position.y, c2.position.y };

        // Energy from clusters
        float photonEnergy[2] = { c1.energy, c2.energy };

        // Append one row per photon
        for (int j = 0; j < 2; ++j) {
            vEvent.push_back(eventNumber);     
            vPhoton.push_back(j + 1);           // Photon number: 1 or 2
            vCrystal.push_back(matched_id[j]);
            vX.push_back(x[j]);
            vY.push_back(y[j]);
            vPhotonEnergy.push_back(photonEnergy[j]);
        }

        eventNumber++;  // Only increment after valid event

        // Check if the mass is within the two sigma range
        if (mass < lowerTwoSigmaRaw || mass > upperTwoSigmaRaw) {
            continue; // Skip events outside the two sigma range
        }

        // Store two sigma photon events
        for (int j = 0; j < 2; ++j) {
            vEventTwoSigma.push_back(eventNumber);
            vPhotonTwoSigma.push_back(j + 1); // Photon number: 1 or 2
            vCrystalTwoSigma.push_back(matched_id[j]);
            vXTwoSigma.push_back(x[j]);
            vYTwoSigma.push_back(y[j]);
            vPhotonEnergyTwoSigma.push_back(photonEnergy[j]);
        }

        // Apply the correction to the photon energy
        float correctedEnergy1 = applyCorrection(photonEnergy[0], matched_id[0], corrections);
        float correctedEnergy2 = applyCorrection(photonEnergy[1], matched_id[1], corrections);

        // Fill the corrected mass histogram
        Vec4 corrected_q1(correctedEnergy1 * dir1.X(), correctedEnergy1 * dir1.Y(), correctedEnergy1 * dir1.Z(), correctedEnergy1);
        Vec4 corrected_q2(correctedEnergy2 * dir2.X(), correctedEnergy2 * dir2.Y(), correctedEnergy2 * dir2.Z(), correctedEnergy2);
        Vec4 corrected_pi0 = corrected_q1 + corrected_q2;
        double correctedMass = corrected_pi0.M() * 1000; // MeV

        hist_massCorrected->Fill(correctedMass);

        // Fill the histogram of crystal positions
        hxy_assoc->Fill(c1.position.x, c1.position.y);
        hxy_assoc->Fill(c2.position.x, c2.position.y);
    }

    // Write the CSV file with all photon events
    WritePhotonEventCSV("photon_output_18x275_25.04.1_run21-22.csv", vEvent, vPhoton, vCrystal, vX, vY, vPhotonEnergy);
    WritePhotonEventCSV("photon_output_18x275_TwoSigma_25.04.1_run21-22.csv", vEventTwoSigma, vPhotonTwoSigma, vCrystalTwoSigma, vXTwoSigma, vYTwoSigma, vPhotonEnergyTwoSigma);

    TCanvas* c = new TCanvas("c", "Pi0 Invariant Mass", 800, 600);
    hist_mass->Draw();
    // No statistics box
    hist_mass->SetStats(kFALSE);
    hist_mass->SetLineColor(kBlack);
    hist_massCorrected->Draw("same");
    hist_massCorrected->SetLineColor(kBlue);


    // Draw a vertical line at the pi0 mass in MeV
    TLine* line = new TLine(massPDG, 0, massPDG, hist_massCorrected->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineStyle(2); // Dashed line
    line->SetLineWidth(2);
    line->Draw("same");

    // Draw vertical lines at the two sigma range
    TLine* lineLowerTwoSigma = new TLine(lowerTwoSigmaRaw, 0, lowerTwoSigmaRaw, hist_mass->GetMaximum());
    TLine* lineUpperTwoSigma = new TLine(upperTwoSigmaRaw, 0, upperTwoSigmaRaw, hist_mass->GetMaximum());
    lineLowerTwoSigma->SetLineColor(kGreen+2);
    lineLowerTwoSigma->SetLineStyle(2); // Dashed line
    lineLowerTwoSigma->SetLineWidth(2);
    lineLowerTwoSigma->Draw("same");
    lineUpperTwoSigma->SetLineColor(kGreen+2);
    lineUpperTwoSigma->SetLineStyle(2); // Dashed line
    lineUpperTwoSigma->SetLineWidth(2);
    lineUpperTwoSigma->Draw("same");

    // Fit with background
    TF1 *fitFunc = new TF1("fitFunc", "gaus(0) + pol1(3)", lowerThreeSigma, upperThreeSigma);
    // Set initial parameters, p1 = massPDG, p2 = sigma = 0.00437, p3 = background slope = 326.566
    fitFunc->SetParameters(288.939, massPDG, 4.32389, 419.771, -2.53081); //1000, 
    hist_massCorrected->Fit("fitFunc", "R");
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("same");
    double fitMean = fitFunc->GetParameter(1);
    double fitMeanError = fitFunc->GetParError(1);
    double fitSigma = fitFunc->GetParameter(2);
    double fitSigmaError = fitFunc->GetParError(2);

    // Add a legend to differentiate histograms
    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(hist_mass, "Raw m_{#pi^{0}} Distribution", "l");
    legend->AddEntry(hist_massCorrected, "Calibrated m_{#pi^{0}} Distribution", "l");
    legend->AddEntry(line, "m_{#pi^{0}} (PDG)", "l");
    legend->AddEntry(fitFunc, "Fit", "l");
    legend->SetBorderSize(0);
    legend->SetFillColor(0); // Transparent background
    legend->Draw();
    
    //Show peak position and sigma on the plot, also errors
    TLatex* fitText1 = new TLatex(0.15, 0.85,
    Form("m_{#pi^{0}} = %.3f#pm %.3f MeV", fitMean, fitMeanError));
    fitText1->SetNDC();
    fitText1->SetTextColor(kRed);
    fitText1->SetTextSize(0.03);
    fitText1->Draw();

    TLatex* fitText2 = new TLatex(0.15, 0.80,
        Form("#sigma_{m} = %.3f#pm %.3f MeV", fitSigma, fitSigmaError));
    fitText2->SetNDC();
    fitText2->SetTextColor(kRed);
    fitText2->SetTextSize(0.03);
    fitText2->Draw();

    TCanvas* c2 = new TCanvas("c2", "Crystals Contributing to Pion Mass Reconstruction", 1000, 1000);
    gStyle->SetPalette(55);
    hxy_assoc->GetXaxis()->SetTitle("X [mm]");
    hxy_assoc->GetYaxis()->SetTitle("Y [mm]");
    hxy_assoc->GetXaxis()->SetLabelSize(0.025);
    hxy_assoc->GetYaxis()->SetLabelSize(0.025);
    hxy_assoc->GetXaxis()->SetTitleSize(0.03);
    hxy_assoc->GetYaxis()->SetTitleSize(0.03);
    hxy_assoc->GetZaxis()->SetLabelSize(0.025);
    hxy_assoc->GetZaxis()->SetTitleSize(0.03);
    hxy_assoc->GetXaxis()->SetTitleOffset(1.2);
    hxy_assoc->GetYaxis()->SetTitleOffset(1.2);
    hxy_assoc->GetZaxis()->SetTitleOffset(0.7);
    hxy_assoc->Draw("COLZ");
    hxy_assoc->SetMinimum(1);
    gPad->SetLogz();

    // Cross section for SIDIS at 18x275 GeV // 5.36877564e+10 picobarns
    // root -l -b -q root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/SIDIS/pythia6-eic/1.0.0/18x275/q2_0to1/pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run9.ab.hepmc3.tree.root -e 'hepmc3_tree->Scan("hepmc3_event.attribute_string", "hepmc3_event.attribute_id==0 && hepmc3_event.attribute_name ==\"GenCrossSection\"", "colsize=30", 1, hepmc3_tree->GetEntries()-1)'
    float instLum = 1.0e33; // cm^-2 s^-1 // Expected instantaneous luminosity for the experiment
    //float SIDIS_cross_section_18x275 = 5.36877564e+10; // pb // Cross section for SIDIS at 18x275 GeV // Data Set 9
    float SIDIS_cross_section_18x275 = 7.42519148e5; // pb // Cross section for SIDIS at 18x275 GeV needs reference!!!
    float SIDIS_cross_section_18x275_cm2 = SIDIS_cross_section_18x275 * 1.0e-36; // Convert to cm^2
    float NumberOfEvents = nEntries; // Total number of events processed

    float intLum = NumberOfEvents / SIDIS_cross_section_18x275_cm2; // cm^-2
    float t0 = intLum / instLum; // seconds

    // Output results
    std::cout << "Cross Section for SIDIS at " << collisionEnergy << " GeV: " << SIDIS_cross_section_18x275_cm2 << " cm^2" << std::endl;
    std::cout << "Integrated Luminosity: " << intLum << " cm^-2" << std::endl;
    std::cout << "Instantaneous Luminosity: " << instLum << " cm^-2 s^-1" << std::endl;
    std::cout << "Total number of DIS events: " << NumberOfEvents << std::endl;
    std::cout << "Time to record data set: " << t0 << " seconds" << std::endl;

    // Skip not crystal bins
    std::set<std::pair<double, double>> emptyBinCoords;

    // --- Load CSV file ---
    std::ifstream inFile("empty_bins.csv");
    std::string line_read;
    std::getline(inFile, line_read);  // Skip header

    while (std::getline(inFile, line_read)) {
        std::stringstream ss(line_read);
        std::string x_str, y_str;
        std::getline(ss, x_str, ',');
        std::getline(ss, y_str, ',');

        double x = std::stod(x_str);
        double y = std::stod(y_str);
        emptyBinCoords.insert({x, y});
    }
    inFile.close();

    double tolerance = 1e-2;  // tolerance in mm

    auto isInEmptyList = [&](double x, double y) {
        for (const auto& [ex, ey] : emptyBinCoords) {
            if (std::abs(ex - x) < tolerance && std::abs(ey - y) < tolerance)
                return true;
        }
        return false;
    };

    // Total number of crytals
    //int totalCrystals = 59 * 59 - emptyBinCoords.size(); // Total number of crystals minus empty bins
    int totalCrystals = 2740; // Total number of crystals in the EcalEndcapN detector 
    //std::cout << "Total number of crystals: " << totalCrystals << std::endl;

    // Empty bins
    int emptyBins = emptyBinCoords.size();

    // Compute empty crystals
    int emptyCrystals = 0;

    // Create a time color map using the hxy_assoc histogram, which contains the positions of clusters associated with electrons
    TH2F *timeMap = (TH2F*)hxy_assoc->Clone("timeMap");
    timeMap->SetTitle("Time in Hours to Calibrate Crystals");
    timeMap->Reset();  // Clear contents

    for (int ix = 1; ix <= hxy_assoc->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hxy_assoc->GetNbinsY(); ++iy) {
            float binEntries = hxy_assoc->GetBinContent(ix, iy);
            double binx = hxy_assoc->GetXaxis()->GetBinCenter(ix);
            double biny = hxy_assoc->GetYaxis()->GetBinCenter(iy);

            float t;
            if (binEntries >= 10) {
                t = (1000.0 * t0) / binEntries / 3600.0;  // in hours
            } else if (isInEmptyList(binx, biny)) {
                t = 0.0;  // Empty bin, set to 0.0
            } else {
                t = 0.0;  // No data or less than 10 entries
                emptyCrystals++;
            }

            timeMap->SetBinContent(ix, iy, t);
        }
    }

    std::cout << "Percentage of calibrated crystals: " 
              << (100.0 * (totalCrystals - emptyCrystals) / totalCrystals) << "%" << std::endl;

    TCanvas *c4 = new TCanvas("c4", "Time Map", 1000, 1000);
    timeMap->GetXaxis()->SetTitle("X [mm]");
    timeMap->GetYaxis()->SetTitle("Y [mm]");
    timeMap->GetXaxis()->SetLabelSize(0.025);
    timeMap->GetYaxis()->SetLabelSize(0.025);
    timeMap->GetXaxis()->SetTitleSize(0.03);
    timeMap->GetYaxis()->SetTitleSize(0.03);
    timeMap->GetZaxis()->SetLabelSize(0.025);
    timeMap->GetZaxis()->SetTitleSize(0.03);
    timeMap->GetXaxis()->SetTitleOffset(1.2);
    timeMap->GetYaxis()->SetTitleOffset(1.2);
    timeMap->GetZaxis()->SetTitleOffset(0.7);
    timeMap->Draw("COLZ");
    gPad->SetLogz();
    // TLatex timeLatex;
    // timeLatex.SetNDC(); // Normalized device coordinates (0–1)
    // timeLatex.SetTextAlign(22); // Centered alignment
    // timeLatex.SetTextSize(0.030); // Adjust size as needed
    // timeLatex.DrawLatex(0.5, 0.93, "Time in Hours to Calibrate Crystals");

    timeMap->SetMinimum(1e-1); // Set minimum for log scale
    //timeMap->SetMaximum(1e6); // Set maximum for log scale
    timeMap->SetStats(false);

    timer.Stop();
    timer.Print();
}
