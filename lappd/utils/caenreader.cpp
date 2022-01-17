#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "TROOT.h"
#include "ROOT/RVec.hxx"

// I don't know how to write c++, so
// this code is probably horrendous.
// Compile with .L caenreader.cpp+ in ROOT command line,
// or with g++ but I don't understand makefiles

namespace CAENReader
{

    struct CAENWave
    {
        std::vector<float> header;
        std::vector<float> wave;
        float baseline;
    };

    inline float mv_to_adc(float millivolts, int bitR)
    {
        int resolution = pow(bitR, 2) - 1;
        return resolution * millivolts / 1000.0;
    }

    inline float adc_to_mv(float adc, int bitR)
    {
        int resolution = pow(bitR, 2) - 1;
        return adc / resolution * 1000;
    }

    std::vector<float> readOne(std::ifstream &infile)
    {
        unsigned int hBuffer = 0;
        std::vector<unsigned int> header;
        for (int i = 0; i < 6; i++)
        {
            infile.read((char *)&hBuffer, sizeof(unsigned int));
            header.push_back(hBuffer);
        }
        float wBuffer = 0;
        std::vector<float> wave;
        for (int i = 0; i < 1024; i++)
        {
            infile.read((char *)&wBuffer, sizeof(float));
            wave.push_back(wBuffer);
        }
        return wave;
    }

    std::vector<float> readOne(std::string filename, int eventNo)
    {
        int byteStart = eventNo * sizeof(float) * 1030;
        std::ifstream infile(filename, std::ios::binary);
        infile.seekg(byteStart);
        auto wave = readOne(infile);
        return wave;
    }

    inline std::vector<float> subtractBaseline(std::vector<float> wave, float baseline)
    {
        std::transform(wave.begin(), wave.end(), wave.begin(), [baseline](float &value) -> float
                       { return value - baseline; });
        return wave;
    }

    inline std::vector<float> convertADCWave(std::vector<float> ADCWave, int bitR)
    {
        std::transform(ADCWave.begin(), ADCWave.end(), ADCWave.begin(), [](float &value) -> float
                       { return adc_to_mv(value, 12); });
        return ADCWave;
    }

    inline float calculateBaseline(ROOT::RVec<float> wave)
    {
        float dev = ROOT::VecOps::StdDev(wave); //TMath::StdDev(wave.begin(), wave.end());
        float mean = ROOT::VecOps::Mean(wave);  //TMath::Mean(wave.begin(), wave.end());
        // int count = 0;
        // float total = 0;
        // for (auto value : wave)
        // {
        //     if (TMath::Abs(value) < mean + dev)
        //     {
        //         total += value;
        //         count++;
        //     }
        // }
        auto slicedWave = wave[abs(wave) < (mean + dev)];
        float slicedMean = ROOT::VecOps::Mean(slicedWave); //TMath::Mean(slicedWave.begin(), slicedWave.end());
        // float redMean = total / float(count);
        return slicedMean;
    }

    std::vector<float> processWave(std::vector<float> wave)
    {
        auto baseline = calculateBaseline(wave);
        auto subWave = subtractBaseline(wave, baseline);
        auto voltWave = convertADCWave(subWave, 12);
        return voltWave;
    }

    inline bool passThreshold(std::vector<float> wave, float threshold)
    {
        for (auto &value : wave)
        {
            if (value < threshold)
            {
                return true;
            }
        }
        return false;
    }

    std::vector<std::vector<float>> coarseDarkSearch(std::string filepath, float threshold)
    {
        std::ifstream infile(filepath, std::ios::binary);
        infile.seekg(0, std::ios::end);
        auto size = infile.tellg();
        std::cout << "Number of waves: " << size / (sizeof(float) * 1030) << std::endl; // (1024+6 )* 32 byte * 8 bits/byte
        infile.seekg(0, std::ios::beg);
        std::vector<std::vector<float>> darkWaves;
        int count = 0;
        int passed = 0;
        while ((infile.tellg() < size))
        {
            auto rawWave = readOne(infile);
            auto wave = processWave(rawWave);
            count++;
            if (passThreshold(wave, threshold))
            {
                passed++;
                darkWaves.push_back(wave);
            }
        }
        return darkWaves;
    }

}