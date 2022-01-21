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
        std::vector<unsigned int> header;
        std::vector<float> wave;
        int eventNo;
    };

    struct CAENProcessedWave
    {
        std::vector<unsigned int> header;
        std::vector<float> wave;
        int eventNo;
        float max;
        float min;
        int peakSample;
    };

    int getNumberEntries(std::string filename)
    {
        std::ifstream infile(filename, std::ios::binary);
        infile.seekg(0, std::ios::end);
        auto size = infile.tellg();
        int nWaves = size / (sizeof(float) * 1030);
        return nWaves;
    }

    inline float mv_to_adc(float millivolts, int bitR)
    {
        float resolution = pow(2, bitR) - 1;
        return resolution * millivolts / 1000.0;
    }

    inline float adc_to_mv(float adc, int bitR)
    {
        float resolution = pow(2, bitR) - 1;
        return adc / resolution * 1000.0;
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

    CAENWave readCAENWave(std::ifstream &infile)
    {
        auto count = infile.tellg();
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
        int eventNo = count / (sizeof(float) * 1030);
        return CAENWave{
            header,
            wave,
            eventNo};
    }

    CAENWave readCAENWave(std::string filename, int eventNo)
    {
        int bytestart = eventNo * sizeof(float) * 1030;
        std::ifstream infile(filename, std::ios::binary);
        infile.seekg(bytestart);
        return readCAENWave(infile);
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

    inline void subtractBaseline(std::vector<float> &wave)
    {
        // Subtract baseline using default calculateBaseline (reduced mean) method
        auto baseline = calculateBaseline(wave);
        std::transform(wave.begin(), wave.end(), wave.begin(), [baseline](float &value) -> float
                       { return value - baseline; });
    }

    inline void subtractBaseline(std::vector<float> &wave, float baseline)
    {
        // Subtract a given baseline from a wave
        std::transform(wave.begin(), wave.end(), wave.begin(), [baseline](float &value) -> float
                       { return value - baseline; });
    }

    inline void convertADCWave(std::vector<float> &ADCWave, int bitR)
    {
        std::transform(ADCWave.begin(), ADCWave.end(), ADCWave.begin(), [](float &value) -> float
                       { return adc_to_mv(value, 12); });
    }

    inline void removeLastNSamples(std::vector<float> &wave, int nSamples)
    {
        wave.erase(wave.end() - nSamples, wave.end());
    }

    void preprocessWave(std::vector<float> &wave)
    {
        removeLastNSamples(wave, 10);
        auto baseline = calculateBaseline(wave);
        subtractBaseline(wave, baseline);
        convertADCWave(wave, 12);
    }

    inline bool passMinThreshold(std::vector<float> wave, float threshold)
    {
        for (auto value : wave)
        {
            if (value < threshold)
            {
                return true;
            }
        }
        return false;
    }

    inline bool passMaxThreshold(std::vector<float> wave, float threshold)
    {
        for (auto value : wave)
        {
            if (value > threshold)
            {
                return false;
            }
        }
        return true;
    }

    struct darkOutput
    {
        std::vector<CAENWave> waves;
        int rejectedMin;
        int rejectedMax;
    };

    darkOutput coarseDarkSearch(std::string filepath, float minThreshold, float maxThreshold)
    {
        struct darkOutput output;
        std::vector<CAENWave> darkWaves;
        int count = 0;
        int passed = 0;
        int failedMin = 0;
        int failedMax = 0;
        std::ifstream infile(filepath, std::ios::binary);
        infile.seekg(0, std::ios::end);
        auto size = infile.tellg();
        // std::cout << "Number of waves: " << size / (sizeof(float) * 1030) << std::endl; // (1024+6 )* 32 byte * 8 bits/byte
        infile.seekg(0, std::ios::beg);
        while ((infile.tellg() < size))
        {
            //auto rawWave = readOne(infile);
            auto caenWave = readCAENWave(infile);
            preprocessWave(caenWave.wave);
            count++;
            if (!passMaxThreshold(caenWave.wave, maxThreshold))
            {
                failedMax++;
                continue;
            }
            if (!passMinThreshold(caenWave.wave, minThreshold))
            {
                failedMin++;
                continue;
            }
            passed++;
            darkWaves.push_back(caenWave);
        }
        std::cout << failedMin << " failed low threshold" << std::endl;
        std::cout << failedMax << " failed high threshold" << std::endl;
        output.waves = darkWaves;
        output.rejectedMin = failedMin;
        output.rejectedMax = failedMax;
        return output;
    }

}