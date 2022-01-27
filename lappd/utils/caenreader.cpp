#include <algorithm>
#include <cstdlib>
#include <list>
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
    static const float NS_PER_SAMPLE = 0.2;

    // Container for header as output by CAEN Wavedump
    struct CAENHeader
    {
        unsigned int eventSize;  // Should be 4120 for x1742 models (6*4byte header + 1024*4byte record)
        unsigned int boardID;    // Not used
        unsigned int pattern;    // Not used
        unsigned int channel;    // This is the channel within the group i.e. ch8 is 0 but in group2
        unsigned int eventCount; // Not exact event count as it rolls over whenever the limit is reached
        unsigned int eventTime;  // Time when event is created in digitiser memory, not trigger time
    };

    struct CAENWave
    {
        std::vector<unsigned int> header;
        std::vector<float> wave;
        int eventNo;
    };

    struct CAENProcessedWave
    {
        CAENHeader header;
        std::vector<float> wave;
        int eventNo;
        float max;
        float peakHeight;
        int peakSample;
    };

    struct WavePair
    {
        CAENProcessedWave leftWave;
        CAENProcessedWave rightWave;
    };

    inline CAENHeader processHeader(std::vector<unsigned int> headerVect)
    {
        return CAENHeader{
            headerVect[0],
            headerVect[1],
            headerVect[2],
            headerVect[3],
            headerVect[4],
            headerVect[5],
        };
    }

    float median(std::vector<float> wave)
    {
        assert(!wave.empty());
        if (wave.size() % 2 == 0)
        {
            const auto median_it1 = wave.begin() + wave.size() / 2 - 1;
            const auto median_it2 = wave.begin() + wave.size() / 2;

            std::nth_element(wave.begin(), median_it1, wave.end());
            const auto e1 = *median_it1;

            std::nth_element(wave.begin(), median_it2, wave.end());
            const auto e2 = *median_it2;

            return (e1 + e2) / 2;
        }
        else
        {
            const auto median_it = wave.begin() + wave.size() / 2;
            std::nth_element(wave.begin(), median_it, wave.end());
            return *median_it;
        }
    }

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
        float dev = ROOT::VecOps::StdDev(wave);
        float mean = ROOT::VecOps::Mean(wave);
        auto slicedWave = wave[abs(wave) < (mean + dev)];
        float slicedMean = ROOT::VecOps::Mean(slicedWave);
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
        // Strip last samples, baseline subtract and convert to mV
        removeLastNSamples(wave, 10);
        auto baseline = calculateBaseline(wave);
        subtractBaseline(wave, baseline);
        convertADCWave(wave, 12);
    }

    // Maybe overload this method such that it runs preprocessWave?
    CAENProcessedWave processWave(CAENWave &caenWave)
    {
        ROOT::RVec<float> rWave = caenWave.wave;
        int wavePeak = ROOT::VecOps::ArgMin(rWave);
        float waveMax = ROOT::VecOps::Max(rWave);
        float waveMin = caenWave.wave[wavePeak];
        CAENHeader processedHeader = processHeader(caenWave.header);
        return CAENProcessedWave{
            processedHeader,
            caenWave.wave,
            caenWave.eventNo,
            waveMax,
            waveMin,
            wavePeak};
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

    std::vector<float> sliceAroundPeak(CAENProcessedWave wave, float lookback, float lookforward)
    {
        int backSamples = std::floor(lookback / NS_PER_SAMPLE);
        int forwardSamples = std::floor(lookforward / NS_PER_SAMPLE);
        int windowStart = wave.peakSample - backSamples;
        int windowEnd = wave.peakSample + forwardSamples;
        if (windowStart < 0 || windowEnd > (int)wave.wave.size())
        {
            std::cout << "Invalid slice, returning original waveform" << std::endl;
            return wave.wave;
        }
        std::vector<float> slicedWave = {wave.wave.begin() + wave.peakSample - backSamples, wave.wave.begin() + wave.peakSample + forwardSamples};
        return slicedWave;
    }

    std::vector<float> sliceAroundPeak(CAENWave wave, int peakSample, float lookback, float lookforward)
    {
        int backSamples = std::floor(lookback * NS_PER_SAMPLE);
        int forwardSamples = std::floor(lookforward * NS_PER_SAMPLE);
        std::vector<float> slicedWave = {wave.wave.begin() + peakSample - backSamples, wave.wave.begin() + peakSample + forwardSamples};
        return slicedWave;
    }

    bool correlateAcrossStrip(CAENProcessedWave firstWave, std::string secondFilename, float peakHeightReq = 0.8)
    {
        CAENWave secondWave = readCAENWave(secondFilename, firstWave.eventNo);
        std::vector<float> slicedSecond = sliceAroundPeak(secondWave, firstWave.peakSample, 2, 2);
        return passMinThreshold(slicedSecond, firstWave.peakHeight * peakHeightReq);
    }

    std::vector<CAENProcessedWave> readAllChannels(std::list<std::string> filenames, int eventNo)
    {
        std::vector<CAENProcessedWave> waves;
        for (auto file : filenames)
        {
            CAENWave wave = readCAENWave(file, eventNo);
            waves.push_back(processWave(wave));
        }
        return waves;
    }

    float integrate(std::vector<float> wave)
    {
        float total = 0;
        for (auto value : wave)
        {
            total += value;
        }
        return abs(total);
    }

    float integratedCharge(std::vector<float> wave, float terminationOhms = 50.0)
    {
        float charge = integrate(wave) * NS_PER_SAMPLE / terminationOhms;
        return charge;
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