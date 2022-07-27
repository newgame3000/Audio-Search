#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <math.h>
#include <vector>
#include <complex>
#include <cmath>  


using namespace std;

const double PI = 3.141592653589793238463;
vector<int> RANGE_FREQ {40, 80, 120, 180, 300};
int WINSIZE = 4096;

void WindwoHann(vector<complex<double>> &signal, int size) {
    for (int i = 0; i < size; i++) {
        double multiplier = 0.5 * (1 - cos(2* PI * i / (size- 1)));
        signal[i] = multiplier * signal[i];
    }
}

void Gausse(vector<complex<double>> &signal, int size) {
    for (int i = 0; i < size; i++) {
        int a = (size - 1) / 2;
        double t = (i - a) / (0.5 * a);
        t = t*t;
        signal[i] *= exp(-t/2);
    }
}


int Reverse(int value, int lg_n) { //Реверс последних lg_n битов в числе
    int res = 0;
    while(lg_n > 0) {
        res = res << 1;
        res = res | (value & 1);
        value = value >> 1;
        lg_n -= 1;
    }
    return res;
}

void FFT(vector<complex<double>> & a, const int &n) {
    int lg_n = 0;
    while ((1 << lg_n) < n) { 
        ++lg_n;
    }

    for (int i = 0; i < n; ++i) {
        int fr = Reverse(i, lg_n);
        if (i < fr) {
            swap (a[i], a[fr]);
        }
    }

    for (int k = 1; k <= lg_n; ++k) { //log_n раз
        int len = 1 << k;
        for (int i = 0; i < n; i += len) { //идём по всем нижним блокам

            complex<double> wz(cos(2 * PI / len), sin(2 * PI / len));
            complex<double> w(1, 0);

            for (int j = 0; j < len / 2; j++) {//по верхнему блоку
                complex<double> u = a[i + j];
                complex<double> v = a[i + j + len / 2] * w; //Сливаем
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v; 
                w *= wz;
            }
        }
    }
}


struct TWaveHeader {
    char ChunkId[4];                    // Информация о формате файла (RIFF), Содержит символы “RIFF” в ASCII кодировке;
    uint32_t ChunkSize;                 // Размер без  chunkId[4];
    char Format[4];                     // Формат потоковых данных (WAVE);
    char SubchunkId[4];                 // Описание параметров WAV-файла (fmt-chunk);
    uint32_t SubchunkSize;              // Размер подструктуры  subchunk1  (16 байт);
    uint16_t WFomatTag;                 // Аудио формат (PCM = 1);
    uint16_t NChannels;                 // Количество каналов (Моно = 1, Стерео = 2);
    uint32_t SamplesPerSec;             // Частота дискретизации в Гц;
    uint32_t ByteRate;                  // Кол-во передаваемых байт в секунду воспроизведения;
    uint16_t BlockAlign;                // Размер сэмпла в байтах 16 бит = 2 байта моно, 32 бита = 4 байта стерео (включая все каналы);
    uint16_t BitsPerSample;             // Количество бит в сэмпле. Так называемая “глубина” или точность звучания. 8 бит, 16 бит и т.д. /// битов на отсчет
};

struct Chunk
{
    char ID[4]; //"data" = 0x61746164
    uint32_t size;  //Chunk data bytes
};

void writeGraphicsFFT(FILE *gr_f, vector<vector<complex<double>>> &analysis) {
    vector<double> magnitude;
    vector<double> phase;

    for (uint i = 0; i < analysis.size(); ++i) {
        for (uint j = 0; j < analysis[i].size(); ++j) {
            magnitude.push_back(sqrt(analysis[i][j].real() * analysis[i][j].real() + analysis[i][j].imag() * analysis[i][j].imag()));
            phase.push_back(atan2(analysis[i][j].imag(), analysis[i][j].real()));
        }
    }

    for (uint i = 0; i < magnitude.size(); ++i) {

        string name = to_string(magnitude[i]) + "\n";
        fprintf(gr_f, name.c_str());

        name = to_string(phase[i]) + "\n";
        fprintf(gr_f, name.c_str());
    }
}

vector<vector<complex<double>>> fileWindiowFFT(FILE * file, FILE *gr, FILE *gr_f, bool write_graphics) {

    TWaveHeader header;
    fread(&header, sizeof(TWaveHeader), 1, file);

    Chunk chunk;
    while (true) {
        fread(&chunk, sizeof(chunk), 1, file);
        if (*(unsigned int*)&chunk.ID == 0x61746164) { //data
            break;
        }
        fseek(file, chunk.size, SEEK_CUR);
    }

    short point;
    vector<complex<double>> signal;
    vector<double> graphic;

    int k = 0;
    double mean = 0;
    int window = 0;
    
    vector<vector<complex<double>>> analysis;

    while (chunk.size > 0) {
        k += 1;
        fread(&point, sizeof(uint16_t), 1, file);
        double p = double(point);

        chunk.size -= sizeof(uint16_t);

        if (header.NChannels  == 2) {
            if (k % 2 == 1) {
                mean = p;
                continue;
            }
            p = (p + mean) / 2;
        }

        window += 1;

        signal.push_back(p);       

        if (window == WINSIZE) {
            window = 0;
            Gausse(signal, WINSIZE);
            analysis.push_back(signal);
            signal.clear();
        }

        if (write_graphics) {
            graphic.push_back(p);
            string name = to_string(graphic[graphic.size() - 1]) + "\n";
            fprintf(gr, name.c_str());
        }
    }

    while(analysis[analysis.size() - 1].size() < (uint)WINSIZE) {
        analysis[analysis.size() - 1].push_back(complex<double>(0));
    }   

    for (uint i = 0; i < analysis.size(); ++i) {
        FFT(analysis[i], analysis[i].size());
    }
    
    if(write_graphics) {
        writeGraphicsFFT(gr_f, analysis);
    }    

    return analysis;
}


int getIndex(double curfreq) {

    if (curfreq < 30 || curfreq > RANGE_FREQ[RANGE_FREQ.size() - 1]) {
        return -1;
    }

    for (uint i = 0; i < RANGE_FREQ.size(); ++i) {
        if (curfreq < RANGE_FREQ[i]) {
            return (int)i;
        }
    }

    return -1;
}

vector<vector<double>> analysisOfFile(vector<vector<complex<double>>> analysis, int size) {

    vector<vector<double>> freqs(analysis.size(), vector<double>(RANGE_FREQ.size()));
    for (uint i = 0; i < analysis.size(); ++i) { //Проходимся по всем окнам
        double curfreq = 0;
        vector<double> max(RANGE_FREQ.size());
        for (uint j = 0; j < analysis[i].size() / 2; ++j) { //Нас интересует только половина

            curfreq += 44100 / size;
            int ind = getIndex(curfreq);
            if (ind == -1) {
                continue;
            }

            if (floor(log2(abs(analysis[i][j])) > max[ind])) {
                freqs[i][ind] = floor(log2(abs(analysis[i][j])));
            }
        }
    }

    return freqs;
}


int Search(vector<vector<double>> &freqsSearch, vector<vector<double>> &freqsVariant) {
    int ans = 0;

    for (uint i = 0; i < freqsVariant.size(); ++i) { //По окнам одного
        for (uint k = 0; k < freqsSearch.size(); ++k) { //По окнам второго

            bool flag = true;

            for (uint j = 0; j < freqsVariant[i].size(); ++j) {
            
                if(freqsVariant[i][j] != freqsSearch[k][j]) {
                    flag = false;
                    break;
                }
            }

            if (flag) {
                ans += 1;
            } 
        }
    }

    return ans;
}


int main(int argc, char* argv[]) {
    string name = "/home/kirill/FFT/audio/patterns/";
    name += argv[1];

    FILE *file;
    if ((file = fopen(name.c_str(), "rb")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    
    FILE *gr;
    if ((gr = fopen("gr.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    FILE *gr_f;
    if ((gr_f = fopen("gr_f.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }
    
    vector<vector<complex<double>>> windowsFFTSearch = fileWindiowFFT(file, gr, gr_f, true);
    vector<vector<double>> freqsSearch = analysisOfFile(windowsFFTSearch, WINSIZE);

    fclose(gr);
    fclose(gr_f);
    fclose(file);


    if ((file = fopen(argv[2], "rb")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    vector<string> searchFilenames;
    string filename = "";

    char c;

    while(!feof(file)) {
        fread(&c, sizeof(c), 1, file);
        if(c == '\n' && filename != "") {
            searchFilenames.push_back(filename);
            filename = "";
        } else {
            filename += c;
        }
    }

    fclose(file);

    int rightIndex = -1;
    int maxWin = 0;

    //cout << freqsSearch.size()<<endl;
    string searchPath = "/home/kirill/FFT/audio/searchFiles/";

    for (uint i = 0; i < searchFilenames.size(); ++i) {

        string filePath = searchPath + searchFilenames[i];

        if ((file = fopen(filePath.c_str(), "rb")) == NULL) {
            cout << "Ошибка открытия файла\n";
            return -1;
        }


        cout << "Анализ файла: " << searchFilenames[i] << endl;
        vector<vector<complex<double>>> windowsVariantFFT = fileWindiowFFT(file, gr, gr_f, false);

        if (windowsVariantFFT.size() < windowsFFTSearch.size()) {
            cout << "Количество совпадений: 0" << endl;
            continue;
        }

        vector<vector<double>> freqsVariant = analysisOfFile(windowsVariantFFT, WINSIZE);
        int ans = Search(freqsSearch, freqsVariant);
        cout << "Количество совпадений: " << ans << endl;
        if ((uint)ans > freqsSearch.size() / 2 && ans > maxWin) {
            maxWin = ans;
            rightIndex = i;
        }
        fclose(file);
    }


    if (rightIndex != -1) {
        cout << searchFilenames[rightIndex] << endl;
    } else {
        cout << "Не найдено!\n";
    }
}


