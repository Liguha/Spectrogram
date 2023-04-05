#include <iostream>
#include <bitset>
#include <fstream>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <algorithm>
#include <string>
#include <complex>

using namespace std;

const int INF = 1e9;
const double PI = acos(-1);

using ll = long long;
using uchar = unsigned char;
using cmd = complex<double>;

template <typename T>
void fft_r(vector <T>& p, vector <cmd>& res, cmd w, int l, int r)
{
    if (l == r)
    {
        int k = p.size();
        int n = 0;
        while (k > 1)
        {
            n++;
            k /= 2;
        }
        int c = 0;
        for (int i = 0; i < n; i++)
        {
            c = c | (l & 1);
            c = c << 1;
            l = l >> 1;
        }
        c = c >> 1;
        res[r] = p[c];
        return;
    }
    int m = (l + r) / 2;
    fft_r(p, res, w * w, l, m);
    fft_r(p, res, w * w, m + 1, r);
    int n = r - l + 1;
    int k = n / 2;
    cmd wt = 1;
    vector <cmd> v;
    for (int i = 0; i < n; i++)
    {
        cmd a = res[l + i % k] + wt * res[m + 1 + i % k];
        v.push_back(a);
        wt *= w;
    }
    for (int i = 0; i < n; i++)
        res[l + i] = v[i];
}

template <typename T>
vector <cmd> fft(vector <T>& v)
{
    int n = v.size();
    vector <cmd> out(n);
    fft_r(v, out, polar(1.0, 2 * acos(-1) / n), 0, n - 1);
    return out;
}

vector <double> rfft(vector <cmd>& v)
{
    int n = v.size();
    vector <cmd> tmp(n);
    fft_r(v, tmp, polar(1.0, -2 * acos(-1) / n), 0, n - 1);
    vector <double> out(n);
    for (int i = 0; i < n; i++)
        out[i] = tmp[i].real() / n;
    return out;
}

#pragma pack(push,1)
struct RGB
{
    uchar B;
    uchar G;
    uchar R;

    RGB(int _R = 0, int _G = 0, int _B = 0)
    {
        R = (uchar)_R;
        B = (uchar)_B;
        G = (uchar)_G;
    }
};
struct BMPheader
{
    char type[2] = { 'B','M' };
    __int32 size = -1;
    __int32 reserved = 0;
    __int32 offset = 54;
    __int32 headSize = 40;
    __int32 cols = -1;
    __int32 rows = -1;
    __int16 planes = 1;
    __int16 bitCount = 24;
    __int32 compression = 0;
    __int32 imageSize = -1;
    char otherParams[16];
};
#pragma pack(pop)

struct BMPimage
{
    BMPheader head;
    int w, h;
    vector <vector <RGB>> mat;
};

BMPimage load_image(string name)
{
    ifstream img(name, ios::binary);
    BMPimage out;
    img.read((char*)&out.head, sizeof(BMPheader));
    out.w = out.head.cols;
    out.h = out.head.rows;
    out.mat = vector <vector <RGB>>(out.h, vector<RGB>(out.w));
    for (int y = 0; y < out.h; y++)
    {
        for (int x = 0; x < out.w; x++)
        {
            RGB col;
            img.read((char*)&col, 3);
            out.mat[y][x] = col;
        }
        char del[4];
        img.read(del, out.w % 4);
    }
    img.close();
    return out;
}

BMPimage new_image(int rows, int cols)
{
    BMPimage img;
    img.w = cols;
    img.h = rows;
    img.head.rows = rows;
    img.head.cols = cols;
    img.mat = vector <vector <RGB>>(rows, vector<RGB>(cols));
    return img;
}

void save_image(BMPimage& img, string name)
{
    ofstream file(name, ios::binary);
    file.write((char*)&img.head, sizeof(BMPheader));
    for (int y = 0; y < img.h; y++)
    {
        for (int x = 0; x < img.w; x++)
            file.write((char*)&img.mat[y][x], sizeof(RGB));
        char zero = 0;
        for (int i = 0; i < img.w % 4; i++)
            file.write((char*)&zero, 1);
    }
    file.close();
}

#pragma pack(push,1)
struct WAVfmt
{
    char ID[4] = { 'f', 'm', 't', ' ' };
    __int32 size = 16;
    __int16 compression = 1;
    __int16 channels = 1;
    __int32 sampleRate = 44100;
    __int32 bytesPerSecond = 88200;
    __int16 blockAlign = 2;
    __int16 bitsPerSample = 16;
};
#pragma pack(pop)

struct WAVaudio
{
    WAVfmt fmt;
    int sampleRate;
    int channels;
    vector <vector <int>> audio;
};

WAVaudio load_audio(string name)
{
    char ID[5];
    ID[4] = '\0';
    ifstream file(name, ios::binary);
    WAVaudio out;
    int size;
    while (true)
    {
        file.read(ID, 4);
        file.read((char*)&size, 4);
        if (string(ID) == "RIFF")
        {
            file.read(ID, 4);
            continue;
        }
        if (string(ID) == "fmt ")
        {
            file.read((char*)&out.fmt + 8, sizeof(WAVfmt) - 8);
            out.fmt.size = size;
            out.channels = out.fmt.channels;
            out.sampleRate = out.fmt.sampleRate;
            continue;
        }
        if (string(ID) == "data")
            break;
        file.seekg((int)file.tellg() + size);
    }
    int n = out.fmt.blockAlign / out.fmt.channels;
    size /= out.fmt.blockAlign;
    out.audio = vector <vector <int>>(out.channels, vector<int>(size));
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < out.channels; j++)
            file.read((char*)&out.audio[j][i], n);
    }
    file.close();
    return out;
}

WAVaudio new_audio(int sampleRate, int channels)
{
    WAVaudio out;
    out.sampleRate = sampleRate;
    out.fmt.sampleRate = sampleRate;
    out.channels = channels;
    out.fmt.channels = channels;
    out.fmt.blockAlign = 2 * channels;
    out.fmt.bytesPerSecond = 2 * sampleRate * channels;
    out.audio = vector <vector <int>>(channels);
    return out;
}

void save_audio(WAVaudio& audio, string name)
{
    ofstream file(name, ios::binary);
    char ID[] = "RIFF";
    char WAVE[] = "WAVE";
    int size = 8 + sizeof(WAVfmt) + audio.audio[0].size() * audio.fmt.blockAlign;
    file.write(ID, 4);
    file.write((char*)&size, 4);
    file.write(WAVE, 4);
    file.write((char*)&audio.fmt, sizeof(WAVfmt));
    char IDdata[] = "data";
    size = audio.audio[0].size() * audio.fmt.blockAlign;
    file.write(IDdata, 4);
    file.write((char*)&size, 4);
    size = audio.audio[0].size();
    int n = audio.fmt.blockAlign / audio.channels;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < audio.channels; j++)
            file.write((char*)&audio.audio[j][i], n);
    }
    file.close();
}

int rgb2gray(RGB col)
{
    return 0.2989 * col.R + 0.5870 * col.G + 0.1140 * col.B;
}

// f: [0,1] -> [0,1]
vector <RGB> make_gradient(int n, RGB l, RGB r, double(*f)(double))
{
    vector <RGB> out(n);
    int dR = (int)r.R - l.R;
    int dG = (int)r.G - l.G;
    int dB = (int)r.B - l.B;
    for (int i = 0; i < n; i++)
    {
        double k = f((double)(i + 1) / n);
        out[i] = RGB(l.R + dR * k, l.G + dG * k, l.B + dB * k);
    }
    return out;
}

int nearest(int a)
{
    int res = 1;
    while (a > res)
        res *= 2;
    return res;
}

BMPimage make_spectrogram(WAVaudio& audio, int H, vector <RGB> gradient)
{
    H = 2 * nearest(H);
    int W = audio.audio[0].size() / H;
    vector <vector <double>> mat(H, vector <double>(W));
    double mx = 1;
    BMPimage img = new_image(H / 2, W);
    vector <double> a(H);
    vector <cmd> b(H);
    for (int x = 0; x < W; x++)
    {
        for (int i = 0; i < H; i++)
            a[i] = audio.audio[0][x * H + i];
        b = fft(a);
        for (int y = 0; y < H / 2; y++)
        {
            mat[y][x] = abs(b[y]);
            if (mat[y][x] > mx)
                mx = mat[y][x];
        }
    }

    double k = log(155);
    for (int y = 0; y < img.h; y++)
    {
        for (int x = 0; x < img.w; x++)
        {
            int idx = min((int)(gradient.size() * mat[y][x] / mx), (int)gradient.size() - 1);
            img.mat[y][x] = gradient[idx];
        }
    }
    return img;
}

WAVaudio make_noise(BMPimage& img, int rate)
{
    const double magnitude = 1e5;
    const int contrast = 30;
    const int black = 100;

    int H = 2 * nearest(img.h);
    vector <cmd> a(H);
    vector <double> b(H);
    int mx = 1;
    WAVaudio audio = new_audio(rate, 1);
    for (int x = 0; x < img.w; x++)
    {
        for (int i = 0; i < img.h; i++)
        {
            a[i] = polar(magnitude * (rand() % (contrast * rgb2gray(img.mat[i][x]) + black)), 0.0);
            a[H - i - 1] = a[i];
        }
        b = rfft(a);
        for (int i = 0; i < H; i++)
        {
            audio.audio[0].push_back(b[i]);
            if (b[i] > mx)
                mx = b[i];
        }
    }
    for (int i = 0; i < audio.audio[0].size(); i++)
        audio.audio[0][i] = 32000ll * audio.audio[0][i] / mx;
    return audio;
}

double f(double x)
{
    double ans = 1 + log(x) / log(155);
    return max(0.0, ans);
}

int main()
{
    // пример использования
    // генерация звука, спектрограмма которого повторяет заданную картинку
    BMPimage img1 = load_image("original.bmp"); // загружаем картинку, по которой сгенерируется звук
    WAVaudio audio1 = make_noise(img1, 44100);
    save_audio(audio1, "sound.wav");
    
    // генерация спектрограммы указанного звука
    WAVaudio audio2 = load_audio("sound.wav"); // в данном примере можно было использовать audio1, т.к. audio1 == audio2
    vector <RGB> gradient = make_gradient(255, RGB(0, 0, 0), RGB(0, 200, 200), f);
    BMPimage img2 = make_spectrogram(audio2, img1.h, gradient);
    save_image(img2, "result.bmp");
}