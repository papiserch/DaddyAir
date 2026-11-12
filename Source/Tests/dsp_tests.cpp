// BABY AIR - DSP validation suite. Standalone (no JUCE). Build: see README / Tests/CMakeLists.txt
//   ./baby_air_dsp_tests            run all tests
//   ./baby_air_dsp_tests --calibrate  bisect levelMax per mode to hit the target "excitation" (see docs/CALIBRATION.md)
#include "../DSP/ExciterEngine.h"
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

using namespace babyair;
using cd = std::complex<double>;

static int g_fail = 0;
#define CHECK(cond, ...)                                                    \
    do { if (! (cond)) { ++g_fail; std::printf ("  [FAIL] "); std::printf (__VA_ARGS__); std::printf ("\n"); } \
         else { std::printf ("  [ ok ] "); std::printf (__VA_ARGS__); std::printf ("\n"); } } while (0)

// ----------------------------------------------------------------------------- utilities
static void fft (std::vector<cd>& a)
{
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i)
    {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap (a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1)
    {
        const double ang = -2.0 * kPi / (double) len;
        const cd wl (std::cos (ang), std::sin (ang));
        for (size_t i = 0; i < n; i += len)
        {
            cd w (1);
            for (size_t j = 0; j < len / 2; ++j)
            {
                const cd u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v; a[i + j + len / 2] = u - v; w *= wl;
            }
        }
    }
}

struct Stereo { std::vector<float> l, r; };

static Stereo pinkNoise (size_t n, float rmsDb, unsigned seed)
{
    std::mt19937 rng (seed);
    std::uniform_real_distribution<float> u (-1.0f, 1.0f);
    Stereo s { std::vector<float> (n), std::vector<float> (n) };
    for (auto* ch : { &s.l, &s.r })
    {
        float b0 = 0, b1 = 0, b2 = 0, b3 = 0, b4 = 0, b5 = 0, b6 = 0; double sq = 0;
        for (size_t i = 0; i < n; ++i)
        {
            const float w = u (rng);   // Paul Kellet pink filter
            b0 = 0.99886f * b0 + w * 0.0555179f; b1 = 0.99332f * b1 + w * 0.0750759f;
            b2 = 0.96900f * b2 + w * 0.1538520f; b3 = 0.86650f * b3 + w * 0.3104856f;
            b4 = 0.55000f * b4 + w * 0.5329522f; b5 = -0.7616f * b5 - w * 0.0168980f;
            const float p = b0 + b1 + b2 + b3 + b4 + b5 + b6 + w * 0.5362f; b6 = w * 0.115926f;
            (*ch)[i] = p; sq += (double) p * p;
        }
        const float g = dbToLin (rmsDb) / (float) std::sqrt (sq / (double) n);
        for (auto& v : *ch) v *= g;
    }
    return s;
}

static Stereo sine (size_t n, double fs, double f, float db)
{
    Stereo s { std::vector<float> (n), std::vector<float> (n) };
    const float a = dbToLin (db);
    for (size_t i = 0; i < n; ++i) s.l[i] = s.r[i] = a * (float) std::sin (2.0 * kPi * f * (double) i / fs);
    return s;
}

static Stereo runEngine (const ExciterParams& p, const Stereo& in, double fs, int block = 512, int* latencyOut = nullptr, int nCh = 2)
{
    ExciterEngine e;
    e.setParams (p);
    e.prepare (fs, block, nCh);
    if (latencyOut) *latencyOut = e.getLatencySamples();
    Stereo out = in;
    const size_t n = in.l.size();
    for (size_t pos = 0; pos < n; pos += (size_t) block)
    {
        const int len = (int) std::min ((size_t) block, n - pos);
        float* ch[2] = { out.l.data() + pos, out.r.data() + pos };
        e.process (ch, nCh, len);
    }
    return out;
}

// Mean band power (Welch, Hann, 4096) in [f1,f2] (dB, arbitrary reference), skipping the first `skip` samples
static double bandPowerDb (const std::vector<float>& x, double fs, double f1, double f2, size_t skip)
{
    const size_t N = 4096, hop = 2048;
    double acc = 0; int frames = 0;
    for (size_t p = skip; p + N <= x.size(); p += hop, ++frames)
    {
        std::vector<cd> a (N);
        for (size_t i = 0; i < N; ++i) a[i] = (double) x[p + i] * (0.5 - 0.5 * std::cos (2 * kPi * (double) i / (double) N));
        fft (a);
        for (size_t k = 0; k < N / 2; ++k)
        {
            const double f = (double) k * fs / (double) N;
            if (f >= f1 && f <= f2) acc += std::norm (a[k]);
        }
    }
    return 10.0 * std::log10 (std::max (acc / std::max (frames, 1), 1e-30));
}

// coherent single-bin helper: choose f so that it falls exactly on a bin of an N-point FFT
static double coherentFreq (double fs, size_t N, double f) { return std::round (f / (fs / (double) N)) * fs / (double) N; }

static std::vector<double> spectrumMag (const std::vector<float>& x, size_t start, size_t N)
{
    std::vector<cd> a (N);
    for (size_t i = 0; i < N; ++i) a[i] = (double) x[start + i];
    fft (a);
    std::vector<double> m (N / 2);
    for (size_t k = 0; k < N / 2; ++k) m[k] = std::abs (a[k]) * 2.0 / (double) N;
    return m;
}

static double excitationDb (int mode, double fs, int osLog2, float inDb = -18.0f, double f1 = 5000, double f2 = 16000)
{
    const size_t n = (size_t) (fs * 4);
    const auto in = pinkNoise (n, inDb, 1234);
    ExciterParams p; p.mode = mode; p.midHigh = 1.0f; p.high = 1.0f; p.osLog2 = osLog2;
    const auto out = runEngine (p, in, fs);
    const size_t skip = (size_t) fs;
    return bandPowerDb (out.l, fs, f1, f2, skip) - bandPowerDb (in.l, fs, f1, f2, skip);
}


// THD (%) of a 3 kHz / -20 dBFS tone at 100%/100% (musical-safety metric used by the calibration)
static double toneThdPct (int mode)
{
    const double fs = 48000; const size_t N = 32768, warm = 48000;
    const double f = coherentFreq (fs, N, 3000.0);
    const auto in = sine (warm + N, fs, f, -20.0f);
    ExciterParams p; p.mode = mode; p.midHigh = 1; p.high = 1;
    const auto m = spectrumMag (runEngine (p, in, fs).l, warm, N);
    const size_t k0 = (size_t) std::llround (f / (fs / (double) N));
    double all = 0;
    for (size_t k = 0; k < m.size(); ++k) if (std::abs ((long) k - (long) k0) > 2 && k > 2) all += m[k] * m[k];
    return 100.0 * std::sqrt (all) / m[k0];
}

// Voice-like test signal: harmonic series on 190 Hz with a -10 dB/oct tilt above 500 Hz and slow vibrato, -20 dBFS RMS
static Stereo voiceLike (double fs, size_t n)
{
    Stereo s { std::vector<float> (n, 0.0f), std::vector<float> (n, 0.0f) };
    std::mt19937 rng (21); std::uniform_real_distribution<double> u (0, 2 * kPi);
    double ph[80]; for (auto& p : ph) p = u (rng);
    double acc = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const double t = (double) i / fs, f0 = 190.0 * (1.0 + 0.012 * std::sin (2 * kPi * 5.2 * t));
        double v = 0;
        for (int h = 1; h <= 80 && h * f0 < 0.45 * fs; ++h)
        {
            const double fr = h * f0, a = fr < 500 ? 1.0 : std::pow (fr / 500.0, -10.0 / 6.02);
            v += a * std::sin (ph[h - 1] + 2 * kPi * fr * t);
        }
        s.l[i] = (float) v; acc += v * v;
    }
    const float g = dbToLin (-20.0f) / (float) std::sqrt (acc / (double) n);
    for (size_t i = 0; i < n; ++i) { s.l[i] *= g; s.r[i] = s.l[i]; }
    return s;
}

// ----------------------------------------------------------------------------- tests
static void testOversampler()
{
    std::printf ("\n== Oversampler: latency + transparency (identity shaper) ==\n");
    const double rates[] = { 44100, 48000, 96000, 192000 };
    for (double fs : rates)
        for (int lg = 0; lg <= 4; ++lg)
        {
            Oversampler os; os.prepare (fs, lg);
            const int L = os.latencySamples();
            const size_t n = 16384;
            std::vector<float> x (n), y (n);
            const double fr[3] = { 500.0, 5000.0, std::min (15000.0, 0.34 * fs) };
            for (size_t i = 0; i < n; ++i)
            {
                double v = 0; for (double f : fr) v += 0.25 * std::sin (2 * kPi * f * (double) i / fs);
                x[i] = (float) v;
                y[i] = os.process (x[i], [] (float*, int) {});
            }
            double e = 0, s = 0;
            for (size_t i = 2000; i < n; ++i) { const double d = (double) y[i] - (double) x[i - (size_t) L]; e += d * d; s += (double) x[i] * x[i]; }
            const double errDb = 10 * std::log10 (e / s + 1e-30);
            std::printf ("  fs=%6.0f  %2dx  latency=%3d samples (%.2f ms)  error=%.1f dB\n", fs, 1 << lg, L, 1000.0 * L / fs, errDb);
            if (errDb > -60.0) { ++g_fail; std::printf ("  [FAIL] identity error too high\n"); }
        }
}

static double aliasDb (int lg, bool adaa, double fs)
{
    const size_t N = 32768, warm = 8192;
    const double f0 = coherentFreq (fs, N, 9000.0);
    Oversampler os; os.prepare (fs, lg);
    AdaaShaper sh; sh.setup (AdaaShaper::Type::Tanh, 0.0f, 0.3f, adaa);
    const float drive = 8.0f;
    std::vector<float> y (N + warm);
    for (size_t i = 0; i < N + warm; ++i)
    {
        const float x = 0.9f * (float) std::sin (2 * kPi * f0 * (double) i / fs);
        y[i] = os.process (x, [&] (float* b, int n) { for (int k = 0; k < n; ++k) b[k] = sh.process (drive * b[k]); });
    }
    const auto m = spectrumMag (y, warm, N);
    const size_t k0 = (size_t) std::llround (f0 / (fs / (double) N));
    double harm = m[k0] * m[k0], non = 0;
    for (size_t k = 1; k < m.size(); ++k) if (k % k0 != 0 && k != 0) non += m[k] * m[k];
    return 10 * std::log10 (non / harm + 1e-30);
}

static void testAliasing()
{
    std::printf ("\n== Aliasing (9 kHz sine through tanh drive x8 @44.1k; non-harmonic energy re fundamental) ==\n");
    const double fs = 44100;
    const double a1 = aliasDb (0, false, fs), a1d = aliasDb (0, true, fs), a2d = aliasDb (1, true, fs),
                 a4d = aliasDb (2, true, fs), a8 = aliasDb (3, false, fs), a16 = aliasDb (4, false, fs);
    std::printf ("  1x naive      : %7.1f dB\n  1x ADAA       : %7.1f dB\n  2x ADAA       : %7.1f dB\n  4x ADAA       : %7.1f dB\n  8x naive      : %7.1f dB\n  16x naive     : %7.1f dB\n", a1, a1d, a2d, a4d, a8, a16);
    CHECK (a1d < a1 - 6.0, "ADAA reduces aliasing vs naive at 1x (%.1f dB better)", a1 - a1d);
    CHECK (a4d < a2d, "4x ADAA cleaner than 2x ADAA");
    CHECK (a8 < a1 - 25.0, "8x oversampling >25 dB cleaner than no oversampling (%.1f dB)", a1 - a8);
    CHECK (a16 <= a8 + 3.0, "16x not worse than 8x");
}

static void testBypassAndLatency()
{
    std::printf ("\n== Dry path alignment / bypass / zero-amount ==\n");
    const double fs = 48000;
    const auto in = pinkNoise (48000, -12.0f, 7);
    for (int lg : { 1, 2, 3, 4 })
    {
        ExciterParams p; p.midHigh = 0; p.high = 0; p.osLog2 = lg;
        int L = 0; const auto out = runEngine (p, in, fs, 256, &L);
        double maxd = 0; for (size_t i = (size_t) L; i < in.l.size(); ++i) maxd = std::max (maxd, (double) std::abs (out.l[i] - in.l[i - (size_t) L]));
        CHECK (maxd < 1e-6, "amount=0%% is the latency-aligned dry signal, %2dx, L=%d, max err=%.2e", 1 << lg, L, maxd);
    }
    ExciterParams p; p.mix = 0.0f; int L = 0;
    auto out = runEngine (p, in, fs, 100, &L);
    double maxd = 0; for (size_t i = (size_t) L + 100; i < in.l.size(); ++i) maxd = std::max (maxd, (double) std::abs (out.l[i] - in.l[i - (size_t) L]));
    CHECK (maxd < 1e-6, "MIX=0%% == dry (aligned), max err=%.2e", maxd);
    p = ExciterParams(); p.power = false;
    out = runEngine (p, in, fs, 100, &L);
    maxd = 0; for (size_t i = (size_t) L + 100; i < in.l.size(); ++i) maxd = std::max (maxd, (double) std::abs (out.l[i] - in.l[i - (size_t) L]));
    CHECK (maxd < 1e-6, "Power OFF == true bypass (delayed by PDC), max err=%.2e", maxd);
}

static void testLowEndUntouched()
{
    std::printf ("\n== Frequency response of the fundamental (sine -24 dBFS, 100%%/100%%, gain in dB) ==\n");
    const double fs = 48000; const size_t N = 16384, warm = 48000;
    const double freqs[] = { 100, 500, 1000, 2000, 3000, 5000, 8000, 12000, 16000 };
    std::printf ("  %-11s", "mode\\Hz"); for (double f : freqs) std::printf ("%7.0f", f); std::printf ("\n");
    for (int mode = 0; mode < 4; ++mode)
    {
        std::printf ("  %-11s", getModeConfig (mode).name);
        double lowDev = 0;
        for (double f0 : freqs)
        {
            const double f = coherentFreq (fs, N, f0);
            const auto in = sine (warm + N, fs, f, -24.0f);
            ExciterParams p; p.mode = mode; p.midHigh = 1; p.high = 1;
            int L; const auto out = runEngine (p, in, fs, 512, &L);
            const size_t k = (size_t) std::llround (f / (fs / (double) N));
            const double g = 20 * std::log10 (spectrumMag (out.l, warm, N)[k] / spectrumMag (in.l, warm, N)[k]);
            std::printf ("%7.2f", g);
            if (f0 <= 1000) lowDev = std::max (lowDev, std::abs (g));
        }
        std::printf ("\n");
        CHECK (lowDev < 0.05, "%s: <=1 kHz untouched (max dev %.3f dB)", getModeConfig (mode).name, lowDev);
    }
}

static void testModesDiffer()
{
    std::printf ("\n== Harmonic spectra, 3 kHz sine -20 dBFS, 100%%/100%% (dB re output fundamental) ==\n");
    const double fs = 48000; const size_t N = 32768, warm = 48000;
    const double f = coherentFreq (fs, N, 3000.0);
    const auto in = sine (warm + N, fs, f, -20.0f);
    double H[4][6], thd[4];
    std::printf ("  %-11s %7s %7s %7s %7s %7s | THD    even/odd\n", "mode", "H2", "H3", "H4", "H5", "H6");
    for (int mode = 0; mode < 4; ++mode)
    {
        ExciterParams p; p.mode = mode; p.midHigh = 1; p.high = 1;
        const auto out = runEngine (p, in, fs);
        const auto m = spectrumMag (out.l, warm, N);
        const size_t k0 = (size_t) std::llround (f / (fs / (double) N));
        double ev = 0, od = 0;
        for (int h = 2; h <= 6; ++h)
        {
            const size_t k = k0 * (size_t) h;
            H[mode][h - 1] = k < m.size() ? 20 * std::log10 (m[k] / m[k0] + 1e-12) : -200;
            const double pw = k < m.size() ? (m[k] / m[k0]) * (m[k] / m[k0]) : 0;
            (h % 2 ? od : ev) += pw;
        }
        double all = 0; for (size_t k = k0 + 2; k < m.size(); ++k) if (std::abs ((long) k - (long) k0) > 2) all += m[k] * m[k];
        thd[mode] = 100.0 * std::sqrt (all) / m[k0];
        std::printf ("  %-11s %7.1f %7.1f %7.1f %7.1f %7.1f | %5.2f%%  %5.1f dB\n", getModeConfig (mode).name, H[mode][1], H[mode][2], H[mode][3], H[mode][4], H[mode][5], thd[mode], 10 * std::log10 ((ev + 1e-20) / (od + 1e-20)));
    }
    int diffPairs = 0, pairs = 0;
    for (int a = 0; a < 4; ++a) for (int b = a + 1; b < 4; ++b)
    {
        ++pairs; double d = 0; for (int h = 1; h < 6; ++h) d = std::max (d, std::abs (H[a][h] - H[b][h]));
        if (d > 2.0) ++diffPairs;
    }
    CHECK (diffPairs == pairs, "all %d mode pairs differ by >2 dB in at least one harmonic (%d/%d)", pairs, diffPairs, pairs);
    CHECK (thd[3] > thd[0], "BABY BOO THD (%.2f%%) > CLEAN THD (%.2f%%)", thd[3], thd[0]);
}

// ----------------------------------------------------------------------------- synthetic material
// NOTE: these are SYNTHETIC stand-ins (formant-filtered glottal pulses, noise bursts, ...). They are NOT recordings
// of real voices/instruments and are no substitute for listening tests.
static void normaliseRms (Stereo& s, float db)
{
    double a = 0; for (float v : s.l) a += (double) v * v;
    const float g = dbToLin (db) / (float) std::sqrt (a / (double) s.l.size());
    for (size_t i = 0; i < s.l.size(); ++i) { s.l[i] *= g; s.r[i] = s.l[i]; }
}
static Stereo formantVoice (double fs, size_t n, double f0, double F1, double F2, double F3, bool whisper, unsigned seed)
{
    std::mt19937 rng (seed); std::uniform_real_distribution<float> u (-1, 1);
    Stereo s { std::vector<float> (n, 0.0f), std::vector<float> (n, 0.0f) };
    Biquad bp[3]; const double F[3] = { F1, F2, F3 };
    for (int k = 0; k < 3; ++k)   // constant-skirt band-pass resonators
    {
        const double w = 2 * kPi * F[k] / fs, q = 8.0, al = std::sin (w) / (2 * q), a0 = 1 + al;
        bp[k].setCoeffs ({ al / a0, 0.0, -al / a0, -2 * std::cos (w) / a0, (1 - al) / a0 });
    }
    double ph = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const double t = (double) i / fs, env = 0.55 + 0.45 * std::sin (2 * kPi * 3.5 * t) * std::sin (2 * kPi * 3.5 * t);
        ph += f0 * (1 + 0.01 * std::sin (2 * kPi * 5.0 * t)) / fs; ph -= std::floor (ph);
        const float glottal = whisper ? 0.7f * u (rng) : (float) ((1.0 - 2.0 * ph) + 0.03 * u (rng)); // saw-like pulse + breath
        float y = 0; for (auto& f : bp) y += f.process (glottal);
        s.l[i] = (float) (y * env);
    }
    normaliseRms (s, -20.0f); return s;
}
static Stereo noiseBursts (double fs, size_t n, double hp, double lp, double period, double decayMs, float db, unsigned seed)
{
    std::mt19937 rng (seed); std::uniform_real_distribution<float> u (-1, 1);
    Stereo s { std::vector<float> (n, 0.0f), std::vector<float> (n, 0.0f) };
    Butter4 h, l; h.setHighpass (fs, hp); l.setLowpass (fs, lp);
    for (size_t i = 0; i < n; ++i)
    {
        const double t = std::fmod ((double) i / fs, period);
        const float env = (float) std::exp (-t * 1000.0 / decayMs);
        s.l[i] = l.process (h.process (u (rng))) * env;
    }
    normaliseRms (s, db); return s;
}
static Stereo synthPad (double fs, size_t n)
{
    Stereo s { std::vector<float> (n, 0.0f), std::vector<float> (n, 0.0f) };
    Butter4 lp; lp.setLowpass (fs, 7000);
    const double fr[4] = { 220.0, 277.18, 329.63, 440.0 };
    for (size_t i = 0; i < n; ++i)
    {
        float v = 0; for (double f : fr) for (double dt : { -0.003, 0.003 }) { const double ph = std::fmod ((double) i * f * (1 + dt) / fs, 1.0); v += (float) (2 * ph - 1); }
        s.l[i] = lp.process (v);
    }
    normaliseRms (s, -20.0f); return s;
}
static Stereo mixDown (std::initializer_list<const Stereo*> parts, std::initializer_list<float> gainsDb)
{
    const size_t n = (*parts.begin())->l.size();
    Stereo s { std::vector<float> (n, 0.0f), std::vector<float> (n, 0.0f) };
    auto g = gainsDb.begin();
    for (auto* p : parts) { const float k = dbToLin (*g++); for (size_t i = 0; i < n; ++i) s.l[i] += k * p->l[i]; }
    for (size_t i = 0; i < n; ++i) s.r[i] = s.l[i];
    return s;
}
struct Material { const char* name; Stereo sig; };
static std::vector<Material> makeMaterials (double fs)
{
    const size_t n = (size_t) (fs * 4);
    std::vector<Material> m;
    m.push_back ({ "voice male",   formantVoice (fs, n, 110, 700, 1200, 2600, false, 1) });
    m.push_back ({ "voice female", formantVoice (fs, n, 220, 850, 1400, 3000, false, 2) });
    m.push_back ({ "voice bass",   formantVoice (fs, n,  82, 600, 1000, 2400, false, 3) });
    m.push_back ({ "voice high",   formantVoice (fs, n, 330, 900, 1600, 3200, false, 4) });
    m.push_back ({ "whisper",      formantVoice (fs, n, 150, 750, 1300, 2800, true,  5) });
    auto v = formantVoice (fs, n, 130, 700, 1200, 2600, false, 6);
    auto sib = noiseBursts (fs, n, 5500, 9000, 0.4, 70, -24.0f, 7);
    m.push_back ({ "voice+sibil.", mixDown ({ &v, &sib }, { 0.0f, 0.0f }) });
    m.push_back ({ "hi-hats",      noiseBursts (fs, n, 7000, 16000, 0.25, 55, -20.0f, 8) });
    m.push_back ({ "synth pad",    synthPad (fs, n) });
    auto hat = noiseBursts (fs, n, 7000, 16000, 0.25, 55, -26.0f, 9); auto pad = synthPad (fs, n);
    auto kick = Stereo { std::vector<float> (n), std::vector<float> (n) };
    for (size_t i = 0; i < n; ++i) { const double t = std::fmod ((double) i / fs, 0.5); kick.l[i] = kick.r[i] = (float) (std::sin (2 * kPi * 55.0 * t) * std::exp (-t * 7.0)); }
    m.push_back ({ "full mix",     mixDown ({ &v, &hat, &pad, &kick }, { -2.0f, -2.0f, -10.0f, -4.0f }) });
    return m;
}

// ----------------------------------------------------------------------------- tests (intensity / identity)
static double bandDelta (const Stereo& o, const Stereo& in, double fs, double a, double b)
{
    return bandPowerDb (o.l, fs, a, b, (size_t) fs) - bandPowerDb (in.l, fs, a, b, (size_t) fs);
}
static double rmsDb (const std::vector<float>& x, size_t skip)
{
    double a = 0; for (size_t i = skip; i < x.size(); ++i) a += (double) x[i] * x[i];
    return 10 * std::log10 (a / (double) (x.size() - skip) + 1e-30);
}

static void testIntensity()
{
    std::printf ("\n== Intensity law, voice-like signal (-20 dBFS), both knobs together. dB change in the AIR band 9-20 kHz and PRESENCE band 2.5-9 kHz ==\n");
    const double fs = 48000; const auto voice = voiceLike (fs, (size_t) fs * 4);
    std::printf ("  %-11s | %-34s | %-34s\n", "mode", "AIR 9-20k @25/50/75/100 %", "PRESENCE 2.5-9k @25/50/75/100 %");
    for (int mode = 0; mode < 4; ++mode)
    {
        double air[4], pres[4]; int k = 0;
        for (float a : { 0.25f, 0.5f, 0.75f, 1.0f })
        {
            ExciterParams p; p.mode = mode; p.midHigh = a; p.high = a;
            const auto o = runEngine (p, voice, fs);
            air[k] = bandDelta (o, voice, fs, 9000, 20000); pres[k] = bandDelta (o, voice, fs, 2500, 9000); ++k;
        }
        std::printf ("  %-11s | %+6.1f %+6.1f %+6.1f %+6.1f dB        | %+6.1f %+6.1f %+6.1f %+6.1f dB\n", getModeConfig (mode).name, air[0], air[1], air[2], air[3], pres[0], pres[1], pres[2], pres[3]);
        const double* m = (mode == 2) ? pres : air;    // VINTAGE is a warmth/density mode: judged in the presence band
        bool steps = true; for (int i = 1; i < 4; ++i) steps &= (m[i] - m[i - 1] >= 1.5);
        CHECK (steps && m[3] >= 6.0, "%s: 25/50/75/100%% are >=1.5 dB apart and 100%% reaches >=+6 dB (%+.1f dB)", getModeConfig (mode).name, m[3]);
    }
}

static void testIndependence()
{
    std::printf ("\n== MID HIGH vs HIGH independence (voice-like, 100%% on one knob, 0%% on the other; delta dB) ==\n");
    const double fs = 48000; const auto voice = voiceLike (fs, (size_t) fs * 4);
    std::printf ("  %-11s | MID only: 3-8k / 13-20k | HIGH only: 3-8k / 13-20k\n", "mode");
    for (int mode = 0; mode < 4; ++mode)
    {
        ExciterParams pm; pm.mode = mode; pm.midHigh = 1; pm.high = 0;
        ExciterParams ph = pm; ph.midHigh = 0; ph.high = 1;
        const auto om = runEngine (pm, voice, fs), oh = runEngine (ph, voice, fs);
        const double m1 = bandDelta (om, voice, fs, 3000, 8000), m2 = bandDelta (om, voice, fs, 13000, 20000);
        const double h1 = bandDelta (oh, voice, fs, 3000, 8000), h2 = bandDelta (oh, voice, fs, 13000, 20000);
        std::printf ("  %-11s | %+6.1f / %+6.1f        | %+6.1f / %+6.1f\n", getModeConfig (mode).name, m1, m2, h1, h2);
        if (mode != 2)
            CHECK (m1 > m2 + 6.0 && h2 > h1 + 6.0, "%s: MID HIGH acts in the presence region, HIGH acts in the air region", getModeConfig (mode).name);
        else
            CHECK (m1 > 1.0 && h2 >= m2, "%s: MID HIGH acts in presence, HIGH reaches higher", getModeConfig (mode).name);
    }
}

static void testMaterials()
{
    std::printf ("\n== Synthetic material, 100%%/100%%: AIR 9-20k delta | RMS delta | peak delta | TRIM that would level-match ==\n");
    std::printf ("  (synthetic stand-ins, NOT real recordings; 'level-match' = the TRIM value that cancels the RMS change)\n");
    const double fs = 48000; const auto mats = makeMaterials (fs);
    bool allBounded = true;
    for (const auto& m : mats)
    {
        std::printf ("  %-13s", m.name);
        for (int mode = 0; mode < 4; ++mode)
        {
            double air[4]; double rmsD = 0, pkD = 0; int k = 0;
            for (float a : { 0.25f, 0.5f, 0.75f, 1.0f })
            {
                ExciterParams p; p.mode = mode; p.midHigh = a; p.high = a;
                const auto o = runEngine (p, m.sig, fs);
                air[k++] = bandDelta (o, m.sig, fs, 9000, 20000);
                if (a == 1.0f)
                {
                    rmsD = rmsDb (o.l, (size_t) fs) - rmsDb (m.sig.l, (size_t) fs);
                    float po = 0, pi = 0; for (size_t i = (size_t) fs; i < o.l.size(); ++i) { po = std::max (po, std::abs (o.l[i])); pi = std::max (pi, std::abs (m.sig.l[i])); }
                    pkD = linToDb (po) - linToDb (pi);
                    bool fin = true; for (float v : o.l) fin &= std::isfinite (v);
                    allBounded &= fin && po <= std::max (pi, 0.98f) + 0.02f;
                }
            }
            std::printf ("\n      %-10s air %+5.1f %+5.1f %+5.1f %+5.1f | RMS %+4.1f dB | peak %+4.1f dB | trim %+4.1f dB", getModeConfig (mode).name, air[0], air[1], air[2], air[3], rmsD, pkD, -rmsD);
        }
        std::printf ("\n");
    }
    CHECK (allBounded, "all materials/modes: finite, and peaks never exceed max(input peak, 0.98) thanks to harmonic protection");
}

static void testSibilance()
{
    std::printf ("\n== Sibilance protection (voice + sibilant noise bursts 5.5-9 kHz), BABY BOO 100%% ==\n");
    const double fs = 48000; const size_t n = (size_t) fs * 4;
    auto v = formantVoice (fs, n, 130, 700, 1200, 2600, false, 6), sib = noiseBursts (fs, n, 5500, 9000, 0.4, 70, -24.0f, 7);
    const auto sig = mixDown ({ &v, &sib }, { 0.0f, 0.0f });
    ExciterParams p; p.mode = 3; p.midHigh = 1; p.high = 1;
    const auto o = runEngine (p, sig, fs);
    // sibilance "burst peaks" in 9-20k vs average 9-20k: protection should keep the burst peak from jumping out
    auto burstRatio = [&] (const std::vector<float>& x)
    {   // crest of 9-20 kHz band envelope over 20 ms windows
        Butter4 h; h.setHighpass (fs, 9000); double mx = 0, mean = 0; int cnt = 0;
        const size_t w = (size_t) (fs * 0.02); double acc = 0; size_t k = 0;
        for (size_t i = (size_t) fs; i < x.size(); ++i) { const float y = h.process (x[i]); acc += (double) y * y; if (++k == w) { mx = std::max (mx, acc / w); mean += acc / w; ++cnt; acc = 0; k = 0; } }
        return 10 * std::log10 (mx / (mean / cnt));
    };
    const double c0 = burstRatio (sig.l), c1 = burstRatio (o.l);
    std::printf ("  air-band crest (peak-window / mean-window): dry %.1f dB  ->  processed %.1f dB\n", c0, c1);
    CHECK (c1 <= c0 + 3.0, "sibilant bursts are not exaggerated beyond +3 dB relative to the average air level");
}

static void testIMDAliasingModes()
{
    std::printf ("\n== IMD and aliasing per mode (100%%/100%%, 8x) ==\n");
    const double fs = 48000; const size_t N = 32768, warm = 48000;
    const double f1 = coherentFreq (fs, N, 3000.0), f2 = coherentFreq (fs, N, 3600.0);
    Stereo two = sine (warm + N, fs, f1, -26.0f);
    const auto s2 = sine (warm + N, fs, f2, -26.0f);
    for (size_t i = 0; i < two.l.size(); ++i) two.l[i] = two.r[i] = two.l[i] + s2.l[i];
    std::printf ("  %-11s | IMD 3rd-order products (2f1-f2, 2f2-f1) re tones | audible aliasing (6.1k sine, <20 kHz non-harmonic energy re fund.) 2x / 8x\n", "mode");
    for (int mode = 0; mode < 4; ++mode)
    {
        ExciterParams p; p.mode = mode; p.midHigh = 1; p.high = 1;
        const auto m = spectrumMag (runEngine (p, two, fs).l, warm, N);
        auto bin = [&] (double f) { return (size_t) std::llround (f / (fs / (double) N)); };
        const double tone = std::max (m[bin (f1)], m[bin (f2)]), im = std::max (m[bin (2 * f1 - f2)], m[bin (2 * f2 - f1)]);
        double imLow[2]; int q = 0;
        for (float a : { 0.25f, 0.5f })
        {
            ExciterParams pa = p; pa.midHigh = a; pa.high = a;
            const auto ma = spectrumMag (runEngine (pa, two, fs).l, warm, N);
            imLow[q++] = 20 * std::log10 (std::max (ma[bin (2 * f1 - f2)], ma[bin (2 * f2 - f1)]) / std::max (ma[bin (f1)], ma[bin (f2)]) + 1e-12);
        }
        double al[2]; int k = 0;
        for (int lg : { 1, 3 })
        {
            const double f0 = coherentFreq (fs, N, 6100.0);
            const auto in = sine (warm + N, fs, f0, -18.0f);
            ExciterParams q = p; q.osLog2 = lg;
            const auto mm = spectrumMag (runEngine (q, in, fs).l, warm, N);
            const size_t k0 = bin (f0); double harm = mm[k0] * mm[k0], non = 0;
            for (size_t b = 3; b < bin (20000.0); ++b) if (b % k0 != 0) non += mm[b] * mm[b];   // audible band only (<20 kHz): content at 22-24 kHz is the legitimate 4th harmonic reflected by the half-band decimator
            al[k++] = 10 * std::log10 (non / harm + 1e-30);
        }
        std::printf ("  %-11s | IMD 25%%: %6.1f  50%%: %6.1f  100%%: %6.1f dB | %7.1f dB / %7.1f dB\n", getModeConfig (mode).name, imLow[0], imLow[1], 20 * std::log10 (im / tone + 1e-12), al[0], al[1]);
        CHECK (al[1] <= al[0] + 0.5, "%s: 8x is not worse than 2x for aliasing", getModeConfig (mode).name);
    }
}

static void testProtection()
{
    std::printf ("\n== Headroom protection ==\n");
    const double fs = 48000; const auto in = pinkNoise ((size_t) fs * 3, -3.0f, 77);
    float pin = 0; for (float v : in.l) pin = std::max (pin, std::abs (v));
    for (int mode = 0; mode < 4; ++mode)
    {
        ExciterParams p; p.mode = mode; p.midHigh = 1; p.high = 1;
        const auto o = runEngine (p, in, fs);
        float po = 0; for (size_t i = (size_t) fs; i < o.l.size(); ++i) po = std::max (po, std::abs (o.l[i]));
        std::printf ("  %-11s input peak %.3f -> output peak %.3f\n", getModeConfig (mode).name, pin, po);
        CHECK (po <= std::max (pin, 0.98f) + 0.02f, "%s: the exciter never pushes the peak beyond max(input peak, 0.98)", getModeConfig (mode).name);
    }
    ExciterParams p; p.mode = 3; p.midHigh = 0; p.high = 0;
    const auto o = runEngine (p, in, fs); double d = 0;
    for (size_t i = 200; i < o.l.size(); ++i) d = std::max (d, (double) std::abs (o.l[i] - in.l[i - 47]));
    CHECK (d < 1e-6, "with 0%% the loud dry signal passes untouched (max err %.1e)", d);
}


static void testSampleRatesAndRobustness()
{
    std::printf ("\n== Sample rates, extreme input, mono, block sizes ==\n");
    for (double fs : { 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 })
    {
        const auto in = pinkNoise ((size_t) (fs * 2), -18.0f, 99);
        ExciterParams p; p.mode = 3; p.midHigh = 1; p.high = 1;
        int L; const auto out = runEngine (p, in, fs, 512, &L);
        bool fin = true; float pk = 0; for (float v : out.l) { fin &= std::isfinite (v); pk = std::max (pk, std::abs (v)); }
        std::printf ("  fs=%6.0f  latency=%3d (%.2f ms)  excitation(5-16k)=+%4.1f dB  peak=%.3f\n", fs, L, 1000.0 * L / fs, excitationDb (3, fs, 3), pk);
        CHECK (fin && pk < 2.0f, "finite and bounded at %.0f Hz", fs);
    }
    {   // full-scale white noise + clipped square + NaN/Inf injection
        const double fs = 48000; std::mt19937 rng (5); std::uniform_real_distribution<float> u (-1, 1);
        Stereo in { std::vector<float> (96000), std::vector<float> (96000) };
        for (size_t i = 0; i < in.l.size(); ++i) { in.l[i] = u (rng); in.r[i] = ((i / 40) & 1) ? 1.0f : -1.0f; }
        in.l[5000] = NAN; in.r[6000] = INFINITY;
        for (int mode = 0; mode < 4; ++mode)
        {
            ExciterParams p; p.mode = mode; p.midHigh = 1; p.high = 1; p.trimDb = 18;
            const auto out = runEngine (p, in, fs);
            bool fin = true; float pk = 0; for (float v : out.l) { fin &= std::isfinite (v); pk = std::max (pk, std::abs (v)); }
            for (float v : out.r) fin &= std::isfinite (v);
            CHECK (fin && pk <= 16.0f, "%s: finite & bounded with 0 dBFS noise/square, NaN/Inf and +18 dB trim (peak %.2f)", getModeConfig (mode).name, pk);
        }
    }
    {   // block-size independence
        const auto in = pinkNoise (48000, -14.0f, 3);
        ExciterParams p; p.mode = 3; p.midHigh = 0.7f; p.high = 0.4f;
        const auto a = runEngine (p, in, 48000, 1), b = runEngine (p, in, 48000, 64), c = runEngine (p, in, 48000, 513);
        CHECK (a.l == b.l && b.l == c.l && a.r == c.r, "bit-identical output for block sizes 1 / 64 / 513");
    }
    {   // mono
        const auto in = pinkNoise (48000, -18.0f, 3);
        ExciterParams p; p.mode = 3; p.midHigh = 1; p.high = 1;
        const auto out = runEngine (p, in, 48000, 512, nullptr, 1);
        bool fin = true; for (float v : out.l) fin &= std::isfinite (v);
        CHECK (fin, "mono processing is finite");
    }
    {   // mono compatibility: identical L/R -> identical outputs; and L/R = (M,S) -> sum stays sane
        const auto m = pinkNoise (48000, -18.0f, 3);
        Stereo in { m.l, m.l };
        ExciterParams p; p.mode = 3; p.midHigh = 1; p.high = 1;
        const auto out = runEngine (p, in, 48000);
        CHECK (out.l == out.r, "identical L/R inputs give identical L/R outputs (centre stays centred)");
    }
}

static void testTransitions()
{
    std::printf ("\n== Click tests (mode crossfade, power toggle, parameter jumps) ==\n");
    const double fs = 48000; const size_t n = 96000;
    const auto in = sine (n, fs, 4000.0, -12.0f);
    auto maxJump = [] (const std::vector<float>& y, size_t a, size_t b) { double m = 0; for (size_t i = a + 1; i < b; ++i) m = std::max (m, (double) std::abs (y[i] - y[i - 1])); return m; };

    ExciterEngine e; ExciterParams p; p.mode = 0; p.midHigh = 1; p.high = 1; e.setParams (p); e.prepare (fs, 256, 2);
    Stereo out = in;
    for (size_t pos = 0; pos < n; pos += 256)
    {
        if (pos == 24576) { p.mode = 3; e.setParams (p); }
        if (pos == 49152) { p.mode = 2; p.midHigh = 0.1f; p.high = 1.0f; e.setParams (p); }
        if (pos == 73728) { p.power = false; e.setParams (p); }
        float* ch[2] = { out.l.data() + pos, out.r.data() + pos };
        e.process (ch, 2, (int) std::min<size_t> (256, n - pos));
    }
    // A mode change legitimately changes the steady-state slope (BABY BOO is much brighter than CLEAN), so the
    // reference is the louder-slope of the steady state BEFORE and AFTER each transition.
    const double r1 = std::max (maxJump (out.l, 12000, 24000), maxJump (out.l, 38000, 48000));
    const double r2 = std::max (maxJump (out.l, 38000, 48000), maxJump (out.l, 62000, 72000));
    const double r3 = std::max (maxJump (out.l, 62000, 72000), maxJump (out.l, 84000, 96000));
    std::printf ("  steady-state max |dy|: %.4f / %.4f / %.4f\n", r1, r2, r3);
    CHECK (maxJump (out.l, 24500, 36000) < r1 * 1.25, "CLEAN -> BABY BOO crossfade is click-free (max |dy| %.4f vs ref %.4f)", maxJump (out.l, 24500, 36000), r1);
    CHECK (maxJump (out.l, 49500, 61000) < r2 * 1.25, "BABY BOO -> VINTAGE + knob jump is click-free (max |dy| %.4f vs ref %.4f)", maxJump (out.l, 49500, 61000), r2);
    CHECK (maxJump (out.l, 73700, 82000) < r3 * 1.25, "Power OFF transition is click-free (max |dy| %.4f vs ref %.4f)", maxJump (out.l, 73700, 82000), r3);
}

static void testCpu()
{
    std::printf ("\n== CPU (single core, -O2, stereo, 44.1 kHz, realtime %% = process time / audio time) ==\n");
    const double fs = 44100; const auto in = pinkNoise ((size_t) (fs * 5), -18.0f, 11);
    for (int mode : { 0, 3 })
        for (int lg = 1; lg <= 4; ++lg)
        {
            ExciterParams p; p.mode = mode; p.midHigh = 0.8f; p.high = 0.8f; p.osLog2 = lg;
            const auto t0 = std::chrono::steady_clock::now();
            const auto out = runEngine (p, in, fs);
            const double sec = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
            std::printf ("  %-11s %2dx : %5.1f %% of one core  (%s)\n", getModeConfig (mode).name, 1 << lg, 100.0 * sec / 5.0, out.l[1000] == out.l[1000] ? "ok" : "nan");
        }
}

// ----------------------------------------------------------------------------- calibration
// Calibration target (100%/100%, voice-like signal -20 dBFS, 48 kHz, 8x). Because the added signal is mostly
// harmonics, the metric is the band where the original has almost no energy (see docs/CALIBRATION.md).
struct CalTarget { double lo, hi, dB; };
static const CalTarget kCal[4] = { { 9000, 20000, 12.0 },    // CLEAN   : +12 dB air
                                   { 9000, 20000, 18.0 },    // BRIGHT  : +18 dB air
                                   { 2500, 9000,   8.0 },    // VINTAGE : +8 dB presence/density (warmth mode)
                                   { 9000, 20000, 24.0 } };  // BABY BOO: +24 dB air
static double calMetric (int mode, const Stereo& voice)
{
    ExciterParams p; p.mode = mode; p.midHigh = 1; p.high = 1;
    return bandDelta (runEngine (p, voice, 48000), voice, 48000, kCal[mode].lo, kCal[mode].hi);
}
static void calibrate()
{
    std::printf ("Calibrating levelMax by bisection (voice-like, 100%%/100%%) ...\n");
    const auto voice = voiceLike (48000, 48000 * 4);
    for (int mode = 0; mode < 4; ++mode)
    {
        double lo = 0.05, hi = 400.0;
        for (int it = 0; it < 16; ++it)
        {
            const double mid = std::sqrt (lo * hi);
            for (auto& b : modeConfigForCalibration (mode).band) b.levelMax = (float) mid;
            (calMetric (mode, voice) < kCal[mode].dB ? lo : hi) = mid;
        }
        const double l = std::sqrt (lo * hi);
        for (auto& b : modeConfigForCalibration (mode).band) b.levelMax = (float) l;
        std::printf ("  %-11s levelMax = %.3f -> %+.1f dB (band %.0f-%.0f Hz, target %+.0f dB)  tone THD %.0f%%\n", getModeConfig (mode).name, l,
                     calMetric (mode, voice), kCal[mode].lo, kCal[mode].hi, kCal[mode].dB, toneThdPct (mode));
    }
}


int main (int argc, char** argv)
{
    if (argc > 1 && ! std::strcmp (argv[1], "--calibrate")) { calibrate(); return 0; }
    testOversampler();
    testAliasing();
    testBypassAndLatency();
    testLowEndUntouched();
    testModesDiffer();
    testIntensity();
    testIndependence();
    testMaterials();
    testSibilance();
    testIMDAliasingModes();
    testProtection();
    testSampleRatesAndRobustness();
    testTransitions();
    testCpu();
    std::printf ("\n%s (%d failure%s)\n", g_fail ? "SOME TESTS FAILED" : "ALL TESTS PASSED", g_fail, g_fail == 1 ? "" : "s");
    return g_fail ? 1 : 0;
}
