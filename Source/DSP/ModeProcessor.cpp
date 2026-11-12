#include "ModeProcessor.h"

namespace babyair
{
using T = AdaaShaper::Type;

// ---------------------------------------------------------------------------------------------
// Mode tables (v1.1 redesign: extreme-intensity DSP)
//
// Source bands sit LOWER than the audible target region so that the generated harmonics land in it:
//   MID HIGH : source 1.8-5.5 kHz  -> harmonics 3-11 kHz (presence)           [outLp keeps it out of the air band]
//   HIGH     : source 4.5-11 kHz   -> harmonics 9-20 kHz (air)
// (v1.0 used 7-18 kHz as the HIGH source: its harmonics fell above Nyquist, which made the knob nearly dead.)
//
//   CLEAN  "PURE AIR"      Chebyshev-dominant (clean, exactly-controlled 2nd/3rd) + a touch of tanh glue
//   BRIGHT "CRYSTAL SHINE" odd-dominant tanh + Chebyshev 3rd/4th, strong pre-emphasis, transient emphasis
//   VINTAGE "ANALOG GLOW"  asymmetric Vintage curve + odd tanh + soft rectifier (even), LP on harmonics, leveling
//   BABY BOO "DIGITAL AIR" NO analog saturator: Chebyshev T2/T3/T4 synthesis + soft rectifier, near level-
//                          independent drive (rho .9), heavy pre-emphasis, strong transient emphasis
// levelMax values come from `dsp_tests --calibrate` (air-band target, see docs/CALIBRATION.md)
// ---------------------------------------------------------------------------------------------
static ModeConfig kModes[(int) Mode::Count] =
{
    { "CLEAN",
      { //            inHp  inLp   outHp  outLp   emphDb emphHz dMin  dMax  rho   level  curve1    bias   w1    w2    w3    shelfDb shelfHz c2    c3    c4    trans
        BandConfig { 1800, 5500,  3000,  9000,   4.0f,  3000,  0.8f, 4.0f, 0.60f, 2.00f, T::Tanh,    0.05f, 0.35f, 0.0f, 0.0f, 6.0f,  5000,  0.90f, 0.40f, 0.0f, 0.3f },
        BandConfig { 4500, 11000, 9000,  0,      4.0f,  7000,  0.8f, 4.0f, 0.60f, 1.05f, T::Tanh,    0.05f, 0.35f, 0.0f, 0.0f, 6.0f,  12000, 0.90f, 0.40f, 0.0f, 0.3f } },
      1.0f },
    { "BRIGHT",
      {
        BandConfig { 1800, 5500,  3000,  9000,   7.0f,  3000,  1.0f, 7.0f, 0.70f, 1.79f, T::Tanh,    0.03f, 0.5f, 0.9f, 0.15f, 6.0f, 5000,  0.40f, 0.80f, 0.30f, 0.6f },
        BandConfig { 4500, 11000, 9000,  0,      7.0f,  7000,  1.0f, 7.0f, 0.70f, 1.79f, T::Tanh,    0.03f, 0.5f, 0.9f, 0.15f, 6.0f, 12000, 0.40f, 0.80f, 0.30f, 0.6f } },
      1.5f },
    { "VINTAGE",
      {
        BandConfig { 1800, 5500,  2400,  10000,  0.0f,  3000,  1.0f, 7.0f, 0.35f, 8.22f, T::Vintage, 0.00f, 1.0f, 0.35f, 0.30f, 0.0f, 6000, 0.0f, 0.0f, 0.0f, 0.0f },
        BandConfig { 4500, 11000, 7500,  15000,  0.0f,  7000,  1.0f, 7.0f, 0.35f, 8.22f, T::Vintage, 0.00f, 1.0f, 0.35f, 0.30f, 0.0f, 12000, 0.0f, 0.0f, 0.0f, 0.0f } },
      1.0f },
    { "BABY BOO",
      {
        BandConfig { 1800, 5500,  3500,  10000,  10.0f, 3000,  1.0f, 6.0f, 0.90f, 0.97f, T::Tanh,    0.00f, 0.0f, 0.0f, 0.25f, 9.0f, 5000,  1.00f, 0.90f, 0.60f, 1.0f },
        BandConfig { 4500, 11000, 9000,  0,      10.0f, 7000,  1.0f, 6.0f, 0.90f, 0.97f, T::Tanh,    0.00f, 0.0f, 0.0f, 0.25f, 9.0f, 12000, 1.00f, 0.90f, 0.60f, 1.0f } },
      2.0f }
};

const ModeConfig& getModeConfig (int mode) noexcept
{
    return kModes[std::clamp (mode, 0, (int) Mode::Count - 1)];
}

ModeConfig& modeConfigForCalibration (int mode) noexcept
{
    return kModes[std::clamp (mode, 0, (int) Mode::Count - 1)];
}

// ---------------------------------------------------------------------------------------------
// BandProcessor
// ---------------------------------------------------------------------------------------------
static constexpr float kRefLevel = 0.08f;   // target band RMS the leveler steers towards (~ -22 dBFS)
static constexpr float kGMin = 0.05f, kGMax = 30.0f;

void BandProcessor::prepare (double fs, int log2Os, const BandConfig& cfg)
{
    cfg_ = cfg;
    DescribingFunction::init();
    const bool useAdaa = log2Os <= 2; // ADAA for 1x/2x/4x; 8x/16x rely on oversampling alone (see README)
    os_.prepare (fs, log2Os);
    hp_.setHighpass (fs, cfg.inHpHz);
    lp_.setLowpass (fs, cfg.inLpHz);
    outHp_.setHighpass (fs, cfg.outHpHz);
    outLpOn_ = cfg.outLpHz > 0.0f;
    if (outLpOn_) outLp_.setLowpass (fs, cfg.outLpHz);
    oShelfOn_ = std::abs (cfg.oShelfDb) > 0.01f;
    if (oShelfOn_) oShelf_.setCoeffs (biquad::highShelf (fs, cfg.oShelfHz, cfg.oShelfDb));
    emphOn_ = std::abs (cfg.emphDb) > 0.01f;
    if (emphOn_) emph_.setCoeffs (biquad::highShelf (fs, cfg.emphHz, cfg.emphDb));
    envFast_.prepare (fs, 0.4, 6.0);
    env_.prepare (fs, 8.0, 8.0);   // symmetric => exact RMS of a sine at any frequency (needed for the fundamental correction)
    sh1_.setup (cfg.curve1, cfg.bias1, 0.3f, useAdaa);
    sh2_.setup (T::Tanh, 0.0f, 0.3f, useAdaa);
    sh3_.setup (T::Rectifier, 0.0f, 0.35f, useAdaa);
    gLen_ = std::min (os_.latencySamples() + 1, kMaxGDelay);
    logRef_ = std::log (kRefLevel);
    reset();
}

void BandProcessor::reset() noexcept
{
    hp_.reset(); lp_.reset(); outHp_.reset(); outLp_.reset(); emph_.reset(); oShelf_.reset(); dc_.reset();
    os_.reset();
    sh1_.reset(); sh2_.reset(); sh3_.reset();
    env_.reset (kRefLevel); envFast_.reset (kRefLevel);
    std::fill (std::begin (gDelay_), std::end (gDelay_), 1.0f);
    gPos_ = 0;
    idle_ = true;
}

float BandProcessor::process (float x, float amt, float sibGain) noexcept
{
    if (amt < 1.0e-5f) { idle_ = true; return 0.0f; }        // 0% = no excitation at all
    if (idle_) { reset(); idle_ = false; }

    float b = lp_.process (hp_.process (x));
    if (emphOn_) b = emph_.process (b);

    const float e = env_.process (b) + 1.0e-6f;
    const float g = std::clamp (std::exp (cfg_.rho * (logRef_ - std::log (e))), kGMin, kGMax);

    // Delay line of g, aligned with the oversampler latency, to undo the leveler exactly in time.
    const float gD = gDelay_[gPos_];
    gDelay_[gPos_] = g;
    if (++gPos_ >= gLen_) gPos_ = 0;

    const float d1 = cfg_.driveMin + (cfg_.driveMax - cfg_.driveMin) * amt;
    const float d2 = d1 * 1.6f, d3 = d1 * 0.8f;
    const float i1 = 1.0f / d1, i2 = 1.0f / d2, i3 = 1.0f / d3;

    // exact describing-function gain of tanh at the expected peak level => removes the FUNDAMENTAL from the
    // residual, so what is added is harmonic/intermodulation content only (not a hidden EQ boost or cut)
    const float A = 1.4142f * e * g;
    const float k1 = DescribingFunction::tanhGain (d1 * A), k2 = DescribingFunction::tanhGain (d2 * A);
    const float w1 = cfg_.w1, w2 = cfg_.w2, w3 = cfg_.w3;

    // Chebyshev synthesiser: the band signal is normalised by its own envelope (so the harmonic content does not
    // depend on input level), softly bounded to [-1,1], and T2/T3/T4 turn it into exact 2nd/3rd/4th harmonics.
    // "Density": the 2nd harmonic is always present, the 3rd grows linearly and the 4th quadratically with amount.
    constexpr float kN = 1.6f;
    const float cNorm = kN / (e * g);                 // u (=b*g, upsampled) * cNorm  ==  b * kN / e
    const float cOut = 1.4142f * e * g;               // restores the source-band amplitude, expressed in the g-domain
    const float k2c = cfg_.c2, k3c = cfg_.c3 * amt, k4c = cfg_.c4 * amt * amt;
    const bool useCheb = (cfg_.c2 + cfg_.c3 + cfg_.c4) > 0.0f;

    float y = os_.process (b * g, [&] (float* buf, int n) noexcept
    {
        for (int i = 0; i < n; ++i)
        {
            const float u = buf[i];
            float acc = 0.0f;
            if (w1 > 0.0f) acc += w1 * (sh1_.process (d1 * u) * i1 - k1 * u);
            if (w2 > 0.0f) acc += w2 * (sh2_.process (d2 * u) * i2 - k2 * u);
            if (w3 > 0.0f) acc += w3 * (sh3_.process (d3 * u) * i3);
            if (useCheb)
            {
                const float v = u * cNorm;
                const float x = v / std::sqrt (1.0f + v * v);          // soft bound, near-linear for small v
                const float x2 = x * x;
                const float t2 = 2.0f * x2 - 1.0f, t3 = x * (4.0f * x2 - 3.0f), t4 = 8.0f * x2 * (x2 - 1.0f) + 1.0f;
                acc += cOut * (k2c * t2 + k3c * t3 + k4c * t4);
            }
            buf[i] = acc;
        }
    });

    // transient emphasis of the added harmonics (fast / slow envelope ratio)
    const float eF = envFast_.process (b) + 1.0e-6f;
    y *= 1.0f + cfg_.transAmt * std::clamp (eF / e - 1.0f, 0.0f, 2.5f);

    const float level = cfg_.levelMax * amt;   // linear in amount: 25/50/75/100 % are ~6/3.5/2.5 dB apart where harmonics dominate
    y *= level * sibGain / gD;
    y = outHp_.process (y);
    if (oShelfOn_) y = oShelf_.process (y);
    if (outLpOn_) y = outLp_.process (y);
    y = dc_.process (y);
    return sanitize (y);
}

// ---------------------------------------------------------------------------------------------
// ModeVoice
// ---------------------------------------------------------------------------------------------
void ModeVoice::prepare (double sampleRate, int log2Os)
{
    fs_ = sampleRate;
    log2Os_ = log2Os;
    for (int c = 0; c < kMaxChannels; ++c)
    {
        sibHp_[c].setHighpass (fs_, 5500.0);
        sibLp_[c].setLowpass (fs_, 9000.0);
    }
    sibEnv_.prepare (fs_, 0.3, 30.0);
    fullEnv_.prepare (fs_, 10.0, 80.0);
    sibAtk_ = 1.0f - std::exp (-1.0f / (float) (0.001 * 1.0 * fs_));    // 1 ms attack
    sibRel_ = 1.0f - std::exp (-1.0f / (float) (0.001 * 40.0 * fs_));   // 40 ms release
    setMode (mode_);
}

void ModeVoice::setMode (int mode)
{
    mode_ = std::clamp (mode, 0, (int) Mode::Count - 1);
    const ModeConfig& m = getModeConfig (mode_);
    for (int b = 0; b < 2; ++b)
        for (int c = 0; c < kMaxChannels; ++c)
            band_[b][c].prepare (fs_, log2Os_, m.band[b]);
    reset();
}

void ModeVoice::reset() noexcept
{
    for (int b = 0; b < 2; ++b)
        for (int c = 0; c < kMaxChannels; ++c)
            band_[b][c].reset();
    for (int c = 0; c < kMaxChannels; ++c) { sibHp_[c].reset(); sibLp_[c].reset(); }
    sibEnv_.reset(); fullEnv_.reset();
    sibGain_ = 1.0f;
}

void ModeVoice::process (const float* in, int nCh, float amtMid, float amtHigh, float* out) noexcept
{
    // --- stereo-linked sibilance detector: ratio of 5.5-9 kHz energy to broadband energy ---
    float pSib = 0.0f, pFull = 0.0f;
    for (int c = 0; c < nCh; ++c)
    {
        const float s = sibLp_[c].process (sibHp_[c].process (in[c]));
        pSib += s * s;
        pFull += in[c] * in[c];
    }
    const float eSib = sibEnv_.processPower (pSib / (float) nCh);
    const float eFull = fullEnv_.processPower (pFull / (float) nCh);
    const float ratio = eSib / (eFull + 1.0e-5f);
    const float amount = getModeConfig (mode_).sibAmount;
    const float over = std::max (0.0f, ratio - 0.35f);
    const float target = 1.0f / (1.0f + amount * 8.0f * over);   // amount > 1 = stronger protection
    sibGain_ += (target < sibGain_ ? sibAtk_ : sibRel_) * (target - sibGain_);

    for (int c = 0; c < nCh; ++c)
        out[c] = band_[0][c].process (in[c], amtMid, sibGain_) + band_[1][c].process (in[c], amtHigh, sibGain_);
}
} // namespace babyair
