// BABY AIR - factory presets (exposed to the host as programs)
#pragma once

namespace babyair
{
struct FactoryPreset
{
    const char* name;
    int mode;           // 0 Clean, 1 Bright, 2 Vintage, 3 Baby Boo
    float midHigh, high; // %
    bool link;
    float trimDb, mix;  // dB, %
};

// v1.1: the DSP is far more intense than v1.0 (100% is "extreme" by design), so presets sit in the useful
// 10-50% zone and carry a small negative TRIM to approximate level-matched A/B.
inline constexpr FactoryPreset kFactoryPresets[] =
{
    { "Clean Vocal",          0, 30.0f, 40.0f, false, -1.0f, 100.0f },
    { "Bright Vocal",         1, 25.0f, 35.0f, false, -1.5f, 100.0f },
    { "Trap Vocal",           3, 22.0f, 35.0f, false, -2.0f, 100.0f },
    { "Reggaeton Vocal",      3, 28.0f, 25.0f, false, -2.0f,  90.0f },
    { "Airy Adlibs",          3, 15.0f, 45.0f, false, -2.0f, 100.0f },
    { "Modern Hi-Hats",       1, 15.0f, 50.0f, false, -1.0f,  85.0f },
    { "Vintage Air",          2, 45.0f, 30.0f, false, -1.5f, 100.0f },
    { "Baby Boo Signature",   3, 35.0f, 55.0f, false, -3.0f, 100.0f },
    { "Urban Master Bus",     3, 10.0f, 18.0f, true,  -0.5f,  60.0f },
};
inline constexpr int kNumFactoryPresets = (int) (sizeof (kFactoryPresets) / sizeof (kFactoryPresets[0]));
} // namespace babyair
