# BABY AIR™ — Harmonic Exciter (Baby Demon)

Excitador armónico VST3 / AU / Standalone en C++17 + JUCE 8.

## Estado de verificación (leer primero) — v1.1 "rediseño de intensidad"
| Parte | Estado |
|---|---|
| Motor DSP (`Source/DSP`, sin JUCE) | Rediseñado, compilado con g++ 13 (`-Wall -Wextra`, sin warnings) y **probado**: `./build_tests.sh` → *ALL TESTS PASSED* |
| Capa JUCE y interfaz | Solo cambió el texto «BABY BOOST» → «BABY BOO» y los nombres de preset. **Sin compilar en este entorno** (sin JUCE/CMake/red): la compilación real ocurre en GitHub Actions |
| VST3 | **No existe aún un binario verificado.** Lo genera el workflow al subir el código (ver abajo) |
| Escucha real | **No realizada.** Pruebas con material *sintético*; ver `docs/CALIBRATION.md` |

## Compilar
Requisitos: CMake ≥ 3.22, compilador C++17 (Xcode / VS 2022 / GCC 11+), git (para descargar JUCE).
```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release            # descarga JUCE 8.0.6
# o con JUCE local:  cmake -B build -DJUCE_PATH=/ruta/a/JUCE
cmake --build build --config Release
```
Salidas en `build/BabyAir_artefacts/Release/` (VST3, Standalone; AU solo en macOS). Copia el `.vst3` a
`%COMMONPROGRAMFILES%\VST3` (Windows), `~/Library/Audio/Plug-Ins/VST3` (macOS) o `~/.vst3` (Linux).
AAX no incluido (requiere SDK de Avid). Solo DSP, sin JUCE: `./build_tests.sh` (o `-DBABY_AIR_BUILD_PLUGIN=OFF`).

## Compilar en GitHub (sin instalar nada)
1. Crea un repositorio vacío en GitHub y sube el contenido de esta carpeta:
   ```sh
   git init && git add . && git commit -m "BABY AIR 1.0.0"
   git branch -M main
   git remote add origin https://github.com/TU_USUARIO/baby-air.git
   git push -u origin main
   ```
2. Pestaña **Actions**: el workflow *Build BABY AIR* se ejecuta solo (también con *Run workflow*).
   Primero corre los tests DSP y, si pasan, compila en Windows, macOS (universal arm64+x86_64) y Linux.
3. Descarga los plugins en **Artifacts** de la ejecución (`BABY_AIR-Windows`, `-macOS`, `-Linux`).
4. Para publicar una release con los zips: `git tag v1.0.0 && git push origin v1.0.0`.

Si un job falla, abre el log del paso *Build*: los errores de compilación de la capa JUCE (nunca compilada hasta ahora) aparecerán ahí y son lo primero que hay que corregir.
Los binarios macOS salen sin firmar ni notarizar y los de Windows sin firmar: el sistema operativo mostrará avisos de seguridad.
**Licencia:** JUCE es AGPLv3 o comercial; distribuir BABY AIR requiere cumplir una de las dos. Este repositorio no incluye ningún archivo LICENSE: elige el tuyo.

## Controles
| Control | Función |
|---|---|
| ON | Bypass real con transición suave (la salida sigue retardada por la latencia, en sincronía con la compensación del host) |
| AIR MODE | CLEAN / BRIGHT / VINTAGE / BABY BOO. Crossfade de 40 ms; los valores de los knobs se conservan |
| MID HIGH (0–100 %) | Excitación de la región 2–7.5 kHz (armónicos hasta ~18 kHz). 0 % = ninguna excitación |
| HIGH (0–100 %) | Excitación de la región 7–18 kHz |
| LINK | Al arrastrar un knob, el otro se mueve conservando la diferencia (el knob arrastrado se limita para que nunca cambie). Solo actúa en gestos de usuario, no con automatización del DAW |
| TRIM | Salida ±18 dB (0 dB por defecto) |
| MIX (0–100 %) | Dry/Wet. `wet = dry + armónicos`, ambas rutas alineadas a la muestra |
| OS 2x/4x/8x/16x | Oversampling real (por defecto 8x). Cambiarlo cambia la latencia (se informa al host) |
Medidores IN/OUT: RMS (300 ms) + pico + peak hold (1.5 s) + LED de clip (clic para borrar). Presets de fábrica: menú de programas del DAW.

## Arquitectura
`Parameters` → `PluginProcessor` → `ExciterEngine` (suavizado, bypass, crossfade de modo, MIX/TRIM) → `ModeVoice` (detector de sibilancia + 2 `BandProcessor` por canal) → `Oversampler` + `AdaaShaper`.
- La señal original nunca se filtra: solo se **suma** contenido armónico. Por eso no hay huecos ni picos de cruce.
- Cada banda: aislamiento Butterworth 4º → pre-énfasis → nivelador dependiente de la señal → 3 caminos no lineales en paralelo (tanh asimétrico, tanh impar, rectificador suave) → solo el **residuo no lineal** (se resta el fundamental con la describing function exacta de tanh) → high-pass de salida + tilt + DC blocker.
- Oversampling: half-band FIR de fase lineal en cascada, latencia entera exacta (47 muestras a 8x/44.1 kHz = 1.07 ms). ADAA de 1er orden en 1x/2x/4x; a 8x/16x solo oversampling (ahorra CPU).
- Sin asignaciones, bloqueos ni E/S en el callback; atomics para medidores. Sin limitador: solo una cota de seguridad de ±16 (+24 dBFS) y saneado de NaN/Inf. No compilar con `-ffast-math`.
- Un cambio de modo o de calidad recalcula coeficientes (decenas de µs) en el hilo de audio: sin asignaciones, pero con un coste puntual.

## Los cuatro modos (v1.1)
| Modo | Identidad | Técnica |
|---|---|---|
| CLEAN | PURE AIR | armónicos Chebyshev exactos (2.º/3.º) + algo de tanh, pre-énfasis +4 dB |
| BRIGHT | CRYSTAL SHINE | tanh impar + Chebyshev 3.º/4.º, pre-énfasis +7 dB, énfasis de transitorios |
| VINTAGE | ANALOG GLOW | curva asimétrica (pares) + tanh impar + rectificador suave, filtro de armónicos, compresión suave |
| BABY BOO | DIGITAL AIR | sin saturador analógico: síntesis Chebyshev T2/T3/T4 + rectificador, pre-énfasis +10 dB, transitorios, sibilancia ×2 |

MID HIGH = región de presencia (armónicos 3–10 kHz); HIGH = aire (9–20 kHz). Al 100 % la intensidad es **extrema a propósito** (BABY BOO ≈ +24 dB de aire sobre voz sintética): el RMS sube hasta ~8 dB, compensa con TRIM al comparar.
Los IDs de parámetros y el orden de los modos no cambian: los proyectos de DAW guardados siguen cargando (el modo 4 solo cambia de etiqueta).

## Resultados medidos y diagnóstico
Ver `docs/CALIBRATION.md` (diagnóstico de la v1.0, tablas 25/50/75/100 %, IMD, aliasing, CPU y limitaciones).

## Limitaciones conocidas
- **Sin pruebas auditivas ni material real.** Los objetivos de intensidad son estimaciones; ajustables en un solo número por modo (`kCal` + `--calibrate`).
- Al 100 % el IMD es alto (BABY BOO ≈ −1 dB, CLEAN ≈ −6 dB en dos tonos): intensidad extrema por diseño. A 25–50 % es mucho menor.
- Hi-hats con energía ya en 7–16 kHz cambian poco (+0.1 a +1.3 dB de aire).
- La protección de picos es de picos de muestra, no de true-peak.
- El logo es un dibujo vectorial aproximado. Latencia ≈ 1 ms (el host debe soportar PDC).
- La etiqueta de versión de la interfaz sigue en «v1.0.0» para no tocar la UI.
