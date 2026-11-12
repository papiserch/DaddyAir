# Calibración de BABY AIR (v1.1 — rediseño de intensidad)

Reproducible: `./build_tests.sh --calibrate` (ajusta `levelMax`) y `./build_tests.sh` (suite completa).

## Diagnóstico de la v1.0 (por qué el efecto era "1 de 10")
Medido con señal tipo voz, ambos knobs al 100 %:

| Modo v1.0 | Aire 9–20 kHz | Causa |
|---|---|---|
| CLEAN | +0.5 dB | nivel de armónicos limitado por un techo de THD (10 %) |
| BRIGHT | +1.9 dB | idem (14 %) |
| BABY BOOST | +7.3 dB (MID) / +3.2 dB (HIGH) | idem (40 %) |

Tres causas: (1) techo de THD autoimpuesto en la calibración; (2) la banda fuente de HIGH (7–18 kHz) generaba armónicos que caen por encima de Nyquist, así que el knob casi no hacía nada (+0.0 a +0.3 dB en CLEAN/BRIGHT); (3) ley de nivel `amt^1.25`, con zona muerta hasta ~50 %.

## Rediseño
- Bandas fuente más bajas: MID HIGH 1.8–5.5 kHz → armónicos 3–10 kHz; HIGH 4.5–11 kHz → armónicos 9–20 kHz (filtro de salida de 24 dB/oct en MID para que no invada el aire).
- Nuevo generador **Chebyshev T2/T3/T4** sobre la señal de banda normalizada por su envolvente: armónicos exactos, independientes del nivel de entrada, sin carácter de saturador. "Densidad": 2.º siempre, 3.º crece lineal y 4.º cuadráticamente con el knob.
- Modos: CLEAN = Chebyshev + algo de tanh; BRIGHT = tanh impar + Chebyshev 3.º/4.º + énfasis de transitorios; VINTAGE = curva asimétrica + rectificador suave, sin Chebyshev; **BABY BOO = sin saturador analógico**: Chebyshev + rectificador, nivelador casi independiente del nivel (rho 0.9), pre-énfasis +10 dB, énfasis de transitorios, sibilancia ×2.
- Ley de nivel lineal en el knob (antes `^1.25`). Rango visual 0–100 % intacto.
- **Protección de picos** (`protectHarmonics`): reduce solo los armónicos añadidos cuando dry+armónicos pasaría de ~−1 dBFS (rodilla suave desde 0.9). Nunca modifica el dry; es protección de picos de muestra, no de true-peak.

## Objetivos de calibración (100 %/100 %, voz sintética −20 dBFS, 48 kHz, 8x)
`levelMax` se obtiene por bisección contra estos objetivos (editables en `kCal` de `dsp_tests.cpp`):
CLEAN +12 dB y BRIGHT +18 dB y BABY BOO +24 dB en la banda de aire 9–20 kHz; VINTAGE +8 dB en presencia 2.5–9 kHz (es un modo de densidad/calidez).
Se probaron objetivos aún mayores (+15/+21/+28), pero se moderaron: ya con los actuales un tono puro tiene THD de 146–429 % al 100 %. **Si al escuchar sigue siendo poco (o demasiado), es un único número por modo.**

## Resultados medidos (voz tipo sintética, ambos knobs)
| Modo | Aire 9–20k @25/50/75/100 % | Presencia 2.5–9k @25/50/75/100 % |
|---|---|---|
| CLEAN | +4.1 / +8.7 / +12.1 / +14.8 | +1.8 / +5.1 / +8.1 / +10.7 |
| BRIGHT | +3.2 / +8.7 / +13.7 / +18.0 | +0.7 / +3.6 / +7.7 / +11.7 |
| VINTAGE | +0.1 / +0.8 / +2.5 / +4.8 | +0.3 / +2.1 / +5.1 / +8.0 |
| BABY BOO | +10.2 / +16.3 / +20.4 / +24.0 | +3.5 / +8.3 / +12.4 / +16.1 |

Independencia (un knob al 100 %, el otro a 0): HIGH solo → presencia ≈ +0.0–0.3 dB, aire +3.0 (VINTAGE) a +24 dB (BABY BOO). MID HIGH solo → presencia +8 a +17.5 dB, aire 13–20 kHz +0.1 a +8.8 dB.

## Otras mediciones
- **Fundamental** (seno −24 dBFS, 100 %): ±0.3 dB en CLEAN/BRIGHT/BABY BOO; VINTAGE hasta +3.9 dB en 3 kHz (curva asimétrica genera fundamental propio); 0.00 dB por debajo de 1 kHz en todos.
- **IMD** (3 kHz + 3.6 kHz, −26 dBFS, productos de 3.er orden respecto a los tonos): CLEAN −30/−18/−6 dB, BRIGHT −24/−13/−2 dB, VINTAGE −30/−21/−13 dB, BABY BOO −23/−12/−1 dB (a 25/50/100 %). **Al 100 % CLEAN no es "limpio" en sentido IMD**: es limpio en carácter (armónicos controlados), no en distorsión por intermodulación.
- **Aliasing audible** (<20 kHz, seno 6.1 kHz, 2x → 8x): CLEAN −45 → −71 dB, BRIGHT −36 → −62, VINTAGE −86 → −98, BABY BOO −36 → −56. El 4.º armónico real cae a ~23.6 kHz por la zona de transición del decimador (inaudible).
- **Nivel (voces sintéticas, 100 %)**: el RMS sube entre +1 y +8 dB (BABY BOO en voz aguda +8.3 dB). Es energía armónica real, pero **al comparar A/B hay que compensar con TRIM** (los tests imprimen el valor para igualar RMS; 100 % en BABY BOO pide −4 a −8 dB).
- **Hi-hats sintéticos**: casi no cambian (+0.1 a +1.3 dB de aire) porque ya tienen su energía en 7–16 kHz; las voces y sintes sí cambian mucho.
- **Sibilancia** (voz + ráfagas de ruido 5.5–9 kHz, BABY BOO 100 %): la cresta del aire baja de 10.1 a 4.7 dB: las ráfagas no sobresalen.
- **Picos**: nunca superan max(pico de entrada, 0.98) por culpa del exciter (comprobado en los 4 modos). Con ruido rosa −18 dBFS y BABY BOO al 100 % la protección ya actúa (pico 0.980): el modo es muy intenso.
- **CPU** (estéreo, 44.1 kHz, un núcleo, máquina de desarrollo): CLEAN 4.6 / 6.5 / 9.0 / 12.7 % y BABY BOO 5.3 / 7.2 / 6.0 / 8.2 % a 2x / 4x / 8x / 16x.

## Lo que NO se ha hecho
Ninguna prueba auditiva. Todo el material es sintético (pulsos glóticos filtrados por formantes, ruido, sierras): sirve para medir, no para juzgar si algo suena bien. Faltan voces reales (masculina/femenina/susurro/adlibs), 808 y mezclas reales, y el ajuste por oído de los objetivos de arriba.
