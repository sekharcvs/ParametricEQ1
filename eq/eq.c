#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "options.h"
#include "constants.h"
#include "fft.h"
#include "rom.h"

int main(int argc, char** argv)
{
    FILE *fin, *fout, *feq;
    int Fs;
    short i, j, frameLen, olaLen, winLen, n, b;
#ifdef IO_SHORT
    short inputPCM[MS2NSAMPLES(48000, FRAME_SIZE_MS)], outputPCM[MS2NSAMPLES(48000, FRAME_SIZE_MS)];
#else
    float inputPCM[MS2NSAMPLES(48000, FRAME_SIZE_MS)], outputPCM[MS2NSAMPLES(48000, FRAME_SIZE_MS)];
#endif
    float fftBuffReal[MS2NSAMPLES(48000, WIN_SIZE_MS)/2 + 1], fftBuffImag[MS2NSAMPLES(48000, WIN_SIZE_MS)/2 + 1];
    float eqGains[N_BANDS], temp;
    short bandEdgesFs[N_BANDS + 1], n_bands;

#ifdef DEBUGGING
    FILE *f1, *f2, *f3;
    float winBuff[MS2NSAMPLES(48000, WIN_SIZE_MS)];
    long cnt = 0;
#endif

    float inpMem[MS2NSAMPLES(48000, OLA_SIZE_MS)], outMem[MS2NSAMPLES(48000, OLA_SIZE_MS)], *win;

    kiss_fft_cfg cfg, cfg1;
    kiss_fft_cpx *out, *inp;

    if (argc != 5)
    {
        fprintf(stderr, "please check number of input arguements\n");
        exit(-1);
    }

    /* Parse Args */
    /* Args order: exe_name inputStream.extension outputStream.extension eqGainsStream.extension Fs */
    /* eqGainsStream.extension is a binary file with N_BANDS floating point values (32-bit per sample) corresponding to each frame */

    fin = fopen(argv[1], "rb");
    fout = fopen(argv[2], "wb");
    feq = fopen(argv[3], "rb");

#ifdef DEBUGGING
    f1 = fopen("windowedSignal.pcm", "wb");
    f2 = fopen("fftSignalReal.pcm", "wb");
    f3 = fopen("fftSignalImag.pcm", "wb");
#endif

    Fs = (int)atoi(argv[4]);

    for (i = 0; i < N_BANDS; i++)
    {
        eqGains[i] = 0;
    }

    switch (Fs)
    {
    case 48000:
        win = hannWin48k;
        break;
    case 44100:
        win = hannWin44k1;
        break;
    case 32000:
        win = hannWin32k;
        break;
    case 16000:
        win = hannWin16k;
        break;
    case 8000:
        win = hannWin8k;
        break;
    default:
        fprintf(stderr, "Not an allowed sample rate\n");
        exit(-1);
        break;
    }

    n_bands = 0;
    while (1)
    {

    }

    /* Calculate band edges in terms of DFT bins at the current sampling rate */
    bandEdgesFs[0] = 0;
    for (i = 1; i < N_BANDS + 1; i++)
    {
        bandEdgesFs[i] = (short)(roundf(bandEdges[i] * Fs/48000.0f));
        bandEdgesFs[i] = (bandEdgesFs[i] > bandEdgesFs[i - 1]) ? (bandEdgesFs[i]) : (bandEdgesFs[i] + 1);
    }

    frameLen = MS2NSAMPLES(Fs, FRAME_SIZE_MS);
    olaLen   = MS2NSAMPLES(Fs, OLA_SIZE_MS);
    winLen = frameLen + olaLen;

    cfg  = kiss_fft_alloc(winLen, 0, NULL, 0);
    cfg1 = kiss_fft_alloc(winLen, 1, NULL, 0);

    inp = (kiss_fft_cpx *)malloc(winLen * sizeof(kiss_fft_cpx));
    out = (kiss_fft_cpx *)malloc(winLen * sizeof(kiss_fft_cpx));

    if (inp == NULL || out == NULL || cfg == NULL)
    {
        fprintf(stderr, "malloc did not succeed.. exiting\n");
        exit(-1);
    }

    for (i = 0; i < olaLen; i++)
    {
        inpMem[i] = 0;
        outMem[i] = 0;
    }

    while (1)
    {
#ifdef IO_SHORT
        n = (short)fread(inputPCM, sizeof(short), frameLen, fin);
#else
        n = (short)fread(inputPCM, sizeof(float), frameLen, fin);
#endif
        if (n < frameLen)
        {
            break;
        }

        if (!feof(feq))
        {
            /* In case feq does not specify eq gains for all frames, or rather is just a single frame's entry, just use the last valid set of eq gains for all frames */
            n = fread(eqGains, sizeof(float), N_BANDS, feq);
        }
        if (n != N_BANDS)
        {
            fprintf(stderr, "Number of eq gains provided needs to be a multiple of %d.. exiting\n", N_BANDS);
            exit(-1);
        }

        for (i = 0, j = 0; i < olaLen; i++, j++)
        {
            (inp[i]).r = win[j] * inpMem[i]; (inp[i]).i = 0;
        }

        for (i = olaLen; i < frameLen; i++)
        {
            (inp[i]).r = (float)(inputPCM[i - olaLen]); (inp[i]).i = 0;
        }

        for (i = frameLen, j = 0; i < winLen; i++, j++)
        {
            (inp[i]).r = sqrtf(1 - win[j] * win[j]) * inputPCM[i - olaLen]; (inp[i]).i = 0;
            inpMem[i - frameLen] = inputPCM[i - olaLen];
        }

        kiss_fft(cfg, inp, out);

        /* EQ Processing goes here */

        for (i = 0; i < winLen / 2 + 1; i++)
        {
            fftBuffReal[i] = out[i].r;
            fftBuffImag[i] = out[i].i;
        }

        for (b = 0; b < N_BANDS; b++)
        {
            /* Band loop */
            temp = powf(10.0f, eqGains[b] / 20.0f);
            for (i = bandEdgesFs[b]; i < bandEdgesFs[b + 1]; i++)
            {
                /* Bin loop */
                fftBuffReal[i] *= temp;
                fftBuffImag[i] *= temp;
            }
        }

        /* Real signal fft buffer generation */
        out[0].r = fftBuffReal[0];
        out[0].i = 0;

        out[winLen / 2].r = fftBuffReal[winLen / 2];
        out[winLen / 2].i = 0;

        for (i = 1; i < MS2NSAMPLES(Fs, WIN_SIZE_MS) / 2; i++)
        {
            out[i].r = fftBuffReal[i];
            out[i].i = fftBuffImag[i];
        }

        for (i = winLen - 1, j = 1; i > winLen / 2; i--, j++)
        {
            out[i].r = fftBuffReal[j];
            out[i].i = -fftBuffImag[j];
        }


#ifdef DEBUGGING
        fwrite(winBuff, sizeof(float), winLen, f1);
        fwrite(fftBuffReal, sizeof(float), winLen, f2);
        fwrite(fftBuffImag, sizeof(float), winLen, f3);
#endif
        /***************************/

        kiss_fft(cfg1, out, inp);

        for (i = 0, j = 0; i < olaLen; i++, j++)
        {
#ifdef IO_SHORT
            outputPCM[i] = (short)(roundf(win[j] * (inp[i]).r + (sqrtf(1 - win[j] * win[j])) * outMem[i])/winLen);
#else
            outputPCM[i] = (win[j] * (inp[i]).r + (sqrtf(1 - win[j] * win[j])) * outMem[i])/winLen;
#endif
        }

        for (i = olaLen; i < frameLen; i++)
        {
#ifdef IO_SHORT
            outputPCM[i] = (short)(roundf((inp[i]).r)/winLen);
#else
            outputPCM[i] = ((inp[i]).r)/winLen;
#endif
        }

        for (i = frameLen; i < winLen; i++)
        {
            outMem[i - frameLen] = (inp[i]).r;
        }

#ifdef IO_SHORT
        fwrite(outputPCM, sizeof(short), frameLen, fout);
#else
        fwrite(outputPCM, sizeof(float), frameLen, fout);
#endif

#ifdef DEBUGGING
        fprintf(stdout, "%-8ld\b\b\b\b\b\b\b\b", cnt);
        cnt++;
#endif
    }

    free(cfg);
    free(cfg1);
    free(out);
    free(inp);

    fclose(fin);
    fclose(fout);
    fclose(feq);

#ifdef DEBUGGING
    fclose(f1);
    fclose(f2);
    fclose(f3);
#endif



    return 0;
}


