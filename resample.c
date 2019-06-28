/*
 * FILE: resample.c
 * Sampling-rate-conversion main program (command line usage)
 */

static char resampleVersion[]
	= "\n\tresample version 0.1 (Feb. 1, 2006 - jos@ccrma.stanford.edu)\n\n\
Copyright 1982-2006 by Julius Smith.\n\
Interlaved buffer by philippe_44@outlook.com.\n\
This is free software. See the Lesser GNU Public License (LGPL) for copying conditions.\n\
There is NO warranty;  not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\
Modification for interleaved buffer & misc philippe_44@outlook.comn\
";

#define USAGE "\
\n\
USAGE: One of the following:\n\
\n\
	  resample -to srate [-noFilterInterp] [-quality <b|l|q>] [-f filterFile] [-terse] input output\n\
  	  \tb: linear interpolation, l:low 13 taps, m: medium 21 taps\n\
	  resample -by factor [options as above] input output\n\
	  resample -version\n\
	  input & output must be using wav format\
\n\
Options can be abbreviated.\n\n\
Report bugs to <bug-resample@w3k.org>.\n\n\
"

#include "resample16.h"
#include "stdefs.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>

static int trace = 1;		/* controls verbosity of output */

static char comment[256] = "";

static void fail(char *s)
{
    fprintf(stderr,"\n*** resample: %s \n",s); /* Display error message  */
    fprintf(stderr,USAGE);
    exit(1);			/* Exit, indicating error */
}

static void fails(char *s, char *s2)
{
    printf("resample: ");           /* Display error message  */
    printf(s,s2);
    printf("\n\n");
    exit(1);                        /* Exit, indicating error */
}

typedef unsigned char u8_t;

static struct wave_header_s {
	u8_t 	chunk_id[4];
	u8_t	chunk_size[4];
	u8_t	format[4];
	u8_t	subchunk1_id[4];
	u8_t	subchunk1_size[4];
	u8_t	audio_format[2];
	u8_t	channels[2];
	u8_t	sample_rate[4];
	u8_t    byte_rate[4];
	u8_t	block_align[2];
	u8_t	bits_per_sample[2];
	u8_t	subchunk2_id[4];
	u8_t	subchunk2_size[4];
} wave_in, wave_out = {
		{ 'R', 'I', 'F', 'F' },
		{ 0x24, 0xff, 0xff, 0xff },
		{ 'W', 'A', 'V', 'E' },
		{ 'f','m','t',' ' },
		{ 16, 0, 0, 0 },
		{ 1, 0 },
		{ 2, 0 },
		{ 0x44, 0xac, 0x00, 0x00  },
		{ 0x10, 0xb1, 0x02, 0x00 },
		{ 4, 0 },
		{ 16, 0 },
		{ 'd', 'a', 't', 'a' },
		{ 0x00, 0xff, 0xff, 0xff },
	};

#define BUFSIZE (2*4096)

int main(int argc, char *argv[])
{
	double factor = 0;	/* factor = Sndout/Sndin */
	double insrate, newsrate=0;
	BOOL interpFilt = TRUE;	/* TRUE means interpolate filter coeffs */
	int quality = RESAMPLE16_BASIC;
	BOOL knowFactor = FALSE;	/* Used to detect insufficient command-line spec */
	int inCount, outCount, outCountReal;
	int nChans, result;
	int inType, inFormat;
	int outType, outFormat;
	FILE *in, *out;
	struct resample16_s *r;
	HWORD *ibuf = malloc(BUFSIZE);
	HWORD *obuf = malloc(BUFSIZE);
	int rate;

	struct stat statbuf;
	char *insfname, *outsfname, *argv0;
	char filterFile[512] = "";

	if (argc == 1) {
	fprintf(stderr, USAGE);
	exit(1);
	}

	argv0 = argv[0];
	while (--argc && **(++argv)=='-') {
	++(argv[0]); /* skip over '-' */
	switch (*argv[0]) {
	case 'b':			       /* -by factor */
		if (--argc)
		  sscanf(*(++argv),"%lf",&factor);
		if (trace)
		  printf("Sampling-rate conversion factor set to %f.\n",factor);
		knowFactor = TRUE;
		break;
	case 'f':  			       /* -filter filterFile */
		if (--argc)
		strcpy(filterFile, *++argv);
		else
		fail("Need to specify filter file name");
		if (trace)
		  printf("Filter file set to %s.\n",filterFile);
		break;
	case 'n':			       /* -noFilterInterpolation */
		interpFilt = FALSE;
		if (trace)
		  printf("Filter-table interpolation disabled.\n");
		break;
	case 'q':
		if (--argc)
		{
		  char qual;
		  sscanf(*(++argv),"%c",&qual);
		  if (qual == 'b') quality = RESAMPLE16_BASIC;
		  else if (qual == 'l') quality = RESAMPLE16_LOW;
		  else if (qual == 'm') quality = RESAMPLE16_MED;
		}
		if (trace)
		  printf("Quality %d", quality);
		break;
	case 't':
		if (*(argv[0]+1) == 'e') { 		/* -terse */
		trace = 0;
		break;
		}
		if (--argc)				/* -to srate */
		  sscanf(*(++argv),"%lf",&newsrate);
		if (trace)
		  printf("Target sampling-rate set to %f.\n",newsrate);
		knowFactor = TRUE;
		break;
	case 'v':			       /* -version */
		printf(resampleVersion);
		if (argc == 1)
		exit(0);
		break;
	default:
		fprintf(stderr,"Unknown switch -%s\n",argv[0]);
		fprintf(stderr,USAGE);
		exit(1);
	}
	}

	if (!knowFactor)
	  fail("Must specify sampling-rate conversion factor via '-to' or '-by' option");

	if (argc < 1)
	  fail("Need to specify input soundfile");
	insfname = *argv;

	if (argc < 2) {
	fprintf(stderr, USAGE);
	exit(1);
	}
	else
	  outsfname = *++argv;

	in = fopen(insfname, "rb");
	out = fopen(outsfname, "wb");

	fread(&wave_in, sizeof(wave_in), 1, in);
	insrate = wave_in.sample_rate[2] << 16 | wave_in.sample_rate[1] << 8 | wave_in.sample_rate[0];
	rate = factor ? rate * factor : newsrate;
	wave_out.sample_rate[0] = rate;
	wave_out.sample_rate[1] = rate >> 8;
	wave_out.sample_rate[2] = rate >> 16;
	fwrite(&wave_out, sizeof(wave_out), 1, out);

	r = resample16_create(newsrate / insrate, quality, NULL, interpFilt);

	while (1) {
		int n = fread(ibuf, 4, BUFSIZE/4, in);
		if (!n) break;
		n = resample16(r, ibuf, n, obuf);
		fwrite(obuf, 4, n, out);
	}

	resample16_delete(r);

	fclose(in);
	fclose(out);

	return(0);
}
