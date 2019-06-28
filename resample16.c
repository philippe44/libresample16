/* resamplesubs.c - sampling rate conversion subroutines */
// Altered version

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "resample16.h"
#include "stdefs.h"
#include "lowfilter.h"
#include "F21T8.h"

/* Check for illegal constants */
#if (Np >= 16)
#error "Error: Np>=16"
#endif
#if (Nb+Nhg+NLpScl >= 32)
#error "Error: Nb+Nhg+NLpScl>=32"
#endif
#if (Nh+Nb > 32)
#error "Error: Nh+Nb>32"
#endif

typedef int sample_t;

typedef struct resample16_s {
	float factor;              /* factor = Sndout/Sndin */
	BOOL fastMode;

	UHWORD LpScl;               /* Unity-gain scale factor */
	UHWORD Nwing;               /* Filter table size */
	UHWORD Nmult;               /* Filter length for up-conversions */
	const HWORD *Imp;           /* Filter coefficients */
	const HWORD *ImpD;          /* ImpD[n] = Imp[n+1]-Imp[n] */
	BOOL Interp;

	HWORD	*Xbuf;
	UHWORD	Ncreep, Xoff;
} resample16_t;

static WORD FilterUp(const HWORD Imp[], const HWORD ImpD[], UHWORD Nwing, BOOL Interp,
		  HWORD *Xp, HWORD Ph, HWORD Inc);

static WORD FilterUD(const HWORD Imp[], const HWORD ImpD[], UHWORD Nwing, BOOL Interp,
		  HWORD *Xp, HWORD Ph, HWORD Inc, UHWORD dhb);

/* fixed point alignment to 16 bits signed values
 */
static INLINE HWORD WordToHword(WORD v, int scl)
{
	WORD llsb = (1<<(scl-1));
	v += llsb;          /* round */
	v >>= scl;
	if (v>MAX_HWORD) {
		return MAX_HWORD;
	} else if (v < MIN_HWORD) {
		return MIN_HWORD;
	}
	return (HWORD) v;
}

#define INTERP(B) 					\
	while (Time < endTime)			\
	{                               \
		iconst = Time & Pmask;		\
		Xp = &B[2*(Time>>Np)];		\
		x1 = *(Xp);					\
		x2 = *(Xp+2);   			\
		x1 *= ((1<<Np)-iconst);		\
		x2 *= iconst;				\
		v = x1 + x2;                \
		*Y++ = WordToHword(v,Np); 	\
		x1 = *(Xp+1);               \
		x2 = *(Xp+3);   		   	\
		x1 *= ((1<<Np)-iconst);		\
		x2 *= iconst;				\
		v = x1 + x2;                \
		*Y++ = WordToHword(v,Np);	\
		Time += dtb;    			\
	}

/* Sampling rate conversion using linear interpolation for maximum speed.
 */
static int resampleFast(resample16_t *r, HWORD *X, int inCount, HWORD *Y)
{
	HWORD iconst;
	HWORD *Xp, *Ystart = Y;
	WORD v,x1,x2;

	float dt;
	UWORD dtb;
	UWORD Time, endTime;

	dt = 1.0/r->factor;            /* Output sampling period */
	dtb = dt*(1<<Np) + 0.5;     /* Fixed-point representation */

	// copy beginning of new block (up to 2*Xoff+creep)
	memcpy(r->Xbuf + (2*r->Xoff - r->Ncreep)*2, X, (2*r->Xoff + r->Ncreep)*2*2);

	// always start @ Xoff and do 2*Xoff data
	Time = (r->Xoff<<Np);
	endTime = Time + ((2*r->Xoff)<<Np);

	INTERP(r->Xbuf);

	// start @ Ncreep in main block
	Time = (r->Xoff + r->Ncreep)<<Np;
	endTime = (inCount - r->Xoff)<<Np;

	INTERP(X);

	// memorize end of block, up to 2*Xoff
	r->Ncreep = r->Xoff - (inCount - (Time>>Np));
	memcpy(r->Xbuf, X + (inCount - 2*r->Xoff + r->Ncreep)*2, (2*r->Xoff - r->Ncreep)*2*2);

	return (Y - Ystart) / 2;
}

#define FILTERUP(F,B)						\
	while (Time < endTime)      		\
	{                              		\
		Xp = &B[2*(Time>>Np)];     		\
		/* Perform left-wing inner product */       					\
		v = F(r->Imp, r->ImpD, r->Nwing, r->Interp, Xp,					\
					(HWORD)(Time&Pmask), -1);							\
		/* Perform right-wing inner product */                          \
		v += F(r->Imp, r->ImpD, r->Nwing, r->Interp, Xp+2,				\
					(HWORD)((((Time)^Pmask)+1)&Pmask), 1);				\
		v >>= Nhg; 														\
		v *= r->LpScl;             										\
		*Y++ = WordToHword(v,NLpScl);   								\
		/* Perform left-wing inner product */       					\
		v = F(r->Imp, r->ImpD, r->Nwing, r->Interp, Xp+1,				\
					(HWORD)(Time&Pmask), -1);							\
		/* Perform right-wing inner product */                          \
		v += F(r->Imp, r->ImpD, r->Nwing, r->Interp, Xp+3,				\
					(HWORD)((((Time)^Pmask)+1)&Pmask), 1);				\
		v >>= Nhg; 														\
		v *= r->LpScl;             										\
		*Y++ = WordToHword(v,NLpScl);   								\
		Time += dtb;    												\
	}

/* Sampling rate up-conversion only subroutine;
 * Slightly faster than down-conversion;
 */

static int resampleUp(resample16_t *r, HWORD X[], int inCount, HWORD Y[])
{
	HWORD *Xp, *Ystart = Y;
	WORD v;

	float dt;                  /* Step through input signal */
	UWORD Time, endTime;
	UWORD dtb;             /* Fixed-point versions of Dh,Dt */

	dt = 1.0/r->factor;            /* Output sampling period */
	dtb = dt*(1<<Np) + 0.5;     /* Fixed-point representation */

	// copy beginning of new block (up to 2*Xoff+creep)
	memcpy(r->Xbuf + (2*r->Xoff - r->Ncreep)*2, X, (2*r->Xoff + r->Ncreep)*2*2);

	// always start @ Xoff and do 2*Xoff data
	Time = (r->Xoff<<Np);
	endTime = Time + ((2*r->Xoff)<<Np);

	FILTERUP(FilterUp, r->Xbuf);

	// start @ Ncreep in main block
	Time = (r->Xoff + r->Ncreep)<<Np;
	endTime = (inCount - r->Xoff)<<Np;

	FILTERUP(FilterUp, X);

	// memorize end of block, up to 2*Xoff
	r->Ncreep = r->Xoff - (inCount - (Time>>Np));
	memcpy(r->Xbuf, X + (inCount - 2*r->Xoff + r->Ncreep)*2, (2*r->Xoff - r->Ncreep)*2*2);

	return (Y - Ystart) / 2;
}

#define FILTERUD(F,B, ...)						\
	while (Time < endTime)      		\
	{                              		\
		Xp = &B[2*(Time>>Np)];     		\
		/* Perform left-wing inner product */       					\
		v = F(r->Imp, r->ImpD, r->Nwing, r->Interp, Xp,					\
					(HWORD)(Time&Pmask), -1, dhb);						\
		/* Perform right-wing inner product */                          \
		v += F(r->Imp, r->ImpD, r->Nwing, r->Interp, Xp+2,				\
					(HWORD)((((Time)^Pmask)+1)&Pmask), 1, dhb);			\
		v >>= Nhg; 														\
		v *= r->LpScl;             										\
		*Y++ = WordToHword(v,NLpScl);   								\
		/* Perform left-wing inner product */       					\
		v = F(r->Imp, r->ImpD, r->Nwing, r->Interp, Xp+1,				\
					(HWORD)(Time&Pmask), -1, dhb);						\
		/* Perform right-wing inner product */                          \
		v += F(r->Imp, r->ImpD, r->Nwing, r->Interp, Xp+3,				\
					(HWORD)((((Time)^Pmask)+1)&Pmask), 1, dhb);			\
		v >>= Nhg; 														\
		v *= r->LpScl;             										\
		*Y++ = WordToHword(v,NLpScl);   								\
		Time += dtb;    												\
	}


/*
	Sampling rate conversion subroutine */
static int resampleUD(resample16_t *r, HWORD X[], int inCount, HWORD Y[])
{
	HWORD *Xp, *Ystart = Y;
	WORD v;

	float dh;                  /* Step through filter impulse response */
	float dt;                  /* Step through input signal */
	UWORD Time, endTime;
	UWORD dhb, dtb;             /* Fixed-point versions of Dh,Dt */

	dt = 1.0/r->factor;            /* Output sampling period */
	dtb = dt*(1<<Np) + 0.5;     /* Fixed-point representation */

	dh = MIN(Npc, r->factor*Npc);  /* Filter sampling period */
	dhb = dh*(1<<Na) + 0.5;     /* Fixed-point representation */

	// copy beginning of new block (up to 2*Xoff+creep)
	memcpy(r->Xbuf + (2*r->Xoff - r->Ncreep)*2, X, (2*r->Xoff + r->Ncreep)*2*2);

	// always start @ Xoff and do 2*Xoff data
	Time = (r->Xoff<<Np);
	endTime = Time + ((2*r->Xoff)<<Np);

	FILTERUD(FilterUD, r->Xbuf);

	// start @ Ncreep in main block
	Time = (r->Xoff + r->Ncreep)<<Np;
	endTime = (inCount - r->Xoff)<<Np;

	FILTERUD(FilterUD, X);

	// memorize end of block, up to 2*Xoff
	r->Ncreep = r->Xoff - (inCount - (Time>>Np));
	memcpy(r->Xbuf, X + (inCount - 2*r->Xoff + r->Ncreep)*2, (2*r->Xoff - r->Ncreep)*2*2);

	return (Y - Ystart) / 2;
}

static WORD FilterUD( const HWORD Imp[], const HWORD ImpD[],
		     UHWORD Nwing, BOOL Interp,
		     HWORD *Xp, HWORD Ph, HWORD Inc, UHWORD dhb)
{
    HWORD a;
	const HWORD *Hp, *Hdp, *End;
    WORD v, t;
    UWORD Ho;
    
    v=0;
    Ho = (Ph*(UWORD)dhb)>>Np;
    End = &Imp[Nwing];
    if (Inc == 1)		/* If doing right wing...              */
    {				/* ...drop extra coeff, so when Ph is  */
	End--;			/*    0.5, we don't do too many mult's */
	if (Ph == 0)		/* If the phase is zero...           */
	  Ho += dhb;		/* ...then we've already skipped the */
    }				/*    first sample, so we must also  */
				/*    skip ahead in Imp[] and ImpD[] */
    if (Interp)
      while ((Hp = &Imp[Ho>>Na]) < End) {
	  t = *Hp;		/* Get IR sample */
	  Hdp = &ImpD[Ho>>Na];  /* get interp (lower Na) bits from diff table*/
	  a = Ho & Amask;	/* a is logically between 0 and 1 */
	  t += (((WORD)*Hdp)*a)>>Na; /* t is now interp'd filter coeff */
	  t *= *Xp;		/* Mult coeff by input sample */
	  if (t & 1<<(Nhxn-1))	/* Round, if needed */
	    t += 1<<(Nhxn-1);
	  t >>= Nhxn;		/* Leave some guard bits, but come back some */
	  v += t;			/* The filter output */
	  Ho += dhb;		/* IR step */
	  Xp += 2*Inc;		/* Input signal step. NO CHECK ON BOUNDS */
      }
    else 
      while ((Hp = &Imp[Ho>>Na]) < End) {
	  t = *Hp;		/* Get IR sample */
	  t *= *Xp;		/* Mult coeff by input sample */
	  if (t & 1<<(Nhxn-1))	/* Round, if needed */
	    t += 1<<(Nhxn-1);
	  t >>= Nhxn;		/* Leave some guard bits, but come back some */
	  v += t;			/* The filter output */
	  Ho += dhb;		/* IR step */
	  Xp += 2*Inc;		/* Input signal step. NO CHECK ON BOUNDS */
      }
    return(v);
}

static WORD FilterUp(const HWORD Imp[], const HWORD ImpD[],
		     UHWORD Nwing, BOOL Interp,
		     HWORD *Xp, HWORD Ph, HWORD Inc)
{
	const HWORD *Hp, *Hdp = NULL, *End;
	HWORD a = 0;
	WORD v, t;
	int count = 0;

	v=0;
	Hp = &Imp[Ph>>Na];
	End = &Imp[Nwing];
	if (Interp) {
	Hdp = &ImpD[Ph>>Na];
	a = Ph & Amask;
	}
	if (Inc == 1)		/* If doing right wing...              */
	{				/* ...drop extra coeff, so when Ph is  */
	End--;			/*    0.5, we don't do too many mult's */
	if (Ph == 0)		/* If the phase is zero...           */
	{			/* ...then we've already skipped the */
		Hp += Npc;		/*    first sample, so we must also  */
		Hdp += Npc;		/*    skip ahead in Imp[] and ImpD[] */
	}
	}
	if (Interp)
      while (Hp < End) {
	  t = *Hp;		/* Get filter coeff */
	  t += (((WORD)*Hdp)*a)>>Na; /* t is now interp'd filter coeff */
	  Hdp += Npc;		/* Filter coeff differences step */
	  t *= *Xp;		/* Mult coeff by input sample */
	  if (t & (1<<(Nhxn-1)))  /* Round, if needed */
	    t += (1<<(Nhxn-1));
	  t >>= Nhxn;		/* Leave some guard bits, but come back some */
	  v += t;			/* The filter output */
	  Hp += Npc;		/* Filter coeff step */
	  Xp += 2*Inc;		/* Input signal step. NO CHECK ON BOUNDS */
      } 
    else 
      while (Hp < End) {
	  t = *Hp;		/* Get filter coeff */
	  t *= *Xp;		/* Mult coeff by input sample */
	  if (t & (1<<(Nhxn-1)))  /* Round, if needed */
	    t += (1<<(Nhxn-1));
	  t >>= Nhxn;		/* Leave some guard bits, but come back some */
	  v += t;			/* The filter output */
	  Hp += Npc;		/* Filter coeff step */
	  Xp += 2*Inc;		/* Input signal step. NO CHECK ON BOUNDS */
      }
    return(v);
}

int resample16(struct resample16_s *r, HWORD X[], int inCount, HWORD Y[])
{
	if (r->fastMode)
		return resampleFast(r, X, inCount, Y);
	else if (r->factor >= 1)
		return resampleUp(r, X, inCount, Y);
	else
		return resampleUD(r, X, inCount, Y);
}

struct resample16_s*  resample16_create(float factor, resample16_filter_e filter, resample16_filter_t *custom, BOOL interp)
{
	resample16_t *r = malloc(sizeof(resample16_t));

	r->Ncreep = 0;
	r->factor = factor;

	if (filter == RESAMPLE16_BASIC) {
		r->fastMode = 1;
		r->Xoff = 10;
	} else {
		r->fastMode = 0;
		r->Interp = interp;
		if (custom) {
			r->Nmult = custom->Nmult;
			r->Imp = custom->Imp;
			r->ImpD = custom->ImpD;
			r->LpScl = custom->LpScl;
			r->Nwing = custom->Nwing;
		} else if (filter == RESAMPLE16_LOW) {
			r->Nmult = LOW_FILTER_NMULT;
			r->Imp = LOW_FILTER_IMP;
			r->ImpD = LOW_FILTER_IMPD;
			r->LpScl = LOW_FILTER_SCALE;
			r->Nwing = LOW_FILTER_NWING;
		} else if (filter == RESAMPLE16_MED) {

			r->Nmult = F21T8_NMULT;
			r->Imp = F21T8_IMP;
			r->ImpD = F21T8_IMPD;
			r->LpScl = F21T8_SCALE;
			r->Nwing = F21T8_NWING;
		} else {
			return NULL;
		}

		// Calc reach of LP filter wing & give some creeping room
		r->Xoff = ((r->Nmult+1)/2.0) * MAX(1.0,1.0/r->factor) + 10;

		// reduce clipping probability
		r->LpScl *= 0.95f;
		if (r->factor < 1) r->LpScl = r->LpScl*r->factor + 0.5;
	}

	// 4*Xoff frames of 2*2 bytes each
	r->Xbuf = calloc(r->Xoff * 4, 4);

	return r;
};

void resample16_delete(struct resample16_s *r) {
	r->Ncreep = 0;
	if (r->Xbuf) free(r->Xbuf);
	free(r);
}

void resample16_flush(struct resample16_s *r) {
	r->Ncreep = 0;
	memset(r->Xbuf, 0, r->Xoff * 4 * 4);
}

