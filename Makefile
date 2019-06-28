EXAMPLE = resampler
LIBRARY = libresample16.a
FILTER  = buildfilter

BASE		= ..
CFLAGS		+= 
SRC			= .
DEFINES 	=
OBJ			= .

vpath %.c $(SRC)

INCLUDE = -I. \
		  -I$(BASE)
 
SOURCES_EXP = resample.c 
SOURCES_LIB = resample16.c
SOURCES_FLT	= windowfilter.c filterkit.c
		
OBJECTS_EXP	= $(patsubst %.c,$(OBJ)/%.o,$(filter %.c,$(SOURCES_EXP))) 
OBJECTS_LIB	= $(patsubst %.c,$(OBJ)/%.o,$(filter %.c,$(SOURCES_LIB))) 
OBJECTS_FLT	= $(patsubst %.c,$(OBJ)/%.o,$(filter %.c,$(SOURCES_FLT))) 

all: $(LIBRARY) $(EXAMPLE) $(FILTER)

$(EXAMPLE): $(OBJECTS_EXP)
	$(CC) $(OBJECTS_EXP) $(LIBRARY) -o $@
	
$(FILTER): $(OBJECTS_FLT)
	$(CC) $(OBJECTS_FLT) -lrt -lm -o $@

$(LIBRARY): $(OBJECTS_LIB)	
	$(AR) -rcs $(LIBRARY) $(OBJECTS_LIB)

$(OBJ)/%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) $< -c -o $@

clean:
	rm -f *.o *.a

