# The program requires stacksize unlimited
# bash: ulimit -s unlimited
# csh:  limit  stacksize unlimited
C    = gcc
OPT  =	-g -O3  -fexpensive-optimizations
OBJS = AssignStructure.o \
       AssignNumberInStructure.o\
       EffectiveCounts.o \
       MI_Compute.o      \
       Correction_APC.o  \
       MaxEnt.o          \
       Pseudoinverse.o   \
       ReadMSA.o         \
       UseColumns.o      \
       svdcmp.o          \
       sort2.o           \
       readpdb.o         \
       NearestContact.o  \
       contact_map.o     \
       RankOrder.o       \
       Results_MI_DI.o   \
       MakeFileNames.o   \
       Calculate_MI_DI.o \
       GetSeedSequence.o \
       raDI.o

LIBS =  -lm -fstack-protector-all 
  
.SUFFIXES: .c .o
.c.o:
	$(C) -c -O3  -fexpensive-optimizations $*.c 

raDI: $(OBJS) 
	$(C) -o raDI  $(OBJS)  $(LIBS) $(OPT)
clean:
	\rm *.o 

