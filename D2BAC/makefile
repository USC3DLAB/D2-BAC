SYSTEM     = ultrasparc64_8_6.2
LIBFORMAT  = static_pic_mt
YANG_D2BAC =/home/semihatakan/D2BAC

#CCC = /opt/SUNWspro/SC5.0/bin/CC
CCC = /opt/SUNWspro/bin/CC
CC  = /opt/SUNWspro/WS6/bin/cc
#CCC = /usr/local/bin/gcc

#------------------------------------------------------------
# Libraries used with examples
#------------------------------------------------------------

LIBS = -lcplex -lm -lsocket -lnsl
#For parallel
#LIBS = -lparcplex -lm -lsocket -lnsl -mt -lmtsk

#------------------------------------------------------------
# Compiler options
#------------------------------------------------------------

OPT    =  -xtarget=ultra -xarch=v9 -KPIC
#For parallel
#OPT    = -xtarget=ultra -xarch=v9 -KPIC

CPPOPT = -g  -xildoff -xtarget=ultra -xarch=v9  #-pto -KPIC -DNDEBUG

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

# The first two definitions are for compilation within the cplex70
#    directory tree.
# CPLEXDIR   = ../../..
# CONCERTDIR = ../../../../concert10
CPLEXDIR   = /usr/bin/cplex/cplex100
CONCERTDIR = /usr/bin/cplex/concert22

CPLEXINCDIR   = $(CPLEXDIR)/include
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CONCERTINCDIR = $(CONCERTDIR)/include
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

EXDIR         = $(CPLEXDIR)/examples
EXSRC         = $(EXDIR)/src
EXINC         = $(EXDIR)/include
EXDATA        = $(EXDIR)/data

CFLAGS  = $(OPT) -I$(CPLEXINCDIR)
LDFLAGS = -L$(CPLEXLIBDIR) $(LIBS)

CPPCFLAGS  =  $(CPPOPT) -I$(CONCERTINCDIR) -I$(CPLEXINCDIR)
CPPLDFLAGS = -L$(CONCERTLIBDIR) -L$(CPLEXLIBDIR) -lilocplex -lconcert $(LIBS)
#CPPLDFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lconcert $(LIBS)

# ------------------------------------------------------------

clean : 
	/bin/rm -rf *.o *~

realclean : 
	/bin/rm -rf *.o *~ 
	/bin/rm -rf *.ord *.sos *.lp *.sav *.net *.msg
# -------------------------------------------------------------


d2bac: d2bacFuncs.o  d2bacMain.o
	$(CCC)   $(CPPCFLAGS) d2bacFuncs.o d2bacMain.o -o d2bac   $(CPPLDFLAGS) 
d2bacMain.o: $(YANG_D2BAC)/ d2bacMain.c
	$(CCC) -c $(CPPCFLAGS) $(YANG_D2BAC)/d2bacMain.c -o d2bacMain.o
d2bacFuncs.o: $(YANG_D2BAC)/ d2bacFuncs.c
	$(CCC) -c $(CPPCFLAGS) $(YANG_D2BAC)/d2bacFuncs.c -o d2bacFuncs.o
d2bac.o: $(YANG_D2BAC)/d2bacFuncs.c d2bacMain.c
	$(CCC) -c $(CPPCFLAGS) $(YANG_D2BAC)/d2bacFuncs.c d2bacMain.c -o d2bacFuncs.o d2bacMain.o
to compile, type at command line prompt: make d2bac

