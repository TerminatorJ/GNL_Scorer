AM_CPPFLAGS = -DPKGDATADIR=\"$(pkgdatadir)/\"

bin_PROGRAMS = hybrid-min hybrid-ss-min ct-energy ct2rnaml

hybrid_min_SOURCES = energy.c hybrid-min.c energy.h getopt.h constraints.h options.h util.h xmalloc.h
hybrid_min_CPPFLAGS = -DENERGY=int -DINTEGER -DPRECISION=1000 $(AM_CPPFLAGS)
hybrid_min_LDADD = getopt.o getopt1.o xmalloc.o -lm

hybrid_ss_min_SOURCES = energy.c hybrid-ss-min.c energy.h getopt.h constraints.h options.h util.h xmalloc.h
hybrid_ss_min_CPPFLAGS = -DENERGY=int -DINTEGER -DPRECISION=1000 $(AM_CPPFLAGS)
hybrid_ss_min_LDADD = getopt.o getopt1.o xmalloc.o -lm

ct_energy_SOURCES = getopt.c getopt1.c energy.c xmalloc.c ct-energy.c energy.h getopt.h util.h xmalloc.h
ct_energy_LDADD = -lm

ct2rnaml_SOURCES = ct2rnaml.c xmalloc.c xmalloc.h
