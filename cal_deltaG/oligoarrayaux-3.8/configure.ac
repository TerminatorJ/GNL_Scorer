AC_PREREQ(2.57)
AC_INIT(OligoArrayAux, 3.8, markhn@rpi.edu)
AC_CONFIG_SRCDIR(src/energy.c)
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h:config.in)
AM_INIT_AUTOMAKE
AC_CANONICAL_BUILD
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PATH_PROG(PERL, perl)
if test -z "$PERL"; then
   AC_MSG_ERROR(Perl not found)
fi
AC_PROG_PERL_VERSION(5.6.1, , AC_MSG_ERROR(Perl 5.6.1 or better is required))

AC_MSG_CHECKING(whether compiler needs -mieee)
case $build/$CC in
  (alpha*/gcc*)
    CFLAGS="$CFLAGS -mieee"
    AC_MSG_RESULT(yes)
  ;;
  (*)
    AC_MSG_RESULT(no)
  ;;
esac

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS(fcntl.h ieeefp.h libintl.h limits.h stdlib.h string.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.

AC_MSG_CHECKING(whether scanf() understands "inf")
AC_RUN_IFELSE(AC_LANG_PROGRAM([[#include <stdio.h>]], [[double d; sscanf("inf", "%lg", &d); return finite(d);]]),
              [INFINITY="inf"; AC_MSG_RESULT(yes)],
	      [INFINITY="999999"; AC_MSG_RESULT(no)],
	      [INFINITY="999999"; AC_MSG_RESULT(can't tell, assuming no)])
AC_SUBST(INFINITY)

AC_ARG_ENABLE(force,
	      AC_HELP_STRING([--enable-force], [enable support for force constraints (off by default)]),
	      [if test "X$enableval" = "Xyes"; then
	       enable_force=1
	      else
	       enable_force=0
	      fi],
	      enable_force=0)
AC_DEFINE_UNQUOTED(ENABLE_FORCE, $enable_force, [Enable force constraints])

# Variables for man pages
AC_SUBST(DATE, 'May 2010')
AC_SUBST_FILE(REFERENCES)
REFERENCES=$srcdir/doc/references

AC_CONFIG_FILES(Makefile src/Makefile scripts/Makefile data/Makefile doc/Makefile)
AC_CONFIG_FILES(scripts/melt.pl, [chmod +x scripts/melt.pl])
AC_CONFIG_FILES(scripts/h-num.pl, [chmod +x scripts/h-num.pl])
AC_CONFIG_FILES(scripts/ss-count.pl, [chmod +x scripts/ss-count.pl])
AC_CONFIG_FILES(doc/hybrid-min.1 doc/hybrid-ss-min.1 doc/ct-energy.1 doc/melt.pl.1 doc/h-num.pl.1 doc/ss-count.pl.1 doc/ct2rnaml.1)
AC_OUTPUT
