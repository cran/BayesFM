AC_PREREQ(2.69)
AC_INIT([DESCRIPTION])
AC_PROG_GREP
AC_PROG_SED

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

FC=`"${R_HOME}/bin/R" CMD config FC`

if test -z "${FC}"; then
  AC_MSG_WARN(R could not find a Fortran compiler. Now trying to detect one.)
  AC_PROG_FC
fi

AC_MSG_CHECKING(for Fortran compiler)

case "${FC}" in

  *gfortran*)
    dnl check if GNU Fortran is actually used
    GNU_FORTRAN=`${FC} --version | ${GREP} 'GNU Fortran'`
    if test -z "${GNU_FORTRAN}"; then
      AC_MSG_RESULT(no)
      AC_MSG_ERROR(Fortran compiler (${FC}) does not seem to be GNU Fortran.)
    fi
    dnl minimum version of GNU Fortran that supports code
    MINI_GFORTRAN_VERSION="4.6.3"
    dnl check version of GNU Fortran
    FC_VERSION=`echo "__GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__" | \
                 ${FC} -E -P - | ${SED} -e 's/ //g'`
    AC_MSG_RESULT(gfortran ${FC_VERSION})
    AS_VERSION_COMPARE(${FC_VERSION}, ${MINI_GFORTRAN_VERSION},
      [AC_MSG_ERROR([this version of GNU Fortran cannot be used to compile this package. Please try to upgrade gfortran (>= ${MINI_GFORTRAN_VERSION}) or to use a different Fortran compiler.])], [], [])
    ;;

  *g95*)
    AC_MSG_RESULT(${FC})
    AC_MSG_ERROR([G95 cannot compile the Fortran code contained in this package. Please try to use a different Fortran compiler.])
    ;;

  *)
    AC_MSG_RESULT(${FC})
    AC_MSG_WARN([the Fortran code contained in this package will be compiled with ${FC}. Please contact package maintainer if compilation fails.])
    ;;

esac

AC_OUTPUT
