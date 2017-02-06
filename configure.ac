AC_PREREQ(2.69)
AC_INIT([DESCRIPTION])
AC_PROG_GREP
AC_PROG_AWK

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

FC=`"${R_HOME}/bin/R" CMD config FC`

if test -z "${FC}"; then
  AC_PROG_FC()
fi

AC_MSG_CHECKING(for gfortran)
if test "${FC}" = gfortran; then
  FC_VERSION=`${FC} -v 2>&1 | grep " version " | awk '{print $3;}'`
  AS_VERSION_COMPARE("${FC_VERSION}", "4.8.4",
    [AC_MSG_RESULT(version not supported)
     AC_MSG_ERROR(only gfortran version 4.8.4 or greater can be used)],
    [AC_MSG_RESULT($FC $FC_VERSION)],
    [AC_MSG_RESULT($FC $FC_VERSION)])
else
  AC_MSG_RESULT(no)
  AC_MSG_ERROR(only gfortran version 4.8.4 or greater can be used)
fi
