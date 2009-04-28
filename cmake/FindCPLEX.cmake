FIND_PATH(CPLEX_INCLUDE_DIR
  ilcplex/cplex.h
  PATHS "C:/ILOG/CPLEX91/include")

FIND_LIBRARY(CPLEX_LIBRARY
  NAMES cplex91
  PATHS "C:/ILOG/CPLEX91/lib/msvc7/stat_mda")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

FIND_PATH(CPLEX_BIN_DIR
  cplex91.dll
  PATHS "C:/ILOG/CPLEX91/bin/x86_win32")

IF(CPLEX_FOUND)
  SET(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
  SET(CPLEX_LIBRARIES ${CPLEX_LIBRARY})
ENDIF(CPLEX_FOUND)

MARK_AS_ADVANCED(CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_BIN_DIR)

IF(CPLEX_FOUND)
  SET(LEMON_HAVE_LP TRUE)
  SET(LEMON_HAVE_MIP TRUE)
  SET(LEMON_HAVE_CPLEX TRUE)
ENDIF(CPLEX_FOUND)
