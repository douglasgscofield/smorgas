// smorgas_global.h  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå
// University
//
// This file contains macros, debug-based includes, and atomic typedefs used
// throughout smorgas.  By convention, STL headers are not included here, nor
// are global classes, rather they are included by the header file for each
// file where they are needed (even if that ends up being every file).
//
// Things here are things you want available "without thinking about it".  This
// file should be fully self-contained, includable by every other header file
// in the project, and none of what it defines should be redundant.


#ifndef _SMORGAS_GLOBAL_H_
#define _SMORGAS_GLOBAL_H_

// Std C/C++ includes
//
// #define NDEBUG  // uncomment to remove assert() code
#include <assert.h>

#ifdef _WITH_DEBUG
#define IF_DEBUG(__lvl__) if (opt_debug >= __lvl__)
#define _DEBUG(__lvl__) if (opt_debug >= __lvl__)
#define DEBUG(__lvl__) (opt_debug >= __lvl__)
#else
#define IF_DEBUG(__lvl__) if (0)
#define _DEBUG(__lvl__) if (0)
#define DEBUG(__lvl__) (false)
#endif

#endif // _SMORGAS_GLOBAL_H_

