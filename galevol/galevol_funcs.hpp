#ifndef GALEVOL_FUNCS_HPP_
#define GALEVOL_FUNCS_HPP_

#define GEFUNCS galevol_funcs::

#include <vector>
#include <iostream>

using std::string;
using std::vector;

class galevol_funcs
{

public:

    static void   trim       (char*, int*);
    static int    sizestring (const char*);
    static void   split      (const char*, const char*, char***, int*);
    static void   readfile   (const char*, char***, int*);
    static double convert_string (string);

  template <typename T> static string ge_tostring (T);

    template <typename T> static void ge_malloc( T**   , size_t                );
    template <typename T> static void ge_malloc( T***  , size_t, size_t           );
    template <typename T> static void ge_malloc( T**** , size_t, size_t, size_t      );
    template <typename T> static void ge_malloc( T*****, size_t, size_t, size_t, size_t );
    template <typename T> static void ge_free  ( T*                         );
    template <typename T> static void ge_free  ( T**   , size_t                );
    template <typename T> static void ge_free  ( T***  , size_t, size_t           );
    template <typename T> static void ge_free  ( T**** , size_t, size_t, size_t      );
};

#endif
