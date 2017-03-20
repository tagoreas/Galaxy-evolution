#ifndef GALEVOL_FUNCS_TEMPLATES_CPP_
#define GALEVOL_FUNCS_TEMPLATES_CPP_

#include "galevol_funcs.hpp"
#include <unistd.h>
#include <stdlib.h>
#include <sstream>
#include <string>

template <typename T>
std::string galevol_funcs::ge_tostring (T a)
{
  std::string Result2;
  std::ostringstream convert2;
  convert2 << a;
  Result2 = convert2.str();
  return Result2;
}

template <typename T>
void galevol_funcs::ge_malloc( T **ptr, size_t dim1 )
{
    size_t total = (size_t)dim1;

    *ptr = (T*)NULL;

    if(!total/* || total > CONSTANT uintmax/sizeof(T) */)
    {
        return;
    }

    do
    {
        *ptr = (T*)malloc( dim1 * sizeof(T) );

        if( !(*ptr) )
        {
            std::cerr << "pixsrc :: WARNING: trying to malloc "
                      << 1.0e-9*(double)(dim1 *sizeof(T)) << " Gbytes\n"
                "waiting for memory to be available .." << std::endl;
            usleep (10000000ll);
        }
    }
    while( !(*ptr) );
}

template <typename T>
void galevol_funcs::ge_malloc( T ***ptr, size_t dim1 , size_t dim2 )
{
    size_t total = (size_t)(dim1*dim2);

    *ptr = (T**)NULL;

    if(!total/* || total > CONSTANT uintmax/sizeof(T) */)//if( !total || total > CONSTANT uintmax/sizeof(T*) )
    {
        return;
    }

    do
    {
        *ptr = (T**)malloc( dim1 * sizeof(T*) );

        if( !(*ptr) )
        {
            std::cerr << "pixsrc :: WARNING: trying to malloc "
                      << 1.0e-9*(dim1 *sizeof(T*)) << " Gbytes\n"
                "waiting for memory to be available .." << std::endl;
            usleep (10000000ll);
        }

    }
    while( !(*ptr) );

    for(size_t j=0; j<dim1; ++j)
    {
        ge_malloc<T>( &((*ptr)[j]), dim2 );
    }
}

template <typename T>
void galevol_funcs::ge_malloc( T ****ptr, size_t dim1 , size_t dim2 , size_t dim3 )
{
    size_t total = (size_t)(dim1*dim2*dim3);

    *ptr = (T***)NULL;

    if(!total/* || total > CONSTANT uintmax/sizeof(T) */)//if( !total || total > CONSTANT uintmax/sizeof(T**) )
    {
        return;
    }

    do
    {
        *ptr = (T***)malloc( dim1 * sizeof(T**) );

        if( !(*ptr) )
        {
            std::cerr << "pixsrc :: WARNING: trying to malloc "
                      << 1.0e-9*(dim1 *sizeof(T**)) << " Gbytes\n"
                "waiting for memory to be available .." << std::endl;
            usleep (10000000ll);
        }

    }
    while( !(*ptr) );

    for (size_t j=0; j<dim1; ++j)
    {
        ge_malloc<T>( &((*ptr)[j]), dim2, dim3 );
    }
}

template <typename T>
void galevol_funcs::ge_malloc( T *****ptr, size_t dim1 , size_t dim2 , size_t dim3, size_t dim4 )
{
    size_t total = (size_t)(dim1*dim2*dim3*dim4);

    *ptr = (T****)NULL;

    if(!total/* || total > CONSTANT uintmax/sizeof(T) */)//if( !total || total > CONSTANT uintmax/sizeof(T***) )
    {
        return;
    }

    do
    {
        *ptr = (T****)malloc( dim1 * sizeof(T***) );

        if( !(*ptr) )
        {
            std::cerr << "pixsrc :: WARNING: trying to malloc "
                      << 1.0e-9*(dim1 *sizeof(T***)) << " Gbytes\n"
                "waiting for memory to be available .." << std::endl;
            usleep (10000000ll);
        }
    }
    while( !(*ptr) );

    for (size_t j=0; j<dim1; ++j)
    {
        ge_malloc<T>( &((*ptr)[j]), dim2, dim3, dim4 );
    }
}

template <typename T>
void galevol_funcs::ge_free( T *ptr )
{
    free( ptr );
    ptr = (T*)NULL;
}

template <typename T>
void galevol_funcs::ge_free( T **ptr, size_t dim1 )
{
    for (size_t j=0; j<dim1; ++j)
        ge_free<T>( ptr[j] );
    free( ptr );
    ptr = (T**)NULL;
}

template <typename T>
void galevol_funcs::ge_free( T ***ptr, size_t dim1, size_t dim2 )
{
    for (size_t j=0; j<dim1; ++j)
        ge_free<T>( ptr[j], dim2 );
    free( ptr );
    ptr = (T***)NULL;
}

template <typename T>
void galevol_funcs::ge_free( T ****ptr, size_t dim1, size_t dim2, size_t dim3 )
{
    for (size_t j=0; j<dim1; ++j)
        ge_free<T>( ptr[j], dim2, dim3 );
    free( ptr );
    ptr = (T****)NULL;
}

#endif
