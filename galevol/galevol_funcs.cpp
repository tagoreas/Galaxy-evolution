#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include "galevol_funcs_templates.cpp"


using std::ifstream;                            // reading in file
using std::stringstream;

void galevol_funcs::readfile(const char *filename, char ***vec, int *numlines)
{
    // removes comments signaled with a '#' and whitespace

    // open file
    ifstream file;
    file.open(filename);
    if(file.fail())
    {
        GEFUNCS ge_malloc( &(*vec), 0 );
        *numlines = 0;
        return;
    }

    char *buffer, *ptr;
    int size;
    GEFUNCS ge_malloc( &buffer, 100000 );

    // find out how many "good" lines in file
    *numlines = 0;
    while(!file.eof())
    {
        file.getline( buffer, 100000 );

        // can't just call sizestring because need to exclude comments
        size = 0;
        ptr = buffer;

        // remove comments from string
        while( *ptr )
        {
            if( *ptr=='#' )
            {
                *ptr = 0;
                break;
            }

            ++size;
            ++ptr;
        }

        // remove white space
        GEFUNCS trim( buffer, &size );

        if( size )
        {
            ++*numlines;
        }
    }
    file.close();

    if( !*numlines )
    {
        GEFUNCS ge_free( buffer );
        return;
    }

    // alloctae memory for file read in
    GEFUNCS ge_malloc( &(*vec), *numlines );

    // rewind file buffer to beginning of file
    // file.seekg( 0, std::ios::beg );
    file.open( filename );

    // actually read in file to memory
    int position = 0;
    while(!file.eof())
    {
        file.getline( buffer, 100000 );

        // can't just call sizestring because need to exclude comments
        size = 0;
        ptr = buffer;

        // remove comments from string
        while( *ptr )
        {
            if( *ptr=='#' )
            {
                *ptr = 0;
                break;
            }

            ++size;
            ++ptr;
        }

        // remove white space
        GEFUNCS trim( buffer, &size );

        // write to memory
        if( size )
        {
            GEFUNCS ge_malloc( &((*vec)[position]), size+1 );
            strcpy( (*vec)[position], buffer );
            ++position;
        }
    }
    file.close();

    GEFUNCS ge_free( buffer );
}
void galevol_funcs::trim( char *line, int *size0 )
{
    // size does NOT include terminating null character

    int size = GEFUNCS sizestring( line );

    // remove whitespace at end of string
    while( size && ( line[size-1]=='\n' || line[size-1]=='\r' ||
                     line[size-1]==' '  || line[size-1]=='\t'    ) )
    {
        line[size-1] = 0;
        --size;
    }

    // remove whitespace at beginning of string
    while( size && ( line[  0   ]=='\n' || line[  0   ]=='\r' ||
                     line[  0   ]==' '  || line[  0   ]=='\t'    ) )
    {
        // shift entire string over, not copying first character
        char *ptr1 = line;
        char *ptr2 = ptr1+1;
        for( int j=0; j<size; ++j)
        {
            *ptr1 = *ptr2;
            ptr1 = ptr2;
            ++ptr2;
        }
        --size;
    }

    if( size0 )
    {
        *size0 = size;
    }
}
int galevol_funcs::sizestring( const char *str )
{
    // size does NOT include terminating character

    int size = 0;
    const char *ptr = str;
    while( *ptr )
    {
        ++size;
        ++ptr;
    }
    return size;
}
void galevol_funcs::split( const char *line0, const char *delimiter0, char ***out, int *num )
{
    // this function will trim whitespace from ends of string before splitting it
    // a dleimiter of " " will match spaces and tabs
    const char *delimiter = strcmp (delimiter0, " ") ? delimiter0 : " \t";

    // copy string to split
    int sizestr = GEFUNCS sizestring( line0 );
    char line[sizestr+1];
    strcpy( line, line0 );

    // trim whitespace
    GEFUNCS trim( line, &sizestr );

    // find number of strings result after splitting
    char *ptr;
    ptr = strtok( line, delimiter );
    int num0 = 0;
    while( ptr )
    {
        ++num0;
        ptr = strtok( NULL, delimiter );
    }
    if( num )
        *num = num0;

    GEFUNCS ge_malloc( &(*out), num0 );

    // split strings

    // re-copy because strtok destroys string
    strcpy( line, line0 );
    ptr = strtok( line, delimiter );
    int position = 0;
    while( ptr )
    {
        GEFUNCS ge_malloc( &((*out)[position]), GEFUNCS sizestring(ptr)+1 );
        strcpy( (*out)[position], ptr );
        ++position;

        ptr = strtok( NULL, delimiter );
    }
}

double galevol_funcs::convert_string (string a)
{
    double i;
    std::stringstream out(a);
    out >> i;
    return i;
}
