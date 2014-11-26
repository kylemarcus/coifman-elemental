#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <vector>
#include <typeinfo>
#include <time.h>
#include <limits>
#include <cfloat>
#include <cmath>
#include <ctgmath>
#include "math.h"



/* configuration */

#define DEBUG_MATRIX_ENABLED 0
#define DEBUG_PRINT_ENABLED 1
#define USE_DOUBLE

//#define SCALE 6
#define SCALE 1
#define NDV_BUF 1024



/* macros */

// floating point size used within the program.
#ifdef USE_DOUBLE
    #define FP double
#else
    #define FP float
#endif

// debug matrix by calling a function to write matrix to a file
//#define DEBUG_MATRIX(args...) do { if (DEBUG_MATRIX_ENABLED) printMatrix(##args); } while (0)
#define DEBUG_MATRIX 

// debug output to stderr
#define DEBUG(x) do { if (DEBUG_PRINT_ENABLED) std::cerr << x << std::endl; } while (0)
#define DEBUG_VAR(x) do { if (DEBUG_PRINT_ENABLED) std::cerr << #x << ": " << x << std::endl; } while (0)

// print macros for stdout
#define PRINT(x) do { std::cout << x << std::endl; } while (0)
#define PRINT_VAR(x) do { std::cout << #x << ": " << x << std::endl; } while (0)

// creates a string
#define STR2(x) #x
#define STR(y) STR2(y)

// Uses Box Muller's algo to emulate MATLAB's randn function.  This is a modified
// version of the one found at the address below.  I modified it to use the
// thread save version of rand function so I could use it with OpenMP.
// http://www.developpez.net/forums/d544518/c-cpp/c/equivalent-randn-matlab-c/
#define TWOPI (6.2831853071795864769252867665590057683943387987502)
// RAND is a macro All which returns a pseudo-random number from a uniform
// distribution on the interval [0 1]
#define RAND (rand_r(&seed)) / ((double) RAND_MAX)
// randn macro returns a pseudo-random numbers from a normal
// distribution with mean zero and standard deviation one.
//  This macro uses Box Muller's algorithm
#define randn (sqrt (-2.0 * log (RAND)) * cos (TWOPI * RAND))



/* debug functions */
 
// Mainly used for debug, this function prints a 
// col-major ordered matrix to a csv file
void printMatrix(const FP * matrix, int rows, int cols, char * fname) {
    
    std::ofstream csv (fname);
    for( int i=0; i<rows; i++ ) {
        for( int j=0; j<cols; j++ ) {
            //std::cout << matrix[ i + (j*rows)] << " ";
            csv << matrix[ i + (j*rows)] << ",";
        }
        //std::cout << std::endl;
        csv << std::endl;
    }
    csv.close();
    
}