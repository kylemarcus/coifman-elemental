/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_IDENTITY_INC
#include ELEM_MAX_INC
#include "coifman.h"
using namespace elem;
using namespace std;

int 
main( int argc, char* argv[] )
{
    
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    
    /*
     * Make sure there are the correct number of arguments.
     * file - input file, 1+ are allowed
     */
    if (argc < 2)
    {
        PRINT("Incorrect Arguments\nUsage: " << argv[0] << " <file> ...");
        return 1;
    }
    
    DEBUG("Scale: " << STR(SCALE));
    DEBUG("Using: " << STR(FP));
    DEBUG("Threads: " << getenv("OMP_NUM_THREADS") << std::endl);

    // input data values from files
    FP * data = NULL;
    
    // distance matrix
    FP * distM = NULL;
    
    // affinity matrix
    FP * affinity = NULL;
    
    // rows and cols of input matrix
    int rows = 0, cols = 0;

    // loop through input files and put them into the data matrix
    DEBUG("There are " << argc-1 << " input files");
    for (int i=1; i<argc; i++ )
    {
        
        std::string line, tmp;
        
        // no data values in input files must be kept track of
        int ndvIndexArr[NDV_BUF];
        int ndvCounter = 0;

        int nx;     // nCols
        int ny;     // nRows
        int xllc;   // nxllcorner
        int yllc;   // nyllcorner
        int dx;     // cellSize
        FP ndv;  // noDataValue

        DEBUG("Working on input file " << argv[i]);

        // open file
        std::ifstream inputf (argv[i]);

        if (inputf.is_open())
        {
            
            DEBUG("File " << argv[i] << " is open");

            // read in header from input file
            inputf >> tmp >> nx;
            inputf >> tmp >> ny;
            inputf >> tmp >> xllc;
            inputf >> tmp >> yllc;
            inputf >> tmp >> dx;
            inputf >> tmp >> ndv;

            rows = ny;
            cols = nx;
            
            DEBUG_VAR(nx);
            DEBUG_VAR(ny);
            DEBUG_VAR(xllc);
            DEBUG_VAR(yllc);
            DEBUG_VAR(dx);
            DEBUG_VAR(ndv);

            // create data array
            if (data == NULL) {
                data = new FP [ nx * ny * (argc-1) ];
                DEBUG("Creating data[" << nx * ny * (argc-1) << "]");
            }

            // read in matrix and find min and max
            FP min, max;
            for (int r=0; r<ny; ++r) {
                
                for (int c=0; c<nx; ++c) {

                    // store column-major order
                    int k = (nx * r) + c + ((i-1)*nx*ny);
                    inputf >> data[k];

                    // check max and min
                    if (r==0 && c==0) {
                        max = min = data[k];
                    }
                    else if ( (data[k] < min) && (data[k] != ndv) ) {
                        min = data[k];
                    }
                    else if ( (data[k] > max) && (data[k] != ndv) ) {
                        max = data[k];
                    }

                    // keep track of ndv index
                    if (data[k] == ndv) {
                        if (ndvCounter < NDV_BUF) {
                            ndvIndexArr[ndvCounter++] = k;
                        } else {
                            PRINT("NDV buffer size exceded!");
                            return 1;
                        }
                    }

                }
                
            }  // done reading in matrix from file
            
            DEBUG_MATRIX(data,cols,rows,"inputTranspose_matrix.csv");

            DEBUG_VAR(min);
            DEBUG_VAR(max);
            DEBUG_VAR(ndvCounter);

            // set all ndv's to min
            for (int k=0; k<ndvCounter; k++) {
                data[ ndvIndexArr[k] ] = min;
            }

            // change numbers to go from 0 to 1
            for (int r=0; r<(nx*ny); r++) {
                data[ r + ((i-1)*nx*ny) ] -= min;
                data[ r + ((i-1)*nx*ny) ] /= max-min;
            }
            
            DEBUG_MATRIX(data, (nx*ny), (argc-1), "data_matrix.csv");

            inputf.close();
            DEBUG("File " << argv[i] << " is closed");

        }

        else
        {
            PRINT("Unable to open file " << argv[i]);
            return 1;
        }

    }
    
    DEBUG("Done with input file");
    
    // create Euclidean distance matrix
    int n = rows * cols;
    DEBUG_VAR(n);
    
    /*
    // in MATLAB, D is distM
    distM = new FP [size_t(n)*size_t(n)]();
    DEBUG("Creating distance matrix (" << STR(FP) << ") " << n << "x" << n);
    
    double tic = omp_get_wtime();
    
    int r;  //rows in dm
    int c;  //cols in dm
    int cd;  //cols in data
    #pragma omp parallel for shared(data, distM) private(r,c,cd) collapse(2)
    for( r=0; r<n; r++ ) { 
        for( c=0; c<n; c++ ) { 
            for( cd=0; cd<(argc-1); cd++ ) { 
                distM[ r + (c*n) ] += std::pow( data[r + (cd*n) ] - data[c + (cd*n)], FP(2) );
            }
            distM[ r + (c*n) ] = std::pow( distM[ r + (c*n) ], FP(0.5) );
            distM[ r + (c*n) ] = std::fabs( distM[ r + (c*n) ] );
        }   
    }
    
    DEBUG("Done creating distance matrix");
    */

    try
    {
        
        //const int n = 10;
        const bool print = true;
        int rank = mpi::WorldRank();

        if( mpi::WorldRank() == 0 )
        {
            const Int commSize = mpi::Size( comm );
            std::cout << "Will create matrices distributed over " 
                      << commSize << " processes" << std::endl;
        }

        const Grid grid( comm );

        // Local buffers
        {
            // Allocate local data
            const Int gridHeight = grid.Height();
            const Int gridWidth = grid.Width();
            const Int gridRow = grid.Row();
            const Int gridCol = grid.Col();
            const Int localHeight = Length( n, gridRow, gridHeight );
            const Int localWidth = Length( n, gridCol, gridWidth );
            std::vector<double> localData( localHeight*localWidth );

            // Fill local data for identity
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                // Form global column index from local column index
                const Int j = gridCol + jLoc*gridWidth;  // global row
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                {
                    // Form global row index from local row index
                    const Int i = gridRow + iLoc*gridHeight;  // global col
                    //if (rank == 0) {
                    //    cout << "j: " << j << " i: " << i << endl;
                    //}
                    
                    // If diagonal entry, set to one, otherwise zero
                    //if( i == j )
                    //    localData[iLoc+jLoc*localHeight] = 1.;
                    //else
                    //    localData[iLoc+jLoc*localHeight] = 0.;
                    
                    int index = iLoc+jLoc*localHeight;
                    int cd;
                    for( cd=0; cd<(argc-1); cd++ ) { 
                        //distM[ r + (c*n) ] += std::pow( data[r + (cd*n) ] - data[c + (cd*n)], FP(2) );
                        localData[index] += std::pow( data[j + (cd*n) ] - data[i + (cd*n)], 2.0 );
                    }
                    //distM[ r + (c*n) ] = std::pow( distM[ r + (c*n) ], FP(0.5) );
                    localData[index] = std::pow( localData[index], 0.5 );
                    //distM[ r + (c*n) ] = std::fabs( distM[ r + (c*n) ] );
                    localData[index] = std::fabs( localData[index] );
                    
                    //localData[iLoc+jLoc*localHeight] = rank;
                }
            }

            DistMatrix<double> distM(grid);
            distM.Attach( n, n, grid, 0, 0, localData.data(), localHeight );
            ValueIntPair<double> max = Max(distM);
            
            double diam_D = max.value;
            double d_k = 2 * std::pow( diam_D/2.0, 2.0 );

            if( print ) {
                if (rank == 0) {
                    cout << "diam_D: " << diam_D << endl;
                    cout << "d_k: " << d_k << endl;
                }
                Write( distM, "distance_matrix", ASCII);
            }
            
            DistMatrix<double> affinity(grid);
            
            for (int scale=0; scale<SCALE; scale++) {
        
                DEBUG_VAR(scale);
        
                // MATLAB: epsilon_s=d_k/(4^s);
                double epsilon_s = d_k / double(std::pow(4,scale));
                DEBUG_VAR(epsilon_s);

                // MATLAB: Dist=(D.^2)./((epsilon_s)^2);
                // MATLAB: affinity = exp(-Dist)
                /*
                int r,c,i;
                DEBUG("Creating affinity matrix");
                #pragma omp parallel for shared(distM, affinity) private(r,c,i) collapse(2)
                for (r=0; r<n; r++) {
                    for (c=0; c<n; c++) {
                        i = r + (c*n);
                        affinity[i] = std::exp( -(std::pow( distM[i], FP(2) ) / std::pow( epsilon_s, FP(2) )) );
                    }
                }
                */
                
                std::vector<double> localDataAffinity( localHeight*localWidth );
                
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                {
                    // Form global column index from local column index
                    const Int r = gridCol + jLoc*gridWidth;  // global row
                    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    {
                        // Form global row index from local row index
                        const Int c = gridRow + iLoc*gridHeight;  // global col
                        
                        int i = iLoc+jLoc*localHeight;
                        
                        localDataAffinity[i] = std::exp( -(std::pow( localData[i], 2.0 ) / std::pow( epsilon_s, 2.0 )) );
                        
                    }
                    
                }
                affinity.Attach( n, n, grid, 0, 0, localDataAffinity.data(), localHeight );
                if (print) Write( affinity, "affinity_matrix", ASCII);
                
                
                DEBUG("Done creating affinity matrix");
            } //end scale for loop
            
            
            
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
