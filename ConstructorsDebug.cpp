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
using namespace elem;
using namespace std;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        
        const int n = 10;
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
                    if (rank == 0) {
                        cout << "j: " << j << " i: " << i << endl;
                    }
                    // If diagonal entry, set to one, otherwise zero
                    //if( i == j )
                    //    localData[iLoc+jLoc*localHeight] = 1.;
                    //else
                    //    localData[iLoc+jLoc*localHeight] = 0.;
                    localData[iLoc+jLoc*localHeight] = rank;
                }
            }

            DistMatrix<double> X(grid);
            X.Attach( n, n, grid, 0, 0, localData.data(), localHeight );
            if( print )
                Print( X, "Identity constructed from local buffers" );
            
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
