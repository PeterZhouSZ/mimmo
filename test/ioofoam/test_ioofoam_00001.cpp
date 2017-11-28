/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
 \ *---------------------------------------------------------------------------*/

#include "IOOFOAM.hpp"
#include <exception>

// =================================================================================== //
int test1() {
    
    std::cout<<"Waiting for a proper test. I do nothing for now"<<std::endl;
    return 0;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {
    
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
    
    #if ENABLE_MPI==1
    MPI::Init(argc, argv);
    
    {
        #endif
        /**<Calling mimmo Test routines*/
        
        int val = test1() ;
        
        #if ENABLE_MPI==1
    }
    
    MPI::Finalize();
    #endif
    
    return val;
}
