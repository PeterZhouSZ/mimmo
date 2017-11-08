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

#include "mimmo_iogeneric.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;



// =================================================================================== //
/*!
 * Reading & managing a quad planar mesh as a degenerate volume mesh.
 */
int test4() {
	
    MimmoGeometry * reader = new MimmoGeometry();
    reader->setIOMode(2);
    reader->setReadDir("geodata");
    reader->setReadFilename("naca0012sq");
    reader->setReadFileType(FileType::SQVTUPLANAR);
    reader->setWriteDir(".");
    reader->setWriteFilename("out-mesh-naca");
    reader->setWriteFileType(FileType::SQVTUPLANAR);
    reader->setCodex();
    reader->execute();
    
    auto obj = reader->getGeometry();
    
    obj->getPatch()->buildAdjacencies();
    obj->getPatch()->buildInterfaces();
    obj->getPatch()->update();
    
    std::cout<<obj->getNCells()<<std::endl;
    
//ROCCO 19/10/2017
//estraggo un bordo dalla mesh SQVTUPLANAR naca e mi porto dietro le sue interfacce
//In un oggetto MimmoObject di typo 4- Curve3D, le normali di qualsiasi genere sono 
//sempre NON definite , i numeri che vengono fuori sono spazzatura.
// L'unica soluzione e': estrarre il bordo dalla mesh originale e salvarlo in un MimmoObject di tipo 4 e 
// contestualmente salvare le normali in un piercedVector chiedendole direttamente alla mesh originale e non quella di bordo.
// Qui sotto un esempio:
    
    //puntatore alla mesh madre obj
    //puntatore alla patchkernel dentro obj: obj->getPatch();
    
    //Extracting free borders of the mesh. 
    bitpit::PiercedVector<Interface> interfaceMeshBorder;
    //instanzio il pierced vector delle normali che saranno riferite a ciascuna cella LINE del bordo objB
    // ma che in realta' per la mesh madre sono delle normali all'interfaccia.
    bitpit::PiercedVector< std::array<double,3> > objB_normals;
    for (const bitpit::Interface & it : obj->getPatch()->getInterfaces())
    {      //cicla sulle interfacce se sono di bordo le inserisce nel piercedvector interfaceMeshBorder
        if(it.isBorder()){
            interfaceMeshBorder.insert(it.getId(), it);
            //valuto normale all'interfaccia 
            darray3E interface_normal = (static_cast<bitpit::VolumeKernel*>(obj->getPatch()))->evalInterfaceNormal(it.getId());
            //... e la salvo nel pierced corrispondente
            objB_normals.insert(it.getId(), interface_normal);
        }
    }

    //instanzio un MimmoObject ti tipo curve3D e lo riempio
    MimmoObject * objB = new MimmoObject(4);
    //visito le interfacce per prendere i vertici e le LINE che le connettono
    for (bitpit::Interface & it0 :interfaceMeshBorder){
        int n = it0.getConnectSize();
        long* con = it0.getConnect(); 
        livector1D cc(n);
        std::set<long int> idcheck;                     //salvo gli id dei vertici in una struttura set per non riperterli
        for(int i=0; i<n; i++){                     //ciclo sulla connettività di ciascuna interfaccia e se il vertice non è ancora stato inserito lo inserisco  
            if(idcheck.count(con[i]) == 0){
                objB->addVertex(obj->getPatch()->getVertexCoords(con[i]), con[i]);
                idcheck.insert(con[i]);
            }
            cc[i] = con[i];
        }
        //aggiungo la cella connessa al mimmoobject curva 3D
        objB->addConnectedCell(cc, bitpit::ElementType::LINE, it0.getId());
    }
    objB->buildAdjacencies();

    //Ora hai un oggetto objB che contiene i bordi, connessi con celle di tipo linea, della tua mesh originale
    // ed un pierced vector objB_normals che contiene tutte le normali ad ogni cella. Puoi entrare nel PiercedVector
    // con l'id della cella dato da objB e ottenere la normale a quella cella. Tieni conto che le normali sono ora definite:
    //
    //                  normale
    //                  ^
    //                  |
    //                  |
    //   x------------------------------x
    //  V1            cella id         V2 
    //
    // adesso per trovare le normali sui vertici di objB devo fare questo:
    // dato l'id di un vertice di objB, trovo il ring di celle attorno al vertice
    // chiedendolo a objB. Hai una serie di id di cella che formano il ring. Recuperi
    // le normali delle celle nel ring visitando il pierced vector, e fai una media tra le
    // normali cella per ottenere la normale al vertice. 
    // Di seguito viene riportato un esempio per recuperare le normali ai vertici e posizionarle
    // in un pierced vector riferito ai vertici di objB.
    
    bitpit::PiercedVector< std::array<double,3> > objB_vertexNormals;
    
    for(const bitpit::Cell & cell : objB->getCells()){ //ciclo su ogni cella
        long idC = cell.getId();
        const long * conn = cell.getConnect();
        std::size_t sizeConn = cell.getVertexCount();
        
        for(std::size_t i =0; i<sizeConn; ++i){
        
            if(!objB_vertexNormals.exists(conn[i])){ // mi accerto che non ho gia' calcolato la normale al vertice
                
                //valuto il ring di celle attorno al vertice conn[i], che appartiene alla cella idC in posizione i;
                std::vector<long> listRing = objB->getPatch()->findCellVertexOneRing(idC, i);
                std::array<double,3> vNormal = {{0,0,0}};
                //medio la normale sul ring.
                for(const long & idRing : listRing){
                    vNormal += objB_normals[idRing];
                }
                vNormal /= double(listRing.size()); //media
                
                //salvo la normale nel Pierced Vector objB_vertexNormals. 
                objB_vertexNormals.insert(conn[i], vNormal);
            }//endif
        }//end loop sul ring
    }//end loop su tutte le celle
    
    // output purpose only
    std::vector<std::array<double,3> > datanormals;
    datanormals.reserve(objB_normals.size());
    for(const std::array<double,3> & normal : objB_normals ){
        datanormals.push_back(normal);
    }
    
    std::vector<std::array<double,3> > datavertexnormals;
    datavertexnormals.reserve(objB_vertexNormals.size());
    for(const std::array<double,3> & vnormal : objB_vertexNormals ){
        datavertexnormals.push_back(vnormal);
    }
    
    objB->getPatch()->getVTK().addData("cellNormals", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::CELL, datanormals);
    objB->getPatch()->getVTK().addData("vertexNormals", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, datavertexnormals);
    objB->getPatch()->write("BorderOfNaca");
    objB->getPatch()->getVTK().removeData("cellNormals");
    objB->getPatch()->getVTK().removeData("vertexNormals");
        
    delete reader;
    bool check = true;
    return int(!check);
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

        int val = test4() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
