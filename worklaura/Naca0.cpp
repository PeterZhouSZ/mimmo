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
#include "mimmo_core.hpp"
#include "myHeader.hpp"
#include "Operators.hpp"
#include "bitpit.hpp"

#include "mimmo_manipulators.hpp"
#include "mimmo_utils.hpp"
#include "mimmo_geohandlers.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;

// =================================================================================== //
 

 /*Voglio creare una funzione che effettua una deformazione con l'RBF
 */
void evalDef (MimmoGeometry*  geom ){
  
//     geom->setIOMode(IOMode::WRITE);
//     geom->setWriteDir(".");
//     geom->setWriteFileType(FileType::SQVTUPLANAR);
//     geom->setWriteFilename("laura_func_def_00001.0000");

    /* Create MRBF block.
     */
    dvecarr3E rbfNodes(1,{{1.0,0.0,5.3522e-08}});
    dvecarr3E displ(1,{{0.0,0.0,-0.05}});
    
   MRBF* mrbf = new MRBF();
    mrbf->setMode(MRBFSol::NONE);
    mrbf->setSupportRadius(0.0003);
    mrbf->setPlotInExecution(true);
    mrbf->setNode(rbfNodes);
    mrbf->setFunction(bitpit::RBFBasisFunction::LINEAR);
    mrbf->setDisplacements(displ);
    
    /* Create applier block.
     */
    Apply* applier = new Apply();
        
    addPin(geom, mrbf, PortType::M_GEOM, PortType::M_GEOM);
    addPin(geom, applier, PortType::M_GEOM, PortType::M_GEOM);
    addPin(mrbf, applier, PortType::M_GDISPLS, PortType::M_GDISPLS);
     
    Chain ch0;
    ch0.addObject(geom);
    ch0.addObject(mrbf);
    ch0.addObject(applier);
    ch0.exec(true);
    
    delete applier;
    delete mrbf;
    
    applier = NULL;
    mrbf = NULL;
    
    return;  
}



// =================================================================================== //


void test1(){
    MimmoGeometry * mimmo0 = new MimmoGeometry();
    mimmo0->setIOMode(IOMode::READ);
    mimmo0->setReadDir("./geodata");
    mimmo0->setReadFilename("naca0012sq");
    mimmo0->setReadFileType(FileType::SQVTUPLANAR);
    mimmo0->execute();
    
    MimmoGeometry * mimmo0_Copy = new MimmoGeometry();
    mimmo0_Copy->setIOMode(IOMode::READ);
    mimmo0_Copy->setReadDir("./geodata");
    mimmo0_Copy->setReadFilename("naca0012sq");
    mimmo0_Copy->setReadFileType(FileType::SQVTUPLANAR);
    mimmo0_Copy->execute();
    
    MimmoObject* obj_u = mimmo0->getGeometry();
    obj_u->getPatch()->write("mimmoPrima");
    obj_u->buildAdjacencies();					//costruisco le adiacenze  
    obj_u->buildInterfaces();					//costruisco le interfacce
  
    evalDef(mimmo0_Copy);
    MimmoObject* obj_d = mimmo0_Copy->getGeometry();
    obj_d->getPatch()->write("mimmoDef");
    obj_d->buildAdjacencies();					//costruisco le adiacenze  
    obj_d->buildInterfaces();					//costruisco le interfacce
  
       //plotting quaterniorns
    bitpit::PiercedVector<double> ortRbf = evalCellOrth(obj_d);
    dvector1D tempOrtRbf = setCounter(ortRbf);
    plotDataQUAD("./", "OrthRBF", obj_d, tempOrtRbf, "Orthogonality");
    
    exit(1);
    /*Devo etrarre i bordi del corpo deformato dovrei estrarre solo il profilo
     */
    MimmoObject * objBord = new MimmoObject(4,2);
    PatchKernel* m_patch = obj_u->getPatch();
    extractBorders (m_patch, objBord);					//Estraendo ora i bordi del corpo non deformato
    objBord->getPatch()->write("objBord");//Bordi del corpo non deformato
    
    PatchKernel* m_patch1 = obj_d->getPatch();
 
  bitpit::PiercedVector<bitpit::Interface>  it_0 = getBorderInterfaces(obj_d->getPatch());				//Estraggo le interfacce di bordo del MimmoObject deformato
  bitpit::PiercedVector<bitpit::Interface>  it_1;									//Sto inizializzando questo pierced vector per memorizzare le interfacce della cella che ho eliminato
  for (bitpit::Interface &i : it_0){											//Ciclo sulle interfacce di bordo
    int n= i.getConnectSize();												//Restituisce la size della connettivitàdell'elemento
    long int* con = i.getConnect();											//E' un puntatore alla connettività di un elemento 
    bool check = false;
    int count = 0;				//Inizializzo un contatore per poter aggiungere l'elemento(interfaccia) se tutti i suoi vertici (elementi della connettività) rispettano contemporaneamente dei limiti che andrò ad imporre 
    for(long j=0; j<n; j++){												//Ciclo sulla connettività
      std::array<double,3> a ;												//salvo in a le coordinate del vertice di posto j della connettività
      a = obj_u->getPatch()->getVertexCoords(con[j]);
      if( a[0]>-1.5 && a[0]<1.5 && a[2]>-0.5 && a[2]<0.5)			//Impongo ch esia x che y appartengono all'insieme [6,15]
	count++;
	if (count == n){												//Se ciascun vertice che appartiene all'interfaccia rispetta queste limitazioni
	  it_1.insert(i.getId(), i);											//Allora l'interfaccia viene inserita
	}
    }
  }   

  MimmoObject * obj_round = new MimmoObject(4,2);			//Inizializzo un nuovo MimmoObject che rappresenterà la porzione di bordo ruotato
  dvecarr3E vertices(4, {{0,0,0}});
  livector2D c(4, livector1D(2,0));
  for (bitpit::Interface & it :it_1){ 					//Ciclo sul PiercedVector che contiene le interfacce del quadrato che ho eliminato
    int n= it.getConnectSize();
    long* con = it.getConnect(); 
    livector1D cc(n);							//Vector cove andrò a salvare la connettività di ciascuna interfaccia
    std::set<long int> idcheck;						//salvo gli id dei vertici in una struttura set per non riperterli
    for(int i=0; i<n; i++){						//ciclo sulla connettività di ciascuna interfaccia e se il vertice non è ancora stato inserito lo inserisco  
      if(idcheck.find(con[i])== idcheck.end()){
	obj_round->addVertex(m_patch1->getVertexCoords(con[i]), con[i]);
	idcheck.insert(con[i]);
      }
      cc[i] = con[i];
    }
   obj_round->addConnectedCell(cc, bitpit::ElementType::LINE, it.getId());//Aggiungo l'elemento con connettività cc
  }
  obj_round->buildAdjacencies();  
  
//   obj_round->getPatch()->write("obj_round");//E' l'ala deformata
  completeBorders(objBord,obj_round);
  obj_round->getPatch()->write("obj_roundCompleteBord");//E' l'ala deformata

    bitpit::PiercedVector<Interface> interfaceMeshBorder;
    bitpit::PiercedVector< std::array<double,3> > objB_normals;    //instanzio il pierced vector delle normali che saranno riferite a ciascuna cella LINE del bordo objB per la mesh madre sono normali all'interfaccia
   
   for (const bitpit::Interface & it : obj_u->getPatch()->getInterfaces())
    {      //cicla sulle interfacce se sono di bordo le inserisce nel piercedvector interfaceMeshBorder
        if(it.isBorder()){
            interfaceMeshBorder.insert(it.getId(), it);
            //valuto normale all'interfaccia 
            darray3E interface_normal = (static_cast<bitpit::VolumeKernel*>(obj_u->getPatch()))->evalInterfaceNormal(it.getId());
            //... e la salvo nel pierced corrispondente
            objB_normals.insert(it.getId(), interface_normal);
        }
    }
    
 matchConnectivity(objBord, obj_round);					//Mi assicuro che la connettività sia uguale   

 //Calcolo le normali ai vertici facendo una media sulle normali alle interfacce
        bitpit::PiercedVector< std::array<double,3> > objB_vertexNormals;

    for(const bitpit::Cell & cell : objBord->getCells()){ //ciclo su ogni cella//////////////////////////////////////////////////////////////////////////////objB-->objBord
        long idC = cell.getId();
        const long * conn = cell.getConnect();
        std::size_t sizeConn = cell.getVertexCount();
        
        for(std::size_t i =0; i<sizeConn; ++i){
        
            if(!objB_vertexNormals.exists(conn[i])){ // mi accerto che non ho gia' calcolato la normale al vertice
                
                //valuto il ring di celle attorno al vertice conn[i], che appartiene alla cella idC in posizione i;
                std::vector<long> listRing = objBord->getPatch()->findCellVertexOneRing(idC, i);//////////////////////////////////////////////////////////////objB-->objBord
                std::array<double,3> vNormal = {{0,0,0}};
                //medio la normale sul ring.
                for(const long & idRing : listRing){
                    vNormal += objB_normals[idRing];
                }
                vNormal /= double(listRing.size()); //media
		vNormal /= -norm2(vNormal);///////////////////////////////////////////////////////////ho aggiunto un meno per invertre il verso delle normali
                
                //salvo la normale nel Pierced Vector objB_vertexNormals. 
                objB_vertexNormals.insert(conn[i], vNormal);
            }//endif
        }//end loop sul ring
    }//end loop su tutte le celle    
 
//    {//FOR DEBUG ONLY visualizzo le normali ai vertici su objBord
//       std::vector<std::array<double,3> > datanormals;
//       datanormals.reserve(objB_vertexNormals.size());
//       
//       for(const std::array<double,3> & normal : objB_vertexNormals ){
// 	datanormals.push_back(normal);
//       }
//   
//       objBord->getPatch()->getVTK().addData("vertexNormals", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, datanormals); 
//       objBord->getPatch()->write("./objBord_visualizzazioneNormali");
//       objBord->getPatch()->getVTK().removeData("vertexNormals");
//    }

  
   PatchKernel* m_patchB = objBord->getPatch();
   PatchKernel* m_patchR = obj_round->getPatch();
   bitpit::PiercedVector< std::array<double,3> > objR_normals;    //instanzio il pierced vector delle normali che saranno riferite a ciascuna cella LINE del bordo objB per la mesh madre sono normali all'interfaccia
    
    for(const bitpit::Cell & cell : objBord->getCells()){ //ciclo su ogni cella//////////////////////////////////////////////////////////////////////////////objB-->objBord
      long idC = cell.getId();
        const long * conn = cell.getConnect();
	std::array < double,3 > normal_def;  
	std::array < double,3 > n_u = objB_normals[idC];
	//Differenza di coordinate tra i due vertici che compongono la cella considerata indeformata e poi normalizzo il vettore ottenuto
	std::array<double, 3> v_u = m_patchB->getVertexCoords(conn[1])-m_patchB->getVertexCoords(conn[0]);
	double normv_u = norm2(v_u);
	if (normv_u > 1.E-18){
	    v_u /= normv_u;
	} 
	//Differenza di coordinate tra i due vertici che compongono la cella considerata deformata e poi normalizzo il vettore ottenuto
	std::array<double, 3> v_d = m_patchR->getVertexCoords(conn[1])-m_patchR->getVertexCoords(conn[0]);
	double normv_d = norm2(v_d);
	if ( normv_d > 1.E-18){
	    v_d /= normv_d;
	} 
	
	//Calcolo il coseno tra due vettori normalizzati
	double Dx = dotProduct ( v_d, v_u );
	
	double Dz = dotProduct(v_d, n_u);
	
	double a = Dx*n_u[0] + Dz*n_u[2];
	double b = n_u[1];
	double c = -1.0*Dz*n_u[0] + Dx*n_u[2];
	
	std::array < double,3 > n_d = {{ a, b, c }};
	objR_normals.insert(idC, n_d);
    }
    
//Calcolo le normali ai vertici facendo una media sulle normali alle interfacce
        bitpit::PiercedVector< std::array<double,3> > objR_vertexNormals;

    for(const bitpit::Cell & cell : obj_round->getCells()){ //ciclo su ogni cella//////////////////////////////////////////////////////////////////////////////objB-->objBord
        long idC = cell.getId();
        const long * conn = cell.getConnect();
        std::size_t sizeConn = cell.getVertexCount();
        
        for(std::size_t i =0; i<sizeConn; ++i){
        
            if(!objR_vertexNormals.exists(conn[i])){ // mi accerto che non ho gia' calcolato la normale al vertice
                
                //valuto il ring di celle attorno al vertice conn[i], che appartiene alla cella idC in posizione i;
                std::vector<long> listRing = obj_round->getPatch()->findCellVertexOneRing(idC, i);//////////////////////////////////////////////////////////////objB-->objBord
                std::array<double,3> vNormal = {{0,0,0}};
                //medio la normale sul ring.
                for(const long & idRing : listRing){
                    vNormal += objR_normals[idRing];
                }
                vNormal /= double(listRing.size()); //media
                vNormal /= -norm2(vNormal);///////////////////////////////////////////////////////////ho aggiunto un meno per invertre il verso delle normali
                //salvo la normale nel Pierced Vector objB_vertexNormals. 
                objR_vertexNormals.insert(conn[i], vNormal);
            }//endif
        }//end loop sul ring
    }//end loop su tutte le celle    

 /* 
   {//FOR DEBUG ONLY 
      std::vector<std::array<double,3> > datanormals_def;
      datanormals_def.reserve(objR_vertexNormals.size());      
      for(const std::array<double,3> & normal : objR_vertexNormals ){
	datanormals_def.push_back(normal);
      }
  
      obj_round->getPatch()->getVTK().addData("vertexNormals", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, datanormals_def); 
      obj_round->getPatch()->write("./obj_round_visualizzazioneNormali");
      obj_round->getPatch()->getVTK().removeData("vertexNormals");
   }
   */
  
    bitpit::PiercedVector<Quaternion> Q, Q_mod;
    bitpit::PiercedVector<Quaternion>  T, T_mod;
    std::array<double, 3> COR{0.0,0.0,0.0};					//Voglio imporre come centro di rotazione il punto (1, 0, 1)

    int flag0 = 0;
    int flag1 = 1;
    int flag2 = 2;
    evalQuaternions( objBord , obj_round, objB_vertexNormals, objR_vertexNormals, COR, Q, T);				//Valuto i quaternioni secondo il centro di rotazione specificato 
    initializeAllQuaternions(obj_u, Q, T);  
//     plotQuaternionField ("./", "Quats_fNaca0", obj_u, Q, T);	//Plotto idue campi di quaternioni inizializzati su obj 
   
//    Contiene gli id dei vertici di bordo perchè lì non voglio modificare il quaternioni perchè rappresentano le C.I.

    std::vector<long int> idB = objBord->getVertices().getIds();		//vector che contiene solo id di obj2 che coincide con gli id dei vertici di bordo   
    bitpit::PiercedVector< double > err_trasl0, err_trasl1;
    bitpit::PiercedVector< double > err_rot0, err_rot1;
    
    for ( long & i: Q.getIds() ){
      Q_mod.insert( i, Q[i]);
      T_mod.insert( i,T[i]);
    }
    //Smoothing considerado le lunghezze modificate
    
//     std::unordered_map< long int , std::vector<long> >   mapOneRing_interfaces;
//     std::unordered_map< long int , std::vector<long> > mapOneRing = evalVertexOneRingAdj( obj_u , mapOneRing_interfaces);
    
    std::unordered_map< long int , std::vector<long> > mapOneRing =evalVertexOneRing(obj_u);
    std::unordered_map< long int , std::vector<double> > coefficientiModificati = mapCoeffModDumping( obj_u, objBord, mapOneRing ,20, 4.0, flag1);
    int n = 2;

       
  for ( long & vv :obj_u->getVertices().getIds() ){
    double err = 1;
      err_trasl1.insert( vv , err );
      err_rot1.insert( vv , err );
    }
    //IL piercedvector err_trasl  e err_rot1 devono essere imnizializzati 
    
    smoothing ( obj_u, idB, T_mod , n, coefficientiModificati, mapOneRing, err_trasl1, flag0, flag1 );			//Smoothing del campo traslazione con lerp, flag1 
    smoothing ( obj_u, idB, Q_mod , n, coefficientiModificati, mapOneRing, err_rot1, flag2,flag1);			//Smoothing del campo rotazione con LERP, flag1

    plotQuaternionField ("./", "QuatsfNaca_n=50", obj_u, Q_mod, T_mod);	//Plotto idue campi di quaternioni   

    //plotting quaterniorns
    bitpit::PiercedVector<double> Ort_mod = evalCellOrth(obj_u);
    dvector1D tempOrth_mod = setCounter(Ort_mod);
    plotDataQUAD("./", "Orth0Naca_err_n=50", obj_u, tempOrth_mod, "Orthogonality");

    
    //Deforming obj_u
    for (bitpit::Vertex &vv : obj_u->getVertices()){
      Quaternion P(vv.getCoords());
      darray3E vv_round = sumQuats(RotQuats(P, Q_mod[vv.getId()],COR),T_mod[vv.getId()]).getVector();
      obj_u->modifyVertex(vv_round, vv.getId());
    }
    obj_u->getPatch()->write("Deformation_err_n=50");
   
    delete objBord;
    delete obj_round;
    delete mimmo0;
    delete mimmo0_Copy;
}

int test() {
    int check = 1;
    
    return check;
    
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
	
#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling your test*/
	test1();
	
	int val = test() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
