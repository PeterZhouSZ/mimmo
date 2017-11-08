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

using namespace std;
using namespace bitpit;
using namespace mimmo;

// =================================================================================== //
/*Date due configurazioni diverse di uno stesso oggetto(deformato e non)
 *trovare i quaternioni che rappresentano questa deformazione.
 */
void test1(){
 
  MimmoObject * geom0 = new MimmoObject(4,2);
  MimmoObject * geom0_rot = new MimmoObject(4,2);
  
  std::array<double, 3> COR{0.0,0.0,0.0};					//Voglio imporre come centro di rotazione il pt centrale del quadrato  

  dvecarr3E vertices(4, {{0,0,0}});
  dvecarr3E verticesROT(4, {{0,0,0}});
  livector2D conns(4, livector1D(2,0));
  vertices[1][0] = 1.0;
  vertices[2][0] = 1.0; vertices[2][1] = 1.0;
  vertices[3][1] = 1.0;
  
  conns[0] = {{0,1}};
  conns[1] = {{1,2}};
  conns[2] = {{2,3}};
  conns[3] = {{3,0}};

  verticesROT[0][0] =  0.0;		 verticesROT[0][1] = 0.0;
  verticesROT[1][0] =  0.5;		 verticesROT[1][1] = 0.866025; 
  verticesROT[2][0] =  -0.366025;	 verticesROT[2][1] = 1.36603;
  verticesROT[3][0] =  -0.866025; 	 verticesROT[3][1] = 0.5;
  //Aggiungo a geom0 (obj_undeformed) vertici e la connettività
  for(darray3E & vv : vertices){
    geom0->addVertex(vv);
  }
  long int i=0;
  for(livector1D & cc : conns){    
    geom0->addConnectedCell(cc, bitpit::ElementType::LINE,  i);
    i++;
  } 
  //Aggiungo a geom0_rot (obj_deformed) vertici e la connettività (sarà la stessa del MimmoObject non-deformato)
  for(darray3E & vv : verticesROT){
    geom0_rot->addVertex(vv);
  }
  long int j=0;
  for(livector1D  & cc : conns)
  {    
    geom0_rot->addConnectedCell(cc, bitpit::ElementType::LINE,  j);
    j++;
  } 
  geom0->getPatch()->buildAdjacencies();
  geom0->getPatch()->buildInterfaces();
  
  geom0_rot->getPatch()->buildAdjacencies();
  geom0_rot->getPatch()->buildInterfaces();
  
  MimmoGeometry * writer = new MimmoGeometry();
  writer->setIOMode(1);
  writer->setDir(".");
  writer->setFilename("planarquads_Laura.0000");
  writer->setFileType(FileType::CURVEVTU);
  writer->setGeometry(geom0);
  writer->execute();  
  
  writer->setIOMode(1);
  writer->setDir(".");
  writer->setFilename("planarquads_Laura.0001");
  writer->setFileType(FileType::CURVEVTU);
  writer->setGeometry(geom0_rot);
  writer->execute();

  //Calcolo i quaternioni 
  bitpit::PiercedVector<Quaternion> Q;
  bitpit::PiercedVector<Quaternion> T; 
   
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
/////////////////////////////    //HO COMMENTATO TUTTO PERCHè STO INSERENDO NELLA FUNZIONE EVAL QUATERNIONS I PIERCED VECTOR CONTENENTI LE NORMALI ///////////
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  
//   evalQuaternions(geom0,geom0_rot, COR, Q, T);
  
  delete geom0;
  delete geom0_rot;
  delete writer;
  
}

// =================================================================================== //
/*Esempio di più rotazioni con vari quaternioni 
 */
void test2(){
  //Definisco una griglia 7*7 ed elimino una cella interna 
  double h=21, l=21;
  int Nh= 7, Nl= 7;
  double dx = h/Nh;
  double dy = l/Nl;
  std::vector<std::array<double,3>> coord ((Nh+1)*(Nl+1),{0.0,0.0,0.0});
  for(int j=0; j<=Nl; j++){
    for(int i=0; i<=Nh;i++){
      int nG=i+(Nh+1)*j;
      coord[nG][0] = i*dx;
      coord[nG][1] = j*dx;
    }
  }
  livector2D con(Nh*Nl, livector1D(4,0));
  for ( int jc=0; jc<Nl; jc++){
    for(int ic=0; ic<Nh; ic++){
      int nC1 = ic+Nh*jc;
      con[nC1][0] = ic+(Nh+1)*jc;
      con[nC1][1] = (ic+1)+(Nh+1)*jc;
      con[nC1][2] = (ic+1)+(Nh+1)*(jc+1);
      con[nC1][3] = ic+(Nh+1)*(jc+1);
    }
  }
  MimmoObject * obj = new MimmoObject(2,2);
  for(darray3E &vv : coord){
    obj->addVertex(vv);
  }
  long int i=0;
  for(livector1D & cc : con){    
    obj->addConnectedCell(cc, bitpit::ElementType::QUAD, i);
    i++;
  } 
  obj->getPatch()->buildAdjacencies(); 
  obj->getPatch()->deleteCell(24);
  obj->getPatch()->buildInterfaces();
  
//    Prova di modifyCOR
//    obj->getPatch()->write("Modify_000");
//    darray3E vv_0 = {{10.5, 10.5, 0.0}};
//    modifyCOR(obj, vv_0, 0) ;
//    obj->getPatch()->write("Modify_001");
//    modifyCOR(obj, vv_0, 1);  
//    obj->getPatch()->write("Modify_002");
  std::array<double, 3> COR{10.5,10.5,0.0};
  std::array<double, 3> O{0.0, 0.0,0.0};
  
  bitpit::PiercedVector<bitpit::Interface>  it_0 = getBorderInterfaces(obj->getPatch());	//Estraggo le interfacce di bordo
  bitpit::PiercedVector<bitpit::Interface>  it_1;						//Salverò le interfacce relative al quadrato eliminato
  for (bitpit::Interface &i : it_0){					//Ciclo sulle interfacce di bordo ed estraggo quelle della cella che ho eliminato e le salvo nel piercedvector it_1  
    int n= i.getConnectSize();
    long int* con = i.getConnect();
    bool check = false;
    int count = 0;
    for(int j=0; j<n; j++){
      std::array<double,3> a ;
      a = obj->getPatch()->getVertexCoords(con[j]);
      if( a[0]>6 && a[0]<15 && a[1]>6 && a[1]<15) 			//Per COSTRUZIONE di griglia so che le coordinate dei vertici cercate appartengono a questi intervalli
	count++;
	if (count == n){ 						//L'interfaccia viene inserita nel piercedvector se tutti i vertici di essa appartengono agli intervalli specificati (n=size della connettività) 
	  it_1.insert(i.getId(), i);
	}
    }
  }

  PatchKernel* m_patch = obj->getPatch();
  MimmoObject * obj1 = new MimmoObject(4,2);
  dvecarr3E vertices(4, {{0,0,0}});
  livector2D c(4, livector1D(2,0));
  for (bitpit::Interface & it :it_1){ 					//Ciclo sulle interfacce del quadrato che ho eliminato
    int n= it.getConnectSize();
    long* con = it.getConnect(); 
    livector1D cc(n);
    std::set<long int> idcheck;						//salvo gli id dei vertici in una struttura set per non riperterli
    for(int i=0; i<n; i++){						//ciclo sulla connettività di ciascuna interfaccia e se il vertice non è ancora stato inserito lo inserisco  
      if(idcheck.find(con[i])== idcheck.end()){
	obj1->addVertex(m_patch->getVertexCoords(con[i]), con[i]);
	idcheck.insert(con[i]);
      }
      cc[i] = con[i];
    }
   obj1->addConnectedCell(cc, bitpit::ElementType::LINE, it.getId());
  }
  obj1->buildAdjacencies();
  
  std::array<double, 3> b{0.0,0.0,1.0};					//Definisco l'asse di rotezione b sarà uguale per ogni rotazione nel seguito
  Quaternion q((90*M_PI)/180,b),  qu0;
     
/*Faccio una copia del MimmoObject obj1 e modifico le coordinate dei vertici assegnando dei nuovi valori alle coordinate ottenendo obj2
 */
  MimmoObject * obj2 = new MimmoObject(4,2);
  obj2->setHARDCopy(obj1);
  for (bitpit::Vertex &vv : obj2->getVertices())
    {
      darray3E v_0 = {{9.0, 9.0, 0.0}};
      darray3E v_1 = {{12.0, 9.0, 0.0}};
      darray3E v_2 = {12.0, 12.0, 0.0};
      darray3E v_3 = {9.0, 12.0, 0.0};
      if (vv.getCoords() == v_0){
	darray3E vv_mod2 = {{-9.0, 9.0, 0.0}};
	obj2->modifyVertex(vv_mod2, vv.getId());
      }
      if (vv.getCoords() == v_1){
	darray3E vv_mod2 = {{-9.0, 12.0, 0.0}};
	obj2->modifyVertex(vv_mod2, vv.getId());
      }
      if (vv.getCoords() == v_2){
	darray3E vv_mod2 = {{-12.0, 12.0, 0.0}};
	obj2->modifyVertex(vv_mod2, vv.getId());
      }
      if (vv.getCoords() == v_3){
	darray3E vv_mod2 = {{-12.0, 9.0, 0.0}};
	obj2->modifyVertex(vv_mod2, vv.getId());
      }
     }  
      

/*Faccio un'altra copia di obj1 modificando le coordinate ruotando i punti con il quaternione q1((M_PI),b) ottenendo obj3
 */
  std::array<double, 3> COR3{-10.5,10.5,0.0};

  MimmoObject * obj3 = new MimmoObject(4,2);
  Quaternion q1((M_PI),b);
  obj3->setHARDCopy(obj1);
  for (bitpit::Vertex &vv : obj3->getVertices())
    {
      Quaternion P(vv.getCoords());
      darray3E vv_mod2 = RotQuats(P, q1,COR).getVector();
      obj3->modifyVertex(vv_mod2, vv.getId());
    }  
    
/*Faccio una copia dell'ultimo MimmoObject: obj3 e modifico le coordinate dei vertici dell'ultimo MimmoObject con il quaternione  q4((60*M_PI)/180,b) ottenendo obj4
 */   

  MimmoObject * obj4 = new MimmoObject(4,2);
  Quaternion q4((60*M_PI)/180,b);
  obj4->setHARDCopy(obj3);
  for (bitpit::Vertex &vv : obj4->getVertices())
    {
      Quaternion P(vv.getCoords());
      darray3E vv_mod2 = RotQuats(P, q4,COR).getVector();
      obj4->modifyVertex(vv_mod2, vv.getId());
    }
    
/*Faccio una copia dell'ultimo MimmoObject: obj4 e vado a memorizzare i vertici ruotati con q5((120*M_PI)/180,b) ottenendo obj5
 */   
  MimmoObject * obj5 = new MimmoObject(4,2);
  Quaternion q5((120*M_PI)/180,b);
  obj5->setHARDCopy(obj4);
  for (bitpit::Vertex &vv : obj5->getVertices())
    { 
      darray3E v_0 = {{9.0, 9.0, 0.0}};
      Quaternion P(vv.getCoords());
      darray3E vv_mod2 = RotQuats(P, q5,COR).getVector()+v_0;
      obj5->modifyVertex(vv_mod2, vv.getId());
    }
  
  MimmoGeometry * writer0 = new MimmoGeometry();
  writer0->setIOMode(1);
  writer0->setDir(".");
  writer0->setFilename("ProvaQuats_Laura.0000");
  writer0->setFileType(FileType::CURVEVTU);
  writer0->setGeometry(obj1);
  writer0->execute();  
  
  writer0->setIOMode(1);
  writer0->setDir(".");
  writer0->setFilename("ProvaQuats_Laura.0001");
  writer0->setFileType(FileType::CURVEVTU);
  writer0->setGeometry(obj2);
  writer0->execute();
  
  writer0->setIOMode(1);
  writer0->setDir(".");
  writer0->setFilename("ProvaQuats_Laura.0002");
  writer0->setFileType(FileType::CURVEVTU);
  writer0->setGeometry(obj3);
  writer0->execute();
  
  writer0->setIOMode(1);
  writer0->setDir(".");
  writer0->setFilename("ProvaQuats_Laura.0003");
  writer0->setFileType(FileType::CURVEVTU);
  writer0->setGeometry(obj4);
  writer0->execute();
  
  writer0->setIOMode(1);
  writer0->setDir(".");
  writer0->setFilename("ProvaQuats_Laura.0004");
  writer0->setFileType(FileType::CURVEVTU);
  writer0->setGeometry(obj5);
  writer0->execute();
  
/*Calcolo tra gli ultimi due MimmoObject i quaternioni rotazione e traslazione e plotto il campo
 */
  bitpit::PiercedVector<Quaternion> Q;
  bitpit::PiercedVector<Quaternion> T;

   
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
/////////////////////////////    //HO COMMENTATO TUTTO PERCHè STO INSERENDO NELLA FUNZIONE EVAL QUATERNIONS I PIERCED VECTOR CONTENENTI LE NORMALI ///////////
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  
//   evalQuaternions(obj4, obj5, O, Q, T);
//   plotQuaternionField("./", "Test2-ExampleQuaternions", obj5, Q, T);

  Quaternion  q0, Q_2((90*M_PI/180)*2,b), Q_slerp, Q_3((270*M_PI/180)*2,b), Q_4((180*M_PI/180)*2,b),qq2((90*M_PI/180),b), Q_7((60*M_PI/180),b);	//Inizializzo due quaternioni con cui voglio andare a fare l'interpolazione e un quaternione vuoto dove andremo a memorizzare il quaternione dopo aver fatto l'interpolazione sferica 
 
  delete obj;
  delete obj1;
  delete obj2;
  delete obj3;
  delete obj4;
  delete obj5;
  delete writer0;
  
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
// 	test1();
	test2();
	
	int val = test() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
