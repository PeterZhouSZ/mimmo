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

//#include <unordered_map>
using namespace std;
using namespace bitpit;
using namespace mimmo;

// =================================================================================== //
void test1(){
 MimmoObject * obj = new MimmoObject(4,2);
  MimmoObject * obj_rot = new MimmoObject(4,2);
  
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

  for(darray3E & vv : vertices){
    obj->addVertex(vv);
  }
  long int i=0;
  for(livector1D & cc : conns){    
    obj->addConnectedCell(cc, bitpit::ElementType::LINE,  i);
    i++;
  } 

  obj->getPatch()->buildAdjacencies();
  obj->getPatch()->buildInterfaces();
  
  obj_rot->setHARDCopy(obj);
  std::array<double, 3> O{0.0,0.0,0.0};
  std::array<double, 3> b{0.0,0.0,1.0};					//Definisco l'asse di rotezione b sarà uguale per ogni rotazione nel seguito
  Quaternion q1((30*M_PI/180),b);
  for (bitpit::Vertex &vv : obj_rot->getVertices())
    {
      Quaternion P(vv.getCoords());
      darray3E vv_mod2 = RotQuats(P, q1, O).getVector();
      obj_rot->modifyVertex(vv_mod2, vv.getId());
    }  
    
  MimmoGeometry * writer = new MimmoGeometry();
  writer->setIOMode(1);
  writer->setDir(".");
  writer->setFilename("DiffCoR_Laura.0000");
  writer->setFileType(FileType::CURVEVTU);
  writer->setGeometry(obj);
  writer->execute();  
  
  writer->setIOMode(1);
  writer->setDir(".");
  writer->setFilename("DiffCoR_Laura.0001");
  writer->setFileType(FileType::CURVEVTU);
  writer->setGeometry(obj_rot);
  writer->execute();

  int flag0 = 0;
  int flag1 = 1;
  int flag2 = 2;
  
  darray3E vv_0 = {{21, 21, 0.0}};
  std::array<double, 3> COR{10.5,10.5,0.0};					//Voglio imporre come centro di rotazione il pt centrale del quadrato  
  
  //Calcolo i quaternioni 
  bitpit::PiercedVector<Quaternion> Q;
  bitpit::PiercedVector<Quaternion> T; 
  bitpit::PiercedVector< std::array<double, 3> >   pv = evalPV_VertexNormal(obj);			//memorizzo in due piercedvector le normali ai vertici dei MimmoObject uno rappresenta quello undeformato e il secondo deformato
  bitpit::PiercedVector< std::array<double, 3> >   pv_rot = evalPV_VertexNormal(obj_rot);
  evalQuaternions(obj,obj_rot, pv, pv_rot, COR, Q, T);

  delete obj;
  delete obj_rot;
  delete writer;
  }

void test2(){
  /*Genero una griglia 7x7 con (x,y) appartiene a [0,21]x[0,21]
   * aggiungo i vertici e la connettività ad un MimmoObject obj estraggo una sottoporzione di bordo obj2
   */
  double h=21, l=21;
  int Nh= 51, Nl= 51;
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
  obj->getPatch()->buildAdjacencies();					//costruisco le adiacenze  
  obj->getPatch()->deleteCell(1300);					//Elimino la cella 24 che coincide in questo specifico caso a quella centrale
  obj->getPatch()->buildInterfaces();					//costruisco le interfacce
  
  bitpit::PiercedVector<bitpit::Interface>  it_0 = getBorderInterfaces(obj->getPatch());	//Estraggo le interfacce di bordo del MimmoObject
  bitpit::PiercedVector<bitpit::Interface>  it_1;			//Sto inizializzando questo pierced vector per memorizzare le interfacce della cella che ho eliminato
  for (bitpit::Interface &i : it_0){					//Ciclo sulle interfacce di bordo
    int n= i.getConnectSize();						//Restituisce la size della connettivitàdell'elemento
    long int* con = i.getConnect();					//E' un puntatore alla connettività di un elemento 
    bool check = false;
    int count = 0;			//Inizializzo un contatore per poter aggiungere l'elemento(interfaccia) se tutti i suoi vertici (elementi della connettività) rispettano contemporaneamente dei limiti che andrò ad imporre 
    for(long j=0; j<n; j++){						//Ciclo sulla connettività
      std::array<double,3> a ;						//salvo in a le coordinatedel vertice di posto j della connettività
      a = obj->getPatch()->getVertexCoords(con[j]);
      if( a[0]>6 && a[0]<15 && a[1]>6 && a[1]<15)			//Impongo ch esia x che y appartengono all'insieme [6,15]
	count++;
	if (count == n){						//Se ciascun vertice che appartiene all'interfaccia rispetta queste limitazioni
	  it_1.insert(i.getId(), i);					//Allora l'interfaccia viene inserita
	}
    }
  }   
 
  PatchKernel* m_patch = obj->getPatch();
  MimmoObject * obj2 = new MimmoObject(4,2);			//Inizializzo un nuovo MimmoObject che rappresenterà la porzione di bordo ruotato
  dvecarr3E vertices(4, {{0,0,0}});
  livector2D c(4, livector1D(2,0));
  for (bitpit::Interface & it :it_1){ 					//Ciclo sul PiercedVector che contiene le interfacce del quadrato che ho eliminato
    int n= it.getConnectSize();
    long* con = it.getConnect(); 
    livector1D cc(n);							//Vector cove andrò a salvare la connettività di ciascuna interfaccia
    std::set<long int> idcheck;						//salvo gli id dei vertici in una struttura set per non riperterli
    for(int i=0; i<n; i++){						//ciclo sulla connettività di ciascuna interfaccia e se il vertice non è ancora stato inserito lo inserisco  
      if(idcheck.find(con[i])== idcheck.end()){
	obj2->addVertex(m_patch->getVertexCoords(con[i]), con[i]);
	idcheck.insert(con[i]);
      }
      cc[i] = con[i];
    }
   obj2->addConnectedCell(cc, bitpit::ElementType::LINE, it.getId());//Aggiungo l'elemento con connettività cc
  }
  obj2->buildAdjacencies();
  
  int flag0 = 0;
  int flag1 = 1;
  int flag2 = 2; 
  std::array<double, 3> COR{10.5,10.5,0.0};					//Voglio imporre come centro di rotazione il pt centrale del quadrato  
  darray3E COR0 = {{0.0, 0.0, 0.0}}; 
  darray3E COR1 = {{3.0, 3.0, 0.0}};
  darray3E COR2 = {{9.0, 9.0, 0.0}};
  darray3E COR3 = {{0.0, 10.5, 0.0}};
  darray3E COR4 = {{6.0, 6.0, 0.0}};


  std::array<double, 3> b{0.0,0.0,1.0};					//Asse di rotazione serve per modificare l'obj2 ruotandolo attorno a b e con q1
  Quaternion q1((30*M_PI/180),b);	//per una prova era stato messo come angolo 0 =>rotazione nulla; 

/*Definisco un campo scalare e un campo vettoriale in ogni punto del bordo interno
 */    
  for (bitpit::Vertex &vv : obj2->getVertices())			//Ciclo sui vertici di obj2
    {
      std::array<double, 3> COR{10.5,10.5,0.0};					//Voglio imporre come centro di rotazione il pt centrale del quadrato  
      Quaternion P(vv.getCoords());					//ogni punto può essere visto come un quaternione con parte scalare nulla 
      std::array<double, 3> b{0.0,0.0,0.0};
      Quaternion t(b);
      darray3E vv_mod2 = RotQuats(P, q1,COR).getVector();			//ruoto ciascun vertice P con il quaternione q1 COME COR STO SCEGLIENDO IL CENTRO DEL QUADRATINO{10.5,10.5,0.0}
      obj2->modifyVertex(vv_mod2, vv.getId());				//modifico le coordinate dei vertici con  i valori trovati
      
     darray3E vv_mod3 = vv.getCoords() + darray3E({0.5, 0.0, 0.0});
     obj2->modifyVertex(vv_mod3, vv.getId());				//modifico le coordinate dei vertici dopo aver introdotto una traslazione
    }       
      
      
      bitpit::PiercedVector<Quaternion> Q, Q_mod;
      bitpit::PiercedVector<Quaternion>  K, K_mod;

    getQuaternionC_I(obj, obj2, COR, Q, K);

//     Voglio partire dalla stessa condizione iniziale e effettuuare due smoothing diversi
    getQuaternionC_I(obj, obj2, COR, Q_mod, K_mod);

//     plotQuaternionField ("./", "QuatsFieldBorderCOR", obj, Q, K);	//Plotto idue campi di quaternioni inizializzati su obj

  std::vector<long int> idB = obj2->getVertices().getIds();		//vector che contiene solo id di obj2 che coincide con gli id dei vertici di bordo
// Ho commentato PERCHè STO INSERENDO NELLO SMOOTHING IL bVtREE
  /*
  //Sto prendendo in considerazione solo il caso in cui la rotazione la tratto con l'algebra di Lie e non più con l'interpolazione sferica
  int n =2;						//n numero di iterazioni per lo smoothing
  bitpit::PiercedVector< double > err_trasl0, err_trasl1;
  bitpit::PiercedVector< double > err_rot0, err_rot1;
  
  smoothing(obj, idB, K, n, err_trasl0, flag0,flag0);			//Smoothing del campo traslazione con lerp,flag0
  smoothing(obj, idB, Q, n, err_rot0, flag2,flag0);			//Smoothing del campo rotazione con Lie ,flag0 

  smoothing(obj, idB, K_mod, n, err_trasl1, flag0,flag1);			//Smoothing del campo traslazione con lerp,flag1 PESI MODIFICATI
  smoothing(obj, idB, Q_mod, n, err_rot1, flag2,flag1);			//Smoothing del campo rotazione con Lie ,flag1 PESI MODIFICATI 

  bitpit::PiercedVector<double> Delta_def;
  
//   plotQuaternionField ("./", "QuatsFieldCORmodifyNO(n=2)", obj, Q, K);	//Plotto i due campi quaternioni   
     
  MimmoObject * obj_r = new MimmoObject(2,2);
  obj_r->setHARDCopy(obj);
  for (bitpit::Vertex &vv : obj_r->getVertices())
    {
      Quaternion P(vv.getCoords());
      darray3E vv_round = sumQuats(RotQuats(P, Q[vv.getId()],COR),K[vv.getId()]).getVector();
      obj_r->modifyVertex(vv_round, vv.getId());
    }   
    PatchKernel* m_patch_r = obj_r->getPatch();
    m_patch_r->write("DeformationCORmodifyNO(n=2)");

//   const std::array< double, 3 > & bitpit::PatchKernel::getVertexCoords 
       
  MimmoObject * obj_rMOD = new MimmoObject(2,2);
  obj_rMOD->setHARDCopy(obj);
  for (bitpit::Vertex &vv : obj_rMOD->getVertices())
    {
      Quaternion P(vv.getCoords());
      darray3E vv_round = sumQuats(RotQuats(P, Q_mod[vv.getId()],COR),K_mod[vv.getId()]).getVector();
      obj_rMOD->modifyVertex(vv_round, vv.getId());
      double scarto=norm2(vv_round-m_patch_r->getVertexCoords(vv.getId()));
      Delta_def.insert(vv.getId(), scarto);

    }
    obj_rMOD->getPatch()->write("DeformationCORmodify(n=2)");
  
  dvector1D tempDelta = setCounter(Delta_def);
  bitpit::VTKElementType cellType = bitpit::VTKElementType::QUAD;
  
  bitpit::VTKUnstructuredGrid output("./","Delta_def",cellType);
  dvecarr3E vert = obj->getVertexCoords();
  ivector2D connc = obj->getCompactConnectivity();
  output.setGeomData( bitpit::VTKUnstructuredField::POINTS, vert );
  output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connc );
  output.setDimensions(connc.size(), vert.size());
  output.addData("Delta_def" , bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, tempDelta);
  output.setCodex(bitpit::VTKFormat::APPENDED);
  output.write();
  
//      bitpit::PiercedVector<double> Ort_def = evalCellOrth(obj_r);
//      dvector1D tempOrth_def = setCounter(Ort_def);
//      plotDataQUAD("./", "OrthogonalityCORmodifyNO(n=2)", obj_r, tempOrth_def, "Orthogonality");
//   
//   delete objB;*/
  delete obj;
  delete obj2;
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
