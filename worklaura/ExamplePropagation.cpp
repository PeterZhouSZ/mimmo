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

/*Nel seguente test si definisce prima le geometrie e i campi vettoriali 
 * Test che prende due MimmoObject: obj e obj2; il sucondo è una sottoporzione di bordo del primo che viene ruotata con q1((90*M_PI/180),b)
 *definisco un campo scalare e due campi vettoriali (uno di questi calcolando i quaternioni e considerando la parte vettoriale) 
 *uso lo smoothing e plotto il campo vettoriale 
 */
void test1(){
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
 

  //Sto definendo sia obj2 che obj_round per poi poter rappresentare con writer0 la geometria prima e dopo la deformazione 
  PatchKernel* m_patch = obj->getPatch();
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
	obj_round->addVertex(m_patch->getVertexCoords(con[i]), con[i]);
	idcheck.insert(con[i]);
      }
      cc[i] = con[i];
    }
   obj_round->addConnectedCell(cc, bitpit::ElementType::LINE, it.getId());//Aggiungo l'elemento con connettività cc
  }
  obj_round->buildAdjacencies();
  
  
   
  int flag0 = 0;
  int flag1 = 1;
  int flag2 = 2;

  bitpit::PiercedVector<double>  pv_campoScalare;		//Inizializzo un campo scalare: assegnerò agli elementi del bordo interno e agli elementi del bordo di coordinate [x,0] un valore !=0 altrinmenti verrà assegnato 0
  bitpit::PiercedVector< std::array<double, 3> > campo_vett;	//Inizializzo un campo vettoriale: assegnerò agli elementi del bordo interno {{1.0,1.0,0.0}} a tutti i restanti vertici assegnerò  {{0.0,0.0,0.0}}
  bitpit::PiercedVector< std::array<double, 3> > campo_vettorialeT;	//Inizializzo un campo vettoriale:dove salverò la parte vettoriale dei quaternioni traslazione 
  bitpit::PiercedVector<Quaternion> Q;
  bitpit::PiercedVector<Quaternion>  K;
  bitpit::PiercedVector<Quaternion>  Tt;
  std::array<double, 3> b{0.0,0.0,8.0};					//Asse di rotazione serve per modificare l'obj2 ruotandolo attorno a b e con q1
  Quaternion q1((30*M_PI/180),b);	//per una prova era stato messo come angolo 0 =>rotazione nulla; 
  std::array<double, 3> COR{10.5,10.5,0.0};					//Voglio imporre come centro di rotazione il pt centrale del quadrato  
//   std::array<double, 3> O{0.0,0.0,0.0};					//Voglio imporre come centro di rotazione il pt centrale del quadrato  

  
  MimmoObject * obj2 = new MimmoObject(4,2); 				 //obj2 contiene solo elementi di bordo
  obj2->setHARDCopy(obj_round);
/*Definisco un campo scalare e un campo vettoriale in ogni punto del bordo interno
 */    
  for (bitpit::Vertex &vv : obj2->getVertices())			//Ciclo sui vertici di obj2
    {
      Quaternion P(vv.getCoords());					//ogni punto può essere visto come un quaternione con parte scalare nulla 
        std::array<double, 3> b{0.0,0.0,0.0};
      Quaternion t(b);
      darray3E vv_mod2 = RotQuats(P, q1,COR).getVector();			//ruoto ciascun vertice P con il quaternione q1
      obj2->modifyVertex(vv_mod2, vv.getId());				//modifico le coordinate dei vertici con  i valori trovati

      darray3E vv_mod3 = vv.getCoords() + darray3E({0.5, 0.0, 0.0});
      obj2->modifyVertex(vv_mod3, vv.getId());				//modifico le coordinate dei vertici dopo aver introdotto una traslazione
      
      
      if (!pv_campoScalare.exists(vv.getId()))
	  pv_campoScalare.insert(vv.getId(), 30);			//Definisco il campo scalare e vettoriale in obj2
	  campo_vett.insert(vv.getId(), {{1.0,10.0,0.0}});
      Tt.insert(vv.getId(), t);        
    }  
    
     Quaternion  q0,  Q_2((0*M_PI/180),b), Q_slerp, Q_3((270*M_PI/180)*2,b), Q_4((180*M_PI/180)*2,b),qq2((360*M_PI/180),b);	//Inizializzo due quaternioni con cui voglio andare a fare l'interpolazione e un quaternione vuoto dove andremo a memorizzare il quaternione dopo aver fatto l'interpolazione sferica 
     
     
   for (const long &id: obj->getPatch()->getVertices().getIds()){	//Ciclo sui vertici del obj generale e se l'id non esiste gli assegna un valore 
    if (!pv_campoScalare.exists(id)) { 
      if (id <8)
	pv_campoScalare.insert(id, 30);
      else
	pv_campoScalare.insert(id, 0);
    }
    if(!campo_vett.exists(id))
      campo_vett.insert(id, {{0.0,0.0,0.0}});
    
    if(!Tt.exists(id))
      Tt.insert(id, q0);
  }

  MimmoObject * objB = new MimmoObject(4,2);
  extractBorders (m_patch, objB);					//Estraendo ora i bordi
  completeBorders(objB, obj2);						//Aggiungo ai bordi di obj2 i bordi mancanti di objB
  matchConnectivity(objB, obj2);					//Mi assicuro che la connettività sia uguale   
  
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
/////////////////////////////    //HO COMMENTATO TUTTO PERCHè STO INSERENDO NELLA FUNZIONE EVAL QUATERNIONS I PIERCED VECTOR CONTENENTI LE NORMALI ///////////
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  
//   evalQuaternions(objB,obj2,COR,Q,K);
//   initializeAllQuaternions(obj, Q, K);
// 
//   
// //   MimmoGeometry * writer0 = new MimmoGeometry();
// //   writer0->setIOMode(1);
// //   writer0->setDir(".");
// //   writer0->setFilename("ProvaProp_Laura.0000");
// //   writer0->setFileType(FileType::CURVEVTU);
// //   writer0->setGeometry(obj_round);
// //   writer0->execute();  
// //   
// //   writer0->setIOMode(1);
// //   writer0->setDir(".");
// //   writer0->setFilename("ProvaProp_Laura.0001");
// //   writer0->setFileType(FileType::CURVEVTU);
// //   writer0->setGeometry(obj2);
// //   writer0->execute();
// 
//   std::vector<long int> idB = obj2->getVertices().getIds();		//vector che contiene solo id di obj2 che coincide con gli id dei vertici di bordo
// 
//   int n =10;						//n numero di iterazioni per lo smoothing
//   smoothing(obj, idB, K, n, flag0,flag0);			//Smoothing del campo traslazione con lerp, flag0
//   smoothing(obj, idB, Q, n, flag2,flag0);			//Smoothing del campo rotazione con LERP, flag0
// //   smoothing(obj, idB, Q, n, flag1);			//Smoothing del campo rotazione con slerp
//   
// //   plotQuaternionField ("./", "Quatsfield L-S_n=10", obj, Q, K);	//Plotto idue campi di quaternioni   
//   plotQuaternionField ("./", "QuatsField(R+T)[L-L](n=10)", obj, Q, K);	//Plotto i due campi quaternioni   
//      
//     
//   MimmoObject * obj_r = new MimmoObject(2,2);
//   obj_r->setHARDCopy(obj);
//   for (bitpit::Vertex &vv : obj_r->getVertices())
//     {
//       Quaternion P(vv.getCoords());
//       darray3E vv_round = sumQuats(RotQuats(P, Q[vv.getId()],COR),K[vv.getId()]).getVector();
//       obj_r->modifyVertex(vv_round, vv.getId());
//     }
//     obj_r->getPatch()->write("Deformation(R+T)[L-L](n=10)");
// //      obj_r->getPatch()->write("DeformationL-S_n=10");
// 
//      bitpit::PiercedVector<double> Ort_def = evalCellOrth(obj_r);
//      dvector1D tempOrth_def = setCounter(Ort_def);
//      plotDataQUAD("./", "Orthogonality(R+T)[L-L](n=10)", obj_r, tempOrth_def, "Orthogonality");

  
  delete objB;
  delete obj;
//   delete obj_r; l'ho commentato perche in eval quaternions sto aggiungendo i pierced vector con le normali 
  delete obj2;
  delete obj_round;
//   delete writer0;
  
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
