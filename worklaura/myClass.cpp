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
\*---------------------------------------------------------------------------*/
#include "myHeader.hpp"
#include "customOperators.hpp"
#include "mimmo_iogeneric.hpp"
#include "mimmo_core.hpp"
#include "Operators.hpp"
#include "MeshSelection.hpp"
//#include "bitpit.hpp"

// #include <iostream>
// any other includes you want not already declared in the header
 using namespace std;
 using namespace bitpit;
 using namespace mimmo;
//using any other namespace
 
 //Define implementation of all your classes methods and namespace methods here.
 //Class1::ciccio(){
    //something bla bla
    //something bla bla
 //}
 /*Devo pensare alla normalizzazione
  */
 Quaternion::Quaternion()
 {
   m_scalar = 0.0;
   m_vector[0] = 0.0;
   m_vector[1] = 0.0;
   m_vector[2] = 0.0;
 }
 
/*Ad ogni punto dllo spazio si può associare un quaternione avente parte scalare nulla
 */
 Quaternion::Quaternion(const std::array<double, 3> &vector)
 {
   m_scalar = 0.0;
   m_vector = vector;
 }
 
/*Voglio creare un costruttore tale che 
 * date le parti scalari e vettoriali (angolo e assi)
 * valuti il quaternione che rappresenta la rotazione
 * Gliangoli devono essere dati in radianti
 * La norm2 è sempre =1
 */
 Quaternion::Quaternion(const double &scalar, const std::array<double, 3> &vector)
 {
     
      double n = norm2(vector);
      m_scalar = cos(scalar/2);
      if (n ==1)
	{
	for(int i=0;i<3;i++)
	  m_vector[i] = sin(scalar/2)*vector[i];
	}
      else 
	if (n > 1.0e-12)
	{
	  for(int i=0;i<3;i++)
	    m_vector[i] = sin(scalar/2)*(1/n)*vector[i];
	}
      else{
	m_vector.fill(0.0);
      }

 }
 
  bool Quaternion::operator==(const Quaternion &other)
  {
      for (int i = 0; i < 3; ++i) 
      {
	if(!utils::DoubleFloatingEqual()(m_scalar, other.m_scalar))
          if (!utils::DoubleFloatingEqual()(m_vector[i], other.m_vector[i])) {
              return false;
          }
      }
      return true;
  }
  
 Quaternion & Quaternion::operator = (const Quaternion & other){
 m_scalar = other.m_scalar;
 m_vector = other.m_vector; 
 return(*this);
  };

 void Quaternion::setScalar(const double &scalar){
   m_scalar = scalar;
 }
 
 void Quaternion::setVector(const std::array<double, 3> &vector){
   m_vector = vector;
 }
 
 double Quaternion::getScalar() const{
   return m_scalar;
 }
 
 double Quaternion::getAngle() const{
   return 2*acos(m_scalar);
 }
 
 
 std::array<double, 3> Quaternion::getVector() {
   return m_vector;
 }
  
 std::array<double, 3> Quaternion::getAxis(){
   if(std::abs(m_scalar)!=1){
     double d = sin(acos(m_scalar));
     m_vector[0] = m_vector[0]/d;
     m_vector[1] = m_vector[1]/d;
     m_vector[2] = m_vector[2]/d;
     return m_vector;
  }
  else {return m_vector;}
 } 
  
//  Quaternion::Quaternion(const Quaternion & other)
//  {
//    *this = other;
//  }

 double Quaternion::getNorm2(const Quaternion & other)  
 {
   double sum = pow(other.m_scalar,2);
   for (int i=0;i<3;i++){
      sum+=pow(other.m_vector[i],2);
  }
   return sqrt(sum);
 }

 Quaternion Quaternion::getConjugate (const Quaternion & other) const
 {
   Quaternion q;
   for(int i=0;i<3;i++)
   {
     if(other.m_vector[i]!= 0) 
       q.m_vector[i] = -other.m_vector[i];
      else
	q.m_vector[i] = 0;
   }
   q.m_scalar = other.m_scalar;
   return (q);
 }
 
 bool Quaternion::isUnitary(Quaternion & other) 
 { 
   if (other.getNorm2(other)==1)
   {
   return true;
   }
   else return false;
 }
 

 /*Devo definire le operazioni sui quaternioni
  */
 
 /* Voglio definire il logaritmo di un quaternione, Sto considerando solo quaternioni unitari
  * lnQ=[0,u*theta/2] 
  */
 
 Quaternion lnQuat (const Quaternion & x)
   {
   Quaternion z;
   z.m_scalar =0;
   if(std::abs(x.m_scalar)!=1)//se m_scalar=1 <=> theta=0 <=> sin(theta=0) =>non posso dividere per una quantità nulla 
     z.m_vector = (x.m_vector/sin(acos(x.m_scalar)))*acos(x.m_scalar);
   else//x.m_vector={0.0,0.0,0.0} Sto assegnando in caso di rotazione nulla il vettore nullo m_vector={0.0,0.0,0.0} 
     z.m_vector = x.m_vector*acos(x.m_scalar);
   return z;
   }  
   

/*Voglio scrivere una funzione chemvaluti l'esponenziale di un quaternione, Q è unitario
 * |lnQ|=theta/2 
 * u_x=(lnQ)_x/(theta/2), u_y=(lnQ)_y/(theta/2), u_z=(lnQ)_z/(theta/2)
 */

Quaternion expQuat ( Quaternion & x)
{
  double a = x.getNorm2(x);
  std::array<double,3> u;
  if(a!=0)//se la norma è diversa da zero; la norma andrà e rappresentare theta/2 con cui poi dovrò andare a dividere
    u = x.m_vector/a;
  else//Se |lnQ|==0 <=> (lnQ)_x=(lnQ)_y=(lnQ)_z=0 <=> u=x.m_vector={0.0,0.0,0.0} In effetti dovrebbe essere una forma indeterminata
    u = x.m_vector;
  
  Quaternion z(2*a,u);				//Moltiplico per 2 perchè il costruttore prende in ingresso l'angolo theta e non theta/2
  return z;  
}  

  
   Quaternion sumQuats (const Quaternion & x, const Quaternion & y)
   {
   Quaternion z;
   z.m_scalar = x.m_scalar + y.m_scalar;
   z.setVector(x.m_vector +y.m_vector);
   return z;
   }  
   
   Quaternion minQuats (const Quaternion & x, const Quaternion & y)
   {
   Quaternion z;
   z.m_scalar = x.m_scalar - y.m_scalar;
   z.setVector(x.m_vector -y.m_vector);
   return z;
   }  
   
//Prodotto scalare*quaternione
   Quaternion prodQuat (double s, const Quaternion & Q)
   {
   Quaternion sQ;
   sQ.m_scalar = s*Q.m_scalar;
   sQ.m_vector = s*Q.m_vector ;

   return sQ;
   } 
   
   
//Q_1*Q_2=x*y   
   Quaternion prodQuats (const Quaternion & Q_1, const Quaternion & Q_2)
   {
   Quaternion z;
   z.m_scalar = Q_1.m_scalar*Q_2.m_scalar-dotProduct(Q_1.m_vector, Q_2.m_vector);
  
   z.m_vector[0] = Q_1.m_scalar*Q_2.m_vector[0] + Q_2.getScalar()*Q_1.m_vector[0] + crossProduct(Q_1.m_vector,Q_2.m_vector)[0] ;
   z.m_vector[1] = Q_1.m_scalar*Q_2.m_vector[1] + Q_2.getScalar()*Q_1.m_vector[1] + crossProduct(Q_1.m_vector,Q_2.m_vector)[1];
   z.m_vector[2] = Q_1.m_scalar*Q_2.m_vector[2] + Q_2.getScalar()*Q_1.m_vector[2] + crossProduct(Q_1.m_vector,Q_2.m_vector)[2];
   return z;
   }  

   /*P_round=Q*P*Q^(-1) nel caso di quaternione unitario Q^(-1)=Q*
    *P_round=q*P*q.getConjugate(q)
    * Introduco tra i parametri il centro di rotazione 
    */
   Quaternion RotQuats (const Quaternion & P,const Quaternion & Q, darray3E & COR)//
   {
   std::array<double, 3> p ;
   p[0]= P.m_vector[0]-COR[0];
   p[1]= P.m_vector[1]-COR[1];
   p[2]= P.m_vector[2]-COR[2];
   Quaternion Pcor(p);
   std::array<double, 3> r = prodQuats(Q,prodQuats(Pcor,Q.getConjugate(Q))).getVector()+COR;//
   Quaternion P_round(r);
   return P_round;
   }  
   
   /*La funzione prende due punti ne considera il quaternione associato 
    *Date le due normali nei punti calcola asse  e angolo di rotazione e definisco il quaternione rotazione Q
    * Valuto il quaternione traslazione dalla formula X_d =QX_uQ* +T
    */
   void evalQuaternions (std::array<double, 3> &P_u, std::array<double, 3> &P_d, std::array<double, 3>& n_u, std::array<double, 3>& n_d, darray3E & COR, Quaternion &Q, Quaternion &T){
   Quaternion P_U(P_u), P_D(P_d);					//Associo ai punti i rispettivi quaternioni
   n_u /=norm2(n_u);							//Essendo mormali al vertice posso anche non normalizzarle
   n_d /=norm2(n_d);
   double a = dotProduct(n_u, n_d);
   double angle;
   if ( a > 1 && a < 1 +1.0e-8) a=1;
   if ( a < -1 && a > -1 -1.0e-8) a=-1;  
   //Voglio imporre che se il prodotto scalare è 1 allora l'angolo è zero
   if(a == 1)
     angle = 0;
   else
     angle = std::acos(a);				//angle=cos^-1(n_u*n_d)   
   std::array<double, 3> u = crossProduct(n_u, n_d);			//u = n_u x n_d
   if (norm2(u) > 1.0e-18) u /=norm2(u);				//Normalizzo solo se la norma e' diversa da zero
//   Q = Quaternion(angle,u);
   Q = Quaternion(0.0,{0.0, 0.0, 0.0});
   T = minQuats(P_D, RotQuats(P_U, Q, COR));

   }

   
/*Funzione che dato un MimmoObject restituisce un PiercedVector contenente le normali ai vertici 
 *per ogni cella chiedo la connettività e per ogni elemento della connettività calcolo la normale
 */
 bitpit::PiercedVector< std::array<double, 3> >  evalPV_VertexNormal(MimmoObject* geo){
 bitpit::SurfaceKernel * nor = static_cast <bitpit::SurfaceKernel * >(geo->getPatch());

 bitpit::PiercedVector< std::array<double, 3> > PVnor;
 for(bitpit::Cell &cell : nor->getCells()){				//Ciclo sulle celle
   int n = cell.getVertexCount();                             
   long* con = cell.getConnect();					//per ogni cella memorizzo la connettività
   for (int i=0; i<n;i++){						//ciclo sui vertici della connettività di ogni cella
     if(!PVnor.exists(con[i])){						//se l'id globale del vertice non esiste nel PiercedVector 
       PVnor.insert(con[i], nor->evalVertexNormal(cell.getId(),i));	//si inserisce un nuovo elemento (vettore normale e id associato)
    }
   }  
 }
 return PVnor; 
}

/*Voglio scrivere una funzione che prende due mimmo obj e calcola i quaternioni
 *Devo ciclare sugli id dei vertici e valutare i quats
 * se il primo mimmo ha più elementi a quei vertici assegno il quaternione q ;
 * Così che i due MimmoObject possono contenere un numero diverso di vertici.
 * Si suppone che il secondo sia contenuto nel primo 
 */
 void evalQuaternions (MimmoObject* geom0, MimmoObject* geom0_rot, bitpit::PiercedVector< std::array<double, 3> >  & pv, bitpit::PiercedVector< std::array<double, 3> >  & pv_rot, darray3E & COR, bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T){ 
//    bitpit::PiercedVector< std::array<double, 3> >   pv = evalPV_VertexNormal(geom0);			//memorizzo in due piercedvector le normali ai vertici dei MimmoObject uno rappresenta quello undeformato e il secondo deformato
//    bitpit::PiercedVector< std::array<double, 3> >   pv_rot = evalPV_VertexNormal(geom0_rot);
   bitpit::PiercedVector< bitpit::Vertex> pv_f = geom0_rot->getVertices();
   Quaternion q(0,{0.0,0.0,0.0});											//inizializzo un quaternione nullo
   for (long & i : geom0->getVertices().getIds()){							//ciclo sugli elementi del primo MimmoObject
     if (pv_f.exists(i)){								//se l'id del vertice esiste anche nel secondo mimmo calcolo i quaternioni
       Quaternion Q_rot, Q_t;
       darray3E g0 = geom0->getVertexCoords(i);
       darray3E g0_r = geom0_rot->getVertexCoords(i);
       std::array<double, 3> pv_i = pv[i];
       std::array<double, 3> pv_rot_i = pv_rot[i];
       evalQuaternions(g0, g0_r, pv_i, pv_rot_i, COR, Q_rot, Q_t);
       Q.insert(i,Q_rot);
       T.insert(i,Q_t);
     }
     else{												//altrimenti assegno dei quaternioni nulli 
       Q.insert(i, q);
       T.insert(i, q);
    }
   }
 
 }
 
// =================================================================================== //
/*Scrivo una funzione che estrae le interfacce di bordo e le salva in un PiercedVector
 */
bitpit::PiercedVector<Interface> getBorderInterfaces(PatchKernel* m_patch){
  bitpit::PiercedVector<Interface> itBorder;
  for (bitpit::Interface & it : m_patch->getInterfaces()){ 		//cicla sulle interfacce se sono di bordo le inserisce nel piercedvector itBorder
     if(it.isBorder())
       itBorder.insert(it.getId(), it);
  }
  return itBorder;
}
/*Funzione che estrae un MimmoObject objB contenente solo i bordi di obj (MimmoObject di partenza) 
 */
void extractBorders (PatchKernel* m_patch, MimmoObject* objB){
//  PatchKernel* m_patch = obj->getPatch();
  bitpit::PiercedVector<bitpit::Vertex> & objBVert = objB->getVertices();
  for (bitpit::Interface & it0 :getBorderInterfaces(m_patch)){ 		//Ciclo sulle interfacce di bordo 
    int n= it0.getConnectSize();
    long* con = it0.getConnect(); 
    livector1D cc(n);
//     std::set<long int> idcheck;						//salvo gli id dei vertici in una struttura set per non riperterli
    for(int i=0; i<n; i++){						//ciclo sulla connettività di ciascuna interfaccia e se il vertice non è ancora stato inserito lo inserisco  
      if(!objBVert.exists(con[i])){
	objB->addVertex(m_patch->getVertexCoords(con[i]), con[i]);
// 	idcheck.insert(con[i]);
      }
      cc[i] = con[i];
    }
   objB->addConnectedCell(cc, bitpit::ElementType::LINE, it0.getId());
  }
  if(!objB->areAdjacenciesBuilt()) objB->buildAdjacencies();
//   objB->getPatch()->buildInterfaces();
}

/*Funzione che dato un MimmoObject obj (contenente solo dei bordi) e considerato un'altro MimmoObject obj1 (contenente un sottoinsieme del primo MimmoObject) aggiunge a quest'ultimo le celle mancanti di obj
 *aggiungo a obj1 i bordi mancanti dal confronto con obj
 */
void completeBorders(MimmoObject* obj, MimmoObject* obj1){
  std::set<long int> ids1;
   for(bitpit::Cell &cc : obj1->getCells()){		//sto salvando in un set gli id delle celle del MimmoObject obj1 al quale inserirò altri elementi 
     ids1.insert(cc.getId());
   }   
   for (bitpit::Cell & cell : obj->getPatch()->getCells()){
//ciclo sulle celle del primo MimmoObject obj e se quell'id di cella non lo trovo recupero la connettività della cella ciclo su di essa e aggiungo i vertici e la cella 
     long int idcell = cell.getId();
     if(ids1.find(idcell)== ids1.end()){ 
       int n= cell.getConnectSize();
       long* con = cell.getConnect();
       livector1D cc(n);
       for(int i=0; i<n; i++){
	 if(!obj1->getVertices().exists(con[i])){
	  obj1->addVertex(obj->getPatch()->getVertexCoords(con[i]), con[i]);
	 }
	 cc[i] = con[i];
       }
       obj1->addConnectedCell(cc , bitpit::ElementType::LINE, idcell);
     }
   }
   if(!obj1->areAdjacenciesBuilt()) obj1->buildAdjacencies();  
}


/*Funzione che assegna ai vertici non di bordo un quaternione nullo 
 *Dai quaternioni al bordo inizializzo i quaternioni non al bordo a zero
 */
void initializeAllQuaternions (MimmoObject* obj, bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T){
  Quaternion q(0,{0.0,0.0,0.0});
  for (const long &id: obj->getPatch()->getVertices().getIds()){			//Ciclos sui vertici del MimmoObject
    if (!T.exists(id)) {								//se l'id non esiste nel PiercedVector T (contiene gli elementi con stesso id di Q ) aggiungo un quaternione di default
      T.insert(id, q);
      Q.insert(id, q);
    }
  }
}
// =================================================================================== //
// /*Scrivo una funzione che prende una patch indeformata un MimmoObject che rappresenta il bordo del corpo deformato 
//  *e due pv di quats per la rotazione e per la traslazione
//  *questi due MimmoObject li uso nella funzione
//  *void quats (MimmoObject* geom0, MimmoObject* geom0_rot,Q, bitpit::PiercedVector<Quaternion> &T)
//  */ 
//  void evalBorderQuaternions (PatchKernel* m_patch, MimmoObject* geo, bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T){  
//    MimmoObject * objB = new MimmoObject(4,2);
//    extractBorders(m_patch, objB);
//    completeBorders(objB, geo);
//    quats (objB, geo, Q, T);
// }  


// /*Creare una funzione che presa una patch scorre tutte le celle e dalla connettività percorre tutti i vertici della cella
//  *e per ogni vertice salva in una UNORDERED_MAP il vettore dei VERTICIONERINGAdiacenti
//  * si serve della funzione findVertexIdOneRing
//  *sto costruendo le interfacce 
//  * riempie la mappa per le interfaacee 
//  */
// 
std::unordered_map< long int , std::vector<long> > 
evalVertexOneRingAdj (MimmoObject* obj, 
		      std::unordered_map< long int , std::vector<long> >  &mapOneRing_interfaces)
{
  std::unordered_map< long int , std::vector<long> > mapOneRing;
  
  if(!obj->areInterfacesBuilt())	obj->buildInterfaces();
  
  int type = obj->getType();
  std::map<std::pair<long, long>, bool> visitedge;
  long id_interf;
  
  for (bitpit::Cell & cell :obj->getCells()){
    int edgecount;
    if (type == 1 || obj->getPatch()->getDimension() == 2){
      edgecount = cell.getFaceCount();
    }else{
      edgecount = cell.getEdgeCount();
    }
  
   for (int ie=0; ie<edgecount; ie++){
      bitpit::ConstProxyVector<long> econn;
      
      if (type == 1 || obj->getPatch()->getDimension() == 2){
	econn = cell.getFaceConnect(ie);
	id_interf=cell.getInterface(ie);//salvo l'id globale dell'interfaccia con indice locale ie
      }else{
	econn = cell.getEdgeConnect(ie);	
	id_interf=cell.getInterface(ie);//salvo l'id globale dell'interfaccia con indice locale ie
      }
      std::vector<long int> m_conn;
      if (!visitedge[std::pair<long, long>(econn[0],econn[1])]){	      
	visitedge[std::pair<long, long>(econn[0],econn[1])] = true;
	visitedge[std::pair<long, long>(econn[1],econn[0])] = true;
	mapOneRing[econn[0]].push_back(econn[1]);
	mapOneRing[econn[1]].push_back(econn[0]);
	mapOneRing_interfaces[econn[0]].push_back(id_interf);//Aggiungo al vettore associato al vertice econn[0] l'interfaccia con id: id_interf
	mapOneRing_interfaces[econn[1]].push_back(id_interf);//Aggiungo al vettore associato al vertice econn[1] l'interfaccia con id: id_interf
      }
    }
  }
  return mapOneRing;
}

/*Creare una funzione che presa una patch scorre tutte le celle e dalla connettività percorre tutti i vertici della cella
 *e per ogni vertice salva in una UNORDERED_MAP il vettore dei VERTICIONERING
 * si serve della funzione findVertexIdOneRing
 */
std::unordered_map< long int , std::vector<long> > evalVertexOneRing (MimmoObject* obj){
  std::unordered_map< long int , std::vector<long> > mapOneRing;
  std::set<long int> idcheck;
  PatchKernel* m_patch = obj->getPatch();
  for (bitpit::Cell &cell: obj->getCells()){
        long c = cell.getId();
        int n= cell.getConnectSize();
        long* con = cell.getConnect();
        for(long j=0; j<n; j++){
            if(idcheck.find(con[j]) == idcheck.end()){  
                std::vector<long> vOneRing = findVertexIdOneRing(c,j,m_patch);
                idcheck.insert(con[j]);
                std::pair< long int , std::vector<long> > vertOneRing (con[j], vOneRing );
                mapOneRing.insert(vertOneRing);
            }    
        }
  }
  return mapOneRing;
}


/*Voglio creare una funzione che preso l'id globale della cella, l'id locale del vertice e la geometria 
 * restituisca un vettore di id dei verticiOneRing
 */
std::vector< long > findVertexIdOneRing(const long &id, const int &vertex, PatchKernel* m_patch){

  std::vector< long > cellVertexOneRing = m_patch->findCellVertexOneRing(id,vertex);		//Salvo l'ID DELLE CELLE OneRing del vertice con id locale: vertex e id globale della cella: id
  // Ciclo sugli id del vettore findCellVertexOneRing(c,j) per ogni id della cella Ciclo sugli elementi della rispettiva connettivita
  std::vector<long int> verticesVertexOneRing;
  std::set<long> vtemp;
  long idV = m_patch->getCell(id).getConnect()[vertex];       //Connettività della cella contenente il vertice per cui vogliamo trovare il vettore dei verticiOneRing
   
  for (long &x : cellVertexOneRing){				//Ciclo sulle celle OneRing
        int n= m_patch->getCell(x).getConnectSize();
        long * v = m_patch->getCell(x).getConnect();		//Connettività delle celle OneRing
        for (int i=0; i<n; i++ ){					//Ciclo sulla connettività della cella OneRing 
            vtemp.insert(v[i]);
        }
  } 
  vtemp.erase(idV);
  verticesVertexOneRing.insert(verticesVertexOneRing.end(), vtemp.begin(), vtemp.end());
  return verticesVertexOneRing;
}


/*Voglio creare una funzione che preso l'id globale della cella e quello locale dei vertici e la geometria 
 * restituisca la distanza con tutti i vertici OneRing 
 */
std::vector<double> findDistanceOneRing(const long &id, const int &vertex, PatchKernel* m_patch){
  std::vector<double> dist;
  std::array< double, 3 > cc = m_patch->getVertex(m_patch->getCell(id).getConnect()[vertex]).getCoords();
  std::vector<long> verts = findVertexIdOneRing(id,vertex, m_patch);
  for(const long &v: verts){
    std::array< double, 3 > coords = m_patch->getVertex(v).getCoords();
    double d_i = norm2(coords-cc);
//    if(d_i > 0)
      dist.push_back(d_i);
  }
  return dist;
}
/*Voglio creare una funzione che presi due id di due vertici restituisca la distanza tra di essi
 */
double evalDistanceVertices(const long & vv, const long &x, PatchKernel* m_patch){
  double l_i = norm2(m_patch->getVertex(vv).getCoords()-m_patch->getVertex(x).getCoords());
  return l_i;
}

/*Funzione che confronta le connettività del primo e del secondo MimmoObject e se sono diverse la connettività del secondo viene modificata 
 *imponendola uguale a quella del primo MimmoObject.Inrwealtà l'andremo a cambiare e besta senza controllare e poi in caso cambiare. 
 * Uguaglia la connettivitàdel secondo MimmoObject al primo
 */
void matchConnectivity(MimmoObject * obj1, MimmoObject * obj2){
//   livector2D Conn = obj1->getConnectivity();
  PatchKernel* m_patch =obj1->getPatch();
  
  for (bitpit::Cell &cell : obj2->getCells()){			//Ciclo sulle celle del secondo MimmoObject
    int n= cell.getConnectSize();
    long* con = cell.getConnect();				//Connettività delle celle del secondo MimmoObject
    long id = cell.getId();
    bitpit::Cell c = m_patch->getCell(id);
    long* con1 = c.getConnect();	//Connettivitàdelle celle del primo mimmo 
    for(long j=0; j<n; j++){					//Ciclo sugli elementi della connettività 
      con[j] = con1[j];						//Impongo l'elemento j-esimo della connettività del secondo mimmo = all'elemento j-esimo della connettività del primo mimmo 
    }
  }
}
/*LERP di un singolo vertice con campo scalare 
 * Ho bisogno di passargli la patch per valutare la distanza con i vertici vicini 
 * l'id del vertice preso in considerazione 
 * mapOneRing[id_vertex]
 * evo immettere anche i valori del campo in corrispondenza dei vertici 
 */
double lerp(PatchKernel* m_patch, long int &id_vertex,  std::vector<long> & mapOneRing_id_vertex, bitpit::PiercedVector<double> &pv_campoScalare){
    std::vector<double> fOneRing;						//vettore che conterrà i valori del campo scalare dei verticiOneRing
    std::vector<double> lOneRing;						//vettore dove salverò tutti i reciproci delle distanze del vertice vv con i Vertici OneRing
    double s, f_new;
    for (long &x:  mapOneRing_id_vertex){
        lOneRing.push_back(std::pow(1/evalDistanceVertices(id_vertex, x, m_patch), 3.0));
      fOneRing.push_back(pv_campoScalare[x]);
    }  
    sum(lOneRing,s);									//somma delle componenti del vettore l 
    f_new = dotProduct(fOneRing,lOneRing)/s;
    return f_new;
}

 /*LERP di un singolo vertice con campo quaternioni
  * restituisce il quaternione associato al vertice id_vertex dopo l'interpolazione 
  * Per definire il quaternione Quaternion T(t) con t[i] = dotProduct(fOneRing[i],l)/s
  * sto usando il costruttore per cui impongo solo la parte vettoriale 
  * 
  * Inserisco una flag per poter fare l'interlolazione lineare con la rotazione facendo attenzione che il risultato in questo caso è il quaternione logaritmo
  * flag==0 =>lerp 
  * flag==2 =>lerp con Lie algebra
  * 
  * Devo modificere per poter aggiungere il caso in cui voglio considerare un cofficiente di diffusività non costante 
  * e voglio considerare una dumping function 
  * Quindi inserisco una variabile booleana: modify 
  *se essa è =1 uso le lughezze modificate 
  *se è =!1 uso le lunghezze e i pesi non modificati
  */
Quaternion lerp(PatchKernel* m_patch, long int &id_vertex,  std::vector<long> & ringVert, std::vector<double> &weights, bitpit::PiercedVector< Quaternion > &Q, int &flag, int &modify){//, bool modify
    std::array<double, 3> t, Q_vec;						//per memorizzare le coordinate dopo l'interpolazione
    Quaternion q;
    double s;

    bitpit::PiercedVector< double > pesi;
    std::vector<std::vector<double> > fOneRing(3 );//Tutti e 3 gli elementisono inizializzati a 0
    std::vector<double> l;
    l.reserve(ringVert.size());
    
    int locC = 0;
    for (const long & idV:  ringVert)
    {
      if(modify==1){
         l.push_back( weights[locC] /evalDistanceVertices(id_vertex, idV, m_patch) );			//Valuto le lunghezze 
      }
      else{
         l.push_back(1.0/evalDistanceVertices(id_vertex, idV, m_patch));			//Valuto le lunghezze 
      }
      
	//Se devo fare l'interpolazione di un campo quaternioni che rappresentano la traslazione salvo in Q_vec la parte vettoriale del quaternione Q[x]
	if(flag==0){
	  Q_vec = Q[idV].getVector();
	}
      
	//Se devo fare l'interpolazione di un campo quaternioni che rappresentano la rotezione usando l'algebra di Lie salvo in Q_vec la parte vettoriale del QUATERNIONE LOGARITMO lnQuat(Q[x])
	if(flag==2){
	  q = lnQuat(Q[idV]);
	  Q_vec = q.getVector();
	}
	
	for(int i=0; i<3; ++i){
	  fOneRing[i].push_back(Q_vec[i]);			//per ogni vvOneRing mi salvo il vampo vettoriale
//           fOneRing[i] = Q_vec[i];          //per ogni vvOneRing mi salvo il vampo vettoriale

	}
	++locC;
        
    }
      
      sum(l,s);
      
      for (int i=0; i< 3 ; i++){						//ciclo sulle 3 componenti e interpolo    
	t[i] = dotProduct(fOneRing[i],l)/s;
      }
      Quaternion T(t);
      
      if(flag==2){ 
	T=expQuat(T);
      }
      return T;
  }

/*SLERP di un singolo vertice 
 * restituisce il quaternione associato al vertice id_vertex dopo l'interpolazione 
 */
Quaternion slerp(PatchKernel* m_patch, long int &id_vertex,  std::vector<long> & mapOneRing_id_vertex,  bitpit::PiercedVector<Quaternion> &Q){
  Quaternion Q_slerp;//Inizializzo un primo quaternione 
//   std::random_shuffle(mapOneRing_id_vertex.begin(), mapOneRing_id_vertex.end());
  Q_slerp = Q[mapOneRing_id_vertex[0]]; 							//Assegno a Q_slerp il valore del quaternione corrispondente al primo elemento di mapOneRing_id_vertex
  double l1 = evalDistanceVertices(id_vertex, mapOneRing_id_vertex[0] , m_patch);		//Salvo la distanza del primo vertice OneRing con il vertice fissato id_vertex
  double dist = 1/l1;										//Salvo il reciproco della distanza del primo vertice OneRing con il vertice fissato id_vertex 
  for (long &x: mapOneRing_id_vertex){								//Ciclo sul vettore dei vertici OneRing del vertice  id_vertex
    double d = std::pow(evalDistanceVertices(id_vertex, x, m_patch), 1.0);  					// Valuto la distanza tra il vertice id_vertex e un vertice oneRing x
    dist +=1/d;
    double s=(1/d)/dist;
    /*Definisco l'angolo tra i due quaternioni unitari come se ogni quaternione fosse un vettore a 4 dimensioni
     * calcolo cos^-1(dotProduct(Q_slerp,Q[x])) (LE NORME di Q_slerp e Q[x] =1) 
     */
    double Omega = acos(Q_slerp.getScalar()*Q[x].getScalar()
		      +Q_slerp.getVector()[0]*Q[x].getVector()[0] 
		      +Q_slerp.getVector()[1]*Q[x].getVector()[1]
		      +Q_slerp.getVector()[2]*Q[x].getVector()[2]);
    /*Voglio calcolare Il quaternione dato dall'interpolazione sferica di Q_slerp e Q[x] tramite la formula 
     * Q_slerp = (sin((1-alpha)Omega)/sin(Omega)) * Q_1 + (sin(alpha*Omega)/sin(Omega)) * Q_2 == slerp(Q_1,Q_2,alpha). Con alpha funzione delle distanze.
     *Se Omega=0 l'angolo tra i quaternioni (visti come vettori a 4 dimensioni) è 0 => sen(0)=0 non posso dividere per 0 => Q_slerp := Q[x]
     */

    if(abs(Omega) > 1.0e-12){
      Q_slerp = sumQuats(prodQuat(sin((1-s)*Omega)/sin(Omega), Q_slerp),prodQuat(sin(s*Omega)/sin(Omega), Q[x]));
//      l_slerp2=sin((1-s)*Omega)/sin(Omega)*l1+sin(s*Omega)/sin(Omega)*d;//Per Omega piccoli il seno è approssimabile al suo angolo 
//      l_slerp2=(1-s)*l1+s*d;     
    }
    else{
      Q_slerp = Q[x];
    }
  }
 return Q_slerp;
}


/*Funzione che mi fa l'interpolazione lineare avendo un campo Scalare 
 * io ho un campo scalare e ciclando sugli id dei vertici
 * std::unordered_map ad ogni vertice associa il vettore dei verticiOneRing (obj)
 * Il vettore idB contiene gli id dei vertici di bordo (Il campo su questi punti rappresenta la C.I. quindi il valore del campo non deve essere modificato)
 */
void interpolation (MimmoObject* obj, 
                    std::vector<long int> idB,
                    std::unordered_map< long int ,
                    std::vector<long> > & mapOneRing,
                    bitpit::PiercedVector<double> &pv_campoScalare, 
                    bitpit::PiercedVector<double> &pvNew)
{
  pvNew.clear();								//Mi assicuro che il piercedvector dove andrò a memorizzare i valori del campo non contenga nulla  
  PatchKernel* m_patch = obj->getPatch();
  for (long & vv : obj->getVertices().getIds()){				//Ciclo sui vertici del MimmoObjectma voglio fare l'interpolazione solo sui vertici che non sono di bordo
    if(find(idB.begin(),idB.end(), vv)==idB.end()){				//Per tutti i vertici che non sono contenuti nel vettore idB(identifico così in tal modo i vertici interni) faccio l'interpolazione
      pvNew.insert(vv, lerp(m_patch, vv , mapOneRing[vv], pv_campoScalare));
  }
  else 
    pvNew.insert(vv, pv_campoScalare[vv]);
  }
}

/*Funzione di smoothing che preso in input un MimmoObject, un campoScalare e il numero di iterazioni che si vogliono fare chiama n volte la funzione lerp
 */
void smoothing (MimmoObject* obj,
                std::vector<long int> idB, 
                bitpit::PiercedVector<double> &pv_campoScalare,
                int &n)
{
  std::unordered_map< long int , std::vector<long> > mapOneRing_interfaces;
  std::unordered_map< long int , std::vector<long> > mapOneRing = evalVertexOneRingAdj (obj,mapOneRing_interfaces);
  for (int i=0; i<n; i++){
    bitpit::PiercedVector<double> pv_interm ;
    interpolation(obj, idB, mapOneRing, pv_campoScalare, pv_interm);
    pv_campoScalare = pv_interm;
  }
}

/*Funzione che mi fa l'interpolazione lineare dato un campo di quaternioni
 * ho messo l'if sull'essere di bordo perche' se e' di bordo il valore non deve essere modificato
 * devo specificare che tipo di interpolazione voglio fare 
 * flag==0 => lerp
 * flag==1 => slerp
 *Inserisco std::unordered_map< long int , std::vector<long> > &mapOneRing_interfaces perchèdevo identificare l'interfaccia comune a due vertici per identificare i centro cella delle celle che condividono quell'interfaccia
 * Voglio inserire 
 * flag==2 =>lerp con Lie algebra
 */
void interpolation (MimmoObject* obj,
                    std::vector<long int> idB,
                    std::unordered_map< long int ,
                    std::vector<long> > &mapOneRing,
                    std::unordered_map< long int , std::vector<double> > &mapCoeffModificati,
                    bitpit::PiercedVector< Quaternion > &Q, 
                    bitpit::PiercedVector< Quaternion > &pvNew, 
                    int &flag, int &modify)
{
//   pvNew.clear();
  PatchKernel* m_patch = obj->getPatch();
  for (long & vv : obj->getVertices().getIds()){
        if(find(idB.begin(),idB.end(), vv)==idB.end()){   //Se l'id non esiste nel vettore che contiene gli id di bordo
            //Se flag==0 allora farò un'interpolazione lineare su un campo traslazioni se flag==2 allora voglio effettuare un'interpolazione lineare su un campo rotazione usando l'algebra di Lie
            if (flag==0 || flag==2){
                pvNew[vv] = lerp(m_patch, vv , mapOneRing[vv], mapCoeffModificati[vv] , Q, flag, modify);/////////////////////////////////////////////////////////////////////////////////////
            }
            //Se flag==1 farò una slerp
            if(flag==1){
                pvNew[vv] = slerp(m_patch, vv , mapOneRing[vv], Q);
            } 
        }
    else{ //Se l'id è di bordo non voglio in quei punti modificare il valore del campo 
        pvNew[vv] = Q[vv];
    }
  }
} 

/*Funzione che valuti icoefficienti per i pesi modificati 
 * se voglio i pesi modificati
 */
std::unordered_map< long int , std::vector<double> >  mapCoeffModificati (MimmoObject* obj,MimmoObject * objB,  std::unordered_map <long int , std::vector<long>> & mapOneRing, std::unordered_map <long int , std::vector<long>> &mapOneRing_interfaces, int &modify ){  

  std::unordered_map< long int , std::vector<double> > mapCoeffMod; 
  bitpit::PiercedVector< double > eta = evalEta( obj, objB );
  bitpit::PiercedVector< double > tau = evalTau( obj, modify ); 

  for ( long & id_vertex : obj->getVertices().getIds() ){
    double tau_i = tau[id_vertex];
    double eta_i = eta[id_vertex];
    PatchKernel* m_patch = obj->getPatch();
    
    bitpit::PiercedVector< double > pesi;
    std::vector<std::vector<double> > fOneRing(3);
    
    std::vector<long> it = mapOneRing_interfaces[id_vertex];//salvo in vettore di long gli id delle interfacce relative al vertice id_vertex
    std::vector<long> vertOneRing = mapOneRing[id_vertex];//salvo in un vettore di long gli id dei vertici oneRing di id_vertex
    
    for (long & x : vertOneRing) {	//Ciclo sul vettore contenete gli id dei vertici contenuti nel OneRing
      
      if(modify==1){			//Se voglio calcolare le lunghezze modificate
	std::vector<long> m_x = mapOneRing_interfaces[x];//vettore contenente le interfacce del OneRing del vertice x appartenente al oneRingdi id_vertex
	double tau_ij;
	double Aij;
	std::array< long, 2 > cell_ON;
	tau_ij=(tau_i+tau[x])/2;//Calcolo la media dei due valori dei vertici
      
	/*Devo calcolare l'area di interfaccia A_ij
	* Devo passare dall'edgecount
	* Per ogni coppia di vertici devo passare ai centri cella delle celle che contengono entrambi i vertici tramite getOwnerNeigh mi faccio dare le celle
	* poi dovrò ricavare i centri cella e calcolare la distanza tra i due centri cella e quella distanza sarà A_ij
	*/
	int counter=0;
	for ( long & id_it : it ){
	  for ( long &i : m_x ) { //elementi delle interfacce collegate al nodo xhe appartiene al oneRing di id_vertex
	    if ( id_it == i ) {
           if(m_patch->getInterface(i).isBorder()){
                bitpit::VolUnstructured * temp = static_cast<VolUnstructured*>(m_patch);
                Aij = temp->evalCellVolume(temp->getInterface(i).getOwner())/temp->evalInterfaceArea(i);           
	       }else{ 
                cell_ON = m_patch->getInterface(i).getOwnerNeigh();//Sto salvando in questo vettore gli id delle due celle che hanno in comune questa interfaccia condivisa da entrambi i vettori mapOneRing_interfaces
                Aij = norm2( m_patch->evalCellCentroid(cell_ON[0])-m_patch->evalCellCentroid(cell_ON[1]) );//Calcolo la distanza tra i due centri cella(trovati tramite evalCellCentroid)  
	      }	
	    }
	  }
	}//chiude il ciclo sulle interfacceOneRing
        //mapCoeffMod[id_vertex].push_back( std::pow ( tau_ij*Aij, eta_i ) );			//Valuto i coefficienti per generare poi i pesi modificati
        mapCoeffMod[id_vertex].push_back( std::pow ( tau_ij*Aij, 1.0 ) );         //Valuto i coefficienti per generare poi i pesi modificati
          
    }//chiude if modify
      
    }//chiude il ciclo sui vertOneRing
  }//chiude il ciclo sui vertici
  return mapCoeffMod;
}//chiude la funzione


/*Funzione che valuti icoefficienti per i pesi modificati utilizzando la funzione di dumping optimad. 
 */
std::unordered_map< long int , std::vector<double> >  
mapCoeffModDumping (MimmoObject* obj,
                    MimmoObject * objB,  
                    std::unordered_map <long int , std::vector<long>> & mapOneRing, 
                    double soglia, 
                    double expfactor,
                    int &modify ){  

  std::unordered_map< long int , std::vector<double> > mapCoeffMod; 
  bitpit::PiercedVector< double > dump = evalDumpOptimad( obj, objB, soglia, expfactor );

  for ( long & id_vertex : obj->getVertices().getIds() ){
    PatchKernel* m_patch = obj->getPatch();
    std::vector<long> vertOneRing = mapOneRing[id_vertex];//salvo in un vettore di long gli id dei vertici oneRing di id_vertex
    for (long & x : vertOneRing) {  //Ciclo sul vettore contenete gli id dei vertici contenuti nel OneRing
        if(modify==1){            //Se voglio calcolare le lunghezze modificate
            mapCoeffMod[id_vertex].push_back( dump[x] );         //Valuto i coefficienti per generare poi i pesi modificati
        }else{
            mapCoeffMod[id_vertex].push_back( 1.0 );         //Valuto i coefficienti per generare poi i pesi modificati
        }
    }//chiude il ciclo sui vertOneRing
  }//chiude il ciclo sui vertici
  return mapCoeffMod;
}//chiude la funzione
      
      
/*Funzione che calcola gli eta  per ogni vertice cioè quell'esponente da utilizzare nella formula del tau (vedi evalTau) 
 * passo puntatore a mesh intera obj ed puntatore al suo bordo completo objB
 */
bitpit::PiercedVector< double > evalEta (MimmoObject* obj, MimmoObject * objB ){  
  //La distanza che si deve calcolare la voglio calcolare dall'oggetto composto solo dalle celle (elementi linea in questo caso) che si deformano
  //Devo quindi ricavare un nuovo mimmoobject che contenga solo le celle di bordo che si stanno deformando
  bitpit::PiercedVector< double > dist;						//Definisco un PiercedVector in cui memorizzo per ogni vertice la distanza che questo ha dal corpo che si deforma 
  bitpit::PiercedVector< double > eta; 
  
  //MimmoObject * objB = new MimmoObject(4,2);
  //PatchKernel* m_patch_d = obj->getPatch();
  //extractBorders (m_patch_d, objB);						//Estraendo ora i bordi del corpo deformato
 
  //Metodo per estrazione della sola parte del  bordo che mi interessa (profilo)
  SelectionByBox * selector = new SelectionByBox();
  selector->setOrigin( {{ 0.5, 0.0, 0.0}} );
  selector->setSpan( {{ 2.0, 0.2 , 0.6}} );
  selector->setGeometry(objB);
  selector->execute();
  
  MimmoObject * objB_Def = selector->getPatch(); 
//   objB_Def->getPatch()->write("provaObj_def");
//   exit(1);
  
  if(!objB_Def->areAdjacenciesBuilt())  objB_Def->buildAdjacencies();

  if(!objB_Def->isBvTreeSync()){
      objB_Def->buildBvTree();
    }  
    BvTree * tree = objB_Def->getBvTree();
    
    long id;
    std::array < double, 3 > coord; 
    double d_max = 200;
    //Sto ciclando sui vertici di obj per calcolare la distanza dal bordo deformato per calcolare la dumping function
    //Devo salvare in un piercedvector per ogni vertice il valore della distanza
    for ( const bitpit::Vertex & vv : obj->getVertices()){
      double r = 1.E+18;
      coord = vv.getCoords();
      double d = bvTreeUtils::distance( &coord, tree, id, r);
      //Fisso una distanza massima dopo la quale voglio che le lunghezze modifcate coincidano con quelle normali 
      //e normalizzo per quella distanza per le lunghezze maggiori la dumping function sarà settata a zero
      dist.insert( vv.getId(), d/d_max );
    }//chiude il for sui vertici per il calcolo delle distanze
     //Date le distanze devo calcolare in ogni vertice il valore di eta
     double soglia = 5e-05;
     for ( long & id : obj->getVertices().getIds() ){
       if ( dist[ id ] < soglia )
	 eta.insert( id , 1 );
       else
	 if ( dist[ id ] >= soglia && dist[ id ] <= 1){
	   double Eta = std::pow( soglia ,2)/std::pow( (dist[ id ] ) , 2 );
	   eta.insert( id , Eta);
	}
	 else
	  if(dist[id]>1)
	    eta.insert( id , 0 );
    }	

    delete selector;
    
        dvector1D disttemp = setCounter(dist);
     obj->getPatch()->getVTK().addData("dist", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, disttemp);
     obj->getPatch()->write("plotdist");
     obj->getPatch()->getVTK().removeData("dist");


     dvector1D etatemp = setCounter(eta);
     obj->getPatch()->getVTK().addData("eta", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, etatemp);
     obj->getPatch()->write("ploteta");
     obj->getPatch()->getVTK().removeData("eta");
    
    return (eta); 
}

/*Funzione di dumping Optimad sulla mesh obj che utilizza la distanza di level set di ogni punto di mesh dal corpo objB 
 * passo puntatore a mesh intera obj ed puntatore al suo bordo completo objB
 * soglia ed expfactor sono due parametri da definire da me. vedi appunti da Edoardo.
 */
bitpit::PiercedVector< double > evalDumpOptimad (MimmoObject* obj, MimmoObject * objB, double soglia, double expfactor ){  
  //La distanza che si deve calcolare la voglio calcolare dall'oggetto composto solo dalle celle (elementi linea in questo caso) che si deformano
  //Devo quindi ricavare un nuovo mimmoobject che contenga solo le celle di bordo che si stanno deformando
  bitpit::PiercedVector< double > dist;                     //Definisco un PiercedVector in cui memorizzo per ogni vertice la distanza che questo ha dal corpo che si deforma 
  bitpit::PiercedVector< double > dumping; 
  
  //Metodo per estrazione della sola parte del  bordo che mi interessa (profilo)
  SelectionByBox * selector = new SelectionByBox();
  selector->setOrigin( {{ 0.5, 0.0, 0.0}} );
  selector->setSpan( {{ 2.0, 0.2 , 0.6}} );
  selector->setGeometry(objB);
  selector->execute();
  
  MimmoObject * objB_Def = selector->getPatch(); 
  if(!objB_Def->areAdjacenciesBuilt())  objB_Def->buildAdjacencies();

  if(!objB_Def->isBvTreeSync()){
      objB_Def->buildBvTree();
  }  
  BvTree * tree = objB_Def->getBvTree();
    
    long id;
    std::array < double, 3 > coord; 
    //Sto ciclando sui vertici di obj per calcolare la distanza dal bordo deformato per calcolare la dumping function
    //Devo salvare in un piercedvector per ogni vertice il valore della distanza
    for ( const bitpit::Vertex & vv : obj->getVertices()){
        double r = 1.E+18;
        coord = vv.getCoords();
        double d = bvTreeUtils::distance( &coord, tree, id, r);
        dist.insert( vv.getId(), d);
    }
    
     for ( long & id : obj->getVertices().getIds() ){
        double distnorm = soglia/std::max(1.0e-8, dist[id]); 
        double Eta = std::max(1.0, std::pow(distnorm, expfactor));
        dumping.insert( id , Eta);
     }
    
    delete selector;
    
//       dvector1D disttemp = setCounter(dist);
//      obj->getPatch()->getVTK().addData("dist", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, disttemp);
//      obj->getPatch()->write("plotdist");
//      obj->getPatch()->getVTK().removeData("dist");
// 
// 
//      dvector1D dumptemp = setCounter(dumping);
//      obj->getPatch()->getVTK().addData("optimad dumping", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, dumptemp);
//      obj->getPatch()->write("plotdumpingOptimad");
//      obj->getPatch()->getVTK().removeData("optimad dumping");
    
    return (dumping); 
}



bitpit::PiercedVector< double > evalTau (MimmoObject* obj,  int &modify ){  
  
  bitpit::PiercedVector< double > tau;   //Definisco un PiercedVector in cui memorizzo per ogni vertice il tau associato
  bitpit::PiercedVector< double > tau_cell;   //Definisco un PiercedVector in cui memorizzo per ogni cella il tau associato ad essa
  bitpit::PiercedVector< double > volume_cell;   //Definisco un PiercedVector in cui memorizzo per ogni cella il volume associato ad essa

 if (modify==1){
    bitpit::VolUnstructured * vol = static_cast <bitpit::VolUnstructured * >(obj->getPatch());
    for(bitpit::Cell &cell : vol->getCells())
    {
      double volume = vol->evalCellVolume(cell.getId());
      volume_cell.insert(cell.getId(), volume);
    }
    double massimo, minimo;
    evalMinMaxPV(volume_cell, minimo, massimo);//Per calcolare il massimo e minimo globale 
     
    /*Fino a questo momento ho calcolato questo fattore di diffusività nella cella devo ora assegnare tale valore ai vertici 
     * Voglio ciclare sui vertici e per ogni vertice con il OneRing andare a mediare sui valori di tau_cell associati alle celle del OneRing
     */
    std::set<long int> idcheck;//Per non rivalutare il tau nello stesso vertice 
    std::vector< long > OneRing;
  
    for(bitpit::Cell &cell : obj->getCells()){//ciclo sulle celle per definire tau
        long c = cell.getId();
        int n= cell.getConnectSize();
        long* con = cell.getConnect();
        for(long j=0; j<n; j++){                  //ciclo sulla connettività delle celle
            if(idcheck.count(con[j]) < 1){      //Se non si è mai passati dal vertice gli chiedo il one ring di celle e ciclo su esso
                OneRing = obj->getPatch()->findCellVertexOneRing(c, j);
                double avgvolume = 0.0;
                double Tau = 0.0;
                for (long &id : OneRing ){                //Ciclo sugli id dele celle OneRing del vertice con[j]
                    avgvolume+=volume_cell[id];                 //Voglio fare una media dei valori tau delle celle quindi sommo tutti i valori che tau assume in ogni cella 
                }
                idcheck.insert(con[j]);             //Inserisco l'id del vertice nel vettore i modo da valutare una sola volta il valore di tau[con[j]]     

                avgvolume /= double(OneRing.size());                  //Devo dividere per il numero di elementi contenuti nel OneRing
                Tau = 1.0+ ( massimo - minimo ) / avgvolume;
                tau.insert(con[j],Tau);
            }
        }
    }//chiudo il cclo sulle celle per definire tau
  
  }else{
        for (long & vv : obj->getVertices().getIds()){ 
            tau.insert( vv, 1);
        }
  }  
/*
  dvector1D tautemp = setCounter(tau);
  obj->getPatch()->getVTK().addData("tau", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, tautemp);
  obj->getPatch()->write("plottau");
  obj->getPatch()->getVTK().removeData("tau");
 */
  return (tau); 
}//chiude evalTau

/*Funzione che preso in input un MimmoObject, un campo Quaternion, il numero di iterazioni che si vogliono fare e il tipo di interpolazione che si vuole fare 
 * chiama n volte la funzione lerp;
 * PRIMA FLAG =>TIPO DI INTERPOLAZIONE che voglio fare:
 * flag==0 => lerp, flag==1 => slerp, flag==2 => lie algebra
 * 
 * Voglio introdurre un fattore di diffusività locale,
 * SECONDA FLAG =>PESI MODIFICATI O NO 
 * La variabile modify permette di scegliere se usare i pesi modificati o meno
 * modify!=1 lunghezze non modificate 
 * modify==1 lunnghezze modificate
 * Salvo l'errore a ogni iterazione in un pierced vector che dò io vuoto 
 * 
 * iL PIERCED VECTOR DEVE ESSERE GIA' RIEMPITO
 */
void smoothing (MimmoObject* obj,
                std::vector<long int> idB, 
                bitpit::PiercedVector< Quaternion > &Q, 
                int &n, 
                std::unordered_map< long int , std::vector<double> > &mapCoeffModificati,
                std::unordered_map <long int , std::vector<long>> & mapOneRing,
                bitpit::PiercedVector< double > & err_iter, 
                int &flag,  
                int &modify) 
{

  double errmax = 1.E+18;
  int itMAX = 200000;
  double err_threshold = 1.0E-20;
  int itLOC = 0;
  
  //inizializzo un piercedvector dove vado a salvare la norma dei quaternioni Q_old che sono contenuti nel piercedvector Q
  bitpit::PiercedVector< double > norm_old ;  
  bitpit::PiercedVector< Quaternion > Q_new ;
  bitpit::PiercedVector< double > N ;
  double min, max, min_old , max_old ;
  {
    for ( long & i : Q.getIds() )
    {
        double compT = Q[i].getNorm2(Q[i]);
        norm_old.insert( i, compT );
        N.insert( i, 0);
    } 
  } 
      evalMinMaxPV( norm_old , min_old, max_old );

   
  for ( long & vv :obj->getVertices().getIds() ){
      Quaternion x;
      Q_new.insert( vv , x );
    }
    
  while (itLOC < itMAX && errmax >= err_threshold){
    interpolation(obj, idB, mapOneRing, mapCoeffModificati , Q, Q_new, flag, modify);
    //Per il calcolo dell'errore devo valutare max(Q_new-Q_old)/max(Q_old)= max ( Q_new - Q ) / max( Q )
    //Ciclo sugli id e valuto il quaternione differenza 
    for ( long & i : Q.getIds() ){
      Quaternion q = minQuats ( Q_new[i], Q[i] );
      double nor = q.getNorm2(q);
      N[i] =  nor ;
    }
    //Valuto l'errore servendomi del massimo e minimo delle norme dei piercedvector delle norme dei quaternioni vecchi e dei quaternioni differenza
    //Valuto il massimo e il minimo del piercedvector contenente le norme dei quaternioni Q_old
      

    //Valuto il massimo e il minimo del piercedvector contenente le norme dell a differenza tra Q_new - Q_old
    evalMinMaxPV( N, min, max ); 
    errmax = max/max_old;
    //Inseisco nel PiercedVector all'id che rappresenta il numero dell'iterazione il valore max associato all'iterazione n-esima
    err_iter[ itLOC ] = errmax ; 
    //Assegno al piercedvector norm_old i valori aggiornati
//     for ( long & i : Q_new.getIds() ){
//         double compT = Q_new[i].getNorm2( Q_new[i] );
//         norm_old[i] = compT;
//     } 
    Q = Q_new;
    std::cout << "iterazione " << itLOC <<  "     max " << max << "    errmax  "  << errmax << std::endl;
    itLOC++;
  }
  std::cout<<"ALL DONE SMOOTHING"<<std::endl; 
}

 /*Voglio scrivere una funzione che modifica il centro di rotazione
  * potrei mettere una flag per definire se da {0.0,0.0,0.0} passo a COR o viceversa
  * flag==0 se dall'origine voglio passare ad un'altro COR
  * flag==1 se da unCOR generico voglio passare all'origine
  * NON LA STO USANDO
  */
void modifyCOR (MimmoObject* obj, darray3E & COR, int flag){
  for (bitpit::Vertex &vv : obj->getVertices())
  { 
    darray3E vv_new ;
    if(flag==0)//Impongo che il vertice COR corrisponda con l'origine modificando le coordinate di tutti i punti darray3E
       vv_new =vv.getCoords()-COR;
    if(flag==1)//Per ripotare tutto l'obj alle condizioni iniziali quindi per ogni vertice modifico le coordinate aggiungendo COR così da essere l'inversa del caso flag==0darray3E 
      vv_new =vv.getCoords()+COR;
    
    obj->modifyVertex(vv_new, vv.getId());
  }
}
/*Funzione che plotta il campo scalare rappresentante l'angolo della rotazione e il campo vettoriale rappresentante la traslazione 
 * in input si deve dare la stringa che indica la directory , stringa che indica il nome che verà assegnato al file , il MimmoObject e i due PiercedVector di quaternioni 
 */

void plotQuaternionField (std::string dir, std::string name, MimmoObject * geom0_rot, bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T){
  dvector1D temp_Q(Q.size()); 
  dvecarr3E temp_T(T.size()); 
  {
    int counter = 0;
    for(bitpit::Vertex v : geom0_rot->getPatch()->getVertices()){
      temp_Q[counter] = Q[v.getId()].getAngle();
      temp_T[counter] = T[v.getId()].getVector();
      ++counter;
    }
  }
  
  ///////bitpit::VTKElementType::LINE   
  bitpit::VTKElementType cellType = bitpit::VTKElementType::QUAD;
  
  bitpit::VTKUnstructuredGrid output(dir,name,cellType);
  dvecarr3E vert = geom0_rot->getVertexCoords();
  ivector2D connc = geom0_rot->getCompactConnectivity();
  output.setGeomData( bitpit::VTKUnstructuredField::POINTS, vert );
  output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connc );
  output.setDimensions(connc.size(), vert.size());
  output.addData("AngleQuaternion" , bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, temp_Q);
  output.addData("Translation" , bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, temp_T);
  output.setCodex(bitpit::VTKFormat::ASCII);
  output.write(); 
}


/*Funzione che plotta un campo vettoriale rappresentante la traslazione 
 * in input si deve dare la stringa che indica la directory , stringa che indica il nome che verà assegnato al file , il MimmoObject
 * e un PiercedVector<std::array<double,3>>rappreenta il campo vettoriale 
 */

void plotTranslationField (std::string dir, std::string name, MimmoObject * obj, bitpit::PiercedVector< std::array<double, 3> > &campo_vettoriale){
  dvecarr3E temp_T(campo_vettoriale.size());
  {
    int counter = 0;
    for(long v : campo_vettoriale.getIds()){
      temp_T[counter] = campo_vettoriale[v];
      ++counter;
    }
  }
  bitpit::VTKElementType cellType = bitpit::VTKElementType::QUAD;
  
  bitpit::VTKUnstructuredGrid output(dir,name,cellType);
  dvecarr3E vert = obj->getVertexCoords();
  ivector2D connc = obj->getCompactConnectivity();
  output.setGeomData( bitpit::VTKUnstructuredField::POINTS, vert );
  output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connc );
  output.setDimensions(connc.size(), vert.size());
  output.addData("Translation" , bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, temp_T);
  output.setCodex(bitpit::VTKFormat::APPENDED);
  output.write(); 
}

/*Voglio creare una funzione che mi restituisca un pierced vector con l'ortogonalità 
 *che prenda un MimmoObject
 * bitpit::PiercedVector<double>
 */
bitpit::PiercedVector<double> evalCellOrth (MimmoObject * geom ){
    bitpit::PiercedVector< Interface >::iterator itInterfaces; 
    bitpit::PiercedVector<double>::iterator itInterfaces1;
    bitpit::PiercedVector<double> PVCosQuadAlpha;
    bitpit::PiercedVector<double> PVOrt;
    double cosquad;
    darray3E FirstVec;
    darray3E SecondVec;
    double Orth;
    PatchKernel* m_patch = geom->getPatch();
    if(!geom->areAdjacenciesBuilt())    geom->buildAdjacencies();
    if(!geom->areInterfacesBuilt())     geom->buildInterfaces(); 
    bitpit::VolUnstructured * nor = static_cast <bitpit::VolUnstructured * >(m_patch);

   /*Ciclo sulle interfacce calcolo i due vettori e il coseno quadro formato dai due vettori
    *Tengo in considerazione delle interfacce al bordo imponendo CosQuadAlpha=2
    */
    for(itInterfaces = m_patch->interfaceBegin();   itInterfaces != m_patch->interfaceEnd();    itInterfaces++){
        if (itInterfaces->isBorder()) {
            cosquad = 2;  
            PVCosQuadAlpha.insert(itInterfaces.getId(), cosquad);
        }else{
            FirstVec = m_patch->evalCellCentroid(itInterfaces->getOwnerNeigh()[0])-m_patch->evalCellCentroid(itInterfaces->getOwnerNeigh()[1]);
            SecondVec = nor->evalInterfaceNormal(itInterfaces.getId());
            cosquad = pow( (dotProduct(FirstVec,SecondVec)), 2.0)/ (pow(norm2(FirstVec), 2.0)* pow(norm2(SecondVec), 2.0) );  
            PVCosQuadAlpha.insert(itInterfaces.getId(), cosquad);    
        } 
    }

    for(bitpit::Cell &cell : m_patch->getCells()){
      int dim = cell.getInterfaceCount();
      long * CellInter = cell.getInterfaces();
      double minimo = PVCosQuadAlpha[CellInter[0]]; 
      
      for(int j=1; j<dim;j++){
            minimo = std::min(minimo, PVCosQuadAlpha[cell.getInterfaces()[j]]);
      };
      Orth = minimo;
      if (Orth > 1.0) {
        Orth = 2;
      }
      PVOrt.insert(cell.getId(), Orth);
    }
    
    return  PVOrt;
}

/*Creo una funzione che calcoli Ort_ratio= Orth_deformed/Orth_undeformed = Orth_obj1/Orth_obj0
 */
bitpit::PiercedVector<double> evalOrth_ratio (MimmoObject * obj0 ,MimmoObject * obj1){
  bitpit::PiercedVector<double> pvOrth_ratio;
  bitpit::PiercedVector<double> temp_0 = evalCellOrth(obj0);
  bitpit::PiercedVector<double> temp_1 = evalCellOrth(obj1);
  for(bitpit::Cell &cell : obj0->getPatch()->getCells()){ 
    double a = temp_1[cell.getId()]/temp_0[cell.getId()];
    pvOrth_ratio.insert(cell.getId() , a);
  }

   return(pvOrth_ratio);
}



dvector1D setCounter (bitpit::PiercedVector<double> pv){
   dvector1D temp_pv(pv.size()); 
  {
    int counter = 0;
    for(const double & val : pv){
      temp_pv[counter] = val;
      ++counter;
    }
  }
return(temp_pv);
}


/*Dovrei in plotDataQUAD discriminare il caso in cui VTKLocation::POINT oppure VTKLocation::CELL
 */
/*
void plotDataLine (std::string dir, std::string name, MimmoObject * obj, dvector1D temp_pv, std::string NameProp ){
 
 bitpit::VTKElementType cellType = bitpit::VTKElementType::LINE;
 
 bitpit::VTKUnstructuredGrid output(dir,name,cellType);
 dvecarr3E vert = obj->getVertexCoords();
 ivector2D connc = obj->getCompactConnectivity();
 output.setGeomData( bitpit::VTKUnstructuredField::POINTS, vert );
 output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connc );
 output.setDimensions(connc.size(), vert.size());
 output.addData(NameProp , bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, temp_pv);
 output.setCodex(bitpit::VTKFormat::APPENDED);
 output.write();
 }
 */

void plotDataQUAD (std::string dir, std::string name, MimmoObject * obj, dvector1D temp_pv, std::string NameProp ){
 
 bitpit::VTKElementType cellType = bitpit::VTKElementType::QUAD;
 
 bitpit::VTKUnstructuredGrid output(dir,name,cellType);
 dvecarr3E vert = obj->getVertexCoords();
 ivector2D connc = obj->getCompactConnectivity();
 output.setGeomData( bitpit::VTKUnstructuredField::POINTS, vert );
 output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connc );
 output.setDimensions(connc.size(), vert.size());
 output.addData(NameProp , bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, temp_pv);
 output.setCodex(bitpit::VTKFormat::APPENDED);
 output.write();
 }
 
/*Voglio scrivere una funzione che calcoli il MINIMO di un PiercedVector
 */
 void evalMinMaxPV (bitpit::PiercedVector<double>  &PV, double &minimo, double &massimo){
   
     double a = pow (10.0,10.0);
     double b = 0;
     for(double &val : PV){
       a = min(a,val);
       b = max(b,val);
     }
     minimo= a;
     massimo=b;
 }
 
 /*Funzione che prende un MimmoObject che rappresenta la griglia non deformata e un MimmoObject che rappresenta solo il corpo deformato
  *valuta tutti i quaternioni
  * Lenormali non sono giuste non la sto usando la funzione
  *devo specificare anche il centro di rotazione 
  */
 void getQuaternionC_I(MimmoObject * obj, MimmoObject * obj_r, darray3E & COR , bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T){
   MimmoObject * objB = new MimmoObject(4,2);
   
   if(!obj->areAdjacenciesBuilt())  obj->buildAdjacencies();					//costruisco le adiacenze  
   if(!obj->areInterfacesBuilt())   obj->buildInterfaces();					//costruisco le interfacce
//     obj_r->buildAdjacencies();
//     obj_r->getPatch()->buildInterfaces();					//costruisco le interfacce

   PatchKernel* m_patch = obj->getPatch();
   extractBorders (m_patch, objB);					//Estraendo ora i bordi
   completeBorders(objB, obj_r);						//Aggiungo ai bordi di obj2 i bordi mancanti di objB
//      objB->buildAdjacencies();
//      objB->getPatch()->buildInterfaces();					//costruisco le interfacce

   matchConnectivity(objB, obj_r);					//Mi assicuro che la connettività sia uguale   
   
   bitpit::PiercedVector< std::array<double,3> > pv;
   bitpit::PiercedVector< std::array<double,3> > pv_rot;
   
   evalQuaternions( objB, obj_r, pv, pv_rot, COR, Q, T);				//Valuto i quaternioni secondo il centro di rotazione specificato

   initializeAllQuaternions(obj, Q, T);
   delete objB;
 }

 /*Funzione per calcolare le normale su una giometria le cui celle sono delle linee
  */
/*
 bitpit::PiercedVector< std::array<double,3> > evalNormal (MimmoObject * obj){
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
//Calcolo le normali ai vertici facendo una media sulle normali alle interfacce
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
    
    
    }//endfunction
 
 */
