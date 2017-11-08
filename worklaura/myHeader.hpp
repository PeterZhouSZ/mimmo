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
#ifndef __MYHEADER_HPP__
#define __MYHEADER_HPP__

#include "bitpit.hpp"
//#include "MimmoNamespace.hpp"

#include "mimmo_iogeneric.hpp"
#include "mimmo_core.hpp"
#include "Operators.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;
//#include "BaseManipulation.hpp"
// any other include you want
//define prototypes of your classes, namespaces, enum etc.. here 


class Quaternion {
protected:
  double m_scalar;
  std::array<double, 3> m_vector;
  double m_magnitude;
public:
  Quaternion();
  Quaternion(const std::array<double, 3> &vector);
  Quaternion(const double &scalar, const std::array<double, 3> &vector);
  
  Quaternion(const Quaternion & other) = default;
  Quaternion(Quaternion &&other) = default;
  Quaternion & operator = (const Quaternion & other);
  Quaternion & operator = (Quaternion &&other) = default;
  bool operator == (const Quaternion &other);
  void setVector(const std::array<double, 3> &vector);    
  void setScalar(const double &scalar);
  
  double getScalar() const;
  double getAngle() const;
  std::array<double, 3> getVector();
  std::array<double, 3> getAxis();
  double getNorm2 (const Quaternion & other);
  Quaternion getConjugate ( const Quaternion & other) const;
  bool isUnitary (Quaternion & other) ; 
  friend Quaternion sumQuats (const Quaternion & x, const Quaternion & y);
  friend Quaternion minQuats (const Quaternion & x, const Quaternion & y);
  friend Quaternion lnQuat (const Quaternion & x);
  friend Quaternion expQuat ( Quaternion & x);
  
  friend Quaternion prodQuat (double s, const Quaternion & Q);
  friend Quaternion prodQuats (const Quaternion & Q_1, const Quaternion & Q_2);
  friend Quaternion RotQuats (const Quaternion & P, const Quaternion & Q, darray3E & COR);
};

   bitpit::PiercedVector< std::array<double, 3> > evalDeformedNormals(MimmoObject*obj_undef, MimmoObject* obj_def);//, bitpit::PiercedVector< std::array<double, 3> > & normal_undef
   
   bitpit::PiercedVector< std::array<double, 3> >  evalPV_VertexNormal(MimmoObject* geo);
   void evalQuaternions (std::array<double, 3> &P_u, std::array<double, 3> &P_d, std::array<double, 3>& n_u, std::array<double, 3> &n_d, darray3E & COR , Quaternion &Q, Quaternion &T);
//    void evalQuaternions (MimmoObject* geom0, MimmoObject* geom0_rot, darray3E & COR, bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T);
   void evalQuaternions (MimmoObject* geom0, MimmoObject* geom0_rot, bitpit::PiercedVector< std::array<double, 3> >  & pv, bitpit::PiercedVector< std::array<double, 3> >  & pv_rot, darray3E & COR, bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T);

   
   bitpit::PiercedVector<Interface> getBorderInterfaces(PatchKernel* m_patch) ;  
   void extractBorders (PatchKernel* m_patch, MimmoObject* objB);
   void completeBorders(MimmoObject* obj, MimmoObject* obj1);
   void initializeAllQuaternions (MimmoObject* obj, bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T);
   
   std::vector< long > findVertexIdOneRing(const long &id, const int &vertex, PatchKernel* m_patch);
   std::vector<double> findDistanceOneRing(const long &id, const int &vertex, PatchKernel* m_patch);
   std::unordered_map< long int , std::vector<long> > evalVertexOneRingAdj (MimmoObject* obj, std::unordered_map< long int , std::vector<long> >  &mapOneRing_interfaces);
   std::unordered_map< long int , std::vector<long> > evalVertexOneRing (MimmoObject* obj);
   double evalDistanceVertices(const long & vv, const long &x, PatchKernel* m_patch);
   void matchConnectivity(MimmoObject * obj1, MimmoObject * obj2);

   double lerp(PatchKernel* m_patch, long int &id_vertex,  std::vector<long> & mapOneRing_id_vertex, bitpit::PiercedVector<double> &pv_campoScalare);
//    Quaternion lerp(PatchKernel* m_patch, long int &id_vertex,  std::vector<long> & mapOneRing_id_vertex,  bitpit::PiercedVector< double > &tau, double &eta_id, bitpit::PiercedVector< Quaternion > &Q, int &flag, int &modify);
   Quaternion lerp(PatchKernel* m_patch, long int &id_vertex,  std::vector<long> & ringVert, std::vector<double> & weights, bitpit::PiercedVector< Quaternion > &Q, int &flag, int &modify);

   
   Quaternion slerp(PatchKernel* m_patch, long int &id_vertex,  std::vector<long> & mapOneRing_id_vertex,  bitpit::PiercedVector<Quaternion> &Q);
   
   bitpit::PiercedVector< double > evalEta (MimmoObject* obj, MimmoObject * objB);
   bitpit::PiercedVector< double > evalTau (MimmoObject* obj,  int &modify );
   
   bitpit::PiercedVector< double > evalDumpOptimad (MimmoObject* obj, MimmoObject * objB, double soglia, double expfactor);
   
   std::unordered_map< long int , std::vector<double> > mapCoeffModificati (MimmoObject* obj,MimmoObject * objB,  std::unordered_map <long int , std::vector<long>> & mapOneRing, std::unordered_map <long int , std::vector<long> > &mapOneRing_interfaces, int &modify );
   
   std::unordered_map< long int , std::vector<double> > mapCoeffModDumping (MimmoObject* obj,
                                                                            MimmoObject * objB,  
                                                                            std::unordered_map <long int , std::vector<long>> & mapOneRing,
                                                                            double soglia, 
                                                                            double expfactor,
                                                                            int &modify );   
   
   
   void interpolation (MimmoObject* obj, std::vector<long int> idB, std::unordered_map< long int , std::vector<long> > & mapOneRing, bitpit::PiercedVector<double> &pv_campoScalare, bitpit::PiercedVector<double> &pvNew);
//    void interpolation (MimmoObject* obj, std::vector<long int> idB, std::unordered_map< long int , std::vector<long> > &mapOneRing, std::unordered_map< long int, std::pair<long, long>  > &mapOneRing_interfaces,  bitpit::PiercedVector< double > &tau,  bitpit::PiercedVector< double > &eta,  bitpit::PiercedVector<Quaternion> &Q,  bitpit::PiercedVector<Quaternion> &pvNew, int &flag, int &modify);
  void interpolation (MimmoObject* obj, std::vector<long int> idB, std::unordered_map< long int , std::vector<long> > &mapOneRing, std::unordered_map< long int , std::vector<double> > &mapCoeffModificati, bitpit::PiercedVector< Quaternion > &Q, bitpit::PiercedVector< Quaternion > &pvNew, int &flag, int &modify);

   
   void smoothing (MimmoObject* obj, std::vector<long int> idB, bitpit::PiercedVector<double> &pv_campoScalare, int &n); 
//    void smoothing (MimmoObject* obj, std::vector<long int> idB, bitpit::PiercedVector<Quaternion> &Q, int &n, int &flag); //, int & modify
   void smoothing (MimmoObject* obj,
                   std::vector<long int> idB, 
                   bitpit::PiercedVector< Quaternion > &Q,  
                   int &n, 
                   std::unordered_map< long int , std::vector<double> > &mapCoeffModificati ,
                   std::unordered_map <long int , std::vector<long>> & mapOneRing,
                   bitpit::PiercedVector< double > & err_iter, 
                   int &flag,  
                   int &modify);
  
   void modifyCOR (MimmoObject* obj, darray3E & COR, int flag);
   
   void plotQuaternionField (std::string dir, std::string name, MimmoObject * geom0_rot, bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T);
   void plotTranslationField (std::string dir, std::string name, MimmoObject * obj, bitpit::PiercedVector< std::array<double, 3> > &campo_vettoriale);
   
   //FUNZIONI commentate
   //   void evalBorderQuaternions (PatchKernel* m_patch, MimmoObject* geo, bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T);
   
   /*Funzioni in CalcoloOrtoOrto
    * 
    */
   bitpit::PiercedVector<double> evalCellOrth (MimmoObject * geom );
   bitpit::PiercedVector<double> evalOrth_ratio (MimmoObject * obj0 ,MimmoObject * obj1);
   dvector1D setCounter (bitpit::PiercedVector<double> pv);
   void plotDataQUAD (std::string dir, std::string name, MimmoObject * obj, dvector1D temp_pv, std::string NameProp );
//     double evalMinPV (bitpit::PiercedVector<double> * PV);
   void evalMinMaxPV (bitpit::PiercedVector<double>  &PV,  double &minimo, double &massimo);
   void getQuaternionC_I(MimmoObject * obj, MimmoObject * obj_r, darray3E & COR , bitpit::PiercedVector<Quaternion> &Q, bitpit::PiercedVector<Quaternion> &T);
//  bitpit::PiercedVector< std::array<double,3> > evalNormal (MimmoObject * obj);
   
   #endif /* __MYHEADER_HPP__ */
   
   
