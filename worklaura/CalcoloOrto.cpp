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
#include "myHeader.hpp"

#include "mimmo_manipulators.hpp"
#include "mimmo_utils.hpp"
#include "bitpit.hpp"
#include "mimmo_geohandlers.hpp"
#include "Operators.hpp"
#include "mimmo_iogeneric.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;

/*Voglio creare una funzione che mi restituisca un piersed vector con l'ortogonalità 
 *che prenda un MimmoGeometry
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
    m_patch->buildAdjacencies();
    m_patch->buildInterfaces(); 
    bitpit::VolUnstructured * nor = static_cast <bitpit::VolUnstructured * >(m_patch);

   /*Ciclo sulle interfacce calcolo i due vettori e il coseno quadro formato dai due vettori
    *Tengo in considerazione delle interfacce al bordo imponendo CosQuadAlpha=2
    */
    for(itInterfaces = m_patch->interfaceBegin();itInterfaces != m_patch->interfaceEnd();itInterfaces++){
    if (itInterfaces->isBorder()) {
    cosquad = 2;  
    PVCosQuadAlpha.insert(itInterfaces.getId(), cosquad);
    }
    else{
    FirstVec = m_patch->evalCellCentroid(itInterfaces->getOwnerNeigh()[0])-m_patch->evalCellCentroid(itInterfaces->getOwnerNeigh()[1]);
    SecondVec = nor->evalInterfaceNormal(itInterfaces.getId());
    cosquad = pow( (dotProduct(FirstVec,SecondVec)), 2.0)/ pow(norm2(FirstVec), 2.0)* pow(norm2(SecondVec), 2.0);  
    PVCosQuadAlpha.insert(itInterfaces.getId(), cosquad);    
    } 
}
    for(auto &cell : m_patch->getCells()){
      int dim = cell.getInterfaceCount();
      long * CellInter = cell.getInterfaces();
      double minimo = PVCosQuadAlpha[CellInter[0]]; 
      
      for(int j=1; j<dim;j++){
	minimo = min(minimo, PVCosQuadAlpha[cell.getInterfaces()[j]]);
	};
      Orth = minimo;
      if (Orth > 1.0) {
	Orth = 2;
	}
      PVOrt.insert(cell.getId(), Orth);
    }
    return(PVOrt);
}

/*Creo una funzione che calcoli Ort_ratio= Orth_deformed/Orth_undeformed = Orth_obj1/Orth_obj0
 */
bitpit::PiercedVector<double> evalOrth_ratio (MimmoObject * obj0 ,MimmoObject * obj1){
  bitpit::PiercedVector<double> pvOrth_ratio;
  bitpit::PiercedVector<double> temp_0 = evalCellOrth(obj0);
  bitpit::PiercedVector<double> temp_1 = evalCellOrth(obj1);
  for(auto &cell : obj0->getPatch()->getCells()){ 
    double a = temp_1[cell.getId()]/temp_0[cell.getId()];
    pvOrth_ratio.insert(cell.getId() , a);
  }

   return(pvOrth_ratio);
}



dvector1D setCounter (bitpit::PiercedVector<double> pv){
   dvector1D temp_pv(pv.size()); 
  {
    int counter = 0;
    for(auto val : pv){
      temp_pv[counter] = val;
      ++counter;
    }
  }
return(temp_pv);
}


void plotDataQUAD (std::string dir, std::string name, MimmoObject * obj, dvector1D temp_pv, std::string NameProp ){
 
 bitpit::VTKElementType cellType = bitpit::VTKElementType::QUAD;
 
 bitpit::VTKUnstructuredGrid output(dir,name,cellType);
 auto vert = obj->getVertexCoords();
 auto connc = obj->getCompactConnectivity();
 output.setGeomData( bitpit::VTKUnstructuredField::POINTS, vert );
 output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connc );
 output.setDimensions(connc.size(), vert.size());
 output.addData(NameProp , bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, temp_pv);
 output.setCodex(bitpit::VTKFormat::APPENDED);
 output.write();
 }

/*Voglio creare una funzione che prende un MimmoObject e restituisce un PiercedVector contenente l'area di ogni cella
 */
bitpit::PiercedVector<double> evalCellsArea (MimmoObject * obj ){
   bitpit::PiercedVector<double> PVevalCellsArea;
   obj->getPatch()->buildAdjacencies();
   obj->getPatch()->buildInterfaces(); 
   
   bitpit::VolUnstructured * vol = static_cast <bitpit::VolUnstructured * >(obj->getPatch());
   for(auto &cell : vol->getCells()){
     auto volume = vol->evalCellVolume(cell.getId());
     PVevalCellsArea.insert(cell.getId(), volume);
   }
   
dvector1D temp_pvarea(PVevalCellsArea.size()); 
{
  int counter = 0;
  for(auto val : PVevalCellsArea){
    temp_pvarea[counter] = val;
    ++counter;
  }
}
 std::string dir = ".";
 std::string name = "area";
 bitpit::VTKElementType cellType = bitpit::VTKElementType::QUAD;
 
 bitpit::VTKUnstructuredGrid output(dir,name,cellType);
 auto vert = obj->getVertexCoords();
 auto connc = obj->getCompactConnectivity();
 output.setGeomData( bitpit::VTKUnstructuredField::POINTS, vert );
 output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connc );
 output.setDimensions(connc.size(), vert.size());
 output.addData("CellArea", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, temp_pvarea);
 output.setCodex(bitpit::VTKFormat::APPENDED);
 output.write();
 
  return(PVevalCellsArea); 
}

bitpit::PiercedVector<double> evalArea_ratio (MimmoObject * obj0 ,MimmoObject * obj1){
  bitpit::PiercedVector<double> pvArea_ratio;
  bitpit::PiercedVector<double> temp_0 = evalCellsArea(obj0);
  bitpit::PiercedVector<double> temp_1 = evalCellsArea(obj1);
  for(auto &cell : obj0->getPatch()->getCells()){ 
    double a = (temp_1[cell.getId()])/temp_0[cell.getId()];
    pvArea_ratio.insert(cell.getId() , a);
  }

   return(pvArea_ratio);
}


 

 /*Voglio creare una funzione che effettua una deformazione
 */
void evalDef (MimmoGeometry*  geom ){
  
    geom->setIOMode(IOMode::WRITE);
    geom->setWriteDir(".");
    geom->setWriteFileType(FileType::SQVTUPLANAR);
    geom->setWriteFilename("laura_func_def_00001.0000");
 

    /* Create MRBF block.
     */
    dvecarr3E rbfNodes(1,{{1.0,0.0,5.3522e-08}});
    dvecarr3E displ(1,{{0.0,0.0,-0.2}});
   
   MRBF* mrbf = new MRBF();
    mrbf->setMode(MRBFSol::NONE);
    mrbf->setSupportRadius(0.0003);
    mrbf->setPlotInExecution(true);
    mrbf->setNode(rbfNodes);
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
/*
 * 
 */
void test4() {

    MimmoGeometry * mimmo0 = new MimmoGeometry();
    mimmo0->setIOMode(IOMode::READ);
    mimmo0->setReadDir("./");
    mimmo0->setReadFilename("naca0012sq");
    mimmo0->setReadFileType(FileType::SQVTUPLANAR);
    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::SQVTUPLANAR);
    mimmo0->setWriteFilename("laura_output_00001.0000");
    mimmo0->execute();
    
    MimmoGeometry * mimmo0_HARDCopy = new MimmoGeometry();
    mimmo0_HARDCopy->setHARDCopy(mimmo0);
    
    auto obj_undef = mimmo0->getGeometry();
    evalDef(mimmo0_HARDCopy);
    auto obj_def = mimmo0_HARDCopy->getGeometry();
       
    bitpit::PiercedVector<double> Ort_undef = evalCellOrth(obj_undef);
    dvector1D tempOrth_undef = setCounter(Ort_undef);
    plotDataQUAD("./", "OrthogonalityUndef", obj_undef, tempOrth_undef, "Orthogonality");
    
    bitpit::PiercedVector<double> Area_undef = evalCellsArea(obj_undef);
    dvector1D tempArea_undef = setCounter(Area_undef);
    plotDataQUAD("./", "Area_undef", obj_undef, tempArea_undef, "Area");
    
    bitpit::PiercedVector<double> Ort_def = evalCellOrth(obj_def);
    dvector1D temp_def = setCounter(Ort_def);
    plotDataQUAD("./", "OrthogonalityDef", obj_def, temp_def, "Orthogonality");
        
    bitpit::PiercedVector<double> Area_def = evalCellsArea(obj_def);
    dvector1D tempArea_def = setCounter(Area_def);
    plotDataQUAD("./", "Area_def", obj_def, tempArea_def, "Area");
    
    bitpit::PiercedVector<double> Ort = evalOrth_ratio(obj_undef,obj_def);
    dvector1D temp_Orth_ratio = setCounter(Ort);
    plotDataQUAD("./", "Orthogonality_Ratio", obj_def, temp_Orth_ratio, "Orthogonality_ratio");
    
    bitpit::PiercedVector<double> Ar = evalArea_ratio(obj_undef,obj_def);
    dvector1D temp_Area_ratio = setCounter(Ar);
    plotDataQUAD("./", "Area_Ratio", obj_def, temp_Area_ratio, "Area_ratio");
    
    delete mimmo0;
    delete mimmo0_HARDCopy;
     
    mimmo0 = NULL;
    mimmo0_HARDCopy = NULL;
     
    return;  
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

        test4() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
//	return val;
	
}
