#include <vtkPolyData.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyhedron.h>
#include <vtkPolyLine.h>
#include <vtkFeatureEdges.h>
#include <vtkDiskSource.h>
#include <vtkFloatArray.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <sstream>
#include <map>
#include <vector>
#include <vtkProperty.h>

// first => p1 p2 from a edge; second
typedef std::pair<std::vector<vtkIdType>,std::vector<vtkIdType>> PairType;

std::map<vtkIdType, std::vector<vtkIdType>> mapPrimalPointByDualPolygon;
std::map<std::string, vtkIdType> mapDualPolygonByPrimalEdge;
std::map<vtkIdType, std::vector<vtkIdType>> mapDualCellsByDualVertex;

#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

//#define DEBUG

void WriteMeshToVTK(vtkSmartPointer<vtkUnstructuredGrid> polyData, std::string filename);

vtkSmartPointer<vtkUnstructuredGrid> ReadMeshFromVTK(std::string filename);

void findDualPoints(int count, vtkIdType cellCounter, vtkSmartPointer<vtkUnstructuredGrid> &mesh,
		    vtkSmartPointer<vtkPoints> &intermediaryMeshPoints, vtkSmartPointer<vtkIdList> &idListPoints);

vtkIdType createDualPolygon(vtkSmartPointer<vtkUnstructuredGrid> &mesh, vtkSmartPointer<vtkUnstructuredGrid> & lambdaMesh, vtkIdList * idsCells);

void getEdgeCells (vtkSmartPointer<vtkUnstructuredGrid> & _mesh, PairType* element_mapEdge, vtkCellArray * cells, vtkIdList * idsCells);

bool compareCellsByFaces(vtkSmartPointer<vtkUnstructuredGrid> & _mesh, vtkIdType cellId1, vtkIdType cellId2);

void compareCellsByFaces(vtkSmartPointer<vtkUnstructuredGrid> & _mesh, vtkIdType cellId,
			 vtkIdList * listIdsCellsToCompare, vtkIdList * listIdsCellsNeighbors);

void fillMapDualCellsByDualVertex(vtkSmartPointer<vtkUnstructuredGrid> & _mesh, std::map<vtkIdType, std::vector<vtkIdType>> & mapPolygonByDualVertex,
				  std::map<vtkIdType, std::vector<vtkIdType>> & mapPolygonByPrimalVertex,
				  vtkSmartPointer<vtkUnstructuredGrid> & meshToWrite);

void WriteMeshToPolyVTK(vtkSmartPointer<vtkPolyData> polyData, std::string filename);



void WriteMeshToVTK(vtkSmartPointer<vtkUnstructuredGrid> polyData,
		    std::string filename)
{
  //Write the mesh.
  vtkSmartPointer<vtkUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  std::stringstream ssFilename;
  ssFilename << filename << ".vtk";
  writer->SetFileName(ssFilename.str().c_str());
  writer->SetInputData(polyData);
  writer->SetFileTypeToASCII();
  writer->Update();
}

vtkSmartPointer<vtkUnstructuredGrid> ReadMeshFromVTK(std::string filename)
{
  vtkSmartPointer<vtkUnstructuredGridReader> reader =
    vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName ( filename.c_str() );
  reader->Update();

  //Obtain the polydata.
  vtkSmartPointer<vtkUnstructuredGrid> polyData = reader->GetOutput();
  return polyData;
}

void findDualPoints(int count, vtkIdType cellCounter, vtkSmartPointer<vtkUnstructuredGrid> &mesh, vtkSmartPointer<vtkPoints> &dualMesh,
		    vtkSmartPointer<vtkIdList> &idListPoints)
{
  double p[3];
  double x=0;
  double y=0;
  double z=0;

  for (unsigned int i = 0; i < idListPoints->GetNumberOfIds(); ++i)
    {
      vtkIdType v = idListPoints->GetId(i);
      mesh->GetPoint(v,p);
      x += p[0];
      y += p[1];
      z += p[2];
    }
  x /= (double)idListPoints->GetNumberOfIds();
  y /= (double)idListPoints->GetNumberOfIds();
  z /= (double)idListPoints->GetNumberOfIds();

  if(count < 10) //print the 10th first duals points
    {
      std::cout << "point " << cellCounter <<" : ";
      std::cout << x << " ";
      std::cout << y << " ";
      std::cout << z;
      std::cout << std::endl;
    }

  dualMesh->InsertNextPoint(x,y,z);
}

vtkIdType createDualPolygon(vtkSmartPointer<vtkUnstructuredGrid> &mesh,
			    vtkSmartPointer<vtkUnstructuredGrid> & lambdaMesh, vtkIdList * idsCells)
{

  vtkSmartPointer<vtkIdList> polygon = vtkSmartPointer<vtkIdList>::New();

  /*std::cout<<cells->GetNextCell(idListPoints1)<<" ";
    std::cout << idsCells->GetNumberOfIds() << std::endl;*/

  vtkIdType idPreviousCell = -1;
  vtkIdType idCurrentCell = idsCells->GetId(0);


  while (polygon->GetNumberOfIds() < idsCells->GetNumberOfIds())
  {
    vtkSmartPointer<vtkIdList> adjacentCells = vtkSmartPointer<vtkIdList>::New();

    compareCellsByFaces(mesh, idCurrentCell, idsCells, adjacentCells);

    std:cout <<"number of Id from adjacentCells"<< adjacentCells->GetNumberOfIds() << std::endl;
    if (adjacentCells->GetNumberOfIds() == 2)
		{
		  polygon->InsertNextId(idCurrentCell);

		  if (idPreviousCell == -1)
		  {
		    idPreviousCell = idCurrentCell;
		    idCurrentCell = adjacentCells->GetId(0);
		  }
		  else
		  {
		    if(idPreviousCell == adjacentCells->GetId(0))
				{
				  idPreviousCell = idCurrentCell;
				  idCurrentCell = adjacentCells->GetId(1);
				}
		    else
				{
				  idPreviousCell = idCurrentCell;
				  idCurrentCell = adjacentCells->GetId(0);
				}
		  }

		}
	  else
			break;

	}

  vtkIdType idNewPolygon = -1;
  if (polygon->GetNumberOfIds() == idsCells->GetNumberOfIds()) {
    idNewPolygon = lambdaMesh->InsertNextCell(VTK_POLYGON, polygon);
    //mapPrimalPointByDualPolygon->insert(std::make_pair(idNewPolygon, ));
  }

  return idNewPolygon;
}

void getEdgeCells (vtkSmartPointer<vtkUnstructuredGrid> & _mesh, PairType* element_mapEdge/*vtkIdType cellId*/,
		   vtkCellArray * cells, vtkIdList * idsCells) {
  /**
   * Get the neightbors cells at an edge
   */

  // All points ids from a cell
  // two points ids !defined in tough! forms a segment of current cell
  vtkSmartPointer<vtkIdList> edgeByTwoPtIds = vtkSmartPointer<vtkIdList>::New();
  edgeByTwoPtIds->InsertNextId(element_mapEdge->first[0]);
  edgeByTwoPtIds->InsertNextId(element_mapEdge->first[1]);

  vtkSmartPointer<vtkIdList> cellIdsNeighborsFromEdge = vtkSmartPointer<vtkIdList>::New();

  _mesh->GetCellNeighbors(element_mapEdge->second[0],edgeByTwoPtIds,  cellIdsNeighborsFromEdge);


  // Stock all the points from Current Cell in vtkCellArray
  idsCells->InsertNextId(element_mapEdge->second[0]);
  cells->InsertNextCell(_mesh->GetCell(element_mapEdge->second[0]));

  for (int counterCurrentCell = 0; counterCurrentCell < cellIdsNeighborsFromEdge->GetNumberOfIds(); counterCurrentCell++) {
    vtkIdType idCurrentCell = cellIdsNeighborsFromEdge->GetId(counterCurrentCell);
    cells->InsertNextCell(_mesh->GetCell(idCurrentCell));
    idsCells->InsertNextId(idCurrentCell);
  }

}

bool compareCellsByFaces(vtkSmartPointer<vtkUnstructuredGrid> & _mesh, vtkIdType cellId1, vtkIdType cellId2) {
  /**
   *	compare two cell by their	faces, determine if they are adjacent.
   */

  // Recovere all points ids from each cell
  vtkSmartPointer<vtkIdList> ptIdsFromCell1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> ptIdsFromCell2 = vtkSmartPointer<vtkIdList>::New();
  _mesh->GetCellPoints(cellId1, ptIdsFromCell1);
  _mesh->GetCellPoints(cellId2, ptIdsFromCell2);

  // Course all points ids from each cell and add the points ids common in list
  vtkSmartPointer<vtkIdList> resultPointsIdsCommonCompareFromCells = vtkSmartPointer<vtkIdList>::New();
  for (int counterPointIdFromCell1 = 0; counterPointIdFromCell1 < ptIdsFromCell1->GetNumberOfIds(); counterPointIdFromCell1++) {
    for (int counterPointIdFromCell2 = 0; counterPointIdFromCell2 < ptIdsFromCell2->GetNumberOfIds(); counterPointIdFromCell2++) {
      if (ptIdsFromCell1->GetId(counterPointIdFromCell1) == ptIdsFromCell2->GetId(counterPointIdFromCell2)) {
				resultPointsIdsCommonCompareFromCells->InsertUniqueId(ptIdsFromCell1->GetId(counterPointIdFromCell1));
      }
    }
  }

  // true if face common
  if (resultPointsIdsCommonCompareFromCells->GetNumberOfIds() >= 3)
    return true;

  return false;
}

void compareCellsByFaces(vtkSmartPointer<vtkUnstructuredGrid> & _mesh, vtkIdType cellId,
			 vtkIdList * listIdsCellsToCompare, vtkIdList * listIdsCellsNeighbors) {
  /**
   *	compare cell with list cells by their faces, determine if they are adjacent.
   */

  // Recovere all points ids from cell
  vtkSmartPointer<vtkIdList> ptIdsFromCell = vtkSmartPointer<vtkIdList>::New();
  _mesh->GetCellPoints(cellId, ptIdsFromCell);

  // Course all ids from list cells
  for (int counterCurrentCellOfList = 0;
       counterCurrentCellOfList < listIdsCellsToCompare->GetNumberOfIds();
       counterCurrentCellOfList++)
  {
    if (listIdsCellsToCompare->GetId(counterCurrentCellOfList) != cellId)
		{
		  int counterCommonPoint = 0;
		  // Stock all points from current cell
	  	vtkSmartPointer<vtkIdList> ptIdsFromCellToCompare = vtkSmartPointer<vtkIdList>::New();
	  	_mesh->GetCellPoints(listIdsCellsToCompare->GetId(counterCurrentCellOfList),ptIdsFromCellToCompare);
	  	// Course all points ids from current cell

	  	for (int counterPointIdFromCellToCompare = 0;
	       counterPointIdFromCellToCompare < ptIdsFromCellToCompare->GetNumberOfIds();
	       counterPointIdFromCellToCompare++)
	    {
	      // Course all points ids from each cell
	      for (int counterPointIdFromCell1 = 0;
		   			counterPointIdFromCell1 < ptIdsFromCell->GetNumberOfIds();
		   			counterPointIdFromCell1++)
				{
		  		if (ptIdsFromCell->GetId(counterPointIdFromCell1) == ptIdsFromCellToCompare
		      		->GetId(counterPointIdFromCellToCompare))
		    		counterCommonPoint++;
				}
	    	if (counterCommonPoint >= 3)
					listIdsCellsNeighbors->InsertUniqueId(listIdsCellsToCompare->GetId(counterCurrentCellOfList));
	  	}
		}
  }
}

void fillMapDualCellsByDualVertex(vtkSmartPointer<vtkUnstructuredGrid> & _mesh, std::map<vtkIdType, std::vector<vtkIdType>> & mapPolygonByDualVertex,
				  std::map<vtkIdType, std::vector<vtkIdType>> & mapPolygonByPrimalVertex,
				  vtkSmartPointer<vtkUnstructuredGrid> & meshToWrite) {
  //vtkSmartPointer<vtkCellArray> primalCells = mesh->GetCells();
  /* --- Notes
   * std::map<vtkIdType, std::vector<vtkIdType>> mapPrimalPointByDualPolygon;
   *	std::map<std::string, std::vector<vtkIdType>> mapDualPolygonByPrimalEdge;
   *	std::map<vtkIdType, std::vector<vtkIdType>> mapDualCellsByDualVertex;
   **/

  // Course all point dual in map
  for (std::map<vtkIdType, std::vector<vtkIdType>>::iterator it_mapPolygonByDualVertex = mapPolygonByDualVertex.begin();
       it_mapPolygonByDualVertex != mapPolygonByDualVertex.end(); it_mapPolygonByDualVertex++)
	{
    if(it_mapPolygonByDualVertex->second.size() >= 3)
    {
      std::vector<vtkIdType> vectorDualCells;
      // Course element id polygon in vector from map Polygon By Dual Vertex
      std::cout << "Dual Vertex [ " << it_mapPolygonByDualVertex->first;
      for (std::vector<vtkIdType>::iterator it_elementVectorFromMapPolygonByDualVertex = it_mapPolygonByDualVertex->second.begin();
	   it_elementVectorFromMapPolygonByDualVertex != it_mapPolygonByDualVertex->second.end(); it_elementVectorFromMapPolygonByDualVertex++)
      {
				vectorDualCells.insert(vectorDualCells.end(), mapPrimalPointByDualPolygon[*it_elementVectorFromMapPolygonByDualVertex].begin(),
		       mapPrimalPointByDualPolygon[*it_elementVectorFromMapPolygonByDualVertex].end());

      }

      std::sort(vectorDualCells.begin(), vectorDualCells.end(), std::greater<int>());
      auto last  = std::unique(vectorDualCells.begin(), vectorDualCells.end());
      vectorDualCells.erase(last, vectorDualCells.end());

			std::vector<vtkIdType> vectorDualCellsCheck;

      std::cout << " ] Cells [ ";
      for (std::vector<vtkIdType>::iterator iter = vectorDualCells.begin();iter != vectorDualCells.end();iter++)
      {
				if (mapPolygonByPrimalVertex[*iter].size() > 1)
				  {
	    		//vectorDualCells.push_back(*it_elementVectorFromMapPolygonByDualVertex);
	    		std::cout << *iter << ", ";
					vectorDualCellsCheck.push_back(*iter);
					}
      }

      std::cout << " ]" << endl;

      mapDualCellsByDualVertex.insert(std::make_pair(it_mapPolygonByDualVertex->first, vectorDualCellsCheck));
    }
  }

	// ----------------Add intormations

	vtkSmartPointer<vtkPoints> primalPoints = _mesh->GetPoints();
	meshToWrite->SetPoints(primalPoints);

	for(std::map<vtkIdType,std::vector<vtkIdType>>::iterator it_mapDualCellsByDualVertex = mapDualCellsByDualVertex.begin();
	    it_mapDualCellsByDualVertex != mapDualCellsByDualVertex.end(); it_mapDualCellsByDualVertex ++ )
	{
			//std::cout << it_mapDualCellsByDualVertex->second.size() << endl;
			if ( it_mapDualCellsByDualVertex->second.size() > 3)
			{
				  meshToWrite->InsertNextCell(VTK_TETRA, (_mesh->GetCell(it_mapDualCellsByDualVertex->first)->GetPointIds()));
			}
	}


}

void WriteMeshToPolyVTK(vtkSmartPointer<vtkPolyData> polyData,
			std::string filename)
{
  //Write the mesh.
  vtkSmartPointer<vtkPolyDataWriter> writer =
    vtkSmartPointer<vtkPolyDataWriter>::New();
  std::stringstream ssFilename;
  ssFilename << filename << ".vtk";
  writer->SetFileName(ssFilename.str().c_str());
  writer->SetInputData(polyData);
  writer->SetFileTypeToASCII();
  writer->Update();
}


int main ( int argc, char *argv[] )
{
  if(argc != 2)
  {
      std::cout << "Usage: " << argv[0] << "  Filename(.vtk)" << " <Number of iterations>" << std::endl;
      return EXIT_FAILURE;
  }

  vtkSmartPointer<vtkUnstructuredGrid> intermediaryMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> intermediaryMeshPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> intermediaryMeshCells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> intermediaryMeshLines = vtkSmartPointer<vtkCellArray>::New();
  intermediaryMesh->SetPoints(intermediaryMeshPoints);

  vtkSmartPointer<vtkUnstructuredGrid> dualMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> dualMeshPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> dualMeshCells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> dualMeshLines = vtkSmartPointer<vtkCellArray>::New();
  dualMesh->SetPoints(dualMeshPoints);


  std::string inputFilename = argv[1];
  vtkSmartPointer<vtkUnstructuredGrid> mesh = ReadMeshFromVTK(inputFilename);

  mesh->BuildLinks();
  mesh->ComputeBounds();

  vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
  unsigned int numberOfPoints = points->GetNumberOfPoints();

  //contains ids Edges in key associated with the edge and the cells that compose it
  std::map<vtkIdType, PairType*> *mapEdges = new std::map<vtkIdType,PairType*>();

  vtkIdType idEdgeCounter = 0;

  for (unsigned int pointCounter = 0; pointCounter < numberOfPoints;
       ++pointCounter)
  {
      double* point = points->GetPoint(pointCounter);

      vtkSmartPointer<vtkIdList> adjacentCells = vtkSmartPointer<vtkIdList>::New();
      mesh->GetPointCells(pointCounter, adjacentCells);

      if (adjacentCells->GetNumberOfIds() > 1)
			{
	  		std::vector<vtkIdType> neighbouringVertices;
	  		// This loop extracts the neighbourhood of the current point.
	  		for (unsigned int i = 0; i < adjacentCells->GetNumberOfIds(); ++i)
	    	{
	      	vtkSmartPointer<vtkIdList> vertices = vtkSmartPointer<vtkIdList>::New();

	      	mesh->GetCellPoints(adjacentCells->GetId(i), vertices);

	      	for (unsigned int k = 0; k < vertices->GetNumberOfIds(); ++k)
					{
					  if (vertices->GetId(k) != pointCounter)
					    neighbouringVertices.push_back(vertices->GetId(k));
					}
	    	}

	  		std::sort(neighbouringVertices.begin(), neighbouringVertices.end(), std::greater<int>());
	  		auto last  = std::unique(neighbouringVertices.begin(), neighbouringVertices.end());
	  		neighbouringVertices.erase(last, neighbouringVertices.end());


	  		for (unsigned int i = 0; i < neighbouringVertices.size(); ++i)
	    	{
	      //====== test : Etraire les arêtes voisines d'un sommet i ========

	      	vtkSmartPointer<vtkIdList> line = vtkSmartPointer<vtkIdList>::New();
	      	std::vector<vtkIdType> edge;

	      	vtkIdType p1 = pointCounter;
	      	vtkIdType p2 = neighbouringVertices[i];

	      	if(p1 < p2)
					{
					  //mapDualPolygonByPrimalEdge
					  //line->InsertNextId(p1);
					  //line->InsertNextId(p2);

		  			edge.push_back(p1);
		  			edge.push_back(p2);

		  			vtkSmartPointer <vtkIdList> adjacentCells1 = vtkSmartPointer<vtkIdList>::New();
		  			mesh->GetPointCells(p1, adjacentCells1);

		  			vtkSmartPointer <vtkIdList> adjacentCells2 = vtkSmartPointer<vtkIdList>::New();
		  			mesh->GetPointCells(p2, adjacentCells2);

		  			//vtkSmartPointer<vtkIdList> cellsBelongEdge = vtkSmartPointer<vtkIdList>::New();
		  			std::vector<vtkIdType> cellsBelongEdge;

		  			for(int i= 0; i< adjacentCells1->GetNumberOfIds(); i++ )
		    		{
		      		for(int j= 0; j< adjacentCells2->GetNumberOfIds(); j++ )
							{
			  				if(adjacentCells1->GetId(i) == adjacentCells2->GetId(j))
			    			{
			      			//cellsBelongEdge->InsertUniqueId(adjacentCells1->GetId(i));
			      			cellsBelongEdge.push_back(adjacentCells1->GetId(i));
			    			}
							}

		    		}
		  			std::sort(cellsBelongEdge.begin(), cellsBelongEdge.end(), std::greater<vtkIdType>());
		  			auto last  = std::unique(cellsBelongEdge.begin(), cellsBelongEdge.end());
		  			cellsBelongEdge.erase(last, cellsBelongEdge.end());

		  			//print informations about an edge E, the twice vertex id of E and all the neigbors cell from the edge E
		  			std::cout << "Id arete : " << idEdgeCounter << " -> id sommets : " << p1 << " " << p2 << "-> id cells : " ;

		  			for(int i=0; i<cellsBelongEdge.size(); i++)
		    			std::cout << cellsBelongEdge[i] <<  " ";
		  			std::cout << endl;

		  			//mapEdges->insert(std::pair<vtkIdType, vtkIdList * >(idEdgeCounter , line));
		  			PairType *pairLineCells = new PairType(edge,cellsBelongEdge);
		  			//std::cout << " number vertex on cells :  " <<pairLineCells->first.size() <<endl;
		  			//std::cout <<"Id arete : " <<idEdgeCounter <<" number of cells neigbors :  " <<pairLineCells->second.size() <<endl;

		  			mapEdges->insert(std::pair<vtkIdType, PairType*>(idEdgeCounter , pairLineCells));

		  			//std::cout << mapEdges->at(idEdgeCounter)->first <<endl;

		  			//std::cout << "Identifiant ligne : " << idEdgeCounter << "id sommets : " << p1 << " " << p2 <<endl;
		  			idEdgeCounter ++;

					}

	    	}

			}

  }

  std::map<vtkIdType,std::vector<vtkIdType>> mapPolygonByPrimalVertex;
  std::map<vtkIdType,std::vector<vtkIdType>> mapPolygonByDualVertex;


  unsigned int nEdges = 0;
  for(std::map<vtkIdType,PairType*>::iterator it_mapEdges = mapEdges->begin(); it_mapEdges != mapEdges->end();it_mapEdges ++)
  {
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkIdList> idsCells = vtkSmartPointer<vtkIdList>::New();
    getEdgeCells(mesh, it_mapEdges->second , cells, idsCells);

    nEdges++;
    vtkIdType idPolygon = createDualPolygon(mesh, dualMesh,idsCells);

		if (idPolygon >= 0)
		{
	  	std::stringstream ssEdge;
	  	std::stringstream ssEdgeReverse;
	  	ssEdge << it_mapEdges->second->first[0] << '_' << it_mapEdges->second->first[1];
	  	ssEdgeReverse << it_mapEdges->second->first[1] << '_' << it_mapEdges->second->first[0];
	  	mapDualPolygonByPrimalEdge.insert(std::make_pair(ssEdge.str(), idPolygon));
	  	mapDualPolygonByPrimalEdge.insert(std::make_pair(ssEdgeReverse.str(), idPolygon));

	  	std::vector<vtkIdType> vectorIdsEdge;
	  	vectorIdsEdge.push_back(it_mapEdges->second->first[0]);
	  	vectorIdsEdge.push_back(it_mapEdges->second->first[1]);
	  	mapPrimalPointByDualPolygon.insert(std::make_pair(idPolygon , vectorIdsEdge));
	  	//std::cout<< " test : "<<(*it_mapEdges).second->first[0] <<endl;

	  	// --------------------- Polygon By Primal Vertex ------------------------
	  	std::map<vtkIdType,std::vector<vtkIdType>>::iterator point1Exist = mapPolygonByPrimalVertex.find((*it_mapEdges).second->first[0]);
	  	if (point1Exist != mapPolygonByPrimalVertex.end())
	    	mapPolygonByPrimalVertex[(*it_mapEdges).second->first[0]].push_back(idPolygon);
	  	else
	    {
	      std::vector<vtkIdType> polygons;
	      polygons.push_back(idPolygon);
	      mapPolygonByPrimalVertex.insert(std::make_pair((*it_mapEdges).second->first[0],polygons));
	    }

	  	std::map<vtkIdType,std::vector<vtkIdType>>::iterator point2Exist = mapPolygonByPrimalVertex.find((*it_mapEdges).second->first[1]);
	  	if (point2Exist != mapPolygonByPrimalVertex.end())
	    	mapPolygonByPrimalVertex[(*it_mapEdges).second->first[1]].push_back(idPolygon);
	  	else
	    {
	      std::vector<vtkIdType> polygons;
	      polygons.push_back(idPolygon);
	      mapPolygonByPrimalVertex.insert(std::make_pair((*it_mapEdges).second->first[1],polygons));
	    }
	  	// --------------------------------------------------------------------------
	  	for (int counterIdsCells = 0; counterIdsCells < idsCells->GetNumberOfIds(); counterIdsCells++)
			{
	    	std::map<vtkIdType,std::vector<vtkIdType>>::iterator dualPointExist = mapPolygonByDualVertex.find(idsCells->GetId(counterIdsCells));
	    	if (dualPointExist != mapPolygonByDualVertex.end())
	      	mapPolygonByDualVertex[idsCells->GetId(counterIdsCells)].push_back(idPolygon);
	    	else
	      {
					std::vector<vtkIdType> polygons;
					polygons.push_back(idPolygon);
					mapPolygonByDualVertex.insert(std::make_pair(idsCells->GetId(counterIdsCells), polygons));
	      }
	  	}
		}
      //	mapPolygonByPrimalVertex.insert((*it_mapEdges).second->first[1],idPolygon);
  }


  //std::cout<< "number of edges : " << nEdges << "map number element of Polygon : "<<mapPolygonByPrimalVertex.size() << " / " << mapPolygonByPrimalVertex[12].size()<<endl;
  for (unsigned int i = 0;i < mapPolygonByPrimalVertex.size();i++)
    std::cout<< "Nombre de facettes pour le sommet primal : " << i  << " = " << mapPolygonByPrimalVertex[i].size()<<endl;
  for (unsigned int i = 0;i < mapPolygonByDualVertex.size();i++)
    std::cout<< "Nombre de facettes pour le sommet dual : " << i  << " = " << mapPolygonByDualVertex[i].size()<<endl;


  std::cout << "Number of Cells in dualMesh : " << dualMesh->GetCells()->GetNumberOfCells() << endl;
  fillMapDualCellsByDualVertex(mesh, mapPolygonByDualVertex, mapPolygonByPrimalVertex, intermediaryMesh);

  //std::cout << idsCells->GetNumberOfIds()<<endl;

  // ------ TODO :
  //recup tout les points du dual qui deviendrons primal
  //getpoints
  //setpoints

  //--------------------------------------
  // Course all point dual in map
  /*for (std::map<vtkIdType, std::vector<vtkIdType>>::iterator it_mapPolygonByPrimalVertex
    = mapPolygonByPrimalVertex.begin(); it_mapPolygonByPrimalVertex != mapPolygonByPrimalVertex.end();
    it_mapPolygonByPrimalVertex++)
    {
    //std::cout << it_mapPolygonByDualVertex->second.size();
    if(it_mapPolygonByPrimalVertex->second.size() >= 3)
    {
    //std::cout << it_mapPolygonByDualVertex <<endl;
    vtkIdType npts;
    vtkSmartPointer<vtkIdList> idsCells = vtkSmartPointer<vtkIdList>::New();
    // course all polygon from map iterator
    for (std::vector<vtkIdType>::iterator it_mapIdPolygonFromVector = it_mapPolygonByPrimalVertex->second.begin(); it_mapIdPolygonFromVector != it_mapPolygonByPrimalVertex->second.end(); it_mapIdPolygonFromVector++) {
    vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
    //std::cout<<*it_mapIdPolygonFromVector <<endl;
    dualMesh->GetCellPoints (*it_mapIdPolygonFromVector, pts);
    for (int i = 0; i < pts->GetNumberOfIds(); i++)
    idsCells->InsertNextId(pts->GetId(i));
    }

    //createDualPolygon(mesh, intermediaryMesh, idsCells);

    //intermediaryMesh->InsertNextCell
    }
    }*/
  //intermediaryMesh->SetPoints(dualMeshPoints);

  //--------------------------------------

  // --------- Set DualMesh --------------
  // TODO : Test
  // Meth 1 :
  //intermediaryMesh->SetCells(VTK_TETRA, cells);

  // Meth 2 : => Functional
  //intermediaryMesh->SetPoints(points);
  //intermediaryMesh->SetCells(VTK_TETRA, cells);

  // Meth 3 :
  /*for (int counterCurrentCell = 0; counterCurrentCell < cells->GetNumberOfCells(); counterCurrentCell++) {
    vtkSmartPointer<vtkIdList> ptsIdsFromCurrentCell = cells->;
    intermediaryMesh->InsertNextCell(VTK_POLYGON, ptsIdsFromCurrentCell);
    }*/
  // -------------------------------------
  std::cout << "compare cells 0 and 1 by faces" << compareCellsByFaces(mesh, 0, 3) << endl;
  /**
     Create the cell data to store the index to the dual vertex of the cell.
     in the primal triangular cells.
  */
  vtkSmartPointer<vtkIdTypeArray> cellData =
    vtkSmartPointer<vtkIdTypeArray>::New();
  mesh->GetCellData()->SetScalars(cellData);

  /*
   * Contains information related to the cells.
   * - Index of the primal vertex at the origin of the dual cell.
   * - Value attached to the primal vertex (curvature, indication value, etc).
   *
   */
  vtkSmartPointer<vtkFloatArray> arrayData =
    vtkSmartPointer<vtkFloatArray>::New();
  arrayData->SetNumberOfComponents(5);
  arrayData->SetName ("ArrayData");
  mesh->GetCellData()->AddArray(arrayData);

  /**
   * Traverse all the 2-dimensional cells in the polyData in order to create
   * the dual vertices.
   * This loop also detects the primal cells (faces) that are only connected to
   * one other face.
   *
   */
  vtkSmartPointer<vtkCellArray> primalCells = mesh->GetCells();
  vtkSmartPointer<vtkIdList> idListPoints = vtkSmartPointer<vtkIdList>::New();
  vtkIdType cellCounter = 0;

  int count = 0;

  for (cellCounter = 0, primalCells->InitTraversal();
       primalCells->GetNextCell(idListPoints) ;
       ++cellCounter)
  {
      //std::cout << "Les ids des points de la cellule " << cellCounter << " sont ";
      std::vector<vtkIdType> CellPoints;

      //======= set dual point  ========
      findDualPoints(count,cellCounter,mesh,dualMeshPoints,idListPoints);

      count ++;

      /*
			Adding data to the primal and dual meshes. We have to store the coordinates
			and the index of the barycenter in each cell of the primal mesh.
      */

      double* cellCounterData = new double(cellCounter);

      vtkIdType position = mesh->GetCellData()->GetScalars()
													->InsertNextTuple(cellCounterData);
      delete cellCounterData;

      double* cellData = new double[5];
      cellData[0] = 1;
      cellData[1] = 1;
      cellData[2] = -1.0;
      cellData[3] = -1.0;
      cellData[4] = 0.0;
      vtkIdType idNewCellData = arrayData->InsertNextTuple(cellData);
      delete[] cellData;
    }

  dualMesh->SetPoints(dualMeshPoints);


  //Write the dual dual PolyData.
  std::stringstream ssDualDual;
  ssDualDual << "ProcessedMesh";
  WriteMeshToVTK(mesh, ssDualDual.str().c_str());

  std::stringstream ssDualDual2;
  ssDualDual2 << "ProcessedintermediaryMesh";
  WriteMeshToVTK(intermediaryMesh, ssDualDual2.str().c_str());

  std::stringstream ssDualDual3;
  ssDualDual3 << "ProcessedDualMesh";
  WriteMeshToVTK(dualMesh, ssDualDual3.str().c_str());

  return EXIT_SUCCESS;
}
