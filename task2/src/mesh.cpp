#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;

void make_node_count_tag(Mesh *m)
{
	Tag tagNodeCount = m->CreateTag("Node_Count", DATA_INTEGER, CELL, NONE, 1);
	for (Mesh::iteratorCell icell = m->BeginCell(); icell != m->EndCell(); icell++) {
		Cell cell = icell->getAsCell();
		cell.Integer(tagNodeCount) = cell.getNodes().size(); // obviously 3 for out triangular meshes
	}
}

void make_area_tag(Mesh* m)
{
	Tag tagArea = m->CreateTag("Area", DATA_REAL, CELL, NONE, 1);
	for (Mesh::iteratorCell icell = m->BeginCell(); icell != m->EndCell(); icell++) {
		Cell cell = icell->getAsCell();
		cell.Real(tagArea) = cell.Volume(); // 'Volume' is volume for 3D meshes and area for 2D meshes
	}
}

void make_coord_tag(Mesh* m)
{
	Tag tagCoord = m->CreateTag("Coords", DATA_REAL, NODE, NONE, 3);
	for (Mesh::iteratorNode inode = m->BeginNode(); inode != m->EndNode(); inode++) {
		Node node = inode->getAsNode();
		node.RealArray(tagCoord)[0] = node.Coords()[0];
		node.RealArray(tagCoord)[1] = node.Coords()[1];
		node.RealArray(tagCoord)[2] = 0;
	}
}

int main(int argc, char *argv[])
{
	if (argc != 2) {
		cout << "Usage: mesh <mesh.vtk>" << endl;
		return 1;
	}
    Mesh *m = new Mesh;
	m->Load(argv[1]);
	make_node_count_tag(m);
	make_area_tag(m);
	make_coord_tag(m);
	m->Save("res.vtk");
	delete m;
	cout << "Success!";
	return 0;
}
