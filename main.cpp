#include <QCoreApplication>
#include <iostream>
#include "lsd.h"

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
using namespace vcg;
class CFace;
class CVertex;
class CEdge;
struct MyUsedTypes : public vcg::UsedTypes<	Use<CVertex>::AsVertexType,	Use<CEdge>::AsEdgeType,Use<CFace>::AsFaceType>{};

/// compositing wanted proprieties
class CVertex : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d,vcg::vertex::VFAdj, vcg::vertex::Color4b, vcg::vertex::TexCoord2d,vcg::vertex::BitFlags,vcg::vertex::VEAdj, vcg::vertex::Qualityd,vcg::vertex::Mark>{};
class CEdge   : public vcg::Edge<MyUsedTypes, vcg::edge::VertexRef, vcg::edge::BitFlags,vcg::edge::EEAdj,vcg::edge::EFAdj,vcg::edge::VEAdj>{};
class CFace   : public vcg::Face<  MyUsedTypes, vcg::face::VertexRef, vcg::face::FFAdj,vcg::face::VFAdj, vcg::face::Normal3d, vcg::face::Color4b, vcg::face::BitFlags,vcg::face::Mark,vcg::face::WedgeTexCoord2d> {};
class CMesh   : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace>,std::vector<CEdge>> {};

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    CMesh m;
    vcg::tri::io::Importer<CMesh>::Open(m,"../../Deformation10/models/reduced100Old.ply");
    UpdateNormal<CMesh>::PerVertexNormalizedPerFaceNormalized(m);
    //vector<size_t> roi={1,2,3,4};
    //vector<size_t> handle={1,2};
    vector<size_t> roi,handle;
    vector<CMesh::CoordType> newpos;
    CMesh::ScalarType scale=30;
    for(CMesh::VertexIterator vi=m.vert.begin();vi!=m.vert.end();vi++)
        if(vi->P().Y()<230 && vi->P().Y()>10 ){
            roi.push_back(vcg::tri::Index(m,&*vi));
            if(vi->P().Y()<130 && vi->P().Y()>110 ){
                handle.push_back(vcg::tri::Index(m,&*vi));
                newpos.push_back(vi->P()+vi->N()*scale);
            }

        }
    //tri::io::ExporterPLY<CMesh>::Save(m,"antes.ply",false,vcg::tri::io:);
    cout<<"comeza "<<endl;
    LaplacianSurfaceDeformation<CMesh> lsd(m,roi,handle);
    lsd.prepareDeform();
    lsd.doDeform(newpos);
    tri::io::ExporterPLY<CMesh>::Save(m,"deforma.ply");
    return 0;

}
