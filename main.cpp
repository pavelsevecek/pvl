#include "Cloud.hpp"
#include "Graph.hpp"
#include "KdTree.hpp"
#include "Kernels.hpp"
#include "PlyReader.hpp"
#include "PlyWriter.hpp"
#include "Refinement.hpp"
#include "Simplification.hpp"
#include "Svd.hpp"
#include "UniformGrid.hpp"
#include <iostream>
#include <sstream>

#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp
                          // file
#include "catch.hpp"

using namespace Pvl;
/*std::cout << "Test" << std::endl;
Pvl::KdTree<Pvl::Vec3f> tree;
std::vector<Pvl::Vec3f> points = {
    { 0.f, 0.f, 0.f },
    { 0.f, 0.f, 4.f },
    { 1.f, 0.f, 0.f },
};
tree.build(points);

// std::vector<uint32_t> neighs;
std::array<uint32_t, 3> neighs;
int n = tree.rangeQuery(Pvl::Vec3f(0.f, 0.f, 0.f), 1.5f, neighs.begin());

std::cout << "Found " << n << " neighs" << std::endl;
for (auto& n : neighs) {
    std::cout << " - " << n << std::endl;
}

std::stringstream ss;
Pvl::PlyWriter ply(ss);
for (auto& p : points) {
    ply << p;
}
ply.close();
std::cout << ss.str() << "\n";

std::stringstream iss(ss.str());
Pvl::PlyReader rd(iss);
std::vector<Pvl::Vec3f> ps = rd.read();
std::cout << "Read points \n";
for (auto& p : ps) {
    std::cout << p[0] << "," << p[1] << "," << p[2] << "\n";
}

std::cout << "Grid\n";
Pvl::UniformGrid<float, 3> grid(Pvl::Vec3i(2, 3, 2));
int cntr = 0;
for (auto& f : grid) {
    std::cout << cntr++ << "-" << f << "\n";
}*/

/* Pvl::Graph graph;
 Pvl::VertexHandle a = graph.addVertex();
 Pvl::VertexHandle b = graph.addVertex();
 Pvl::VertexHandle c = graph.addVertex();
 Pvl::VertexHandle d = graph.addVertex();
 Pvl::VertexHandle e = graph.addVertex();
 Pvl::VertexHandle f = graph.addVertex();

 Pvl::FaceHandle f1 = graph.addFace(a, b, c);
 Pvl::FaceHandle f2 = graph.addFace(e, f, a);
 Pvl::FaceHandle f3 = graph.addFace(d, e, a);
 Pvl::FaceHandle f4 = graph.addFace(a, c, d);

 std::cout << "Faces = " << graph.numFaces() << std::endl;
 std::cout << "Vertices = " << graph.numVertices() << std::endl;

 for (Pvl::VertexHandle vh : graph.vertexRing(a)) {
     std::cout << " - " << vh.index() << std::endl;
 }*/

TEST_CASE("edge circulator inner", "[circulators]") {
    Graph graph;
    VertexHandle center = graph.addVertex();
    VertexHandle v0 = graph.addVertex();
    VertexHandle v1 = graph.addVertex();
    VertexHandle v2 = graph.addVertex();
    VertexHandle v3 = graph.addVertex();
    graph.addFace(v0, v1, center);
    graph.addFace(v1, v2, center);
    graph.addFace(v2, v3, center);
    graph.addFace(v3, v0, center);

    Vertex::EdgeRange edgeRing = graph.edgeRing(center);
    for (EdgeHandle eh : edgeRing) {
        std::cout << graph.from(graph.halfEdge(eh)) << "-" << graph.to(graph.halfEdge(eh))
                  << std::endl;
    }
}

TEST_CASE("edge circulator boundary", "[circulators]") {
    Graph graph;
    VertexHandle center = graph.addVertex();
    VertexHandle v0 = graph.addVertex();
    VertexHandle v1 = graph.addVertex();
    VertexHandle v2 = graph.addVertex();
    VertexHandle v3 = graph.addVertex();
    graph.addFace(center, v0, v1);
    graph.addFace(center, v1, v2);
    graph.addFace(center, v2, v3);

    Vertex::EdgeRange edgeRing = graph.edgeRing(center);
    for (EdgeHandle eh : edgeRing) {
        std::cout << graph.from(graph.halfEdge(eh)) << "-" << graph.to(graph.halfEdge(eh))
                  << std::endl;
    }
}

TEST_CASE("collapse allowed", "[simplify]") {
    Graph graph;
    VertexHandle vA = graph.addVertex();
    VertexHandle vB = graph.addVertex();
    VertexHandle v0 = graph.addVertex();
    VertexHandle v1 = graph.addVertex();

    graph.addFace(vA, vB, v0);
    graph.addFace(vA, v0, v1);
    graph.addFace(vB, v1, v0);

    HalfEdgeHandle eh = graph.halfEdge(vA, vB);
    REQUIRE_FALSE(graph.removed(eh));
    REQUIRE_FALSE(graph.collapseAllowed(eh));
}

VertexHandle operator"" _vh(unsigned long long int i) {
    return VertexHandle(i);
}

/* Creates a following graph
 *
 *   1 - 2 - 3
 *  / \ / \ / \
 * 0 - A - B - 4
 *  \ / \ / \ /
 *   7 - 6 - 5
 */
struct SimpleGraphFixture {
    Graph graph;
    VertexHandle vA;
    VertexHandle vB;
    std::array<VertexHandle, 8> vs;

    SimpleGraphFixture() {
        for (int i = 0; i < 8; ++i) {
            vs[i] = graph.addVertex();
        }
        vA = graph.addVertex();
        vB = graph.addVertex();

        graph.addFace(vs[0], vA, vs[1]);
        graph.addFace(vA, vs[2], vs[1]);
        graph.addFace(vA, vB, vs[2]);
        graph.addFace(vB, vs[3], vs[2]);
        graph.addFace(vB, vs[4], vs[3]);
        graph.addFace(vs[5], vs[4], vB);
        graph.addFace(vs[6], vs[5], vB);
        graph.addFace(vs[7], vs[6], vA);
        graph.addFace(vs[6], vB, vA);
        graph.addFace(vs[0], vs[7], vA);
    }
};

TEST_CASE_METHOD(SimpleGraphFixture, "halfedge", "[graph]") {
    HalfEdgeHandle eh = graph.halfEdge(vA, vB);
    REQUIRE(graph.valid(eh));
    REQUIRE(graph.from(eh) == vA);
    REQUIRE(graph.to(eh) == vB);
}

TEST_CASE_METHOD(SimpleGraphFixture, "edge range", "[graph]") {
    Graph::EdgeRange edges = graph.edgeRange();
    REQUIRE(std::distance(edges.begin(), edges.end()) == 19);
    REQUIRE(std::all_of(
        edges.begin(), edges.end(), [&](EdgeHandle eh) { return graph.valid(eh); }));

    HalfEdgeHandle heh = graph.halfEdge(vA, vB);
    int visited = std::count_if(edges.begin(), edges.end(), [this, heh](EdgeHandle eh) {
        return heh == graph.halfEdge(eh);
    });
    REQUIRE(visited == 1);
}

TEST_CASE_METHOD(SimpleGraphFixture, "collapse simple", "[simplify]") {
    Graph::EdgeRange edges = graph.edgeRange();
    HalfEdgeHandle collapsible = graph.halfEdge(vA, vB);

    REQUIRE(std::all_of(edges.begin(), edges.end(), [&](EdgeHandle eh) {
        HalfEdgeHandle heh = graph.halfEdge(eh);
        if (heh == collapsible) {
            return graph.collapseAllowed(heh);
        } else {
            return !graph.collapseAllowed(heh);
        }
    }));

    graph.collapse(collapsible);

    Graph::EdgeRange edges2 = graph.edgeRange();
    REQUIRE(std::all_of(edges2.begin(), edges2.end(), [&](EdgeHandle eh) {
        if (!graph.valid(eh)) {
            // removed
            return true;
        } else {
            HalfEdgeHandle heh = graph.halfEdge(eh);
            return !graph.collapseAllowed(heh);
        }
    }));
}

TEST_CASE("collapse simple 2", "[simplify]") {
    Graph graph;
    std::array<Pvl::VertexHandle, 8> vs;
    for (int i = 0; i < 8; ++i) {
        vs[i] = graph.addVertex();
    }
    Pvl::VertexHandle vA = graph.addVertex();
    Pvl::VertexHandle vB = graph.addVertex();

    graph.addFace(vs[0], vA, vs[1]);
    graph.addFace(vA, vs[2], vs[1]);
    graph.addFace(vA, vB, vs[2]);
    graph.addFace(vB, vs[3], vs[2]);
    graph.addFace(vB, vs[4], vs[3]);
    graph.addFace(vs[5], vs[4], vB);
    graph.addFace(vs[6], vs[5], vB);
    graph.addFace(vs[7], vs[6], vA);
    graph.addFace(vs[6], vB, vA);
    graph.addFace(vs[0], vs[7], vA);

    std::cout << "precollapse" << std::endl;
    for (FaceHandle fh : graph.faceRange()) {
        std::cout << "Face " << fh << " - " << std::flush;
        for (VertexHandle vh : graph.vertexRing(fh)) {
            std::cout << vh << " ";
        }
        std::cout << std::endl;
    }

    HalfEdgeHandle eh = graph.halfEdge(vA, vB);
    graph.collapse(eh);

    REQUIRE(graph.valid(vA));
    REQUIRE_FALSE(graph.valid(vB));
    auto range = graph.vertexRing(vA);
    std::set<VertexHandle> ring(range.begin(), range.end());
    REQUIRE(
        ring == std::set<VertexHandle>({ 0_vh, 1_vh, 2_vh, 3_vh, 4_vh, 5_vh, 6_vh, 7_vh }));
    std::cout << "postcollapse " << std::endl;
    for (FaceHandle fh : graph.faceRange()) {
        std::cout << "Face " << fh << " - " << std::flush;
        if (!graph.valid(fh)) {
            std::cout << " removed!" << std::endl;
            continue;
        }
        for (VertexHandle vh : graph.vertexRing(fh)) {
            REQUIRE(graph.valid(vh));
            std::cout << vh << " ";
        }
        for (HalfEdgeHandle heh : graph.halfEdgeRing(fh)) {
            REQUIRE(graph.valid(heh));
        }
        for (FaceHandle f : graph.faceRing(fh)) {
            REQUIRE(graph.valid(f));
        }
        std::cout << std::endl;
    }
}

TEST_CASE("simplify simple", "[simplify]") {
    Pvl::TriangleMesh<Pvl::Vec3f, std::size_t> mesh;
    std::array<Pvl::VertexHandle, 8> vs;
    for (int i = 0; i < 8; ++i) {
        vs[i] = mesh.addVertex();
    }
    Pvl::VertexHandle vA = mesh.addVertex();
    Pvl::VertexHandle vB = mesh.addVertex();
    mesh.points.resize(mesh.numVertices());
    mesh.points[vA] = Pvl::Vec3f(0, 0, 0);
    mesh.points[vB] = Pvl::Vec3f(1, 0, 0);
    mesh.points[vs[0]] = Pvl::Vec3f(-0.5, 0, 0);
    mesh.points[vs[1]] = Pvl::Vec3f(-0.1, 0.5, 0);
    mesh.points[vs[2]] = Pvl::Vec3f(0.5, 0.6, 0);
    mesh.points[vs[3]] = Pvl::Vec3f(1.3, 0.3, 0);
    mesh.points[vs[4]] = Pvl::Vec3f(1.6, 0, 0);
    mesh.points[vs[5]] = Pvl::Vec3f(1.3, -0.3, 0);
    mesh.points[vs[6]] = Pvl::Vec3f(1, -0.8, 0);
    mesh.points[vs[7]] = Pvl::Vec3f(0, -0.9, 0);

    mesh.addFace(vs[0], vA, vs[1]);
    mesh.addFace(vA, vs[2], vs[1]);
    mesh.addFace(vA, vB, vs[2]);
    mesh.addFace(vB, vs[3], vs[2]);
    mesh.addFace(vB, vs[4], vs[3]);
    mesh.addFace(vs[5], vs[4], vB);
    mesh.addFace(vs[6], vs[5], vB);
    mesh.addFace(vs[7], vs[6], vA);
    mesh.addFace(vs[6], vB, vA);
    mesh.addFace(vs[0], vs[7], vA);

    Pvl::HalfEdgeHandle eh = mesh.halfEdge(vA, vB);
    std::cout << "A-B edge = " << eh << std::endl;
    {
        std::ofstream ofs("base.ply");
        Pvl::PlyWriter writer(ofs);
        writer << mesh;
    }
    std::cout << "A ring:" << std::endl;
    for (Pvl::VertexHandle vh : mesh.vertexRing(vA)) {
        std::cout << vh << ",";
    }
    std::cout << std::endl;
    std::cout << "B ring:" << std::endl;
    for (Pvl::VertexHandle vh : mesh.vertexRing(vB)) {
        std::cout << vh << ",";
    }
    std::cout << std::endl;

    mesh.collapse(eh, Pvl::Vec3f(0.5, 0., 0.));
    std::cout << "After collapse A ring:" << std::endl;
    std::cout << "A ring:" << std::endl;
    for (Pvl::VertexHandle vh : mesh.vertexRing(vA)) {
        std::cout << vh << ",";
    }
    std::cout << std::endl;

    {
        std::ofstream ofs("simplified.ply");
        Pvl::PlyWriter writer(ofs);
        writer << mesh;
    }
}

TEST_CASE("Test bunny", "[mesh]") {
    std::ifstream ifs("/home/pavel/projects/pvl/data/bunny-fixed.ply");
    Pvl::PlyReader reader(ifs);
    auto bunny = reader.readMesh();

    auto vertices = bunny.vertexRange();
    auto faces = bunny.faceRange();
    auto halfEdges = bunny.halfEdgeRange();
    auto edges = bunny.edgeRange();

    auto valid = [&](auto h) { return bunny.valid(h); };
    REQUIRE(std::all_of(vertices.begin(), vertices.end(), valid));
    REQUIRE(std::all_of(faces.begin(), faces.end(), valid));
    REQUIRE(std::all_of(halfEdges.begin(), halfEdges.end(), valid));
    REQUIRE(std::all_of(edges.begin(), edges.end(), valid));
}

TEST_CASE("simplify bunny", "[simplify]") {
    std::ifstream ifs("/home/pavel/projects/pvl/data/bunny-fixed.ply");
    Pvl::PlyReader reader(ifs);
    auto bunny = reader.readMesh();

    std::cout << "Simplifying bunny " << std::endl;
    Pvl::simplify(bunny,
        Pvl::EdgeLengthCost<Pvl::Vec3f>{},
        Pvl::MidpointPlacement<Pvl::Vec3f>{},
        Pvl::EdgeCountStop{ 33000 });
    {
        std::ofstream ofs("simplified.ply");
        Pvl::PlyWriter writer(ofs);
        writer << bunny;
    }
}
/*   std::ifstream ifs("/home/pavel/projects/pvl/data/pc.ply");
   Pvl::PlyReader reader(ifs);
   std::vector<Pvl::Vec3f> points;
   {
       auto all = reader.readCloud();
       for (std::size_t i = 0; i < all.size(); i += 10) {
           points.push_back(all[i]);
       }
   }
   std::cout << "Loaded cloud with " << points.size() << " points " << std::endl;
   std::cout << "Estimating normals" << std::endl;
   std::vector<Pvl::Vec3f> normals = Pvl::estimateNormals(points);

   {
       std::ofstream ofs("pc-raw.ply");
       Pvl::PlyWriter writer(ofs);
       writer.write(points, normals);
   }

   std::cout << "Orienting normals" << std::endl;
   Pvl::orientNormals(points, normals);
   {
       std::ofstream ofs("pc-oriented.ply");
       Pvl::PlyWriter writer(ofs);
       writer.write(points, normals);
   }
*/
/*std::ifstream
ifs("/home/pavel/projects/random/mesh/data/bunny/reconstruction/bun_zipper.ply");
Pvl::PlyReader reader(ifs);
auto mesh = reader.readMesh();*/

/* for (int i = 0; i < 200; ++i) {
     std::cout << "iter = " << i << std::endl;
     Pvl::laplacianSmoothing<Pvl::ParallelTag>(mesh);
 }
 std::ofstream ofs("bunny.ply");
 Pvl::PlyWriter writer(ofs);
 writer << mesh;

 std::ofstream ofs2("boundary.ply");
 Pvl::PlyWriter boun(ofs2);
 for (auto vh : mesh.vertexRange()) {
     if (mesh.boundary(vh)) {
         boun << mesh.points[vh];
     }
 }*/

/*std::vector<Pvl::Vec3f> normals = Pvl::estimateNormals(mesh.points);

{
    std::ofstream ofs("bunny-raw.ply");
    Pvl::PlyWriter writer(ofs);
    writer.write(mesh.points, normals);
}

Pvl::orientNormals(mesh.points, normals);
{
    std::ofstream ofs("bunny-oriented.ply");
    Pvl::PlyWriter writer(ofs);
    writer.write(mesh.points, normals);
*/
