#pragma once

#include "TriangleMesh.hpp"
#include "Vector.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace Pvl {

class PlyReader {
    std::istream& in_;

public:
    PlyReader(std::istream& in)
        : in_(in) {}

    std::vector<Vec3f> readCloud() {
        std::string line;
        // std::string countTag("element vertex ");
        while (std::getline(in_, line)) {
            /*   if (line.size() > countTag.size() && line.substr(0, countTag.size()) ==
               countTag) {

               }*/
            if (line == "end_header") {
                break;
            }
        }
        std::vector<Vec3f> points;
        while (std::getline(in_, line)) {
            std::stringstream ss(line);
            Vec3f p;
            ss >> p[0] >> p[1] >> p[2];
            points.push_back(p);
        }
        return points;
    }

    TriangleMesh<Vec3f, int> readMesh() {
        std::string line;
        std::size_t numVertices = 0;
        std::size_t numFaces = 0;
        while (std::getline(in_, line)) {
            sscanf(line.c_str(), "element vertex %zu", &numVertices);
            sscanf(line.c_str(), "element face %zu", &numFaces);
            if (line == "end_header") {
                break;
            }
        }
        TriangleMesh<Vec3f, int> mesh;
        std::cout << "Loading mesh with " << numVertices << " vertices and " << numFaces
                  << " faces" << std::endl;
        for (std::size_t i = 0; i < numVertices; ++i) {
            std::getline(in_, line);
            std::stringstream ss(line);
            Vec3f p;
            ss >> p[0] >> p[1] >> p[2];
            mesh.addVertex();
            mesh.points.push_back(p);
        }
        for (std::size_t i = 0; i < numFaces; ++i) {
            std::getline(in_, line);
            std::stringstream ss(line);
            int dummy, a, b, c;
            ss >> dummy >> a >> b >> c;
            mesh.addFace(VertexHandle(a), VertexHandle(b), VertexHandle(c));
        }
        return mesh;
    }
};

} // namespace Pvl
