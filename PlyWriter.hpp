#pragma once

#include "TriangleMesh.hpp"
#include "Vector.hpp"
#include <fstream>

namespace Pvl {
class PlyWriter {
    std::ostream& out_;
    std::size_t vertexPos_;
    std::size_t facePos_;

    bool closed_ = false;
    std::size_t vertexCnt_ = 0;
    std::size_t faceCnt_ = 0;

public:
    PlyWriter(std::ostream& out)
        : out_(out) {
        out << "ply\n";
        out << "format ascii 1.0\n";
        out << "comment Created by PVL library\n";
        out << "element vertex ";
        vertexPos_ = out.tellp();
        out << "              \n";
        out << "property float x\n";
        out << "property float y\n";
        out << "property float z\n";
        out << "property float nx\n";
        out << "property float ny\n";
        out << "property float nz\n"; /// \todo
        out << "element face ";
        facePos_ = out.tellp();
        out << "              \n";
        out << "property list uchar int vertex_index\n";
        out << "end_header\n";
    }

    template <typename T, int Dim>
    PlyWriter& operator<<(const Vector<T, Dim>& p) {
        for (int i = 0; i < Dim; ++i) {
            out_ << p[i] << " ";
        }
        if (Dim == 3) {
            out_ << "0 0 0\n";
        }
        ++vertexCnt_;
        return *this;
    }

    template <typename T, int Dim>
    PlyWriter& operator<<(const std::vector<Vector<T, Dim>>& cloud) {
        for (auto& p : cloud) {
            *this << p;
        }
        return *this;
    }

    template <typename T, int Dim>
    PlyWriter& write(const std::vector<Vector<T, Dim>>& cloud,
        const std::vector<Vector<T, Dim>>& normals) {
        for (std::size_t i = 0; i < cloud.size(); ++i) {
            const auto& p = cloud[i];
            const auto& n = normals[i];
            out_ << p[0] << " " << p[1] << " " << p[2] << " " << n[0] << " " << n[1] << " "
                 << n[2] << "\n";
        }
        vertexCnt_ += cloud.size();
        return *this;
    }

    template <typename Vec, typename Index>
    PlyWriter& operator<<(const TriangleMesh<Vec, Index>& mesh) {
        for (Index i = 0; i < Index(mesh.numVertices()); ++i) {
            const Vec& p = mesh.points[i];
            out_ << p[0] << " " << p[1] << " " << p[2] << " 0 0 0\n";
        }
        std::size_t validCnt = 0;
        for (Index i = 0; i < Index(mesh.numFaces()); ++i) {
            if (!mesh.valid(FaceHandle(i))) {
                continue;
            }
            auto f = mesh.faceVertices(FaceHandle(i));
            out_ << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
            validCnt++;
        }
        vertexCnt_ += mesh.numVertices();
        faceCnt_ += validCnt;
        return *this;
    }

    void close() {
        if (!closed_) {
            out_.seekp(vertexPos_);
            out_ << vertexCnt_;
            out_.seekp(facePos_);
            out_ << faceCnt_;
        }
        closed_ = true;
    }

    ~PlyWriter() {
        close();
    }
};
} // namespace Pvl
