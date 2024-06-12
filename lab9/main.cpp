#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "Vector.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <omp.h>
#include <chrono>
using namespace std;

#define PRINT_VEC(v) (printf("%s: (%f %f %f)\n", #v, (v)[0], (v)[1], (v)[2]))

class TriangleIndices
{
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group){};
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;    // indices within the uv coordinates array
    int ni, nj, nk;       // indices within the normals array
    int group;            // face group
};

class TriangleMesh
{
public:
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;

    TriangleMesh(){};
    ~TriangleMesh() {}

    void readOBJ(const char *obj)
    {

        char matfile[255];
        char grp[255];

        FILE *f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f))
        {
            char line[255];
            if (!fgets(line, 255, f))
                break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's')
            {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ')
            {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
                {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
                }
                else
                {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n')
            {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't')
            {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f')
            {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char *consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    if (k0 < 0)
                        t.ni = normals.size() + k0;
                    else
                        t.ni = k0 - 1;
                    if (k1 < 0)
                        t.nj = normals.size() + k1;
                    else
                        t.nj = k1 - 1;
                    if (k2 < 0)
                        t.nk = normals.size() + k2;
                    else
                        t.nk = k2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (j0 < 0)
                            t.uvi = uvs.size() + j0;
                        else
                            t.uvi = j0 - 1;
                        if (j1 < 0)
                            t.uvj = uvs.size() + j1;
                        else
                            t.uvj = j1 - 1;
                        if (j2 < 0)
                            t.uvk = uvs.size() + j2;
                        else
                            t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3)
                        {
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            if (k0 < 0)
                                t.ni = normals.size() + k0;
                            else
                                t.ni = k0 - 1;
                            if (k1 < 0)
                                t.nj = normals.size() + k1;
                            else
                                t.nj = k1 - 1;
                            if (k2 < 0)
                                t.nk = normals.size() + k2;
                            else
                                t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true)
                {
                    if (consumedline[0] == '\n')
                        break;
                    if (consumedline[0] == '\0')
                        break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        if (k0 < 0)
                            t2.ni = normals.size() + k0;
                        else
                            t2.ni = k0 - 1;
                        if (k2 < 0)
                            t2.nj = normals.size() + k2;
                        else
                            t2.nj = k2 - 1;
                        if (k3 < 0)
                            t2.nk = normals.size() + k3;
                        else
                            t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (j0 < 0)
                                t2.uvi = uvs.size() + j0;
                            else
                                t2.uvi = j0 - 1;
                            if (j2 < 0)
                                t2.uvj = uvs.size() + j2;
                            else
                                t2.uvj = j2 - 1;
                            if (j3 < 0)
                                t2.uvk = uvs.size() + j3;
                            else
                                t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                if (k0 < 0)
                                    t2.ni = normals.size() + k0;
                                else
                                    t2.ni = k0 - 1;
                                if (k2 < 0)
                                    t2.nj = normals.size() + k2;
                                else
                                    t2.nj = k2 - 1;
                                if (k3 < 0)
                                    t2.nk = normals.size() + k3;
                                else
                                    t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1)
                                {
                                    if (i0 < 0)
                                        t2.vtxi = vertices.size() + i0;
                                    else
                                        t2.vtxi = i0 - 1;
                                    if (i2 < 0)
                                        t2.vtxj = vertices.size() + i2;
                                    else
                                        t2.vtxj = i2 - 1;
                                    if (i3 < 0)
                                        t2.vtxk = vertices.size() + i3;
                                    else
                                        t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else
                                {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(f);
    }
};

pair<int, int> construct_pair(int a, int b)
{

    return a < b ? make_pair(a, b) : make_pair(b, a);
}

class Polygon
{
public:
    vector<Vector> vertices;
    Polygon(){};
};

void save_svg_w_points(const vector<Polygon> &polygons, const vector<Vector> &points, const string &filename, const string &fillcol = "none", const std::string &pointcol = "red");

int main(int argc, char **argv)
{
    TriangleMesh mesh = TriangleMesh();
    mesh.readOBJ("goethe.obj");
    printf("Vertices size: %d\n", mesh.vertices.size());

    auto hash_pair = [](const pair<int, int> a)
    {
        return hash<int>{}(a.first) ^ hash<int>{}(a.second);
    };

    unordered_map<pair<int, int>, int, decltype(hash_pair)> edge_count(0, hash_pair);

    for (auto tid : mesh.indices)
    {
        pair<int, int> A, B, C;
        A = construct_pair(tid.vtxi, tid.vtxj);
        B = construct_pair(tid.vtxj, tid.vtxk);
        C = construct_pair(tid.vtxi, tid.vtxk);
        edge_count[A] += 1;
        edge_count[B] += 1;
        edge_count[C] += 1;

        // printf("%d %d %d\n", edge_count[A], edge_count[B], edge_count[C]);
        // break;
    }
    // printf("edge count: %d\n", edge_count.size());

    vector<pair<int, int>> bounding_egdes;
    int starting = -1, current = -1, next;
    unordered_map<int, vector<int>> graph;
    int n = 0;
    unordered_map<int, vector<int>> adj_list;
    for (auto itr = edge_count.begin(); itr != edge_count.end(); itr++)
    {
        pair<int, int> edge = itr->first;
        if (itr->second == 1)
        {
            // cout << "here\n";
            if (starting < 0)
            {
                starting = edge.first;
            }
            n += 1;
            // bounding_egdes.push_back(itr->first);
            graph[edge.first].push_back(edge.second);
            graph[edge.second].push_back(edge.first);
        }
        adj_list[edge.first].push_back(edge.second);
        adj_list[edge.second].push_back(edge.first);
    }

    double s = 0.;

    vector<int> boundaries;
    current = starting;
    cout << current << " " << graph[current].size() << endl;
    next = graph[current][0];
    while (next != starting)
    {
        boundaries.push_back(current);
        s += (mesh.vertices[next] - mesh.vertices[current]).norm();
        int tmp = next;
        // cout << graph[next].size() << endl;
        if (graph[next][0] != current)
        {
            next = graph[next][0];
        }
        else
        {
            next = graph[next][1];
        }
        current = tmp;
    }
    s += (mesh.vertices[next] - mesh.vertices[current]).norm();

    double cs = 0.;
    vector<Vector> v = mesh.vertices;
    double theta;
    vector<int> interior;
    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        interior.push_back(i);
    }

    Polygon disk;
    vector<Vector> vertices = mesh.vertices;
    for (int i = 0; i < boundaries.size(); i++)
    {
        theta = 2 * M_PI * cs / s;
        int next = boundaries[i < boundaries.size() - 1 ? i + 1 : 0];
        int cur = boundaries[i];
        for (auto itr = interior.begin(); itr < interior.end(); itr++)
        {
            if (*itr == cur)
            {
                interior.erase(itr);
            }
        }

        vertices[cur] = Vector(cos(theta), sin(theta), 0);
        disk.vertices.push_back(vertices[cur]);
        cs += (mesh.vertices[next] - mesh.vertices[cur]).norm();
    }
    auto start = std::chrono::high_resolution_clock::now();
    int num_itr = 20000;
    for (int i = 0; i < num_itr; i++)
    {
        vector<Vector> new_v = vertices;

        // #pragma omp parallel for schedule(dynamic, 1)
        for (int idx = 0; idx < interior.size(); idx++)
        {
            int j = interior[idx];
            Vector sum(0, 0, 0);
            if (!adj_list[j].size())
            {
                printf("zero\n");
                continue;
            }

            for (auto k : adj_list[j])
            {
                sum = sum + vertices[k];
            }
            new_v[j] = sum / adj_list[j].size();
        }
        // if (i % 1000 == 0)
        // {
        //     PRINT_VEC(new_v[0]);
        // }
        vertices = new_v;
    }

    vector<Polygon> scene;
    for (auto tid : mesh.indices)
    {
        Polygon t;
        t.vertices = {vertices[tid.vtxi], vertices[tid.vtxj], vertices[tid.vtxk]};
        scene.push_back(t);
    }

    scene.push_back(disk);

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken for embedding: " << elapsed.count() << " seconds" << std::endl;

    save_svg_w_points(scene, disk.vertices, "image_1.svg");
}

void save_svg_w_points(const vector<Polygon> &polygons, const vector<Vector> &points, const string &filename, const string &fillcol, const std::string &pointcol)
{
    FILE *f = fopen(filename.c_str(), "w+");
    if (!f)
    {
        fprintf(stderr, "Error opening file %s for writing\n", filename.c_str());
        return;
    }

    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");

    // Draw polygons
    for (const auto &polygon : polygons)
    {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (const auto &vertex : polygon.vertices)
        {
            fprintf(f, "%3.3f, %3.3f ", (vertex[0] * 500 + 500), (500 - vertex[1] * 500));
        }
        fprintf(f, "\" fill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }

    // Draw points
    for (const auto &point : points)
    {
        fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r = \"5\" fill = \"%s\" />\n", (point[0] * 500 + 500), (500 - point[1] * 500), pointcol.c_str());
    }

    fprintf(f, "</svg>\n");
    fclose(f);
}