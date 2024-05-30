#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "Vector.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <queue>

using namespace std;

class Polygon
{
public:
    vector<Vector> vertices;
    Polygon()
    {
        this->vertices = {};
    }

    Polygon(vector<Vector> &vertices)
    {
        this->vertices = vertices;
    }
};

class Line
{
public:
    Vector u, N;
    Line(Vector u, Vector N)
    {
        N.normalize();
        this->u = u;
        this->N = N;
    }
};

void save_svg(const vector<Polygon> &polygons, string filename, string fillcol = "none")
{
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i = 0; i < polygons.size(); i++)
    {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++)
        {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

void save_svg_w_points(const vector<Polygon> &polygons, const vector<Vector> &points, const string &filename, const string &fillcol = "none", const std::string &pointcol = "red")
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
            fprintf(f, "%3.3f, %3.3f ", (vertex[0] * 1000), (1000 - vertex[1] * 1000));
        }
        fprintf(f, "\" fill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }

    // Draw points
    for (const auto &point : points)
    {
        fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r = \"5\" fill = \"%s\" />\n", (point[0] * 1000), (1000 - point[1] * 1000), pointcol.c_str());
    }

    fprintf(f, "</svg>\n");
    fclose(f);
}
// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const vector<Polygon> &polygons, string filename, int frameid, int nbframes)
{
    FILE *f;
    if (frameid == 0)
    {
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
    }
    else
    {
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (int i = 0; i < polygons.size(); i++)
    {
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++)
        {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
    }
    fprintf(f, "<animate\n");
    fprintf(f, "    id = \"frame%u\"\n", frameid);
    fprintf(f, "    attributeName = \"display\"\n");
    fprintf(f, "    values = \"");
    for (int j = 0; j < nbframes; j++)
    {
        if (frameid == j)
        {
            fprintf(f, "inline");
        }
        else
        {
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n    keyTimes = \"");
    for (int j = 0; j < nbframes; j++)
    {
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n   dur = \"5s\"\n");
    fprintf(f, "    begin = \"0s\"\n");
    fprintf(f, "    repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1)
    {
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}

Vector edgeLineIntersect(Vector &prevVertex, Vector &curVertex, Line clipEdge)
{
    double t = dot(clipEdge.u - prevVertex, clipEdge.N) / dot(curVertex - prevVertex, clipEdge.N);
    if (t < 0 || t > 1)
    {
        // cerr << "No intersection" << endl;
        return Vector(INFINITY, INFINITY, INFINITY);
    }
    return prevVertex + t * (curVertex - prevVertex);
}

Polygon polygonClip(Polygon subjectPolygon, Polygon clipPolygon)
{
    Polygon outPolygon;
    for (int i = 0; i < clipPolygon.vertices.size(); i++)
    {
        outPolygon = Polygon();
        Vector dir, normal;
        if (i < clipPolygon.vertices.size() - 1)
            dir = clipPolygon.vertices[i + 1] - clipPolygon.vertices[i];
        else
            dir = clipPolygon.vertices[0] - clipPolygon.vertices[i];
        normal = Vector();
        normal[0] = dir[1];
        normal[1] = -dir[0];
        normal.normalize();
        Line clipEdge = Line(clipPolygon.vertices[i], normal);

        for (int i = 0; i < subjectPolygon.vertices.size(); i++)
        {
            Vector curVertex = subjectPolygon.vertices[i];
            Vector prevVertex = subjectPolygon.vertices[i > 0 ? i - 1 : subjectPolygon.vertices.size() - 1];
            Vector intersection = edgeLineIntersect(prevVertex, curVertex, clipEdge);
            if (dot(curVertex - clipEdge.u, clipEdge.N) <= 0)
            {
                if (dot(prevVertex - clipEdge.u, clipEdge.N) > 0)
                {
                    outPolygon.vertices.push_back(intersection);
                }
                outPolygon.vertices.push_back(curVertex);
            }
            else if (dot(prevVertex - clipEdge.u, clipEdge.N) <= 0)
            {
                outPolygon.vertices.push_back(intersection);
            }
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
}

Polygon VoronoiPolygon(vector<Vector> &pointSet, int index)
{
    vector<Vector> boundVertices = {Vector(0, 0, 0),
                                    Vector(0, 1, 0),
                                    Vector(1, 1, 0),
                                    Vector(1, 0, 0)};
    Polygon subjectPolygon = Polygon(boundVertices);
    Polygon outPolygon;

    auto maxDistComp = [pointSet, index](Vector a, Vector b)
    {
        return (a - pointSet[index]).norm2() > (b - pointSet[index]).norm2();
    };

    for (int i = 0; i < pointSet.size(); i++)
    {
        if (i == index)
        {
            continue;
        }

        outPolygon = Polygon();
        Vector M = (pointSet[i] + pointSet[index]) / 2.;
        Line clipEdge = Line(M, pointSet[i] - pointSet[index]);

        for (int i = 0; i < subjectPolygon.vertices.size(); i++)
        {
            Vector curVertex = subjectPolygon.vertices[i];
            Vector prevVertex = subjectPolygon.vertices[i > 0 ? i - 1 : subjectPolygon.vertices.size() - 1];
            Vector intersection = edgeLineIntersect(prevVertex, curVertex, clipEdge);
            if (dot(curVertex - clipEdge.u, clipEdge.N) < 0)
            {
                if (dot(prevVertex - clipEdge.u, clipEdge.N) >= 0)
                {
                    outPolygon.vertices.push_back(intersection);
                }
                outPolygon.vertices.push_back(curVertex);
            }
            else if (dot(prevVertex - clipEdge.u, clipEdge.N) < 0)
            {
                outPolygon.vertices.push_back(intersection);
            }
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
}

int main()
{
    vector<Vector> sPol_vertices = {
        Vector(60. / 200, 20. / 200, 0),
        Vector(120. / 200, 30. / 200, 0),
        Vector(140. / 200, 80. / 200, 0),
        Vector(110. / 200, 150. / 200, 0),
        Vector(50. / 200, 130. / 200, 0),
        Vector(30. / 200, 60. / 200, 0),
        Vector(70. / 200, 100. / 200, 0)};
    // vector<Vector> sPol_vertices = {
    //     Vector(100. / 200, 190. / 200, 0),
    //     Vector(120. / 200, 130. / 200, 0),
    //     Vector(180. / 200, 130. / 200, 0),
    //     Vector(140. / 200, 90. / 200, 0),
    //     Vector(160. / 200, 30. / 200, 0),
    //     Vector(100. / 200, 70. / 200, 0),
    //     Vector(40. / 200, 30. / 200, 0),
    //     Vector(60. / 200, 90. / 200, 0),
    //     Vector(20. / 200, 130. / 200, 0),
    //     Vector(80. / 200, 130. / 200, 0)};

    vector<Vector> clipPol_vertices = {
        Vector(40. / 200, 40. / 200, 0),
        Vector(140. / 200, 40. / 200, 0),
        Vector(140. / 200, 140. / 200, 0),
        Vector(40. / 200, 140. / 200, 0)};

    Polygon subject = Polygon(sPol_vertices);
    // Polygon clip = Polygon(clipPol_vertices);
    // save_svg({subject, clip}, "test.svg");
    // Polygon outPolygon = polygonClip(subject, clip);

    vector<Polygon> voronoi;
    for (int i = 0; i < sPol_vertices.size(); i++)
    {
        voronoi.push_back(VoronoiPolygon(sPol_vertices, i));
    }
    // save_svg({subject, clip}, "test.svg");
    save_svg_w_points(voronoi, sPol_vertices, "voronoi.svg");
}