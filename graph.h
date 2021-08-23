#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <queue>
#include <random>
#include <ctime>
#include <set>
using namespace std;

struct Color {
    double red = 0, green = 0, blue = 0;
    explicit Color(double red, double green, double blue) : red(red), green(green), blue(blue){

    }
};

const long long MAX_COST = 1000000000;
const long long MIN_COST = -MAX_COST;
const double STANDART_RADIUS = 15;
class Vertex {
private:
    double x;
    double y;
    double radius;
    long long distanse;
    int num = 0;
    bool is_active = false;
    bool visible = true;
public:
    Color color;
    Vertex(double x, double y, double radius = STANDART_RADIUS);
    Vertex(double x, double y, double radius, Color color);
    void moveTo(double x, double y);
    void moveByVector(double x, double y);
    void activate();
    void deactivate();
    void setDist(long long dist);
    void setVisible(bool vis);
    bool operator== (const Vertex& other) const;
    bool operator!= (const Vertex& other) const;
    double getX() const;
    double getY() const;
    int getNum() const;
    long long getDist() const;
    double getRadius() const;
    bool isActive() const;
    bool isVisible() const;
    bool intersect(double x, double y) const;
    void setNum(int nm);
};


class Graph {
private:
    bool numbered = true;
    bool oriented = false;
    void resize(int size);
    void renumerate();
    void dfsAlgo(Vertex* X, vector<bool>& used, vector<pair<Vertex*, Vertex*> >& r);
    void bfsAlgo(Vertex* X, vector<bool>& used, vector<pair<Vertex*, Vertex*> >& r);
    void setCost(int x, int y, long long cost);
    void addEdge(int XNum, int YNum, long long cost = 0);
    int getSnmPr(int X, vector<int>& pr);
    Vertex* randVertex(int width, int height) const;
    long long rnd() const;
public:
    int numeration = 0;
    Graph(int size);
    Graph(int size, bool oriented);
    vector<Vertex*> vertexes;
    vector<vector<int> > edges;
    vector<vector<long long> > costs;
    void setNumbered(bool numbered);
    void setNumeration(int numeration);
    void deleteVertex(int num);
    void deleteVertex(const Vertex* X);
    void addVertex(Vertex* X);
    void addEdge(const Vertex* X, const Vertex* Y, long long cost = 0);
    void deleteEdge(int x, int y);
    void deleteEdge(const Vertex* X, const Vertex* Y);
    void deactivateAllVertex();
    void setVisible(bool vis);
    void setCost(const Vertex* X, const Vertex* Y, long long cost);
    int outDegree(const Vertex* X) const;
    int inDegree(const Vertex* X)const;
    int getNumeration() const;
    bool isNumerated() const;
    bool isOriented() const;
    bool areConnected(int x, int y) const;
    bool areConnected(const Vertex* X, const Vertex* Y) const;
    int getCost(int x, int y) const;
    int getCost(const Vertex* X, const Vertex* Y) const;
    void generateComponent(int size, int width, int height, double probality, int maxCost = 0);
    vector<pair<Vertex*, Vertex*> > dfs(Vertex* X);
    vector<pair<Vertex*, Vertex*> > bfs(Vertex* X);
    vector<pair<pair<Vertex*, Vertex*>, long long> > kruskal(Vertex* X);
    vector<pair<pair<Vertex*, Vertex*>, long long> > prima(Vertex* X);
    vector<pair<pair<Vertex*, Vertex*>, long long> > dijkstra(Vertex* X);
    vector<pair<pair<Vertex*, Vertex*>, long long> > bellman_ford(Vertex* X);
    vector<pair<pair<Vertex*, Vertex*>, long long> > floyd_warshall(Vertex* X);
    vector<pair<pair<Vertex*, Vertex*>, long long> > jonson(Vertex* X);
};
