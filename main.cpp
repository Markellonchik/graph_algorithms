#include <GL/glut.h>
#include "graph.h"
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <utility>
#include <map>

using namespace std;

const long double eps = 0.000000001;

struct Point {
    double x, y;
};
struct Polygon {
    vector<Point> points;
    bool isActive = true;
};
double dist(Point X, Point Y) {
    return sqrt(pow(X.x - Y.x, 2) + pow(X.y - Y.y, 2));
}
double trianglearea(double x1, double y1, double x2, double y2, double x3, double y3) {
    return abs( x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2) ) / 2;
}
double area(const Polygon& X)
{
    double r = 0;
    for(size_t i = 0; i < X.points.size(); ++i) {
        const Point& f = X.points[i];
        const Point& s = (i + 1 < X.points.size()) ? X.points[i + 1] : X.points[0];

        r += f.x * s.y - f.y * s.x;
    }
    r = abs(r / 2);
    return r;
}
double area(const Polygon& X, const Point& Y) {
    double r = 0;
    for(size_t i = 0; i < X.points.size(); ++i) {
        const Point& f = X.points[i];
        const Point& s = (i + 1 < X.points.size()) ? X.points[i + 1] : X.points[0];
        r += trianglearea(f.x, f.y, s.x, s.y, Y.x, Y.y);
    }
    return r;
}
bool intersect(const Polygon& X, const Point& Y) {
    return abs(area(X) - area(X, Y)) < eps;
}

const string DFS = "Run depth-first search";
const string BFS = "Run breadth-first search";
const string KRUSKAL = "Run Kruskal's algorithm";
const string PRIMA = "Run Prim's algorithm";
const string DIJKSTRA = "Run Dijkstra's algorithm";
const string BELLMAN_FORD = "Run Bellman-Ford algorithm";
const string FLOYD_WARSHALL = "Run Floyd–Warshall algorithm";
Vertex *lastAlgoX, *lastAlgoY;
vector<long long> lastDist;
vector<pair<Vertex*, Vertex*> > algo;
vector<pair<pair<Vertex*, Vertex*>, long long> > algoSnmDist;
map<pair<Vertex*, Vertex*>, bool> shownEd;
string algoName, algoDescription;
const Color ALGO_DESCRIPTION_COLOR(0.0, 0.0, 0.0);
size_t algoPos;

const double WINDOW_WIDTH = 1300, WINDOW_HEIGHT = 700;
const int TIMER_DELAY = 1000;
const int CIRCLE_ACCUR = 40;
const Color edgesColor(0.8, 0, 0.8);
const Color numbersColor(1.0, 1.0, 1.0);

const Color DISTANSE_COLOR(0.0, 0.0, 1.0);
const Color ALGORITHM_SHOW_EDGE(1337.0 / 1488, 322.0 / 1000, 228.0 / 255);
const Color ALGORITHM_SHOW_VERTEX(0.0, 1.0, 0.0);
const Color activeVertexColor(1.0, 0.0, 0.0);
const double ACTIVE_CIRCLE_WIDTH = 5;

const Color EDGES_COST_COLOR(0.0, 0.0, 0.0);

const Color interectedVertexColor(1.0, 1.0, 0.0);
const double INTERSECTED_CIRCLE_WIDTH = 5;
int mouse_button_clicked;
int lastX, lastY;
Polygon* active_polygon;
const Color ACTIVE_INFORMATION_LINES(1.0, 0.0, 0.0);
const Color INFORMATION_LINES(1.0, 1.0, 1.0);

const double INFORMATION_WIDTH = 250, INFORMATION_HEIGHT = 600;
const Color INFORMATION_COLOR(1.0, 0.5, 0.0);
vector<pair<Polygon, string> > Information;
Vertex* activeVertex;
Vertex* intersectVertex;
Vertex *lastEdgeX, *lastEdgeY;

bool distanse_shown = false;
bool is_algorithm_shown = false;
bool isWeighed = true;

# define M_PI           3.14159265358979323846

Graph graph(1, false);

void setColor(const Color& color) {
    glColor3f(color.red, color.green, color.blue);
}

void PrintFilledCircle(double x, double y, double r) {
    glBegin(GL_POLYGON);
    for(double i = 0; i < 2 * M_PI; i += M_PI / CIRCLE_ACCUR) {
        glVertex2f(x + cos(i) * r, y + sin(i) * r);
    }
    glEnd();
}
void PrintCircle(double x, double y, double r) {
    glBegin(GL_LINES);
    for(double i = 0; i < 2 * M_PI; i += M_PI / CIRCLE_ACCUR) {
        glVertex2f(x + cos(i) * r, y + sin(i) * r);
        glVertex2f(x + cos(i + M_PI / CIRCLE_ACCUR) * r, y + sin(i + M_PI / CIRCLE_ACCUR) * r);
    }
    glEnd();
}

void RenderNumbers(double x, double y, void *font, string s, Color color = numbersColor)
{
    setColor(color);
    glRasterPos2f(x, y);

    for(char c : s)
        glutBitmapCharacter(font, c);
}

void PrintAlgorithmDescription() {
    int textLength = algoDescription.size();
    RenderNumbers(-textLength * 6, WINDOW_HEIGHT / 2 - 50, GLUT_BITMAP_TIMES_ROMAN_24, algoDescription, ALGO_DESCRIPTION_COLOR);
}

void PrintPolygonLines(const Polygon& X) {
    glBegin(GL_LINES);
    for(size_t i = 0; i < X.points.size(); ++i) {
        const Point& A = X.points[i];
        const Point& B = i + 1 < X.points.size() ? X.points[i + 1] : X.points[0];
        glVertex2f(A.x, A.y);
        glVertex2f(B.x, B.y);
    }
    glEnd();
}


void PrintActiveVertexInformation() {
    if(activeVertex && !is_algorithm_shown) {
        setColor(INFORMATION_COLOR);
        glBegin(GL_POLYGON);
            glVertex2f(WINDOW_WIDTH / 2, WINDOW_HEIGHT / 2);
            glVertex2f(WINDOW_WIDTH / 2, WINDOW_HEIGHT / 2 - INFORMATION_HEIGHT);
            glVertex2f(WINDOW_WIDTH / 2 - INFORMATION_WIDTH, WINDOW_HEIGHT / 2 - INFORMATION_HEIGHT);
            glVertex2f(WINDOW_WIDTH / 2 - INFORMATION_WIDTH, WINDOW_HEIGHT / 2);
        glEnd();

        double x = WINDOW_WIDTH / 2 - INFORMATION_WIDTH / 2;
        double y = WINDOW_HEIGHT / 2 - activeVertex->getRadius() * 2;
        setColor(activeVertexColor);
            PrintFilledCircle(x, y, activeVertex->getRadius() + ACTIVE_CIRCLE_WIDTH);
        glColor3f(activeVertex->color.red, activeVertex->color.green, activeVertex->color.blue);

        PrintFilledCircle(x, y, activeVertex->getRadius());
        stringstream out;
        out << activeVertex->getNum() + graph.getNumeration();
        string text = out.str();
            RenderNumbers(x - 4 - static_cast<int>(text.size()) * 3, y - 5 - static_cast<int>(text.size()) * 2, GLUT_BITMAP_TIMES_ROMAN_24, out.str());

        setColor(INFORMATION_LINES);
        for(const pair<Polygon, string>& X : Information) {
            if(&X.first == active_polygon) setColor(ACTIVE_INFORMATION_LINES);
            PrintPolygonLines(X.first);
            if(&X.first == active_polygon) setColor(INFORMATION_LINES);
            RenderNumbers(X.first.points[0].x + 5, X.first.points[0].y / 2 + X.first.points[1].y / 2 - 5, GLUT_BITMAP_HELVETICA_18, X.second);
        }
    }

}

bool shownEdge(Vertex* X, Vertex* Y) {
    if(algoName == DFS || algoName == BFS) return (X == lastAlgoX && Y == lastAlgoY) ||
                                                  (Y == lastAlgoX && X == lastAlgoY);
    if(algoName == KRUSKAL || algoName == PRIMA) {
        return shownEd[{X, Y}] || shownEd[{Y, X}];
    }
    if(algoName == DIJKSTRA || algoName == BELLMAN_FORD) return (X == lastAlgoX && Y == lastAlgoY) ||
                                    (Y == lastAlgoX && X == lastAlgoY);
    return false;
}

void printEdges() {
    setColor(edgesColor);
    glLineWidth(3.0f);
    glBegin(GL_LINES);
    for(size_t i = 0; i < graph.edges.size(); ++i) {
        if(graph.vertexes[i]->isVisible())
        for(size_t k = 0; k < graph.edges[i].size(); ++k) {
            size_t j = graph.edges[i][k];
            if(graph.vertexes[j]->isVisible()) {
                Vertex& X = *graph.vertexes[i];
                Vertex& Y = *graph.vertexes[j];
                if(is_algorithm_shown && shownEdge(&X, &Y)) setColor(ALGORITHM_SHOW_EDGE);
                glVertex2f(X.getX(), X.getY());
                glVertex2f(Y.getX(), Y.getY());
                if(is_algorithm_shown && shownEdge(&X, &Y)) setColor(edgesColor);
            }
        }
    }
    glEnd();
    if(isWeighed) {
        double edgeLenX, edgeLenY;
        for(size_t i = 0; i < graph.edges.size(); ++i) {
            if(graph.vertexes[i]->isVisible())
            for(size_t k = 0; k < graph.edges[i].size(); ++k) {
                size_t j = graph.edges[i][k];
                int cost = graph.costs[i][k];
                if(graph.vertexes[j]->isVisible()) {
                    Vertex& X = *graph.vertexes[i];
                    Vertex& Y = *graph.vertexes[j];
                    edgeLenX = abs(X.getX() - Y.getX());
                    edgeLenY = abs(X.getY() - Y.getY());
                    stringstream out;
                    out << cost;
                    int textLen = out.str().size();
                    if(edgeLenX < textLen * 14 + STANDART_RADIUS * 2 && edgeLenY < 30 + STANDART_RADIUS * 2)continue;
                    RenderNumbers((X.getX() + Y.getX()) / 2 - textLen * 6, (X.getY() + Y.getY()) / 2 - 6, GLUT_BITMAP_TIMES_ROMAN_24, out.str(), EDGES_COST_COLOR);
                }
            }
        }
    }
}

bool shownVertex(Vertex* X) {
    return X == lastAlgoX;
}
bool shownVertexY(Vertex* Y) {
    return Y == lastAlgoY;
}

void printVertexes() {
    for(Vertex* V : graph.vertexes) {
        if(!V->isVisible())continue;
        if(intersectVertex && *V == *intersectVertex) {
            setColor(interectedVertexColor);
            PrintFilledCircle(V->getX(), V->getY(), V->getRadius() + INTERSECTED_CIRCLE_WIDTH);
        }
        if(V->isActive()) {
            setColor(activeVertexColor);
            PrintFilledCircle(V->getX(), V->getY(), V->getRadius() + ACTIVE_CIRCLE_WIDTH);
        }
        if(is_algorithm_shown && shownVertex(V)) {
            setColor(ALGORITHM_SHOW_VERTEX);
            PrintFilledCircle(V->getX(), V->getY(), V->getRadius() + INTERSECTED_CIRCLE_WIDTH);
        }

        glColor3f(V->color.red, V->color.green, V->color.blue);

        PrintFilledCircle(V->getX(), V->getY(), V->getRadius());

        if(graph.isNumerated()) {
            stringstream out;
            out << V->getNum() + graph.getNumeration();
            string text = out.str();
            RenderNumbers(V->getX() - 4 - static_cast<int>(text.size()) * 3, V->getY() - 5 - static_cast<int>(text.size()) * 2, GLUT_BITMAP_TIMES_ROMAN_24, out.str());
        }
        if(is_algorithm_shown && distanse_shown) {
            stringstream out;
            if(V->getDist() >= MAX_COST * algoSnmDist.size()) out << "inf";
            else out << V->getDist();
            string text = out.str();
            RenderNumbers(V->getX() - 4 - static_cast<int>(text.size()) * 3, V->getY() + V->getRadius(), GLUT_BITMAP_TIMES_ROMAN_24, text, shownVertexY(V) ? ALGORITHM_SHOW_EDGE : DISTANSE_COLOR);
        }
    }
}

void Display() {
    glClear(GL_COLOR_BUFFER_BIT);
    printEdges();
    printVertexes();
    if(is_algorithm_shown) {
        PrintAlgorithmDescription();
    }
    PrintActiveVertexInformation();
    glutSwapBuffers();
}

void Timer(int deep) {
    glutPostRedisplay();
    glutTimerFunc(TIMER_DELAY, Timer, 0);
}

void deactivateAllVertexes() {
    activeVertex = nullptr;
    graph.deactivateAllVertex();
}

const double POLYGON_DISTANSE = 5;

void setActive(Vertex& act) {
    activeVertex = &act;
    Information.clear();

    double x = WINDOW_WIDTH / 2 - INFORMATION_WIDTH;
    double y = WINDOW_HEIGHT / 2 - activeVertex->getRadius() * 2;
    stringstream out;
    y -= INFORMATION_HEIGHT / 13 + 10;

    out.str(string());

    if(!graph.isOriented()) {
        out << "Vertex degree: " << graph.outDegree(activeVertex);
        Polygon pol;
        pol.isActive = false;
        pol.points = {{x, y}, {x, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y}};
        Information.push_back({pol, out.str()});
        y -= INFORMATION_HEIGHT / 10 + POLYGON_DISTANSE;
        out.str(string());
    } else {
        out << "Out degree: " << graph.outDegree(activeVertex);
        Polygon pol;
        pol.isActive = false;
        pol.points = {{x, y}, {x, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y}};
        Information.push_back({pol, out.str()});
        y -= INFORMATION_HEIGHT / 10 + POLYGON_DISTANSE;
        out.str(string());
        pol.points = {{x, y}, {x, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y}};
        Information.push_back({pol, out.str()});
        y -= INFORMATION_HEIGHT / 10 + POLYGON_DISTANSE;
        out.str(string());
    }

    Polygon pol;

    pol.points = {{x, y}, {x, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y}};
    Information.push_back({pol, DFS});
    y -= INFORMATION_HEIGHT / 10 + POLYGON_DISTANSE;

    pol.points = {{x, y}, {x, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y}};
    Information.push_back({pol, BFS});
    y -= INFORMATION_HEIGHT / 10 + POLYGON_DISTANSE;

    if(isWeighed) {
        pol.points = {{x, y}, {x, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y}};
        Information.push_back({pol, KRUSKAL});
        y -= INFORMATION_HEIGHT / 10 + POLYGON_DISTANSE;

        pol.points = {{x, y}, {x, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y}};
        Information.push_back({pol, PRIMA});
        y -= INFORMATION_HEIGHT / 10 + POLYGON_DISTANSE;

        pol.points = {{x, y}, {x, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y}};
        Information.push_back({pol, BELLMAN_FORD});
        y -= INFORMATION_HEIGHT / 10 + POLYGON_DISTANSE;

        pol.points = {{x, y}, {x, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y}};
        Information.push_back({pol, DIJKSTRA});
        y -= INFORMATION_HEIGHT / 10 + POLYGON_DISTANSE;

        pol.points = {{x, y}, {x, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y - INFORMATION_HEIGHT / 10}, {x + INFORMATION_WIDTH, y}};
        Information.push_back({pol, FLOYD_WARSHALL});
        y -= INFORMATION_HEIGHT / 10 + POLYGON_DISTANSE;
    }
}

void setAlgorithmShow(bool show) {
    lastAlgoX = lastAlgoY = nullptr;

    if(!show)shownEd.clear(), distanse_shown = false, lastDist.clear();
    is_algorithm_shown = show;
    graph.setVisible(!show);
}
void nextAlgorithmStep() {
    if(algoName == DFS || algoName == BFS) {
        if(algoPos < algo.size()) {
            lastAlgoX = algo[algoPos].first;
            lastAlgoY = algo[algoPos].second;
            algo[algoPos++].second->setVisible(true);
            stringstream out;
            if(algoName == DFS)
                out << "DFS algorithm. Component consists of " << algoPos << " vertexes";
            else out << "BFS algorithm. Component consists of " << algoPos << " vertexes";
            algoDescription = out.str();
        } else {
            if(algoPos == algo.size()) {
                stringstream out;
                if(algoName == DFS)out << "DFS algorithm ended. Finally, component consists of " << algo.size() << " vertexes";
                else out << "BFS algorithm ended. Finally, component consists of " << algo.size() << " vertexes";
                algoDescription = out.str();
                algoPos++;
            } else {
                setAlgorithmShow(false);
            }
        }
    } else if(algoName == KRUSKAL || algoName == PRIMA) {
        if(algoPos < algoSnmDist.size()) {
            lastAlgoX = algoSnmDist[algoPos].first.first;
            lastAlgoY = algoSnmDist[algoPos].first.second;
            shownEd[{lastAlgoX, lastAlgoY}] = true;
            stringstream out;
            if(algoName == KRUSKAL)out << "Kruskal's algorithm. Cost of the minimum-spanning-tree is " <<  algoSnmDist[algoPos++].second;
            else                 out << "Prim's algorithm. Cost of the minimum-spanning-tree is " <<  algoSnmDist[algoPos++].second;
            algoDescription = out.str();
        } else {
            if(algoPos == algoSnmDist.size()) {
                stringstream out;
                if(algoName == KRUSKAL) out << "Kruskal's algorithm ended. Finally, cost of the minimum-spanning-tree is " << algoSnmDist.back().second;
                else                  out << "Prim's algorithm ended. Finally, cost of the minimum-spanning-tree is " << algoSnmDist.back().second;
                algoDescription = out.str();
                algoPos++;
            } else {
                setAlgorithmShow(false);
            }
        }
    } else if(algoName == DIJKSTRA || algoName == BELLMAN_FORD || algoName == FLOYD_WARSHALL) {
        if(algoPos < algoSnmDist.size()) {
            lastAlgoX = algoSnmDist[algoPos].first.first;
            lastAlgoY = algoSnmDist[algoPos].first.second;
            lastDist.push_back(lastAlgoY->getDist());
            lastAlgoY->setDist(algoSnmDist[algoPos++].second);

            stringstream out;
            if(algoName == DIJKSTRA)out << "Dijkstra's algorithm.";
            else if(algoName == BELLMAN_FORD)out << "Bellman-Ford algorithm.";
            else if(algoName == FLOYD_WARSHALL)out << "Floyd-Warshall algorithm.";

            algoDescription = out.str();
        } else {
            if(algoPos == algoSnmDist.size()) {
                stringstream out;
                if(algoName == DIJKSTRA) out << "Dijkstra's algorithm ended.";
                else if(algoName == BELLMAN_FORD)out << "Bellman-Ford algorithm ended.";
                else if(algoName == FLOYD_WARSHALL)out << "Floyd-Warshall algorithm ended.";
                algoDescription = out.str();
                algoPos++;
            } else {
                setAlgorithmShow(false);
            }
        }
    }
    glutPostRedisplay();
}

void prevAlgorithmStep() {
    if(algoName == DFS || algoName == BFS) {
        if(algoPos > 1) {
            if(algoPos <= algo.size()) {
                algo[--algoPos].second->setVisible(false);
                lastAlgoX = algo[algoPos - 1].first;
                lastAlgoY = algo[algoPos - 1].second;
            } else algoPos--;
            stringstream out;
            if(algoName == DFS)
                out << "DFS algorithm. Component consists of " << algoPos << " vertexes";
            else out << "BFS algorithm. Component consists of " << algoPos << " vertexes";
            algoDescription = out.str();
        }
    } else if(algoName == KRUSKAL || algoName == PRIMA) {
        if(algoPos > 1) {
            stringstream out;
            if(algoName == KRUSKAL) {
                out << "Kruskal's algorithm. Cost of the minimum-spanning-tree is " <<  algoSnmDist[--algoPos - 1].second;
            } else {
                out << "Prim's algorithm. Cost of the minimum-spanning-tree is " <<  algoSnmDist[--algoPos - 1].second;
            }
            algoDescription = out.str();

            shownEd[{algoSnmDist[algoPos].first.first, algoSnmDist[algoPos].first.second}] = false;
            lastAlgoX = algoSnmDist[algoPos - 1].first.first;
            lastAlgoY = algoSnmDist[algoPos - 1].first.second;
        }
    } else if(algoName == DIJKSTRA || algoName == BELLMAN_FORD || algoName == FLOYD_WARSHALL) {
        if(algoPos > 1) {
            stringstream out;
            if(algoName == DIJKSTRA) {
                out << "Dijkstra algorithm.";
            } else if(algoName == BELLMAN_FORD) {
                out << "Bellman-Ford algorithm.";
            } else if(algoName == FLOYD_WARSHALL) {
                out << "Floyd-Warshall algorithm.";
            }
            algoDescription = out.str();
            --algoPos;
            if(algoPos < algoSnmDist.size()) {
                algoSnmDist[algoPos].first.second->setDist(lastDist.back());
                lastDist.pop_back();
            }
            lastAlgoX = algoSnmDist[algoPos - 1].first.first;
            lastAlgoY = algoSnmDist[algoPos - 1].first.second;
        }
    }
    glutPostRedisplay();
}

void runAlgorithm(string algoName) {
    ::algoName = algoName;
    algoPos = 0;
    if(algoName == DFS) {
        algo = graph.dfs(activeVertex);
        setAlgorithmShow(true);
        nextAlgorithmStep();
    }
    else if(algoName == BFS) {
        algo = graph.bfs(activeVertex);
        setAlgorithmShow(true);
        nextAlgorithmStep();
    }
    else if(algoName == KRUSKAL || algoName == PRIMA) {
        if(algoName == KRUSKAL) algoSnmDist = graph.kruskal(activeVertex);
        else                    algoSnmDist = graph.prima(activeVertex);
        setAlgorithmShow(true);
        for(size_t i = 0; i < algoSnmDist.size(); ++i) {
            algoSnmDist[i].first.first->setVisible(true);
            algoSnmDist[i].first.second->setVisible(true);
        }
        if(algoName == KRUSKAL) algoDescription = "Kruskal's algorithm. Cost of the minimum-spanning-tree is 0";
        else                    algoDescription = "Prima's algorithm. Cost of the minimum-spanning-tree is 0";
        nextAlgorithmStep();
    }else if(algoName == DIJKSTRA || algoName == BELLMAN_FORD || algoName == FLOYD_WARSHALL) {
        setAlgorithmShow(true);
        if(algoName == DIJKSTRA) algoSnmDist = graph.dijkstra(activeVertex);
        else if(algoName == BELLMAN_FORD) algoSnmDist = graph.bellman_ford(activeVertex);
        else if(algoName == FLOYD_WARSHALL) algoSnmDist = graph.floyd_warshall(activeVertex);

        for(size_t i = 0; i < algoSnmDist.size(); ++i) {
            algoSnmDist[i].first.first->setVisible(true);
            algoSnmDist[i].first.second->setVisible(true);
            algoSnmDist[i].first.second->setDist(MAX_COST * algoSnmDist.size());
        }
        distanse_shown = true;
        nextAlgorithmStep();
    }
}



void MouseMove(int mouseX, int mouseY) {
    double x = mouseX - WINDOW_WIDTH / 2;
    double y = -mouseY + WINDOW_HEIGHT / 2;

    active_polygon = nullptr;
    for(pair<Polygon, string>& but : Information) {
        if(intersect(but.first, {x, y}) && but.first.isActive) {
            active_polygon = &but.first;
            break;
        }
    }

    lastX = x;
    lastY = y;
    glutPostRedisplay();
}

void MousePressed(int button, int state, int mouseX, int mouseY) {
    double x = mouseX - WINDOW_WIDTH / 2;
    double y = -mouseY + WINDOW_HEIGHT / 2;
    if(is_algorithm_shown) return;
    if(state == GLUT_DOWN) mouse_button_clicked = button;
    switch(state) {
    case GLUT_DOWN :
        if(button == GLUT_LEFT_BUTTON && activeVertex)
        for(const pair<Polygon, string>& but : Information) {
            if(intersect(but.first, {x, y}) && but.first.isActive) {
                lastEdgeX = lastEdgeY = nullptr;
                runAlgorithm(but.second);
                return;
            }
        }
        deactivateAllVertexes();
        for(size_t i = graph.vertexes.size(); i > 0; i--) {
            Vertex& V = *graph.vertexes[i - 1];
            if(V.intersect(x, y)) {
                if(&V != lastEdgeX && &V != lastEdgeY) lastEdgeX = lastEdgeY = nullptr;
                V.activate();
                setActive(V);
                break;
            }
        }
        break;
    case GLUT_UP :
        intersectVertex = nullptr;
        if(button == GLUT_RIGHT_BUTTON && activeVertex)
        for(size_t i = graph.vertexes.size(); i > 0; i--) {
            Vertex& V = *graph.vertexes[i - 1];
            if(V.intersect(x, y) && V != *activeVertex) {
                lastEdgeX = lastEdgeY = nullptr;
                if(graph.areConnected(activeVertex, &V))
                    graph.deleteEdge(activeVertex, &V);
                else {
                    graph.addEdge(activeVertex, &V);
                    lastEdgeX = activeVertex;
                    lastEdgeY = &V;
                }
                setActive(*activeVertex);
                break;
            }
        }
        break;
    }
    lastX = x;
    lastY = y;
    glutPostRedisplay();
}

void MousePressedMove(int mouseX, int mouseY) {
    double x = mouseX - WINDOW_WIDTH / 2;
    double y = -mouseY + WINDOW_HEIGHT / 2;
    if(is_algorithm_shown) return;
    if(activeVertex) {
        if(mouse_button_clicked == GLUT_LEFT_BUTTON) {
            activeVertex->moveByVector(x - lastX, y - lastY);
        }
        if(mouse_button_clicked == GLUT_RIGHT_BUTTON) {
            intersectVertex = nullptr;
            for(size_t i = graph.vertexes.size(); i > 0; i--) {
                Vertex& V = *graph.vertexes[i - 1];
                if(V.intersect(x, y) && V != *activeVertex) {
                    intersectVertex = &V;
                    break;
                }
            }
        }
    }

    lastX = x;
    lastY = y;
    glutPostRedisplay();
}

void KeyPressed(unsigned char key, int mouseX, int mouseY) {
    double x = mouseX - WINDOW_WIDTH / 2;
    double y = -mouseY + WINDOW_HEIGHT / 2;
    if(lastEdgeX != nullptr && (!(key >= '0' && key <= '9') && !(key == '\b') && !(key == '-')))lastEdgeX = lastEdgeY = nullptr;
    if(lastEdgeX != nullptr) {
        long long cost = graph.getCost(lastEdgeX, lastEdgeY);
        if(key >= '0' && key <= '9')
            graph.setCost(lastEdgeX, lastEdgeY, cost * 10 + ((cost < 0) ? -(key - '0') : (key - '0')));
        else if(key == '\b')
            graph.setCost(lastEdgeX, lastEdgeY, cost / 10);
        else if(key == '-')
            graph.setCost(lastEdgeX, lastEdgeY, -cost);
        else lastEdgeX = lastEdgeY = nullptr;
    } else if(key == 27) {
        setAlgorithmShow(false);
        deactivateAllVertexes();
    } else if(key == 'i') {
        setAlgorithmShow(true);
    } else if(!is_algorithm_shown) {
        switch(key)
        {
        case 32 :
            break;
        case 'n' :
            graph.addVertex(new Vertex(x, y));
            break;
        case 242 :
            graph.addVertex(new Vertex(x, y));
            break;
        case '-' :
            graph.setNumbered(false);
            break;
        case '0' :
            graph.setNumeration(0);
            graph.setNumbered(true);
            break;
        case '1' :
            graph.setNumeration(1);
            graph.setNumbered(true);
            break;
        case 127:
            if(activeVertex) {
                graph.deleteVertex(activeVertex);
                deactivateAllVertexes();
            }
            break;
        }
    }
    glutPostRedisplay();
}
void SpecialKeyPressed(int key, int mouseX, int mouseY) {
    double x = mouseX - WINDOW_WIDTH / 2;
    double y = -mouseY + WINDOW_HEIGHT / 2;
    if(is_algorithm_shown) {
        if(key == GLUT_KEY_RIGHT) {
            nextAlgorithmStep();
        } else if(key == GLUT_KEY_LEFT) {
            prevAlgorithmStep();
        }

    } else {
        switch(key) {
        }
    }
    glutPostRedisplay();
}

void Initialize(){
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-WINDOW_WIDTH / 2, WINDOW_WIDTH / 2, -WINDOW_HEIGHT / 2, WINDOW_HEIGHT / 2, -200, 200);

    glutReshapeFunc([](int,int){
    glutReshapeWindow(WINDOW_WIDTH, WINDOW_HEIGHT);
    });

    glutDisplayFunc(Display);
//    glutTimerFunc(TIMER_DELAY, Timer, 0);

    glutPassiveMotionFunc(MouseMove);
    glutMouseFunc(MousePressed);
    glutMotionFunc(MousePressedMove);

    glutKeyboardFunc(KeyPressed);
    glutSpecialFunc(SpecialKeyPressed);
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("Grhap");

    graph.generateComponent(10, WINDOW_WIDTH / 2 - INFORMATION_WIDTH, WINDOW_HEIGHT / 2, 0.2, 1000);

    Initialize();
    glutMainLoop();

    return 0;
}
