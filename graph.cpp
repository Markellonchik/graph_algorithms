#include "graph.h"

Vertex::Vertex(double x, double y, double radius) : x(x), y(y), radius(radius), color(0, 0, 0){
}
Vertex::Vertex(double x, double y, double radius, Color color) : Vertex(x, y, radius){
    this->color = color;
}
bool Vertex::operator== (const Vertex& other) const{
    return num == other.getNum();
}
bool Vertex::operator!= (const Vertex& other) const{
    return !(*this == other);
}
double Vertex::getX() const{
    return x;
}
double Vertex::getY() const{
    return y;
}
int Vertex::getNum() const{
    return num;
}
double Vertex::getRadius() const{
    return radius;
}
long long Vertex::getDist() const {
    return distanse;
}
bool Vertex::isActive() const {
    return is_active;
}
bool Vertex::isVisible() const {
    return visible;
}

bool Vertex::intersect(double xa, double ya) const {
    return sqrt( pow(x - xa, 2) + pow(y - ya, 2) ) < radius;
}

void Vertex::moveTo(double xA, double yA) {
    x = xA;
    y = yA;
}
void Vertex::moveByVector(double addX, double addY) {
    x += addX;
    y += addY;
}
void Vertex::activate() {
    is_active = true;
}
void Vertex::deactivate() {
    is_active = false;
}
void Vertex::setNum(int nm) {
    num = nm;
}
void Vertex::setVisible(bool vis) {
    visible = vis;
}
void Vertex::setDist(long long dist) {
    distanse = dist;
}

void Graph::renumerate() {
    for(size_t i = 0; i < vertexes.size(); ++i)
        vertexes[i]->setNum(i);
}
void Graph::resize(int size) {
    edges.assign(size, vector<int>());
    costs.assign(size, vector<long long>());
    vertexes.assign(size, nullptr);
    for(size_t i = 0; i < size; ++i)vertexes[i] = new Vertex(0, 0, STANDART_RADIUS);
    renumerate();
}
Graph::Graph(int size) {
    srand(time(0));
    resize(size);
}
Graph::Graph(int size, bool oriented) : Graph(size) {
    this->oriented = oriented;
}
bool Graph::isOriented() const {
    return oriented;
}

void Graph::setNumbered(bool yes) {
    numbered = yes;
}
void Graph::setNumeration(int numFirst) {
    numeration = numFirst;
}

void Graph::deleteVertex(int num) {
    for(size_t i = 0; i < edges.size(); ++i) {
        for(size_t j = 0; j < edges[i].size(); ++j) {
            if(edges[i][j] == num) {
                edges[i].erase(edges[i].begin() + j);
                costs[i].erase(costs[i].begin() + j);
                break;
            }
        }
    }
    for(size_t i = num; i + 1 < vertexes.size(); ++i) {
        vertexes[i] = vertexes[i + 1];
        costs[i] = costs[i + 1];
        edges[i] = edges[i + 1];
    }

    vertexes.pop_back();
    costs.pop_back();
    edges.pop_back();
    for(size_t i = 0; i < edges.size(); ++i) {
        for(int& j : edges[i])
            if(j >= num) j--;
    }

    renumerate();
}

void Graph::deleteVertex(const Vertex* X) {
    deleteVertex(X->getNum());
}
void Graph::addVertex(Vertex* X) {
    vertexes.push_back(X);
    edges.push_back({});
    costs.push_back({});
    renumerate();
}
int Graph::getNumeration() const {
    return numeration;
}
bool Graph::isNumerated() const {
    return numbered;
}
bool Graph::areConnected(int x, int y) const{
    return count(edges[x].begin(), edges[x].end(), y);
}
bool Graph::areConnected(const Vertex* X, const Vertex* Y) const{
    return areConnected(X->getNum(), Y->getNum());
}
int Graph::outDegree(const Vertex* X) const {
    return edges[X->getNum()].size();
}
int Graph::inDegree(const Vertex* X) const {
    int res = 0;
    for(size_t i = 0; i < edges.size(); ++i) {
        res += count(edges[i].begin(), edges[i].end(), X->getNum());
    }
    return res;
}

void Graph::setCost(int x, int y, long long cost) {
    cost = min(cost, MAX_COST);
    cost = max(cost, MIN_COST);
    for(size_t i = 0; i < edges[x].size(); ++i) {
        if(edges[x][i] == y) costs[x][i] = cost;
    }
    for(size_t i = 0; i < edges[y].size(); ++i) {
        if(edges[y][i] == x) costs[y][i] = cost;
    }
}

void Graph::setCost(const Vertex* X, const Vertex* Y, long long cost) {
    setCost(X->getNum(), Y->getNum(), cost);
}

int Graph::getCost(int x, int y) const {
    for(size_t i = 0; i < edges[x].size(); ++i) {
        if(edges[x][i] == y) return costs[x][i];
    }
    return -1;
}

int Graph::getCost(const Vertex* X, const Vertex* Y) const{
    return getCost(X->getNum(), Y->getNum());
}

void Graph::deleteEdge(int x, int y) {
    for(size_t i = 0; i < edges[x].size(); ++i)
        if(edges[x][i] == y) {
            edges[x].erase(edges[x].begin() + i);
            costs[x].erase(costs[x].begin() + i);
        }
    if(!oriented) {
        for(size_t i = 0; i < edges[y].size(); ++i)
            if(edges[y][i] == x) {
                edges[y].erase(edges[y].begin() + i);
                costs[y].erase(costs[y].begin() + i);
            }
    }
}
void Graph::deleteEdge(const Vertex* X, const Vertex* Y) {
    deleteEdge(X->getNum(), Y->getNum());
}
void Graph::addEdge(int XNum, int YNum, long long cost) {
    cost = min(cost, MAX_COST);
    cost = max(cost, MIN_COST);
    for(const int& nm : edges[XNum]) {
        if(nm == YNum) return;
    }
    edges[XNum].push_back(YNum);
    costs[XNum].push_back(cost);

    if(!oriented) {
        edges[YNum].push_back(XNum);
        costs[YNum].push_back(cost);
    }
}

void Graph::addEdge(const Vertex* X, const Vertex* Y, long long cost) {
    addEdge(X->getNum(), Y->getNum(), cost);
}
void Graph::deactivateAllVertex() {
    for(Vertex* V : vertexes) V->deactivate();
}
void Graph::setVisible(bool vis) {
    for(Vertex* V : vertexes) V->setVisible(vis);
}
void Graph::dfsAlgo(Vertex* X, vector<bool>& used, vector<pair<Vertex*, Vertex*> >& r){
    used[X->getNum()] = true;
    for(const int& j : edges[X->getNum()]) {
        if(!used[j]) {
            r.push_back({X, vertexes[j]});
            dfsAlgo(vertexes[j], used, r);
        }
    }
}
vector<pair<Vertex*, Vertex*> > Graph::dfs(Vertex* X){
    vector<pair<Vertex*, Vertex*> > r;
    vector<bool> used(vertexes.size(), false);
    r.push_back({nullptr, X});
    dfsAlgo(X, used, r);
    return r;
}
void Graph::bfsAlgo(Vertex* X, vector<bool>& used, vector<pair<Vertex*, Vertex*> >& r) {
    queue<size_t> q;
    vector<int> pr(vertexes.size(), 0);
    used[X->getNum()] = true;
    q.push(X->getNum());

    while(!q.empty()) {
        size_t& i = q.front();
        q.pop();
        if((int)i != X->getNum())
        r.push_back({vertexes[pr[i]], vertexes[i]});
        for(int& j : edges[i]) {
            if(!used[j]) {
                pr[j] = i;
                used[j] = true;
                q.push(j);
            }
        }
    }
}
vector<pair<Vertex*, Vertex*> > Graph::bfs(Vertex* X){
    vector<pair<Vertex*, Vertex*> > r;
    vector<bool> used(vertexes.size(), false);
    r.push_back({nullptr, X});
    bfsAlgo(X, used, r);
    return r;
}

int Graph::getSnmPr(int X, vector<int>& pr) {
    if(pr[X] == X) return X;
    return pr[X] = getSnmPr(pr[X], pr);
}
vector<pair<pair<Vertex*, Vertex*>, long long> > Graph::kruskal(Vertex* X) {
    int N = vertexes.size();
    vector<pair<Vertex*, Vertex*> > component;
    vector<bool> used(N, false);
    component.push_back({nullptr, X});
    dfsAlgo(X, used, component);

    vector<pair<long long, pair<int, int> > > componentEdges;

    for(size_t k = 0; k < component.size(); ++k) {
        int XNum = component[k].second->getNum();
        for(size_t i = 0; i < edges[XNum].size(); ++i) {
            int j = edges[XNum][i];
            if(used[j]) {
                componentEdges.push_back({costs[XNum][i], {XNum, j}});
            }
        }
    }

    vector<int> pr(N, 0), sz(N, 1);
    for(int i = 0; i < N; ++i)pr[i] = i;
    sort(begin(componentEdges), end(componentEdges));

    vector<pair<pair<Vertex*, Vertex*>, long long> > res;
    res.push_back({{X, X}, 0});
    long long sum = 0;
    for(pair<long long, pair<int, int> >& edge : componentEdges) {
        int X = getSnmPr(edge.second.first, pr);
        int Y = getSnmPr(edge.second.second, pr);
        if(X == Y) continue;

        if(sz[X] < sz[Y]) swap(X, Y), swap(edge.second.first, edge.second.second);
        sz[X] += sz[Y];
        pr[Y] = X;

        sum += edge.first;
        res.push_back({{vertexes[edge.second.first], vertexes[edge.second.second]}, sum});
    }
    return res;
}

vector<pair<pair<Vertex*, Vertex*>, long long> > Graph::prima(Vertex* X) {
    vector<pair<pair<Vertex*, Vertex*>, long long> > res;
    res.push_back({{X, X}, 0});
    set<pair<long long, pair<int,int> > > ed;
    vector<bool> used(vertexes.size(), false);

    int XNum = X->getNum(), YNum;
    used[XNum] = true;
    for(size_t i = 0; i < edges[XNum].size(); ++i) {
        ed.insert({costs[XNum][i], {XNum, edges[XNum][i]}});
    }

    long long sum = 0;
    while(ed.size() > 0) {
        XNum = ed.begin()->second.first;
        YNum = ed.begin()->second.second;
        used[YNum] = true;

        sum += ed.begin()->first;
        ed.erase(ed.begin());
        res.push_back({{vertexes[XNum], vertexes[YNum]}, sum});

        for(size_t i = 0; i < edges[YNum].size(); ++i) {
            int j = edges[YNum][i];
            if(used[j]) {
                ed.erase({costs[YNum][i], {j, YNum}});
            } else {
                ed.insert({costs[YNum][i], {YNum, j}});
            }
        }

    }

    return res;
}

vector<pair<pair<Vertex*, Vertex*>, long long> > Graph::dijkstra(Vertex* X) {
    vector<pair<pair<Vertex*, Vertex*>, long long> > res;
    res.push_back({{X, X}, 0});
    int N = vertexes.size();

    vector<long long> dist(N, N * MAX_COST);
    set<pair<long long, int > > st;

    int XNum = X->getNum(), YNum;
    dist[XNum] = 0;

    st.insert({dist[XNum], XNum});
    while(st.size() > 0) {
        XNum = st.begin()->second;
        st.erase(st.begin());
        for(size_t i = 0; i < edges[XNum].size(); ++i) {
            int j = edges[XNum][i];
            int cost = costs[XNum][i];
            if(dist[j] > dist[XNum] + cost) {
                st.erase({dist[j], j});
                dist[j] = dist[XNum] + cost;
                res.push_back({{vertexes[XNum], vertexes[j]}, dist[j]});
                st.insert({dist[j], j});
            }
        }
    }

    return res;
}

vector<pair<pair<Vertex*, Vertex*>, long long> > Graph::bellman_ford(Vertex* X) {
    vector<pair<pair<Vertex*, Vertex*>, long long> > res;
    res.push_back({{X, X}, 0});

    int N = vertexes.size();
    vector<long long> dist(N, N * MAX_COST);

    int XNum = X->getNum(), YNum;
    dist[XNum] = 0;

    for(int iter = 0; iter < N - 1; ++iter) {
        for(size_t i = 0; i < vertexes.size(); ++i) {
            for(size_t j = 0; j < edges[i].size(); ++j) {
                int to = edges[i][j];
                int cost = costs[i][j];
                if(dist[to] > dist[i] + cost) {
                    dist[to] = dist[i] + cost;
                    res.push_back({{vertexes[i], vertexes[to]}, dist[to]});
                }
            }
        }
    }

    return res;
}

vector<pair<pair<Vertex*, Vertex*>, long long> > Graph::floyd_warshall(Vertex* X) {
    vector<pair<pair<Vertex*, Vertex*>, long long> > res;
    int N = vertexes.size();
    vector<vector<long long> > dist(N, vector<long long>(N, N * MAX_COST));
    vector<bool> used(N, false);

    vector<pair<Vertex*, Vertex*> > G;
    dfsAlgo(X, used, G);

    for(size_t i = 0; i < N; ++i) {
        for(size_t j = 0; j < edges[i].size(); ++j) {
            dist[i][edges[i][j]] = costs[i][j];
        }
    }
    for(size_t k = 0; k < N; ++k) {
        for(size_t i = 0; i < N; ++i) {
            for(size_t j = 0; j < N; ++j) {
                if(dist[i][j] > dist[i][k] + dist[k][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
    res.push_back({{X, X}, 0});
    for(size_t i = 0; i < N; ++i) {
        if(i != X->getNum() && used[i])
        res.push_back({{X, vertexes[i]}, dist[X->getNum()][i]});
    }
    return res;
}

vector<pair<pair<Vertex*, Vertex*>, long long> > Graph::jonson(Vertex* X) {

}

Vertex* Graph::randVertex(int width, int height) const{
    double x = rand() % (width + width + 1) - width;
    double y = rand() % (height + height + 1) - height;
    return new Vertex(x, y);
}

long long Graph::rnd() const {
    long long r = rand();
    return (r << 15) + rand();
}

void Graph::generateComponent(int size, int width, int height, double probality, int maxCost) {
    vector<Vertex*> component(size);
    for(size_t i = 0; i < component.size(); ++i) {
        component[i] = randVertex(width - STANDART_RADIUS * 2, height - STANDART_RADIUS * 2);
        addVertex(component[i]);
    }
    for(size_t i = 0; i < component.size(); ++i) {
        for(size_t j = 0; j < component.size(); ++j) {
            if(i == j)continue;
            if(i > j && !oriented)continue;
            if(rand() % 10001 < probality * 10000) {
                addEdge(component[i], component[j], rnd() % (maxCost + 1));
            }
        }
    }
}

