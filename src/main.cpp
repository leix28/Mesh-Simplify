#include <cstdio>
#include <string>
#include <vector>
#include <set>
#include <cassert>
#include <cmath>
#include <queue>
#include <algorithm>
#define BUFFER_SIZE 1024
#define INFD 1e8
#define EPS 1e-8
#define TOLERATE 2.0
using std::min;
using std::max;
using std::make_pair;

typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;
typedef std::pair<int, int> Edge;

void printVector(const Vector &v) {
  for (auto x : v)
    printf("%.4lf\t", x);
  printf("\n");
}

void printMatrx(const Matrix &m) {
  for (auto v : m)
    printVector(v);
}

double norm(const Vector &v) {
  double t = 0;
  for (auto x : v) t += x * x;
  return sqrt(t);
}

Vector crossProduct(const Vector &a, const Vector &b) {
  assert(a.size() == 3 && b.size() == 3);
  Vector c(3);
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}

double innerProduct(const Vector &a, const Vector &b) {
  assert(a.size() == b.size());
  double c = 0;
  for (int i = 0; i < a.size(); i++)
    c += a[i] * b[i];
  return c;
}

Matrix outerProduct(const Vector &a, const Vector &b) {
  Matrix c(a.size(), Vector(b.size(), 0));
  for (int i = 0; i < a.size(); i++)
    for (int j = 0; j < b.size(); j++)
      c[i][j] = a[i] * b[j];
  return c;
}

Vector innerProduct(const Vector &a, const Matrix &b) {
  assert(a.size() == b.size());
  if (a.size() == 0) return Vector();
  Vector c(b[0].size(), 0);
  for (int i = 0; i < b.size(); i++)
    for (int j = 0; j < b[0].size(); j++)
      c[j] += a[i] * b[i][j];
  return c;
}

Vector operator + (const Vector &a, const Vector &b) {
  assert(a.size() == b.size());
  Vector c(a.size());
  for (int i = 0; i < a.size(); i++)
    c[i] = a[i] + b[i];
  return c;
}

Matrix operator + (const Matrix &a, const Matrix &b) {
  assert(a.size() == b.size());
  Matrix c(a.size());
  for (int i = 0; i < a.size(); i++)
    c[i] = a[i] + b[i];
  return c;
}

Vector operator - (const Vector &a, const Vector &b) {
  assert(a.size() == b.size());
  Vector c(a.size());
  for (int i = 0; i < a.size(); i++)
    c[i] = a[i] - b[i];
  return c;
}

Vector operator * (const double &a, const Vector &b) {
  Vector c(b.size());
  for (int i = 0; i < b.size(); i++)
    c[i] = a * b[i];
  return c;
}

Vector operator / (const Vector &a, const double &b) {
  assert(b != 0);
  Vector c(a.size());
  for (int i = 0; i < a.size(); i++)
    c[i] = a[i] / b;
  return c;
}

Vector solveEquation(Matrix m, int n) {
  assert(m.size() >= n);
  if (m.size() == 0) return Vector();
  assert(m[0].size() > n);
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++)
      if (fabs(m[i][i]) < fabs(m[j][i])) m[i].swap(m[j]);
    if (fabs(m[i][i]) < EPS) throw 200;
    m[i] = m[i] / m[i][i];
    for (int j = i + 1; j < n; j++)
      m[j] = m[j] - m[j][i] * m[i];
  }
  Vector v(n);
  for (int i = n - 1; i >= 0; i--) {
    assert(fabs(m[i][i] - 1) < EPS);
    v[i] = -m[i][n];
    for (int j = i + 1; j < n; j++) {
      v[i] -= m[i][j] * v[j];
    }
  }

  return v;
}


class Model {
  std::vector< Vector > vertex;
  std::vector<bool> removed;
  std::vector< std::set<Edge> > face;
  std::set<Edge> edge;
  std::priority_queue< std::pair<double, Edge> > heap;
  size_t faceN;

  double edgeLen(const Edge &e) {
    return norm(vertex[e.first] - vertex[e.second]);
  }

public:
  void clear() {
    vertex.clear();
    removed.clear();
    face.clear();
    edge.clear();
    faceN = 0;
  }

  size_t getEdgeN() {
    return edge.size();
  }

  size_t getVertexN() {
    return vertex.size();
  }

  size_t getFaceN() {
    return faceN;
  }

  void selfCheck() {
    std::set<Edge> ss;

    size_t vertexN = getVertexN();
    for (int i = 0; i < vertexN; i++) {
      if (removed[i]) {
        assert(face[i].size() == 0);
      } else {
        for (const auto &x : face[i]) {
          assert(!removed[x.first] && !removed[x.second]);
          assert(face[x.first].find(make_pair(x.second, i)) != face[x.first].end());
          assert(face[x.second].find(make_pair(i, x.first)) != face[x.second].end());
          ss.insert(make_pair(min(x.first, x.second), max(x.first, x.second)));
        }
      }
    }

    for (const auto &x : ss) {
      assert(x.first < x.second);
      assert(edge.find(x) != edge.end());
    }
    //there may be some e in edge which is actually absent.
    //(for models which are not closed)
  }

  void loadFromFile(std::string filename) {
    clear();


    char buffer[BUFFER_SIZE];
    FILE *file = fopen(filename.c_str(), "r");
    std::vector<std::string> vertexIn;
    std::vector<std::string> faceIn;
    while (fgets(buffer, BUFFER_SIZE, file) != NULL) {
      int ptr = 0;
      while (buffer[ptr] != 0 && buffer[ptr] != 'v' && buffer[ptr] != 'f' && buffer[ptr] != '#') ptr++;
      if (buffer[ptr] == 'v') vertexIn.push_back(std::string(buffer));
      if (buffer[ptr] == 'f') faceIn.push_back(std::string(buffer));
    }
    fclose(file);


    size_t vertexN = vertexIn.size();
    vertex.resize(vertexN, Vector(3, 0));
    removed.resize(vertexN, false);
    face.resize(vertexN);
    faceN = faceIn.size();

    for (int i = 0; i < vertexN; i++) {
      sscanf(vertexIn[i].c_str(), "%*s%lf%lf%lf", &vertex[i][0], &vertex[i][1], &vertex[i][2]);
    }


    for (const auto &f : faceIn) {
      int v[3];
      sscanf(f.c_str(), "%*s%d%d%d", v, v + 1, v + 2);
      v[0] --; v[1] --; v[2] --;
      face[v[0]].insert(make_pair(v[1], v[2]));
      face[v[1]].insert(make_pair(v[2], v[0]));
      face[v[2]].insert(make_pair(v[0], v[1]));
      std::sort(v, v + 3);
      assert(0 <= v[0] && v[0] < v[1] && v[1] < v[2] && v[2] < vertexN);
      edge.insert(make_pair(v[0], v[1]));
      edge.insert(make_pair(v[1], v[2]));
      edge.insert(make_pair(v[0], v[2]));
    }
  }

  void saveToFile(std::string filename) {
    FILE *file = fopen(filename.c_str(), "w");
    size_t vertexN = vertex.size();
    std::vector<int> vertexID(vertexN, 0);
    int vertexReal = 0;

    for (int i = 0; i < vertexN; i++) {
      if (removed[i]) continue;
      vertexID[i] = ++vertexReal;
      fprintf(file, "v %.8lf %.8lf %.8lf\n", vertex[i][0], vertex[i][1], vertex[i][2]);
    }

    for (int i = 0; i < vertexN; i++) {
      if (removed[i]) continue;
      for (const auto &f : face[i]) {
        assert(!removed[f.first] && !removed[f.second]);
        assert(vertexID[f.first] && vertexID[f.second] && vertexID[i]);
        if (i < f.first && i < f.second) {
          fprintf(file, "f %d %d %d\n", vertexID[i], vertexID[f.first], vertexID[f.second]);
        }
      }
    }
  }

  std::pair<Vector, double> getPosition(const Edge &e) {
    // std::set<int> neighbor1;
    // for (const auto &f : face[e.first]) {
    //   neighbor1.insert(f.first);
    //   neighbor1.insert(f.second);
    // }
    // std::set<int> neighbor2;
    // for (const auto &f : face[e.second]) {
    //   neighbor2.insert(f.first);
    //   neighbor2.insert(f.second);
    // }
    // int cnt = 0;
    // for (auto x : neighbor1) {
    //   if (neighbor2.find(x) != neighbor2.end()) cnt++;
    // }
    // if (cnt > 2) {
    //   return make_pair((vertex[e.first] + vertex[e.second]) / 2, -1);
    // }
    // if (cnt < 2) {
    //   return make_pair((vertex[e.first] + vertex[e.second]) / 2, INFD);
    // }

    Matrix q(4, Vector(4, 0));
    for (const auto &f : face[e.first]) {
      auto n = crossProduct(vertex[f.first] - vertex[e.first], vertex[f.second] - vertex[e.first]);
      n = n / norm(n);
      n.push_back(-innerProduct(vertex[e.first], n));
      q = q + outerProduct(n, n);
    }
    for (const auto &f : face[e.second]) {
      auto n = crossProduct(vertex[f.first] - vertex[e.second], vertex[f.second] - vertex[e.second]);
      n = n / norm(n);
      n.push_back(-innerProduct(vertex[e.second], n));
      q = q + outerProduct(n, n);
    }

    Vector v;
    try {
      v = solveEquation(q, 3);
    } catch(...) {
      v = (vertex[e.first] + vertex[e.second]) / 2;
    }
    if (norm(v - vertex[e.first]) + norm(v - vertex[e.second]) > TOLERATE * norm(vertex[e.first] - vertex[e.second])) {
      v = (vertex[e.first] + vertex[e.second]) / 2;
    }
    v.push_back(1);
    double cost = innerProduct(innerProduct(v, q), v);
    assert(cost > -EPS);
    v.pop_back();
    return make_pair(v, cost);
  }

  std::pair<Edge, Vector> selectEdge(double threshold) {
    Edge idx = make_pair(-1, -1);
    Vector pos;
    std::pair<double, Edge> tmp;
    while (!heap.empty()) {
      tmp = heap.top();
      heap.pop();
      if (edge.find(tmp.second) == edge.end()) continue;
      if (removed[tmp.second.first] || removed[tmp.second.second]) continue;
      if (edgeLen(tmp.second) > threshold) continue;
      auto act = getPosition(tmp.second);
      if (fabs(act.second + tmp.first) > EPS) continue;
      idx = tmp.second;
      pos = act.first;
      break;
    }
    printf("%lf %d %d", -tmp.first, idx.first, idx.second);
    return std::make_pair(idx, pos);
  }

  bool faceReverse(const Edge &e, const Vector &v1, const Vector &v2) {
    const auto &x = vertex[e.first];
    const auto &y = vertex[e.second];
    return innerProduct(crossProduct(x - v1, y - v1), crossProduct(x - v2, y - v2)) < 0;
    return 0;
  }

  void addToHeap(const Edge &e, double threshold) {
    if (edgeLen(e) > threshold) return;
    auto pos = getPosition(e);
    heap.push(make_pair(-pos.second, e));
  }

  void updateNeighborEdge(int v, double threshold) {
    std::set<int> neighbor;
    for (const auto &f : face[v]) {
      neighbor.insert(f.first);
      neighbor.insert(f.second);
    }
    for (auto x : neighbor) {
      addToHeap(make_pair(min(x, v), max(x, v)), threshold);
    }
  }

  void removeEdge(const Edge &e, const Vector &v, double threshold) {
    std::vector<Edge> toRev;
    for (const auto &f : face[e.first]) {
      if (f.first == e.second || f.second == e.second) continue;
      auto reverse = faceReverse(f, vertex[e.first], v);
      if (!reverse) continue;
      toRev.push_back(f);
      assert(face[f.second].find(make_pair(e.first, f.first)) != face[f.second].end());
      face[f.second].erase(make_pair(e.first, f.first));
      face[f.second].insert(make_pair(f.first, e.first));

      assert(face[f.first].find(make_pair(f.second, e.first)) != face[f.first].end());
      face[f.first].erase(make_pair(f.second, e.first));
      face[f.first].insert(make_pair(e.first, f.second));
    }
    for (const auto &f : toRev) {
      face[e.first].erase(f);
      face[e.first].insert(make_pair(f.second, f.first));
    }


    for (const auto &f : face[e.second]) {
      assert(face[f.second].find(make_pair(e.second, f.first)) != face[f.second].end());
      face[f.second].erase(make_pair(e.second, f.first));
      auto reverse = faceReverse(f, vertex[e.second], v);
      if (f.first != e.first && f.second != e.first) {
        if (reverse) {
          face[f.second].insert(make_pair(f.first, e.first));
        } else {
          face[f.second].insert(make_pair(e.first, f.first));
        }
      }

      assert(face[f.first].find(make_pair(f.second, e.second)) != face[f.first].end());
      face[f.first].erase(make_pair(f.second, e.second));
      if (f.first != e.first && f.second != e.first) {
        if (reverse) {
          face[f.first].insert(make_pair(e.first, f.second));
        } else {
          face[f.first].insert(make_pair(f.second, e.first));
        }
      }

      if (f.first == e.first || f.second == e.first)
        faceN--;
      else {
        if (reverse) {
          face[e.first].insert(make_pair(f.second, f.first));
        } else {
          face[e.first].insert(f);
        }
      }

      auto tmp = make_pair(min(e.second, f.first), max(e.second, f.first));
      if (edge.find(tmp) != edge.end())
        edge.erase(tmp);
      tmp = make_pair(min(e.second, f.second), max(e.second, f.second));
      if (edge.find(tmp) != edge.end())
        edge.erase(tmp);
      if (f.first != e.first && f.second != e.first) {
        edge.insert(make_pair(min(e.first, f.first), max(e.first, f.first)));
        edge.insert(make_pair(min(e.first, f.second), max(e.first, f.second)));
      }
    }

    edge.erase(e);
    vertex[e.first] = v;
    vertex[e.second].clear();
    removed[e.second] = true;
    face[e.second].clear();

    std::set<int> neighbor;
    for (const auto &f: face[e.first]) {
      neighbor.insert(f.first);
      neighbor.insert(f.second);
    }
    for (auto nb : neighbor) {
      updateNeighborEdge(nb, threshold);
    }
  }

  void buildHeap(double threshold) {
    while (!heap.empty()) heap.pop();
    for (const auto &e : edge) {
      addToHeap(e, threshold);
    }
  }

  void simplify(size_t target, double threshold) {
    buildHeap(threshold);
    while (faceN > target) {
      printf("%c%zu ", 13, faceN);
      auto e = selectEdge(threshold);
      if (e.first != make_pair(-1, -1))
        removeEdge(e.first, e.second, threshold);
      else {
        printf("%cERROR: No enough edges under threshold.\n", 13);
        printf("Warning: Current result will be save.\n");
        return;
      }
      // selfCheck();
      fflush(stdout);
    }
  }

} model;

int main(int argc, char **argv) {
  if (argc < 4) {
    printf("Usage:\n ./main [Input Object] [Output Object] [Simplify Rate] [Threshold Value]");
    return 0;
  }
  std::string inputModelFileName(argv[1]);
  std::string outputModelFileName(argv[2]);
  double simplifyRate = atof(argv[3]);
  double threshold;
  if (argc == 5) {
    threshold = atof(argv[4]);
  } else {
    printf("Warning: use threshold = INF (default)\n");
    threshold = INFD;
  }

  printf("inputModelFileName: %s\n", inputModelFileName.c_str());
  printf("outputModelFileName: %s\n", outputModelFileName.c_str());
  printf("simplifyRate: %.4lf\n", simplifyRate);
  printf("threshold: %.4lf\n", threshold);
  printf("------------------------------------\n");


  time_t start = time(0);

  model.loadFromFile(inputModelFileName);

  size_t all = model.getFaceN();
  size_t simple = all * simplifyRate;

  printf("vertex: %zu\n", model.getVertexN());
  printf("edge: %zu\n", model.getEdgeN());
  printf("simple / all = %zu / %zu\n", simple, all);
  model.simplify(simple, threshold);

  model.saveToFile(outputModelFileName);
  // model.selfCheck();
  time_t end = time(0);
  printf("%cSave to [%s] successfully. Time %ld sec.\n", 13, outputModelFileName.c_str(), end - start);
  return 0;
}
