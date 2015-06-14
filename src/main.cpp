#include <cstdio>
#include <string>
#include <vector>
#include <set>
#include <cassert>
#include <cmath>
#include <algorithm>
#define BUFFER_SIZE 1024
#define INFD 1e8
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

Vector operator - (const Vector &a, const Vector &b) {
  assert(a.size() == b.size());
  Vector c(a.size());
  for (int i = 0; i < a.size(); i++)
    c[i] = a[i] - b[i];
  return c;
}

Vector operator / (const Vector &a, const double &b) {
  assert(b != 0);
  Vector c(a.size());
  for (int i = 0; i < a.size(); i++)
    c[i] = a[i] / b;
  return c;
}




class Model {
  std::vector< Vector > vertex;
  std::vector<bool> removed;
  std::vector< std::set<Edge> > face;
  std::set<Edge> edge;
  int faceN;

  double edgeLen(Edge e) {
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

  int getEdgeN() {
    return edge.size();
  }

  int getVertexN() {
    return vertex.size();
  }

  int getFaceN() {
    return faceN;
  }

  void selfCheck() {
    std::set<Edge> ss;

    int vertexN = getVertexN();
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

    for (const auto &x : ss)
      assert(edge.find(x) != edge.end());
    for (const auto &x : edge)
      assert(ss.find(x) != ss.end());
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


    int vertexN = vertexIn.size();
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
    int vertexN = vertex.size();
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

  double getCost(int vid, Vector vpos) {
    auto v = vertex[vid] - vpos;
    // printVector(v);
    double cost = 0;
    for (const auto &f : face[vid]) {
      auto n = crossProduct(vertex[f.first] - vertex[vid], vertex[f.second] - vertex[vid]);
      n = n / norm(n);
      // printVector(n);
      cost += innerProduct(v, n) * innerProduct(v, n);
    }
    // printf("%.4lf\n", cost);
    return cost;
  }

  Edge selectEdge(double threshold) {
    Edge idx = make_pair(-1, -1);
    double best = INFD;
    for (const auto &e : edge) {
      if (edgeLen(e) > threshold) continue;
      auto v = (vertex[e.first] + vertex[e.second]) / 2.0;
      double c = getCost(e.first, v) + getCost(e.second, v);
      // printf("%.4lf\n", c);
      if (c < best) {
        best = c;
        idx = e;
      }
    }
    assert(idx != make_pair(-1, -1));
    printf("%lf %d %d", best, idx.first, idx.second);
    return idx;
  }

  void removeEdge(Edge e) {
    auto v = (vertex[e.first] + vertex[e.second]) / 2.0;
    edge.erase(e);
    vertex[e.first] = v;
    vertex[e.second].clear();
    removed[e.second] = true;

    for (const auto &f : face[e.second]) {
      assert(face[f.second].find(make_pair(e.second, f.first)) != face[f.second].end());
      face[f.second].erase(make_pair(e.second, f.first));
      if (f.first != e.first && f.second != e.first) {
        face[f.second].insert(make_pair(e.first, f.first));
      }

      assert(face[f.first].find(make_pair(f.second, e.second)) != face[f.first].end());
      face[f.first].erase(make_pair(f.second, e.second));
      if (f.first != e.first && f.second != e.first) {
        face[f.first].insert(make_pair(f.second, e.first));
      }

      if (f.first == e.first || f.second == e.first)
        faceN--;
      else
        face[e.first].insert(f);

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
    face[e.second].clear();
  }

  void simplify(int target, double threshold) {
    while (faceN > target) {
      printf("%c%d\t", 13, faceN);
      auto e = selectEdge(threshold);
      removeEdge(e);
      selfCheck();
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

  model.loadFromFile(inputModelFileName);

  int all = model.getFaceN();
  int simple = all * simplifyRate;

  printf("vertex: %d\n", model.getVertexN());
  printf("edge: %d\n", model.getEdgeN());
  printf("simple / all = %d / %d\n", simple, all);
  model.simplify(simple, threshold);

  model.saveToFile(outputModelFileName);
  model.selfCheck();
  return 0;
}
