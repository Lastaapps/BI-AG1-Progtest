#include <system_error>
#ifndef __PROGTEST__
#include <cassert>
#include <iostream>
#include <memory>
#include <limits>
#include <optional>
#include <algorithm>
#include <bitset>
#include <list>
#include <array>
#include <vector>
#include <deque>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <queue>

using Place = size_t;

struct Map {
  size_t places;
  Place start, end;
  std::vector<std::pair<Place, Place>> connections;
  std::vector<std::vector<Place>> items;
};

template < typename F, typename S >
struct std::hash<std::pair<F, S>> {
  std::size_t operator () (const std::pair<F, S> &p) const noexcept {
    // something like boost::combine would be much better
    return std::hash<F>()(p.first) ^ (std::hash<S>()(p.second) << 1);
  }
};

#endif
using namespace std;
using Point = Place;
using BackMap = unordered_map<Point, unordered_map<Point, vector<Point>>>;

struct Graph {
  unordered_multimap<Point, Point> edges;
};

inline unordered_set<Point> valuableRooms(const Map& map);
inline Graph constructGraph(const Map& map);
inline unordered_map<Point, vector<Point>> bfs(const Graph& graph, const Place start, const unordered_set<Place>& targets);
inline vector<Point> reconstructPath(const unordered_map<Point, Point> parents, const Point end);
inline vector<Point> findShortest(const BackMap& map, const vector<vector<Point>>& items, const Point start, const Point end);
inline list<Point> resolveBacktracking(const vector<Point>& path, const BackMap& map);

template<typename K, typename V>
ostream& operator <<(ostream &os, const pair<K, V>& data) {
  os << "<" << data.first << ", " << data.second << ">";
  return os;
}
ostream& operator <<(ostream &os, const vector<uint8_t>& data) {
  os << "{";
  for (const auto& i0 : data)
    os << (uint32_t)i0 << ", ";
  os << "}";
  return os;
}
template<typename T>
ostream& operator <<(ostream &os, const vector<T>& data) {
  os << "{";
  for (const auto& i0 : data)
    os << i0 << ", ";
  os << "}";
  return os;
}
template<typename T>
ostream& operator <<(ostream &os, const list<T>& data) {
  os << "[";
  for (const auto& i0 : data)
    os << i0 << ", ";
  os << "]";
  return os;
}
template<typename T>
ostream& operator <<(ostream &os, const queue<T>& data) {
  os << "[";
  for (const auto& i0 : data)
    os << i0 << ", ";
  os << "]";
  return os;
}
template<typename K, typename V>
ostream& operator <<(ostream &os, const unordered_map<K, V>& data) {
  os << "{";
  for (const auto& i0 : data)
    os << "(" << i0.first << ", " << i0.second << "), ";
  os << "}";
  return os;
}

list<Point> find_path(const Map &map) {
  const unordered_set<Point> valueable = valuableRooms(map);
  const Graph graph = constructGraph(map);

  //for (const auto& edge : graph.edges) cout << edge.first << " -> " << edge.second << endl;

  BackMap filtered;
  for (const auto p : valueable) {
    const auto res = bfs(graph, p, valueable);
    filtered.emplace(make_pair(p, res));
  }

  //for (const auto& i0 : filtered) {
  //  cout << i0.first << " -> " << endl;
  //  for (const auto& i1 : i0.second) {
  //    cout << " - " << i1.first << " = " << i1.second << endl;
  //  }
  //}
  
  const vector<Point> shortestPaht = findShortest(filtered, map.items, map.start, map.end);
  return resolveBacktracking(shortestPaht, filtered);
}

inline unordered_set<Point> valuableRooms(const Map& map) {
  unordered_set<Point> result;
  result.emplace(map.start);
  result.emplace(map.end);
  for (const auto& data : map.items)
    for (const auto item : data)
      result.emplace(item);
  return result;
}

inline Graph constructGraph(const Map& map) {
  Graph graph;
  for (const auto& edge : map.connections) {
    graph.edges.insert(make_pair(edge.first, edge.second));
    graph.edges.insert(make_pair(edge.second, edge.first));
  }
  return graph;
}

inline unordered_map<Point, vector<Point>> bfs(const Graph& graph, const Place start, const unordered_set<Place>& targets) {
  unordered_map<Point, Point> parents;
  unordered_set<Point> visited;
  unordered_map<Point, vector<Point>> result;
  size_t targesFound = 0;

  queue<Point> q;
  q.emplace(start);

  while (!q.empty()) {
    Point p = q.front();
    q.pop();
    visited.emplace(p);

    if (targets.find(p) != targets.end()) {
      result.emplace(make_pair(p, reconstructPath(parents, p)));
      ++targesFound;
    }

    if (targesFound == targets.size())
      break;

    const auto range = graph.edges.equal_range(p);
    for (auto itr = range.first; itr != range.second; ++itr) {
      Point to = itr -> second;
      if (visited.find(to) != visited.end())
        continue;

      q.emplace(to);
      parents.emplace(to, p);
    }
  }
  return result;
}

inline vector<Point> reconstructPath(const unordered_map<Point, Point> parents, Point end) {
  vector<Point> out;
  while(true) {
    out.emplace_back(end);

    auto next = parents.find(end);
    if (next == parents.end())
      break;

    end = next -> second;
  }
  return out;
}




using LevelPoint = pair<Point, uint16_t>;
inline LevelPoint levelPoint(const Point p, const uint16_t level) {
  return make_pair(p, level);
}

struct ValuedGraph {
  unordered_multimap<LevelPoint, pair<LevelPoint, size_t>> edges;
};
inline ValuedGraph constructValuedGraph(const BackMap& filtered, const vector<vector<Point>>& items, const Point start, const Point end, const vector<uint16_t> order);
inline void buildConnections(ValuedGraph& graph, const BackMap& map, const vector<Point>& src, const vector<Point>& dest, const uint16_t index);
inline vector<Point> dijkstra(const ValuedGraph& graph, const LevelPoint start, const LevelPoint end, const size_t shortest);
inline vector<Point> reconstructPath(const unordered_map<LevelPoint, LevelPoint> parents, LevelPoint end);

inline vector<Point> findShortest(const BackMap& map, const vector<vector<Point>>& items, const Point s, const Point e) {
  vector<uint16_t> indexes;
  for (uint16_t i = 0; i < items.size(); ++i)
    indexes.emplace_back(i);
  LevelPoint start = levelPoint(s, 0);
  LevelPoint end = levelPoint(e, items.size() + 1);

  vector<Point> shortest;
  do {
    const size_t shortestDistance = shortest.size() == 0 ? SIZE_MAX : shortest.size();
    const ValuedGraph graph = constructValuedGraph(map, items, s, e, indexes);

    //for (const auto& i0 : graph.edges) cout << i0.first << " -> (" << i0.second.first << ", " << i0.second.second << ")" << endl;

    const vector<Point> path = dijkstra(graph, start, end, shortestDistance);
    if (path.size() != 0) 
      shortest = move(path);
  } while(next_permutation(indexes.begin(), indexes.end()));
  return shortest;
}

inline ValuedGraph constructValuedGraph(const BackMap& map, const vector<vector<Point>>& items, const Point start, const Point end, const vector<uint16_t> order) {
  ValuedGraph graph;

  const vector<Point> starting = {start};
  vector<Point> const * prev = &starting;

  //cout << "Order: " << order << endl;
  for (size_t i = 0; i < order.size(); ++i) {
    const size_t index = order[i];
    const vector<Point>& current = items[index];
    //cout << "Current: " << current << endl;
    //cout << "Prev:    " << *prev << endl;
    buildConnections(graph, map, *prev, current, i);
    prev = &current;
  }
  //cout << "Current: " << "{" << end << ", }" << endl;
  //cout << "Prev:    " << *prev << endl;
  buildConnections(graph, map, *prev, {end}, order.size());

  return graph;
}

inline void buildConnections(ValuedGraph& graph, const BackMap& map, const vector<Point>& src, const vector<Point>& dest, const uint16_t index) {
  for (const auto from : src) {
    const auto& destinations = map.find(from);
    if (destinations == map.end())
      continue;
    const auto& points = destinations->second;

    for (const auto to : dest) {
      const auto& edge = points.find(to);
      if (edge == points.end())
        continue;
      const size_t length = edge->second.size() - 1;
      graph.edges.emplace(make_pair(levelPoint(from, index), make_pair(levelPoint(to, index + 1), length)));
    }
  }
}

inline vector<Point> dijkstra(const ValuedGraph& graph, const LevelPoint start, const LevelPoint end, const size_t shortest) {
  unordered_map<LevelPoint, size_t> evaluation;
  unordered_map<LevelPoint, LevelPoint> parents;
  queue<pair<Point, uint16_t>> q;
  bool endReached = false;

  q.emplace(start);
  evaluation.emplace(make_pair(start, 0));

  while (!q.empty()) {
    const auto p = q.front();
    q.pop();
    if (p == end) endReached = true;
    size_t length = evaluation.find(p) -> second;

    const auto range = graph.edges.equal_range(p);
    for (auto itr = range.first; itr != range.second; ++itr) {
      const LevelPoint current = itr->second.first;
      const size_t edgeLength = itr->second.second;
      const auto currItr = evaluation.find(current);
      const size_t sum = length + edgeLength;

      if (sum >= shortest)
        continue;
      
      if (currItr == evaluation.end()) {
        q.emplace(current);
        evaluation.emplace(make_pair(current, sum));
        parents.emplace(make_pair(current, p));
      } else {
        const size_t pathLength = currItr->second;
        if (length + edgeLength < pathLength)  {
          evaluation.find(current) -> second = sum;
          parents.find(current) -> second = p;
        }
      }
    }
  }
  //cout << "Evalul:  " << evaluation << endl;
  //cout << "Parents: " << parents << endl;

  return endReached ? reconstructPath(parents, end) : vector<Point>{};
}

inline vector<Point> reconstructPath(const unordered_map<LevelPoint, LevelPoint> parents, LevelPoint end) {
  vector<Point> out;
  if (parents.empty())
    return out;

  Point latest = -1;
  while(true) {
    if (end.first != latest)
      out.emplace_back(end.first);
    latest = end.first;

    auto next = parents.find(end);
    if (next == parents.end())
      break;

    end = next -> second;
  }
  //cout << "Out " << out << endl;
  return out;
}


inline list<Point> resolveBacktracking(const vector<Point>& path, const BackMap& map) {
  list<Point> result;
  if (path.empty()) return result;
  result.emplace_back(path[path.size() - 1]);

  for (auto pathItr = ++path.rbegin(); pathItr != path.rend(); ++ pathItr) {
    const vector<Point>& rooms = map.find(*pathItr)->second.find(*(pathItr - 1))->second;
    for (auto roomItr = ++rooms.begin(); roomItr != rooms.end(); ++roomItr)
      result.emplace_back(*roomItr);
  }

  return result;
}

#ifndef __PROGTEST__

using TestCase = std::pair<size_t, Map>;

const vector<TestCase> examples = {
  // 0
  TestCase{ 1, Map{ 2, 0, 0,
    { { 0, 1 } },
    { { 0 } }
  }},
  // 1
  TestCase{ 3, Map{ 2, 0, 0,
    { { 0, 1 } },
    { { 1 } }
  }},
  // 2
  TestCase{ 3, Map{ 4, 0, 1,
    { { 0, 2 }, { 2, 3 }, { 0, 3 }, { 3, 1 } },
    {}
  }},
  // 3
  TestCase{ 4, Map{ 4, 0, 1,
    { { 0, 2 }, { 2, 3 }, { 0, 3 }, { 3, 1 } },
    { { 2 } }
  }},
  // 4
  TestCase{ 0, Map{ 4, 0, 1,
    { { 0, 2 }, { 2, 3 }, { 0, 3 }, { 3, 1 } },
    { { 2 }, {} }
  }},
  // 5
  TestCase{ 0, Map{ 4, 0, 1,
    { { 0, 2 }, { 2, 3 }, { 0, 3 }, { 3, 1 } },
    { { 2 }, {} }
  }},
  // 6
  TestCase{ 0, Map{ 2, 0, 1,
    {},
    {}
  }},
  // 7
  TestCase{ 0, Map{ 3, 0, 1,
    {{0, 1}},
    {{2}}
  }},
  // 8
  TestCase{ 0, Map{ 3, 0, 1,
    {{0, 2}},
    {{2}}
  }},
  // 9
  TestCase{ 3, Map{ 3, 0, 2,
    {{0, 1}, {1, 2}},
    {}
  }},
  // 10
  TestCase{ 4, Map{ 2, 0, 3,
    {{0, 2}, {0, 1}, {1, 2}, {1, 3}, {2, 3}},
    {{2}, {1}}
  }},
  // 11
  TestCase{ 3, Map{ 2, 0, 3,
    {{0, 2}, {0, 1}, {1, 2}, {1, 3}, {2, 3}},
    {{2, 1}}
  }},
  // 12
  TestCase{ 1, Map{ 1, 0, 0,
    {},
    {}
  }},
  // 13
  TestCase{ 0, Map{ 2, 0, 0,
    {},
    {{1}}
  }},
  // 14
  TestCase{ 4, Map{ 7, 0, 1,
    {{0, 1}, {0, 3}, {0, 3}, {1, 2}, {1, 6}, {2, 3}, {4, 5}, {5, 6}},
    {{0, 1}, {3, 5}, {2, 4}}
  }},
  // 0
  //TestCase{ 0, Map{ 2, 0, 1,
  //  {},
  //  {}
  //}},
};

void memoryTest() {
  size_t max = 1000;
  size_t itemsMax = 4;
  vector<Point> rooms;
  vector<pair<Point, Point>> edges;
  vector<vector<Point>> items;

  for (size_t i = 0; i < max; ++i)
    rooms.emplace_back(i);
  for (size_t i = 0; i < 10 * max; ++i)
    edges.emplace_back(make_pair(rand() % 1000, rand() % 1000));
  for (size_t i = 0; i < itemsMax; ++i) {
    items.emplace_back(vector<Point>());
    for (size_t j = 0; j < max / 10; ++j)
      items[i].emplace_back(rand() % 1000);
  }
}

int main() {
  int fail = 0;
  for (size_t i = 0; i < examples.size(); i++) {
    std::cout << "Running test for map " << i << std::endl;

    auto sol = find_path(examples[i].second);
    cout << "Solution: " << sol << endl;

    if (sol.size() != examples[i].first) {
      std::cout << "Wrong answer for map " << i << std::endl;
      fail++;
    }
    cout << endl;
  }

  if (fail) std::cout << "Failed " << fail << " tests" << std::endl;
  else std::cout << "All tests completed" << std::endl;

  for (size_t i = 0; i < 1000; i++)
    memoryTest();

  return 0;
}

#endif
