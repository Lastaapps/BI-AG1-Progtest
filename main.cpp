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
using Mask = bitset<12>;
using LookupMap = unordered_map<Point, Mask>;
using StatePoint = pair<Point, Mask>;
using ParentsMap = unordered_map<StatePoint, StatePoint>;

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

struct Graph {
  unordered_multimap<Point, Point> edges;
};



inline LookupMap buildItemsLookup(const Map& params) {
  const vector<vector<Point>>& src = params.items;
  unordered_map<Point, vector<Point>> map; 

  for (uint8_t i = 0; i < src.size(); ++i)
    for (const auto point : src[i]) {
      const auto found = map.find(point);
      if (found == map.end())
        map.emplace(point, vector<Point>{i});
      else
        found -> second.emplace_back(i);
    }
  LookupMap lookup;
  for (auto const& entry: map) {
    Mask m = 0;
    for (const auto item : entry.second)
      m |= 1 << item;
    lookup.emplace(entry.first, m);
  }
  if (lookup.count(params.start) == 0)
    lookup.emplace(params.start, 0);
  if (lookup.count(params.end) == 0)
    lookup.emplace(params.end, 0);
  return lookup;
}


inline Mask createEndMask(const uint8_t items) {
  uint16_t m = 0;
  for (uint8_t i = 0; i < items; ++i) {
      m |= 1 << i;
  }
  return m;
}

inline vector<Point> backtrack(const ParentsMap& parents, const StatePoint end) {
  StatePoint latest = make_pair(-1, 0);
  StatePoint from = end;
  vector<Point> out;

  while(true) {
    const auto next = parents.find(from) -> second;
    if (next == latest) break;
    out.emplace_back(from.first);
    latest = from;
    from = next;
  }
  return out;
}

inline vector<Point> smartBfs(const Graph& graph, const LookupMap& lookup, const Map& params) {
  queue<StatePoint> q;
  ParentsMap parents;

  const auto startMaskSearch = lookup.find(params.start);
  const Mask startMask = startMaskSearch == lookup.end() ? 0 : startMaskSearch->second;
  const Mask endMask = createEndMask(params.items.size());

  const StatePoint start = make_pair(params.start, startMask);
  const StatePoint end = make_pair(params.end, endMask);
  //cout << "Starting BFS from " << start << " to " << end << endl;

  q.emplace(start);
  parents.emplace(start, start);

  while (!q.empty()) {
    const StatePoint p = q.front();
    q.pop();

    if (p == end) return backtrack(parents, p);

    const auto neighbors = graph.edges.equal_range(p.first);
    for (auto itr = neighbors.first; itr != neighbors.second; ++itr) {
      const auto target = itr->second;
      const auto items = lookup.find(target);
      const Mask nextMask = p.second | (items == lookup.end() ? 0 : items->second);
      const StatePoint nextState = make_pair(target, nextMask);
      //cout << "Child: " << p << " -> " << target << endl;
      if (parents.count(nextState) == 0) {
        parents.emplace(nextState, p);
        q.emplace(nextState);
      }
    }
  }
  return vector<Point>();
}


inline Graph buildGraph(const Map& map) {
  Graph graph;
  for (const auto& edge : map.connections) {
    graph.edges.insert(make_pair(edge.first, edge.second));
    graph.edges.insert(make_pair(edge.second, edge.first));
  }
  return graph;
}

inline std::list<Place> find_short_path(const Map &map) {
  const Graph graph = buildGraph(map);
  const LookupMap lookup = buildItemsLookup(map);
  //for (const auto& edge : graph.edges) cout << edge.first << " -> " << edge.second << endl;
  
  const vector<Point> points = smartBfs(graph, lookup, map);
  list<Point> result;
  for (auto itr = points.rbegin(); itr != points.rend(); ++itr)
    result.emplace_back(*itr);
  return result;
}



// --- LONG PATH --------------------------------------------------------------

using ValuedGraph = unordered_multimap<Point, pair<Point, size_t>>;

inline unordered_map<Point, size_t> bfsDistance(const Graph& graph, const Place start, const unordered_set<Place>& targets) {
  unordered_set<Point> visited;
  unordered_map<Point, size_t> result;
  size_t targesFound = 0;

  queue<pair<Point, size_t>> q;
  q.emplace(make_pair(start, 0));

  while (!q.empty()) {
    const auto p = q.front();
    q.pop();

    if (targets.find(p.first) != targets.end()) {
      result.emplace(p);
      ++targesFound;
    }

    if (targesFound >= targets.size())
      break;

    const auto range = graph.edges.equal_range(p.first);
    for (auto itr = range.first; itr != range.second; ++itr) {
      const Point to = itr -> second;
      if (visited.find(to) != visited.end())
        continue;

      visited.emplace(to);
      q.emplace(make_pair(to, p.second + 1));
    }
  }
  return result;
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

ValuedGraph buildValuedGraph(const Graph& graph, const LookupMap& lookup, const Map& params) {
  auto valueable = valuableRooms(params);
  ValuedGraph distances;

  auto itr = valueable.begin();
  while (itr != valueable.end()) {
    const auto p = *itr;
    auto tmpItr = itr++;
    valueable.erase(tmpItr);

    auto current = valueable;

    const Mask groups = lookup.find(p)->second;

    for (uint8_t i = 0; i < params.items.size(); ++i) {
      if ((groups & Mask(1 << i)) == 0) continue; 
      for (const auto neighbor : params.items[i])
        current.erase(neighbor);
    }

    const auto res = bfsDistance(graph, p, current);
    for (auto const& entry : res) {
      distances.emplace(p, entry);
      distances.emplace(entry.first, make_pair(p, entry.second));
    }
  }
  return distances;
}

ValuedGraph buildValuedGraphOld(const Graph& graph, const LookupMap& lookup, const Map& params) {
  const auto valueable = valuableRooms(params);
  ValuedGraph distances;
  for (const auto p : valueable) {
    const auto res = bfsDistance(graph, p, valueable);
    for (auto const& entry : res)
      distances.emplace(p, entry);
  }
  return distances;
}

vector<Point> smartDijikstra(const ValuedGraph& graph, const LookupMap& lookup, const Map& params) {
  unordered_map<StatePoint, size_t> evaluation;
  unordered_map<StatePoint, StatePoint> parents;
  queue<StatePoint> q;
  bool endReached = false;

  const auto startMaskSearch = lookup.find(params.start);
  const Mask startMask = startMaskSearch == lookup.end() ? 0 : startMaskSearch->second;
  const Mask endMask = createEndMask(params.items.size());

  const StatePoint start = make_pair(params.start, startMask);
  const StatePoint end = make_pair(params.end, endMask);

  q.emplace(start);
  evaluation.emplace(start, 0);
  parents.emplace(start, start);

  while (!q.empty()) {
    const auto p = q.front();
    q.pop();
    if (p == end) endReached = true;
    const size_t base = evaluation.find(p) -> second;

    const auto range = graph.equal_range(p.first);
    for (auto itr = range.first; itr != range.second; ++itr) {
      const Point targetPoint = itr->second.first;
      const size_t edgeLength = itr->second.second;

      const auto items = lookup.find(targetPoint);
      const Mask nextMask = p.second | items->second;

      const StatePoint target = make_pair(targetPoint, nextMask);
      const auto currItr = evaluation.find(target);
      const size_t sum = base + edgeLength;

      if (currItr == evaluation.end()) {
        q.emplace(target);
        evaluation.emplace(target, sum);
        parents.emplace(target, p);
      } else {
        const size_t pathLength = currItr->second;
        if (base + edgeLength < pathLength)  {
          evaluation.find(target) -> second = sum;
          parents.find(target) -> second = p;
        }
      }
    }
  }
  //cout << "Evalul:  " << evaluation << endl;
  //cout << "Parents: " << parents << endl;
  //cout << "End:     " << end <<  ", reached: " << endReached << endl;

  return endReached ? backtrack(parents, end) : vector<Point>{};
}

inline vector<Point> reconstructPath(const unordered_map<Point, Point> parents, Point end) {
  vector<Point> out;
  Point prev = -1;
  //cout << "Parents: " << parents << endl;
  //cout << "End:     " << end << endl;
  while(true) {
    if (end == prev) break;

    out.emplace_back(end);

    auto next = parents.find(end);
    if (next == parents.end())
      break;

    prev = end;
    end = next -> second;
  }
  return out;
}

inline vector<Point> bfsPredecesors(const Graph& graph, const Place start, const Place end) {
  unordered_map<Point, Point> parents;

  queue<Point> q;
  q.emplace(start);
  parents.emplace(start, start);

  while (!q.empty()) {
    const Point p = q.front();
    q.pop();

    if (p == end)
      return reconstructPath(parents, p);

    const auto range = graph.edges.equal_range(p);
    for (auto itr = range.first; itr != range.second; ++itr) {
      const Point to = itr -> second;
      if (parents.find(to) != parents.end())
        continue;

      q.emplace(to);
      parents.emplace(to, p);
    }
  }
  return vector<Point>();
}

inline list<Point> resolveBacktracking(const Graph& graph, const vector<Point>& path) {
  list<Point> result;
  if (path.empty()) return result;
  result.emplace_back(path[path.size() - 1]);

  for (auto pathItr = ++path.rbegin(); pathItr != path.rend(); ++ pathItr) {
    Point from = *(pathItr - 1);
    Point to = *pathItr;
    const auto rooms = bfsPredecesors(graph, from, to);
    //cout << "Rooms: " << rooms << endl;

    for (auto roomItr = ++rooms.rbegin(); roomItr != rooms.rend(); ++roomItr)
      result.emplace_back(*roomItr);
  }
  return result;
}

inline std::list<Place> find_long_path(const Map &map) {
  //cout << "Dij: Init" << endl;
  const Graph graph = buildGraph(map);
  const LookupMap lookup = buildItemsLookup(map);
  //cout << "Dij: Evaulation" << endl;
  const ValuedGraph valued = buildValuedGraph(graph, lookup, map);
  //const ValuedGraph valued = buildValuedGraphOld(graph, lookup, map);
  //for (const auto& i0 : valued) cout << i0.first << " -> (" << i0.second.first << ", " << i0.second.second << ")" << endl;

  //cout << "Dij: Dijkstra" << endl;
  const vector<Point> points = smartDijikstra(valued, lookup, map);
  //cout << "Dij: Resolving" << endl;
  return resolveBacktracking(graph, points);
}

std::list<Place> find_path(const Map &map) {
  return (map.places < 4098 && map.items.size() <= 6) ? find_short_path(map) : find_long_path(map);
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
  cout << "Mem: Init" << endl;
  size_t rooms = 10000;
  size_t itemsMax = 12;
  size_t maxInRoom = 5;
  vector<pair<Point, Point>> edges;
  vector<vector<Point>> items;

  for (size_t i = 0; i < 10 * rooms; ++i)
    edges.emplace_back(make_pair(rand() % rooms, rand() % rooms));
  for (size_t i = 0; i < itemsMax; ++i) {
    items.emplace_back(vector<Point>());
    for (size_t j = 0; j < maxInRoom; ++j)
      items[i].emplace_back(rand() % rooms);
  }
  const auto map = Map{rooms, rand() % rooms, rand() % rooms, edges, items};
  cout << "Mem: Run" << endl;
  find_path(map);
  cout << "Mem: Done" << endl;
}

int main() {
  int fail = 0;
  for (size_t i = 0; i < examples.size(); i++) {
    std::cout << "OLD: Running test for map " << i << std::endl;
    auto solOld = find_short_path(examples[i].second);
    cout << "OLD: Solution: " << solOld << endl;
    if (solOld.size() != examples[i].first) {
      std::cout << "OLD: Wrong answer for map " << i << std::endl;
      fail++;
      cout << endl;
      continue;
    }
    
    std::cout << "NEW: Running test for map " << i << std::endl;
    auto solNew = find_long_path(examples[i].second);
    cout << "NEW: Solution: " << solNew << endl;
    if (solNew.size() != examples[i].first) {
      std::cout << "NEW: Wrong answer for map " << i << std::endl;
      fail++;
      cout << endl;
      continue;
    }

    if (solOld != solNew) {
      std::cout << "Wrong answer for map " << i << std::endl;
      fail++;
    }
    cout << endl;
  }

  if (fail) std::cout << "Failed " << fail << " tests" << std::endl;
  else std::cout << "All tests completed" << std::endl;

  //cout << "Starting memory test" << endl;
  //for (size_t i = 0; i < 1; i++) memoryTest();

  return 0;
}

#endif
