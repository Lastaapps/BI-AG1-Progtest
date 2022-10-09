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



LookupMap buildItemsLookup(const vector<vector<Point>>& src) {
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
  return lookup;
}


inline Mask createEndMask(const uint8_t items) {
  uint16_t m = 0;
  for (uint8_t i = 0; i < items; ++i) {
      m |= 1 << i;
  }
  return m;
}

vector<Point> backtrack(const ParentsMap& parents, const StatePoint end) {
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

vector<Point> smartBfs(const Graph& graph, const LookupMap& lookup, const Map& params) {
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

std::list<Place> find_path(const Map &map) {
  const Graph graph = buildGraph(map);
  const LookupMap lookup = buildItemsLookup(map.items);
  //for (const auto& edge : graph.edges) cout << edge.first << " -> " << edge.second << endl;
  
  const vector<Point> points = smartBfs(graph, lookup, map);
  list<Point> result;
  for (auto itr = points.rbegin(); itr != points.rend(); ++itr)
    result.emplace_back(*itr);
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
  size_t max = 10000;
  size_t itemsMax = 8;
  size_t maxInRoom = 5;
  vector<pair<Point, Point>> edges;
  vector<vector<Point>> items;

  for (size_t i = 0; i < 10 * max; ++i)
    edges.emplace_back(make_pair(rand() % 1000, rand() % 1000));
  for (size_t i = 0; i < itemsMax; ++i) {
    items.emplace_back(vector<Point>());
    for (size_t j = 0; j < maxInRoom; ++j)
      items[i].emplace_back(rand() % 1000);
  }
  const auto map = Map{max, 0, 1, edges, items};
  find_path(map);
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

  //cout << "Starting memory test" << endl;
  //for (size_t i = 0; i < 1; i++) memoryTest();

  return 0;
}

#endif
